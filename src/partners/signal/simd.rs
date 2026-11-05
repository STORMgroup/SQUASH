//! AVX2 kernels of the weighted edit distance's rows: the same values, directions, starts and choices as the scalar
//! loops, which stay the reference (randomized tests below) and run on CPUs without AVX2. A row's horizontal steps form
//! a min-plus prefix scan with a constant step cost: out[t] = min(v[t], out[t - 1] + gap), ties to v; per block of 8
//! cells a log-step scan (shifts 1, 2, 4) and the previous block's last value carried in.

static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();

/// AVX2 present.
pub fn on() -> bool {
    *ON.get_or_init(|| {
        #[cfg(target_arch = "x86_64")]
        {
            is_x86_feature_detected!("avx2")
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            false
        }
    })
}

#[cfg(test)]
/// `tracked`'s second pass with a constant horizontal cost (scalar reference): out[t] = min(v[t], out[t - 1]
/// + gap) with out[-1] = inf, ties to v; r[t] = 2 where the horizontal step wins; (first cheapest value, its
/// index; usize::MAX when none is below inf).
pub fn hscan_scalar(v: &[i32], gap: i32, inf: i32, out: &mut [i32], r: &mut [u8]) -> (i32, usize) {
    let (mut best_v, mut bk) = (inf, usize::MAX);
    let mut left = inf;
    for (t, ((o, &vv), rr)) in out.iter_mut().zip(v.iter()).zip(r.iter_mut()).enumerate() {
        let h = gap + left;
        let val = if h < vv {
            *rr = 2;
            h
        } else {
            vv
        };
        *o = val;
        left = val;
        if val < best_v {
            best_v = val;
            bk = t;
        }
    }
    (best_v, bk)
}

/// A semi-global row's horizontal pass with start tracking (`edit_seed_band`; scalar reference): out[j] = min(v[j],
/// out[j - 1] + gap) with out[-1] = left0 (ties to v); st[j] = st[j - 1] (st[-1] = st0) where the horizontal
/// step wins, else sv[j].
pub fn hscan_st_scalar(v: &[i32], sv: &[u32], gap: i32, left0: i32, st0: u32, out: &mut [i32], st: &mut [u32]) {
    let (mut left, mut ls) = (left0, st0);
    for j in 0..v.len() {
        let h = left + gap;
        if h < v[j] {
            out[j] = h;
            st[j] = ls;
        } else {
            out[j] = v[j];
            st[j] = sv[j];
        }
        left = out[j];
        ls = st[j];
    }
}

/// Byte masks of an 8-bit lane mask (0xff per set lane).
#[cfg(target_arch = "x86_64")]
static BYTE_MASKS: [u64; 256] = {
    let mut t = [0u64; 256];
    let mut m = 0;
    while m < 256 {
        let mut b = 0;
        let mut v = 0u64;
        while b < 8 {
            if (m >> b) & 1 == 1 {
                v |= 0xff << (8 * b);
            }
            b += 1;
        }
        t[m] = v;
        m += 1;
    }
    t
};

/// `hscan_scalar` with AVX2. Safety: AVX2 present; `out` and `r` as long as `v`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub unsafe fn hscan_avx2(v: &[i32], gap: i32, inf: i32, out: &mut [i32], r: &mut [u8]) -> (i32, usize) {
    use std::arch::x86_64::*;
    let n = v.len();
    assert!(out.len() >= n && r.len() >= n);
    let big = _mm256_set1_epi32(i32::MAX / 2);
    let g1 = _mm256_set1_epi32(gap);
    let g2 = _mm256_set1_epi32(2 * gap);
    let g4 = _mm256_set1_epi32(4 * gap);
    let ramp = _mm256_mullo_epi32(_mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8), g1);
    let idx1 = _mm256_setr_epi32(0, 0, 1, 2, 3, 4, 5, 6);
    let idx2 = _mm256_setr_epi32(0, 0, 0, 1, 2, 3, 4, 5);
    let idx4 = _mm256_setr_epi32(0, 0, 0, 0, 0, 1, 2, 3);
    let eight = _mm256_set1_epi32(8);
    let mut carry = inf;
    let mut bv = _mm256_set1_epi32(inf);
    let mut bi = _mm256_set1_epi32(i32::MAX);
    let mut tv = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    let mut t = 0usize;
    while t + 8 <= n {
        let vv = _mm256_loadu_si256(v.as_ptr().add(t) as *const __m256i);
        let mut x = vv;
        let s = _mm256_blend_epi32::<0b0000_0001>(_mm256_permutevar8x32_epi32(x, idx1), big);
        x = _mm256_min_epi32(x, _mm256_add_epi32(s, g1));
        let s = _mm256_blend_epi32::<0b0000_0011>(_mm256_permutevar8x32_epi32(x, idx2), big);
        x = _mm256_min_epi32(x, _mm256_add_epi32(s, g2));
        let s = _mm256_blend_epi32::<0b0000_1111>(_mm256_permutevar8x32_epi32(x, idx4), big);
        x = _mm256_min_epi32(x, _mm256_add_epi32(s, g4));
        x = _mm256_min_epi32(x, _mm256_add_epi32(_mm256_set1_epi32(carry), ramp));
        _mm256_storeu_si256(out.as_mut_ptr().add(t) as *mut __m256i, x);
        carry = _mm256_extract_epi32::<7>(x);
        let m = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(vv, x))) as usize;
        let p = r.as_mut_ptr().add(t) as *mut u64;
        let bm = BYTE_MASKS[m];
        p.write_unaligned((p.read_unaligned() & !bm) | (0x0202_0202_0202_0202 & bm));
        let lt = _mm256_cmpgt_epi32(bv, x);
        bv = _mm256_blendv_epi8(bv, x, lt);
        bi = _mm256_blendv_epi8(bi, tv, lt);
        tv = _mm256_add_epi32(tv, eight);
        t += 8;
    }
    let mut lv = [0i32; 8];
    let mut li = [0i32; 8];
    _mm256_storeu_si256(lv.as_mut_ptr() as *mut __m256i, bv);
    _mm256_storeu_si256(li.as_mut_ptr() as *mut __m256i, bi);
    let (mut best_v, mut bk) = (inf, usize::MAX);
    for l in 0..8 {
        if lv[l] < best_v || (lv[l] == best_v && lv[l] < inf && (li[l] as usize) < bk) {
            best_v = lv[l];
            bk = li[l] as usize;
        }
    }
    let mut left = carry;
    while t < n {
        let h = gap + left;
        let val = if h < v[t] {
            r[t] = 2;
            h
        } else {
            v[t]
        };
        out[t] = val;
        left = val;
        if val < best_v {
            best_v = val;
            bk = t;
        }
        t += 1;
    }
    (best_v, bk)
}

/// `hscan_st_scalar` with AVX2 (starts by a segmented scan: a lane whose horizontal step wins takes its left
/// neighbour's start). Safety: AVX2 present; `sv`, `out`, `st` as long as `v`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub unsafe fn hscan_st_avx2(v: &[i32], sv: &[u32], gap: i32, left0: i32, st0: u32, out: &mut [i32], st: &mut [u32]) {
    use std::arch::x86_64::*;
    let n = v.len();
    assert!(sv.len() >= n && out.len() >= n && st.len() >= n);
    let big = _mm256_set1_epi32(i32::MAX / 2);
    let ones = _mm256_set1_epi32(-1);
    let g1 = _mm256_set1_epi32(gap);
    let g2 = _mm256_set1_epi32(2 * gap);
    let g4 = _mm256_set1_epi32(4 * gap);
    let ramp = _mm256_mullo_epi32(_mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8), g1);
    let idx1 = _mm256_setr_epi32(0, 0, 1, 2, 3, 4, 5, 6);
    let idx2 = _mm256_setr_epi32(0, 0, 0, 1, 2, 3, 4, 5);
    let idx4 = _mm256_setr_epi32(0, 0, 0, 0, 0, 1, 2, 3);
    let (mut carry, mut cst) = (left0, st0);
    let mut j = 0usize;
    while j + 8 <= n {
        let vv = _mm256_loadu_si256(v.as_ptr().add(j) as *const __m256i);
        let mut x = vv;
        let s = _mm256_blend_epi32::<0b0000_0001>(_mm256_permutevar8x32_epi32(x, idx1), big);
        x = _mm256_min_epi32(x, _mm256_add_epi32(s, g1));
        let s = _mm256_blend_epi32::<0b0000_0011>(_mm256_permutevar8x32_epi32(x, idx2), big);
        x = _mm256_min_epi32(x, _mm256_add_epi32(s, g2));
        let s = _mm256_blend_epi32::<0b0000_1111>(_mm256_permutevar8x32_epi32(x, idx4), big);
        x = _mm256_min_epi32(x, _mm256_add_epi32(s, g4));
        x = _mm256_min_epi32(x, _mm256_add_epi32(_mm256_set1_epi32(carry), ramp));
        _mm256_storeu_si256(out.as_mut_ptr().add(j) as *mut __m256i, x);
        carry = _mm256_extract_epi32::<7>(x);
        let mut f = _mm256_cmpgt_epi32(vv, x);
        let mut s0 = _mm256_loadu_si256(sv.as_ptr().add(j) as *const __m256i);
        {
            let sh = _mm256_blend_epi32::<0b0000_0001>(_mm256_permutevar8x32_epi32(s0, idx1), ones);
            let fh = _mm256_blend_epi32::<0b0000_0001>(_mm256_permutevar8x32_epi32(f, idx1), ones);
            s0 = _mm256_blendv_epi8(s0, sh, f);
            f = _mm256_and_si256(f, fh);
            let sh = _mm256_blend_epi32::<0b0000_0011>(_mm256_permutevar8x32_epi32(s0, idx2), ones);
            let fh = _mm256_blend_epi32::<0b0000_0011>(_mm256_permutevar8x32_epi32(f, idx2), ones);
            s0 = _mm256_blendv_epi8(s0, sh, f);
            f = _mm256_and_si256(f, fh);
            let sh = _mm256_blend_epi32::<0b0000_1111>(_mm256_permutevar8x32_epi32(s0, idx4), ones);
            let fh = _mm256_blend_epi32::<0b0000_1111>(_mm256_permutevar8x32_epi32(f, idx4), ones);
            s0 = _mm256_blendv_epi8(s0, sh, f);
            f = _mm256_and_si256(f, fh);
            s0 = _mm256_blendv_epi8(s0, _mm256_set1_epi32(cst as i32), f);
        }
        _mm256_storeu_si256(st.as_mut_ptr().add(j) as *mut __m256i, s0);
        cst = _mm256_extract_epi32::<7>(s0) as u32;
        j += 8;
    }
    let (mut left, mut ls) = (carry, cst);
    while j < n {
        let h = left + gap;
        if h < v[j] {
            out[j] = h;
            st[j] = ls;
        } else {
            out[j] = v[j];
            st[j] = sv[j];
        }
        left = out[j];
        ls = st[j];
        j += 1;
    }
}

/// `hscan_st_scalar`, with AVX2 when `on()`.
pub fn hscan_st(v: &[i32], sv: &[u32], gap: i32, left0: i32, st0: u32, out: &mut [i32], st: &mut [u32]) {
    #[cfg(target_arch = "x86_64")]
    {
        if on() {
            // SAFETY: on() checked AVX2; the slices are as long as v (asserted inside)
            return unsafe { hscan_st_avx2(v, sv, gap, left0, st0, out, st) };
        }
    }
    hscan_st_scalar(v, sv, gap, left0, st0, out, st)
}
