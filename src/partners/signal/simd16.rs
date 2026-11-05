//! 16-bit lanes for the weighted edit distance's path alignment: `tracked` on q-symbols with 16 cells per AVX2
//! instruction instead of 8, each row in one fused pass (the diagonal / vertical candidates, the horizontal
//! min-plus scan, the directions and the cheapest cell; the 32-bit kernel makes two passes). Values are kept below
//! `INF` by subtracting the row's cheapest reachable value from the whole row whenever its largest reachable value
//! passes `rebase_at` (one constant per row keeps every comparison, so the same directions, cheapest cells and band
//! centres); unreachable cells stay pinned at `INF` by saturating arithmetic. A row whose reachable values then
//! still span `lim` or more returns None and the pair is aligned again with the 32-bit kernel, so every end and
//! path, hence every verdict, is the 32-bit kernel's. The directions of unreachable cells may differ (no path
//! crosses them). On with AVX2.

use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
/// alignments run with 16-bit lanes, and those that fell back to the 32-bit kernel
pub static RUNS: AtomicU64 = AtomicU64::new(0);
pub static FALLBACKS: AtomicU64 = AtomicU64::new(0);
pub const INF: i16 = i16::MAX;
/// a rebased row's reachable values must span less than this
pub const LIM: i16 = 16384;
/// a row is rebased when its largest reachable value passes this
pub const REBASE_AT: i16 = 24000;

/// AVX2 present (`simd::on()`).
pub fn on() -> bool {
    super::simd::on()
}

/// The alignments counted since the last call (16-bit runs, fallbacks).
pub fn take() -> (u64, u64) {
    (RUNS.swap(0, Relaxed), FALLBACKS.swap(0, Relaxed))
}

thread_local! {
    static BUF: std::cell::RefCell<(Vec<i16>, Vec<i16>)> = Default::default();
}

/// Frees this thread's 16-bit rows (`release_thread_buffers`).
pub fn release() {
    BUF.with(|b| *b.borrow_mut() = Default::default());
}

/// `tracked` (q-symbols, weighted edit distance with `gap`) with 16-bit lanes: Some(end) as the 32-bit kernel,
/// None to fall back (a rebased row's reachable values span `lim` or more, or the band too wide for 16-bit
/// values). Safety: AVX2 present.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
pub unsafe fn tracked16_avx2(x: &[u8], y: &[u8], band: i64, shift: i64, gap: i16, dir: &mut Vec<u8>,
                             lefts: &mut Vec<i64>, lim: i16, rebase_at: i16) -> Option<Option<i64>> {
    let width = (2 * band + 1) as usize;
    // a row adds at most 7 (substitution) + gap x width (horizontal steps) to its predecessors' values
    if width > 8000 || rebase_at as i64 + 7 + (gap as i64) * (width as i64 + 32) >= INF as i64 - 64 {
        return None;
    }
    RUNS.fetch_add(1, Relaxed);
    let r = BUF.with(|b| {
        let mut b = b.borrow_mut();
        let (pc, cc) = &mut *b;
        // SAFETY: AVX2 present (this function's contract)
        unsafe { body(x, y, band, shift, gap, dir, lefts, pc, cc, lim, rebase_at) }
    });
    if r.is_none() {
        FALLBACKS.fetch_add(1, Relaxed);
    }
    r
}

/// `tracked_body` (gap set) on 16-bit rows, one fused pass per row (`row16_avx2`).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
unsafe fn body(x: &[u8], y: &[u8], band: i64, shift: i64, gap: i16, dir: &mut Vec<u8>, lefts: &mut Vec<i64>,
               pc: &mut Vec<i16>, cc: &mut Vec<i16>, lim: i16, rebase_at: i16) -> Option<Option<i64>> {
    let n = x.len();
    let m = y.len() as i64;
    let width = (2 * band + 1) as usize;
    let pad = (shift + 2) as usize;
    dir.clear();
    dir.resize(n * width, 255);
    lefts.clear();
    pc.clear();
    pc.resize(width + 2 * pad, INF);
    cc.clear();
    cc.resize(width + 2 * pad, INF);
    let mut centre = 0i64;
    let mut pl = -band;
    lefts.push(pl);
    pc[pad + band as usize] = (x[0] as i16 - y[0] as i16).abs();
    dir[band as usize] = 3;
    let mut best_k = band as usize;
    for i in 1..n {
        let want = pl + best_k as i64 + 1;
        centre = want.clamp(centre + 1 - shift, centre + 1 + shift);
        let l = centre - band;
        lefts.push(l);
        let klo = (-l).max(0);
        let khi = (m - 1 - l).min(width as i64 - 1);
        if khi < klo {
            return Some(None);
        }
        let (klo, khi) = (klo as usize, khi as usize);
        for v in cc[..pad + klo].iter_mut() {
            *v = INF;
        }
        for v in cc[pad + khi + 1..].iter_mut() {
            *v = INF;
        }
        let off = (pad as i64 + l - pl) as usize;
        let yrow = &y[(l + klo as i64) as usize..=(l + khi as i64) as usize];
        let len = yrow.len();
        let row = &mut dir[i * width + klo..i * width + klo + len];
        let up = &pc[klo + off..klo + off + len];
        let dg = &pc[klo + off - 1..klo + off - 1 + len];
        let out = &mut cc[pad + klo..pad + klo + len];
        let (best_v, t, max_v) = row16_avx2(x[i] as i16, yrow, up, dg, gap, out, row);
        if !(best_v < INF) {
            return Some(None);
        }
        if max_v > rebase_at {
            if max_v - best_v >= lim {
                return None;
            }
            for v in out.iter_mut() {
                if *v != INF {
                    *v -= best_v;
                }
            }
        }
        best_k = klo + t;
        std::mem::swap(pc, cc);
        pl = l;
    }
    Some(Some(pl + best_k as i64))
}

/// One row of `tracked_body` (gap set) on 16-bit values with saturating steps, the reference of `row16_avx2`:
/// diagonal (substitution |x - y|) vs vertical (gap) candidates, ties diagonal (dir 1 / 0); then out[t] =
/// min(candidate, out[t - 1] + gap) with out[-1] = INF, ties to the candidate (dir 2 where the horizontal step wins).
/// Returns the first cheapest value and its index (usize::MAX when none is below INF) and the largest value below
/// INF (-1 when none).
#[allow(dead_code)] // the tests' reference
pub fn row16_scalar(xi: i16, yrow: &[u8], up: &[i16], dg: &[i16], gap: i16, out: &mut [i16], dir: &mut [u8])
                    -> (i16, usize, i16) {
    let (mut best_v, mut bk, mut max_v) = (INF, usize::MAX, -1i16);
    let mut left = INF;
    for t in 0..yrow.len() {
        let dv = dg[t].saturating_add((xi - yrow[t] as i16).abs());
        let uv = up[t].saturating_add(gap);
        let d = !(uv < dv);
        let v = if d { dv } else { uv };
        dir[t] = d as u8;
        let h = left.saturating_add(gap);
        let val = if h < v {
            dir[t] = 2;
            h
        } else {
            v
        };
        out[t] = val;
        left = val;
        if val < best_v {
            best_v = val;
            bk = t;
        }
        if val != INF && val > max_v {
            max_v = val;
        }
    }
    (best_v, bk, max_v)
}

/// The 16 bytes of a 16-lane mask (0 / -1 per lane), in lane order. Safety: AVX2 present.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn mask_bytes(m: std::arch::x86_64::__m256i) -> std::arch::x86_64::__m128i {
    use std::arch::x86_64::*;
    _mm256_castsi256_si128(_mm256_permute4x64_epi64::<0b11_01_10_00>(_mm256_packs_epi16(m, m)))
}

/// The constants of `row16_avx2`'s vector step.
#[cfg(target_arch = "x86_64")]
#[derive(Clone, Copy)]
struct C16 {
    inf: std::arch::x86_64::__m256i,
    g1: std::arch::x86_64::__m256i,
    g2: std::arch::x86_64::__m256i,
    g4: std::arch::x86_64::__m256i,
    g8: std::arch::x86_64::__m256i,
    ramp: std::arch::x86_64::__m256i,
    minus1: std::arch::x86_64::__m256i,
    sixteen: std::arch::x86_64::__m256i,
    xv: std::arch::x86_64::__m256i,
    ones8: std::arch::x86_64::__m128i,
    twos8: std::arch::x86_64::__m128i,
}

/// The state `row16_avx2` carries from vector to vector: the last cell's value, the cheapest value and its index
/// by lane, the largest reachable value by lane, the lanes' indices.
#[cfg(target_arch = "x86_64")]
struct S16 {
    carry: i16,
    bv: std::arch::x86_64::__m256i,
    bi: std::arch::x86_64::__m256i,
    mx: std::arch::x86_64::__m256i,
    tv: std::arch::x86_64::__m256i,
}

/// One vector (16 cells) of `row16_avx2` from y bytes `yp`, the previous row's `upp` / `dgp`, into `outp` / `dirp`;
/// only the lanes set in `valid` count for the cheapest cell and the largest value (the padded last vector).
/// Safety: AVX2 present; 16 readable / writable elements at each pointer.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn step16(c: &C16, st: &mut S16, yp: *const u8, upp: *const i16, dgp: *const i16, outp: *mut i16, dirp: *mut u8,
                 valid: std::arch::x86_64::__m256i) {
    use std::arch::x86_64::*;
    let yv = _mm256_cvtepu8_epi16(_mm_loadu_si128(yp as *const __m128i));
    let cost = _mm256_abs_epi16(_mm256_sub_epi16(c.xv, yv));
    let dv = _mm256_adds_epi16(_mm256_loadu_si256(dgp as *const __m256i), cost);
    let uv = _mm256_adds_epi16(_mm256_loadu_si256(upp as *const __m256i), c.g1);
    // the vertical step only where it is strictly cheaper (ties diagonal)
    let lt = _mm256_cmpgt_epi16(dv, uv);
    let vv = _mm256_blendv_epi8(dv, uv, lt);
    // x shifted up by k cells across the 128-bit halves, INF shifted in: alignr over [INF | low half]
    let mut x = vv;
    let lo = _mm256_permute2x128_si256::<0x02>(x, c.inf);
    x = _mm256_min_epi16(x, _mm256_adds_epi16(_mm256_alignr_epi8::<14>(x, lo), c.g1));
    let lo = _mm256_permute2x128_si256::<0x02>(x, c.inf);
    x = _mm256_min_epi16(x, _mm256_adds_epi16(_mm256_alignr_epi8::<12>(x, lo), c.g2));
    let lo = _mm256_permute2x128_si256::<0x02>(x, c.inf);
    x = _mm256_min_epi16(x, _mm256_adds_epi16(_mm256_alignr_epi8::<8>(x, lo), c.g4));
    x = _mm256_min_epi16(x, _mm256_adds_epi16(_mm256_permute2x128_si256::<0x02>(x, c.inf), c.g8));
    x = _mm256_min_epi16(x, _mm256_adds_epi16(_mm256_set1_epi16(st.carry), c.ramp));
    _mm256_storeu_si256(outp as *mut __m256i, x);
    st.carry = _mm256_extract_epi16::<15>(x) as i16;
    // directions: 2 where the horizontal step won (out < candidate), else 0 (vertical) / 1 (diagonal)
    let hz = _mm256_cmpgt_epi16(vv, x);
    let db = _mm_blendv_epi8(_mm_andnot_si128(mask_bytes(lt), c.ones8), c.twos8, mask_bytes(hz));
    _mm_storeu_si128(dirp as *mut __m128i, db);
    let xs = _mm256_blendv_epi8(c.inf, x, valid);
    let less = _mm256_cmpgt_epi16(st.bv, xs);
    st.bv = _mm256_blendv_epi8(st.bv, xs, less);
    st.bi = _mm256_blendv_epi8(st.bi, st.tv, less);
    st.mx = _mm256_max_epi16(st.mx, _mm256_blendv_epi8(xs, c.minus1, _mm256_cmpeq_epi16(xs, c.inf)));
    st.tv = _mm256_add_epi16(st.tv, c.sixteen);
}

/// Horizontal min / max of 16 lanes. Safety: AVX2 present.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn hmin16(v: std::arch::x86_64::__m256i) -> i16 {
    use std::arch::x86_64::*;
    let mut m = _mm_min_epi16(_mm256_castsi256_si128(v), _mm256_extracti128_si256::<1>(v));
    m = _mm_min_epi16(m, _mm_srli_si128::<8>(m));
    m = _mm_min_epi16(m, _mm_srli_si128::<4>(m));
    m = _mm_min_epi16(m, _mm_srli_si128::<2>(m));
    _mm_extract_epi16::<0>(m) as i16
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn hmax16(v: std::arch::x86_64::__m256i) -> i16 {
    use std::arch::x86_64::*;
    let mut m = _mm_max_epi16(_mm256_castsi256_si128(v), _mm256_extracti128_si256::<1>(v));
    m = _mm_max_epi16(m, _mm_srli_si128::<8>(m));
    m = _mm_max_epi16(m, _mm_srli_si128::<4>(m));
    m = _mm_max_epi16(m, _mm_srli_si128::<2>(m));
    _mm_extract_epi16::<0>(m) as i16
}

/// `row16_scalar` with AVX2, 16 cells per vector (`step16`): the candidates, a log-step min-plus prefix scan within
/// each vector (INF shifted in) with the carry from the previous vector, the directions as bytes, the cheapest
/// cell and the largest reachable value by lane, reduced at the end. The last partial vector runs on INF-padded
/// copies (the padding lanes come after the row's cells, so they change none of them and are never cheapest).
/// Safety: AVX2 present; `up`, `dg`, `out`, `dir` as long as `yrow`; yrow.len() < 32767.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
pub unsafe fn row16_avx2(xi: i16, yrow: &[u8], up: &[i16], dg: &[i16], gap: i16, out: &mut [i16], dir: &mut [u8])
                         -> (i16, usize, i16) {
    use std::arch::x86_64::*;
    let n = yrow.len();
    assert!(up.len() >= n && dg.len() >= n && out.len() >= n && dir.len() >= n && n < i16::MAX as usize);
    let g1 = _mm256_set1_epi16(gap);
    let c = C16 {
        inf: _mm256_set1_epi16(INF),
        g1,
        g2: _mm256_set1_epi16(2 * gap),
        g4: _mm256_set1_epi16(4 * gap),
        g8: _mm256_set1_epi16(8 * gap),
        ramp: _mm256_mullo_epi16(_mm256_setr_epi16(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), g1),
        minus1: _mm256_set1_epi16(-1),
        sixteen: _mm256_set1_epi16(16),
        xv: _mm256_set1_epi16(xi),
        ones8: _mm_set1_epi8(1),
        twos8: _mm_set1_epi8(2),
    };
    let mut st = S16 {
        carry: INF,
        bv: c.inf,
        bi: _mm256_set1_epi16(i16::MAX),
        mx: c.minus1,
        tv: _mm256_setr_epi16(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
    };
    let mut t = 0usize;
    while t + 16 <= n {
        step16(&c, &mut st, yrow.as_ptr().add(t), up.as_ptr().add(t), dg.as_ptr().add(t), out.as_mut_ptr().add(t),
               dir.as_mut_ptr().add(t), c.minus1);
        t += 16;
    }
    if t < n {
        let rem = n - t;
        let (mut yb, mut ub, mut db) = ([0u8; 16], [INF; 16], [INF; 16]);
        let (mut ob, mut rb) = ([INF; 16], [0u8; 16]);
        yb[..rem].copy_from_slice(&yrow[t..]);
        ub[..rem].copy_from_slice(&up[t..n]);
        db[..rem].copy_from_slice(&dg[t..n]);
        // the first rem lanes
        let lane = _mm256_setr_epi16(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
        let valid = _mm256_cmpgt_epi16(_mm256_set1_epi16(rem as i16), lane);
        step16(&c, &mut st, yb.as_ptr(), ub.as_ptr(), db.as_ptr(), ob.as_mut_ptr(), rb.as_mut_ptr(), valid);
        out[t..n].copy_from_slice(&ob[..rem]);
        dir[t..n].copy_from_slice(&rb[..rem]);
    }
    let max_v = hmax16(st.mx);
    let best_v = hmin16(st.bv);
    if best_v == INF {
        return (INF, usize::MAX, max_v);
    }
    let at = _mm256_blendv_epi8(_mm256_set1_epi16(i16::MAX), st.bi, _mm256_cmpeq_epi16(st.bv, _mm256_set1_epi16(best_v)));
    (best_v, hmin16(at) as usize, max_v)
}
