//! The token alignment of a donor window to a target window: the overlap DP (+2 / -2 / -3, free end gaps within a
//! margin), its banded form, and the per-thread DP buffers.

const GAP: i64 = -3;

/// Reusable DP buffers, one set per thread. Every trace cell an alignment reads is written by
/// that alignment first, so the buffers are never cleared.
#[derive(Default)]
struct DpScratch {
    trace: Vec<u8>,
    row: Vec<i32>,
    next: Vec<i32>,
    last_col: Vec<i32>,
}

thread_local! {
    static DP_SCRATCH: std::cell::RefCell<DpScratch> = std::cell::RefCell::new(DpScratch::default());
}

/// Overlap alignment of two token windows: +2/-2/-3, free end gaps (within `margin`
/// tokens of each window end when given), diag > target-only > donor-only.
///
/// Each row is filled in two passes with the recurrence, scores and tie-breaks of the one-pass
/// loop (`align_overlap_reference` in the tests): first max(diag, up), diag winning ties, which
/// has no dependency along the row; then the donor-only moves as a running max, taken only when
/// strictly better. Scores are i32 (|score| <= 3 (n + m)) and the buffers are reused per thread.
pub fn align_overlap(target: &[u32], donor: &[u32], margin: Option<usize>) -> Vec<(Option<usize>, Option<usize>)> {
    DP_SCRATCH.with(|s| align_full(target, donor, margin, &mut s.borrow_mut()))
}

fn align_full(target: &[u32], donor: &[u32], margin: Option<usize>, s: &mut DpScratch) -> Vec<(Option<usize>, Option<usize>)> {
    let g = GAP as i32;
    let (n, m) = (target.len(), donor.len());
    assert!(n + m < 1 << 28, "alignment too long for i32 scores");
    let free_n = margin.map_or(n, |x| x.min(n));
    let free_m = margin.map_or(m, |x| x.min(m));
    if s.trace.len() < n * m {
        let before = s.trace.capacity();
        s.trace.resize(n * m, 0);
        crate::partners::mem::grew(before, s.trace.capacity());
    }
    let DpScratch { trace, row, next, last_col } = s;
    // Trace of cell (i, j), i, j >= 1, at (i - 1) * m + j - 1: 0 diag, 1 target-only (up), 2 donor-only (left).
    row.clear();
    row.extend((0..=m).map(|j| if j > free_m { (j - free_m) as i32 * g } else { 0 }));
    next.clear();
    next.resize(m + 1, 0);
    last_col.clear();
    last_col.resize(n + 1, 0);
    last_col[0] = row[m];
    for i in 1..=n {
        next[0] = if i > free_n { row[0] + g } else { 0 };
        let left = next[0];
        dp_row(target[i - 1], g, &row[..m], &row[1..=m], &donor[..m], left, &mut next[1..=m], &mut trace[(i - 1) * m..i * m]);
        std::mem::swap(row, next);
        last_col[i] = row[m];
    }
    // Best end on the last row or column: (score, i + j, i); keys are unique.
    let mut end: Option<((i32, usize, usize), usize, usize)> = None;
    let mut consider = |key: (i32, usize, usize), i: usize, j: usize| {
        if end.map_or(true, |(k, _, _)| key > k) {
            end = Some((key, i, j));
        }
    };
    for j in 0..=m {
        let s = row[j] + g * m.saturating_sub(j + free_m) as i32;
        consider((s, n + j, n), n, j);
    }
    for i in 0..n {
        let s = last_col[i] + g * n.saturating_sub(i + free_n) as i32;
        consider((s, i + m, i), i, m);
    }
    let (_, i, j) = end.unwrap();
    let path = traceback(n, m, i, j, |i, j| trace[(i - 1) * m + j - 1]);
    release_trace(trace);
    path
}

/// One DP row: cells j = 1..=len of `out` from the previous row's cells j - 1 (`diag_in`) and j (`up_in`), the donor
/// tokens j - 1, and the cell before the first (`left`). First max(diag, up), diag winning ties; then the donor-only
/// moves as a running max, taken only when strictly better (trace 0 diag, 1 up, 2 left). The slices of one length
/// let the first pass compile to SIMD (AVX2 when the CPU has it).
#[inline]
fn dp_row(t: u32, g: i32, diag_in: &[i32], up_in: &[i32], donor: &[u32], left: i32, out: &mut [i32], tr: &mut [u8]) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx2") {
            // SAFETY: the CPU supports AVX2 (checked above).
            return unsafe { dp_row_avx2(t, g, diag_in, up_in, donor, left, out, tr) };
        }
    }
    dp_row_body(t, g, diag_in, up_in, donor, left, out, tr)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn dp_row_avx2(t: u32, g: i32, diag_in: &[i32], up_in: &[i32], donor: &[u32], left: i32, out: &mut [i32], tr: &mut [u8]) {
    dp_row_body(t, g, diag_in, up_in, donor, left, out, tr)
}

#[inline(always)]
fn dp_row_body(t: u32, g: i32, diag_in: &[i32], up_in: &[i32], donor: &[u32], left: i32, out: &mut [i32], tr: &mut [u8]) {
    let len = out.len();
    let (diag_in, up_in, donor, tr) = (&diag_in[..len], &up_in[..len], &donor[..len], &mut tr[..len]);
    for k in 0..len {
        let diag = diag_in[k] + if t == donor[k] { 2 } else { -2 };
        let up = up_in[k] + g;
        let take_up = up > diag;
        out[k] = if take_up { up } else { diag };
        tr[k] = take_up as u8;
    }
    let mut left = left;
    for k in 0..len {
        let l = left + g;
        if l > out[k] {
            out[k] = l;
            tr[k] = 2;
        }
        left = out[k];
    }
}

/// A thread keeps its DP trace buffer for the next alignment, but not one larger than `TRACE_KEEP` bytes (a rare full
/// DP of two long windows), which is freed after use.
const TRACE_KEEP: usize = 64 << 20;

fn release_trace(trace: &mut Vec<u8>) {
    if trace.capacity() > TRACE_KEEP {
        let before = trace.capacity();
        *trace = Vec::new();
        crate::partners::mem::grew(before, 0);
    }
}

/// Free every pool thread's DP buffers once a round's replay is done; the partner search does not use them, and the
/// next replay grows them again.
pub(super) fn release_dp_scratch(pool: &rayon::ThreadPool) {
    pool.broadcast(|_| {
        DP_SCRATCH.with(|s| {
            let mut s = s.borrow_mut();
            crate::partners::mem::grew(s.trace.capacity(), 0);
            *s = DpScratch::default();
        })
    });
}

/// The alignment from the chosen end (i, j): overhang after it, the traced path, overhang before it.
fn traceback(n: usize, m: usize, mut i: usize, mut j: usize, op: impl Fn(usize, usize) -> u8) -> Vec<(Option<usize>, Option<usize>)> {
    let mut tail: Vec<(Option<usize>, Option<usize>)> = (i..n).map(|a| (Some(a), None)).collect();
    tail.extend((j..m).map(|b| (None, Some(b))));
    let mut result = Vec::with_capacity(n + m);
    while i > 0 && j > 0 {
        match op(i, j) {
            0 => {
                result.push((Some(i - 1), Some(j - 1)));
                i -= 1;
                j -= 1;
            }
            1 => {
                result.push((Some(i - 1), None));
                i -= 1;
            }
            _ => {
                result.push((None, Some(j - 1)));
                j -= 1;
            }
        }
    }
    result.extend((0..i).rev().map(|a| (Some(a), None)));
    result.extend((0..j).rev().map(|b| (None, Some(b))));
    result.reverse();
    result.extend(tail);
    result
}

/// Centre of row i's band: the window's corner-to-corner line, round(i m / n).
fn band_center(i: usize, n: usize, m: usize) -> usize {
    (i * m + n / 2) / n
}

/// `rounds --band W`: `align_overlap` restricted to cells within W of the window's
/// corner-to-corner line (`band_center`), cells outside unreachable, with the same scores and
/// tie-breaks inside. When the banded path comes within W/4 of the band's edge, the band does not
/// connect the window's corners (a window much wider than tall), or the band covers the donor
/// anyway, the full DP is used. The result can differ from the full DP's when the full DP's best
/// path leaves the band.
pub fn align_overlap_banded(target: &[u32], donor: &[u32], margin: Option<usize>, band: usize) -> Vec<(Option<usize>, Option<usize>)> {
    let (n, m) = (target.len(), donor.len());
    if band == 0 || n == 0 || 2 * band >= m {
        return align_overlap(target, donor, margin);
    }
    let banded = match DP_SCRATCH.with(|s| align_band(target, donor, margin, band, &mut s.borrow_mut())) {
        Some(x) => x,
        None => return align_overlap(target, donor, margin),
    };
    let limit = band - band / 4;
    let deviation = banded
        .iter()
        .filter_map(|&(a, b)| Some((a?, b?)))
        .map(|(x, y)| {
            let c = band_center(x + 1, n, m);
            if y + 1 > c { y + 1 - c } else { c - (y + 1) }
        })
        .max()
        .unwrap_or(0);
    if deviation > limit {
        align_overlap(target, donor, margin)
    } else {
        banded
    }
}

/// The banded DP; None when no end cell is reachable inside the band.
fn align_band(target: &[u32], donor: &[u32], margin: Option<usize>, w: usize, s: &mut DpScratch) -> Option<Vec<(Option<usize>, Option<usize>)>> {
    const NEG: i32 = i32::MIN / 4;
    let g = GAP as i32;
    let (n, m) = (target.len(), donor.len());
    assert!(n + m < 1 << 28, "alignment too long for i32 scores");
    let free_n = margin.map_or(n, |x| x.min(n));
    let free_m = margin.map_or(m, |x| x.min(m));
    let width = 2 * w + 1;
    let bounds = |i: usize| {
        let c = band_center(i, n, m);
        (c.saturating_sub(w), (c + w).min(m))
    };
    if s.trace.len() < n * width {
        let before = s.trace.capacity();
        s.trace.resize(n * width, 0);
        crate::partners::mem::grew(before, s.trace.capacity());
    }
    let DpScratch { trace, row, next, last_col } = s;
    // Row i holds its band [lo, hi]; the cells just outside it that the next row reads are NEG.
    // Trace of cell (i, j) at (i - 1) * width + j - lo(i).
    let (_, hi0) = bounds(0);
    row.clear();
    row.extend((0..=m).map(|j| if j > hi0 { NEG } else if j > free_m { (j - free_m) as i32 * g } else { 0 }));
    next.clear();
    next.resize(m + 1, NEG);
    last_col.clear();
    last_col.resize(n + 1, NEG);
    last_col[0] = row[m];
    for i in 1..=n {
        let (lo, hi) = bounds(i);
        if lo == 0 {
            next[0] = if i > free_n { row[0] + g } else { 0 };
        } else {
            next[lo - 1] = NEG;
        }
        let first = lo.max(1);
        if first <= hi {
            let left = next[first - 1];
            let tr = &mut trace[(i - 1) * width..i * width];
            dp_row(target[i - 1], g, &row[first - 1..hi], &row[first..=hi], &donor[first - 1..hi], left, &mut next[first..=hi],
                   &mut tr[first - lo..=hi - lo]);
        }
        if i < n {
            let (_, hi_next) = bounds(i + 1);
            for v in next.iter_mut().take(hi_next + 1).skip(hi + 1) {
                *v = NEG;
            }
        }
        std::mem::swap(row, next);
        last_col[i] = if hi == m { row[m] } else { NEG };
    }
    let mut end: Option<((i32, usize, usize), usize, usize)> = None;
    let mut consider = |key: (i32, usize, usize), i: usize, j: usize| {
        if key.0 > NEG / 2 && end.map_or(true, |(k, _, _)| key > k) {
            end = Some((key, i, j));
        }
    };
    let (lo_n, _) = bounds(n);
    for j in lo_n..=m {
        let s = row[j] + g * m.saturating_sub(j + free_m) as i32;
        consider((s, n + j, n), n, j);
    }
    for i in 0..n {
        let s = last_col[i] + g * n.saturating_sub(i + free_n) as i32;
        consider((s, i + m, i), i, m);
    }
    let (_, i, j) = match end {
        Some(x) => x,
        None => {
            release_trace(trace);
            return None;
        }
    };
    let path = traceback(n, m, i, j, |i, j| trace[(i - 1) * width + j - bounds(i).0]);
    release_trace(trace);
    Some(path)
}

/// `aligned_core`: the alignment between its first and last paired column.
pub fn aligned_core(alignment: &[(Option<usize>, Option<usize>)]) -> &[(Option<usize>, Option<usize>)] {
    let paired = |x: &&(Option<usize>, Option<usize>)| x.0.is_some() && x.1.is_some();
    match (alignment.iter().position(|x| paired(&x)), alignment.iter().rposition(|x| paired(&x))) {
        (Some(a), Some(b)) => &alignment[a..=b],
        _ => &[],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn core(target: &[u32], donor: &[u32], margin: Option<usize>) -> Vec<(Option<usize>, Option<usize>)> {
        aligned_core(&align_overlap(target, donor, margin)).to_vec()
    }

    /// A one-pass `align_overlap`: the reference of the two-pass DP.
    fn align_overlap_reference(target: &[u32], donor: &[u32], margin: Option<usize>) -> Vec<(Option<usize>, Option<usize>)> {
        let (n, m) = (target.len(), donor.len());
        let free_n = margin.map_or(n, |x| x.min(n));
        let free_m = margin.map_or(m, |x| x.min(m));
        let w = m + 1;
        // Trace: 0 diag, 1 target-only (up), 2 donor-only (left).
        let mut trace = vec![0u8; (n + 1) * w];
        let mut row: Vec<i64> = (0..=m).map(|j| if j > free_m { (j - free_m) as i64 * GAP } else { 0 }).collect();
        for j in 1..=m {
            trace[j] = 2;
        }
        let mut last_col = vec![0i64; n + 1];
        last_col[0] = row[m];
        let mut next = vec![0i64; w];
        for i in 1..=n {
            trace[i * w] = 1;
            next[0] = if i > free_n { row[0] + GAP } else { 0 };
            let t = target[i - 1];
            for j in 1..=m {
                let diag = row[j - 1] + if t == donor[j - 1] { 2 } else { -2 };
                let up = row[j] + GAP;
                let left = next[j - 1] + GAP;
                let (mut best, mut op) = (diag, 0u8);
                if up > best {
                    best = up;
                    op = 1;
                }
                if left > best {
                    best = left;
                    op = 2;
                }
                next[j] = best;
                trace[i * w + j] = op;
            }
            std::mem::swap(&mut row, &mut next);
            last_col[i] = row[m];
        }
        // Best end on the last row or column: (score, i + j, i); keys are unique.
        let mut end: Option<((i64, usize, usize), usize, usize)> = None;
        let mut consider = |key: (i64, usize, usize), i: usize, j: usize| {
            if end.map_or(true, |(k, _, _)| key > k) {
                end = Some((key, i, j));
            }
        };
        for j in 0..=m {
            let s = row[j] + GAP * m.saturating_sub(j + free_m) as i64;
            consider((s, n + j, n), n, j);
        }
        for i in 0..n {
            let s = last_col[i] + GAP * n.saturating_sub(i + free_n) as i64;
            consider((s, i + m, i), i, m);
        }
        let (_, mut i, mut j) = end.unwrap();
        let mut tail: Vec<(Option<usize>, Option<usize>)> = (i..n).map(|a| (Some(a), None)).collect();
        tail.extend((j..m).map(|b| (None, Some(b))));
        let mut result = Vec::with_capacity(n + m);
        while i > 0 && j > 0 {
            match trace[i * w + j] {
                0 => {
                    result.push((Some(i - 1), Some(j - 1)));
                    i -= 1;
                    j -= 1;
                }
                1 => {
                    result.push((Some(i - 1), None));
                    i -= 1;
                }
                _ => {
                    result.push((None, Some(j - 1)));
                    j -= 1;
                }
            }
        }
        result.extend((0..i).rev().map(|a| (Some(a), None)));
        result.extend((0..j).rev().map(|b| (None, Some(b))));
        result.reverse();
        result.extend(tail);
        result
    }

    /// Deterministic pairs: a random target over a small alphabet (many ties) and a donor made from
    /// it with substitutions, insertions, deletions and optional end overhangs.
    fn related_pairs(count: usize, seed: u64) -> Vec<(Vec<u32>, Vec<u32>)> {
        let mut x = seed;
        let mut next = move |k: u64| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (x >> 33) % k
        };
        (0..count)
            .map(|c| {
                let n = next(if c % 7 == 0 { 3 } else { 400 }) as usize;
                let alphabet = 2 + next(6);
                let target: Vec<u32> = (0..n).map(|_| next(alphabet) as u32).collect();
                let mut donor: Vec<u32> = (0..next(30)).map(|_| next(alphabet) as u32).collect();
                for &t in &target {
                    match next(20) {
                        0 => {}
                        1 => donor.extend([t, next(alphabet) as u32]),
                        2 => donor.push(next(alphabet) as u32),
                        _ => donor.push(t),
                    }
                }
                donor.extend((0..next(30)).map(|_| next(alphabet) as u32));
                (target, donor)
            })
            .collect()
    }

    #[test]
    fn two_pass_dp_matches_the_one_pass_reference() {
        let mut cases = related_pairs(600, 7);
        cases.extend([(vec![], vec![1, 2]), (vec![1, 2], vec![]), (vec![], vec![]), (vec![5], vec![5])]);
        for (case, (target, donor)) in cases.iter().enumerate() {
            for margin in [None, Some(0), Some(2), Some(16)] {
                assert_eq!(
                    align_overlap(target, donor, margin),
                    align_overlap_reference(target, donor, margin),
                    "case {} margin {:?}",
                    case,
                    margin
                );
            }
        }
    }

    #[test]
    fn banded_dp_matches_the_full_dp_near_the_diagonal_and_falls_back_off_it() {
        let cases = related_pairs(400, 11);
        let mut banded_used = 0;
        for (target, donor) in &cases {
            let full = align_overlap(target, donor, Some(16));
            // Band 0 and bands that cover the donor are the full DP.
            assert_eq!(align_overlap_banded(target, donor, Some(16), 0), full);
            assert_eq!(align_overlap_banded(target, donor, Some(16), donor.len()), full);
            let banded = align_overlap_banded(target, donor, Some(16), 64);
            if donor.len() > 128 {
                banded_used += 1;
            }
            assert_eq!(banded, full);
        }
        assert!(banded_used > 200, "{}", banded_used);
        // A donor shifted by 150 tokens: the banded path would run along the band's edge, so the full DP is used.
        let target: Vec<u32> = (0..600u64).map(|i| (i.wrapping_mul(2654435761) >> 7) as u32 % 100_000).collect();
        let donor: Vec<u32> = (0..150).map(|i| 200_000 + i).chain(target[..450].iter().copied()).collect();
        let full = align_overlap(&target, &donor, None);
        assert_eq!(align_overlap_banded(&target, &donor, None, 32), full);
        let raw = DP_SCRATCH.with(|s| align_band(&target, &donor, None, 32, &mut s.borrow_mut()));
        assert_ne!(raw, Some(full));
        // A window far wider than tall: the band cannot connect its corners, so the full DP is used.
        let (short, wide): (Vec<u32>, Vec<u32>) = ((0..5).collect(), (0..200).map(|x| x % 7).collect());
        assert_eq!(DP_SCRATCH.with(|s| align_band(&short, &wide, Some(16), 4, &mut s.borrow_mut())), None);
        assert_eq!(align_overlap_banded(&short, &wide, Some(16), 4), align_overlap(&short, &wide, Some(16)));
    }

    #[test]
    fn overlap_alignment_leaves_overhangs_unaligned() {
        let full = align_overlap(&[1, 2, 3, 4, 5], &[3, 4, 5, 6, 7], None);
        assert_eq!(aligned_core(&full), &[(Some(2), Some(0)), (Some(3), Some(1)), (Some(4), Some(2))]);
        // The alignment holds every token of both sequences: overhangs before and after the core.
        assert_eq!(
            full,
            vec![(Some(0), None), (Some(1), None), (Some(2), Some(0)), (Some(3), Some(1)), (Some(4), Some(2)), (None, Some(3)), (None, Some(4))]
        );
    }

    #[test]
    fn anchored_alignment_keeps_a_noisy_certified_overlap() {
        let (target, donor) = ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [1, 20, 30, 4, 50, 60, 7, 80, 90, 10]);
        assert!(core(&target, &donor, None).is_empty());
        assert_eq!(core(&target, &donor, Some(2)).len(), 10);
        let expected: Vec<_> = (0..5).map(|k| (Some(k), Some(k))).collect();
        assert_eq!(core(&[1, 2, 3, 4, 5, 6, 7, 8], &[1, 2, 3, 4, 5], Some(4)), expected);
    }
}
