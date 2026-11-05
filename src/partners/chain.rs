//! Chaining a donor's anchors: the colinear chaining DP (score ring, one-byte back links), its members, and the
//! rounding of the chain's offset.

use super::*;

/// An anchor (target token index, donor token index): (i64, i64), or in the search's per-thread buffer (u32, u32)
/// (8 bytes) or, when the target and the donor have at most 65,536 tokens, (u16, u16) (4 bytes; chosen per donor in
/// `target_chains`). Ordered as the pair.
pub trait Anchor: Copy + Ord {
    fn t(&self) -> i64;
    fn d(&self) -> i64;
}

impl Anchor for (i64, i64) {
    fn t(&self) -> i64 {
        self.0
    }
    fn d(&self) -> i64 {
        self.1
    }
}

impl Anchor for (u32, u32) {
    fn t(&self) -> i64 {
        self.0 as i64
    }
    fn d(&self) -> i64 {
        self.1 as i64
    }
}

impl Anchor for (u16, u16) {
    fn t(&self) -> i64 {
        self.0 as i64
    }
    fn d(&self) -> i64 {
        self.1 as i64
    }
}


/// `chain_anchors`: best forward colinear chain -> (anchors, score, t0, t1, d0, d1).
pub fn chain_anchors(mut anchors: Vec<(i64, i64)>, p: &Params) -> Option<(usize, f64, i64, i64, i64, i64)> {
    anchors.sort_unstable();
    anchors.dedup();
    chain_sorted(&anchors, p)
}

thread_local! {
    /// Per-thread scratch of `chain_sorted` (score, back), reused from one donor to the next, so that chaining a donor
    /// allocates nothing. There is no per-anchor count (a chain's anchor count is its back links plus one); the back
    /// links are one byte (`back8`) when the look-back fits one (`chain_back8`), u32 otherwise, and the one not in use
    /// is empty (a target-donor pair inside a tandem array can have tens of millions of anchors). The scratch keeps the
    /// last group's DP for `chain_members_after_sorted` (its anchors' address and length, and the chain's end).
    pub(super) static CHAIN_SCRATCH: std::cell::RefCell<ChainScratch> = std::cell::RefCell::new(ChainScratch::default());
}

/// `CHAIN_SCRATCH`: the scores (a ring, `chain_dp_ring`), the back links (u32 or u8) and the last group chained.
#[derive(Default)]
pub(super) struct ChainScratch {
    pub(super) score: Vec<f64>,
    pub(super) back: Vec<u32>,
    pub(super) back8: Vec<u8>,
    pub(super) last: (usize, usize, usize),
}

impl ChainScratch {
    pub(super) fn bytes(&self) -> usize {
        self.score.capacity() * 8 + self.back.capacity() * 4 + self.back8.capacity()
    }
}

/// `chain_anchors` on anchors already sorted and distinct.
pub(super) fn chain_sorted<A: Anchor>(anchors: &[A], p: &Params) -> Option<(usize, f64, i64, i64, i64, i64)> {
    let n = anchors.len();
    if n == 0 {
        return None;
    }
    CHAIN_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        let s = &mut *scratch;
        let (out, end) = if chain_back8(p) {
            s.back.clear();
            chain_into(anchors, p, &mut s.score, &mut s.back8)
        } else {
            s.back8.clear();
            chain_into(anchors, p, &mut s.score, &mut s.back)
        };
        if memtrack::on() {
            memtrack::set(memtrack::CHAIN, s.bytes());
        }
        s.last = (anchors.as_ptr() as usize, n, end);
        out
    })
}

/// `chain_sorted`'s DP on its scratch: the scores in a ring (`chain_dp_ring`), `back` sized to the group.
fn chain_into<A: Anchor, B: BackLink>(
    anchors: &[A],
    p: &Params,
    score: &mut Vec<f64>,
    back: &mut Vec<B>,
) -> (Option<(usize, f64, i64, i64, i64, i64)>, usize) {
    let n = anchors.len();
    back.clear();
    back.resize(n, B::NONE);
    score.clear();
    // the scores in a ring of a power of two >= lookback + 1 (`chain_dp_ring`)
    score.resize((p.lookback + 1).next_power_of_two().min(n.next_power_of_two()), 1.0f64);
    chain_dp_ring(anchors, p, score, back)
}

/// `chain_dp`'s back link of an anchor that starts its chain.
const NO_BACK: u32 = u32::MAX;

/// A back link of `chain_dp` (the anchor before `j` in its best chain, or none): u32 = that anchor's index (`NO_BACK`
/// for none); u8 = how far back it is, 0 for none. A predecessor is never more than `lookback` anchors back (every
/// look-back starts at `j - lookback`), so with a look-back up to 255 one byte holds it: 1 byte of chaining scratch an
/// anchor instead of 4 (with the score ring), the same chains and members (`chain_ring_equals_the_full_dp`).
trait BackLink: Copy + PartialEq {
    const NONE: Self;
    /// The link of anchor `j` to anchor `i < j`.
    fn link(j: usize, i: usize) -> Self;
    /// The anchor this link of anchor `j` points to (not `NONE`).
    fn target(self, j: usize) -> usize;
}

impl BackLink for u32 {
    const NONE: u32 = NO_BACK;
    fn link(_j: usize, i: usize) -> u32 {
        i as u32
    }
    fn target(self, _j: usize) -> usize {
        self as usize
    }
}

impl BackLink for u8 {
    const NONE: u8 = 0;
    fn link(j: usize, i: usize) -> u8 {
        debug_assert!(i < j && j - i <= 255);
        (j - i) as u8
    }
    fn target(self, j: usize) -> usize {
        j - self as usize
    }
}

/// `chain_sorted` keeps one-byte back links when the look-back fits one (u32 links otherwise).
fn chain_back8(p: &Params) -> bool {
    p.lookback <= 255
}

/// The anchors of the chain ending at `end`, by its back links, in chain order.
fn chain_trace<A: Anchor, B: BackLink>(anchors: &[A], back: &[B], end: usize) -> Vec<(i64, i64)> {
    let mut k = end;
    let mut out = vec![(anchors[k].t(), anchors[k].d())];
    while back[k] != B::NONE {
        k = back[k].target(k);
        out.push((anchors[k].t(), anchors[k].d()));
    }
    out.reverse();
    out
}

/// `chain_members` of the group `chain_sorted` chained last on this thread (sorted and distinct), from the back
/// links it kept instead of a second DP; None when the thread's last group is another.
pub(super) fn chain_members_after_sorted<A: Anchor>(anchors: &[A]) -> Option<Vec<(i64, i64)>> {
    CHAIN_SCRATCH.with(|cell| {
        let s = cell.borrow();
        let n = anchors.len();
        if n == 0 || s.last != (anchors.as_ptr() as usize, n, s.last.2) {
            return None;
        }
        if s.back8.len() == n {
            Some(chain_trace(anchors, &s.back8, s.last.2))
        } else if s.back.len() == n {
            Some(chain_trace(anchors, &s.back, s.last.2))
        } else {
            None
        }
    })
}

/// The chaining DP with one score per anchor (`chain_members`' DP; `chain_sorted` runs it with its scores in a ring,
/// `chain_dp_ring`), on its initialised arrays. Anchor j's predecessor is chosen among the
/// `lookback` anchors before it: one strictly before it on both reads, at most `max_gap` tokens away on each, whose
/// drift |dd - dt| is at most max_drift + drift_fraction x max(dt, dd); the chain's score grows by 1 - 0.1 x drift
/// (an anchor alone scores 1).
///
/// The anchors before `j` at `j`'s own target position (a block right before it, the anchors being sorted) can never
/// precede it (dt = 0), so the look-back starts below them; the window's lower end stays `j - lookback`. In a repeat a
/// target position matches many donor copies and that block fills the window: the same scores and back links, without
/// visiting it.
fn chain_dp<A: Anchor, B: BackLink>(
    anchors: &[A],
    p: &Params,
    score: &mut [f64],
    back: &mut [B],
) -> (Option<(usize, f64, i64, i64, i64, i64)>, usize) {
    let n = anchors.len();
    let mut block = 0usize;   // first anchor at the current target position
    for j in 0..n {
        let (tj, dj) = (anchors[j].t(), anchors[j].d());
        if j > 0 && anchors[j - 1].t() != tj {
            block = j;
        }
        let lo = j.saturating_sub(p.lookback);
        let mut i = block.max(lo);
        while i > lo {
            i -= 1;
            let (ti, di) = (anchors[i].t(), anchors[i].d());
            let (dt, dd) = (tj - ti, dj - di);
            if dt <= 0 || dd <= 0 {
                continue;
            }
            if dt > p.max_gap || dd > p.max_gap {
                if dt > p.max_gap {
                    break;
                }
                continue;
            }
            let drift = (dd - dt).abs();
            if (drift as f64) > p.max_drift as f64 + p.drift_fraction * (dt.max(dd) as f64) {
                continue;
            }
            let s = score[i] + 1.0 - 0.1 * drift as f64;
            if s > score[j] {
                score[j] = s;
                back[j] = B::link(j, i);
            }
        }
    }
    chain_end(anchors, score, back)
}

/// Bytes of chaining scratch per anchor of the group chained (`heavy`'s need): the back links (one byte when they fit
/// one; the scores are in a ring).
pub(super) fn chain_bytes_per_anchor(p: &Params) -> usize {
    if chain_back8(p) { 1 } else { 4 }
}

/// `chain_sorted`'s DP with its scores in a ring instead of one per anchor. An anchor's look-back reads only the
/// `lookback` scores before it (60 by default), so a ring of a power of two above that holds every score it reads; the
/// best end is tracked as the scores are made (the highest, the first index among equals, as `chain_end`). The same
/// arithmetic in the same order as the scalar `chain_dp`, so the same chain, score and back links
/// (`chain_ring_equals_the_full_dp`); the per-anchor scratch is then the back links alone, which matters for
/// target-donor pairs with hundreds of millions of anchors. `ring` is a power of two at least min(lookback + 1, n)
/// long.
fn chain_dp_ring<A: Anchor, B: BackLink>(
    anchors: &[A],
    p: &Params,
    ring: &mut [f64],
    back: &mut [B],
) -> (Option<(usize, f64, i64, i64, i64, i64)>, usize) {
    let n = anchors.len();
    let mask = ring.len() - 1;
    debug_assert!(ring.len().is_power_of_two() && ring.len() >= (p.lookback + 1).min(n));
    let (mut block, mut end, mut best) = (0usize, 0usize, f64::NEG_INFINITY);
    for j in 0..n {
        let (tj, dj) = (anchors[j].t(), anchors[j].d());
        if j > 0 && anchors[j - 1].t() != tj {
            block = j;
        }
        let lo = j.saturating_sub(p.lookback);
        let mut i = block.max(lo);
        let mut sj = 1.0f64;
        while i > lo {
            i -= 1;
            let (ti, di) = (anchors[i].t(), anchors[i].d());
            let (dt, dd) = (tj - ti, dj - di);
            if dt <= 0 || dd <= 0 {
                continue;
            }
            if dt > p.max_gap || dd > p.max_gap {
                if dt > p.max_gap {
                    break;
                }
                continue;
            }
            let drift = (dd - dt).abs();
            if (drift as f64) > p.max_drift as f64 + p.drift_fraction * (dt.max(dd) as f64) {
                continue;
            }
            let s = ring[i & mask] + 1.0 - 0.1 * drift as f64;
            if s > sj {
                sj = s;
                back[j] = B::link(j, i);
            }
        }
        ring[j & mask] = sj;
        if j == 0 || sj > best {
            best = sj;
            end = j;
        }
    }
    let (mut start, mut count) = (end, 1usize);
    while back[start] != B::NONE {
        start = back[start].target(start);
        count += 1;
    }
    (Some((count, best, anchors[start].t(), anchors[end].t(), anchors[start].d(), anchors[end].d())), end)
}

/// The end of `chain_dp`: the best chain's end, anchor count, score and ends.
fn chain_end<A: Anchor, B: BackLink>(anchors: &[A], score: &[f64], back: &[B]) -> (Option<(usize, f64, i64, i64, i64, i64)>, usize) {
    let n = anchors.len();
    // max over (score, -k): the highest score, the first index among equals
    let mut end = 0;
    for k in 1..n {
        if score[k] > score[end] {
            end = k;
        }
    }
    // the anchor count from the back links (one more anchor per link)
    let (mut start, mut count) = (end, 1usize);
    while back[start] != B::NONE {
        start = back[start].target(start);
        count += 1;
    }
    (Some((count, score[end], anchors[start].t(), anchors[end].t(), anchors[start].d(), anchors[end].d())), end)
}

/// The anchors of the chain `chain_sorted` / `chain_anchors` choose for `anchors` (sorted and distinct, or not), as
/// (target index, donor index) in chain order (for the q-symbol check's seed): the same DP, then its back links.
pub(super) fn chain_members<A: Anchor>(anchors: &[A], p: &Params) -> Vec<(i64, i64)> {
    let mut list: Vec<(i64, i64)> = anchors.iter().map(|a| (a.t(), a.d())).collect();
    if !list.windows(2).all(|w| w[0] < w[1]) {
        list.sort_unstable();
        list.dedup();
    }
    let n = list.len();
    if n == 0 {
        return Vec::new();
    }
    let (mut score, mut back) = (vec![1.0f64; n], vec![NO_BACK; n]);
    let (_, end) = chain_dp(&list, p, &mut score, &mut back);
    let mut out = vec![list[end]];
    let mut k = end;
    while back[k] != NO_BACK {
        k = back[k] as usize;
        out.push(list[k]);
    }
    out.reverse();
    out
}

/// s / 2 rounded to the nearest integer for an integer s, halves to the even neighbour (as Python's `round`).
pub fn round_half_even_div2(s: i64) -> i64 {
    let q = s.div_euclid(2);
    if s.rem_euclid(2) == 0 {
        q
    } else if q % 2 == 0 {
        q
    } else {
        q + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chain_scratch_reuse_gives_fresh_results() {
        // A long noisy anchor set, then a short one, then the long one again on the same thread; each must match the
        // chain computed on a new thread (fresh scratch).
        let p = Params::default();
        let long: Vec<(i64, i64)> = (0..400).map(|i| (i * 3, i * 3 + 7 + (i % 5))).chain((0..50).map(|i| (i * 17 + 1, i * 2))).collect();
        let short: Vec<(i64, i64)> = vec![(0, 10), (4, 14), (9, 19), (30, 2), (31, 41)];
        let fresh = |a: Vec<(i64, i64)>| std::thread::spawn(move || chain_anchors(a, &Params::default())).join().unwrap();
        let (l1, s1, l2) = (chain_anchors(long.clone(), &p), chain_anchors(short.clone(), &p), chain_anchors(long.clone(), &p));
        assert_eq!(l1, fresh(long.clone()));
        assert_eq!(s1, fresh(short));
        assert_eq!(l1, l2);
        assert!(l1.unwrap().0 > 100);
    }

    /// The chain members traced from `chain_sorted`'s kept back links equal `chain_members`' (a second DP), and another
    /// group chained in between is noticed (None).
    #[test]
    fn chain_members_after_sorted_equal_a_second_dp() {
        let (_, tokens) = repeat_reads(10, 5);
        let p = Params::default();
        let mut checked = 0;
        for r in 0..tokens.len() {
            for d in 0..tokens.len() {
                let mut a: Vec<(u32, u32)> = tokens[r].iter().enumerate().flat_map(|(i, x)| {
                    tokens[d].iter().enumerate().filter(move |(_, y)| y.value == x.value).map(move |(j, _)| (i as u32, j as u32))
                }).collect();
                a.sort_unstable();
                a.dedup();
                if a.len() < 3 || r == d {
                    continue;
                }
                let c = chain_sorted(&a, &p).unwrap();
                let m = chain_members_after_sorted(&a).expect("the last group");
                assert_eq!(m, chain_members(&a, &p));
                assert_eq!((m.len(), m[0].0, m[m.len() - 1].0), (c.0, c.2, c.3));
                let other: Vec<(u32, u32)> = a[..a.len() - 1].to_vec();
                chain_sorted(&other, &p);
                assert!(chain_members_after_sorted(&a).is_none());
                checked += 1;
            }
        }
        assert!(checked > 100, "{}", checked);
    }

    /// The chaining DP with its scores in a ring and one-byte back links when the look-back fits one (`chain_sorted`)
    /// gives the chain, score and members of the DP with every score kept and u32 links (`chain_dp`), at look-backs
    /// below and above 255.
    #[test]
    fn chain_ring_equals_the_full_dp() {
        let mut state = 29u64;
        let mut next = move || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as i64
        };
        let mut linked = 0;
        for (lookback, max_gap, max_drift, frac) in [(60usize, 120i64, 6i64, 0.25f64), (1, 120, 6, 0.25), (5, 30, 2, 0.1),
                                                      (60, 8, 0, 0.0), (200, 400, 20, 0.5), (63, 120, 6, 0.25), (64, 120, 6, 0.25),
                                                      (255, 400, 20, 0.5), (256, 400, 20, 0.5), (400, 2000, 40, 0.5)] {
            let p = Params { lookback, max_gap, max_drift, drift_fraction: frac, ..Params::default() };
            for _ in 0..300 {
                let mut a: Vec<(u32, u32)> = Vec::new();
                let positions = 1 + next() % 150;
                let mut t = 0i64;
                for _ in 0..positions {
                    t += 1 + next() % (if next() % 7 == 0 { 300 } else { 6 });
                    let copies = 1 + next() % (if next() % 5 == 0 { 40 } else { 3 });
                    for _ in 0..copies {
                        a.push((t as u32, (next() % 2000) as u32));
                    }
                }
                a.sort_unstable();
                a.dedup();
                let ring = (chain_sorted(&a, &p).map(|c| (c.0, c.1.to_bits(), c.2, c.3, c.4, c.5)), chain_members_after_sorted(&a));
                let (mut score, mut back) = (vec![1.0f64; a.len()], vec![NO_BACK; a.len()]);
                let (out, end) = chain_dp(&a, &p, &mut score, &mut back);
                let full = (out.map(|c| (c.0, c.1.to_bits(), c.2, c.3, c.4, c.5)), Some(chain_trace(&a, &back, end)));
                assert_eq!(ring, full, "lookback {} n {}", lookback, a.len());
                linked += full.1.map_or(0, |m| m.len());
            }
        }
        assert!(linked > 1000, "{}", linked);
    }



    #[test]
    fn halves_round_to_even_like_python() {
        let cases = [(5, 2), (7, 4), (-5, -2), (-7, -4), (6, 3), (-6, -3), (1, 0), (-1, 0), (3, 2), (-3, -2)];
        for (s, want) in cases {
            assert_eq!(round_half_even_div2(s), want, "round({}/2)", s);
        }
    }

    #[test]
    fn chain_is_forward_colinear() {
        let p = Params::default();
        let best = chain_anchors(vec![(0, 10), (5, 15), (9, 19), (12, 2), (20, 30)], &p).unwrap();
        assert_eq!(best.0, 4);
        assert_eq!((best.2, best.3, best.4, best.5), (0, 20, 10, 30));
        assert_eq!(chain_anchors(vec![(0, 30), (10, 20), (20, 10)], &p).unwrap().0, 1);
    }
}
