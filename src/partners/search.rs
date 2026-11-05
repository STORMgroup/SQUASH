//! One target's chains over the index (`target_chains`): counting and placing each donor's anchors (in batches
//! for repeat-rich targets), the donor order and the prune, and the search's per-thread scratch.

use super::*;

/// A kept chain of one donor: its anchor count and score, its first and last target tokens (t0, t1), and the row's
/// bounds on the target [ts, te) and on the donor [ds, de) (q-symbols; the implied overlap).
#[derive(Clone, Debug)]
pub struct Chain {
    pub donor: usize,
    pub anchors: usize,
    pub score: f64,
    pub t0: usize,
    pub t1: usize,
    pub ts: i64,
    pub te: i64,
    pub ds: i64,
    pub de: i64,
    /// The lazy check (`LAZY_CHECK`): the q-symbol alignment check this chain still needs, made by the coverage cap
    /// (`target_rows`) only if the chain reaches it with a position still open; None: checked already, or no check.
    pub check: Option<Box<Deferred>>,
}

/// What `signal::Check::verdict` needs for a chain whose check `DonorChainer::chain` deferred: the chain's line
/// (target q of the first anchor, donor q, target q of the last, donor q), the implied overlap [ts, te) on the target
/// and the chain's anchors as (target q, donor q).
#[derive(Clone, Debug)]
pub struct Deferred {
    pub line: (i64, i64, i64, i64),
    pub ts: i64,
    pub te: i64,
    pub q: Vec<(i64, i64)>,
}


/// The most tokens a read may have for its anchors to be (u16, u16).
const NARROW_TOKENS: usize = 1 << 16;

#[cfg(not(test))]
fn narrow_tokens() -> usize {
    NARROW_TOKENS
}

#[cfg(test)]
thread_local! {
    static TEST_NARROW_TOKENS: std::cell::Cell<usize> = std::cell::Cell::new(NARROW_TOKENS);
    static TEST_NARROW_TAKEN: std::cell::Cell<(usize, usize)> = std::cell::Cell::new((0, 0));
}

/// Tests: a lower limit, so that some test reads count as long.
#[cfg(test)]
fn narrow_tokens() -> usize {
    TEST_NARROW_TOKENS.with(|c| c.get())
}

/// How many occurrences ahead `target_chains`' counting and placing passes prefetch the per-donor counters, write
/// positions and anchor slots they will touch (random accesses into per-read arrays and the anchor buffers, most of the
/// search's time on long reads). Prefetching reads nothing the passes do not read, so the anchors, their order and the
/// chains do not depend on it.
const PREFETCH_DISTANCE: usize = 16;

/// The most anchors a donor's count holds (u32::MAX; tests set a smaller cap to reach the saturated path).
#[cfg(not(test))]
#[inline(always)]
fn count_cap() -> u32 {
    u32::MAX
}

#[cfg(test)]
thread_local! {
    static TEST_COUNT_CAP: std::cell::Cell<u32> = std::cell::Cell::new(u32::MAX);
}

#[cfg(test)]
fn count_cap() -> u32 {
    TEST_COUNT_CAP.with(|c| c.get())
}

/// A hint to load the cache line at `p` (x86-64 `prefetcht0`, baseline SSE; nothing elsewhere). Never faults.
#[inline(always)]
fn prefetch<T>(p: *const T) {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        std::arch::x86_64::_mm_prefetch(p as *const i8, std::arch::x86_64::_MM_HINT_T0);
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = p;
}


/// Per-thread buffers of `target_chains`. `count` is zero for every read between targets.
pub struct Scratch {
    pub(super) count: Vec<u32>,
    start: Vec<u32>,
    touched: Vec<u32>,
    hits: Vec<std::ops::Range<usize>>,
    anchors: Vec<(u32, u32)>,
    /// The (u16, u16) anchor buffer (`NARROW_TOKENS`).
    anchors16: Vec<(u16, u16)>,
    /// The target's (value, position) pairs and (first position, multiplicity, index range) per distinct value.
    values: Vec<(u32, u32)>,
    firsts: Vec<(u32, u32, std::ops::Range<usize>)>,
    /// The reads with more than `narrow_tokens()` tokens as bits, computed on first use for the token array (its
    /// address and length) and the limit.
    long: Option<((usize, usize, usize), Vec<u64>)>,
}

impl Scratch {
    pub fn new(reads: usize) -> Scratch {
        Scratch { count: vec![0; reads], start: vec![0; reads], touched: Vec::new(), hits: Vec::new(), anchors: Vec::new(),
                  anchors16: Vec::new(), values: Vec::new(), firsts: Vec::new(), long: None }
    }
    /// After a target, frees the buffers grown past `SCRATCH_KEEP` bytes by an unusually repetitive target (each thread
    /// keeps its buffers for the next target otherwise).
    pub fn trim(&mut self) {
        fn keep<T>(v: &mut Vec<T>) {
            if v.capacity() * std::mem::size_of::<T>() > SCRATCH_KEEP {
                *v = Vec::new();
            }
        }
        keep(&mut self.anchors);
        keep(&mut self.anchors16);
        keep(&mut self.hits);
        CHAIN_SCRATCH.with(|cell| {
            let mut s = cell.borrow_mut();
            keep(&mut s.score);
            keep(&mut s.back);
            keep(&mut s.back8);
            s.last = (0, 0, 0);
            if memtrack::on() {
                memtrack::set(memtrack::CHAIN, s.bytes());
            }
        });
        if memtrack::on() {
            memtrack::set(memtrack::ANCHORS, self.anchors.capacity() * 8 + self.anchors16.capacity() * 4);
            memtrack::set(memtrack::ORDER, self.touched.capacity() * 4 + self.hits.capacity() * 16);
        }
    }
}

/// The most bytes a per-thread buffer keeps between targets (`Scratch::trim`; each thread could otherwise keep several
/// large buffers).
pub(super) const SCRATCH_KEEP: usize = 4 << 20;

/// `target_chains`: the most anchors placed at once (16 M, 128 MB) before a target's donors go in batches.
pub(super) const BATCH_ANCHORS: usize = 16 << 20;

/// The batch size, `BATCH_ANCHORS` (each batch rescans the target's hits; the same chains at any size).
#[cfg(not(test))]
fn batch_anchors() -> usize {
    BATCH_ANCHORS
}

#[cfg(test)]
thread_local! {
    pub(super) static TEST_BATCH_ANCHORS: std::cell::Cell<usize> = std::cell::Cell::new(BATCH_ANCHORS);
}

/// Tests: a smaller batch, so that test targets go in batches.
#[cfg(test)]
fn batch_anchors() -> usize {
    TEST_BATCH_ANCHORS.with(|c| c.get())
}


/// A kept chain (score, t0, t1) in `target_chains`' prune heap, ordered by score (finite).
pub(super) struct ByScore(pub(super) f64, pub(super) usize, pub(super) usize);

impl PartialEq for ByScore {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}

impl Eq for ByScore {}

impl PartialOrd for ByScore {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ByScore {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0).then(self.1.cmp(&other.1)).then(self.2.cmp(&other.2))
    }
}

/// The first open position at or after `x` (`target_chains`' prune; `next[x] == x` = open, the last entry a sentinel
/// that is always open), compressing the links it follows.
pub(super) fn open_from(next: &mut [u32], x: usize) -> usize {
    let mut root = x;
    while next[root] as usize != root {
        root = next[root] as usize;
    }
    let mut y = x;
    while next[y] as usize != root {
        let z = next[y] as usize;
        next[y] = root as u32;
        y = z;
    }
    root
}

#[allow(clippy::too_many_arguments)]
/// The (target, donor) pairs dropped because the donor's anchor count saturated u32 (see `target_chains`); reported and
/// reset by `search_targets`.
pub(super) static SATURATED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// The lazy check: the q-symbol alignment check deferred from chaining to the coverage cap (`DonorChainer::chain`,
/// `target_rows_with`), so it is made only for chains the cap can still keep, in the cap's order (score). Within a
/// search the rows are those of checking every chain at chaining (given the same verdicts); a pair's verdict is cached
/// from the first search whose cap reaches it, computed from that search's chain, and reused by later searches. Counts:
/// chains deferred, checked at the cap, passed over unchecked.
pub(super) static LAZY_CHECK: [std::sync::atomic::AtomicU64; 3] = [const { std::sync::atomic::AtomicU64::new(0) }; 3];

/// The lazy check's counts, reported and reset (both search paths).
pub(super) fn report_lazy_check() {
    use std::sync::atomic::Ordering::Relaxed;
    let c: Vec<u64> = LAZY_CHECK.iter().map(|x| x.swap(0, Relaxed)).collect();
    eprintln!("[partners] lazy check: {} chains deferred to the cap, {} checked there, {} passed over unchecked",
              c[0], c[1], c[2]);
}

/// One target's chains (the search's loop body). `guard`: the uncorrected tokens, their `kept_values` and the rule. The
/// anchors are grouped by donor in two passes over the index (count, then place); each donor's anchors come out in
/// (target index, donor index) order, as `chain_anchors` would sort them.
///
/// `prune` (`partner_rows` passes `Some(p.coverage)`): donors are chained in decreasing anchor count (a bound on the
/// chain's score), and a donor is skipped when every target position its anchors span is already covered by `coverage`
/// kept chains scoring more than its count. Its chain could only score that much and would then fail `target_rows`'s
/// coverage cap (every chain ranked before it either adds depth there or found the depth full), so the rows are the
/// same; only chains that could never become rows are left out. `None`: every donor is chained. A donor with fewer
/// than `min_anchors` anchors is not placed or chained.
pub fn target_chains(
    r: usize,
    lengths: &[i64],
    tokens: &[Vec<Token>],
    index: &Index,
    p: &Params,
    guard: Option<(&[GuardRead], &HashSet<u32>, &Guard)>,
    scratch: &mut Scratch,
    prune: Option<usize>,
) -> Vec<Chain> {
    let Scratch { count, start, touched, hits, anchors, anchors16, values, firsts, long } = scratch;
    let prof = profile::on();
    let clock = std::time::Instant::now();
    touched.clear();
    hits.clear();
    // The counting pass by distinct value: a value at m target positions scans its occurrences once instead of m times
    // (tandem arrays), adding m to each donor's count. Donors are first touched at a value's first position, in target
    // order. In the placing pass only donors with at least `min_anchors` anchors are placed (the others cannot give a
    // chain).
    // A donor's anchor count is at most the target's hits (each occurrence adds to one donor), so below u32::MAX hits
    // no count can overflow. A target with more (two low-complexity reads that each hold one token value tens of
    // thousands of times give the pair billions of anchors) counts with saturation, and a saturated donor is dropped
    // below.
    let cap = count_cap();
    // Each distinct value once, in the order of its first position: its occurrences add its multiplicity m to the
    // donor's anchors.
    let pd = PREFETCH_DISTANCE;
    let n_t = tokens[r].len();
    values.clear();
    values.extend(tokens[r].iter().enumerate().map(|(i, t)| (t.value, i as u32)));
    values.sort_unstable();
    hits.resize(n_t, 0..0);
    firsts.clear();
    let mut k = 0;
    while k < n_t {
        let mut e = k + 1;
        while e < n_t && values[e].0 == values[k].0 {
            e += 1;
        }
        let range = index.range(values[k].0);
        for &(_, i) in &values[k..e] {
            hits[i as usize] = range.clone();
        }
        if !range.is_empty() {
            firsts.push((values[k].1, (e - k) as u32, range));
        }
        k = e;
    }
    firsts.sort_unstable_by_key(|f| f.0);
    let sat = firsts.iter().map(|f| f.1 as u64 * f.2.len() as u64).sum::<u64>() >= cap as u64;
    for (_, m, range) in firsts.iter() {
        for k in range.clone() {
            if k + pd < range.end {
                prefetch(count.as_ptr().wrapping_add(index.occ.read(k + pd) as usize));
            }
            let d = index.occ.read(k);
            let du = d as usize;
            if du != r {
                if count[du] == 0 {
                    touched.push(d);
                }
                count[du] = if sat { count[du].saturating_add(*m).min(cap) } else { count[du] + m };
            }
        }
    }
    if sat {
        // the donors whose count saturated hold no anchors (so placed nowhere and not chained)
        let mut dropped = 0u64;
        for &d in touched.iter() {
            let d = d as usize;
            if count[d] >= cap {
                count[d] = 0;
                dropped += 1;
            }
        }
        if dropped > 0 {
            SATURATED.fetch_add(dropped, std::sync::atomic::Ordering::Relaxed);
        }
    }
    const UNPLACED: u32 = u32::MAX;
    // A target with more than `batch_anchors()` anchors to place (a repeat-rich read: hundreds of millions, GB at once)
    // places and chains its donors in batches of consecutive donors of the chaining order, each batch's anchors placed
    // by another pass over the target's hits. The same anchors per donor, chaining order and pruning state, so the same
    // chains.
    let placeable: usize = touched.iter().map(|&d| count[d as usize] as usize).filter(|&c| c >= p.min_anchors).sum();
    let batch_limit = batch_anchors();
    let batched = placeable > batch_limit;
    // On a heavy thread (`heavy`), the bytes of the largest placement (the anchors placed at once, 8 bytes each at
    // most, and the chaining scratch of the largest donor, `chain_bytes_per_anchor`) taken before placing anything.
    if heavy::active() {
        let largest = touched.iter().map(|&d| count[d as usize] as usize).filter(|&c| c >= p.min_anchors).max().unwrap_or(0);
        let at_once = if batched { batch_limit.max(largest) } else { placeable };
        heavy::take(at_once.min(placeable) * 8 + largest * chain_bytes_per_anchor(p));
    }
    // Each donor's anchors as (u16, u16) in `anchors16` (4 bytes) when the target and the donor have at most
    // `narrow_tokens()` tokens, else as (u32, u32) in `anchors` (the few longer reads; every donor of a longer target):
    // the same anchors per donor, order, batches and chains either way.
    let limit = narrow_tokens();
    let all_wide = tokens[r].len() > limit;
    if !all_wide {
        let key = (tokens.as_ptr() as usize, tokens.len(), limit);
        if long.as_ref().map_or(true, |l| l.0 != key) {
            let mut bits = vec![0u64; tokens.len().div_ceil(64)];
            for (d, t) in tokens.iter().enumerate() {
                if t.len() > limit {
                    bits[d / 64] |= 1 << (d % 64);
                }
            }
            *long = Some((key, bits));
        }
    }
    let long_bits: &[u64] = if all_wide { &[] } else { &long.as_ref().expect("the long reads").1 };
    let wide_d = |d: usize| all_wide || (long_bits[d >> 6] >> (d & 63)) & 1 == 1;
    // Places the anchors of the donors whose write position `count` is set (their position in their buffer; the donors
    // not placed now have count UNPLACED, so that a placed anchor reads one per-read array, not two).
    let pd = PREFETCH_DISTANCE;
    let place = |count: &mut [u32], nb: &mut [(u16, u16)], wb: &mut [(u32, u32)]| {
        for (i, range) in hits.iter().enumerate() {
            for k in range.clone() {
                // `pd` occurrences ahead the donor's write position; `pd / 2` ahead, by then cached, the anchor slot
                // it will write (only for a placed donor, whose position is within its buffer)
                if k + pd < range.end {
                    prefetch(count.as_ptr().wrapping_add(index.occ.read(k + pd) as usize));
                }
                if k + pd / 2 < range.end {
                    let d2 = index.occ.read(k + pd / 2) as usize;
                    if count[d2] != UNPLACED {
                        let at = count[d2] as usize;
                        if wide_d(d2) {
                            prefetch(wb.as_ptr().wrapping_add(at));
                        } else {
                            prefetch(nb.as_ptr().wrapping_add(at));
                        }
                    }
                }
                let d = index.occ.read(k) as usize;
                if d != r && count[d] != UNPLACED {
                    let j = index.occ.pos(k);
                    if wide_d(d) {
                        wb[count[d] as usize] = (i as u32, j);
                    } else {
                        nb[count[d] as usize] = (i as u16, j as u16);
                    }
                    count[d] += 1;
                }
            }
        }
    };
    let total: u32;
    if !batched {
        let (mut total_n, mut total_w) = (0u32, 0u32);
        for &d in touched.iter() {
            let d = d as usize;
            if (count[d] as usize) < p.min_anchors {
                start[d] = UNPLACED; // too few anchors for a chain
                count[d] = UNPLACED;
            } else {
                let at = if wide_d(d) { &mut total_w } else { &mut total_n };
                start[d] = *at;
                *at += count[d];
                count[d] = start[d]; // now the write position
            }
        }
        total = total_n + total_w;
        anchors16.clear();
        anchors16.resize(total_n as usize, (0, 0));
        anchors.clear();
        anchors.resize(total_w as usize, (0, 0));
        memtrack::set(memtrack::ANCHORS, anchors16.capacity() * 4 + anchors.capacity() * 8);
        place(count, &mut anchors16[..], &mut anchors[..]);
    } else {
        total = placeable as u32;
    }
    if prof {
        profile::since(2, clock);
        profile::add(7, 1);
        profile::add(8, total as u64);
        profile::add(9, touched.len() as u64);
    }
    let oclock = std::time::Instant::now();
    // Donor order (donor, first anchor, end, anchors): as gathered, or by decreasing anchor count when pruning.
    let mut order: Vec<(usize, usize, usize, usize)> = Vec::with_capacity(touched.len());
    for &d in touched.iter() {
        let d = d as usize;
        if batched {
            // (0, count): placed batch by batch below
            let c = count[d] as usize;
            count[d] = 0;
            if c >= p.min_anchors {
                order.push((d, 0, c, c));
            }
            continue;
        }
        let (first, end) = (start[d] as usize, count[d] as usize);
        count[d] = 0;
        if start[d] == UNPLACED {
            continue;
        }
        if end - first >= p.min_anchors {
            order.push((d, first, end, end - first));
        }
    }
    if prune.is_some() {
        order.sort_by(|a, b| b.3.cmp(&a.3));
    }
    #[cfg(test)]
    TEST_NARROW_TAKEN.with(|c| {
        let (n, w) = c.get();
        c.set((n + order.iter().any(|o| !wide_d(o.0)) as usize, w + order.iter().any(|o| wide_d(o.0)) as usize));
    });
    memtrack::set(memtrack::ORDER, order.capacity() * 32 + touched.capacity() * 4 + hits.capacity() * 16);
    let mut dc = DonorChainer::new(r, lengths, tokens, p, guard, prune, prof);
    if prof {
        profile::since(16, oclock);
    }
    let mut batches: Vec<std::ops::Range<usize>> = Vec::new();
    if batched {
        let (mut lo, mut sum) = (0usize, 0usize);
        for k in 0..order.len() {
            let c = order[k].2 - order[k].1;
            if k > lo && sum + c > batch_limit {
                batches.push(lo..k);
                lo = k;
                sum = 0;
            }
            sum += c;
        }
        batches.push(lo..order.len());
    } else {
        batches.push(0..order.len());
    }
    if prof && batched {
        profile::add(19, 1);
        profile::add(20, batches.len() as u64);
    }
    for batch in batches {
        if batched {
            let bclock = if prof { Some(std::time::Instant::now()) } else { None };
            for &d in touched.iter() {
                start[d as usize] = UNPLACED;
                count[d as usize] = UNPLACED;
            }
            let (mut off_n, mut off_w) = (0usize, 0usize);
            for k in batch.clone() {
                let (d, first, end, b) = order[k];
                let off = if wide_d(d) { &mut off_w } else { &mut off_n };
                start[d] = *off as u32;
                count[d] = *off as u32;
                order[k] = (d, *off, *off + end - first, b);
                *off += end - first;
            }
            anchors16.clear();
            anchors16.resize(off_n, (0, 0));
            anchors.clear();
            anchors.resize(off_w, (0, 0));
            memtrack::set(memtrack::ANCHORS, anchors16.capacity() * 4 + anchors.capacity() * 8);
            place(count, &mut anchors16[..], &mut anchors[..]);
            if let Some(t) = bclock {
                profile::since(21, t);
            }
        }
        for &(d, first, end, b) in &order[batch] {
            if wide_d(d) {
                let group = &anchors[first..end];
                if !dc.pruned(b, group[0].t() as usize, group[group.len() - 1].t() as usize, group.len()) {
                    dc.chain(d, group);
                }
            } else {
                let group = &anchors16[first..end];
                if !dc.pruned(b, group[0].t() as usize, group[group.len() - 1].t() as usize, group.len()) {
                    dc.chain(d, group);
                }
            }
        }
    }
    if batched {
        for &d in touched.iter() {
            count[d as usize] = 0; // every touched donor's (the last batch's donors and the UNPLACED others)
        }
    }
    dc.chains
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The coverage-cap prune leaves every target's rows unchanged (guard or not, any coverage), and does prune.
    #[test]
    fn pruned_chaining_keeps_the_rows() {
        let mut pruned = 0;
        for (seed, noise) in [(1, 10), (2, 30), (6, 5), (7, 20)] {
            let (lengths, tokens) = noisy_reads(noise, seed);
            let (_, raw) = noisy_reads(60, seed);
            let ids: Vec<String> = (0..tokens.len()).map(|r| format!("read{:02}", r)).collect();
            let g = Guard { min_anchors: 2, band: 30, band_fraction: 0.02 };
            let guard_kept = kept_values(&raw, 0, 0.0);
            let guard_view = guard_tokens(&raw, &guard_kept);   // the guard's view
            for coverage in [1, 2, 3, 20] {
                let p = Params { max_reads: 0, min_anchors: 3, min_overlap: 500, min_cover: 0.1, coverage, ..Params::default() };
                let index = build_index(&tokens, 0, 0.0);
                let mut scratch = Scratch::new(tokens.len());
                for r in 0..tokens.len() {
                    for guard in [None, Some((&guard_view[..], &guard_kept, &g))] {
                        let all = target_chains(r, &lengths, &tokens, &index, &p, guard, &mut scratch, None);
                        let some = target_chains(r, &lengths, &tokens, &index, &p, guard, &mut scratch, Some(coverage));
                        pruned += all.len() - some.len();
                        assert_eq!(target_rows(r, &ids, &tokens, some, coverage), target_rows(r, &ids, &tokens, all, coverage),
                                   "seed {} coverage {} target {}", seed, coverage, r);
                    }
                }
            }
        }
        assert!(pruned > 0);
    }

    /// Batched anchor placement (a target with more than `batch_anchors()` anchors) gives the chains of the one-buffer
    /// placement, with and without pruning and guard; the reused scratch is clean after it.
    #[test]
    fn batched_placement_keeps_the_chains() {
        let key = |c: &Chain| (c.donor, c.anchors, c.score.to_bits(), c.t0, c.t1, c.ts, c.te, c.ds, c.de);
        let mut batched_targets = 0;
        for (seed, noise) in [(1, 10), (2, 30), (6, 5)] {
            let (lengths, tokens) = noisy_reads(noise, seed);
            let (_, raw) = noisy_reads(60, seed);
            let g = Guard { min_anchors: 2, band: 30, band_fraction: 0.02 };
            let guard_kept = kept_values(&raw, 0, 0.0);
            let guard_view = guard_tokens(&raw, &guard_kept);   // the guard's view
            let p = Params { max_reads: 0, min_anchors: 3, min_overlap: 500, min_cover: 0.1, coverage: 2, ..Params::default() };
            let index = build_index(&tokens, 0, 0.0);
            let mut scratch = Scratch::new(tokens.len());
            for r in 0..tokens.len() {
                for guard in [None, Some((&guard_view[..], &guard_kept, &g))] {
                    for prune in [None, Some(2)] {
                        let run = |limit: usize, scratch: &mut Scratch| {
                            TEST_BATCH_ANCHORS.with(|c| c.set(limit));
                            let out: Vec<_> = target_chains(r, &lengths, &tokens, &index, &p, guard, scratch, prune)
                                .iter().map(key).collect();
                            TEST_BATCH_ANCHORS.with(|c| c.set(BATCH_ANCHORS));
                            out
                        };
                        let whole = run(BATCH_ANCHORS, &mut scratch);
                        for limit in [1, 7, 40] {
                            assert_eq!(run(limit, &mut scratch), whole, "seed {} target {} limit {}", seed, r, limit);
                        }
                        batched_targets += 1;
                    }
                }
            }
            assert!(scratch.count.iter().all(|&c| c == 0));
        }
        assert!(batched_targets > 0);
    }

    /// With a count cap (u32::MAX; a small test cap here), a donor whose anchor count reaches it is dropped and every
    /// other donor keeps its chain; a cap above every donor's count (the saturated pass taken, nothing saturated) gives
    /// the same chains. Small batches too; the scratch is clean after.
    #[test]
    fn saturated_counts_drop_only_their_donors() {
        let key = |c: &Chain| (c.donor, c.anchors, c.score.to_bits(), c.t0, c.t1, c.ts, c.te, c.ds, c.de);
        let (mut dropped, mut kept) = (0, 0);
        for (seed, noise, repeats) in [(1, 10, true), (3, 20, true), (2, 30, false)] {
            let (lengths, tokens) = if repeats { repeat_reads(noise, seed) } else { noisy_reads(noise, seed) };
            let p = Params { max_reads: 0, min_anchors: 3, min_overlap: 500, min_cover: 0.1, coverage: 2, ..Params::default() };
            let index = build_index(&tokens, 0, 0.0);
            let mut scratch = Scratch::new(tokens.len());
            for r in 0..tokens.len() {
                let counts: Vec<u64> = (0..tokens.len())
                    .map(|d| if d == r { 0 } else {
                        tokens[r].iter().map(|t| tokens[d].iter().filter(|u| u.value == t.value).count() as u64).sum()
                    })
                    .collect();
                for batch in [BATCH_ANCHORS, 40] {
                    let mut run = |cap: u64| {
                        TEST_BATCH_ANCHORS.with(|c| c.set(batch));
                        TEST_COUNT_CAP.with(|c| c.set(cap as u32));
                        let out: Vec<_> = target_chains(r, &lengths, &tokens, &index, &p, None, &mut scratch, None)
                            .iter().map(key).collect();
                        TEST_BATCH_ANCHORS.with(|c| c.set(BATCH_ANCHORS));
                        TEST_COUNT_CAP.with(|c| c.set(u32::MAX));
                        assert!(scratch.count.iter().all(|&c| c == 0), "seed {} target {}", seed, r);
                        out
                    };
                    let all = run(u32::MAX as u64);
                    if all.is_empty() {
                        continue;
                    }
                    let top = counts.iter().copied().max().unwrap();
                    assert_eq!(run(top + 1), all, "seed {} target {} cap above every count", seed, r);
                    let mut chained: Vec<u64> = all.iter().map(|k| counts[k.0]).collect();
                    chained.sort_unstable();
                    let cap = chained[chained.len() / 2];
                    let expected: Vec<_> = all.iter().filter(|k| counts[k.0] < cap).cloned().collect();
                    assert_eq!(run(cap), expected, "seed {} target {} cap {}", seed, r, cap);
                    dropped += all.len() - expected.len();
                    kept += expected.len();
                }
            }
        }
        assert!(dropped > 0 && kept > 0, "dropped {} kept {}", dropped, kept);
    }

    /// (u16, u16) anchors give the chains of (u32, u32) anchors (every read counted long), with a long-read limit that
    /// makes some targets and donors take the (u32, u32) buffer, alone or beside (u16, u16) donors (pruning, guard,
    /// batched placement); the scratch is clean after.
    #[test]
    fn narrow_anchors_keep_the_chains() {
        let key = |c: &Chain| (c.donor, c.anchors, c.score.to_bits(), c.t0, c.t1, c.ts, c.te, c.ds, c.de);
        TEST_NARROW_TAKEN.with(|c| c.set((0, 0)));
        for (seed, noise, repeats) in [(1, 10, true), (2, 30, false), (6, 5, false)] {
            let (lengths, tokens) = if repeats { repeat_reads(noise, seed) } else { noisy_reads(noise, seed) };
            let (_, raw) = noisy_reads(60, seed);
            let g = Guard { min_anchors: 2, band: 30, band_fraction: 0.02 };
            let guard_kept = kept_values(&raw, 0, 0.0);
            let guard_view = guard_tokens(&raw, &guard_kept);
            let p = Params { max_reads: 0, min_anchors: 3, min_overlap: 500, min_cover: 0.1, coverage: 2, ..Params::default() };
            let index = build_index(&tokens, 0, 0.0);
            let mut lens: Vec<usize> = tokens.iter().map(Vec::len).collect();
            lens.sort_unstable();
            let mut scratch = Scratch::new(tokens.len());
            for r in 0..tokens.len() {
                let guards = if repeats { vec![None] } else { vec![None, Some((&guard_view[..], &guard_kept, &g))] };
                for guard in guards {
                    for prune in [None, Some(2)] {
                        for batch in [BATCH_ANCHORS, 40] {
                            let run = |limit: usize, scratch: &mut Scratch| {
                                TEST_NARROW_TOKENS.with(|c| c.set(limit));
                                TEST_BATCH_ANCHORS.with(|c| c.set(batch));
                                let out: Vec<_> = target_chains(r, &lengths, &tokens, &index, &p, guard, scratch, prune)
                                    .iter().map(key).collect();
                                TEST_NARROW_TOKENS.with(|c| c.set(NARROW_TOKENS));
                                TEST_BATCH_ANCHORS.with(|c| c.set(BATCH_ANCHORS));
                                out
                            };
                            let wide = run(0, &mut scratch);
                            for limit in [NARROW_TOKENS, lens[lens.len() / 2], lens[lens.len() * 3 / 4]] {
                                assert_eq!(run(limit, &mut scratch), wide, "seed {} target {} limit {}", seed, r, limit);
                            }
                        }
                    }
                }
            }
            assert!(scratch.count.iter().all(|&c| c == 0));
        }
        let (narrow, wide) = TEST_NARROW_TAKEN.with(|c| c.get());
        assert!(narrow > 0 && wide > 0, "narrow {} wide {}", narrow, wide);
    }
}
