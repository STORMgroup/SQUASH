//! The index of the seed values: the kept values (repeat rule, reads cap) and each one's occurrences.

use super::*;

/// The index's occurrences (read ordinal, token index) in 6 bytes: the ordinal in the low 24 bits of `read`, the token
/// index's bits 16-23 in its high 8 bits and bits 0-15 in `low`; the search's counting pass reads only `read`. Wide
/// (some ordinal or token index needs more than 24 bits): `read` the whole ordinal, `wide` the whole token index.
/// Either layout gives the same (ordinal, token index) at every position.
#[derive(Default)]
pub struct Occurrences {
    read: Vec<u32>,
    low: Vec<u16>,
    wide: Vec<u32>,
    is_wide: bool,
    mask: u32,
}

impl Occurrences {
    pub(super) fn new(tokens: &[Vec<Token>], capacity: usize) -> Occurrences {
        let is_wide = tokens.len() > 1 << 24 || tokens.iter().any(|t| t.len() > 1 << 24);
        Occurrences::with_layout(is_wide, capacity)
    }
    fn with_layout(is_wide: bool, capacity: usize) -> Occurrences {
        Occurrences {
            read: Vec::with_capacity(capacity),
            low: if is_wide { Vec::new() } else { Vec::with_capacity(capacity) },
            wide: if is_wide { Vec::with_capacity(capacity) } else { Vec::new() },
            is_wide,
            mask: if is_wide { u32::MAX } else { (1 << 24) - 1 },
        }
    }
    fn push(&mut self, d: u32, j: u32) {
        if self.is_wide {
            self.read.push(d);
            self.wide.push(j);
        } else {
            debug_assert!(d < 1 << 24 && j < 1 << 24);
            self.read.push(d | ((j >> 16) << 24));
            self.low.push(j as u16);
        }
    }
    /// Occurrence k's read ordinal.
    #[inline]
    pub(super) fn read(&self, k: usize) -> u32 {
        self.read[k] & self.mask
    }
    /// Occurrence k's token index.
    #[inline]
    pub(super) fn pos(&self, k: usize) -> u32 {
        if self.is_wide { self.wide[k] } else { ((self.read[k] >> 24) << 16) | self.low[k] as u32 }
    }
    #[cfg(test)]
    #[inline]
    pub(super) fn get(&self, k: usize) -> (u32, u32) {
        (self.read(k), self.pos(k))
    }
    #[cfg(test)]
    fn len(&self) -> usize {
        self.read.len()
    }
    fn copy_within(&mut self, src: std::ops::Range<usize>, dest: usize) {
        self.read.copy_within(src.clone(), dest);
        if self.is_wide { self.wide.copy_within(src, dest) } else { self.low.copy_within(src, dest) }
    }
    fn copy_one(&mut self, from: usize, to: usize) {
        self.read[to] = self.read[from];
        if self.is_wide { self.wide[to] = self.wide[from] } else { self.low[to] = self.low[from] }
    }
    fn truncate(&mut self, n: usize) {
        self.read.truncate(n);
        self.low.truncate(n);
        self.wide.truncate(n);
    }
    fn shrink_to_fit(&mut self) {
        self.read.shrink_to_fit();
        self.low.shrink_to_fit();
        self.wide.shrink_to_fit();
    }
    pub(super) fn bytes(&self) -> usize {
        self.read.capacity() * 4 + self.low.capacity() * 2 + self.wide.capacity() * 4
    }
}

/// Every occurrence `(value, read ordinal, token index)` grouped by value: the occurrences in 2^`INDEX_PART_BITS` parts
/// by a hash of the value (one part's 12-byte entries collected and sorted at a time, in parallel on the current rayon
/// pool), each part split into runs `(value, first, length, distinct reads)`, and its (read ordinal, token index) pairs
/// appended part after part (`flat`; not kept when false: `kept_values` needs only the runs). A value's occurrences are
/// in read and token order.
fn value_runs_parts(tokens: &[Vec<Token>], flat: bool) -> (Occurrences, Vec<(u32, usize, usize, usize)>) {
    let part = |v: u32| v.wrapping_mul(0x9E37_79B1) >> (32 - INDEX_PART_BITS);
    let total: usize = if flat { tokens.iter().map(Vec::len).sum() } else { 0 };
    let mut pairs = if flat { Occurrences::new(tokens, total) } else { Occurrences::default() };
    let mut runs = Vec::new();
    let mut base = 0usize;
    for p in 0..(1u32 << INDEX_PART_BITS) {
        let mut entries: Vec<(u32, u32, u32)> = tokens
            .par_iter()
            .enumerate()
            .flat_map_iter(|(r, read)| {
                read.iter().enumerate().filter(move |(_, t)| part(t.value) == p).map(move |(i, t)| (t.value, r as u32, i as u32))
            })
            .collect();
        entries.par_sort_unstable();
        let mut first = 0;
        while first < entries.len() {
            let value = entries[first].0;
            let (mut end, mut reads) = (first + 1, 1);
            while end < entries.len() && entries[end].0 == value {
                reads += (entries[end].1 != entries[end - 1].1) as usize;
                end += 1;
            }
            runs.push((value, base + first, end - first, reads));
            first = end;
        }
        base += entries.len();
        if flat {
            for &(_, r, i) in &entries {
                pairs.push(r, i);
            }
        }
    }
    (pairs, runs)
}

/// `value_runs_parts`: 2^4 = 16 parts (one part's entries at a time: 1/16 of the occurrences).
const INDEX_PART_BITS: u32 = 4;

/// Which runs `build_index` keeps: values seen in 2..=max_reads distinct reads
/// (0: no cap), after dropping the int(repeat_fraction x distinct values) values
/// seen in the most reads (ties: smaller value first).
fn kept_runs(runs: &[(u32, usize, usize, usize)], max_reads: usize, repeat_fraction: f64) -> Vec<bool> {
    let mut kept: Vec<bool> =
        runs.iter().map(|&(_, _, _, reads)| reads >= 2 && (max_reads == 0 || reads <= max_reads)).collect();
    if repeat_fraction > 0.0 {
        let mut ranked: Vec<usize> = (0..runs.len()).collect();
        ranked.sort_unstable_by(|&a, &b| runs[b].3.cmp(&runs[a].3).then(runs[a].0.cmp(&runs[b].0)));
        let n_drop = (repeat_fraction * runs.len() as f64) as usize;
        for &k in &ranked[..n_drop.min(ranked.len())] {
            kept[k] = false;
        }
    }
    kept
}

/// The token values `build_index` keeps (see `kept_runs`).
pub fn kept_values(tokens: &[Vec<Token>], max_reads: usize, repeat_fraction: f64) -> HashSet<u32> {
    let runs = value_runs_parts(tokens, false).1;
    let kept = kept_runs(&runs, max_reads, repeat_fraction);
    runs.iter().zip(kept).filter(|(_, k)| *k).map(|(run, _)| run.0).collect()
}

/// `build_index`: the occurrences `(read ordinal, token index)` of each kept
/// value, in read and token order, as one flat array.
pub struct Index {
    slots: HashMap<u32, (usize, usize)>,
    /// The occurrences (`Occurrences`, 6 bytes each).
    pub(super) occ: Occurrences,
}

impl Index {
    pub(super) fn range(&self, value: u32) -> std::ops::Range<usize> {
        self.slots.get(&value).map_or(0..0, |&(first, length)| first..first + length)
    }
    #[cfg(test)]
    pub fn get(&self, value: u32) -> Vec<(u32, u32)> {
        self.range(value).map(|k| self.occ.get(k)).collect()
    }
    #[cfg(test)]
    pub fn values(&self) -> impl Iterator<Item = u32> + '_ {
        self.slots.keys().copied()
    }
}

#[cfg(test)]
/// `build_index` (on the current rayon pool).
pub fn build_index(tokens: &[Vec<Token>], max_reads: usize, repeat_fraction: f64) -> Index {
    build_index_excluding(tokens, max_reads, repeat_fraction, &[])
}

/// `build_index` without the occurrences of the `excluded` reads (by ordinal; empty: none), which then are never
/// donors. The kept values are still counted over every read, as metaMDBG counts its repetitive minimizers before
/// correction.
pub fn build_index_excluding(tokens: &[Vec<Token>], max_reads: usize, repeat_fraction: f64, excluded: &[bool]) -> Index {
    let clock = std::time::Instant::now();
    let (mut flat, runs) = value_runs_parts(tokens, true);
    if profile::on() {
        profile::since(0, clock);
    }
    let clock = std::time::Instant::now();
    let kept = kept_runs(&runs, max_reads, repeat_fraction);
    if profile::on() {
        // Kept values by the number of reads holding them: values, occurrences, and sum of occurrences^2 (anchor
        // pairs).
        let edges = [2usize, 31, 61, 101, 301, 1001, 3001, usize::MAX];
        let mut h = vec![(0u64, 0u64, 0u128); edges.len() - 1];
        for (&(_, _, len, reads), &k) in runs.iter().zip(&kept) {
            if k {
                let b = edges.windows(2).position(|w| reads >= w[0] && reads < w[1]).unwrap();
                h[b].0 += 1;
                h[b].1 += len as u64;
                h[b].2 += (len as u128) * (len as u128);
            }
        }
        for (b, (v, o, sq)) in h.iter().enumerate() {
            eprintln!("[partners-profile] values in {}..{} reads: {} values, {} occurrences, {} occurrence pairs",
                      edges[b], edges[b + 1].saturating_sub(1), v, o, sq);
        }
    }
    let mut slots = HashMap::new();
    // the kept runs moved left in place (in array order a run never moves right), the excluded reads left out
    let mut write = 0usize;
    for (&(value, first, length, _), &keep) in runs.iter().zip(&kept) {
        if keep {
            let at = write;
            if excluded.is_empty() {
                flat.copy_within(first..first + length, write);
                write += length;
            } else {
                for j in first..first + length {
                    if !excluded[flat.read(j) as usize] {
                        flat.copy_one(j, write);
                        write += 1;
                    }
                }
            }
            slots.insert(value, (at, write - at));
        }
    }
    flat.truncate(write);
    flat.shrink_to_fit();
    let occurrences = flat;
    if profile::on() {
        profile::since(1, clock);
    }
    mem::report("partner index built", &[("value_runs", runs.capacity() * 32), ("index_occurrences", occurrences.bytes()),
        ("index_slots", slots.capacity() * 24)]);
    Index { slots, occ: occurrences }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The index gives every kept value its occurrences in read and token order, as a reference built from the reads
    /// (excluded reads left out, the reads-per-value cap and the repeat cut), at exact capacities; the guard's view has
    /// exact capacities too.
    #[test]
    fn index_equals_the_reference() {
        for (seed, noise, repeat_fraction, max_reads) in [(1, 10, 0.01, 0), (2, 30, 0.0, 0), (3, 5, 0.05, 12), (4, 50, 0.0001, 0)] {
            let (_, tokens) = noisy_reads(noise, seed);
            let excluded: Vec<bool> = (0..tokens.len()).map(|r| (r as u64 + seed as u64) % 5 == 0).collect();
            for ex in [&[][..], &excluded[..]] {
                let b = build_index_excluding(&tokens, max_reads, repeat_fraction, ex);
                assert_eq!(kept_values(&tokens, max_reads, repeat_fraction),
                           reference_kept_values(&tokens, max_reads, repeat_fraction), "seed {}", seed);
                let want: HashMap<u32, Vec<(u32, u32)>> = reference_index(&tokens, max_reads, repeat_fraction).into_iter()
                    .map(|(v, occ)| (v, occ.into_iter().filter(|&(r, _)| ex.is_empty() || !ex[r as usize]).collect()))
                    .collect();
                let mut va: Vec<u32> = want.keys().copied().collect();
                let mut vb: Vec<u32> = b.values().collect();
                va.sort_unstable();
                vb.sort_unstable();
                assert_eq!(va, vb, "seed {}", seed);
                assert!(!va.is_empty());
                for v in &va {
                    assert_eq!(&b.get(*v), &want[v], "seed {} value {}", seed, v);
                }
                assert_eq!(b.occ.read.len(), b.occ.read.capacity());
                assert_eq!(b.occ.low.len(), b.occ.low.capacity());
            }
            let kept = kept_values(&tokens, max_reads, repeat_fraction);
            for (read, view) in tokens.iter().zip(guard_tokens(&tokens, &kept)) {
                // The compact view gives back every kept token's value and start, at exact capacities.
                let expected: Vec<GuardTok> = read.iter().filter(|t| kept.contains(&t.value))
                    .map(|t| GuardTok { value: t.value, start: t.s() as i32 }).collect();
                assert_eq!(view.iter().collect::<Vec<_>>(), expected);
                assert_eq!(view.len(), expected.len());
                assert_eq!((view.values.capacity(), view.steps.capacity()), (view.len(), view.len()));
            }
            // Steps that do not fit a u16, or go back, come back through the far starts.
            let odd = vec![Token::new(1, 5, 9), Token::new(2, 70_000, 70_009), Token::new(3, 69_990, 70_000), Token::new(4, 70_100, 70_110)];
            let view = GuardRead::new(odd.iter());
            assert_eq!(view.iter().map(|t| (t.value, t.start)).collect::<Vec<_>>(), vec![(1, 5), (2, 70_000), (3, 69_990), (4, 70_100)]);
        }
    }

    /// The flat index and the two-pass anchor grouping give the chains of the hash-map reference (`reference_chains`),
    /// guard included.
    #[test]
    fn flat_index_keeps_the_chains_of_the_hash_map_index() {
        for (seed, noise, repeat_fraction, max_reads) in [(1, 10, 0.01, 0), (2, 30, 0.0, 0), (3, 5, 0.05, 12), (4, 50, 0.0001, 0)] {
            let (lengths, tokens) = noisy_reads(noise, seed);
            let (_, raw) = noisy_reads(60, seed);
            let p = Params { max_reads, min_anchors: 3, min_overlap: 500, min_cover: 0.1, repeat_fraction, ..Params::default() };
            let g = Guard { min_anchors: 2, band: 30, band_fraction: 0.02 };
            let mut a: Vec<u32> = reference_kept_values(&raw, max_reads, repeat_fraction).into_iter().collect();
            let mut b: Vec<u32> = kept_values(&raw, max_reads, repeat_fraction).into_iter().collect();
            a.sort_unstable();
            b.sort_unstable();
            assert_eq!(a, b);
            let reference = reference_index(&tokens, max_reads, repeat_fraction);
            let index = build_index(&tokens, max_reads, repeat_fraction);
            let mut values: Vec<u32> = index.values().collect();
            values.sort_unstable();
            let mut expected: Vec<u32> = reference.keys().copied().collect();
            expected.sort_unstable();
            assert_eq!(values, expected);
            for v in &values {
                assert_eq!(index.get(*v), reference[v].as_slice());
            }
            let guard_kept = kept_values(&raw, max_reads, repeat_fraction);
            let guard_view = guard_tokens(&raw, &guard_kept);   // the guard's view
            let key = |c: &Chain| (c.donor, c.anchors, c.score.to_bits(), c.t0, c.t1, c.ts, c.te, c.ds, c.de);
            let mut scratch = Scratch::new(tokens.len());
            let mut found = 0;
            for r in 0..tokens.len() {
                for guard in [None, Some((&guard_view[..], &guard_kept, &g))] {
                    let mut x: Vec<_> = reference_chains(r, &lengths, &tokens, &reference, &p, guard).iter().map(key).collect();
                    let mut y: Vec<_> = target_chains(r, &lengths, &tokens, &index, &p, guard, &mut scratch, None).iter().map(key).collect();
                    x.sort();
                    y.sort();
                    assert_eq!(x, y, "seed {} target {}", seed, r);
                    found += y.len();
                }
            }
            assert!(found > 0);
        }
    }

    /// `Occurrences` gives back every (ordinal, token index) it was given, narrow (24-bit fields, near the bit
    /// boundaries) and wide, after copies and truncation.
    #[test]
    fn occurrences_round_trip() {
        let mut state = 11u64;
        let mut next = move || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as u32
        };
        for wide in [false, true] {
            let limit: u32 = if wide { u32::MAX } else { (1 << 24) - 1 };
            let mut pairs: Vec<(u32, u32)> = vec![(0, 0), (limit, limit), (1, 65535), (2, 65536), (3, 65537), ((1 << 24) - 1, 131071)];
            pairs.extend((0..5000).map(|_| (next() % limit.min(1 << 24), next() % (limit.min(1 << 24)))));
            if wide {
                pairs.extend((0..100).map(|_| (next() | (1 << 31), next() | (1 << 30))));
            }
            let mut o = Occurrences::with_layout(wide, 0);
            for &(d, j) in &pairs {
                o.push(d, j);
            }
            assert_eq!((0..o.len()).map(|k| o.get(k)).collect::<Vec<_>>(), pairs, "wide {}", wide);
            o.copy_within(10..20, 0);
            o.copy_one(30, 25);
            let mut expect = pairs.clone();
            expect.copy_within(10..20, 0);
            expect[25] = expect[30];
            o.truncate(4000);
            expect.truncate(4000);
            o.shrink_to_fit();
            assert_eq!((0..o.len()).map(|k| (o.read(k), o.pos(k))).collect::<Vec<_>>(), expect, "wide {}", wide);
        }
    }

    #[test]
    fn repeat_fraction_drops_the_most_frequent_values_and_zero_means_no_cap() {
        let t = |v: u32| Token::new(v, 0, 1);
        // value 7 in all 4 reads, 8 in 3, 9 and 10 in 2 reads; 11 and 12 in one read each.
        let reads = vec![vec![t(7), t(8), t(9), t(11)], vec![t(7), t(8), t(9)], vec![t(7), t(8), t(10)],
                         vec![t(7), t(10), t(12)]];
        let keys = |m: &Index| { let mut k: Vec<u32> = m.values().collect(); k.sort(); k };
        assert_eq!(keys(&build_index(&reads, 0, 0.0)), vec![7, 8, 9, 10]);
        assert_eq!(keys(&build_index(&reads, 3, 0.0)), vec![8, 9, 10]);
        // 6 distinct values: 0.34 x 6 = 2.04 -> the two most frequent (7, 8) are dropped.
        assert_eq!(keys(&build_index(&reads, 0, 0.34)), vec![9, 10]);
    }
}
