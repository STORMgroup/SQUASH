//! One target's chains without the index (`target_chains_free`), for searches restricted to a few donors per target.

use super::*;

/// `target_chains` without the index. The correction rounds' partner searches from `--hops-from` restrict each target's
/// donors to its two-hop candidates (`allowed`, a few dozen), so a whole index of every read's occurrences would cost
/// memory and time for little. Instead `kept` = the values the index would keep (`kept_seed_values`), and each allowed
/// donor's anchors come from looking its tokens up in a small table of the target's kept values. A donor gets the
/// anchors the index gives it ((target index, donor index) of equal kept values, sorted), the donors the index's order
/// (first touched at their first matching target position, then by ordinal; then by decreasing anchor count) and the
/// same prune decisions, and its anchors are placed only when it is chained (no batches). So the same chains
/// (`target_chains_free_equals_the_index`).
#[allow(clippy::too_many_arguments)]
pub fn target_chains_free(
    r: usize,
    lengths: &[i64],
    tokens: &[Vec<Token>],
    kept: &HashSet<u32>,
    excluded: &[bool],
    p: &Params,
    guard: Option<(&[GuardRead], &HashSet<u32>, &Guard)>,
    scratch: &mut FreeScratch,
    prune: Option<usize>,
    allowed: &[u32],
) -> Vec<Chain> {
    let FreeScratch { tv, groups, table, cand, anchors } = scratch;
    let prof = profile::on();
    let clock = std::time::Instant::now();
    // The target's kept tokens by value, then position; one group per value (value, first and last position, range).
    tv.clear();
    tv.extend(tokens[r].iter().enumerate().filter(|(_, t)| kept.contains(&t.value)).map(|(i, t)| (t.value, i as u32)));
    tv.sort_unstable();
    groups.clear();
    let mut k = 0;
    while k < tv.len() {
        let mut e = k + 1;
        while e < tv.len() && tv[e].0 == tv[k].0 {
            e += 1;
        }
        groups.push((tv[k].0, tv[k].1, tv[e - 1].1, k as u32, e as u32));
        k = e;
    }
    // Open addressing over the groups (u32::MAX = empty), at most half full.
    let bits = (2 * groups.len().max(8)).next_power_of_two().trailing_zeros();
    let mask = (1usize << bits) - 1;
    let slot = |v: u32| ((v as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - bits)) as usize;
    table.clear();
    table.resize(mask + 1, u32::MAX);
    for (g, gr) in groups.iter().enumerate() {
        let mut s = slot(gr.0);
        while table[s] != u32::MAX {
            s = (s + 1) & mask;
        }
        table[s] = g as u32;
    }
    let lookup = |v: u32| -> Option<usize> {
        let mut s = slot(v);
        loop {
            let g = table[s];
            if g == u32::MAX {
                return None;
            }
            if groups[g as usize].0 == v {
                return Some(g as usize);
            }
            s = (s + 1) & mask;
        }
    };
    // Each allowed donor's anchors counted: (first, last matching target position, donor, anchors).
    cand.clear();
    for &d in allowed {
        let du = d as usize;
        if du == r || (!excluded.is_empty() && excluded[du]) {
            continue;
        }
        let (mut count, mut first, mut last) = (0u64, u32::MAX, 0u32);
        for t in &tokens[du] {
            if let Some(g) = lookup(t.value) {
                let gr = &groups[g];
                count += (gr.4 - gr.3) as u64;
                first = first.min(gr.1);
                last = last.max(gr.2);
            }
        }
        if count > 0 {
            cand.push((first, d, last, count));
        }
    }
    // The index's donor order: first touch (first matching position, then ordinal), then by decreasing anchor count.
    cand.sort_unstable_by_key(|c| (c.0, c.1));
    if prof {
        let placed: u64 = cand.iter().filter(|c| c.3 as usize >= p.min_anchors).map(|c| c.3).sum();
        profile::since(2, clock);
        profile::add(7, 1);
        profile::add(8, placed);
        profile::add(9, cand.len() as u64);
    }
    let oclock = std::time::Instant::now();
    // (donor, first and last matching target position, anchors (the bound), anchors)
    let mut order: Vec<(usize, usize, usize, usize, usize)> = cand.iter().filter(|c| c.3 as usize >= p.min_anchors)
        .map(|c| (c.1 as usize, c.0 as usize, c.2 as usize, c.3 as usize, c.3 as usize)).collect();
    if prune.is_some() {
        order.sort_by(|a, b| b.3.cmp(&a.3));
    }
    let mut dc = DonorChainer::new(r, lengths, tokens, p, guard, prune, prof);
    if prof {
        profile::since(16, oclock);
    }
    for &(d, lo, hi, b, n) in &order {
        if dc.pruned(b, lo, hi, n) {
            continue;
        }
        let pclock = if prof { Some(std::time::Instant::now()) } else { None };
        anchors.clear();
        for (j, t) in tokens[d].iter().enumerate() {
            if let Some(g) = lookup(t.value) {
                let gr = &groups[g];
                anchors.extend(tv[gr.3 as usize..gr.4 as usize].iter().map(|&(_, i)| (i, j as u32)));
            }
        }
        anchors.sort_unstable();
        if let Some(t) = pclock {
            profile::since(21, t);
        }
        dc.chain(d, anchors);
    }
    dc.chains
}

/// `target_chains_free`'s per-thread buffers (a few per target, nothing per read).
#[derive(Default)]
pub struct FreeScratch {
    /// the target's kept tokens (value, position), sorted
    tv: Vec<(u32, u32)>,
    /// per distinct kept value: (value, first position, last position, range in `tv`)
    groups: Vec<(u32, u32, u32, u32, u32)>,
    table: Vec<u32>,
    /// per donor with an anchor: (first, donor, last, anchors)
    cand: Vec<(u32, u32, u32, u64)>,
    anchors: Vec<(u32, u32)>,
}

impl FreeScratch {
    /// After a target, the buffers grown past `SCRATCH_KEEP` bytes are freed (as `Scratch::trim`).
    pub fn trim(&mut self) {
        fn keep<T>(v: &mut Vec<T>) {
            if v.capacity() * std::mem::size_of::<T>() > SCRATCH_KEEP {
                *v = Vec::new();
            }
        }
        keep(&mut self.tv);
        keep(&mut self.groups);
        keep(&mut self.table);
        keep(&mut self.cand);
        keep(&mut self.anchors);
        CHAIN_SCRATCH.with(|cell| {
            let mut s = cell.borrow_mut();
            keep(&mut s.score);
            keep(&mut s.back);
            keep(&mut s.back8);
            s.last = (0, 0, 0);
        });
    }
}

/// The values `build_index_excluding` keeps as seeds (`kept_runs`), counted over every read (`target_chains_free`).
pub fn kept_seed_values(tokens: &[Vec<Token>], p: &Params) -> HashSet<u32> {
    kept_values(tokens, p.max_reads, p.repeat_fraction)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Without the index (`target_chains_free`) a target restricted to allowed donors gets the chains of the index
    /// search over an index of those donors alone, in the same order: repeat and noisy reads, both prune bounds, with
    /// and without the prune, guard or not, excluded reads, the options that change only the kept values (reads cap,
    /// sparse seeds, repeat fraction), batched placement on the index side; and the packed search gives the plain
    /// search's rows.
    #[test]
    fn target_chains_free_equals_the_index() {
        let key = |c: &Chain| (c.donor, c.anchors, c.score.to_bits(), c.t0, c.t1, c.ts, c.te, c.ds, c.de);
        let mut state = 11u64;
        let mut next = move || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as usize
        };
        let (mut compared, mut pruned_differs) = (0usize, false);
        for (seed, noise, repeats) in [(1, 10, true), (2, 30, true), (3, 5, false), (4, 20, false)] {
            let (lengths, tokens) = if repeats { repeat_reads(noise, seed) } else { noisy_reads(noise, seed) };
            let n = tokens.len();
            let allowed: Vec<Vec<u32>> =
                (0..n).map(|t| (0..n as u32).filter(|&d| d as usize != t && next() % 3 != 0).collect()).collect();
            let excluded: Vec<bool> = (0..n).map(|r| r % 11 == 5).collect();
            let guard_kept = kept_values(&tokens, 0, 0.0);
            let raw = guard_tokens(&tokens, &guard_kept);
            let g = Guard::default();
            let settings = [
                (Params { max_reads: 0, min_anchors: 3, min_overlap: 500, min_cover: 0.1, coverage: 3, ..Params::default() }, false),
                (Params { max_reads: 20, min_anchors: 4, min_overlap: 300, min_cover: 0.2, repeat_fraction: 0.01, coverage: 1,
                          ..Params::default() }, true),
                (Params { max_reads: 0, min_anchors: 3, min_overlap: 300, min_cover: 0.1, coverage: 2, ..Params::default() }, false),
            ];
            for (p, with_excluded) in settings.iter() {
                let excl: &[bool] = if *with_excluded { &excluded } else { &[] };
                // per target, the index of its allowed donors (the others excluded; the same kept values)
                let indexes: Vec<Index> = (0..n).map(|r| {
                    let not_allowed: Vec<bool> = (0..n).map(|d| (!excl.is_empty() && excl[d]) || !allowed[r].contains(&(d as u32)))
                        .collect();
                    build_index_excluding(&tokens, p.max_reads, p.repeat_fraction, &not_allowed)
                }).collect();
                let kept = kept_seed_values(&tokens, p);
                for guard in [None, Some((raw.as_slice(), &guard_kept, &g))] {
                    let (mut s, mut f) = (Scratch::new(n), FreeScratch::default());
                    for r in 0..n {
                        if *with_excluded && excluded[r] {
                            continue;
                        }
                        let mut both = Vec::new();
                        for prune in [None, Some(p.coverage)] {
                            for batch in [BATCH_ANCHORS, 40] {
                                TEST_BATCH_ANCHORS.with(|c| c.set(batch));
                                let a = target_chains(r, &lengths, &tokens, &indexes[r], p, guard, &mut s, prune);
                                TEST_BATCH_ANCHORS.with(|c| c.set(BATCH_ANCHORS));
                                let b = target_chains_free(r, &lengths, &tokens, &kept, excl, p, guard, &mut f, prune, &allowed[r]);
                                assert_eq!(a.iter().map(key).collect::<Vec<_>>(), b.iter().map(key).collect::<Vec<_>>(),
                                           "target {} prune {:?} batch {}", r, prune, batch);
                                compared += a.len();
                                both.push(a.len());
                            }
                        }
                        pruned_differs |= both[0] != both[2];
                    }
                }
                let ids: Vec<String> = (0..n).map(|r| format!("read{:02}", r)).collect();
                let pool = rayon::ThreadPoolBuilder::new().num_threads(2).build().unwrap();
                let guard = Some((raw.as_slice(), &guard_kept, &g));
                let plain = partner_rows_ext(&ids, &lengths, &tokens, guard, p, excl, &pool, Some(&allowed), None);
                assert_eq!(partner_rows_packed(&ids, &lengths, &tokens, guard, p, excl, &pool, Some(&allowed), None).unpack(), plain);
            }
        }
        assert!(compared > 1000, "{} chains compared", compared);
        assert!(pruned_differs, "the prune never dropped a chain: the test would not see a wrong prune order");
    }
}
