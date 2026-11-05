//! The per-donor part of a target's search (`DonorChainer`): the prune state, the chain's overlap and cover tests,
//! the guard and the q-symbol check (cached, or deferred to the coverage cap).

use super::*;

/// The per-donor part of `target_chains` and `target_chains_free`: the prune state, the target's lazily built guard
/// map, and the chains kept.
pub(super) struct DonorChainer<'a> {
    r: usize,
    lengths: &'a [i64],
    tokens: &'a [Vec<Token>],
    p: &'a Params,
    guard: Option<(&'a [GuardRead], &'a HashSet<u32>, &'a Guard)>,
    prune: Option<usize>,
    prof: bool,
    // Pruning state: the chains kept so far (score, t0, t1) not yet counted, in a max-heap by score (each counted once,
    // when the bound first falls below its score); per-position counts of counted chains; and the positions with fewer
    // than `coverage` of them ("open"), as "next open position" links with path compression (a position only ever
    // closes).
    pending: std::collections::BinaryHeap<ByScore>,
    depth: Vec<usize>,
    next_open: Vec<u32>,
    // Built on the first chain that reaches the guard.
    target_starts: Option<HashMap<u32, Vec<i64>>>,
    pub(super) chains: Vec<Chain>,
}

impl<'a> DonorChainer<'a> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(r: usize, lengths: &'a [i64], tokens: &'a [Vec<Token>], p: &'a Params,
           guard: Option<(&'a [GuardRead], &'a HashSet<u32>, &'a Guard)>, prune: Option<usize>, prof: bool) -> DonorChainer<'a> {
        let n_layers = tokens[r].len();
        let (mut depth, mut next_open) = (Vec::new(), Vec::new());
        if let Some(coverage) = prune {
            depth = vec![0; n_layers];
            next_open = (0..=n_layers as u32).map(|x| if coverage == 0 && (x as usize) < n_layers { x + 1 } else { x }).collect();
        }
        DonorChainer { r, lengths, tokens, p, guard, prune, prof, pending: std::collections::BinaryHeap::new(), depth,
                       next_open, target_starts: None, chains: Vec::new() }
    }

    /// The prune (`target_chains`): whether a donor whose chain can score at most `b` and whose anchors span target
    /// positions lo..=hi is skipped, every position there being covered already. `n_anchors`: for the profile.
    pub(super) fn pruned(&mut self, b: usize, lo: usize, hi: usize, n_anchors: usize) -> bool {
        let coverage = match self.prune {
            Some(c) => c,
            None => return false,
        };
        let pclock = if self.prof { Some(std::time::Instant::now()) } else { None };
        let c = b as f64;
        // Count every kept chain scoring more than this donor's bound (bounds only fall from here on).
        while self.pending.peek().map_or(false, |top| top.0 > c) {
            let ByScore(_, t0, t1) = self.pending.pop().unwrap();
            for x in t0..=t1 {
                self.depth[x] += 1;
                if self.depth[x] == coverage {
                    self.next_open[x] = x as u32 + 1;   // closed: link to the next position
                }
            }
        }
        let closed = open_from(&mut self.next_open, lo) > hi;
        if let Some(t) = pclock {
            profile::since(17, t);
        }
        if closed && self.prof {
            profile::add(14, 1);
            profile::add(15, n_anchors as u64);
        }
        closed
    }

    /// Chains donor `d`'s anchors (`group`) and keeps the chain when it passes the overlap and cover tests, the guard
    /// and the q-symbol alignment check (its cached verdict, or the check deferred to the coverage cap: the lazy
    /// check).
    pub(super) fn chain<A: Anchor>(&mut self, d: usize, group: &[A]) {
        let (r, p, prof, tokens, lengths, guard) = (self.r, self.p, self.prof, self.tokens, self.lengths, self.guard);
        // `cover_probe` (a diagnostic): what an exact bound would skip; every donor is still chained.
        let probe = cover_probe::on();
        let min_cover = if p.signal_min_cover >= 0.0 && signal::get().is_some() { p.signal_min_cover } else { p.min_cover };
        let skippable = if probe {
            let t = std::time::Instant::now();
            let b = cover_bound(&tokens[r], &tokens[d], lengths[r], lengths[d], group, min_cover, p.min_overlap);
            cover_probe::add(9, t.elapsed().as_nanos() as u64);
            cover_probe::add(0, 1);
            cover_probe::add(3, group.len() as u64);
            if let Some(k) = b {
                cover_probe::add(1 + k as usize, 1);
                cover_probe::add(4, group.len() as u64);
            }
            b.is_some()
        } else {
            false
        };
        let clock = std::time::Instant::now();
        let sorted = group.windows(2).all(|w| w[0] < w[1]);
        let chained = if sorted {
            chain_sorted(group, p)
        } else {
            chain_anchors(group.iter().map(|a| (a.t(), a.d())).collect(), p)
        };
        if probe {
            let ns = clock.elapsed().as_nanos() as u64;
            cover_probe::add(5, ns);
            if skippable {
                cover_probe::add(6, ns);
            }
        }
        if prof {
            profile::since(3, clock);
            profile::add(10, 1);
        }
        let clock = std::time::Instant::now();
        let best = match chained {
            Some(b) if b.0 >= p.min_anchors => b,
            _ => return,
        };
        let (n, score, t0, t1, d0, d1) = best;
        let (tt, dt) = (&tokens[r], &tokens[d]);
        let (t0, t1, d0, d1) = (t0 as usize, t1 as usize, d0 as usize, d1 as usize);
        let offset = round_half_even_div2((dt[d0].s() - tt[t0].s()) + (dt[d1].s() - tt[t1].s()));
        let ts = 0i64.max(-offset);
        let te = lengths[r].min(lengths[d] - offset);
        if te - ts < p.min_overlap || ((tt[t1].e() - tt[t0].s()) as f64) < min_cover * ((te - ts) as f64) {
            if prof {
                profile::since(4, clock);
            }
            return;
        }
        if probe {
            cover_probe::add(8, 1);
            if skippable {
                cover_probe::add(7, 1);
            }
        }
        if prof {
            profile::since(4, clock);
            profile::add(11, 1);
        }
        if let Some((raw, guard_kept, g)) = guard {
            let clock = std::time::Instant::now();
            let starts = self.target_starts.get_or_insert_with(|| guard_starts(raw[r].iter(), guard_kept));
            let pass = guard_supported(starts, raw[d].iter(), ts, te, offset, g);
            if prof {
                profile::since(5, clock);
            }
            if !pass {
                return;
            }
        }
        // The q-symbol alignment check (`signal`), before the coverage cap: a pair without a cached verdict keeps its
        // chain with the check deferred to the cap (`target_rows`), which makes it only for a chain that still reaches
        // an open position, in the cap's order (the lazy check, `LAZY_CHECK`). Such a chain is not counted by the prune
        // (its verdict is unknown), so the prune only counts checked chains.
        let mut deferred: Option<Box<Deferred>> = None;
        if let Some(sc) = signal::get() {
            let line = (tt[t0].s(), dt[d0].s(), tt[t1].s(), dt[d1].s());
            match sc.cached(r, d) {
                Some(false) => return,
                Some(true) => {}
                None => {
                    // the chain's anchors in q (the start of each matched token on both reads), from the DP
                    // `chain_sorted` just ran on this group when it did (no second DP)
                    let members = if sorted { chain_members_after_sorted(group) } else { None };
                    let q: Vec<(i64, i64)> = members.unwrap_or_else(|| chain_members(group, p)).into_iter()
                        .map(|(t, d)| (tt[t as usize].s(), dt[d as usize].s())).collect();
                    LAZY_CHECK[0].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    deferred = Some(Box::new(Deferred { line, ts, te, q }));
                }
            }
        }
        if prof && deferred.is_none() {
            profile::add(12, 1);
        }
        if self.prune.is_some() && deferred.is_none() {
            self.pending.push(ByScore(score, t0, t1));
        }
        self.chains.push(Chain { donor: d, anchors: n, score, t0, t1, ts, te, ds: ts + offset, de: te + offset, check: deferred });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The open-position links answer "first open position at or after x" as a scan of the depths does, while positions
    /// close in any order.
    #[test]
    fn open_position_links_match_a_scan() {
        let mut state = 11u64;
        let mut next = move || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as usize
        };
        for n in [1usize, 2, 7, 50, 300] {
            let mut links: Vec<u32> = (0..=n as u32).collect();
            let mut open = vec![true; n];
            for _ in 0..4 * n {
                let x = next() % n;
                open[x] = false;
                links[x] = x as u32 + 1;
                for _ in 0..3 {
                    let q = next() % n;
                    let want = (q..n).find(|&y| open[y]).unwrap_or(n);
                    assert_eq!(open_from(&mut links, q), want, "n {} q {}", n, q);
                }
            }
        }
    }
}
