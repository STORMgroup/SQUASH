//! Segment skipping (`--skip-from`): the dirty windows of a read (where its inputs changed since the previous
//! round) and the replay restricted to them.

use super::*;

/// q spans (the target's token axis) of the layers a replay changed (chosen allele not the target's own).
pub(super) fn changed_spans(target: &[Token], chosen: &[u32]) -> Vec<(i64, i64)> {
    target.iter().zip(chosen).filter(|(_, &c)| c != 0).map(|(t, _)| (t.s(), t.e())).collect()
}

/// Sorted, merged intervals.
fn merge_intervals(mut v: Vec<(i64, i64)>) -> Vec<(i64, i64)> {
    v.sort_unstable();
    let mut out: Vec<(i64, i64)> = Vec::with_capacity(v.len());
    for (a, b) in v {
        match out.last_mut() {
            Some(last) if a <= last.1 => last.1 = last.1.max(b),
            _ => out.push((a, b)),
        }
    }
    out
}

/// A read's changed spans merged (`merge_intervals`) and stored as 32-bit pairs, at their exact length; `dirty_windows`
/// gives the same windows from them as from the unmerged spans (it pads and merges the spans, and maps a donor's
/// monotonically).
pub(super) fn pack_spans(changed: Vec<(i64, i64)>) -> Vec<(i32, i32)> {
    let narrow = |x: i64| -> i32 { <i32 as std::convert::TryFrom<i64>>::try_from(x).expect("token span outside i32") };
    let merged = merge_intervals(changed);
    let mut out = Vec::with_capacity(merged.len());
    out.extend(merged.into_iter().map(|(a, b)| (narrow(a), narrow(b))));
    out
}

/// What changed since the previous round, for `dirty_windows`.
pub struct RoundChanges<'a> {
    /// The previous round's donor rows per target.
    pub rows: &'a [Vec<DonorRow>],
    /// q spans each read's previous replay changed (on its own token axis, which corrected tokens keep), merged,
    /// 32-bit (token coordinates fit in i32).
    pub spans: &'a [Vec<(i32, i32)>],
}

/// Dirty windows of target `r` (its token q axis), padded by `pad_q` and merged: where its inputs changed since the
/// previous round. (1) Its own tokens the previous replay changed. (2) A partner's changed tokens inside their overlap,
/// mapped linearly through the donor row. (3) A new partner's whole overlap, a dropped partner's old overlap, and the
/// stretch by which a kept partner's overlap moved at either end. Row bounds are on the source-symbol axis
/// (`symbols`), tokens on the token axis (`lengths`), converted as `replay_read` converts them.
pub fn dirty_windows(
    r: usize,
    rows: &[DonorRow],
    changes: &RoundChanges,
    lengths: &[i64],
    symbols: &[i64],
    pad_q: i64,
) -> Vec<(i64, i64)> {
    let to_q = |x: i64| x * lengths[r] / symbols[r].max(1);
    let mut dirty: Vec<(i64, i64)> = changes.spans[r].iter().map(|&(a, b)| (a as i64, b as i64)).collect();
    let previous: HashMap<usize, &DonorRow> = changes.rows[r].iter().map(|x| (x.donor, x)).collect();
    let mut kept: HashSet<usize> = HashSet::new();
    for row in rows {
        match previous.get(&row.donor) {
            None => dirty.push((to_q(row.target_start), to_q(row.target_end))),
            Some(old) => {
                kept.insert(row.donor);
                if old.target_start != row.target_start {
                    dirty.push((to_q(old.target_start.min(row.target_start)), to_q(old.target_start.max(row.target_start))));
                }
                if old.target_end != row.target_end {
                    dirty.push((to_q(old.target_end.min(row.target_end)), to_q(old.target_end.max(row.target_end))));
                }
            }
        }
        let (d, span_d, span_t) = (row.donor, row.donor_end - row.donor_start, row.target_end - row.target_start);
        if span_d <= 0 {
            continue;
        }
        for &(a, b) in &changes.spans[d] {
            let (a, b) = (a as i64, b as i64);
            // Donor token axis -> donor symbols -> target symbols -> target token axis.
            let (a, b) = (a * symbols[d] / lengths[d].max(1), b * symbols[d] / lengths[d].max(1));
            if b < row.donor_start || a > row.donor_end {
                continue;
            }
            let map = |x: i64| row.target_start + (x.clamp(row.donor_start, row.donor_end) - row.donor_start) * span_t / span_d;
            dirty.push((to_q(map(a)), to_q(map(b))));
        }
    }
    for old in &changes.rows[r] {
        if !kept.contains(&old.donor) {
            dirty.push((to_q(old.target_start), to_q(old.target_end)));
        }
    }
    merge_intervals(dirty.into_iter().map(|(a, b)| (a - pad_q, b + pad_q)).collect())
}

/// One read's replay restricted to its dirty windows (segment skipping). Layers outside the windows keep the
/// target's own token. Each run of dirty layers is replayed with `context` layers on both sides and every donor row
/// clipped to that stretch (donor bounds mapped linearly), so a partner is aligned only there; the chosen alleles of
/// the dirty layers are kept. Returns the corrected tokens, the number of dirty layers changed, the changed q spans and
/// the number of layers replayed.
#[allow(clippy::too_many_arguments)]
pub fn replay_windows(
    r: usize,
    windows: &[(i64, i64)],
    rows: &[DonorRow],
    tokens: &[Vec<Token>],
    lengths: &[i64],
    symbols: &[i64],
    pad: usize,
    band: usize,
    context: usize,
) -> (Vec<Token>, usize, Vec<(i64, i64)>, usize) {
    let target = &tokens[r];
    let n = target.len();
    let own: Vec<Allele> = target.iter().map(|t| Allele { insertion: vec![], site: Some(t.value) }).collect();
    let mut chosen: Vec<Allele> = own.clone();
    // Dirty layers: tokens whose span meets a window (both sorted, one pass).
    let mut dirty = vec![false; n];
    let mut w = 0;
    for (i, t) in target.iter().enumerate() {
        while w < windows.len() && windows[w].1 < t.s() {
            w += 1;
        }
        dirty[i] = w < windows.len() && windows[w].0 < t.e();
    }
    // Runs of dirty layers; runs whose contexts would overlap are replayed together.
    let mut runs: Vec<(usize, usize)> = Vec::new();
    let mut i = 0;
    while i < n {
        if !dirty[i] {
            i += 1;
            continue;
        }
        let mut j = i;
        while j + 1 < n && dirty[j + 1] {
            j += 1;
        }
        match runs.last_mut() {
            Some(last) if i <= last.1 + 2 * context + 1 => last.1 = j,
            _ => runs.push((i, j)),
        }
        i = j + 1;
    }
    let (mut changed, mut replayed) = (0, 0);
    let mut spans: Vec<(i64, i64)> = Vec::new();
    let (len, sym) = (lengths[r], symbols[r].max(1));
    for &(i0, i1) in &runs {
        let (c0, c1) = (i0.saturating_sub(context), (i1 + 1 + context).min(n));
        let sub = &target[c0..c1];
        replayed += c1 - c0;
        // The stretch on the symbol axis, and every row clipped to it.
        let (s0, s1) = (sub[0].s() * sym / len, ceil_ratio(sub[c1 - c0 - 1].e(), sym, len));
        let clipped: Vec<DonorRow> = rows
            .iter()
            .filter_map(|row| {
                let (a, b) = (row.target_start.max(s0), row.target_end.min(s1));
                let span_t = row.target_end - row.target_start;
                if a >= b || span_t <= 0 {
                    return None;
                }
                let span_d = row.donor_end - row.donor_start;
                Some(DonorRow {
                    donor: row.donor,
                    target_start: a,
                    target_end: b,
                    donor_start: row.donor_start + (a - row.target_start) * span_d / span_t,
                    donor_end: row.donor_start + ceil_ratio(b - row.target_start, span_d, span_t),
                })
            })
            .collect();
        let replay = replay_read(sub, len, symbols[r], &clipped, tokens, lengths, symbols, pad, band);
        for i in i0..=i1 {
            let p = i - c0;
            let c = replay.chosen[p];
            chosen[i] = replay.layers[p].alleles[c as usize].0.clone();
            if c != 0 {
                changed += 1;
                spans.push((target[i].s(), target[i].e()));
            }
        }
    }
    let corrected = spanned_from_alleles(target, chosen.iter(), len);
    (corrected, changed, spans, replayed)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A random cohort of noisy copies of one sequence (unit spans), all-against-all rows over whole reads.
    fn noisy_cohort(n: usize, len: usize, seed: u64) -> (Vec<Vec<Token>>, Vec<i64>, Vec<Vec<DonorRow>>) {
        let mut x = seed;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let truth: Vec<u32> = (0..len).map(|_| (next() % 1000) as u32).collect();
        let reads: Vec<Vec<u32>> = (0..n)
            .map(|_| {
                let mut r = Vec::new();
                for &t in &truth {
                    match next() % 20 {
                        0 => {}                                    // deletion
                        1 => r.extend([t, (next() % 1000) as u32]), // insertion
                        2 => r.push((next() % 1000) as u32),       // substitution
                        _ => r.push(t),
                    }
                }
                r
            })
            .collect();
        let tokens: Vec<Vec<Token>> = reads.iter().map(|r| unit(r)).collect();
        let lengths: Vec<i64> = tokens.iter().map(|t| t.len() as i64).collect();
        let rows = (0..n)
            .map(|r| {
                (0..n)
                    .filter(|&d| d != r)
                    .map(|d| DonorRow { donor: d, target_start: 0, target_end: lengths[r], donor_start: 0, donor_end: lengths[d] })
                    .collect()
            })
            .collect();
        (tokens, lengths, rows)
    }

    #[test]
    fn segment_replay_of_a_whole_read_is_the_full_replay() {
        let (tokens, lengths, rows) = noisy_cohort(8, 120, 7);
        for r in 0..tokens.len() {
            let full = replay_read(&tokens[r], lengths[r], lengths[r], &rows[r], &tokens, &lengths, &lengths, 8, 0);
            let (corrected, changed, spans, replayed) =
                replay_windows(r, &[(-1000, 1000)], &rows[r], &tokens, &lengths, &lengths, 8, 0, 16);
            assert_eq!(corrected, spanned_tokens(&tokens[r], &full, lengths[r]));
            assert_eq!(changed, full.changed_layers());
            assert_eq!(spans, changed_spans(&tokens[r], &full.chosen));
            assert_eq!(replayed, tokens[r].len());
            // No window: the read is left as it is.
            let (same, changed, spans, replayed) =
                replay_windows(r, &[], &rows[r], &tokens, &lengths, &lengths, 8, 0, 16);
            assert_eq!((same, changed, spans.len(), replayed), (tokens[r].clone(), 0, 0, 0));
        }
    }

    #[test]
    fn a_window_corrects_its_own_layers_as_the_full_replay_does() {
        // One error in the middle of read 0, the others exact: a window around it gives the full replay's tokens.
        let truth: Vec<u32> = (0..200).map(|i| (i * 37 % 1000) as u32).collect();
        let mut reads: Vec<Vec<u32>> = vec![truth.clone(); 6];
        reads[0][100] = 5000;
        let tokens: Vec<Vec<Token>> = reads.iter().map(|r| unit(r)).collect();
        let lengths: Vec<i64> = tokens.iter().map(|t| t.len() as i64).collect();
        let rows: Vec<DonorRow> = (1..6)
            .map(|d| DonorRow { donor: d, target_start: 0, target_end: 200, donor_start: 0, donor_end: 200 })
            .collect();
        let full = replay_read(&tokens[0], 200, 200, &rows, &tokens, &lengths, &lengths, 8, 0);
        let (corrected, changed, spans, replayed) =
            replay_windows(0, &[(90, 110)], &rows, &tokens, &lengths, &lengths, 8, 0, 16);
        assert_eq!(corrected, spanned_tokens(&tokens[0], &full, 200));
        assert_eq!(corrected.iter().map(|t| t.value).collect::<Vec<_>>(), truth);
        assert_eq!((changed, spans, replayed), (1, vec![(100, 101)], 21 + 32));   // layers 90..=110 and 16 on each side
    }

    #[test]
    fn dirty_windows_follow_changed_inputs() {
        // Read 0 (100 q) with partners 1 and 2 (100 q each); lengths = symbols. Previous round: partner 1 over
        // 0..60, partner 3 over 50..100 (dropped now); now partner 1 over 0..70 and partner 2 (new) over 80..100.
        // Read 0 changed 10..11; partner 1 changed its 30..31 (target 30..31).
        let row = |donor, ts, te, ds, de| DonorRow { donor, target_start: ts, target_end: te, donor_start: ds, donor_end: de };
        let lengths = vec![100i64; 4];
        let previous = vec![vec![row(1, 0, 60, 0, 60), row(3, 50, 100, 0, 50)], vec![], vec![], vec![]];
        let spans = vec![vec![(10, 11)], vec![(30, 31)], vec![], vec![(0, 1)]];
        let changes = RoundChanges { rows: &previous, spans: &spans };
        let now = vec![row(1, 0, 70, 0, 70), row(2, 80, 100, 0, 20)];
        assert_eq!(dirty_windows(0, &now, &changes, &lengths, &lengths, 0),
                   vec![(10, 11), (30, 31), (50, 100)]);
        assert_eq!(dirty_windows(0, &now, &changes, &lengths, &lengths, 5), vec![(5, 16), (25, 36), (45, 105)]);
        // Nothing changed: no window.
        let quiet = vec![vec![], vec![], vec![], vec![]];
        let changes = RoundChanges { rows: &previous, spans: &quiet };
        let same = vec![row(1, 0, 60, 0, 60), row(3, 50, 100, 0, 50)];
        assert!(dirty_windows(0, &same, &changes, &lengths, &lengths, 150).is_empty());
    }

    /// Segment skip keeps each read's changed spans merged; `dirty_windows` gives the same windows as
    /// from one span per changed token (random reads, partners and overlapping token spans; pads 0 and 150).
    #[test]
    fn merged_changed_spans_give_the_same_dirty_windows() {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let mut next = |m: i64| -> i64 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % m.max(1) as u64) as i64
        };
        for _ in 0..300 {
            let n = 6usize;
            let lengths: Vec<i64> = (0..n).map(|_| 200 + next(400)).collect();
            let symbols: Vec<i64> = lengths.iter().map(|&l| l * (1 + next(3)) + next(50)).collect();
            let spans: Vec<Vec<(i64, i64)>> = (0..n)
                .map(|r| (0..next(30)).map(|_| { let s = next(lengths[r]); (s, s + 1 + next(12)) }).collect())
                .collect();
            // one span per token, unmerged, and as segment skip stores them (`pack_spans`)
            let merged: Vec<Vec<(i32, i32)>> = spans.iter().map(|v| pack_spans(v.clone())).collect();
            let spans: Vec<Vec<(i32, i32)>> = spans.iter().map(|v| v.iter().map(|&(a, b)| (a as i32, b as i32)).collect()).collect();
            let rows = |next: &mut dyn FnMut(i64) -> i64| -> Vec<Vec<DonorRow>> {
                let mut out = vec![Vec::new(); n];
                for r in 0..n {
                    for d in 0..n {
                        if d == r || next(3) == 0 {
                            continue;
                        }
                        let (ts, ds) = (next(symbols[r] / 2), next(symbols[d] / 2));
                        let (te, de) = (ts + 1 + next(symbols[r] - ts), ds + 1 + next(symbols[d] - ds));
                        out[r].push(DonorRow { donor: d, target_start: ts, target_end: te, donor_start: ds, donor_end: de });
                    }
                }
                out
            };
            let (previous, now) = (rows(&mut next), rows(&mut next));
            for pad in [0, 150] {
                for r in 0..n {
                    let a = dirty_windows(r, &now[r], &RoundChanges { rows: &previous, spans: &spans }, &lengths, &symbols, pad);
                    let b = dirty_windows(r, &now[r], &RoundChanges { rows: &previous, spans: &merged }, &lengths, &symbols, pad);
                    assert_eq!(a, b);
                }
            }
        }
    }
}
