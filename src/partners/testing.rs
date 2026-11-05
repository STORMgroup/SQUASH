//! Shared helpers of the module's tests: reference implementations and synthetic reads.

use super::*;

/// The token values `build_index` keeps: values seen in 2..=max_reads distinct
/// reads (0: no cap), after dropping the int(repeat_fraction x distinct values)
/// values seen in the most reads (ties: smaller value first).
pub(super) fn reference_kept_values(tokens: &[Vec<Token>], max_reads: usize, repeat_fraction: f64) -> HashSet<u32> {
    // value -> (last read seen, distinct reads); reads are visited in order.
    let mut n_reads: HashMap<u32, (u32, usize)> = HashMap::new();
    for (r, read) in tokens.iter().enumerate() {
        for t in read {
            let entry = n_reads.entry(t.value).or_insert((u32::MAX, 0));
            if entry.0 != r as u32 {
                *entry = (r as u32, entry.1 + 1);
            }
        }
    }
    let mut dropped: HashSet<u32> = HashSet::new();
    if repeat_fraction > 0.0 {
        let mut ranked: Vec<(usize, u32)> = n_reads.iter().map(|(&v, &(_, n))| (n, v)).collect();
        ranked.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        let n_drop = (repeat_fraction * ranked.len() as f64) as usize;
        dropped.extend(ranked[..n_drop.min(ranked.len())].iter().map(|&(_, v)| v));
    }
    n_reads
        .into_iter()
        .filter(|(value, (_, reads))| {
            !dropped.contains(value) && *reads >= 2 && (max_reads == 0 || *reads <= max_reads)
        })
        .map(|(value, _)| value)
        .collect()
}

/// `build_index`: token value -> occurrences (read ordinal, token index) for the
/// `kept_values`.
pub(super) fn reference_index(tokens: &[Vec<Token>], max_reads: usize, repeat_fraction: f64) -> HashMap<u32, Vec<(u32, u32)>> {
    let kept = reference_kept_values(tokens, max_reads, repeat_fraction);
    let mut index: HashMap<u32, Vec<(u32, u32)>> = HashMap::new();
    for (r, read) in tokens.iter().enumerate() {
        for (i, t) in read.iter().enumerate() {
            if kept.contains(&t.value) {
                index.entry(t.value).or_default().push((r as u32, i as u32));
            }
        }
    }
    index
}

/// The reference of `target_chains` (every donor's anchors collected in a hash map, no prune): one target's chains.
/// `guard`: the uncorrected tokens, their `kept_values` and the rule.
pub(super) fn reference_chains(
    r: usize,
    lengths: &[i64],
    tokens: &[Vec<Token>],
    index: &HashMap<u32, Vec<(u32, u32)>>,
    p: &Params,
    guard: Option<(&[GuardRead], &HashSet<u32>, &Guard)>,
) -> Vec<Chain> {
    let mut pairs: HashMap<usize, Vec<(i64, i64)>> = HashMap::new();
    for (i, t) in tokens[r].iter().enumerate() {
        if let Some(occ) = index.get(&t.value) {
            for &(d, j) in occ {
                if d as usize != r {
                    pairs.entry(d as usize).or_default().push((i as i64, j as i64));
                }
            }
        }
    }
    // Built on the first chain that reaches the guard.
    let mut target_starts: Option<HashMap<u32, Vec<i64>>> = None;
    let mut chains = Vec::new();
    for (d, anchors) in pairs {
        if anchors.len() < p.min_anchors {
            continue;
        }
        let best = match chain_anchors(anchors, p) {
            Some(b) if b.0 >= p.min_anchors => b,
            _ => continue,
        };
        let (n, score, t0, t1, d0, d1) = best;
        let (tt, dt) = (&tokens[r], &tokens[d]);
        let (t0, t1, d0, d1) = (t0 as usize, t1 as usize, d0 as usize, d1 as usize);
        let offset = round_half_even_div2((dt[d0].s() - tt[t0].s()) + (dt[d1].s() - tt[t1].s()));
        let ts = 0i64.max(-offset);
        let te = lengths[r].min(lengths[d] - offset);
        if te - ts < p.min_overlap || ((tt[t1].e() - tt[t0].s()) as f64) < p.min_cover * ((te - ts) as f64) {
            continue;
        }
        if let Some((raw, guard_kept, g)) = guard {
            let starts = target_starts.get_or_insert_with(|| guard_starts(raw[r].iter(), guard_kept));
            if !guard_supported(starts, raw[d].iter(), ts, te, offset, g) {
                continue;
            }
        }
        chains.push(Chain { donor: d, anchors: n, score, t0, t1, ts, te, ds: ts + offset, de: te + offset, check: None });
    }
    chains
}

/// Reads of `noisy_reads`' kind over a genome with a tandem array (a 7-token unit, one copy in 15 varied, over 1,000
/// tokens): many anchors per target position and per donor position, as in a centromere.
pub(super) fn repeat_reads(noise: u32, seed: u64) -> (Vec<i64>, Vec<Vec<Token>>) {
    let mut state = seed;
    let mut next = move || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (state >> 33) as u32
    };
    let unit: Vec<u32> = (0..7).map(|_| 400 + next() % 50).collect();
    let genome: Vec<u32> = (0..3000)
        .map(|g| if (1000..2000).contains(&g) && next() % 15 != 0 { unit[g % 7] } else { next() % 400 })
        .collect();
    let (mut lengths, mut reads) = (Vec::new(), Vec::new());
    for k in 0..30usize {
        let (start, len) = ((k * 83) % 2400, 300 + (k * 37) % 300);
        let read: Vec<Token> = (0..len)
            .map(|i| {
                let value = if next() % 100 < noise { next() % 450 } else { genome[start + i] };
                Token::new(value, 10 * i as i64, 10 * i as i64 + 10)
            })
            .collect();
        lengths.push(10 * len as i64);
        reads.push(read);
    }
    (lengths, reads)
}

/// Reads cut from one random token genome, with a fraction of tokens replaced.
pub(super) fn noisy_reads(noise: u32, seed: u64) -> (Vec<i64>, Vec<Vec<Token>>) {
    let mut state = seed;
    let mut next = move || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (state >> 33) as u32
    };
    let genome: Vec<u32> = (0..3000).map(|_| next() % 400).collect();
    let (mut lengths, mut reads) = (Vec::new(), Vec::new());
    for k in 0..30usize {
        let (start, len) = ((k * 83) % 2400, 300 + (k * 37) % 300);
        let read: Vec<Token> = (0..len)
            .map(|i| {
                let value = if next() % 100 < noise { next() % 400 } else { genome[start + i] };
                Token::new(value, 10 * i as i64, 10 * i as i64 + 10)
            })
            .collect();
        lengths.push(10 * len as i64);
        reads.push(read);
    }
    (lengths, reads)
}
