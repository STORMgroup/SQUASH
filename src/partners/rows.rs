//! The coverage cap and the partner table: each target's rows, packed rows, and the table's writers.

use super::*;

pub const COLUMNS: &str = "target_id\tdonor_id\trank\torientation\tanchor_count\tshared_seed_pairs\t\
target_start\ttarget_end\tdonor_start\tdonor_end\tchain_score\tanchor_density";


/// The coverage cap on one target's chains, then its rows (`COLUMNS`); see `target_rows_with`.
pub fn target_rows(r: usize, ids: &[String], tokens: &[Vec<Token>], chains: Vec<Chain>, coverage: usize) -> Vec<String> {
    target_rows_with(r, ids, tokens, chains, coverage, &|c: &Chain| deferred_verdict(r, c))
}

/// A deferred chain's q-symbol alignment check (`Deferred`), made now with the installed check.
fn deferred_verdict(r: usize, c: &Chain) -> bool {
    let (Some(dc), Some(sc)) = (c.check.as_deref(), signal::get()) else { return true };
    let prof = profile::on();
    let clock = std::time::Instant::now();
    let v = sc.verdict(r, c.donor, dc.line, dc.ts, dc.te, &dc.q);
    if prof {
        profile::since(18, clock);
        if v {
            profile::add(12, 1);
        }
    }
    v
}

/// `target_rows` with the verdict of a deferred chain from `verdict` (the lazy check, `LAZY_CHECK`). nanoMDBG's
/// coverage cap: the chains in decreasing score (then donor id), each kept while a token position of its span
/// (t0..=t1 on the target) is covered fewer than `coverage` times by the chains kept before it. A deferred chain that
/// reaches the cap with a position open is checked then and kept only if it passes; one with every position covered is
/// passed over unchecked. A failed chain adds no depth, as if it had been dropped before the cap, and a chain passed
/// over would have been rejected by the cap whatever its verdict: the rows of checking every chain before the cap
/// (`lazy_check_keeps_the_rows`), without the checks of the chains that cannot be kept. Chains without a deferred check
/// are never handed to `verdict`. Rows: rank from 1, orientation `+`, the anchor count (twice), the bounds on target
/// and donor (q-symbols), the chain score and anchors per spanned target token.
pub fn target_rows_with(r: usize, ids: &[String], tokens: &[Vec<Token>], chains: Vec<Chain>, coverage: usize,
                        verdict: &dyn Fn(&Chain) -> bool) -> Vec<String> {
    let mut chains = chains;
    chains.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap().then_with(|| ids[a.donor].cmp(&ids[b.donor])));
    let mut depth = vec![0usize; tokens[r].len()];
    let mut rows = Vec::new();
    for c in chains {
        if (c.t0..=c.t1).any(|i| depth[i] < coverage) {
            if c.check.is_some() {
                LAZY_CHECK[1].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if !verdict(&c) {
                    continue;
                }
            }
            for i in c.t0..=c.t1 {
                depth[i] += 1;
            }
            let rank = rows.len() + 1;
            rows.push(format!("{}\t{}\t{}\t+\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{:.4}", ids[r], ids[c.donor], rank, c.anchors,
                              c.anchors, c.ts, c.te, c.ds, c.de, c.score,
                              c.anchors as f64 / ((c.t1 - c.t0 + 1) as f64)));
        } else if c.check.is_some() {
            LAZY_CHECK[2].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
    }
    rows
}


#[allow(clippy::too_many_arguments)]
/// A round's partner rows (text lines per target) kept compressed for the next round, whose unchanged targets copy them
/// back as they were: per target its lines joined by '\n', zstd level 1 (as strings they grow round after round).
pub struct PackedRows {
    pub(super) blobs: Vec<Vec<u8>>,
}

/// One target's lines, packed (`PackedRows`).
pub(super) fn pack_lines(lines: Vec<String>) -> Vec<u8> {
    if lines.is_empty() {
        return Vec::new();
    }
    let mut out = zstd::bulk::compress(lines.join("\n").as_bytes(), 1).expect("compressing partner rows");
    out.shrink_to_fit();
    out
}

impl PackedRows {
    pub fn pack(rows: Vec<Vec<String>>) -> PackedRows {
        PackedRows { blobs: rows.into_par_iter().map(pack_lines).collect() }
    }

    /// Every target's lines.
    #[cfg(test)]
    pub fn unpack(&self) -> Vec<Vec<String>> {
        (0..self.blobs.len()).into_par_iter().map(|r| self.get(r)).collect()
    }

    /// Target `r`'s lines, as packed.
    pub fn get(&self, r: usize) -> Vec<String> {
        let blob = &self.blobs[r];
        if blob.is_empty() {
            return Vec::new();
        }
        let text = String::from_utf8(zstd::stream::decode_all(&blob[..]).expect("decompressing partner rows"))
            .expect("partner rows are UTF-8");
        text.split('\n').map(str::to_owned).collect()
    }

    /// Heap bytes (for the memory report).
    pub fn bytes(&self) -> usize {
        self.blobs.iter().map(|b| b.capacity()).sum::<usize>() + self.blobs.len() * std::mem::size_of::<Vec<u8>>()
    }
}


/// `write_rows` of packed rows: the same bytes.
pub fn write_packed_rows(path: &Path, rows: &PackedRows) -> Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "{}", COLUMNS)?;
    for r in 0..rows.blobs.len() {
        for row in rows.get(r) {
            writeln!(out, "{}", row)?;
        }
    }
    out.flush()?;
    Ok(())
}


/// The selected-donor table, uncompressed: header and rows.
#[cfg(test)]
pub fn write_rows(path: &Path, rows: &[Vec<String>]) -> Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "{}", COLUMNS)?;
    for target in rows {
        for row in target {
            writeln!(out, "{}", row)?;
        }
    }
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The cap with deferred checks (`target_rows_with`) keeps the rows of the check at chaining (chains failing the
    /// check dropped before the cap), for verdicts of every pattern, on repeat and noisy reads at several coverages,
    /// and passes some chains over unchecked (else the test would not see the saving); with no deferred chain it never
    /// asks.
    #[test]
    fn lazy_check_keeps_the_rows() {
        let (mut asked, mut chains_seen, mut compared) = (0usize, 0usize, 0usize);
        for (seed, noise, repeats) in [(1, 10, true), (2, 30, true), (5, 5, true), (6, 5, false), (7, 20, false)] {
            let (lengths, tokens) = if repeats { repeat_reads(noise, seed) } else { noisy_reads(noise, seed) };
            let n = tokens.len();
            let ids: Vec<String> = (0..n).map(|r| format!("read{:02}", r)).collect();
            for coverage in [1usize, 2, 3, 20] {
                let p = Params { max_reads: 0, min_anchors: 3, min_overlap: 300, min_cover: 0.1, coverage, ..Params::default() };
                let index = build_index(&tokens, 0, 0.0);
                let mut s = Scratch::new(n);
                for r in 0..n {
                    let chains = target_chains(r, &lengths, &tokens, &index, &p, None, &mut s, None);
                    for pattern in 0..4u64 {
                        let pass = |d: usize| -> bool {
                            let h = (r as u64 * 1_000_003 + d as u64 * 7_919 + pattern * 104_729) % 7;
                            match pattern { 0 => true, 1 => h != 0, 2 => h > 3, _ => h == 6 }
                        };
                        let eager: Vec<Chain> = chains.iter().filter(|c| pass(c.donor)).cloned().collect();
                        let want = target_rows(r, &ids, &tokens, eager, coverage);
                        let deferred: Vec<Chain> = chains.iter().cloned().map(|mut c| {
                            c.check = Some(Box::new(Deferred { line: (0, 0, 0, 0), ts: 0, te: 0, q: Vec::new() }));
                            c
                        }).collect();
                        let calls = std::cell::Cell::new(0usize);
                        let got = target_rows_with(r, &ids, &tokens, deferred, coverage, &|c: &Chain| {
                            calls.set(calls.get() + 1);
                            pass(c.donor)
                        });
                        assert_eq!(got, want, "seed {} target {} coverage {} pattern {}", seed, r, coverage, pattern);
                        let plain = target_rows_with(r, &ids, &tokens, chains.clone(), coverage, &|_: &Chain| panic!("asked"));
                        assert_eq!(plain, target_rows(r, &ids, &tokens, chains.clone(), coverage));
                        asked += calls.get();
                        chains_seen += chains.len();
                        compared += 1;
                    }
                }
            }
        }
        assert!(compared > 500, "{} cases", compared);
        assert!(asked < chains_seen, "every chain was checked ({} of {}): the test would not see the saving", asked, chains_seen);
    }
}
