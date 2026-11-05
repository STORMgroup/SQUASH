//! The search of every target (`partner_rows_packed`): each target on the pool (heavy targets apart), its rows
//! packed as found.

use super::*;

#[cfg(test)]
/// Every target's rows (`target_chains`, then `target_rows`), in catalog order. `guard`: the uncorrected tokens, their
/// `kept_values` and the rule. `excluded` (by ordinal; empty: none): reads that are neither corrected (no rows) nor
/// donors (the pipeline's low-complexity reads, `--exclude-reads`).
pub fn partner_rows(
    ids: &[String],
    lengths: &[i64],
    tokens: &[Vec<Token>],
    guard: Option<(&[GuardRead], &HashSet<u32>, &Guard)>,
    p: &Params,
    excluded: &[bool],
    pool: &rayon::ThreadPool,
) -> Vec<Vec<String>> {
    partner_rows_ext(ids, lengths, tokens, guard, p, excluded, pool, None, None)
}


/// While a partner search runs, freed memory is returned to the system every `TRIM_SECS` seconds by a helper thread
/// calling `malloc_trim`: glibc otherwise keeps what the threads free in their arenas, and the search's resident memory
/// creeps up.
const TRIM_SECS: u64 = 30;

struct PeriodicTrim {
    stop: std::sync::Arc<std::sync::atomic::AtomicBool>,
    handle: Option<std::thread::JoinHandle<()>>,
}

impl PeriodicTrim {
    fn start() -> PeriodicTrim {
        let secs = TRIM_SECS;
        let stop = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
        let flag = std::sync::Arc::clone(&stop);
        let handle = std::thread::spawn(move || {
            let mut waited = 0u64;
            while !flag.load(std::sync::atomic::Ordering::Relaxed) {
                std::thread::sleep(std::time::Duration::from_millis(250));
                waited += 250;
                if waited >= secs * 1000 {
                    waited = 0;
                    unsafe {
                        libc::malloc_trim(0);
                    }
                }
            }
        });
        PeriodicTrim { stop, handle: Some(handle) }
    }
}

impl Drop for PeriodicTrim {
    fn drop(&mut self) {
        self.stop.store(true, std::sync::atomic::Ordering::Relaxed);
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

thread_local! {
    /// Each pool thread's `Scratch` (five arrays over all reads), reused by every piece of work the thread takes
    /// instead of made anew for each piece rayon splits off; freed when the search ends. A target leaves the scratch as
    /// it found it, so reuse does not change the chains.
    static SEARCH_SCRATCH: std::cell::RefCell<Option<Scratch>> = std::cell::RefCell::new(None);
    /// `target_chains_free`'s buffers (`FreeScratch`), freed when the search ends.
    static FREE_SCRATCH: std::cell::RefCell<FreeScratch> = std::cell::RefCell::new(FreeScratch::default());
}

#[cfg(test)]
/// `partner_rows` with two restrictions (the correction rounds' `--hops-from` / `--changed-from`): `allowed` (per
/// target, sorted ordinals): only these donors are chained; `previous` (a flag per target and the last table's rows):
/// the targets not flagged keep their previous rows and are not searched.
pub fn partner_rows_ext(
    ids: &[String],
    lengths: &[i64],
    tokens: &[Vec<Token>],
    guard: Option<(&[GuardRead], &HashSet<u32>, &Guard)>,
    p: &Params,
    excluded: &[bool],
    pool: &rayon::ThreadPool,
    allowed: Option<&[Vec<u32>]>,
    previous: Option<(&[bool], &PackedRows)>,
) -> Vec<Vec<String>> {
    let _trim = PeriodicTrim::start();   // stops when the search returns
    let keep = |r: usize| previous.map_or_else(Vec::new, |(_, rows)| rows.get(r));
    search_targets(ids, lengths, tokens, guard, p, excluded, pool, allowed, previous.map(|(redo, _)| redo), &keep, &|lines| lines)
}

/// `partner_rows_ext` with every target's rows packed as soon as they are found (and an unchanged target's packed rows
/// kept as they were), so the round's whole table is never held as text.
pub fn partner_rows_packed(
    ids: &[String],
    lengths: &[i64],
    tokens: &[Vec<Token>],
    guard: Option<(&[GuardRead], &HashSet<u32>, &Guard)>,
    p: &Params,
    excluded: &[bool],
    pool: &rayon::ThreadPool,
    allowed: Option<&[Vec<u32>]>,
    previous: Option<(&[bool], &PackedRows)>,
) -> PackedRows {
    let _trim = PeriodicTrim::start();
    let keep = |r: usize| previous.map_or_else(Vec::new, |(_, rows)| rows.blobs[r].clone());
    let blobs = search_targets(ids, lengths, tokens, guard, p, excluded, pool, allowed, previous.map(|(redo, _)| redo),
                               &keep, &pack_lines);
    let packed = PackedRows { blobs };
    mem::report("partner rows", &[("rows_packed", packed.bytes())]);
    packed
}

/// The search of `partner_rows_ext` / `partner_rows_packed`: each target's rows through `finish`, or `keep(r)` for a
/// target `redo` does not flag; returns them in target order. With `allowed` donors (per target) the search goes
/// without the index (`target_chains_free`).
#[allow(clippy::too_many_arguments)]
fn search_targets<T: Send>(
    ids: &[String],
    lengths: &[i64],
    tokens: &[Vec<Token>],
    guard: Option<(&[GuardRead], &HashSet<u32>, &Guard)>,
    p: &Params,
    excluded: &[bool],
    pool: &rayon::ThreadPool,
    allowed: Option<&[Vec<u32>]>,
    redo: Option<&[bool]>,
    keep: &(dyn Fn(usize) -> T + Sync),
    finish: &(dyn Fn(Vec<String>) -> T + Sync),
) -> Vec<T> {
    let heavy_cfg = heavy::config();   // read on the caller's thread (tests set it there)
    pool.install(|| {
        if let Some(allowed) = allowed {
            let clock = std::time::Instant::now();
            let kept = kept_seed_values(tokens, p);
            eprintln!("[partners] no index: {} kept values in {:.1} s", kept.len(), clock.elapsed().as_secs_f64());
            mem::report("partner kept values (no index)", &[("kept_values", kept.capacity() * 5)]);
            let skip = |r: usize| !excluded.is_empty() && excluded[r];
            let rows: Vec<T> = (0..ids.len())
                .into_par_iter()
                .map(|r| {
                    if skip(r) {
                        return finish(Vec::new());
                    }
                    if let Some(redo) = redo {
                        if !redo[r] {
                            return keep(r);
                        }
                    }
                    let chains = FREE_SCRATCH.with(|cell| {
                        let mut scratch = cell.borrow_mut();
                        let chains = target_chains_free(r, lengths, tokens, &kept, excluded, p, guard, &mut scratch,
                                                        Some(p.coverage), &allowed[r]);
                        scratch.trim();
                        chains
                    });
                    let clock = std::time::Instant::now();
                    let rows = target_rows(r, ids, tokens, chains, p.coverage);
                    if profile::on() {
                        profile::since(6, clock);
                        profile::add(13, rows.len() as u64);
                    }
                    finish(rows)
                })
                .collect();
            if profile::on() {
                profile::report("search");
            }
            cover_probe::report("search");
            report_lazy_check();
            rayon::broadcast(|_| {
                CHAIN_SCRATCH.with(|cell| *cell.borrow_mut() = Default::default());
                FREE_SCRATCH.with(|cell| *cell.borrow_mut() = FreeScratch::default());
                signal::release_thread_buffers();
            });
            return rows;
        }
        let index = build_index_excluding(tokens, p.max_reads, p.repeat_fraction, excluded);
        let skip = |r: usize| !excluded.is_empty() && excluded[r];
        let body = |r: usize| -> T {
            if skip(r) {
                return finish(Vec::new());
            }
            if let Some(redo) = redo {
                if !redo[r] {
                    return keep(r);
                }
            }
            let chains = SEARCH_SCRATCH.with(|cell| {
                let mut cell = cell.borrow_mut();
                if cell.as_ref().map_or(true, |s| s.count.len() != ids.len()) {
                    *cell = Some(Scratch::new(ids.len()));
                    memtrack::set(memtrack::FIXED, ids.len() * 4 * 2);
                }
                let scratch = cell.as_mut().expect("the thread's scratch");
                let chains = target_chains(r, lengths, tokens, &index, p, guard, scratch, Some(p.coverage));
                scratch.trim();
                heavy::release();   // a heavy thread's bytes, once its buffers are trimmed (no-op elsewhere)
                chains
            });
            let clock = std::time::Instant::now();
            let rows = target_rows(r, ids, tokens, chains, p.coverage);
            if profile::on() {
                profile::since(6, clock);
                profile::add(13, rows.len() as u64);
            }
            finish(rows)
        };
        // Heavy targets apart (`heavy`) when every target is searched over every donor (no previous rows).
        let rows: Vec<T> = match heavy_cfg {
            Some(cfg) if redo.is_none() => {
                let estimate = |r: usize| -> u64 {
                    if skip(r) { 0 } else { tokens[r].iter().map(|t| index.range(t.value).len() as u64).sum() }
                };
                heavy::run(ids.len(), cfg, &estimate, &body)
            }
            _ => (0..ids.len()).into_par_iter().map(&body).collect(),
        };
        if profile::on() {
            profile::report("search");
        }
        cover_probe::report("search");
        let saturated = SATURATED.swap(0, std::sync::atomic::Ordering::Relaxed);
        if saturated > 0 {
            eprintln!("[partners] {saturated} target-donor pairs dropped: the donor's anchor count over u32");
        }
        report_lazy_check();
        // Each pool thread frees its chaining scratch (a caller's replay does not use it), its search scratch and the
        // signal check's alignment buffers.
        rayon::broadcast(|_| {
            CHAIN_SCRATCH.with(|cell| *cell.borrow_mut() = Default::default());
            SEARCH_SCRATCH.with(|cell| *cell.borrow_mut() = None);
            signal::release_thread_buffers();
            for c in [memtrack::ANCHORS, memtrack::CHAIN, memtrack::ORDER, memtrack::PATH_KEPT, memtrack::FIXED] {
                memtrack::set(c, 0);
            }
        });
        rows
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn allowed_donors_and_kept_rows_restrict_the_search() {
        let (lengths, tokens) = noisy_reads(10, 4);
        let ids: Vec<String> = (0..tokens.len()).map(|r| format!("read{:02}", r)).collect();
        let p = Params { max_reads: 0, min_anchors: 3, min_overlap: 500, min_cover: 0.1, ..Params::default() };
        let pool = rayon::ThreadPoolBuilder::new().num_threads(2).build().unwrap();
        let rows = partner_rows(&ids, &lengths, &tokens, None, &p, &[], &pool);
        // Every donor allowed: the same rows; no donor allowed: none.
        let every: Vec<Vec<u32>> = (0..tokens.len()).map(|t| (0..tokens.len() as u32).filter(|&d| d as usize != t).collect()).collect();
        assert_eq!(partner_rows_ext(&ids, &lengths, &tokens, None, &p, &[], &pool, Some(&every), None), rows);
        let none: Vec<Vec<u32>> = vec![Vec::new(); tokens.len()];
        assert!(partner_rows_ext(&ids, &lengths, &tokens, None, &p, &[], &pool, Some(&none), None).iter().all(Vec::is_empty));
        // Targets not flagged keep the given rows.
        let marker: Vec<Vec<String>> = (0..tokens.len()).map(|t| vec![format!("kept{}", t)]).collect();
        let redo: Vec<bool> = (0..tokens.len()).map(|t| t % 2 == 0).collect();
        let packed = PackedRows::pack(marker.clone());
        let mixed = partner_rows_ext(&ids, &lengths, &tokens, None, &p, &[], &pool, None, Some((&redo, &packed)));
        // Packed rows come back as they were (several lines, one line, an empty line, none).
        let samples = vec![vec!["a\tb".to_owned(), "c".to_owned()], vec!["x".to_owned()], vec![String::new()], vec![]];
        let back = PackedRows::pack(samples.clone());
        assert_eq!((0..samples.len()).map(|i| back.get(i)).collect::<Vec<_>>(), samples);
        // Rows packed while searching = the plain search's rows packed, with and without previous rows, and the
        // table written from them has the same bytes.
        let plain = partner_rows_ext(&ids, &lengths, &tokens, None, &p, &[], &pool, None, Some((&redo, &packed)));
        let live = partner_rows_packed(&ids, &lengths, &tokens, None, &p, &[], &pool, None, Some((&redo, &packed)));
        assert_eq!(live.unpack(), plain);
        assert_eq!(partner_rows_packed(&ids, &lengths, &tokens, None, &p, &[], &pool, None, None).unpack(), rows);
        let dir = std::env::temp_dir().join(format!("packed_rows_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        write_rows(&dir.join("a.tsv"), &plain).unwrap();
        write_packed_rows(&dir.join("b.tsv"), &live).unwrap();
        assert_eq!(std::fs::read(dir.join("a.tsv")).unwrap(), std::fs::read(dir.join("b.tsv")).unwrap());
        std::fs::remove_dir_all(&dir).unwrap();
        for t in 0..tokens.len() {
            assert_eq!(mixed[t], if redo[t] { rows[t].clone() } else { marker[t].clone() });
        }
    }

    /// An excluded read gets no rows and is no read's donor; the other reads keep their rows unless they named an
    /// excluded donor, and an empty exclusion is the default.
    #[test]
    fn excluded_reads_are_neither_targets_nor_donors() {
        let (lengths, tokens) = noisy_reads(10, 5);
        let ids: Vec<String> = (0..tokens.len()).map(|r| format!("read{:02}", r)).collect();
        let p = Params { max_reads: 0, min_anchors: 3, min_overlap: 500, min_cover: 0.1, ..Params::default() };
        let pool = rayon::ThreadPoolBuilder::new().num_threads(2).build().unwrap();
        let all = partner_rows(&ids, &lengths, &tokens, None, &p, &[], &pool);
        assert_eq!(all, partner_rows(&ids, &lengths, &tokens, None, &p, &vec![false; ids.len()], &pool));
        let mut excluded = vec![false; ids.len()];
        for r in [3, 4, 17] {
            excluded[r] = true;
        }
        let some = partner_rows(&ids, &lengths, &tokens, None, &p, &excluded, &pool);
        let donor = |row: &String| row.split('\t').nth(1).unwrap().to_string();
        let named: HashSet<String> = [3, 4, 17].iter().map(|&r| ids[r].clone()).collect();
        assert!(all.iter().flatten().any(|row| named.contains(&donor(row))));
        for r in 0..ids.len() {
            if excluded[r] {
                assert!(some[r].is_empty());
                continue;
            }
            assert!(some[r].iter().all(|row| !named.contains(&donor(row))));
            if all[r].iter().all(|row| !named.contains(&donor(row))) {
                assert_eq!(some[r], all[r], "target {}", r);
            }
        }
        // The kept values are counted over every read: a value seen in one kept read and one excluded read stays.
        let t = |v: u32| Token::new(v, 0, 1);
        let reads = vec![vec![t(7), t(8)], vec![t(7), t(9)], vec![t(8), t(9)]];
        let index = build_index_excluding(&reads, 0, 0.0, &[false, true, false]);
        assert_eq!(index.get(7), &[(0, 0)]);
        assert_eq!(index.get(9), &[(2, 1)]);
        assert_eq!(index.get(8), &[(0, 1), (2, 0)]);
    }
}
