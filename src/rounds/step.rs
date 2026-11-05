//! The `rounds rounds` step: every correction round in one process.

use super::*;

/// Corrected tokens and the number of changed layers of every read of one round, in catalog order.
fn replay_all(
    tokens: &[Vec<Token>],
    lengths: &[i64],
    symbols: &[i64],
    grouped: &[Vec<DonorRow>],
    pad: usize,
    band: usize,
    pool: &rayon::ThreadPool,
) -> Vec<(Vec<Token>, usize)> {
    pool.install(|| {
        (0..tokens.len())
            .into_par_iter()
            .map(|r| {
                let replay = replay_read(&tokens[r], lengths[r], symbols[r], &grouped[r], tokens, lengths, symbols, pad, band);
                (spanned_tokens(&tokens[r], &replay, lengths[r]), replay.changed_layers())
            })
            .collect()
    })
}

pub fn run() -> Result<()> {
    crate::partners::set_mmap_threshold();   // large blocks are returned to the system when freed
    if crate::argv::args().nth(1).as_deref() == Some("rounds") {
        return run_rounds();
    }
    Err("usage: squash __rounds rounds OPTIONS (the correction rounds)".into())
}

/// Whether two reads have the same token values (spans ignored), as the fixed-point rule compares them.
fn same_values(a: &[Token], b: &[Token]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.value == y.value)
}

/// `SQUASH_ROUNDS_TIMING=1`: the start-up's and each round's parts on stderr (`[rounds-timing]` lines; a diagnostic).
fn rounds_timing_on() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("SQUASH_ROUNDS_TIMING").map_or(false, |v| v == "1"))
}

/// `rounds rounds`: every correction round in one process, with the tokens kept in memory, from round 1. Round 1
/// replays the prepared tokens with the first search's partners (`--selected`); round k >= 2 searches partners on round
/// k-1's corrected tokens, with the raw-evidence guard on the prepared tokens (a partner is kept only when the
/// uncorrected tokens still agree), and replays them. Each round writes its token table to
/// `--scratch-dir/roundK/tokens.tsv` (or through pigz, `--tokens-gzip-dir`) and its partner table (uncompressed) to
/// `--scratch-dir/partnersK.tsv`. Stop rules: a fixed point (the round reproduces its input; it is not to be
/// published), `--stop-below F` (converged when fewer than F of the input's token layers changed) and `--max-rounds`.
/// After each round's table is complete, one line goes to stdout: `ROUND<TAB>k<TAB>status<TAB>directory` with status
/// continue, last, converged or fixed_point; the pipeline publishes the rounds from these lines.
/// `--exclude-reads TABLE`: reads that are never searched targets or donors (the pipeline passes the
/// low-complexity reads; round 1's `--selected` already leaves them out).
/// The other options are described where they are parsed below.
fn run_rounds() -> Result<()> {
    use crate::partners::{kept_values, load_excluded, take_guard, take_params};
    use std::sync::Arc;
    use std::time::Instant;

    let mut map: HashMap<String, String> = HashMap::new();
    let mut argv = crate::argv::args().skip(2);
    while let Some(key) = argv.next() {
        if !key.starts_with("--") {
            return Err(format!("unexpected argument {:?}", key).into());
        }
        let value = argv.next().ok_or_else(|| format!("{} needs a value", key))?;
        map.insert(key[2..].to_string(), value);
    }
    let p = take_params(&mut map)?;
    let g = take_guard(&mut map)?;
    // The q-symbol alignment check in every partner search of the rounds, with the first search's
    // verdicts (`--signal-verdicts`) reused; pairs new in a round are checked then.
    let sig = crate::partners::signal::take(&mut map)?;
    fn need(map: &mut HashMap<String, String>, key: &str) -> Result<String> {
        map.remove(key).ok_or_else(|| format!("rounds: missing --{}", key).into())
    }
    fn number<T: std::str::FromStr>(map: &mut HashMap<String, String>, key: &str, default: T) -> Result<T>
    where
        T::Err: std::fmt::Display,
    {
        match map.remove(key) {
            None => Ok(default),
            Some(v) => v.parse::<T>().map_err(|e| format!("--{} {:?}: {}", key, v, e).into()),
        }
    }
    let reads_dir = PathBuf::from(need(&mut map, "reads-dir")?);
    let symbols_path = PathBuf::from(need(&mut map, "symbols")?);
    let selected_path = PathBuf::from(need(&mut map, "selected")?);
    let scratch = PathBuf::from(need(&mut map, "scratch-dir")?);
    // Each donor aligned over its overlap widened by `pad` tokens, with free end gaps within 2 x `pad` tokens of the
    // window's ends (`replay_read`).
    let pad: usize = number(&mut map, "pad", 0)?;
    // Banded alignment of half-width `band` tokens (align_overlap_banded); 0 = the full DP.
    let band: usize = number(&mut map, "band", 0)?;
    let threads: usize = number(&mut map, "threads", 1)?;
    let max_rounds: usize = number(&mut map, "max-rounds", 12)?;
    let stop_below: Option<f64> = map.remove("stop-below").map(|x| x.parse()).transpose()?;
    // Reads that are neither searched targets nor donors.
    let exclude_path = map.remove("exclude-reads").map(PathBuf::from);
    // Segment skipping: from round `skip-from` (0 = off) a round replays only each read's dirty windows
    // (`dirty_windows`, padded by `dirty-pad` q symbols) with 2 x `pad` layers around them. Needs the previous round,
    // so the first round runs in full.
    let skip_from: usize = number(&mut map, "skip-from", 0)?;
    let dirty_pad: i64 = number(&mut map, "dirty-pad", 150)?;
    let skip_context = 2 * pad;
    // From round R (0 = off; R >= 2, round 2 taking round 1's table as the last table): `hops-from` chains only the
    // donors of the last table's donors (two hops: a target's donors and theirs); `changed-from` searches again only
    // the targets whose tokens changed in the last round, the others keep their last rows.
    let hops_from: usize = number(&mut map, "hops-from", 0)?;
    let changed_from: usize = number(&mut map, "changed-from", 0)?;
    // `tokens-gzip-dir DIR` = each round's table through pigz into DIR (TokensGzip), with `tokens-gzip-threads`
    // (default 1) and `tokens-gzip-level` (1-9, default 1), instead of the uncompressed `roundK/tokens.tsv`.
    let gzip_threads: Option<usize> = map.remove("tokens-gzip-threads").map(|x| x.parse()).transpose()?;
    let gzip_level: Option<u32> = map.remove("tokens-gzip-level").map(|x| x.parse()).transpose()?;
    let tokens_gzip = match map.remove("tokens-gzip-dir") {
        Some(dir) => Some(TokensGzip { dir: PathBuf::from(dir), threads: gzip_threads.unwrap_or(1).max(1), level: gzip_level.unwrap_or(1) }),
        None if gzip_threads.is_some() || gzip_level.is_some() => {
            return Err("rounds: --tokens-gzip-threads / --tokens-gzip-level need --tokens-gzip-dir".into());
        }
        None => None,
    };
    if tokens_gzip.as_ref().map_or(false, |g| !(1..=9).contains(&g.level)) {
        return Err("rounds: --tokens-gzip-level must be 1-9".into());
    }
    // `sparse-last INV:FRACTION` = the stopping round's (last or converged) sparse table made in memory into the
    // tokens-gzip directory (`write_sparse_last`; `sparse-gzip-threads` T, default 1), announced by
    // `SPARSE<TAB>k<TAB>directory` after the round's ROUND line.
    let sparse_threads: usize = number(&mut map, "sparse-gzip-threads", 1)?;
    let sparse_last: Option<(u64, f64, usize)> = match map.remove("sparse-last") {
        None => None,
        Some(spec) => {
            let (inv, fraction) = spec.split_once(':').ok_or("rounds: --sparse-last takes INV:FRACTION")?;
            let (inv, fraction): (u64, f64) = (inv.parse()?, fraction.parse()?);
            if inv < 1 || !(fraction > 0.0 && fraction < 1.0) {
                return Err("rounds: --sparse-last needs INV >= 1 and 0 < FRACTION < 1".into());
            }
            if tokens_gzip.is_none() {
                return Err("rounds: --sparse-last needs --tokens-gzip-dir".into());
            }
            Some((inv, fraction, sparse_threads.max(1)))
        }
    };
    if let Some(key) = map.keys().next() {
        return Err(format!("rounds: unknown option --{}", key).into());
    }
    if max_rounds < 1 {
        return Err("rounds: --max-rounds must be at least 1".into());
    }
    // 2 is allowed: round 2 then takes round 1's table (the input selection) as the last round's table.
    if hops_from == 1 || changed_from == 1 {
        return Err("rounds: --hops-from / --changed-from must be 0 (off) or >= 2".into());
    }
    if skip_from == 1 || dirty_pad < 0 {
        return Err("rounds: --skip-from must be 0 (off) or >= 2, --dirty-pad >= 0".into());
    }
    if !scratch.is_dir() {
        return Err(format!("rounds: missing scratch directory {}", scratch.display()).into());
    }

    let load = Instant::now();
    // `SQUASH_ROUNDS_TIMING=1`: the start-up's parts.
    let mut load_clock = Instant::now();
    let mut load_part = |name: &str| {
        if rounds_timing_on() {
            eprintln!("[rounds-timing] load: {} {:.2} s", name, load_clock.elapsed().as_secs_f64());
        }
        load_clock = Instant::now();
    };
    let (ids, lengths, prepared) = load_prepared(&reads_dir)?;
    load_part("prepared tokens");
    let (symbol_ids, symbols) = symbol_counts(&symbols_path)?;
    load_part("q-symbol counts");
    if symbol_ids != ids {
        return Err("source-q and prepared-token catalogs differ".into());
    }
    if ids.iter().collect::<HashSet<_>>().len() != ids.len() {
        return Err("duplicate read ids".into());
    }
    let excluded = match &exclude_path {
        None => Vec::new(),
        Some(path) => load_excluded(path, &ids)?,
    };
    let ids = Arc::new(ids);
    let lengths = Arc::new(lengths);
    let index: HashMap<&str, usize> = ids.iter().enumerate().map(|(i, x)| (x.as_str(), i)).collect();
    // The donor tables are shared (Arc) by the round, the two-hop candidates and segment skip, not cloned.
    let selected = Arc::new(load_donors(&selected_path, &index, ids.len())?.0);
    load_part("round 1 donors");
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads.max(1)).build()?;
    if let Some(o) = &sig {
        if o.verdicts_out.is_some() {
            return Err("rounds: the signal check takes no --signal-verdicts-out".into());
        }
        crate::partners::signal::install(o, &ids, &pool)?;
        // only the searches read the q-symbols: from a local file, pair by pair (`--signal-store` DIR, else the
        // scratch directory)
        crate::partners::signal::store_on_file(o.store.as_deref().unwrap_or(&scratch))?;
    }
    load_part("q-symbol alignment check");
    // the guard's kept values of the uncorrected tokens
    let guard_kept = kept_values(&prepared, p.max_reads, p.repeat_fraction);
    load_part("guard kept values");
    let prepared = Arc::new(prepared);
    let mut input = prepared.clone();
    // `--changed-from` needs only which reads changed in the last round, so that is all that is kept of it.
    let mut changed_flags: Option<Vec<bool>> = None;
    eprintln!("[rounds] {} reads loaded in {:.1} s", ids.len(), load.elapsed().as_secs_f64());
    {
        use crate::partners::mem;
        mem::report("loaded", &[("prepared_tokens", mem::nested(&prepared)), ("selected_rows", mem::nested(&selected)),
            ("guard_kept_values", guard_kept.capacity() * 8)]);
    }
    // The whole prepared reads are needed only while they are a round's input (round 1); the raw-evidence guard
    // keeps its own view (`guard_tokens`: kept values only, value and start) and the whole copy is dropped as soon as
    // no round reads it (`retire_prepared`).
    let mut prepared: Option<Arc<Vec<Vec<Token>>>> = Some(prepared);
    let mut guard_view: Option<Vec<crate::partners::GuardRead>> = None;
    // The guard's view, once built, is kept in a local file between searches (`spill_guard_view`) and read back before
    // each, instead of in memory through every step.
    let mut guard_spill: Option<PathBuf> = None;
    let retire_prepared = |prepared: &mut Option<Arc<Vec<Vec<Token>>>>, guard_view: &mut Option<Vec<crate::partners::GuardRead>>,
                           input: &Arc<Vec<Vec<Token>>>| {
        if let Some(full) = prepared.as_ref() {
            if guard_view.is_none() {
                *guard_view = Some(pool.install(|| crate::partners::guard_tokens(full, &guard_kept)));
            }
            if !Arc::ptr_eq(full, input) {
                *prepared = None;
            }
        }
    };

    let mut writer: Option<std::thread::JoinHandle<Result<()>>> = None;
    // Segment skipping: the previous round's rows and changed spans (kept only when skipping is on).
    let mut previous: Option<(Arc<Vec<Vec<DonorRow>>>, Vec<Vec<(i32, i32)>>)> = None;
    // `--hops-from` / `--changed-from`: the last round's table, as rows (packed, `PackedRows`) and as loaded donors.
    let mut prev_rows: Option<crate::partners::PackedRows> = None;
    let mut prev_grouped: Option<Arc<Vec<Vec<DonorRow>>>> = None;
    // From round 2: round 1's table, the input selection, is the last round's table of round 2: its donors for
    // --hops-from 2, its rows (each target's lines in file order, as the partner search wrote them) for --changed-from
    // 2.
    if hops_from == 2 || changed_from == 2 {
        if hops_from == 2 {
            prev_grouped = Some(Arc::clone(&selected));
        }
        if changed_from == 2 {
            // The rows are copied into round 2's table, so round 1's table must have the partner table's columns.
            let header = open_text(&selected_path)?.lines().next().ok_or("empty selected table")??;
            if header != crate::partners::COLUMNS {
                return Err("rounds: --changed-from 2 needs round 1's table in the partner table's columns".into());
            }
            // read a line at a time
            let c_t = header.split('\t').position(|c| c == "target_id").ok_or("missing column target_id")?;
            let mut rows: Vec<Vec<String>> = vec![Vec::new(); ids.len()];
            for line in open_text(&selected_path)?.lines().skip(1) {
                let line = line?;
                let target = line.split('\t').nth(c_t).unwrap_or("");
                let t = *index.get(target).ok_or("selected target is outside the token catalog")?;
                rows[t].push(line);
            }
            prev_rows = Some(crate::partners::PackedRows::pack(rows));
        }
    }
    for k in 1..=max_rounds {
        let clock = Instant::now();
        // `SQUASH_ROUNDS_TIMING=1`: the round's parts.
        let mut part_clock = Instant::now();
        let mut part = |name: &str| {
            if rounds_timing_on() {
                eprintln!("[rounds-timing] round {}: {} {:.2} s", k, name, part_clock.elapsed().as_secs_f64());
            }
            part_clock = Instant::now();
        };
        let grouped_rows: Arc<Vec<Vec<DonorRow>>> = if k > 1 {
            retire_prepared(&mut prepared, &mut guard_view, &input);
            part("retire_prepared");
            // the guard's view back from its local file for the search
            if guard_view.is_none() {
                if let Some(path) = &guard_spill {
                    guard_view = Some(crate::partners::load_guard_view(path)?);
                }
            }
            part("load_guard_view");
            let guard_arg = Some((guard_view.as_ref().expect("the guard's view").as_slice(), &guard_kept, &g));
            let path = scratch.join(format!("partners{}.tsv", k));
            let tmp = scratch.join(format!(".partners{}.tsv.tmp", k));
            let hops = hops_from > 0 && k >= hops_from && prev_grouped.is_some();
            let changed = changed_from > 0 && k >= changed_from && prev_rows.is_some() && changed_flags.is_some();
            let allowed: Option<Vec<Vec<u32>>> = prev_grouped.as_ref().filter(|_| hops).map(|g| two_hop_donors(g, &pool));
            part("two_hop_donors");
            let redo: Option<Vec<bool>> = if changed { changed_flags.clone() } else { None };
            if let Some(a) = &allowed {
                let n: usize = a.iter().map(Vec::len).sum();
                eprintln!("[rounds] round {}: two hops, {:.0} candidate donors per target", k, n as f64 / a.len().max(1) as f64);
            }
            if let Some(m) = &redo {
                eprintln!("[rounds] round {}: {} of {} targets changed, searched again", k, m.iter().filter(|&&x| x).count(), m.len());
            }
            let previous = redo.as_deref().map(|m| (m, prev_rows.as_ref().unwrap()));
            // each target's rows packed as found (`partner_rows_packed`), the table written from them
            let packed = crate::partners::partner_rows_packed(&ids, &lengths, &input, guard_arg, &p, &excluded, &pool,
                                                              allowed.as_deref(), previous);
            part("search");
            crate::partners::signal::report(&format!("round {}", k));
            crate::partners::write_packed_rows(&tmp, &packed)?;
            part("write_table");
            if hops_from > 0 || changed_from > 0 {
                // moved, not cloned (nothing else keeps them)
                prev_rows = Some(packed);
            }
            std::fs::rename(&tmp, &path)?;
            let searched = Arc::new(load_donors(&path, &index, ids.len())?.0);
            part("load_donors");
            if hops_from > 0 {
                prev_grouped = Some(Arc::clone(&searched));
            }
            // until the next search, the guard's view only in its local file
            if let Some(view) = guard_view.take() {
                if guard_spill.is_none() {
                    let file = scratch.join("guard_view.bin");
                    crate::partners::spill_guard_view(&file, &view)?;
                    guard_spill = Some(file);
                }
            }
            part("spill_guard_view");
            searched
        } else {
            Arc::clone(&selected)
        };
        let grouped: &[Vec<DonorRow>] = &grouped_rows;
        let search_s = clock.elapsed().as_secs_f64();
        // return the search's freed buffers to the system before the replay
        unsafe {
            libc::malloc_trim(0);
        }
        {
            use crate::partners::mem;
            mem::report(&format!("round {} partners searched", k), &[("input_tokens", if prepared.as_ref().map_or(false, |p| Arc::ptr_eq(&input, p)) { 0 } else { mem::nested(&input) }),
                ("prepared_tokens", prepared.as_ref().map_or(0, |p| mem::nested(p))),
                ("guard_view", guard_view.as_ref().map_or(0, |v| v.iter().map(|x| x.bytes()).sum())),
                ("donor_rows", mem::nested(grouped)), ("prev_rows", prev_rows.as_ref().map_or(0, |r| r.bytes())),
                ("prev_grouped", prev_grouped.as_ref().map_or(0, |g| mem::nested(g)))]);
        }
        let (corrected, stats): (Vec<Vec<Token>>, Vec<usize>) = if skip_from == 0 {
            replay_all(&input, &lengths, &symbols, grouped, pad, band, &pool).into_iter().unzip()
        } else {
            let skip = k >= skip_from && previous.is_some();
            let tracked: Vec<(Vec<Token>, usize, Vec<(i32, i32)>, usize)> = pool.install(|| {
                (0..input.len())
                    .into_par_iter()
                    .map(|r| {
                        // each read's changed spans packed here, as it is replayed
                        let (tokens, stats, changed, n) = if let (true, Some((rows, spans))) = (skip, &previous) {
                            let changes = RoundChanges { rows, spans };
                            let windows = dirty_windows(r, &grouped[r], &changes, &lengths, &symbols, dirty_pad);
                            replay_windows(r, &windows, &grouped[r], &input, &lengths, &symbols, pad, band, skip_context)
                        } else {
                            let replay = replay_read(&input[r], lengths[r], symbols[r], &grouped[r], &input, &lengths,
                                                     &symbols, pad, band);
                            (spanned_tokens(&input[r], &replay, lengths[r]), replay.changed_layers(),
                             changed_spans(&input[r], &replay.chosen), input[r].len())
                        };
                        (tokens, stats, pack_spans(changed), n)
                    })
                    .collect()
            });
            let replayed: usize = tracked.iter().map(|x| x.3).sum();
            let layers: usize = input.iter().map(Vec::len).sum();
            eprintln!("[rounds] round {}: {} replayed {} of {} layers ({:.1}%)", k, if skip { "segment skip" } else { "full" },
                      replayed, layers, 100.0 * replayed as f64 / layers.max(1) as f64);
            let mut spans = Vec::with_capacity(tracked.len());
            let mut out = (Vec::with_capacity(tracked.len()), Vec::with_capacity(tracked.len()));
            for (tokens, stats, changed, _) in tracked {
                out.0.push(tokens);
                out.1.push(stats);
                // the changed spans, merged and 32-bit (`pack_spans`): the same dirty windows as one span per
                // changed token
                spans.push(changed);
            }
            previous = Some((Arc::clone(&grouped_rows), spans));
            out
        };
        let replay_s = clock.elapsed().as_secs_f64() - search_s;
        part("replay");
        {
            use crate::partners::mem;
            mem::report(&format!("round {} replayed", k), &[("corrected_tokens", mem::nested(&corrected)),
                ("segment_skip_previous", previous.as_ref().map_or(0, |(r, s)| mem::nested(r) + mem::nested(s)))]);
        }
        release_dp_scratch(&pool);
        let corrected = Arc::new(corrected);
        let unchanged = (0..ids.len()).all(|r| same_values(&corrected[r], &input[r]));
        let changed: usize = stats.iter().sum();
        let total: usize = input.iter().map(Vec::len).sum();
        let status = if unchanged {
            "fixed_point"
        } else if stop_below.map_or(false, |f| (changed as f64) < f * total as f64) {
            "converged"
        } else if k == max_rounds {
            "last"
        } else {
            "continue"
        };
        if status != "continue" {
            // The stopping round: no partner search follows, so its state goes (the guard's view, segment skip's
            // previous rows and spans, the previous tables, the check's cached verdicts, and its input unless that is
            // the prepared reads), and the previous round's writer is waited for here.
            guard_view = None;
            previous = None;
            prev_rows = None;
            prev_grouped = None;
            changed_flags = None;
            crate::partners::signal::release_cache();
            if let Some(handle) = writer.take() {
                handle.join().map_err(|_| "round writer panicked")??;
            }
            if prepared.as_ref().map_or(true, |p| !Arc::ptr_eq(p, &input)) {
                input = Arc::new(Vec::new());
            }
            unsafe {
                libc::malloc_trim(0);
            }
            crate::partners::mem::report(&format!("round {} search state released", k), &[]);
        }
        eprintln!("[rounds] round {}: search {:.1} s, replay {:.1} s; {} ({} of {} layers changed)",
                  k, search_s, replay_s, status, changed, total);
        part("after replay (stats, status)");
        // One writer at a time keeps the ROUND lines in order; the next round computes meanwhile.
        if let Some(handle) = writer.take() {
            handle.join().map_err(|_| "round writer panicked")??;
        }
        part("wait for the previous round's writer");
        crate::partners::mem::report(&format!("round {}, previous round written", k), &[]);
        // The writer holds the corrected reads only (the next round's input), not this round's input.
        let (w_ids, w_lengths, w_corrected) = (ids.clone(), lengths.clone(), corrected.clone());
        let w_gzip = tokens_gzip.clone();
        // the folder of the round's table: the tokens-gzip folder, or `roundK` (written as `.roundK.tmp`, then renamed)
        let directory = match &tokens_gzip {
            Some(g) => g.dir.clone(),
            None => scratch.join(format!("round{}", k)),
        };
        let stage = scratch.join(format!(".round{}.tmp", k));
        writer = Some(std::thread::spawn(move || -> Result<()> {
            let write_clock = Instant::now();
            match &w_gzip {
                Some(g) => g.write(k, &w_ids, &w_lengths, &w_corrected)?,
                None => {
                    std::fs::create_dir(&stage)?;
                    write_token_table(&mut BufWriter::with_capacity(1 << 20, File::create(stage.join("tokens.tsv"))?),
                                      &w_ids, &w_lengths, &w_corrected)?;
                    std::fs::rename(&stage, &directory)?;
                }
            }
            if rounds_timing_on() {
                eprintln!("[rounds-timing] round {}: writer {:.2} s", k, write_clock.elapsed().as_secs_f64());
            }
            {
                let mut out = std::io::stdout().lock();
                writeln!(out, "ROUND\t{}\t{}\t{}", k, status, directory.display())?;
                out.flush()?;
            }
            // --sparse-last: the pipeline publishes the round meanwhile
            if let (Some((inv, fraction, threads)), Some(g), "last" | "converged") = (sparse_last, w_gzip.as_ref(), status) {
                write_sparse_last(&g.dir, k, &w_ids, &w_lengths, &w_corrected, inv, fraction, threads)?;
                let mut out = std::io::stdout().lock();
                writeln!(out, "SPARSE\t{}\t{}", k, g.dir.display())?;
                out.flush()?;
            }
            Ok(())
        }));
        if status != "continue" {
            break;
        }
        part("writer spawned");
        if changed_from > 0 {
            changed_flags = Some(pool.install(|| (0..input.len()).into_par_iter().map(|r| !same_values(&corrected[r], &input[r])).collect()));
        }
        part("changed flags");
        input = corrected;
        retire_prepared(&mut prepared, &mut guard_view, &input);
        part("retire prepared");
        // Return the round's freed memory to the system (glibc keeps it in its per-thread arenas).
        unsafe {
            libc::malloc_trim(0);
        }
        part("malloc_trim");
    }
    if let Some(handle) = writer.take() {
        handle.join().map_err(|_| "round writer panicked")??;
    }
    // the local copy of the guard's view (the q-symbol file is unlinked once open)
    let _ = std::fs::remove_file(scratch.join("guard_view.bin"));
    Ok(())
}

/// `--hops-from`: per target, the sorted donors of its donors and its own donors (not itself).
fn two_hop_donors(grouped: &[Vec<DonorRow>], pool: &rayon::ThreadPool) -> Vec<Vec<u32>> {
    pool.install(|| {
        (0..grouped.len())
            .into_par_iter()
            .map(|t| {
                let mut set: Vec<u32> = Vec::new();
                for row in &grouped[t] {
                    set.push(row.donor as u32);
                    set.extend(grouped[row.donor].iter().map(|x| x.donor as u32));
                }
                set.sort_unstable();
                set.dedup();
                set.retain(|&d| d as usize != t);
                set
            })
            .collect()
    })
}
