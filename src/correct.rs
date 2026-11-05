//! The correct stage: the first partner search, then the correction rounds (which also make the final round's sparse
//! assembly tables), for one work folder, each step run as a child process through `exec::Runner`. The tokens are
//! corrected in rounds against overlapping reads, as nanoMDBG corrects minimizers.
//!
//! Fixed arguments (not options of `squash`): the partner chains' `--max-reads 0 --max-gap 120 --max-drift 6
//! --drift-fraction 0.25 --lookback 60` (the first search and every partner search of the rounds), the guard's
//! `--guard-min-anchors 2 --guard-band 150 --guard-band-fraction 0.15` (the rounds), `rounds --pad 8`, and the
//! q-symbol check's `--signal-window 1000 --signal-band 100 --signal-path 200 --signal-path-shift 1
//! --signal-edit-gap 2` (large preset). Float arguments are written in shortest round-trip form (`py_float`).
//!
//! Files (in the work folder, `Layout`):
//! - the partner search reads `tokens/{reads.tsv.gz,tokens.u32le,token_spans.u64le}` (the uncorrected tokens), the
//!   low-complexity list `correction/low_complexity.tsv` and, large preset, `signal/qsymbols.tsv`; writes
//!   `correction/partners.tsv.gz` (read by the rounds, `--selected`, and by the contig joins, `--join-partners`) and,
//!   large preset, `correction/qsym_verdicts.tsv` (read by the rounds, `--signal-verdicts`).
//! - the rounds read those, the prepared tokens again (round 1's input and the guard) and the q-symbols (each read's
//!   length; large preset also the check); they write, for the final round N only,
//!   `correction/roundN/tokens.tsv.gz` (the corrected tokens; read by the signal layout) and
//!   `correction/roundN/sparse/tokens.tsv.gz` (1 in 8 tokens by value hash, then metaMDBG's
//!   100-ppm repetitive values dropped; read by the bp scale and the graph). An earlier round's table is removed as
//!   soon as the next round is known, so only round N is kept.
//!
//! Stop rule: the rounds stop after `max_rounds` (status `last`), after a round that changed fewer than `stop_below`
//! of its input token layers (`converged`; that round is the final one), or at a round that reproduces its input
//! (`fixed_point`; that round is not kept and the round before it is the final one). With `Gzip::Pigz` the engine
//! writes the sparse table itself for a `last` / `converged` round (`--sparse-last`); after a fixed point, and always
//! with `Gzip::Builtin`, it is made by `__sparse sparse-repeat` from the final round's table (the same text).
//! Resuming is not supported: a rerun starts from round 1 (with segment skipping, two-hop and
//! changed-only partner searches, a resumed process would replay its first round in full and give other tokens).
//!
//! Compression only changes .gz bytes, never the text: with `Gzip::Pigz` (used when `pigz` is on PATH) the engine
//! streams each round's table through `pigz -1` and writes the sparse table through `pigz -9`; `Gzip::Builtin` needs
//! nothing (the engine writes each round's table uncompressed into the scratch folder, where about two rounds' tables
//! are kept at a time, and SQUASH gzips the kept tables in process, in independent gzip members).

use std::fs::{self, File};
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use crate::exec::Runner;
use crate::io::{READS_TABLE, TOKEN_SPANS, TOKEN_VALUES, py_float};
use crate::pipeline::work;

/// The partner search's chain and cap options, in the order they are passed.
#[derive(Clone, Debug, PartialEq)]
pub struct PartnerSearch {
    pub max_reads: u64,
    pub min_anchors: u64,
    pub max_gap: u64,
    pub max_drift: u64,
    pub drift_fraction: f64,
    pub lookback: u64,
    pub min_overlap: u64,
    pub min_cover: f64,
    pub coverage: u64,
    pub repeat_fraction: f64,
}

/// The guard of the rounds' partner searches: a partner is kept only when the uncorrected tokens still agree.
#[derive(Clone, Debug, PartialEq)]
pub struct Guard {
    pub min_anchors: u64,
    pub band: u64,
    pub band_fraction: f64,
}

/// The q-symbol alignment check (large preset): `thr` per 1,000 alignment columns, the seed window (q-symbols) and
/// band, the path's band and shift, the gap cost of the weighted edit distance.
#[derive(Clone, Debug, PartialEq)]
pub struct SignalCheck {
    pub thr: f64,
    pub window: u64,
    pub band: u64,
    pub path: u64,
    pub path_shift: u64,
    pub edit_gap: u64,
    /// The first search's cover test while the check is on (`--signal-min-cover`); the rounds keep `min_cover`.
    pub first_min_cover: f64,
}

/// The correction rounds' options.
#[derive(Clone, Debug, PartialEq)]
pub struct Rounds {
    pub max_rounds: u64,
    pub pad: u64,
    pub band: u64,
    pub stop_below: f64,
    /// Segment skipping from this round (0 = off), with `dirty_pad` q-symbols around a changed input.
    pub skip_from: u64,
    pub dirty_pad: u64,
    /// Two-hop candidate donors / only changed targets searched, from this round (0 = off).
    pub hops_from: u64,
    pub changed_from: u64,
}

/// The sparse assembly selection of the final round: 1 in `inverse_density` tokens by value hash, then the
/// `repeat_ppm` per million most frequent values dropped.
#[derive(Clone, Debug, PartialEq)]
pub struct Sparse {
    pub inverse_density: u64,
    pub repeat_ppm: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Params {
    pub search: PartnerSearch,
    pub guard: Guard,
    /// `Some` = the q-symbol alignment check (large preset); the partner search always checks lazily, only the chains
    /// the coverage cap keeps.
    pub signal: Option<SignalCheck>,
    pub rounds: Rounds,
    pub sparse: Sparse,
}

impl Params {
    /// The small preset (genome < 10 Mb).
    pub fn small() -> Params {
        Params {
            search: PartnerSearch {
                max_reads: 0,
                min_anchors: 4,
                max_gap: 120,
                max_drift: 6,
                drift_fraction: 0.25,
                lookback: 60,
                min_overlap: 500,
                min_cover: 0.5,
                coverage: 20,
                repeat_fraction: 0.0001,
            },
            guard: Guard { min_anchors: 2, band: 150, band_fraction: 0.15 },
            signal: None,
            rounds: Rounds { max_rounds: 12, pad: 8, band: 256, stop_below: 0.0005, skip_from: 0, dirty_pad: 150, hops_from: 0, changed_from: 0 },
            sparse: Sparse { inverse_density: 8, repeat_ppm: 100 },
        }
    }

    /// The large preset (genome >= 10 Mb): the small preset plus the q-symbol alignment check (`thr` per 1,000
    /// columns), 16 rounds, two-hop partners and changed-only partner searches from round 2, segment skipping from
    /// round 8.
    pub fn large() -> Params {
        let mut p = Params::small();
        p.signal = Some(SignalCheck {
            thr: 667.0,
            window: 1000,
            band: 100,
            path: 200,
            path_shift: 1,
            edit_gap: 2,
            first_min_cover: 0.25,
        });
        p.rounds = Rounds { max_rounds: 16, skip_from: 8, dirty_pad: 150, hops_from: 2, changed_from: 2, ..p.rounds };
        p
    }

    /// The preset's values with the pipeline's options on top (`presets.rs` names; `get` = `Values::get`):
    /// `min-anchors`, `min-overlap`, `min-cover`, `coverage`, `seed-repeat-fraction` (the partner search, first search
    /// and rounds), `rounds`, `stop-below`, `band`, `hops-from`, `changed-from`, `skip-from` (`off` = 0), `skip-pad`
    /// (`--dirty-pad`, q-symbols), `qsym-check` (`off`, or the cost per alignment column: `--signal-thr` is per 1,000
    /// columns, 0.667 -> 667), `density` and `repeat-ppm` (the final round's sparse tables). Everything else is fixed.
    pub fn from_values(large: bool, get: &dyn Fn(&str) -> String) -> Result<Params, String> {
        let mut p = if large { Params::large() } else { Params::small() };
        let value = |name: &str| get(name).trim().to_string();
        let int = |name: &str| -> Result<u64, String> {
            value(name).parse::<u64>().map_err(|_| format!("--{}: '{}' is not a non-negative integer", name, value(name)))
        };
        let float = |name: &str| -> Result<f64, String> {
            value(name).parse::<f64>().ok().filter(|x| x.is_finite()).ok_or_else(|| format!("--{}: '{}' is not a number", name, value(name)))
        };
        let from = |name: &str| -> Result<u64, String> {
            match value(name).as_str() {
                "off" => Ok(0),
                _ => match int(name)? {
                    1 => Err(format!("--{}: 1 is not allowed (off or a round >= 2)", name)),
                    n => Ok(n),
                },
            }
        };
        p.search.min_anchors = int("min-anchors")?;
        p.search.min_overlap = int("min-overlap")?;
        p.search.min_cover = float("min-cover")?;
        p.search.coverage = int("coverage")?;
        p.search.repeat_fraction = float("seed-repeat-fraction")?;
        p.rounds.max_rounds = int("rounds")?;
        if p.rounds.max_rounds == 0 {
            return Err("--rounds: at least 1".into());
        }
        p.rounds.stop_below = float("stop-below")?;
        p.rounds.band = int("band")?;
        p.rounds.hops_from = from("hops-from")?;
        p.rounds.changed_from = from("changed-from")?;
        p.rounds.skip_from = from("skip-from")?;
        p.rounds.dirty_pad = int("skip-pad")?;
        p.sparse.inverse_density = int("density")?;
        p.sparse.repeat_ppm = int("repeat-ppm")?;
        if p.sparse.inverse_density == 0 || p.sparse.repeat_ppm == 0 || p.sparse.repeat_ppm >= 1_000_000 {
            return Err("--density must be >= 1 and --repeat-ppm in 1..999999".into());
        }
        p.signal = match value("qsym-check").as_str() {
            "off" => None,
            v => Some(SignalCheck { thr: per_mille(v)?, ..Params::large().signal.expect("the large preset's check") }),
        };
        Ok(p)
    }
}

/// A cost per alignment column as the check's threshold per 1,000 columns (`0.667` -> 667.0), rounded to 1e-6 so that
/// a decimal input gives its decimal value (`0.0041` -> 4.1, not 4.1000000000000005).
pub fn per_mille(v: &str) -> Result<f64, String> {
    let x: f64 = v.trim().parse().map_err(|_| format!("'{}' is not a number", v))?;
    if !(x.is_finite() && x >= 0.0) {
        return Err(format!("'{}' is not a non-negative number", v));
    }
    Ok((x * 1e9).round() / 1e6)
}

/// How the tables are gzipped (the text is the same either way).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Gzip {
    /// `pigz -n -p THREADS` (round tables streamed by the engine at level 1, the partner and sparse tables at level
    /// 9). Needs `pigz` on PATH.
    Pigz { threads: usize },
    /// No external tool: tables gzipped in process (`gzip_members`), the rounds' tables written uncompressed into the
    /// scratch folder first.
    Builtin { threads: usize },
}

/// The correct stage's files in the work folder (its inputs, and the outputs the later steps read).
#[derive(Clone, Debug)]
pub struct Layout {
    pub root: PathBuf,
    /// The low-complexity list, `--exclude-reads` of the search and the rounds; None = the rule is off (no list).
    pub low_complexity: Option<PathBuf>,
}

impl Layout {
    /// `root` with the presets' low-complexity list (`correction/low_complexity.tsv`).
    pub fn new(root: &Path) -> Layout {
        Layout::with_low_complexity(root, true)
    }
    /// `root` with the low-complexity list, or without one (the rule off).
    pub fn with_low_complexity(root: &Path, on: bool) -> Layout {
        Layout { root: root.to_path_buf(), low_complexity: on.then(|| root.join(work::LOW_COMPLEXITY)) }
    }
    /// The uncorrected tokens (a folder of token arrays).
    pub fn tokens_dir(&self) -> PathBuf {
        self.root.join(work::TOKENS)
    }
    pub fn symbols(&self) -> PathBuf {
        self.root.join(work::QSYMBOLS)
    }
    pub fn correction(&self) -> PathBuf {
        self.root.join("correction")
    }
    /// The first search's partner table (gzip).
    pub fn partners(&self) -> PathBuf {
        self.root.join(work::PARTNERS)
    }
    /// The partner search's uncompressed table, before it is gzipped.
    fn partners_plain(&self) -> PathBuf {
        self.correction().join("partners.tsv")
    }
    /// The first search's q-symbol check verdicts (large preset).
    pub fn verdicts(&self) -> PathBuf {
        self.correction().join("qsym_verdicts.tsv")
    }
    /// The gzipped tables of the rounds as the engine writes them, before the final round's are published.
    fn staging(&self) -> PathBuf {
        self.correction().join(".rounds_gzip")
    }
    pub fn round(&self, n: u64) -> PathBuf {
        self.correction().join(format!("round{}", n))
    }
    /// Round N's corrected tokens.
    pub fn tokens(&self, n: u64) -> PathBuf {
        self.round(n).join("tokens.tsv.gz")
    }
    /// Round N's sparse table's folder, the graph's input.
    pub fn sparse_dir(&self, n: u64) -> PathBuf {
        self.round(n).join("sparse")
    }
    /// The round folders (`roundN`) in the correction folder, as (N, path).
    fn rounds(&self) -> Result<Vec<(u64, PathBuf)>, String> {
        let dir = self.correction();
        let mut out = Vec::new();
        if dir.exists() {
            for e in fs::read_dir(&dir).map_err(io_err(&dir))? {
                let e = e.map_err(io_err(&dir))?;
                if let Some(n) = e.file_name().to_string_lossy().strip_prefix("round").and_then(|x| x.parse::<u64>().ok()) {
                    out.push((n, e.path()));
                }
            }
        }
        out.sort();
        Ok(out)
    }
}

/// What the rounds published.
#[derive(Clone, Debug, PartialEq)]
pub struct Corrected {
    /// The final round N (`correction/roundN`).
    pub round: u64,
    /// Why the rounds stopped: `last`, `converged` or `fixed_point` (then round N + 1 reproduced round N).
    pub status: String,
    pub tokens: PathBuf,
    pub sparse_rf: PathBuf,
}

/// Removes the correct stage's outputs from the work folder (for a rerun of the stage): the partner table (and its
/// uncompressed stage), the verdicts, the rounds' folders and their staging folder. The inputs (tokens/, signal/, the
/// low-complexity list) stay.
pub fn clear_outputs(l: &Layout) -> Result<(), String> {
    let plain = l.partners_plain();
    for f in [l.partners(), plain.with_extension("tmp"), plain, l.verdicts()] {
        if f.exists() {
            fs::remove_file(&f).map_err(io_err(&f))?;
        }
    }
    let staging = l.staging();
    let dirs = l.rounds()?.into_iter().map(|(_, d)| d).chain(staging.exists().then(|| staging));
    for d in dirs {
        fs::remove_dir_all(&d).map_err(io_err(&d))?;
    }
    Ok(())
}

/// The final round of a finished correction: the one round whose corrected tokens and sparse assembly table exist.
/// For the stages after this one and `--resume`.
pub fn final_round(l: &Layout) -> Result<u64, String> {
    let dir = l.correction();
    if !dir.is_dir() {
        return Err(format!("{}: no correction folder", dir.display()));
    }
    let found: Vec<u64> = l.rounds()?.into_iter().map(|(n, _)| n)
        .filter(|&n| l.tokens(n).is_file() && l.sparse_dir(n).join("tokens.tsv.gz").is_file()).collect();
    match found.as_slice() {
        [n] => Ok(*n),
        [] => Err(format!("{}: no finished correction round", dir.display())),
        _ => Err(format!("{}: several rounds with sparse tables: {:?}", dir.display(), found)),
    }
}

// ------------------------------------------------------------------------------------------------ arguments

fn path(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn push(args: &mut Vec<String>, key: &str, value: String) {
    args.push(format!("--{}", key));
    args.push(value);
}

/// The chain and cap options, shared by the first search and the rounds' partner searches.
fn search_options(args: &mut Vec<String>, r: &PartnerSearch) {
    push(args, "max-reads", r.max_reads.to_string());
    push(args, "min-anchors", r.min_anchors.to_string());
    push(args, "max-gap", r.max_gap.to_string());
    push(args, "max-drift", r.max_drift.to_string());
    push(args, "drift-fraction", py_float(r.drift_fraction));
    push(args, "lookback", r.lookback.to_string());
    push(args, "min-overlap", r.min_overlap.to_string());
    push(args, "min-cover", py_float(r.min_cover));
    push(args, "coverage", r.coverage.to_string());
    push(args, "repeat-fraction", py_float(r.repeat_fraction));
}

/// The q-symbol check's options; `first` = the first search's extras (`--signal-min-cover`, `--signal-verdicts-out`,
/// `--signal-store`), else the rounds' `--signal-verdicts`.
fn signal_options(args: &mut Vec<String>, s: &SignalCheck, symbols: &Path, first: Option<(&Path, &Path)>, verdicts_in: Option<&Path>) {
    push(args, "signal-symbols", path(symbols));
    push(args, "signal-thr", py_float(s.thr));
    push(args, "signal-window", s.window.to_string());
    if first.is_some() && s.first_min_cover >= 0.0 {
        push(args, "signal-min-cover", py_float(s.first_min_cover));
    }
    push(args, "signal-band", s.band.to_string());
    push(args, "signal-path", s.path.to_string());
    push(args, "signal-path-shift", s.path_shift.to_string());
    push(args, "signal-edit-gap", s.edit_gap.to_string());
    if let Some(v) = verdicts_in {
        push(args, "signal-verdicts", path(v));
    }
    if let Some((verdicts_out, store)) = first {
        push(args, "signal-verdicts-out", path(verdicts_out));
        push(args, "signal-store", path(store));
    }
}

/// The partner search's arguments: the step writes the uncompressed table to `output`; `store` is the q-symbol
/// check's node-local folder (large preset).
pub fn partner_args(l: &Layout, threads: usize, output: &Path, store: &Path, p: &Params) -> Vec<String> {
    let mut a = Vec::new();
    push(&mut a, "output", path(output));
    push(&mut a, "threads", threads.to_string());
    push(&mut a, "reads-dir", path(&l.tokens_dir()));
    search_options(&mut a, &p.search);
    if let Some(list) = &l.low_complexity {
        push(&mut a, "exclude-reads", path(list));
    }
    if let Some(s) = &p.signal {
        signal_options(&mut a, s, &l.symbols(), Some((&l.verdicts(), store)), None);
    }
    a
}

/// The correction rounds' arguments (run from round 1): `scratch` holds the engine's round files and searched
/// tables; with `Gzip::Pigz`, `staging` receives the gzipped round and sparse tables.
pub fn rounds_args(l: &Layout, threads: usize, scratch: &Path, staging: &Path, p: &Params, gzip: Gzip) -> Vec<String> {
    let r = &p.rounds;
    let mut a = vec!["rounds".to_string()];
    push(&mut a, "reads-dir", path(&l.tokens_dir()));
    push(&mut a, "symbols", path(&l.symbols()));
    push(&mut a, "selected", path(&l.partners()));
    push(&mut a, "scratch-dir", path(scratch));
    push(&mut a, "pad", r.pad.to_string());
    push(&mut a, "threads", threads.to_string());
    push(&mut a, "max-rounds", r.max_rounds.to_string());
    if r.band > 0 {
        push(&mut a, "band", r.band.to_string());
    }
    if r.skip_from > 0 {
        push(&mut a, "skip-from", r.skip_from.to_string());
        push(&mut a, "dirty-pad", r.dirty_pad.to_string());
    }
    if r.hops_from > 0 {
        push(&mut a, "hops-from", r.hops_from.to_string());
    }
    if r.changed_from > 0 {
        push(&mut a, "changed-from", r.changed_from.to_string());
    }
    if let Some(list) = &l.low_complexity {
        push(&mut a, "exclude-reads", path(list));
    }
    search_options(&mut a, &p.search);
    push(&mut a, "guard-min-anchors", p.guard.min_anchors.to_string());
    push(&mut a, "guard-band", p.guard.band.to_string());
    push(&mut a, "guard-band-fraction", py_float(p.guard.band_fraction));
    if let Some(s) = &p.signal {
        signal_options(&mut a, s, &l.symbols(), None, Some(&l.verdicts()));
    }
    if let Gzip::Pigz { threads: t } = gzip {
        push(&mut a, "tokens-gzip-dir", path(staging));
        push(&mut a, "tokens-gzip-threads", t.max(1).to_string());
        push(&mut a, "tokens-gzip-level", "1".to_string());
        push(&mut a, "sparse-last", sparse_spec(&p.sparse));
        push(&mut a, "sparse-gzip-threads", t.max(1).to_string());
    }
    push(&mut a, "stop-below", py_float(r.stop_below));
    a
}

/// `--sparse-last INV:FRACTION`: the density and the repeat fraction (ppm / 1e6, by `py_float`; 100 ppm -> `8:0.0001`).
pub fn sparse_spec(s: &Sparse) -> String {
    format!("{}:{}", s.inverse_density, py_float(s.repeat_ppm as f64 / 1e6))
}

/// The arguments of `__sparse sparse-repeat`, which makes the sparse table `out` (uncompressed) from a round's table.
pub fn select_args(tokens: &Path, s: &Sparse, out: &Path, pigz: usize) -> Vec<String> {
    let mut a = vec!["sparse-repeat".to_string()];
    push(&mut a, "tokens", path(tokens));
    push(&mut a, "inverse-density", s.inverse_density.to_string());
    push(&mut a, "fraction", py_float(s.repeat_ppm as f64 / 1e6));
    push(&mut a, "output", path(out));
    push(&mut a, "pigz", pigz.to_string());
    a
}

// ------------------------------------------------------------------------------------------------ steps

/// The correct stage for a work folder: the first partner search, then the rounds. `scratch` is a node-local folder
/// (the q-symbol store, the rounds' uncompressed files); it is left empty.
pub fn correct(runner: &Runner, l: &Layout, threads: usize, scratch: &Path, p: &Params, gzip: Gzip) -> Result<Corrected, String> {
    partner_search(runner, l, threads, scratch, p, gzip)?;
    correction_rounds(runner, l, threads, scratch, p, gzip)
}

fn io_err(what: &Path) -> impl Fn(io::Error) -> String + '_ {
    move |e| format!("{}: {}", what.display(), e)
}

/// A fresh empty folder of SQUASH's own (an earlier failed run's is removed).
fn fresh_dir(dir: &Path) -> Result<(), String> {
    if dir.exists() {
        fs::remove_dir_all(dir).map_err(io_err(dir))?;
    }
    fs::create_dir_all(dir).map_err(io_err(dir))
}

fn require(files: &[PathBuf]) -> Result<(), String> {
    match files.iter().find(|f| !f.is_file()) {
        Some(f) => Err(format!("missing input {}", f.display())),
        None => Ok(()),
    }
}

/// The first partner search (`__partners`; large preset with the q-symbol check, lazily). Writes the partner table
/// (gzip) and, large preset, the check's verdicts.
pub fn partner_search(runner: &Runner, l: &Layout, threads: usize, scratch: &Path, p: &Params, gzip: Gzip) -> Result<(), String> {
    let mut inputs: Vec<PathBuf> = [READS_TABLE, TOKEN_VALUES, TOKEN_SPANS].iter().map(|s| l.tokens_dir().join(s)).collect();
    inputs.extend(l.low_complexity.clone());
    if p.signal.is_some() {
        inputs.push(l.symbols());
    }
    require(&inputs)?;
    check_gzip(gzip)?;
    let (out, plain) = (l.partners(), l.partners_plain());
    for f in [&out, &l.verdicts()] {
        if f.exists() {
            return Err(format!("refusing to overwrite {}", f.display()));
        }
    }
    fs::create_dir_all(l.correction()).map_err(io_err(&l.correction()))?;
    for stale in [&plain, &plain.with_extension("tmp")] {
        if stale.exists() {
            fs::remove_file(stale).map_err(io_err(stale))?;
        }
    }
    let store = scratch.join("partners");
    fresh_dir(&store)?;
    runner.step("partners", "partners", &partner_args(l, threads, &plain, &store, p))?;
    gzip_file(&plain, &out, 9, gzip)?;
    fs::remove_file(&plain).map_err(io_err(&plain))?;
    fs::remove_dir_all(&store).map_err(io_err(&store))?;
    eprintln!("[correct] partner search -> {}", out.display());
    Ok(())
}

/// What to do with the rounds' files as the engine announces them (`ROUND k status dir`): which round's files can go
/// and, once it stops, which round is final. Round k's files are kept until round k + 1 is known, because a fixed
/// point at round k + 1 makes round k the final one.
#[derive(Debug, Default)]
pub struct RoundTracker {
    last: u64,
    pub final_round: Option<(u64, String)>,
}

impl RoundTracker {
    /// Round k's status: the rounds whose files are no longer needed.
    pub fn on_round(&mut self, k: u64, status: &str) -> Result<Vec<u64>, String> {
        if self.final_round.is_some() {
            return Err(format!("round {} announced after the final round", k));
        }
        if k != self.last + 1 {
            return Err(format!("round {} announced after round {}", k, self.last));
        }
        self.last = k;
        let earlier: Vec<u64> = if k > 1 { vec![k - 1] } else { vec![] };
        match status {
            "continue" => Ok(earlier),
            "last" | "converged" => {
                self.final_round = Some((k, status.to_string()));
                Ok(earlier)
            }
            "fixed_point" => {
                if k == 1 {
                    return Err("round 1 reproduces the uncorrected tokens (fixed point): no round to publish".into());
                }
                self.final_round = Some((k - 1, status.to_string()));
                Ok(vec![k])
            }
            other => Err(format!("round {}: unknown status {:?}", k, other)),
        }
    }
}

/// The correction rounds in one engine process (`__rounds rounds`), each round's files removed once the next round
/// is known; publishes the final round's tokens and sparse table. Needs the partner search's outputs.
pub fn correction_rounds(runner: &Runner, l: &Layout, threads: usize, scratch: &Path, p: &Params, gzip: Gzip) -> Result<Corrected, String> {
    let mut inputs: Vec<PathBuf> = [READS_TABLE, TOKEN_VALUES, TOKEN_SPANS].iter().map(|s| l.tokens_dir().join(s)).collect();
    inputs.extend(vec![l.symbols(), l.partners()]);
    inputs.extend(l.low_complexity.clone());
    if p.signal.is_some() {
        inputs.push(l.verdicts());
    }
    require(&inputs)?;
    check_gzip(gzip)?;
    if !l.rounds()?.is_empty() {
        return Err(format!("refusing to overwrite the rounds in {}", l.correction().display()));
    }
    let work = scratch.join("rounds");
    fresh_dir(&work)?;
    let staging = l.staging();   // beside the published rounds: published by rename
    fresh_dir(&staging)?;
    let pigz = matches!(gzip, Gzip::Pigz { .. });
    let table = |k: u64| if pigz { staging.join(format!("round{}.tokens.tsv.gz", k)) } else { work.join(format!("round{}/tokens.tsv", k)) };
    let mut tracker = RoundTracker::default();
    let mut sparse_from_engine: Option<u64> = None;
    let started = std::time::Instant::now();
    runner.step_lines("rounds", "rounds", &rounds_args(l, threads, &work, &staging, p, gzip), |line| {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "ROUND" if f.len() == 4 => {
                let k: u64 = f[1].parse().map_err(|_| format!("bad ROUND line {:?}", line))?;
                // The round's searched partners (round 2 on) are read back by the engine only within the round.
                let searched = work.join(format!("partners{}.tsv", k));
                if searched.exists() {
                    fs::remove_file(&searched).map_err(io_err(&searched))?;
                }
                eprintln!("[correct] round {}: {} ({:.1} s)", k, f[2], started.elapsed().as_secs_f64());
                for done in tracker.on_round(k, f[2])? {
                    remove_round(&work, done, &table(done))?;
                }
                Ok(())
            }
            "SPARSE" if f.len() == 3 => {
                sparse_from_engine = Some(f[1].parse().map_err(|_| format!("bad SPARSE line {:?}", line))?);
                Ok(())
            }
            _ => Ok(()),
        }
    })?;
    let (n, status) = tracker.final_round.clone().ok_or("the rounds stopped without a final round")?;
    // Publish round N: its tokens, then the sparse table (the engine's, or sparse's from the table).
    let tokens = l.tokens(n);
    let rf_dir = l.sparse_dir(n);
    for d in [tokens.parent().unwrap(), rf_dir.as_path()] {
        fs::create_dir_all(d).map_err(io_err(d))?;
    }
    let source = table(n);
    let sparse = match (sparse_from_engine, gzip) {
        (Some(k), Gzip::Pigz { .. }) if k == n => staging.join(format!("round{}.repeat.tsv.gz", n)),
        (Some(k), _) => return Err(format!("the engine wrote the sparse table of round {}, not of the final round {}", k, n)),
        (None, _) => {
            fresh_dir(&work.join("select"))?;
            let out = work.join("select/repeat.tsv");
            let pigz_threads = match gzip {
                Gzip::Pigz { threads } => threads.max(1),
                Gzip::Builtin { .. } => 0,
            };
            runner.step("sparse-select", "sparse", &select_args(&source, &p.sparse, &out, pigz_threads))?;
            out
        }
    };
    if pigz {
        rename(&source, &tokens)?;
    } else {
        gzip_file(&source, &tokens, 1, gzip)?;
    }
    remove_round(&work, n, &source)?;
    let sparse_rf = rf_dir.join("tokens.tsv.gz");
    if sparse.extension().map_or(false, |x| x == "gz") {
        rename(&sparse, &sparse_rf)?;
    } else {
        gzip_file(&sparse, &sparse_rf, 9, gzip)?;
        fs::remove_file(&sparse).map_err(io_err(&sparse))?;
    }
    // the scratch folders: the engine's leftovers (none but empty folders when it ends normally)
    fs::remove_dir_all(&work).map_err(io_err(&work))?;
    fs::remove_dir_all(&staging).map_err(io_err(&staging))?;
    eprintln!("[correct] rounds: round {} is final ({}) -> {}", n, status, l.round(n).display());
    Ok(Corrected { round: n, status, tokens, sparse_rf })
}

/// Round k's files in the engine's scratch folder and its table (`table`; in the staging folder with pigz).
fn remove_round(work: &Path, k: u64, table: &Path) -> Result<(), String> {
    if table.exists() {
        fs::remove_file(table).map_err(io_err(table))?;
    }
    let dir = work.join(format!("round{}", k));
    if dir.exists() {
        fs::remove_dir_all(&dir).map_err(io_err(&dir))?;
    }
    Ok(())
}

fn rename(from: &Path, to: &Path) -> Result<(), String> {
    fs::rename(from, to).map_err(|e| format!("{} -> {}: {}", from.display(), to.display(), e))
}

// ------------------------------------------------------------------------------------------------ gzip

fn check_gzip(gzip: Gzip) -> Result<(), String> {
    if let Gzip::Pigz { .. } = gzip {
        if std::process::Command::new("pigz").arg("--version").output().is_err() {
            return Err("pigz is not installed (needed by Gzip::Pigz)".into());
        }
    }
    Ok(())
}

/// `src` gzipped into `dst` (written beside it as `.NAME.part`, then renamed) at `level`: `pigz -c -LEVEL -n -p T`
/// or `gzip_members`.
pub fn gzip_file(src: &Path, dst: &Path, level: u32, gzip: Gzip) -> Result<(), String> {
    let name = dst.file_name().ok_or_else(|| format!("{}: no file name", dst.display()))?.to_string_lossy().into_owned();
    let part = dst.with_file_name(format!(".{}.part", name));
    match gzip {
        Gzip::Pigz { threads } => {
            let out = File::create(&part).map_err(io_err(&part))?;
            let status = std::process::Command::new("pigz")
                .args(["-c", &format!("-{}", level), "-n", "-p", &threads.max(1).to_string()])
                .arg(src)
                .stdout(out)
                .status()
                .map_err(|e| format!("pigz: {}", e))?;
            if !status.success() {
                return Err(format!("pigz {}: {}", src.display(), status));
            }
        }
        Gzip::Builtin { threads } => {
            let input = File::open(src).map_err(io_err(src))?;
            let out = File::create(&part).map_err(io_err(&part))?;
            gzip_members(input, out, level, threads, GZIP_BLOCK).map_err(io_err(src))?;
        }
    }
    rename(&part, dst)
}

/// The input block of `gzip_members`.
pub const GZIP_BLOCK: usize = 32 << 20;

/// `input` gzipped into `output` on `threads` threads: blocks of `block` bytes, each compressed as its own gzip member
/// (RFC 1952: a gzip file may hold several members; flate2's `MultiGzDecoder`, which every SQUASH step uses to read
/// these tables, Python's gzip, zcat and pigz read the whole text). An empty input gives one empty member.
pub fn gzip_members(mut input: impl Read, output: impl Write, level: u32, threads: usize, block: usize) -> io::Result<()> {
    use flate2::write::GzEncoder;
    use rayon::prelude::*;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads.max(1)).build()
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    let mut out = BufWriter::with_capacity(1 << 20, output);
    let batch = 2 * threads.max(1);
    let mut wrote = false;
    loop {
        let mut blocks: Vec<Vec<u8>> = Vec::with_capacity(batch);
        let mut end = false;
        while blocks.len() < batch {
            let mut buf = Vec::with_capacity(block);
            (&mut input).take(block as u64).read_to_end(&mut buf)?;
            if buf.len() < block {
                end = true;
            }
            if !buf.is_empty() {
                blocks.push(buf);
            }
            if end {
                break;
            }
        }
        if blocks.is_empty() && wrote {
            break;
        }
        if blocks.is_empty() {
            blocks.push(Vec::new());   // an empty input: one empty member
        }
        let members: Vec<io::Result<Vec<u8>>> = pool.install(|| {
            blocks.par_iter().map(|b| {
                let mut enc = GzEncoder::new(Vec::with_capacity(b.len() / 3 + 64), flate2::Compression::new(level));
                enc.write_all(b)?;
                enc.finish()
            }).collect()
        });
        for m in members {
            out.write_all(&m?)?;
        }
        wrote = true;
        if end {
            break;
        }
    }
    out.flush()
}

// ------------------------------------------------------------------------------------------------ development entry



// ------------------------------------------------------------------------------------------------ tests

#[cfg(test)]
mod tests {
    use super::*;

    fn s(v: &[&str]) -> Vec<String> {
        v.iter().map(|x| x.to_string()).collect()
    }

    #[test]
    fn py_float_writes_as_python() {
        // repr() of these floats in Python 3
        for (x, want) in [(0.25, "0.25"), (0.5, "0.5"), (0.0001, "0.0001"), (0.0005, "0.0005"), (0.15, "0.15"),
                          (667.0, "667.0"), (1.0, "1.0"), (100.0 / 1e6, "0.0001"), (1e-05, "1e-05"), (5e-05, "5e-05"),
                          (1.5e-07, "1.5e-07"), (1e15, "1000000000000000.0"), (1e16, "1e+16"), (1.25e17, "1.25e+17"),
                          (0.1, "0.1"), (2.0 / 3.0, "0.6666666666666666"), (123456.789, "123456.789"), (-0.5, "-0.5"),
                          (-3.0, "-3.0"), (0.0, "0.0"), (-0.0, "-0.0"), (1e-4 * 3.0, "0.00030000000000000003")] {
            assert_eq!(py_float(x), want, "{:?}", x);
        }
    }

    #[test]
    fn presets_differ_in_the_check_and_the_rounds() {
        let (small, large) = (Params::small(), Params::large());
        assert_eq!(small.search, large.search);
        assert_eq!(small.guard, large.guard);
        assert_eq!(small.sparse, large.sparse);
        assert!(small.signal.is_none() && large.signal.is_some());
        assert_eq!(small.rounds, Rounds { max_rounds: 12, pad: 8, band: 256, stop_below: 0.0005, skip_from: 0, dirty_pad: 150,
                                          hops_from: 0, changed_from: 0 });
        assert_eq!(large.rounds, Rounds { max_rounds: 16, pad: 8, band: 256, stop_below: 0.0005, skip_from: 8, dirty_pad: 150,
                                          hops_from: 2, changed_from: 2 });
        assert_eq!(sparse_spec(&small.sparse), "8:0.0001");
    }

    /// The partner search's full argument lists in both presets.
    #[test]
    fn partner_args_in_both_presets() {
        let l = Layout::new(Path::new("/R"));
        let base = s(&["--output", "/R/correction/partners.tsv", "--threads", "64",
                       "--reads-dir", "/R/tokens", "--max-reads", "0", "--min-anchors", "4", "--max-gap", "120",
                       "--max-drift", "6", "--drift-fraction", "0.25", "--lookback", "60", "--min-overlap", "500",
                       "--min-cover", "0.5", "--coverage", "20", "--repeat-fraction", "0.0001",
                       "--exclude-reads", "/R/correction/low_complexity.tsv"]);
        let out = l.partners_plain();
        assert_eq!(partner_args(&l, 64, &out, Path::new("/S/partners"), &Params::small()), base);
        let mut large = base.clone();
        large.extend(s(&["--signal-symbols", "/R/signal/qsymbols.tsv", "--signal-thr", "667.0",
                         "--signal-window", "1000", "--signal-min-cover", "0.25", "--signal-band", "100",
                         "--signal-path", "200", "--signal-path-shift", "1", "--signal-edit-gap", "2",
                         "--signal-verdicts-out", "/R/correction/qsym_verdicts.tsv", "--signal-store", "/S/partners"]));
        assert_eq!(partner_args(&l, 64, &out, Path::new("/S/partners"), &Params::large()), large);
    }

    /// The rounds' full argument lists in both presets, with pigz, and without the gzip options under builtin gzip.
    #[test]
    fn rounds_args_in_both_presets() {
        let l = Layout::new(Path::new("/R"));
        let head = s(&["rounds", "--reads-dir", "/R/tokens", "--symbols", "/R/signal/qsymbols.tsv",
                       "--selected", "/R/correction/partners.tsv.gz", "--scratch-dir", "/S/rounds",
                       "--pad", "8", "--threads", "64"]);
        let search = s(&["--exclude-reads", "/R/correction/low_complexity.tsv", "--max-reads", "0",
                          "--min-anchors", "4", "--max-gap", "120", "--max-drift", "6", "--drift-fraction", "0.25",
                          "--lookback", "60", "--min-overlap", "500", "--min-cover", "0.5", "--coverage", "20",
                          "--repeat-fraction", "0.0001", "--guard-min-anchors", "2", "--guard-band", "150",
                          "--guard-band-fraction", "0.15"]);
        let gz = s(&["--tokens-gzip-dir", "/R/correction/.rounds_gzip", "--tokens-gzip-threads", "64",
                     "--tokens-gzip-level", "1", "--sparse-last", "8:0.0001", "--sparse-gzip-threads", "64"]);
        let (work, staging) = (Path::new("/S/rounds"), Path::new("/R/correction/.rounds_gzip"));
        let mut small = head.clone();
        small.extend(s(&["--max-rounds", "12", "--band", "256"]));
        small.extend(search.clone());
        small.extend(gz.clone());
        small.extend(s(&["--stop-below", "0.0005"]));
        assert_eq!(rounds_args(&l, 64, work, staging, &Params::small(), Gzip::Pigz { threads: 64 }), small);
        let mut large = head.clone();
        large.extend(s(&["--max-rounds", "16", "--band", "256", "--skip-from", "8", "--dirty-pad", "150",
                         "--hops-from", "2", "--changed-from", "2"]));
        large.extend(search);
        large.extend(s(&["--signal-symbols", "/R/signal/qsymbols.tsv", "--signal-thr", "667.0",
                         "--signal-window", "1000", "--signal-band", "100", "--signal-path", "200", "--signal-path-shift", "1",
                         "--signal-edit-gap", "2", "--signal-verdicts", "/R/correction/qsym_verdicts.tsv"]));
        large.extend(gz);
        large.extend(s(&["--stop-below", "0.0005"]));
        assert_eq!(rounds_args(&l, 64, work, staging, &Params::large(), Gzip::Pigz { threads: 64 }), large);
        // Builtin: no gzip or sparse options (the engine writes plain tables; sparse makes the sparse ones)
        let builtin = rounds_args(&l, 64, work, staging, &Params::small(), Gzip::Builtin { threads: 64 });
        assert!(!builtin.iter().any(|a| a.starts_with("--tokens-gzip") || a.starts_with("--sparse")));
        assert_eq!(builtin.len(), small.len() - 10);
    }

    /// The pipeline's options table (presets.rs: each option's value in each preset) gives the presets.
    #[test]
    fn values_of_the_options_table_give_the_presets() {
        let table: &[(&str, &str, &str)] = &[
            ("min-anchors", "4", "4"), ("min-overlap", "500", "500"), ("min-cover", "0.5", "0.5"), ("coverage", "20", "20"),
            ("seed-repeat-fraction", "0.0001", "0.0001"), ("rounds", "12", "16"), ("stop-below", "0.0005", "0.0005"),
            ("band", "256", "256"), ("hops-from", "off", "2"), ("changed-from", "off", "2"), ("skip-from", "off", "8"),
            ("skip-pad", "150", "150"), ("qsym-check", "off", "0.667"), ("density", "8", "8"), ("repeat-ppm", "100", "100"),
        ];
        for large in [false, true] {
            let get = |name: &str| -> String {
                let row = table.iter().find(|r| r.0 == name).unwrap_or_else(|| panic!("asked for {}", name));
                (if large { row.2 } else { row.1 }).to_string()
            };
            let want = if large { Params::large() } else { Params::small() };
            assert_eq!(Params::from_values(large, &get).unwrap(), want);
            let l = Layout::new(Path::new("/R"));
            assert_eq!(partner_args(&l, 8, Path::new("/o"), Path::new("/s"), &Params::from_values(large, &get).unwrap()),
                       partner_args(&l, 8, Path::new("/o"), Path::new("/s"), &want));
        }
        // a user's value on top; the check's threshold per 1,000 columns, decimal-exact
        let get = |name: &str| -> String {
            match name {
                "qsym-check" => "0.581".into(),
                "rounds" => "20".into(),
                "hops-from" | "changed-from" | "skip-from" => "off".into(),
                _ => table.iter().find(|r| r.0 == name).unwrap().1.to_string(),
            }
        };
        let p = Params::from_values(false, &get).unwrap();
        assert_eq!(p.signal.as_ref().unwrap().thr, 581.0);
        assert_eq!(p.rounds.max_rounds, 20);
        assert_eq!(per_mille("0.667").unwrap(), 667.0);
        assert_eq!(per_mille("0.6675").unwrap(), 667.5);
        assert!(per_mille("x").is_err() && per_mille("-1").is_err());
        let bad = |key: &'static str, v: &'static str| {
            move |name: &str| -> String { if name == key { v.into() } else { table.iter().find(|r| r.0 == name).unwrap().2.to_string() } }
        };
        assert!(Params::from_values(true, &bad("hops-from", "1")).is_err());
        assert!(Params::from_values(true, &bad("rounds", "0")).is_err());
        assert!(Params::from_values(true, &bad("min-cover", "half")).is_err());
        assert!(Params::from_values(true, &bad("repeat-ppm", "0")).is_err());
    }

    /// The low-complexity rule off: no --exclude-reads in either step.
    #[test]
    fn no_low_complexity_list_no_exclude_reads() {
        let l = Layout::with_low_complexity(Path::new("/R"), false);
        let p = Params::large();
        assert!(!partner_args(&l, 4, Path::new("/o"), Path::new("/s"), &p).iter().any(|a| a == "--exclude-reads"));
        assert!(!rounds_args(&l, 4, Path::new("/w"), Path::new("/g"), &p, Gzip::Builtin { threads: 4 }).iter().any(|a| a == "--exclude-reads"));
        let l = Layout::with_low_complexity(Path::new("/R"), true);
        assert_eq!(l.low_complexity.as_deref(), Some(Path::new("/R/correction/low_complexity.tsv")));
    }

    #[test]
    fn final_round_is_the_one_with_sparse_tables() {
        let dir = std::env::temp_dir().join(format!("squash_final_round_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        let l = Layout::new(&dir);
        assert!(final_round(&l).is_err());
        let touch = |p: PathBuf| {
            fs::create_dir_all(p.parent().unwrap()).unwrap();
            File::create(p).unwrap();
        };
        for k in 1..=3 {
            touch(l.tokens(k));   // tokens without a sparse table
        }
        touch(l.partners());   // the correction folder's other files are not rounds
        touch(l.root.join(work::LOW_COMPLEXITY));
        assert!(final_round(&l).is_err());
        touch(l.sparse_dir(3).join("tokens.tsv.gz"));
        assert_eq!(final_round(&l), Ok(3));
        touch(l.sparse_dir(2).join("tokens.tsv.gz"));
        assert!(final_round(&l).is_err());
        // a rerun's clearing keeps the inputs and removes the rounds
        clear_outputs(&l).unwrap();
        assert!(l.rounds().unwrap().is_empty() && !l.partners().exists() && l.root.join(work::LOW_COMPLEXITY).exists());
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn select_args_name_the_tables() {
        assert_eq!(select_args(Path::new("/T/tokens.tsv.gz"), &Params::small().sparse, Path::new("/S/select/repeat.tsv"), 4),
                   s(&["sparse-repeat", "--tokens", "/T/tokens.tsv.gz", "--inverse-density", "8", "--fraction", "0.0001",
                       "--output", "/S/select/repeat.tsv", "--pigz", "4"]));
    }

    #[test]
    fn layout_names_the_work_paths() {
        let l = Layout::new(Path::new("/R"));
        assert_eq!(l.tokens(12), Path::new("/R/correction/round12/tokens.tsv.gz"));
        assert_eq!(l.sparse_dir(16), Path::new("/R/correction/round16/sparse"));
        assert_eq!(l.partners(), Path::new("/R/correction/partners.tsv.gz"));
    }

    #[test]
    fn round_tracker_keeps_the_round_a_fixed_point_needs() {
        // 12 rounds to the last
        let mut t = RoundTracker::default();
        assert_eq!(t.on_round(1, "continue").unwrap(), Vec::<u64>::new());
        for k in 2..12 {
            assert_eq!(t.on_round(k, "continue").unwrap(), vec![k - 1]);
        }
        assert_eq!(t.on_round(12, "last").unwrap(), vec![11]);
        assert_eq!(t.final_round, Some((12, "last".to_string())));
        assert!(t.on_round(13, "continue").is_err());
        // converged: that round is final
        let mut t = RoundTracker::default();
        t.on_round(1, "continue").unwrap();
        assert_eq!(t.on_round(2, "converged").unwrap(), vec![1]);
        assert_eq!(t.final_round, Some((2, "converged".to_string())));
        // a fixed point at round 3: round 2 is final, round 3's files go
        let mut t = RoundTracker::default();
        t.on_round(1, "continue").unwrap();
        t.on_round(2, "continue").unwrap();
        assert_eq!(t.on_round(3, "fixed_point").unwrap(), vec![3]);
        assert_eq!(t.final_round, Some((2, "fixed_point".to_string())));
        // a converged or last first round; a fixed point in round 1 has nothing to publish
        let mut t = RoundTracker::default();
        assert_eq!(t.on_round(1, "last").unwrap(), Vec::<u64>::new());
        assert_eq!(t.final_round, Some((1, "last".to_string())));
        assert!(RoundTracker::default().on_round(1, "fixed_point").is_err());
        // out of order or unknown
        let mut t = RoundTracker::default();
        assert!(t.on_round(2, "continue").is_err());
        assert!(RoundTracker::default().on_round(1, "stopped").is_err());
    }

    #[test]
    fn gzip_members_round_trip() {
        use flate2::read::MultiGzDecoder;
        let mut state: u64 = 0x2545_F491_4F6C_DD1D;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let text: Vec<u8> = (0..300_000).map(|i| if i % 61 == 60 { b'\n' } else { b'0' + (next() % 10) as u8 }).collect();
        for &(len, block, threads) in &[(0usize, 1000usize, 1usize), (1, 1000, 2), (999, 1000, 3), (1000, 1000, 2),
                                        (1001, 1000, 2), (300_000, 1000, 4), (300_000, 7, 3), (300_000, 1 << 20, 4),
                                        (64_000, 1000, 1)] {
            let mut gz = Vec::new();
            gzip_members(&text[..len], &mut gz, 1, threads, block).unwrap();
            let mut back = Vec::new();
            MultiGzDecoder::new(&gz[..]).read_to_end(&mut back).unwrap();
            assert_eq!(back, &text[..len], "len {} block {} threads {}", len, block, threads);
            // one member per block (at least one)
            let members = gz.windows(3).filter(|w| w == &[0x1f, 0x8b, 0x08]).count();
            assert!(members >= ((len + block - 1) / block).max(1), "len {} block {}", len, block);
        }
        // the same bytes on any thread count
        let mut a = Vec::new();
        let mut b = Vec::new();
        gzip_members(&text[..], &mut a, 6, 1, 4096).unwrap();
        gzip_members(&text[..], &mut b, 6, 4, 4096).unwrap();
        assert_eq!(a, b);
    }
}

// ------------------------------------------------------------------------------------------------ pipeline stage

/// The pipeline's stage `correct` (pipeline.rs): the partner search and the rounds in the run's work folder with the
/// run's values (a rerun of the stage first removes this stage's outputs). The later stages find the final round
/// with `final_round`.
pub fn run(r: &crate::pipeline::Run) -> Result<(), String> {
    let large = r.preset() == crate::presets::Preset::Large;
    let p = Params::from_values(large, &|name| r.v(name).to_string())?;
    let l = Layout::with_low_complexity(&r.root, crate::prepare::low_complexity_value(r.v("low-complexity"))?.is_some());
    clear_outputs(&l)?;
    // pigz when it is installed, else in process: the same text either way.
    let gzip = match check_gzip(Gzip::Pigz { threads: r.cfg.threads }) {
        Ok(()) => Gzip::Pigz { threads: r.cfg.threads },
        Err(_) => Gzip::Builtin { threads: r.cfg.threads },
    };
    correct(&r.runner, &l, r.cfg.threads, &r.tmp, &p, gzip).map(|_| ())
}

/// The final round's folder (`correction/roundN`) of a run whose `correct` stage is done: the assembly's `round`.
pub fn final_round_dir(r: &crate::pipeline::Run) -> Result<PathBuf, String> {
    let l = Layout::new(&r.root);
    Ok(l.round(final_round(&l)?))
}
