//! The pipeline: the steps in order, each run as `squash __STEP ARGS` in a child process (exec.rs) with fixed
//! arguments built from the resolved values, in the work folder `OUT/work`; then the final files are copied to `OUT`.
//!
//! Work folder: `signal/` (q-symbols and the events' POD5 sample spans), `tokens/` (the uncorrected tokens),
//! `correction/` (low-complexity reads, partners, `roundN/` corrected tokens and the last round's sparse table),
//! `bp_scale/`, `graph/`, `threaded/`, `layout/`, `logs/` (each step's log and command line), `stages/` (one marker
//! per finished stage).
//!
//! Stages (for `--stop-after` and `--resume`): signal (POD5 -> events -> q-symbols -> tokens), correct (partner
//! search and correction rounds), assemble (multi-k graph), thread (read threading and contig joins), layout (the
//! signal layout of each unitig).

use crate::exec::Runner;
use crate::presets::{Preset, Values, PARAMS};
use crate::{assemble, correct, prepare};
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

pub const STAGES: &[&str] = &["signal", "correct", "assemble", "thread", "layout"];

/// The work folder (OUT/work): the files that several stages read, as paths under it. The rest is named by the stage
/// that writes it: `correction/` (correct.rs), `graph/`, `threaded/`, `layout/` and the rest of `bp_scale/`
/// (assemble.rs); `logs/` holds each step's command and log, `stages/` the stages done.
pub mod work {
    /// Every read's q-symbols (the signal step).
    pub const QSYMBOLS: &str = "signal/qsymbols.tsv";
    /// Every event's POD5 samples (the signal step; read by the signal layout).
    pub const SPANS: &str = "signal/spans.bin";
    /// Every read's signal length (read by the bp scale).
    pub const SIGNAL_LENGTHS: &str = "bp_scale/signal_lengths.tsv";
    /// The uncorrected tokens (a folder of token arrays, `io::READS_TABLE` ...).
    pub const TOKENS: &str = "tokens";
    /// The low-complexity reads, left out of correction.
    pub const LOW_COMPLEXITY: &str = "correction/low_complexity.tsv";
    /// The first partner search's table (read by the rounds and the contig joins).
    pub const PARTNERS: &str = "correction/partners.tsv.gz";
}

#[derive(Clone, Debug)]
pub struct Config {
    pub reads: Vec<PathBuf>,
    pub out_dir: PathBuf,
    pub threads: usize,
    pub tmp_dir: Option<PathBuf>,
    pub read_ids: Option<PathBuf>,
    pub genome_size: Option<u64>,
    pub keep_work: bool,
    pub resume: bool,
    pub stop_after: Option<String>,
    pub threading: bool,
    pub tangles: bool,
    pub joins: bool,
    pub values: Values,
    pub command_line: String,
}

/// Paths of a run and the step runner.
pub struct Run<'a> {
    pub cfg: &'a Config,
    /// The work folder `OUT/work` (`work`).
    pub root: PathBuf,
    /// Large temporary files.
    pub tmp: PathBuf,
    /// The POD5 input as one path (a file, or a folder of the inputs in the user's order).
    pub pod5: PathBuf,
    /// `--read-ids` as a plain list, if given.
    pub read_ids: Option<PathBuf>,
    pub runner: Runner,
}

impl<'a> Run<'a> {
    pub fn v(&self, name: &str) -> &str {
        self.cfg.values.get(name)
    }
    pub fn preset(&self) -> Preset {
        self.cfg.values.preset
    }
    pub fn threads(&self) -> String {
        self.cfg.threads.to_string()
    }
}

/// Every resolved value of a run, as written to OUT/params.tsv and printed by `--dry-run`.
pub fn params_tsv(cfg: &Config) -> String {
    let mut s = String::from("name\tvalue\tsource\n");
    s += &format!("version\t{}\tsquash\n", crate::version());
    s += &format!("preset\t{}\tuser\n", cfg.values.preset.name());
    for r in &cfg.reads {
        s += &format!("reads\t{}\tuser\n", r.display());
    }
    if let Some(p) = &cfg.read_ids {
        s += &format!("read-ids\t{}\tuser\n", p.display());
    }
    if let Some(g) = cfg.genome_size {
        s += &format!("genome-size\t{}\tuser\n", g);
    }
    s += &format!("threads\t{}\t{}\n", cfg.threads, "run (no effect on outputs)");
    for p in PARAMS {
        let (v, user) = cfg.values.iter().find(|(k, _, _)| *k == p.name).map(|(_, v, u)| (v.to_string(), u)).unwrap();
        s += &format!("{}\t{}\t{}\n", p.name, v, if user { "user" } else { cfg.values.preset.name() });
    }
    for (name, on) in [("threading", cfg.threading), ("tangles", cfg.tangles), ("joins", cfg.joins)] {
        s += &format!("{}\t{}\t{}\n", name, if on { "on" } else { "off" }, if on { cfg.values.preset.name() } else { "user" });
    }
    s
}

/// The values a stage's outputs depend on, for `--resume`: everything but the thread count and the output folders.
fn fingerprint(cfg: &Config) -> String {
    params_tsv(cfg).lines().filter(|l| !l.starts_with("threads\t")).collect::<Vec<_>>().join("\n")
}

fn stage_marker(root: &Path, stage: &str) -> PathBuf {
    root.join("stages").join(format!("{}.done", stage))
}

pub fn run(cfg: &Config) -> Result<(), String> {
    for r in &cfg.reads {
        if !r.exists() {
            return Err(format!("input {} does not exist", r.display()));
        }
    }
    let out = &cfg.out_dir;
    let root = out.join("work");
    if !cfg.resume {
        // A run without --resume starts over: an earlier work folder (one SQUASH made: it has stages/) and earlier
        // outputs are removed first, so that no file of another run is taken for this one's.
        if root.exists() {
            if !root.join("stages").is_dir() {
                return Err(format!("{} exists and is not a SQUASH work folder: choose another -o", root.display()));
            }
            fs::remove_dir_all(&root).map_err(|e| format!("cannot remove the earlier {}: {}", root.display(), e))?;
        }
        for name in assemble::OUTPUT_NAMES {
            let _ = fs::remove_file(out.join(name));
        }
    }
    fs::create_dir_all(root.join("stages")).map_err(|e| format!("cannot create {}: {}", root.display(), e))?;
    let tmp = cfg.tmp_dir.clone().unwrap_or_else(|| root.join("tmp"));
    fs::create_dir_all(&tmp).map_err(|e| format!("cannot create {}: {}", tmp.display(), e))?;
    fs::write(out.join("params.tsv"), params_tsv(cfg)).map_err(|e| format!("cannot write params.tsv: {}", e))?;
    let mut log = Log::open(&out.join("squash.log"))?;
    log.line(&format!("squash {} | {}", crate::version(), cfg.command_line));
    let runner = Runner::new(&root.join("logs")).map_err(|e| format!("cannot create the logs folder: {}", e))?;
    let pod5 = input_path(cfg, &root)?;
    let read_ids = match &cfg.read_ids {
        Some(p) => Some(plain_list(p, &root.join("input_read_ids.txt"))?),
        None => None,
    };
    let r = Run { cfg, root: root.clone(), tmp, pod5, read_ids, runner };

    let print = fingerprint(cfg);
    let mut rerun = !cfg.resume;
    for &stage in STAGES {
        if stage == "thread" && !cfg.threading {
            continue;
        }
        let marker = stage_marker(&root, stage);
        let done = !rerun && fs::read_to_string(&marker).map(|s| s == print).unwrap_or(false);
        if done {
            log.line(&format!("stage {}: kept from the earlier run (--resume)", stage));
        } else {
            rerun = true;
            let _ = fs::remove_file(&marker);
            log.line(&format!("stage {}: start", stage));
            let t = std::time::Instant::now();
            match stage {
                "signal" => {
                    signal_stage(&r, &mut log)?;
                    prepare::run(&r)?;
                }
                "correct" => correct::run(&r)?,
                "assemble" => assemble::run(&r)?,
                "thread" => assemble::thread_stage(&r)?,
                "layout" => assemble::layout_stage(&r)?,
                _ => unreachable!(),
            }
            fs::write(&marker, &print).map_err(|e| format!("cannot write {}: {}", marker.display(), e))?;
            log.line(&format!("stage {}: done in {:.1} s", stage, t.elapsed().as_secs_f64()));
        }
        if cfg.stop_after.as_deref() == Some(stage) {
            log.line(&format!("stopped after stage {} (--stop-after); work kept in {}", stage, root.display()));
            return Ok(());
        }
    }
    export(&r, &mut log)?;
    if !cfg.keep_work {
        fs::remove_dir_all(&root).map_err(|e| format!("cannot remove {}: {}", root.display(), e))?;
        if let Some(t) = &cfg.tmp_dir {
            let _ = fs::remove_dir(t);
        }
    }
    log.line("done");
    Ok(())
}

/// The steps take one POD5 path: a file, or a folder whose *.pod5 files are read in name order. One input is passed
/// as it is; several become a folder of links (`input_pod5/`) whose names keep the user's order.
fn input_path(cfg: &Config, root: &Path) -> Result<PathBuf, String> {
    if cfg.reads.len() == 1 {
        return Ok(cfg.reads[0].clone());
    }
    let dir = root.join("input_pod5");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).map_err(|e| format!("cannot create {}: {}", dir.display(), e))?;
    let mut n = 0usize;
    for r in &cfg.reads {
        let files: Vec<PathBuf> = if r.is_dir() {
            let mut v: Vec<PathBuf> = fs::read_dir(r)
                .map_err(|e| format!("cannot read {}: {}", r.display(), e))?
                .filter_map(|e| e.ok().map(|e| e.path()))
                .filter(|p| p.is_file() && p.extension().and_then(|s| s.to_str()) == Some("pod5"))
                .collect();
            v.sort();
            v
        } else {
            vec![r.clone()]
        };
        for f in files {
            let abs = fs::canonicalize(&f).map_err(|e| format!("{}: {}", f.display(), e))?;
            let link = dir.join(format!("{:06}_{}", n, f.file_name().unwrap().to_string_lossy()));
            std::os::unix::fs::symlink(&abs, &link).map_err(|e| format!("cannot link {}: {}", link.display(), e))?;
            n += 1;
        }
    }
    if n == 0 {
        return Err("no POD5 files in the inputs".into());
    }
    Ok(dir)
}

/// `--read-ids` as the plain list the signal step reads (a gzipped list is decompressed).
fn plain_list(p: &Path, dst: &Path) -> Result<PathBuf, String> {
    let f = fs::File::open(p).map_err(|e| format!("--read-ids {}: {}", p.display(), e))?;
    let mut head = [0u8; 2];
    let gz = {
        use std::io::Read;
        let mut g = fs::File::open(p).map_err(|e| format!("--read-ids {}: {}", p.display(), e))?;
        g.read(&mut head).map(|n| n == 2 && head == [0x1f, 0x8b]).unwrap_or(false)
    };
    if !gz {
        return Ok(p.to_path_buf());
    }
    let mut out = std::io::BufWriter::new(fs::File::create(dst).map_err(|e| format!("{}: {}", dst.display(), e))?);
    for line in BufReader::new(flate2::read::MultiGzDecoder::new(f)).lines() {
        let line = line.map_err(|e| format!("--read-ids {}: {}", p.display(), e))?;
        writeln!(out, "{}", line).map_err(|e| e.to_string())?;
    }
    out.flush().map_err(|e| e.to_string())?;
    Ok(dst.to_path_buf())
}

/// The signal step: POD5 -> events -> q-symbols in one pass (`signal/qsymbols.tsv` and the event spans
/// `signal/spans.bin`; reads in blocks of 40,000 of file order, each block sorted by read id), then each read's signal
/// length (`bp_scale/signal_lengths.tsv`, read by the bp scale) and the coverage estimate.
fn signal_stage(r: &Run, log: &mut Log) -> Result<(), String> {
    fs::create_dir_all(r.root.join("signal")).map_err(|e| e.to_string())?;
    let (w1, w2) = split2(r.v("evdet-win"));
    let (t1, t2) = split2(r.v("evdet-thr"));
    let mut args: Vec<String> = vec![
        r.pod5.display().to_string(),
        "--evdet-win1".into(),
        w1,
        "--evdet-win2".into(),
        w2,
        "--evdet-thr1".into(),
        t1,
        "--evdet-thr2".into(),
        t2,
        "--evdet-peak".into(),
        r.v("evdet-peak").into(),
        "--evdet-event-merge".into(),
        r.v("event-merge").into(),
        "--quant-bit".into(),
        r.v("quant-bits").into(),
        "--threads".into(),
        r.threads(),
        "--signal-q-symbols".into(),
        r.root.join(work::QSYMBOLS).display().to_string(),
    ];
    if let Some(ids) = &r.read_ids {
        args.push("--keep-reads".into());
        args.push(ids.display().to_string());
    }
    args.extend([
        // each event's POD5 samples, for the signal layout
        "--dump-event-spans".into(),
        r.root.join(work::SPANS).display().to_string(),
    ]);
    r.runner.step("signal", "signal", &args)?;
    let lengths = r.root.join(work::SIGNAL_LENGTHS);
    fs::create_dir_all(r.root.join("bp_scale")).map_err(|e| e.to_string())?;
    r.runner.step(
        "bp-scale-lengths",
        "bp-scale",
        &["lengths".into(), "--pod5".into(), r.pod5.display().to_string(), "--output".into(), lengths.display().to_string()],
    )?;
    coverage_report(r, &lengths, log)
}

fn split2(v: &str) -> (String, String) {
    let mut it = v.split(',');
    (it.next().unwrap_or("").to_string(), it.next().unwrap_or("").to_string())
}

/// Signal-estimated bases of the input (samples x bases per second / sample rate, the reads assembled only), and the
/// coverage per strand when the genome size is known.
fn coverage_report(r: &Run, lengths: &Path, log: &mut Log) -> Result<(), String> {
    let bps: f64 = r.v("bases-per-second").parse().unwrap();
    let keep: Option<std::collections::HashSet<String>> = match &r.read_ids {
        Some(p) => Some(fs::read_to_string(p).map_err(|e| e.to_string())?.lines().map(|l| l.trim().to_string()).collect()),
        None => None,
    };
    let f = fs::File::open(lengths).map_err(|e| format!("{}: {}", lengths.display(), e))?;
    let (mut reads, mut bases) = (0u64, 0f64);
    for line in BufReader::new(f).lines().skip(1) {
        let line = line.map_err(|e| e.to_string())?;
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 3 || keep.as_ref().map(|k| !k.contains(c[0])).unwrap_or(false) {
            continue;
        }
        let (s, rate): (f64, f64) = (c[1].parse().unwrap_or(0.0), c[2].parse().unwrap_or(0.0));
        if rate > 0.0 {
            reads += 1;
            bases += s * bps / rate;
        }
    }
    let mut msg = format!("input: {} reads, {:.1} Mb signal-estimated", reads, bases / 1e6);
    if let Some(g) = r.cfg.genome_size {
        let per_strand = bases / (2.0 * g as f64);
        msg += &format!(", about {:.1}x per strand of {:.2} Mb", per_strand, g as f64 / 1e6);
        log.line(&msg);
        if per_strand < 20.0 {
            log.warn("below 20x per strand: SQUASH needs about 20x per strand (each strand assembles on its own)");
        }
        let large = g >= 10_000_000;
        if large != (r.preset() == Preset::Large) {
            log.warn(&format!(
                "a {:.1}-Mb genome is {} 10 Mb, where the paper's rule picks -x {}",
                g as f64 / 1e6,
                if large { "over" } else { "under" },
                if large { "large" } else { "small" }
            ));
        }
    } else {
        log.line(&msg);
    }
    Ok(())
}

/// The files a user gets, copied from the work folder.
fn export(r: &Run, log: &mut Log) -> Result<(), String> {
    let files = assemble::final_files(r)?;
    for (src, name) in files {
        let dst = r.cfg.out_dir.join(&name);
        let tmp = r.cfg.out_dir.join(format!(".{}.part", name));
        fs::copy(&src, &tmp).map_err(|e| format!("cannot copy {} to {}: {}", src.display(), tmp.display(), e))?;
        fs::rename(&tmp, &dst).map_err(|e| format!("cannot write {}: {}", dst.display(), e))?;
        log.line(&format!("output {}", dst.display()));
    }
    Ok(())
}

/// squash.log: each line also goes to the terminal.
pub struct Log {
    f: fs::File,
}

impl Log {
    fn open(p: &Path) -> Result<Log, String> {
        let f = fs::OpenOptions::new().create(true).append(true).open(p).map_err(|e| format!("{}: {}", p.display(), e))?;
        Ok(Log { f })
    }
    pub fn line(&mut self, s: &str) {
        let stamp = chrono_like_now();
        eprintln!("[squash {}] {}", stamp, s);
        let _ = writeln!(self.f, "[{}] {}", stamp, s);
    }
    pub fn warn(&mut self, s: &str) {
        self.line(&format!("WARNING: {}", s));
    }
}

/// Local wall-clock time as YYYY-MM-DD HH:MM:SS (no date crate; from libc's localtime).
fn chrono_like_now() -> String {
    unsafe {
        let t = libc::time(std::ptr::null_mut());
        let mut tm: libc::tm = std::mem::zeroed();
        libc::localtime_r(&t, &mut tm);
        format!("{:04}-{:02}-{:02} {:02}:{:02}:{:02}", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec)
    }
}
