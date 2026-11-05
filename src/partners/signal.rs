//! The q-symbol alignment check of candidate partners (off unless `--signal-symbols` is given), nanoMDBG's overhang
//! rule on the signal. A candidate partner (a chain that passed every other test) is aligned to its target on the
//! reads' q-symbols (one per detected event) over the overlap its chain implies, by a weighted edit distance
//! (substitution = the level difference, insertion / deletion `--signal-edit-gap`).
//!
//! The seed: the `--signal-window`-symbol stretch of the overlap at the chain's middle, placed on the partner from the
//! chain anchor nearest to it (exact token matches give both reads' q positions) and refined by a semi-global search in
//! a diagonal band of +- `--signal-band` around that start (the start kept even when the band's edge holds the best).
//! From the seed's start cell one path runs forward to the overlap's end and one backward to its start, each in a band
//! that follows the path (centred one column after the previous row's cheapest cell, moving at most
//! `--signal-path-shift` columns per row beyond that; free end): first +- B / 2 symbols, and +- B (`--signal-path B`)
//! for a pair that fails in the narrower band. The partner passes when the whole path costs at most `--signal-thr` per
//! 1,000 alignment columns (one number, as nanoMDBG's identity). Verdicts are kept per directed pair (target, donor),
//! each computed from that target's own chain, so they do not depend on which side reaches a pair first;
//! `--signal-verdicts` loads the verdicts of an earlier search, which are reused.

use super::Result;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
use std::sync::{Mutex, OnceLock};

const SHARDS: usize = 256;

// `SQUASH_QSYM_PROFILE=1` (a diagnostic; the verdicts do not change): thread-ns and counts of the check's phases,
// printed by `report`: seed (ns, calls, found inside the band), attempt in half the band (ns), in the full band
// (ns), `tracked` (ns, cells), `trace` (ns), `tracked` calls.
const PROF_ZERO: AtomicU64 = AtomicU64::new(0);
static PROF: [AtomicU64; 9] = [PROF_ZERO; 9];
static PROF_ON: OnceLock<bool> = OnceLock::new();

fn prof_on() -> bool {
    *PROF_ON.get_or_init(|| std::env::var("SQUASH_QSYM_PROFILE").is_ok_and(|v| v == "1"))
}

fn prof_add(i: usize, v: u64) {
    PROF[i].fetch_add(v, Relaxed);
}

/// One read's q-symbols (0-7), 3 bits each (8 symbols in 3 bytes). `unpack_into` restores them for the alignment.
#[derive(Clone, Default)]
pub struct Syms {
    len: usize,
    data: Vec<u8>,
}

impl Syms {
    pub fn new(v: Vec<u8>) -> std::result::Result<Syms, String> {
        if let Some(&s) = v.iter().find(|&&s| s > 7) {
            return Err(format!("q-symbol {} does not fit 3 bits", s));
        }
        let mut data = vec![0u8; (3 * v.len() + 7) / 8 + 1];
        for (i, &s) in v.iter().enumerate() {
            let bit = 3 * i;
            let w = (s as u16) << (bit % 8);
            data[bit / 8] |= w as u8;
            data[bit / 8 + 1] |= (w >> 8) as u8;
        }
        data.shrink_to_fit();
        Ok(Syms { len: v.len(), data })
    }
    pub fn len(&self) -> usize {
        self.len
    }
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
    /// The symbols, unpacked into `out` (cleared first): 8 at a time from each 3 bytes.
    pub fn unpack_into(&self, out: &mut Vec<u8>) {
        unpack_bytes(self.len, &self.data, out);
    }
}

/// `Syms::unpack_into` on a read's stored bytes (`len` symbols, 3 bits each).
fn unpack_bytes(len: usize, data: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.reserve(len);
    let full = len / 8;
    for c in 0..full {
        let b = &data[3 * c..3 * c + 3];
        let w = b[0] as u32 | (b[1] as u32) << 8 | (b[2] as u32) << 16;
        for j in 0..8 {
            out.push(((w >> (3 * j)) & 7) as u8);
        }
    }
    for i in 8 * full..len {
        let bit = 3 * i;
        let w = data[bit / 8] as u16 | (data[bit / 8 + 1] as u16) << 8;
        out.push(((w >> (bit % 8)) & 7) as u8);
    }
}

/// Where the check's q-symbols are: in memory, or written once to a node-local file (`store_on_file`) from which each
/// checked pair's two reads are read (`pread`) into the thread's buffers, so the search does not hold them in memory.
/// The file is unlinked once open, so nothing is left behind. Either way the same symbols, so the same verdicts.
pub enum SymStore {
    Memory(Vec<Syms>),
    File { file: std::fs::File, entries: Vec<SymEntry> },
}

/// One read's place in the file: offset, symbols, stored bytes.
pub struct SymEntry {
    offset: u64,
    len: u32,
    bytes: u32,
}

impl SymStore {
    fn is_empty(&self) -> bool {
        match self {
            SymStore::Memory(v) => v.is_empty(),
            SymStore::File { entries, .. } => entries.is_empty(),
        }
    }

    /// The symbols written to DIR/.check_qsymbols.PID.bin (then opened and unlinked).
    fn to_file(symbols: &[Syms], dir: &Path) -> Result<SymStore> {
        let path = dir.join(format!(".check_qsymbols.{}.bin", std::process::id()));
        let mut out = BufWriter::with_capacity(1 << 24, std::fs::File::create(&path)?);
        let (mut entries, mut offset) = (Vec::with_capacity(symbols.len()), 0u64);
        for s in symbols {
            out.write_all(&s.data)?;
            let len = <u32 as std::convert::TryFrom<usize>>::try_from(s.len).map_err(|_| "a read with 2^32 q-symbols or more")?;
            let bytes = <u32 as std::convert::TryFrom<usize>>::try_from(s.data.len()).map_err(|_| "a read with 2^32 q-symbol bytes or more")?;
            entries.push(SymEntry { offset, len, bytes });
            offset += s.data.len() as u64;
        }
        out.flush()?;
        drop(out);
        let file = std::fs::File::open(&path)?;
        std::fs::remove_file(&path)?;
        Ok(SymStore::File { file, entries })
    }

    /// Read `i`'s symbols, unpacked into `out` (`raw`: this thread's byte buffer).
    fn unpacked(&self, i: usize, raw: &mut Vec<u8>, out: &mut Vec<u8>) {
        match self {
            SymStore::Memory(v) => v[i].unpack_into(out),
            SymStore::File { file, entries } => {
                use std::os::unix::fs::FileExt;
                let e = &entries[i];
                raw.resize(e.bytes as usize, 0);
                file.read_exact_at(raw, e.offset).expect("reading the check's q-symbols from their local file");
                unpack_bytes(e.len as usize, raw, out);
            }
        }
    }
}

thread_local! {
    static UNPACKED: std::cell::RefCell<(Vec<u8>, Vec<u8>)> = std::cell::RefCell::new((Vec::new(), Vec::new()));
    /// A read's stored bytes, read from the file (`SymStore::File`).
    static RAW_BYTES: std::cell::RefCell<Vec<u8>> = std::cell::RefCell::new(Vec::new());
}

pub struct Check {
    // the reads' q-symbols, behind a lock so that they can be moved to a local file and read from it pair by pair
    // (`SymStore`, `store_on_file`)
    symbols: std::sync::RwLock<SymStore>,
    window: usize,
    thr: f32,
    // the seed's band half-width
    band: i64,
    // the path's band half-width and largest extra shift per row
    path: i64,
    path_shift: i64,
    // the gap cost of the weighted edit distance
    edit_gap: i32,
    shards: Vec<Mutex<HashMap<u64, bool>>>,
    loaded: usize,
    // pairs checked, passed, verdicts reused, thread-ns, rows aligned, pairs aligned again in the full band
    stats: [AtomicU64; 6],
}

static CHECK: OnceLock<Check> = OnceLock::new();

/// The installed check, if any.
pub fn get() -> Option<&'static Check> {
    CHECK.get()
}

pub struct Options {
    symbols: PathBuf,
    verdicts: Option<PathBuf>,
    pub verdicts_out: Option<PathBuf>,
    /// `--signal-store DIR`: the q-symbols moved to a file in DIR after loading (`store_on_file`)
    pub store: Option<PathBuf>,
    window: usize,
    thr: f32,
    band: i64,
    path: i64,
    path_shift: i64,
    edit_gap: i64,
}

/// The `--signal-*` options, removed from MAP; None without `--signal-symbols` (the q-symbols table).
/// `--signal-thr` (per 1,000 alignment columns), `--signal-band`, `--signal-path` and `--signal-edit-gap` are
/// required with it.
pub fn take(map: &mut HashMap<String, String>) -> Result<Option<Options>> {
    let symbols = map.remove("signal-symbols").map(PathBuf::from);
    let verdicts = map.remove("signal-verdicts").map(PathBuf::from);
    let verdicts_out = map.remove("signal-verdicts-out").map(PathBuf::from);
    let store = map.remove("signal-store").map(PathBuf::from);
    if let Some(d) = &store {
        if !d.is_dir() {
            return Err(format!("--signal-store {}: not a directory", d.display()).into());
        }
    }
    let num = |map: &mut HashMap<String, String>, key: &str, d: f64| -> Result<f64> {
        match map.remove(key) {
            None => Ok(d),
            Some(v) => v.parse::<f64>().map_err(|e| format!("--{} {:?}: {}", key, v, e).into()),
        }
    };
    let window = num(map, "signal-window", 1000.0)?;
    let thr = num(map, "signal-thr", f64::NAN)?;
    let band = num(map, "signal-band", 0.0)?;
    let path = num(map, "signal-path", 0.0)?;
    let path_shift = num(map, "signal-path-shift", 1.0)?;
    let edit_gap = num(map, "signal-edit-gap", 0.0)?;
    let symbols = match symbols {
        None => {
            if verdicts.is_some() || verdicts_out.is_some() {
                return Err("--signal-verdicts / --signal-verdicts-out need --signal-symbols".into());
            }
            return Ok(None);
        }
        Some(q) => q,
    };
    if !(edit_gap >= 1.0) || edit_gap.fract() != 0.0 || edit_gap > 1000.0 {
        return Err("--signal-edit-gap: an integer from 1 to 1000".into());
    }
    if !(window >= 10.0) || !(thr > 0.0) || !(band >= 1.0) {
        return Err("--signal-window >= 10, --signal-band >= 1, --signal-thr > 0".into());
    }
    if !(path >= 1.0) || !(path_shift >= 1.0) || path.fract() != 0.0 || path_shift.fract() != 0.0 {
        return Err("--signal-path >= 1 and --signal-path-shift >= 1, integers".into());
    }
    if let Some(p) = &verdicts_out {
        if p.exists() {
            return Err(format!("refusing to overwrite {}", p.display()).into());
        }
    }
    Ok(Some(Options { symbols, verdicts, verdicts_out, store, window: window as usize, thr: thr as f32, band: band as i64,
                      path: path as i64, path_shift: path_shift as i64, edit_gap: edit_gap as i64 }))
}

/// The key of the directed pair (target, donor).
fn directed_key(target: usize, donor: usize) -> u64 {
    ((target as u64) << 32) | donor as u64
}

fn shard(k: u64) -> usize {
    (k.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 56) as usize % SHARDS
}

/// Loads the q-symbols of the catalog's reads (the q-symbols table: `#` header lines, then read_id and
/// comma-separated symbols; reads absent from the file get none and always pass) and the earlier verdicts, and
/// installs the check.
pub fn install(o: &Options, ids: &[String], pool: &rayon::ThreadPool) -> Result<()> {
    let clock = std::time::Instant::now();
    let ordinal: HashMap<&str, usize> = ids.iter().enumerate().map(|(i, x)| (x.as_str(), i)).collect();
    let mut symbols: Vec<Syms> = vec![Syms::default(); ids.len()];
    let reader = BufReader::with_capacity(1 << 24, std::fs::File::open(&o.symbols)
        .map_err(|e| format!("{}: {}", o.symbols.display(), e))?);
    let mut batch: Vec<String> = Vec::with_capacity(256);
    let flush = |batch: &mut Vec<String>, symbols: &mut Vec<Syms>| -> Result<()> {
        let parsed: Vec<Result<Option<(usize, Syms)>>> = pool.install(|| batch.par_iter().map(|line| {
            if line.starts_with('#') {
                return Ok(None);
            }
            let (id, rest) = line.split_once('\t').unwrap_or((line.as_str(), ""));
            match ordinal.get(id) {
                None => Ok(None),
                Some(&r) => {
                    let v: std::result::Result<Vec<u8>, _> = rest.trim_end().split(',').filter(|x| !x.is_empty())
                        .map(|x| x.parse::<u8>()).collect();
                    let v = v.map_err(|e| format!("q-symbols of {}: {}", id, e))?;
                    Ok(Some((r, Syms::new(v).map_err(|e| format!("q-symbols of {}: {}", id, e))?)))
                }
            }
        }).collect());
        for p in parsed {
            if let Some((r, v)) = p? {
                symbols[r] = v;
            }
        }
        batch.clear();
        Ok(())
    };
    for line in reader.lines() {
        batch.push(line?);
        if batch.len() == 256 {
            flush(&mut batch, &mut symbols)?;
        }
    }
    flush(&mut batch, &mut symbols)?;
    let (with, n_symbols) = (symbols.iter().filter(|e| !e.is_empty()).count(), symbols.iter().map(Syms::len).sum::<usize>());
    let shards: Vec<Mutex<HashMap<u64, bool>>> = (0..SHARDS).map(|_| Mutex::new(HashMap::new())).collect();
    let mut loaded = 0;
    if let Some(path) = &o.verdicts {
        for line in BufReader::new(std::fs::File::open(path).map_err(|e| format!("{}: {}", path.display(), e))?).lines() {
            let line = line?;
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() != 3 || f[0] == "target" {
                continue;
            }
            let (a, b) = match (ordinal.get(f[0]), ordinal.get(f[1])) {
                (Some(&a), Some(&b)) => (a, b),
                _ => return Err(format!("{}: read pair {} {} not in the catalog", path.display(), f[0], f[1]).into()),
            };
            let k = directed_key(a, b);
            shards[shard(k)].lock().unwrap().insert(k, f[2] == "1");
            loaded += 1;
        }
    }
    eprintln!("[signal] q-symbols of {} of {} reads ({} values) and {} verdicts loaded in {:.1} s; window {}, thr {}, \
               seed band {}, path band {} shift {}, weighted edit distance, gap {}", with, ids.len(), n_symbols, loaded,
              clock.elapsed().as_secs_f64(), o.window, o.thr, o.band, o.path, o.path_shift, o.edit_gap);
    CHECK.set(Check { symbols: std::sync::RwLock::new(SymStore::Memory(symbols)), window: o.window, thr: o.thr,
                      band: o.band, path: o.path, path_shift: o.path_shift, edit_gap: o.edit_gap as i32, shards, loaded,
                      stats: Default::default() })
        .map_err(|_| "the signal check is installed twice".into())
}

/// Prints and resets the counters.
pub fn report(label: &str) {
    if let Some(c) = get() {
        let s: Vec<u64> = c.stats.iter().map(|x| x.swap(0, Relaxed)).collect();
        eprintln!("[signal] {}: {} pairs checked, {} passed ({:.1}%), {} verdicts reused, {} rows aligned, {} aligned again \
                   in the full band, {:.1} thread-s", label, s[0], s[1], 100.0 * s[1] as f64 / s[0].max(1) as f64, s[2], s[4],
                  s[5], s[3] as f64 / 1e9);
        let (runs, falls) = simd16::take();
        if runs > 0 {
            eprintln!("[signal] {}: {} path alignments with 16-bit lanes, {} of them again with 32-bit lanes", label,
                      runs, falls);
        }
        if prof_on() {
            let p: Vec<u64> = PROF.iter().map(|x| x.swap(0, Relaxed)).collect();
            let ts = |ns: u64| ns as f64 / 1e9;
            eprintln!("[signal-profile] {}: seed {:.1} thread-s ({} calls, {} inside the band); attempt in half the band \
                       {:.1} thread-s, in the full band {:.1} thread-s; of which tracked {:.1} thread-s ({} calls, {:.3e} \
                       cells, {:.2} ns per cell), trace {:.1} thread-s", label, ts(p[0]), p[1], p[2], ts(p[3]), ts(p[4]),
                      ts(p[5]), p[8], p[6] as f64, p[5] as f64 / p[6].max(1) as f64, ts(p[7]));
        }
    }
}

/// Frees this thread's alignment buffers (called on every pool thread at the end of a search; the next search grows
/// them again).
pub fn release_thread_buffers() {
    PBUF.with(|pb| *pb.borrow_mut() = PathScratch::default());
    UNPACKED.with(|cell| *cell.borrow_mut() = (Vec::new(), Vec::new()));
    RAW_BYTES.with(|cell| *cell.borrow_mut() = Vec::new());
    simd16::release();
}

/// Moves the check's q-symbols to a file in DIR (node-local scratch); each pair's reads are then read from it when
/// checked (`SymStore`), so the searches do not hold them in memory. No-op without a check or when already on file.
pub fn store_on_file(dir: &Path) -> Result<()> {
    let c = match get() {
        Some(c) => c,
        None => return Ok(()),
    };
    let mut store = c.symbols.write().unwrap();
    if let SymStore::Memory(v) = &*store {
        let clock = std::time::Instant::now();
        let on_file = SymStore::to_file(v, dir)?;
        *store = on_file;
        unsafe {
            libc::malloc_trim(0);
        }
        eprintln!("[signal] q-symbols moved to a local file in {} ({:.1} s)", dir.display(), clock.elapsed().as_secs_f64());
    }
    Ok(())
}

/// Forgets every cached verdict (the correction rounds' last round, once no search follows), to free their memory.
pub fn release_cache() {
    if let Some(c) = get() {
        for s in &c.shards {
            let mut m = s.lock().unwrap();
            m.clear();
            m.shrink_to_fit();
        }
    }
}

/// Every verdict known (loaded and computed), sorted: target, donor (catalog order), pass (1 / 0).
pub fn write_verdicts(path: &Path, ids: &[String]) -> Result<()> {
    let c = get().ok_or("no signal check installed")?;
    let mut all: Vec<(u64, bool)> = Vec::new();
    for s in &c.shards {
        all.extend(s.lock().unwrap().iter().map(|(&k, &v)| (k, v)));
    }
    all.sort_unstable();
    let tmp = path.with_extension("tmp");
    let mut out = BufWriter::new(std::fs::File::create(&tmp)?);
    writeln!(out, "target\tdonor\tpass")?;
    for (k, v) in all {
        writeln!(out, "{}\t{}\t{}", ids[(k >> 32) as usize], ids[(k & 0xffff_ffff) as usize], v as u8)?;
    }
    out.flush()?;
    drop(out);
    std::fs::rename(&tmp, path)?;
    eprintln!("[signal] {} verdicts written ({} loaded before)", c.shards.iter().map(|s| s.lock().unwrap().len()).sum::<usize>(), c.loaded);
    Ok(())
}

impl Check {
    /// Tests: the check with these q-symbols.
    #[cfg(test)]
    fn set_symbols(&mut self, v: Vec<Syms>) {
        *self.symbols.get_mut().unwrap() = SymStore::Memory(v);
    }

    /// The verdict of the directed pair from the cache, without computing (None: not cached), counted as reused
    /// when found, so the caller can skip building the chain's anchors for a pair already checked. Only the
    /// target's own thread checks a directed pair, so a miss here is a miss in `verdict` too.
    pub fn cached(&self, r: usize, d: usize) -> Option<bool> {
        let k = directed_key(r, d);
        let v = self.shards[shard(k)].lock().unwrap().get(&k).copied();
        if v.is_some() {
            self.stats[2].fetch_add(1, Relaxed);
        }
        v
    }

    /// The directed pair's verdict (cached, or computed from the target's own chain): chain `line` = (target q of
    /// the first anchor, donor q, target q of the last, donor q), its anchors as (target q, donor q) in chain order
    /// (at least one), the implied overlap [ts, te) on the target.
    pub fn verdict(&self, r: usize, d: usize, line: (i64, i64, i64, i64), ts: i64, te: i64, anchors: &[(i64, i64)]) -> bool {
        let k = directed_key(r, d);
        let sh = &self.shards[shard(k)];
        if let Some(&v) = sh.lock().unwrap().get(&k) {
            self.stats[2].fetch_add(1, Relaxed);
            return v;
        }
        let clock = std::time::Instant::now();
        let v = self.path_compute(r, d, line, ts, te, anchors);
        self.stats[0].fetch_add(1, Relaxed);
        self.stats[1].fetch_add(v as u64, Relaxed);
        self.stats[3].fetch_add(clock.elapsed().as_nanos() as u64, Relaxed);
        sh.lock().unwrap().insert(k, v);
        v
    }

    fn path_compute(&self, a: usize, b: usize, line: (i64, i64, i64, i64), lo: i64, hi: i64, anchors: &[(i64, i64)]) -> bool {
        let pass = self.path_compute_inner(a, b, line, lo, hi, anchors);
        trim_path_buffers();
        pass
    }

    fn path_compute_inner(&self, a: usize, b: usize, line: (i64, i64, i64, i64), lo: i64, hi: i64, anchors: &[(i64, i64)]) -> bool {
        let store = self.symbols.read().unwrap();
        assert!(!store.is_empty(), "the check's q-symbols are missing");
        // both reads unpacked into this thread's buffers (linear in their lengths; the alignment is not), read from
        // the file first when they are there
        UNPACKED.with(|cell| {
            let (x, y) = &mut *cell.borrow_mut();
            RAW_BYTES.with(|raw| {
                let raw = &mut *raw.borrow_mut();
                store.unpacked(a, raw, x);
                store.unpacked(b, raw, y);
            });
            self.path_core(x, y, line, lo, hi, anchors)
        })
    }

    fn path_core(&self, x: &[u8], y: &[u8], line: (i64, i64, i64, i64), lo: i64, hi: i64, anchors: &[(i64, i64)]) -> bool {
        if x.is_empty() || y.is_empty() {
            return true;
        }
        let w = self.window as i64;
        let lo = lo.max(0);
        let hi = hi.min(x.len() as i64);
        let n = hi - lo;
        if n < w {
            return true;
        }
        let (a0, _, a1, _) = line;
        let ylen = y.len() as i64;
        let mid = ((a0 + a1) / 2 - lo).clamp(0, n - 1);
        let sm = (mid - w / 2).clamp(0, n - w);
        let xw = &x[(lo + sm) as usize..(lo + sm + w) as usize];
        // the seed: the chain anchor nearest to the window's centre (the first of equals) maps the window's start onto
        // y; a band search around it refines the start (kept even at the band's edge)
        let centre = lo + sm + w / 2;
        let k = anchors.partition_point(|a| a.0 < centre);
        let near = if k == 0 {
            anchors[0]
        } else if k == anchors.len() || centre - anchors[k - 1].0 <= anchors[k].0 - centre {
            anchors[k - 1]
        } else {
            anchors[k]
        };
        let ps = near.1 + (lo + sm - near.0);
        let clock = std::time::Instant::now();
        let seed = edit_seed_band(xw, y, ps, self.band, self.edit_gap);
        if prof_on() {
            prof_add(0, clock.elapsed().as_nanos() as u64);
            prof_add(1, 1);
            prof_add(2, seed.is_some() as u64);
        }
        let j0 = seed.map_or(ps.clamp(0, ylen - 1), |(_, st)| st);
        self.stats[4].fetch_add(n as u64, Relaxed);
        // first in half the band; a pair that fails is aligned again in the full band
        let half = (self.path / 2).max(1);
        let thr = self.thr as f64;
        PBUF.with(|pb| {
            let pb = &mut *pb.borrow_mut();
            let prof = prof_on();
            let clock = std::time::Instant::now();
            let first = path_attempt(x, y, lo as usize, hi as usize, sm as usize, j0 as usize, thr, half, self.path_shift,
                                     self.edit_gap, pb);
            if prof {
                prof_add(3, clock.elapsed().as_nanos() as u64);
            }
            if first {
                return true;
            }
            if half == self.path {
                return false;
            }
            self.stats[5].fetch_add(1, Relaxed);
            let clock = std::time::Instant::now();
            let second = path_attempt(x, y, lo as usize, hi as usize, sm as usize, j0 as usize, thr, self.path,
                                      self.path_shift, self.edit_gap, pb);
            if prof {
                prof_add(4, clock.elapsed().as_nanos() as u64);
            }
            second
        })
    }
}

/// The check's per-thread buffers: traceback directions, each row's band start, the path's cost per overlap row, two
/// DP rows, the row's vertical / diagonal candidates, and reversed copies of the reads for the backward half.
#[derive(Default)]
struct PathScratch {
    dir: Vec<u8>,
    lefts: Vec<i64>,
    rowcost: Vec<f64>,
    pc: Vec<i32>,
    cc: Vec<i32>,
    vd: Vec<i32>,
    xr: Vec<u8>,
    yr: Vec<u8>,
}

thread_local! {
    static PBUF: std::cell::RefCell<PathScratch> = Default::default();
}

/// After a pair, this thread's alignment buffers grown past `PATH_KEEP` bytes by a long overlap are freed: the
/// direction matrix is rows x band bytes (~200 MB for a 500,000-symbol overlap), so threads that kept the largest
/// buffer they had needed would hold several GB between long pairs.
const PATH_KEEP: usize = 4 << 20;

/// The bytes this thread's alignment buffers hold (`memtrack`).
fn path_buffer_bytes() -> usize {
    fn b<T>(v: &Vec<T>) -> usize {
        v.capacity() * std::mem::size_of::<T>()
    }
    let p = PBUF.with(|pb| {
        let pb = &*pb.borrow();
        b(&pb.dir) + b(&pb.lefts) + b(&pb.rowcost) + b(&pb.pc) + b(&pb.cc) + b(&pb.vd) + b(&pb.xr) + b(&pb.yr)
    });
    p + UNPACKED.with(|cell| {
        let (x, y) = &*cell.borrow();
        b(x) + b(y)
    }) + RAW_BYTES.with(|cell| b(&cell.borrow()))
}

fn trim_path_buffers() {
    fn keep<T>(v: &mut Vec<T>) {
        if v.capacity() * std::mem::size_of::<T>() > PATH_KEEP {
            *v = Vec::new();
        }
    }
    let track = super::memtrack::on();
    if track {
        super::memtrack::peak(super::memtrack::PATH, path_buffer_bytes());
    }
    PBUF.with(|pb| {
        let pb = &mut *pb.borrow_mut();
        keep(&mut pb.dir);
        keep(&mut pb.lefts);
        keep(&mut pb.rowcost);
        keep(&mut pb.pc);
        keep(&mut pb.cc);
        keep(&mut pb.vd);
        keep(&mut pb.xr);
        keep(&mut pb.yr);
    });
    UNPACKED.with(|cell| {
        let (x, y) = &mut *cell.borrow_mut();
        keep(x);
        keep(y);
    });
    RAW_BYTES.with(|cell| keep(&mut cell.borrow_mut()));
    if track {
        super::memtrack::set(super::memtrack::PATH_KEPT, path_buffer_bytes());
    }
}

pub(super) mod simd;

pub(super) mod simd16;

/// A cost no reachable cell reaches (the DP rows' unreachable cells).
const INF: i32 = i32::MAX / 4;

/// The substitution cost of two q-symbols: their level difference.
#[inline]
fn cost(a: u8, b: u8) -> i32 {
    (a as i32 - b as i32).abs()
}

/// `tracked`'s second pass by `simd::hscan_avx2` (values, directions, first cheapest cell); None: the scalar loop.
fn hscan(v: &[i32], gap: i32, out: &mut [i32], r: &mut [u8]) -> Option<(i32, usize)> {
    #[cfg(target_arch = "x86_64")]
    {
        if simd::on() {
            // SAFETY: simd::on() checked AVX2; out and r are v's row (asserted inside)
            return Some(unsafe { simd::hscan_avx2(v, gap, INF, out, r) });
        }
    }
    let _ = (v, gap, out, r);
    None
}

/// `tracked` with 16-bit lanes (`simd16`): Some(its result, the same end and path as the 32-bit kernel) when it ran;
/// None when it is off (no AVX2, a gap above 64) or a row's values left the 16-bit range (the caller then runs the
/// 32-bit kernel).
#[allow(clippy::too_many_arguments)]
fn tracked16(x: &[u8], y: &[u8], band: i64, shift: i64, gap: i32, dir: &mut Vec<u8>, lefts: &mut Vec<i64>)
             -> Option<Option<i64>> {
    #[cfg(target_arch = "x86_64")]
    {
        if simd16::on() && (0..=64).contains(&gap) {
            // SAFETY: simd16::on() checked AVX2
            return unsafe { simd16::tracked16_avx2(x, y, band, shift, gap as i16, dir, lefts, simd16::LIM, simd16::REBASE_AT) };
        }
    }
    let _ = (x, y, band, shift, gap, dir, lefts);
    None
}

/// The seed window `xw` placed on `y` by a semi-global weighted edit distance (all of `xw`, free start and end in
/// `y`; ties diagonal, then up, then left) in a diagonal band of +- `band` columns around the expected start `ps`
/// (row i's columns ps + i - band ..= ps + i + band; free start in row 0's band, free end in the last row's):
/// (cost, start in y). None when the band leaves `y` or the best path starts or ends at the band's edge.
fn edit_seed_band(xw: &[u8], y: &[u8], ps: i64, band: i64, gap: i32) -> Option<(f32, i64)> {
    #[cfg(target_arch = "x86_64")]
    {
        if simd::on() {
            // SAFETY: simd::on() checked AVX2
            return unsafe { edit_seed_band_avx2(xw, y, ps, band, gap) };
        }
    }
    edit_seed_band_scalar(xw, y, ps, band, gap)
}

/// `edit_seed_band`'s scalar reference.
fn edit_seed_band_scalar(xw: &[u8], y: &[u8], ps: i64, band: i64, gap: i32) -> Option<(f32, i64)> {
    const INF: i32 = i32::MAX / 4;
    let (n, m) = (xw.len(), y.len() as i64);
    let width = (2 * band + 1) as usize;
    if n == 0 || ps - band < 0 || ps + n as i64 - 1 + band >= m {
        return None;
    }
    let sub = |a: u8, b: u8| (a as i32 - b as i32).abs();
    // cell k of row i is column ps + i - band + k; diagonal = same k of row i - 1, up = k + 1, left = k - 1
    let mut pc: Vec<i32> = (0..width).map(|k| sub(xw[0], y[(ps - band + k as i64) as usize])).collect();
    let mut ps_k: Vec<u32> = (0..width as u32).collect();
    let mut cc = vec![INF; width];
    let mut cs = vec![0u32; width];
    for i in 1..n {
        let base = ps + i as i64 - band;
        for k in 0..width {
            let c = sub(xw[i], y[(base + k as i64) as usize]);
            let (mut b, mut st) = (pc[k] + c, ps_k[k]);
            if k + 1 < width && pc[k + 1] + gap < b {
                b = pc[k + 1] + gap;
                st = ps_k[k + 1];
            }
            if k > 0 && cc[k - 1] + gap < b {
                b = cc[k - 1] + gap;
                st = cs[k - 1];
            }
            cc[k] = b;
            cs[k] = st;
        }
        std::mem::swap(&mut pc, &mut cc);
        std::mem::swap(&mut ps_k, &mut cs);
    }
    let (mut bk, mut bv) = (0usize, INF);
    for (k, &v) in pc.iter().enumerate() {
        if v < bv {
            bv = v;
            bk = k;
        }
    }
    let st = ps_k[bk] as usize;
    if bk == 0 || bk == width - 1 || st == 0 || st == width - 1 {
        return None;
    }
    Some((bv as f32, ps - band + st as i64))
}

/// `edit_seed_band_scalar` with AVX2 (two passes per row: the diagonal / vertical candidates vectorized, then the
/// horizontal steps by `simd::hscan_st`): the same values and starts.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn edit_seed_band_avx2(xw: &[u8], y: &[u8], ps: i64, band: i64, gap: i32) -> Option<(f32, i64)> {
    const INF: i32 = i32::MAX / 4;
    let (n, m) = (xw.len(), y.len() as i64);
    let width = (2 * band + 1) as usize;
    if n == 0 || ps - band < 0 || ps + n as i64 - 1 + band >= m {
        return None;
    }
    let sub = |a: u8, b: u8| (a as i32 - b as i32).abs();
    let mut pc: Vec<i32> = (0..width).map(|k| sub(xw[0], y[(ps - band + k as i64) as usize])).collect();
    let mut ps_k: Vec<u32> = (0..width as u32).collect();
    let mut cc = vec![INF; width];
    let mut cs = vec![0u32; width];
    let mut vv = vec![0i32; width];
    let mut sv = vec![0u32; width];
    for i in 1..n {
        let base = (ps + i as i64 - band) as usize;
        let yrow = &y[base..base + width];
        let xi = xw[i];
        // cell k: diagonal from pc[k], vertical from pc[k + 1]; the last cell has no vertical step (its up
        // neighbour is outside the band)
        for ((((v, s), (&pd, &pu)), (&sd, &su)), &yk) in vv[..width - 1].iter_mut().zip(sv[..width - 1].iter_mut())
            .zip(pc[..width - 1].iter().zip(&pc[1..])).zip(ps_k[..width - 1].iter().zip(&ps_k[1..]))
            .zip(&yrow[..width - 1]) {
            let dv = pd + (xi as i32 - yk as i32).abs();
            let uv = pu + gap;
            let up = uv < dv;
            *v = if up { uv } else { dv };
            *s = if up { su } else { sd };
        }
        vv[width - 1] = pc[width - 1] + sub(xi, yrow[width - 1]);
        sv[width - 1] = ps_k[width - 1];
        cc[0] = vv[0];
        cs[0] = sv[0];
        simd::hscan_st(&vv[1..], &sv[1..], gap, vv[0], sv[0], &mut cc[1..], &mut cs[1..]);
        std::mem::swap(&mut pc, &mut cc);
        std::mem::swap(&mut ps_k, &mut cs);
    }
    let (mut bk, mut bv) = (0usize, INF);
    for (k, &v) in pc.iter().enumerate() {
        if v < bv {
            bv = v;
            bk = k;
        }
    }
    let st = ps_k[bk] as usize;
    if bk == 0 || bk == width - 1 || st == 0 || st == width - 1 {
        return None;
    }
    Some((bv as f32, ps - band + st as i64))
}

/// The weighted edit distance of all of `x` against `y` from the anchored cell (0, 0), in a band of +- `band`
/// columns that follows the path (as nanopolish's adaptive banded alignment): row i is centred one column after row
/// i - 1's cheapest cell, moved at most `shift` from the previous centre + 1. Diagonal = the substitution cost, up /
/// left = `gap`; ties diagonal, then up, then left. Free end: returns the last row's cheapest column, the traceback
/// in `dir` / `lefts` (`trace`); None when a row has no reachable cell.
#[allow(clippy::too_many_arguments)]
fn tracked(x: &[u8], y: &[u8], band: i64, shift: i64, gap: i32, dir: &mut Vec<u8>, lefts: &mut Vec<i64>,
           pc: &mut Vec<i32>, cc: &mut Vec<i32>, vd: &mut Vec<i32>) -> Option<i64> {
    #[cfg(target_arch = "x86_64")]
    {
        if simd::on() {
            if let Some(end) = tracked16(x, y, band, shift, gap, dir, lefts) {
                return end;
            }
            // SAFETY: simd::on() checked AVX2
            return unsafe { tracked_avx2(x, y, band, shift, gap, dir, lefts, pc, cc, vd) };
        }
    }
    tracked_body(x, y, band, shift, gap, dir, lefts, pc, cc, vd, false)
}

/// `tracked_body` compiled with AVX2 (its first pass vectorized 8 wide) and the AVX2 horizontal scan.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
unsafe fn tracked_avx2(x: &[u8], y: &[u8], band: i64, shift: i64, gap: i32, dir: &mut Vec<u8>, lefts: &mut Vec<i64>,
                       pc: &mut Vec<i32>, cc: &mut Vec<i32>, vd: &mut Vec<i32>) -> Option<i64> {
    tracked_body(x, y, band, shift, gap, dir, lefts, pc, cc, vd, true)
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn tracked_body(x: &[u8], y: &[u8], band: i64, shift: i64, gap: i32, dir: &mut Vec<u8>, lefts: &mut Vec<i64>,
                pc: &mut Vec<i32>, cc: &mut Vec<i32>, vd: &mut Vec<i32>, fast: bool) -> Option<i64> {
    let n = x.len();
    let m = y.len() as i64;
    let width = (2 * band + 1) as usize;
    // DP rows padded with INF so the predecessors of every band cell can be read without a test
    let pad = (shift + 2) as usize;
    dir.clear();
    dir.resize(n * width, 255);
    lefts.clear();
    pc.clear();
    pc.resize(width + 2 * pad, INF);
    cc.clear();
    cc.resize(width + 2 * pad, INF);
    vd.clear();
    vd.resize(width, 0);
    let mut centre = 0i64;
    let mut pl = -band;
    lefts.push(pl);
    pc[pad + band as usize] = cost(x[0], y[0]);
    dir[band as usize] = 3;
    let mut best_k = band as usize;
    for i in 1..n {
        let want = pl + best_k as i64 + 1;
        centre = want.clamp(centre + 1 - shift, centre + 1 + shift);
        let l = centre - band;
        lefts.push(l);
        let klo = (-l).max(0);
        let khi = (m - 1 - l).min(width as i64 - 1);
        if khi < klo {
            return None;
        }
        let (klo, khi) = (klo as usize, khi as usize);
        for v in cc[..pad + klo].iter_mut() {
            *v = INF;
        }
        for v in cc[pad + khi + 1..].iter_mut() {
            *v = INF;
        }
        // column l + k of row i is column (l + k) - pl of row i - 1. Two passes: the vertical / diagonal
        // candidates have no dependency along the row (vectorizable), then the horizontal steps (each costing the
        // gap) in order.
        let off = (pad as i64 + l - pl) as usize;
        let xi = x[i];
        let yrow = &y[(l + klo as i64) as usize..=(l + khi as i64) as usize];
        let len = yrow.len();
        let row = &mut dir[i * width + klo..i * width + klo + len];
        let up = &pc[klo + off..klo + off + len];
        let dg = &pc[klo + off - 1..klo + off - 1 + len];
        let vd = &mut vd[..len];
        for ((((v, r), &u), &g), &yj) in vd.iter_mut().zip(row.iter_mut()).zip(up).zip(dg).zip(yrow) {
            let dv = g + cost(xi, yj);
            let uv = u + gap;
            let d = !(uv < dv);
            *v = if d { dv } else { uv };
            *r = d as u8;
        }
        let (mut best_v, mut bk) = (INF, usize::MAX);
        let out = &mut cc[pad + klo..pad + klo + len];
        match if fast { hscan(vd, gap, out, row) } else { None } {
            Some((v, t)) => {
                best_v = v;
                bk = if t == usize::MAX { t } else { klo + t };
            }
            None => {
                let mut left = INF;
                for (t, ((o, &v), r)) in out.iter_mut().zip(vd.iter()).zip(row.iter_mut()).enumerate() {
                    let h = gap + left;
                    let val = if h < v {
                        *r = 2;
                        h
                    } else {
                        v
                    };
                    *o = val;
                    left = val;
                    if val < best_v {
                        best_v = val;
                        bk = klo + t;
                    }
                }
            }
        }
        // no cell of the row is reachable (every cell >= INF)
        if !(best_v < INF) {
            return None;
        }
        best_k = bk;
        std::mem::swap(pc, cc);
        pl = l;
    }
    Some(pl + best_k as i64)
}

/// Walks `tracked`'s path back from the last row's column `end`: each cell is one alignment column, costing the
/// substitution when reached diagonally (or the start cell), the gap when reached from above or the left, added to
/// its row.
fn trace(x: &[u8], y: &[u8], dir: &[u8], lefts: &[i64], band: i64, end: i64, gap: i32, mut add: impl FnMut(usize, f64)) {
    let width = (2 * band + 1) as usize;
    let (mut i, mut j) = (x.len() - 1, end);
    loop {
        let d = dir[i * width + (j - lefts[i]) as usize];
        let c = if d == 0 || d == 2 { gap } else { cost(x[i], y[j as usize]) };
        add(i, c as f64);
        match d {
            3 => break,
            0 => i -= 1,
            1 => {
                i -= 1;
                j -= 1;
            }
            2 => j -= 1,
            _ => unreachable!("path off the band"),
        }
    }
}

/// One path of the overlap's rows [lo, hi) of `x` through the seed cell (row sm of the overlap, column j0 of `y`),
/// forward then backward, in a band of +- `band`; true when the whole path costs at most `thr` per 1,000 alignment
/// columns (path cells).
#[allow(clippy::too_many_arguments)]
fn path_attempt(x: &[u8], y: &[u8], lo: usize, hi: usize, sm: usize, j0: usize, thr: f64, band: i64, shift: i64,
                gap: i32, pb: &mut PathScratch) -> bool {
    let PathScratch { dir, lefts, rowcost, pc, cc, vd, xr, yr } = pb;
    let n = hi - lo;
    rowcost.clear();
    rowcost.resize(n, 0.0);
    // alignment columns (path cells) so far
    let mut cols = 0usize;
    let prof = prof_on();
    let width = (2 * band + 1) as u64;
    // forward from the seed's start cell to the overlap's end
    let xs = &x[lo + sm..hi];
    let ys = &y[j0..];
    let clock = std::time::Instant::now();
    let end = tracked(xs, ys, band, shift, gap, dir, lefts, pc, cc, vd);
    if prof {
        prof_add(5, clock.elapsed().as_nanos() as u64);
        prof_add(6, xs.len() as u64 * width);
        prof_add(8, 1);
    }
    let end = match end {
        None => return false,
        Some(e) => e,
    };
    let clock = std::time::Instant::now();
    trace(xs, ys, dir, lefts, band, end, gap, |i, c| {
        rowcost[sm + i] += c;
        cols += 1;
    });
    if prof {
        prof_add(7, clock.elapsed().as_nanos() as u64);
    }
    // backward from the same cell to the overlap's start
    if sm > 0 {
        xr.clear();
        xr.extend(x[lo..=lo + sm].iter().rev());
        yr.clear();
        yr.extend(y[..=j0].iter().rev());
        let clock = std::time::Instant::now();
        let end = tracked(xr, yr, band, shift, gap, dir, lefts, pc, cc, vd);
        if prof {
            prof_add(5, clock.elapsed().as_nanos() as u64);
            prof_add(6, xr.len() as u64 * width);
            prof_add(8, 1);
        }
        let end = match end {
            None => return false,
            Some(e) => e,
        };
        let clock = std::time::Instant::now();
        trace(xr, yr, dir, lefts, band, end, gap, |i, c| {
            if i > 0 {
                rowcost[sm - i] += c;
                cols += 1;
            }
        });
        if prof {
            prof_add(7, clock.elapsed().as_nanos() as u64);
        }
    }
    let total: f64 = rowcost.iter().sum();
    total * 1000.0 <= thr * cols as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The q-symbols read back from the local file (`SymStore::File`) are the in-memory ones (empty reads, every length
    /// modulo 8).
    #[test]
    fn qsymbol_file_store_reads_the_same_symbols() {
        let mut state = 3u64;
        let mut next = move || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as u8
        };
        let mut syms: Vec<Syms> = (0..40).map(|k| Syms::new((0..k * 7 % 53).map(|_| next() % 8).collect()).unwrap()).collect();
        syms.push(Syms::default());   // a read without q-symbols
        let store = SymStore::to_file(&syms, &std::env::temp_dir()).unwrap();
        let memory = SymStore::Memory(syms.clone());
        let (mut raw, mut x, mut y) = (Vec::new(), Vec::new(), Vec::new());
        for i in (0..syms.len()).rev() {
            store.unpacked(i, &mut raw, &mut x);
            memory.unpacked(i, &mut raw, &mut y);
            assert_eq!(x, y);
            assert_eq!(x.len(), syms[i].len());
        }
    }


    fn xs(s: &mut u64) -> u64 {
        *s ^= *s << 13;
        *s ^= *s >> 7;
        *s ^= *s << 17;
        *s
    }

    /// q-symbols: a random sequence and a noisy copy of it (substitutions, insertions, deletions), or an
    /// unrelated one
    fn qpair(s: &mut u64, n: usize, related: bool) -> (Vec<u8>, Vec<u8>) {
        let x: Vec<u8> = (0..n).map(|_| (xs(s) % 8) as u8).collect();
        if !related {
            return (x, (0..n + n / 10).map(|_| (xs(s) % 8) as u8).collect());
        }
        let mut y = Vec::with_capacity(n + n / 5);
        for &a in &x {
            match xs(s) % 20 {
                0 => {}
                1 => {
                    y.push(a);
                    y.push((xs(s) % 8) as u8);
                }
                2 | 3 => y.push((a as i64 + (xs(s) % 3) as i64 - 1).clamp(0, 7) as u8),
                _ => y.push(a),
            }
        }
        (x, y)
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn simd_scans_equal_the_scalar_loops() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let inf = i32::MAX / 4;
        let mut s = 0x9E37_79B9_7F4A_7C15u64;
        for case in 0..100_000 {
            let n = if case % 3 == 0 { (xs(&mut s) % 40) as usize } else { (xs(&mut s) % 3100) as usize };
            let gap = 1 + (xs(&mut s) % 3) as i32;
            let v: Vec<i32> = (0..n).map(|_| match case % 4 {
                0 => (xs(&mut s) % 40) as i32,
                1 => 1000 + (xs(&mut s) % 9) as i32,
                2 => if xs(&mut s) % 5 == 0 { inf + (xs(&mut s) % 30) as i32 } else { (xs(&mut s) % 20) as i32 },
                _ => inf + (xs(&mut s) % 30) as i32,
            }).collect();
            let r0: Vec<u8> = (0..n).map(|_| (xs(&mut s) % 2) as u8).collect();
            let (mut o1, mut o2, mut r1, mut r2) = (vec![0; n], vec![0; n], r0.clone(), r0);
            let a = simd::hscan_scalar(&v, gap, inf, &mut o1, &mut r1);
            let b = unsafe { simd::hscan_avx2(&v, gap, inf, &mut o2, &mut r2) };
            assert_eq!((o1, r1, a), (o2, r2, b), "hscan, case {}", case);
            let sv: Vec<u32> = (0..n).map(|_| (xs(&mut s) % 5000) as u32).collect();
            let (left0, st0) = ((xs(&mut s) % 50) as i32, (xs(&mut s) % 5000) as u32);
            let (mut o1, mut o2, mut s1, mut s2) = (vec![0; n], vec![0; n], vec![0u32; n], vec![0u32; n]);
            simd::hscan_st_scalar(&v, &sv, gap, left0, st0, &mut o1, &mut s1);
            unsafe { simd::hscan_st_avx2(&v, &sv, gap, left0, st0, &mut o2, &mut s2) };
            assert_eq!((o1, s1), (o2, s2), "hscan_st, case {}", case);
        }
    }


    /// Timing of the scalar and AVX2 `tracked` (cargo test --release tracked_speed -- --ignored --nocapture).
    #[cfg(target_arch = "x86_64")]
    #[test]
    #[ignore]
    fn tracked_speed() {
        let mut s = 0x5555_1234_9876_ABCDu64;
        for (related, band) in [(true, 100i64), (false, 100), (false, 200)] {
            let (x, y) = qpair(&mut s, 8000, related);
            let cells = 8000.0 * (2 * band + 1) as f64 * 10.0;
            let (mut d, mut l, mut pc, mut cc, mut v) = Default::default();
            let t = std::time::Instant::now();
            for _ in 0..10 {
                std::hint::black_box(tracked_body(&x, &y, band, 1, 2, &mut d, &mut l, &mut pc, &mut cc, &mut v, false));
            }
            let a = t.elapsed().as_secs_f64() * 1e9 / cells;
            let t = std::time::Instant::now();
            for _ in 0..10 {
                std::hint::black_box(unsafe { tracked_avx2(&x, &y, band, 1, 2, &mut d, &mut l, &mut pc, &mut cc, &mut v) });
            }
            eprintln!("tracked related {} band {}: scalar {:.2}, avx2 {:.2} ns per cell", related, band, a,
                      t.elapsed().as_secs_f64() * 1e9 / cells);
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn simd_alignments_equal_the_scalar_ones() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let mut s = 0x0123_4567_89AB_CDEFu64;
        for case in 0..300 {
            let n = 200 + (xs(&mut s) % 1800) as usize;
            let (x, y) = qpair(&mut s, n, case % 3 != 0);
            let gap = 1 + (case % 3) as i32;
            // the seed search
            let w = (n / 4).max(20);
            let xw = &x[n / 3..n / 3 + w];
            for band in [5i64, 50, 100] {
                let ps = n as i64 / 3 + (xs(&mut s) % 21) as i64 - 10;
                assert_eq!(edit_seed_band_scalar(xw, &y, ps, band, gap),
                           unsafe { edit_seed_band_avx2(xw, &y, ps, band, gap) }, "edit_seed_band, case {}", case);
            }
            // the path: values, end and every direction of the band
            for band in [20i64, 100, 200] {
                let (mut d1, mut l1, mut p1, mut c1, mut v1) = Default::default();
                let (mut d2, mut l2, mut p2, mut c2, mut v2) = Default::default();
                let e1 = tracked_body(&x, &y, band, 1, gap, &mut d1, &mut l1, &mut p1, &mut c1, &mut v1, false);
                let e2 = unsafe { tracked_avx2(&x, &y, band, 1, gap, &mut d2, &mut l2, &mut p2, &mut c2, &mut v2) };
                assert_eq!(e1, e2, "tracked end, case {} band {}", case, band);
                assert_eq!(l1, l2, "tracked lefts, case {} band {}", case, band);
                assert!(d1 == d2, "tracked directions, case {} band {}", case, band);
            }
        }
    }

    fn wave(n: usize, seed: u64) -> Vec<f32> {
        let mut v = Vec::with_capacity(n);
        let mut s = seed;
        for _ in 0..n {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            v.push(((s >> 33) % 1000) as f32 / 250.0 - 2.0);
        }
        v
    }

    fn quantize(v: &[f32]) -> Vec<u8> {
        v.iter().map(|&a| ((a + 2.0) * 1.75).floor().clamp(0.0, 7.0) as u8).collect()
    }

    /// A check of two reads' q-symbols (window 100, seed band 20, path band 20, gap 2).
    fn qcheck(x: &[f32], y: &[f32], thr: f32) -> Check {
        let mut c = Check { symbols: std::sync::RwLock::new(SymStore::Memory(Vec::new())), window: 100, thr, band: 20,
                            path: 20, path_shift: 1, edit_gap: 2,
                            shards: (0..SHARDS).map(|_| Mutex::new(HashMap::new())).collect(), loaded: 0,
                            stats: Default::default() };
        c.set_symbols(vec![Syms::new(quantize(x)).unwrap(), Syms::new(quantize(y)).unwrap()]);
        c
    }

    /// The donor of the verdict tests: the target with 5% of its symbols doubled (another event rate) and small
    /// noise; `diagonal` gives a target position's place in it (the chain anchors).
    fn doubled(x: &[f32]) -> Vec<f32> {
        let noise = wave(x.len() + x.len() / 5, 11);
        let mut y = Vec::new();
        for (i, &v) in x.iter().enumerate() {
            y.push(v + noise[i] * 0.02);
            if i % 20 == 0 {
                y.push(v + noise[i + 1] * 0.02);
            }
        }
        y
    }

    fn diagonal(i: i64) -> i64 {
        i + (i + 19) / 20
    }

    fn anchors(n: i64) -> Vec<(i64, i64)> {
        (0..n).step_by(100).map(|i| (i, diagonal(i))).collect()
    }

    /// 3-bit q-symbols unpack to the input at every length (partial last group included).
    #[test]
    fn packed_q_symbols_round_trip() {
        let mut s: u64 = 7;
        for n in [0usize, 1, 2, 7, 8, 9, 15, 16, 17, 1000, 1003] {
            let v: Vec<u8> = (0..n)
                .map(|_| {
                    s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    (s >> 61) as u8
                })
                .collect();
            let packed = Syms::new(v.clone()).unwrap();
            assert_eq!(packed.len(), n);
            assert!(packed.data.len() <= (3 * n + 7) / 8 + 1);
            let mut out = vec![9u8; 3];
            packed.unpack_into(&mut out);
            assert_eq!(out, v);
        }
        assert!(Syms::new(vec![3, 8]).is_err());
    }

    #[test]
    fn weighted_edit_counts_columns() {
        // x = y with 3 extra symbols in y and 2 substitutions of one level: cost 3 * gap + 2
        let x: Vec<u8> = (0..300).map(|i| ((i * 5 + i / 7) % 8) as u8).collect();
        let mut y = x.clone();
        y[50] = if y[50] == 7 { 6 } else { y[50] + 1 };
        y[200] = if y[200] == 0 { 1 } else { y[200] - 1 };
        y.insert(120, 3);
        y.insert(121, 3);
        y.insert(250, 6);
        let (mut dir, mut lefts, mut pc, mut cc, mut vd) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let end = tracked(&x, &y, 20, 1, 2, &mut dir, &mut lefts, &mut pc, &mut cc, &mut vd).unwrap();
        let (mut total, mut cols) = (0.0, 0);
        trace(&x, &y, &dir, &lefts, 20, end, 2, |_, c| {
            total += c;
            cols += 1;
        });
        assert!(total <= 8.0, "cost {}", total);
        assert_eq!(cols, 303);
        // the banded seed search finds x inside a longer y from a nearby guess; from a far one it finds a costly
        // placement or gives up at the band's edge
        let mut yy: Vec<u8> = (0..100).map(|i| ((i * 3) % 8) as u8).collect();
        yy.extend(&x);
        assert_eq!(edit_seed_band(&x[..100], &yy, 95, 20, 2), Some((0.0, 100)));
        assert!(edit_seed_band(&x[..100], &yy, 40, 20, 2).map_or(true, |(c, _)| c > 30.0));
        assert_eq!(edit_seed_band(&x[..100], &yy, 10, 20, 2), None);
    }

    #[test]
    fn whole_overlap_rule_and_directed_verdicts() {
        let x = wave(1000, 7);
        let y = doubled(&x);
        let line = (100, diagonal(100), 900, diagonal(900));
        let c = qcheck(&x, &y, 300.0);
        assert!(c.verdict(0, 1, line, 0, 1000, &anchors(1000)));
        assert!(c.verdict(0, 1, line, 0, 1000, &anchors(1000)));   // reused
        assert_eq!((c.stats[0].load(Relaxed), c.stats[2].load(Relaxed)), (1, 1));
        assert_eq!(c.cached(0, 1), Some(true));
        // directed: the reverse pair is computed on its own
        assert_eq!(c.cached(1, 0), None);
        let reverse: Vec<(i64, i64)> = anchors(1000).iter().map(|&(t, d)| (d, t)).collect();
        c.verdict(1, 0, (line.1, line.0, line.3, line.2), 0, y.len() as i64, &reverse);
        assert_eq!(c.stats[0].load(Relaxed), 2);
        // an unrelated stretch of 300 in the middle: the whole overlap costs more
        let junk = wave(300, 99);
        let mut y2 = y.clone();
        y2[400..700].copy_from_slice(&junk);
        assert!(!qcheck(&x, &y2, 300.0).verdict(0, 1, line, 0, 1000, &anchors(1000)));
        // the same pair passes a looser threshold
        assert!(qcheck(&x, &y2, 2000.0).verdict(0, 1, line, 0, 1000, &anchors(1000)));
    }

    #[test]
    fn the_seed_follows_the_chain_anchor() {
        // y = 300 unrelated symbols, then x with 5% of its symbols doubled; the chain's middle anchor places the seed
        let x = wave(3000, 7);
        let mut y = wave(300, 5);
        y.extend(doubled(&x));
        let true_anchor = |i: i64| (i, 300 + diagonal(i));
        let anchors: Vec<(i64, i64)> = [100, 800, 1500, 2200, 2900].iter().map(|&i| true_anchor(i)).collect();
        let verdict = |y: &[f32], anchors: &[(i64, i64)]| {
            qcheck(&x, y, 300.0).verdict(0, 1, (100, 100, 2900, 2900), 0, 3000, anchors)
        };
        assert!(verdict(&y, &anchors));
        // an unrelated donor fails
        let junk = wave(3500, 99);
        assert!(!verdict(&junk, &anchors));
    }

    /// xorshift64* for the randomized 16-bit tests
    fn rng(s: &mut u64) -> u64 {
        *s ^= *s >> 12;
        *s ^= *s << 25;
        *s ^= *s >> 27;
        s.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    /// the cells of `tracked`'s path from its end, as `trace` walks them
    fn path_cells(n: usize, dir: &[u8], lefts: &[i64], band: i64, end: i64) -> Vec<(usize, i64, u8)> {
        let width = (2 * band + 1) as usize;
        let (mut i, mut j) = (n - 1, end);
        let mut out = Vec::new();
        loop {
            let d = dir[i * width + (j - lefts[i]) as usize];
            out.push((i, j, d));
            match d {
                3 => break,
                0 => i -= 1,
                1 => {
                    i -= 1;
                    j -= 1;
                }
                2 => j -= 1,
                _ => panic!("path off the band"),
            }
        }
        out
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn row16_avx2_matches_scalar() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let mut s = 0x9e37_79b9_7f4a_7c15u64;
        for _ in 0..4000 {
            let n = (rng(&mut s) % 140) as usize;
            let gap = 1 + (rng(&mut s) % 4) as i16;
            let xi = (rng(&mut s) % 8) as i16;
            let val = |s: &mut u64| if rng(s) % 6 == 0 { simd16::INF } else { (rng(s) % 30000) as i16 };
            let y: Vec<u8> = (0..n).map(|_| (rng(&mut s) % 8) as u8).collect();
            let up: Vec<i16> = (0..n).map(|_| val(&mut s)).collect();
            let dg: Vec<i16> = (0..n).map(|_| val(&mut s)).collect();
            let (mut o1, mut o2) = (vec![0i16; n], vec![0i16; n]);
            let (mut r1, mut r2) = (vec![9u8; n], vec![9u8; n]);
            let a = simd16::row16_scalar(xi, &y, &up, &dg, gap, &mut o1, &mut r1);
            let b = unsafe { simd16::row16_avx2(xi, &y, &up, &dg, gap, &mut o2, &mut r2) };
            assert_eq!((a, &o1, &r1), (b, &o2, &r2));
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn tracked16_matches_the_32bit_kernel() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let mut s = 0x1234_5678_9abc_def1u64;
        let (mut fell, mut ran) = (0, 0);
        for case in 0..600 {
            let n = 20 + (rng(&mut s) % 1500) as usize;
            let x: Vec<u8> = (0..n).map(|_| (rng(&mut s) % 8) as u8).collect();
            // y: x with substitutions, insertions and deletions (or unrelated), then a tail
            let mut y = Vec::new();
            let related = case % 5 != 0;
            for &v in &x {
                let r = rng(&mut s) % 100;
                if !related || r < 8 {
                    y.push((rng(&mut s) % 8) as u8);
                } else if r < 12 {
                    y.push(v);
                    y.push(v);
                } else if r >= 16 {
                    y.push(v);
                }
            }
            y.extend((0..(rng(&mut s) % 50)).map(|_| (rng(&mut s) % 8) as u8));
            if y.is_empty() {
                continue;
            }
            let band = 3 + (rng(&mut s) % 40) as i64;
            let shift = 1 + (rng(&mut s) % 3) as i64;
            let gap = 1 + (rng(&mut s) % 3) as i32;
            let (mut d32, mut l32) = (Vec::new(), Vec::new());
            let (mut pc, mut cc, mut vd) = (Vec::new(), Vec::new(), Vec::new());
            let e32 = tracked_body(&x, &y, band, shift, gap, &mut d32, &mut l32, &mut pc, &mut cc, &mut vd, false);
            // the production limits, rebasing on almost every row, and a tiny limit that makes many pairs fall back
            for (lim, rebase_at) in [(simd16::LIM, simd16::REBASE_AT), (simd16::LIM, 20), (12, 20)] {
                let (mut d16, mut l16) = (Vec::new(), Vec::new());
                match unsafe { simd16::tracked16_avx2(&x, &y, band, shift, gap as i16, &mut d16, &mut l16, lim, rebase_at) } {
                    None => fell += 1,
                    Some(e16) => {
                        ran += 1;
                        assert_eq!(e16, e32, "case {} lim {} rebase {}", case, lim, rebase_at);
                        if let Some(e) = e32 {
                            assert_eq!(l16, l32);
                            assert_eq!(path_cells(n, &d16, &l16, band, e), path_cells(n, &d32, &l32, band, e));
                        }
                    }
                }
            }
        }
        assert!(fell > 0 && ran >= 1200, "fell back {} ran {}", fell, ran);
    }

    /// Timing of the 32- and 16-bit `tracked` kernels on one related pair (cargo test --release tracked_kernel_speed --
    /// --ignored --nocapture).
    #[cfg(target_arch = "x86_64")]
    #[test]
    #[ignore]
    fn tracked_kernel_speed() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let mut s = 0x5151_2626_3737_4848u64;
        let n = 20_000usize;
        let x: Vec<u8> = (0..n).map(|_| (rng(&mut s) % 8) as u8).collect();
        let mut y = Vec::new();
        for &v in &x {
            let r = rng(&mut s) % 100;
            if r < 20 {
                y.push((rng(&mut s) % 8) as u8);
            } else if r < 25 {
                y.push(v);
                y.push(v);
            } else if r >= 30 {
                y.push(v);
            }
        }
        for band in [100i64, 200] {
            let (mut d, mut l) = (Vec::new(), Vec::new());
            let (mut pc, mut cc, mut vd) = (Vec::new(), Vec::new(), Vec::new());
            let reps = 20;
            let t = std::time::Instant::now();
            let mut e32 = None;
            for _ in 0..reps {
                e32 = unsafe { tracked_avx2(&x, &y, band, 1, 2, &mut d, &mut l, &mut pc, &mut cc, &mut vd) };
            }
            let t32 = t.elapsed().as_secs_f64() / reps as f64;
            let t = std::time::Instant::now();
            let mut e16 = None;
            for _ in 0..reps {
                e16 = unsafe { simd16::tracked16_avx2(&x, &y, band, 1, 2, &mut d, &mut l, simd16::LIM, simd16::REBASE_AT) };
            }
            let t16 = t.elapsed().as_secs_f64() / reps as f64;
            let cells = (n as f64) * (2 * band + 1) as f64;
            eprintln!("band {}: 32-bit {:.3} ns/cell, 16-bit {:.3} ns/cell ({:?} {:?})", band, 1e9 * t32 / cells,
                      1e9 * t16 / cells, e32, e16);
        }
    }

    #[test]
    fn chain_members_follow_the_chain() {
        use super::super::{chain_members, chain_sorted, Params};
        let p = Params::default();
        // a colinear run with one off-diagonal anchor
        let a: Vec<(i64, i64)> = vec![(0, 10), (5, 15), (7, 40), (10, 20), (20, 30), (30, 40)];
        let m = chain_members(&a, &p);
        assert_eq!(m.first(), Some(&(0, 10)));
        assert_eq!(m.last(), Some(&(30, 40)));
        assert!(!m.contains(&(7, 40)));
        let (n, _, t0, t1, d0, d1) = chain_sorted(&a, &p).unwrap();
        assert_eq!((m.len(), m[0], *m.last().unwrap()), (n, (t0, d0), (t1, d1)));
        // unsorted input = chain_anchors
        let mut b = a.clone();
        b.reverse();
        assert_eq!(chain_members(&b, &p), m);
    }

}
