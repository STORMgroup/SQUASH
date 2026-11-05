//! The `partners` step: its usage, entry point and the mmap threshold of the correction steps.

use super::*;

const USAGE: &str = "usage: squash __partners --reads-dir TOKENS --output TABLE.tsv [--threads N] [--max-reads N] \
[--min-anchors N] [--max-gap N] [--max-drift N] [--drift-fraction F] [--lookback N] [--min-overlap N] [--min-cover F] \
[--coverage N] [--repeat-fraction F] [--exclude-reads TABLE] \
[--signal-symbols SYMBOLS.tsv --signal-thr F --signal-band N --signal-path N --signal-edit-gap N [--signal-window N] \
[--signal-path-shift N] [--signal-min-cover F] [--signal-verdicts TABLE] [--signal-verdicts-out TABLE] [--signal-store DIR]]";


/// Fixes glibc's mmap threshold at 4 MB, so that every allocation of 4 MB or more is its own mapping, returned to the
/// system when freed, instead of the dynamic threshold rising with the largest block freed (up to 32 MB) and such
/// blocks staying in the arenas.
pub fn set_mmap_threshold() {
    const MMAP_THRESHOLD: i32 = 4 << 20;
    unsafe {
        libc::mallopt(libc::M_MMAP_THRESHOLD, MMAP_THRESHOLD);
    }
}

/// The `partners` step: loads the tokens (`--reads-dir`), the excluded reads and the q-symbol alignment check
/// (`--signal-*`), searches every target and writes the selected-donor table to `--output` (uncompressed, through a
/// temporary file; an existing output is never overwritten), and with `--signal-verdicts-out` the check's verdicts.
pub fn run() -> Result<()> {
    set_mmap_threshold();
    let _memtrack = memtrack::Sampler::start();   // SQUASH_PARTNERS_MEMTRACK=1 (a diagnostic)
    let mut map: HashMap<String, String> = HashMap::new();
    let mut argv = crate::argv::args().skip(1);
    while let Some(key) = argv.next() {
        if !key.starts_with("--") {
            return Err(format!("unexpected argument {:?}\n{}", key, USAGE).into());
        }
        let value = argv.next().ok_or_else(|| format!("{} needs a value\n{}", key, USAGE))?;
        map.insert(key[2..].to_string(), value);
    }
    let p = take_params(&mut map)?;
    let sig = signal::take(&mut map)?;
    let threads: usize = get(&mut map, "threads", 1usize)?;
    let output = PathBuf::from(map.remove("output").ok_or_else(|| format!("missing --output\n{}", USAGE))?);
    let reads_dir = PathBuf::from(map.remove("reads-dir").ok_or_else(|| format!("missing --reads-dir\n{}", USAGE))?);
    let exclude_path = map.remove("exclude-reads").map(PathBuf::from);
    if let Some(key) = map.keys().next() {
        return Err(format!("unknown option --{}\n{}", key, USAGE).into());
    }
    if output.exists() {
        return Err(format!("refusing to overwrite {}", output.display()).into());
    }
    let (ids, lengths, tokens) = load_prepared(&reads_dir)?;
    let unique: HashSet<&String> = ids.iter().collect();
    if unique.len() != ids.len() {
        return Err("duplicate read ids".into());
    }
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads.max(1)).build()?;
    if let Some(o) = &sig {
        signal::install(o, &ids, &pool)?;
        if let Some(dir) = &o.store {
            signal::store_on_file(dir)?;
        }
    }
    let excluded = match &exclude_path {
        None => Vec::new(),
        Some(path) => load_excluded(path, &ids)?,
    };
    // each target's rows packed as found, the table written from them
    let rows = partner_rows_packed(&ids, &lengths, &tokens, None, &p, &excluded, &pool, None, None);
    if let Some(o) = &sig {
        signal::report("first search");
        if let Some(path) = &o.verdicts_out {
            signal::write_verdicts(path, &ids)?;
        }
    }
    let tmp = output.with_extension("tmp");
    write_packed_rows(&tmp, &rows)?;
    std::fs::rename(&tmp, &output)?;
    Ok(())
}
