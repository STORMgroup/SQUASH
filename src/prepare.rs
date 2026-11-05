//! The rest of the signal stage: from the signal step's q-symbols (`signal/qsymbols.tsv`) to the work folder's
//! uncorrected tokens and its list of low-complexity reads. Two steps, each a child process:
//!
//! - `__graph prepare` (`graph/`) sketches the q-symbols into the uncorrected tokens, `tokens/` (`reads.tsv.gz`,
//!   `tokens.u32le`, `token_spans.u64le`; read by the partner search, the rounds and the layout).
//! - `__low-complexity` lists the reads whose most repeated token value is at least `--min-frac` of their tokens (0.2
//!   in both presets); the partner search and the correction rounds take the list as `--exclude-reads`.

use crate::io::{READS_TABLE, TOKEN_VALUES};
use crate::pipeline::work;
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use structopt::StructOpt;

/// A `low-complexity` value: `off`, empty or zero = off (None), else the fraction as given (the step's `--min-frac`),
/// which must be in (0, 1].
pub fn low_complexity_value(v: &str) -> Result<Option<&str>, String> {
    let v = v.trim();
    if v.is_empty() || v == "off" {
        return Ok(None);
    }
    match v.parse::<f64>() {
        Ok(f) if f == 0.0 => Ok(None),
        Ok(f) if f > 0.0 && f <= 1.0 => Ok(Some(v)),
        _ => Err(format!("low-complexity must be a fraction in (0, 1] or off, not {:?}", v)),
    }
}

/// The rest of the signal stage (`pipeline.rs` runs it after the signal step), with the same arguments in both
/// presets:
///
/// 1. `prepare-tokens` (`squash __graph prepare --q-symbols ROOT/signal/qsymbols.tsv --signal-events E
///    --signal-window W --output-dir ROOT/tokens --threads T`);
/// 2. `low-complexity` (`squash __low-complexity --reads-dir ROOT/tokens --min-frac F
///    --output ROOT/correction/low_complexity.tsv --threads T`), skipped when the rule is off.
pub fn run(r: &crate::pipeline::Run) -> Result<(), String> {
    // a rerun of the stage replaces these steps' outputs (and correction/, which is made from them)
    clear_outputs(&r.root)?;
    let path = |x: &Path| x.display().to_string();
    let tokens = r.root.join(work::TOKENS);
    let prepare: Vec<String> = vec![
        "prepare".into(),
        "--q-symbols".into(),
        path(&r.root.join(work::QSYMBOLS)),
        "--signal-events".into(),
        r.v("events").into(),
        "--signal-window".into(),
        r.v("window").into(),
        "--output-dir".into(),
        path(&tokens),
        "--threads".into(),
        r.threads(),
    ];
    r.runner.step("prepare-tokens", "graph", &prepare)?;
    if let Some(frac) = low_complexity_value(r.v("low-complexity"))? {
        let list = r.root.join(work::LOW_COMPLEXITY);
        fs::create_dir_all(list.parent().unwrap()).map_err(|e| e.to_string())?;
        let low: Vec<String> = vec![
            "--reads-dir".into(),
            path(&tokens),
            "--min-frac".into(),
            frac.to_string(),
            "--output".into(),
            path(&list),
            "--threads".into(),
            r.threads(),
        ];
        r.runner.step("low-complexity", "low-complexity", &low)?;
    }
    Ok(())
}

/// The parts of the work folder made from the q-symbols (top-level names): removed before the stage runs again.
const PARTS: [&str; 2] = ["tokens", "correction"];

/// Removes these steps' outputs from a work folder. Everything in `correction/` goes too: it is made from them.
pub fn clear_outputs(root: &Path) -> Result<(), String> {
    for part in PARTS.iter() {
        let p = root.join(part);
        let gone = if p.is_dir() { fs::remove_dir_all(&p) } else if p.exists() { fs::remove_file(&p) } else { Ok(()) };
        gone.map_err(|e| format!("cannot remove {}: {}", p.display(), e))?;
    }
    Ok(())
}

fn write_text(path: &Path, text: &str) -> io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    w.write_all(text.as_bytes())?;
    w.flush()
}

// ------------------------------------------------------------------------------------------------------------------
// Step `__low-complexity`

#[derive(Debug, StructOpt)]
#[structopt(name = "low-complexity", about = "the reads whose most repeated token value is >= F of their tokens")]
struct LowComplexityOpt {
    /// The uncorrected tokens' folder (`graph prepare`'s reads.tsv.gz and tokens.u32le).
    #[structopt(long, parse(from_os_str))]
    reads_dir: PathBuf,
    #[structopt(long)]
    min_frac: f64,
    #[structopt(long, parse(from_os_str))]
    output: PathBuf,
    /// Threads (the table does not depend on them).
    #[structopt(long, default_value = "1")]
    threads: usize,
}

/// One listed read: id, tokens, the most repeated value's count.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LowRow {
    pub read_id: String,
    pub tokens: u64,
    pub top: u64,
}

fn check_min_frac(f: f64) -> Result<(), String> {
    if f > 0.0 && f <= 1.0 {
        Ok(())
    } else {
        Err("--min-frac must be in (0, 1]".into())
    }
}

/// The count of the most repeated value (0 for no values).
fn top_count(values: &[u32]) -> u64 {
    let mut v = values.to_vec();
    v.sort_unstable();
    let (mut best, mut run) = (0u64, 0u64);
    for i in 0..v.len() {
        run = if i > 0 && v[i] == v[i - 1] { run + 1 } else { 1 };
        best = best.max(run);
    }
    best
}

/// Tokens read per batch (the batch's top counts are computed in parallel).
const LOW_BATCH_TOKENS: usize = 1 << 24;

/// The low-complexity reads: every read (in the reads table's order) with at least one token whose most repeated token
/// value holds `top >= min_frac * tokens` of them. Reads the table by its column names and each read's `token_count`
/// values at `token_offset` of the tokens array (seeking when a read does not follow the last one); a range beyond the
/// array is an error.
pub fn low_complexity(dir: &Path, min_frac: f64, pool: &rayon::ThreadPool, batch_tokens: usize) -> Result<Vec<LowRow>, String> {
    let reads_path = dir.join(READS_TABLE);
    let tokens_path = dir.join(TOKEN_VALUES);
    let fail = |p: &Path, m: String| format!("{}: {}", p.display(), m);
    let file = File::open(&reads_path).map_err(|e| fail(&reads_path, e.to_string()))?;
    let mut lines = BufReader::with_capacity(1 << 20, flate2::read::MultiGzDecoder::new(file)).lines();
    let header = match lines.next() {
        Some(h) => h.map_err(|e| fail(&reads_path, e.to_string()))?,
        None => return Err(fail(&reads_path, "empty reads table".into())),
    };
    let h: Vec<&str> = header.split('\t').collect();
    let col = |name: &str| h.iter().position(|c| *c == name).ok_or_else(|| fail(&reads_path, format!("no column {}", name)));
    let (i_id, i_off, i_n) = (col("read_id")?, col("token_offset")?, col("token_count")?);
    let need = i_id.max(i_off).max(i_n);
    let mut reads: Vec<(String, u64, u64)> = Vec::new();
    for line in lines {
        let line = line.map_err(|e| fail(&reads_path, e.to_string()))?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() <= need {
            return Err(fail(&reads_path, format!("malformed row {:?}", line)));
        }
        let num = |i: usize| f[i].parse::<u64>().map_err(|_| fail(&reads_path, format!("bad number in row {:?}", line)));
        reads.push((f[i_id].to_string(), num(i_off)?, num(i_n)?));
    }
    let token_file = File::open(&tokens_path).map_err(|e| fail(&tokens_path, e.to_string()))?;
    let bytes = token_file.metadata().map_err(|e| fail(&tokens_path, e.to_string()))?.len();
    if bytes % 4 != 0 {
        return Err(fail(&tokens_path, "size is not a multiple of 4".into()));
    }
    let total = bytes / 4;
    let mut tokens = BufReader::with_capacity(1 << 22, token_file);
    let mut at = 0u64; // the reader's position, in tokens
    let mut out = Vec::new();
    let mut start = 0usize;
    while start < reads.len() {
        // A batch: consecutive reads with tokens, up to batch_tokens values (at least one read).
        let mut end = start;
        let mut sum = 0usize;
        while end < reads.len() && (end == start || sum + reads[end].2 as usize <= batch_tokens) {
            sum += reads[end].2 as usize;
            end += 1;
        }
        let mut block: Vec<u32> = Vec::with_capacity(sum);
        let mut spans = Vec::with_capacity(end - start);
        let mut raw = Vec::new();
        for (_, off, n) in &reads[start..end] {
            let (off, n) = (*off, *n);
            if n == 0 {
                spans.push((block.len(), 0));
                continue;
            }
            if off.checked_add(n).map_or(true, |e| e > total) {
                return Err(fail(&tokens_path, format!("token range {}+{} beyond its {} tokens", off, n, total)));
            }
            if off != at {
                tokens.seek(SeekFrom::Start(off * 4)).map_err(|e| fail(&tokens_path, e.to_string()))?;
            }
            raw.resize(n as usize * 4, 0);
            tokens.read_exact(&mut raw).map_err(|e| fail(&tokens_path, e.to_string()))?;
            at = off + n;
            spans.push((block.len(), n as usize));
            block.extend(raw.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])));
        }
        let tops: Vec<u64> =
            pool.install(|| spans.par_iter().map(|&(s, n)| if n == 0 { 0 } else { top_count(&block[s..s + n]) }).collect());
        for ((id, _, n), top) in reads[start..end].iter().zip(tops) {
            if *n == 0 {
                continue;
            }
            if top as f64 >= min_frac * *n as f64 {
                out.push(LowRow { read_id: id.clone(), tokens: *n, top });
            }
        }
        start = end;
    }
    Ok(out)
}

/// The low-complexity table: `read_id, tokens, top_count, top_frac` (`top / tokens` with five decimals: `{:.5}` rounds
/// the exact value of the double, ties to even).
pub fn low_complexity_text(rows: &[LowRow]) -> String {
    let mut out = String::from("read_id\ttokens\ttop_count\ttop_frac\n");
    for r in rows {
        out.push_str(&format!("{}\t{}\t{}\t{:.5}\n", r.read_id, r.tokens, r.top, r.top as f64 / r.tokens as f64));
    }
    out
}

pub fn run_low_complexity() -> Result<(), String> {
    let o = LowComplexityOpt::from_iter_safe(crate::argv::args()).map_err(|e| e.to_string())?;
    check_min_frac(o.min_frac)?;
    if o.threads == 0 {
        return Err("--threads must be at least one".into());
    }
    let pool = rayon::ThreadPoolBuilder::new().num_threads(o.threads).build().map_err(|e| e.to_string())?;
    let rows = low_complexity(&o.reads_dir, o.min_frac, &pool, LOW_BATCH_TOKENS)?;
    // written to .NAME.tmp beside the output, then renamed
    let name = o.output.file_name().ok_or_else(|| format!("{}: not a file name", o.output.display()))?;
    let tmp = o.output.with_file_name(format!(".{}.tmp", name.to_string_lossy()));
    write_text(&tmp, &low_complexity_text(&rows)).map_err(|e| format!("{}: {}", tmp.display(), e))?;
    fs::rename(&tmp, &o.output).map_err(|e| format!("{}: {}", o.output.display(), e))?;
    println!("{} low-complexity reads (most repeated token >= {} of the read) -> {}", rows.len(), o.min_frac,
             o.output.display());
    Ok(())
}

// ------------------------------------------------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static COUNTER: AtomicUsize = AtomicUsize::new(0);

    fn test_dir(label: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!(
            "squash_prepare_{}_{}_{}",
            label,
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        if d.exists() {
            fs::remove_dir_all(&d).unwrap();
        }
        fs::create_dir_all(&d).unwrap();
        d
    }

    /// DIR/reads.tsv.gz and DIR/tokens.u32le as `graph prepare` writes them.
    fn write_prepared(dir: &Path, reads: &[(&str, Vec<u32>)]) {
        fs::create_dir_all(dir).unwrap();
        let gz = File::create(dir.join(READS_TABLE)).unwrap();
        let mut w = flate2::write::GzEncoder::new(gz, flate2::Compression::default());
        writeln!(w, "read_id\tread_ordinal\tsource_symbols\ttoken_offset\ttoken_count").unwrap();
        let (mut off, mut flat) = (0usize, Vec::new());
        for (i, (id, vals)) in reads.iter().enumerate() {
            writeln!(w, "{}\t{}\t{}\t{}\t{}", id, i, 10 * vals.len(), off, vals.len()).unwrap();
            off += vals.len();
            flat.extend(vals.iter().flat_map(|v| v.to_le_bytes()));
        }
        w.finish().unwrap();
        fs::write(dir.join(TOKEN_VALUES), flat).unwrap();
    }

    fn reads_fixture() -> Vec<(&'static str, Vec<u32>)> {
        vec![
            ("normal", (0..100).collect()),                                       // top 1 of 100
            ("stall", std::iter::repeat(7).take(60).chain(100..140).collect()), // top 60 of 100
            ("edge", std::iter::repeat(9).take(20).chain(200..280).collect()),  // exactly 0.2
            ("empty", vec![]),
            ("short", vec![5, 5, 6, 6, 6]), // top 3 of 5
        ]
    }

    #[test]
    fn low_complexity_fraction_and_order() {
        let d = test_dir("low");
        let prefix = d.join("reads");
        write_prepared(&prefix, &reads_fixture());
        let pool = rayon::ThreadPoolBuilder::new().num_threads(2).build().unwrap();
        let got = |f: f64, batch: usize| -> Vec<(String, u64)> {
            low_complexity(&prefix, f, &pool, batch).unwrap().into_iter().map(|r| (r.read_id, r.top)).collect()
        };
        let want: Vec<(String, u64)> = vec![("stall".into(), 60), ("edge".into(), 20), ("short".into(), 3)];
        for batch in [1usize, 7, 100, 1 << 20].iter() {
            assert_eq!(got(0.2, *batch), want);
        }
        assert_eq!(got(0.5, 1).into_iter().map(|x| x.0).collect::<Vec<_>>(), vec!["stall", "short"]);
        assert!(got(0.7, 3).is_empty());
        let rows = low_complexity(&prefix, 0.5, &pool, 64).unwrap();
        let text = low_complexity_text(&rows);
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "read_id\ttokens\ttop_count\ttop_frac");
        assert_eq!(lines[1], "stall\t100\t60\t0.60000");
        assert_eq!(lines[2], "short\t5\t3\t0.60000");
        fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn low_complexity_reads_by_offset() {
        // a reads table whose rows are not in token order: each read's own range is read
        let d = test_dir("offsets");
        let prefix = d.join("reads");
        fs::create_dir_all(&prefix).unwrap();
        let gz = File::create(prefix.join(READS_TABLE)).unwrap();
        let mut w = flate2::write::GzEncoder::new(gz, flate2::Compression::fast());
        write!(w, "token_count\tread_id\ttoken_offset\n4\tb\t5\n3\ta\t0\n").unwrap();
        w.finish().unwrap();
        let tokens: Vec<u8> = [1u32, 1, 2, 8, 8, 3, 3, 3, 4, 9].iter().flat_map(|v| v.to_le_bytes()).collect();
        fs::write(prefix.join(TOKEN_VALUES), tokens).unwrap();
        let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
        let rows = low_complexity(&prefix, 0.6, &pool, 2).unwrap();
        assert_eq!(rows, vec![LowRow { read_id: "b".into(), tokens: 4, top: 3 }, LowRow { read_id: "a".into(), tokens: 3, top: 2 }]);
        // a range beyond the array is refused
        let gz = File::create(prefix.join(READS_TABLE)).unwrap();
        let mut w = flate2::write::GzEncoder::new(gz, flate2::Compression::fast());
        write!(w, "read_id\ttoken_offset\ttoken_count\na\t8\t3\n").unwrap();
        w.finish().unwrap();
        assert!(low_complexity(&prefix, 0.2, &pool, 2).unwrap_err().contains("beyond"));
        fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn five_decimals_as_python() {
        // five decimals of the exact double, ties to even (13/64 = 0.203125 exactly)
        let row = |top: u64, n: u64| low_complexity_text(&[LowRow { read_id: "r".into(), tokens: n, top }]);
        assert!(row(13, 64).ends_with("\t0.20312\n"));
        assert!(row(15, 64).ends_with("\t0.23438\n"));
        assert!(row(1, 3).ends_with("\t0.33333\n"));
        assert!(row(2, 3).ends_with("\t0.66667\n"));
        assert!(row(64, 64).ends_with("\t1.00000\n"));
        assert!(row(1, 5).ends_with("\t0.20000\n"));
    }

    #[test]
    fn top_counts() {
        assert_eq!(top_count(&[]), 0);
        assert_eq!(top_count(&[4]), 1);
        assert_eq!(top_count(&[3, 1, 3, 2, 3, 1]), 3);
        assert_eq!(top_count(&[u32::MAX, u32::MAX, 0]), 2);
    }

    #[test]
    fn low_complexity_values() {
        // 0 or empty = off; else the value as given
        assert_eq!(low_complexity_value("0.2").unwrap(), Some("0.2"));
        assert_eq!(low_complexity_value("1").unwrap(), Some("1"));
        assert_eq!(low_complexity_value("off").unwrap(), None);
        assert_eq!(low_complexity_value("0").unwrap(), None);
        assert_eq!(low_complexity_value("").unwrap(), None);
        assert!(low_complexity_value("1.5").is_err());
        assert!(low_complexity_value("-0.1").is_err());
        assert!(low_complexity_value("x").is_err());
    }

    #[test]
    fn clear_outputs_removes_this_part_only() {
        let d = test_dir("clear");
        let root = d.join("root");
        for p in ["signal", "tokens", "correction/round3", "logs", "bp_scale"] {
            fs::create_dir_all(root.join(p)).unwrap();
        }
        fs::write(root.join(work::QSYMBOLS), "r1\t1\n").unwrap();
        clear_outputs(&root).unwrap();
        let mut left: Vec<String> =
            fs::read_dir(&root).unwrap().map(|e| e.unwrap().file_name().to_string_lossy().into_owned()).collect();
        left.sort();
        assert_eq!(left, vec!["bp_scale", "logs", "signal"]);
        assert!(fs::read(root.join(work::QSYMBOLS)).is_ok(), "the signal step's file stays");
        fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn min_frac_range() {
        assert!(check_min_frac(0.2).is_ok());
        assert!(check_min_frac(1.0).is_ok());
        assert!(check_min_frac(0.0).is_err());
        assert!(check_min_frac(1.0000001).is_err());
        assert!(check_min_frac(f64::NAN).is_err());
    }
}
