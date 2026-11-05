//! A run's bases-per-token scale from the raw signal alone: each read's signal-estimated length (samples x bases
//! per second / sample rate, in bp) over its token count, the median over reads. The assembly's limits given in bp
//! (output length, tips, superbubbles, ...) are converted to tokens with it.
//!
//! - `lengths --pod5 P... --output OUT`: POD5 file(s) or directories -> read_id, samples, sample_rate (sorted by read
//!   id), from the reads table alone (no signal is decompressed).
//! - `scale --lengths L --reads R [--bases-per-second B] --output OUT`: bp per token, the median over reads.
//! - `counts --tokens TOKENS.tsv.gz --output OUT.tsv.gz [--pigz T]`: tokens per read of a token table (first-seen
//!   order); `--pigz T` decompresses the input with `pigz -dc -p T`.
//! - `tokens --scale S --bp N`: a length in bp as tokens (printed).

use crate::io::{open_text, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::collections::{BTreeMap, HashMap};
use std::convert::{TryFrom, TryInto};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};


pub const BASES_PER_SECOND: f64 = 400.0;

/// Collects the `*.pod5` files under `dir`, recursively (symlinked folders are not entered; the caller sorts them).
fn pod5_files(dir: &Path, out: &mut Vec<PathBuf>) -> Result<()> {
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let p = entry.path();
        if entry.file_name().to_string_lossy().ends_with(".pod5") {
            out.push(p.clone());
        }
        if entry.file_type()?.is_dir() {
            pod5_files(&p, out)?;
        }
    }
    Ok(())
}

fn sort_paths(paths: &mut [PathBuf]) {
    paths.sort_by_key(|p| p.components().map(|c| c.as_os_str().to_string_lossy().to_string()).collect::<Vec<_>>());
}

const POD5_SIGNATURE: [u8; 8] = [0x8B, b'P', b'O', b'D', 0x0D, 0x0A, 0x1A, 0x0A];

fn read_at(f: &mut File, pos: u64, len: u64) -> Result<Vec<u8>> {
    let mut buf = vec![0u8; len as usize];
    f.seek(SeekFrom::Start(pos))?;
    f.read_exact(&mut buf)?;
    Ok(buf)
}

/// The POD5 footer's embedded tables (FlatBuffer: Footer.contents = [EmbeddedFile {offset, length, format,
/// content_type}]): content type -> (offset, length).
fn pod5_sections(fb: &[u8]) -> Result<HashMap<i16, (u64, u64)>> {
    let bad = || -> Box<dyn Error + Send + Sync> { "malformed POD5 footer".into() };
    let u32_at = |p: usize| fb.get(p..p + 4).map(|b| u32::from_le_bytes(b.try_into().unwrap())).ok_or_else(bad);
    let field = |table: usize, id: usize| -> Result<Option<usize>> {
        let vt = (table as i64 - i32::from_le_bytes(fb.get(table..table + 4).ok_or_else(bad)?.try_into().unwrap()) as i64) as usize;
        let size = u16::from_le_bytes(fb.get(vt..vt + 2).ok_or_else(bad)?.try_into().unwrap()) as usize;
        let entry = vt + 4 + id * 2;
        if entry + 2 > vt + size {
            return Ok(None);
        }
        let off = u16::from_le_bytes(fb.get(entry..entry + 2).ok_or_else(bad)?.try_into().unwrap()) as usize;
        Ok(if off == 0 { None } else { Some(table + off) })
    };
    let root = u32_at(0)? as usize;
    let mut out = HashMap::new();
    if let Some(pos) = field(root, 3)? {
        let v = pos + u32_at(pos)? as usize;
        for i in 0..u32_at(v)? as usize {
            let at = v + 4 + i * 4;
            let t = at + u32_at(at)? as usize;
            let i64_of = |id| -> Result<i64> {
                Ok(match field(t, id)? { Some(p) => i64::from_le_bytes(fb.get(p..p + 8).ok_or_else(bad)?.try_into().unwrap()), None => 0 })
            };
            let ct = match field(t, 3)? { Some(p) => i16::from_le_bytes(fb.get(p..p + 2).ok_or_else(bad)?.try_into().unwrap()), None => 0 };
            out.insert(ct, (i64_of(0)? as u64, i64_of(1)? as u64));
        }
    }
    Ok(out)
}

fn ipc_batches(f: &mut File, (off, len): (u64, u64)) -> Result<Vec<arrow::record_batch::RecordBatch>> {
    let reader = arrow::ipc::reader::FileReader::try_new(std::io::Cursor::new(read_at(f, off, len)?), None)?;
    Ok(reader.collect::<std::result::Result<Vec<_>, _>>()?)
}

/// (read id, samples, sample rate) of one POD5 file from its Reads and RunInfo tables alone (pod5lib's Reader also
/// indexes every signal batch, which means reads scattered through the whole file); None when the tables are not laid
/// out as expected (the caller then uses pod5lib).
fn pod5_reads_table(path: &Path) -> Result<Option<Vec<(String, u64, u16)>>> {
    use arrow::array::{Array, DictionaryArray, FixedSizeBinaryArray, StringArray, UInt16Array, UInt64Array};
    use arrow::datatypes::Int16Type;
    let mut f = File::open(path).map_err(|e| format!("{}: {}", path.display(), e))?;
    let size = f.seek(SeekFrom::End(0))?;
    if size < 64 {
        return Ok(None);
    }
    // ... [FOOTER\0\0][footer FlatBuffer, padded to 8][footer length i64][section marker 16][signature 8]
    let tail = read_at(&mut f, size - 32, 32)?;
    if tail[24..32] != POD5_SIGNATURE {
        return Err(format!("{}: not a POD5 file", path.display()).into());
    }
    let fb_len = u64::try_from(i64::from_le_bytes(tail[0..8].try_into().unwrap())).map_err(|_| "bad POD5 footer length")?;
    let fb_start = (size - 32).checked_sub((fb_len + 7) & !7).filter(|&s| s >= 8).ok_or("bad POD5 footer length")?;
    if read_at(&mut f, fb_start - 8, 8)? != b"FOOTER\0\0" {
        return Err(format!("{}: POD5 footer magic not found", path.display()).into());
    }
    let sections = pod5_sections(&read_at(&mut f, fb_start, fb_len)?)?;
    let (Some(&reads), Some(&runinfo)) = (sections.get(&0), sections.get(&4)) else { return Ok(None) };
    let mut rate_of: HashMap<String, u16> = HashMap::new();
    for b in ipc_batches(&mut f, runinfo)? {
        let ids = b.column_by_name("acquisition_id").and_then(|c| c.as_any().downcast_ref::<StringArray>());
        let rates = b.column_by_name("sample_rate").and_then(|c| c.as_any().downcast_ref::<UInt16Array>());
        let (Some(ids), Some(rates)) = (ids, rates) else { return Ok(None) };
        for i in 0..b.num_rows() {
            rate_of.insert(ids.value(i).to_string(), rates.value(i));
        }
    }
    let mut out = Vec::new();
    for b in ipc_batches(&mut f, reads)? {
        let rid = b.column_by_name("read_id").and_then(|c| c.as_any().downcast_ref::<FixedSizeBinaryArray>());
        let ns = b.column_by_name("num_samples").and_then(|c| c.as_any().downcast_ref::<UInt64Array>());
        let ri = b.column_by_name("run_info").and_then(|c| c.as_any().downcast_ref::<DictionaryArray<Int16Type>>());
        let (Some(rid), Some(ns), Some(ri)) = (rid, ns, ri) else { return Ok(None) };
        let Some(names) = ri.values().as_any().downcast_ref::<StringArray>() else { return Ok(None) };
        if rid.value_length() != 16 {
            return Ok(None);
        }
        for i in 0..b.num_rows() {
            let key = usize::try_from(ri.keys().value(i)).ok().filter(|&k| k < names.len()).ok_or("run_info key out of range")?;
            let rate = *rate_of.get(names.value(key)).ok_or_else(|| format!("{}: unknown run_info {}", path.display(), names.value(key)))?;
            out.push((pod5lib::uuid_to_str(rid.value(i).try_into().unwrap()), ns.value(i), rate));
        }
    }
    Ok(Some(out))
}

/// read id -> (samples, sample rate) of every read of the POD5 inputs; also the number of files.
pub fn signal_lengths(inputs: &[PathBuf]) -> Result<(BTreeMap<String, (u64, u16)>, usize)> {
    let mut files = Vec::new();
    for p in inputs {
        if p.is_dir() {
            let mut found = Vec::new();
            pod5_files(p, &mut found)?;
            sort_paths(&mut found);
            files.extend(found);
        } else {
            files.push(p.clone());
        }
    }
    let mut rows = BTreeMap::new();
    for f in &files {
        let reads = match pod5_reads_table(f)? {
            Some(reads) => reads,
            None => {
                let reader = pod5lib::Reader::open(f).map_err(|e| format!("{}: {}", f.display(), e))?;
                reader.reads().iter().map(|r| (pod5lib::uuid_to_str(&r.read_id), r.num_samples, r.run_info.sample_rate)).collect()
            }
        };
        for (rid, samples, rate) in reads {
            if rows.contains_key(&rid) {
                return Err(format!("read {} appears twice", rid).into());
            }
            rows.insert(rid, (samples, rate));
        }
    }
    Ok((rows, files.len()))
}

fn tmp_of(path: &Path) -> PathBuf {
    let mut name = path.file_name().unwrap_or_default().to_os_string();
    name.push(".tmp");
    path.with_file_name(name)
}

/// A `lengths` table: read -> (samples, sample rate).
pub fn read_lengths(path: &Path) -> Result<HashMap<String, (i64, i64)>> {
    let mut lines = open_text(path)?.lines();
    let head = lines.next().ok_or_else(|| format!("{}: empty", path.display()))??;
    if head != "read_id\tsamples\tsample_rate" {
        return Err(format!("{}: header {:?}", path.display(), head).into());
    }
    let mut out = HashMap::new();
    for line in lines {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 3 {
            return Err(format!("{}: row {:?}", path.display(), line).into());
        }
        out.insert(f[0].to_string(), (f[1].parse()?, f[2].parse()?));
    }
    Ok(out)
}

/// A table with read_id and token_count columns, in first-seen order (a later row of a read replaces its count).
pub fn read_tokens(path: &Path) -> Result<Vec<(String, i64)>> {
    let mut lines = open_text(path)?.lines();
    let head = lines.next().ok_or_else(|| format!("{}: empty", path.display()))??;
    let head: Vec<&str> = head.split('\t').collect();
    let col = |c: &str| head.iter().position(|h| *h == c).ok_or_else(|| format!("{}: no column {}", path.display(), c));
    let (i, j) = (col("read_id")?, col("token_count")?);
    let (mut out, mut at): (Vec<(String, i64)>, HashMap<String, usize>) = (Vec::new(), HashMap::new());
    for line in lines {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        let (r, n) = (f.get(i).ok_or("short row")?.to_string(), f.get(j).ok_or("short row")?.parse::<i64>()?);
        match at.get(&r) {
            Some(&k) => out[k].1 = n,
            None => {
                at.insert(r.clone(), out.len());
                out.push((r, n));
            }
        }
    }
    Ok(out)
}

pub struct Scale {
    pub bp_per_token: f64,
    pub pooled_bp_per_token: f64, // all signal-estimated bases / all tokens
    pub reads: usize,
    pub reads_without_signal: usize, // reads with tokens but no signal length
    pub sample_rates: String, // the distinct sample rates (Hz), comma-separated
}

/// The median over reads (in both tables, tokens > 0) of signal-estimated bases / tokens.
pub fn scale_from(lengths: &HashMap<String, (i64, i64)>, tokens: &[(String, i64)], bases_per_second: f64) -> Result<Scale> {
    let (mut ratios, mut bases, mut count) = (Vec::new(), 0.0f64, 0i64);
    for (read, n) in tokens {
        if *n <= 0 {
            continue;
        }
        let Some(&(samples, rate)) = lengths.get(read) else { continue };
        if rate <= 0 {
            return Err(format!("{}: sample rate {}", read, rate).into());
        }
        let b = samples as f64 * bases_per_second / rate as f64;
        ratios.push(b / *n as f64);
        bases += b;
        count += n;
    }
    if ratios.is_empty() {
        return Err("no read has both a signal length and tokens".into());
    }
    let reads = ratios.len();
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = reads / 2;
    let median = if reads % 2 == 1 { ratios[m] } else { (ratios[m - 1] + ratios[m]) / 2.0 };
    let mut rates: Vec<i64> = tokens.iter().filter_map(|(r, _)| lengths.get(r).map(|x| x.1)).collect();
    rates.sort_unstable();
    rates.dedup();
    Ok(Scale {
        bp_per_token: median,
        pooled_bp_per_token: bases / count as f64,
        reads,
        reads_without_signal: tokens.iter().filter(|(r, n)| *n > 0 && !lengths.contains_key(r)).count(),
        sample_rates: rates.iter().map(|r| r.to_string()).collect::<Vec<_>>().join(","),
    })
}

/// Tokens per read of a token table (one row per token, with a read_id column), in first-seen order.
pub fn token_counts(path: &Path, pigz: usize) -> Result<Vec<(Vec<u8>, u64)>> {
    let mut child = None;
    let mut reader: Box<dyn BufRead> = if pigz > 0 {
        let mut c = std::process::Command::new("pigz").args(["-dc", "-p", &pigz.to_string()]).arg(path)
            .stdout(std::process::Stdio::piped()).spawn().map_err(|e| format!("pigz: {}", e))?;
        let out = c.stdout.take().ok_or("pigz: no stdout")?;
        child = Some(c);
        Box::new(BufReader::with_capacity(1 << 20, out))
    } else {
        let file = File::open(path).map_err(|e| format!("{}: {}", path.display(), e))?;
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    };
    let mut line = Vec::new();
    reader.read_until(b'\n', &mut line)?;
    let head = String::from_utf8_lossy(&line).trim_end_matches('\n').to_string();
    let i = head.split('\t').position(|h| h == "read_id").ok_or_else(|| format!("{}: no read_id column", path.display()))?;
    let (mut out, mut at): (Vec<(Vec<u8>, u64)>, HashMap<Vec<u8>, usize>) = (Vec::new(), HashMap::new());
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        // column i (when read_id is the last column, the newline stays on it)
        let rid = line.splitn(i + 2, |&b| b == b'\t').nth(i).ok_or_else(|| format!("{}: short row", path.display()))?;
        match at.get(rid) {
            Some(&k) => out[k].1 += 1,
            None => {
                at.insert(rid.to_vec(), out.len());
                out.push((rid.to_vec(), 1));
            }
        }
    }
    drop(reader);
    if let Some(mut c) = child {
        let status = c.wait()?;
        if !status.success() {
            return Err(format!("pigz -dc {}: {}", path.display(), status).into());
        }
    }
    Ok(out)
}

/// The scale of a `scale` table.
pub fn bp_per_token(path: &Path) -> Result<f64> {
    let text = fs::read_to_string(path).map_err(|e| format!("{}: {}", path.display(), e))?;
    let value: f64 = text.lines().skip(1).filter_map(|l| l.split_once('\t')).find(|(k, _)| *k == "bp_per_token")
        .ok_or_else(|| format!("{}: no bp_per_token", path.display()))?.1.trim().parse()?;
    if !(value > 0.0) {
        return Err(format!("{}: bp_per_token {}", path.display(), value).into());
    }
    Ok(value)
}

/// A length in bp as tokens at `scale` bp per token: max(1, round(bp / scale)), rounding half to even.
pub fn bp_to_tokens(bp: i64, scale: f64) -> Result<i64> {
    if bp <= 0 || !(scale > 0.0) {
        return Err(format!("bp {} and scale {} must be positive", bp, scale).into());
    }
    Ok(((bp as f64 / scale).round_ties_even() as i64).max(1))
}

fn resolve(p: &Path) -> Result<String> {
    Ok(fs::canonicalize(p).map_err(|e| format!("{}: {}", p.display(), e))?.display().to_string())
}

const USAGE: &str = "usage: squash __bp-scale lengths --pod5 POD5... --output LENGTHS.tsv\n\
       squash __bp-scale scale --lengths LENGTHS.tsv --reads TABLE [--bases-per-second 400] --output SCALE.tsv\n\
       squash __bp-scale counts --tokens TOKENS.tsv.gz --output COUNTS.tsv.gz [--pigz THREADS]";

pub fn run() -> Result<()> {
    let mut argv = crate::argv::args().skip(1);
    let command = argv.next().ok_or(USAGE)?;
    let mut map: HashMap<String, Vec<String>> = HashMap::new();
    let mut key: Option<String> = None;
    for arg in argv {
        if let Some(k) = arg.strip_prefix("--") {
            map.entry(k.to_string()).or_default();
            key = Some(k.to_string());
        } else {
            let k = key.as_ref().ok_or_else(|| format!("unexpected argument {:?}\n{}", arg, USAGE))?;
            map.get_mut(k).unwrap().push(arg);
        }
    }
    let allowed: &[&str] = match command.as_str() {
        "lengths" => &["pod5", "output"],
        "scale" => &["lengths", "reads", "bases-per-second", "output"],
        "counts" => &["tokens", "output", "pigz"],
        _ => return Err(format!("unknown subcommand {}\n{}", command, USAGE).into()),
    };
    if let Some(k) = map.keys().find(|k| !allowed.contains(&k.as_str())) {
        return Err(format!("unknown option --{}\n{}", k, USAGE).into());
    }
    fn one(map: &mut HashMap<String, Vec<String>>, k: &str) -> Result<String> {
        let v = map.remove(k).ok_or_else(|| format!("missing --{}\n{}", k, USAGE))?;
        if v.len() != 1 {
            return Err(format!("--{} takes one value\n{}", k, USAGE).into());
        }
        Ok(v.into_iter().next().unwrap())
    }
    let output = PathBuf::from(one(&mut map, "output")?);
    if output.exists() {
        return Err(format!("refusing to overwrite {}", output.display()).into());
    }
    let tmp = tmp_of(&output);
    match command.as_str() {
        "lengths" => {
            let inputs: Vec<PathBuf> = map.remove("pod5").ok_or_else(|| format!("missing --pod5\n{}", USAGE))?
                .into_iter().map(PathBuf::from).collect();
            if inputs.is_empty() {
                return Err(format!("--pod5 needs a value\n{}", USAGE).into());
            }
            let (rows, files) = signal_lengths(&inputs)?;
            let mut out = BufWriter::new(File::create(&tmp)?);
            out.write_all(b"read_id\tsamples\tsample_rate\n")?;
            for (rid, (samples, rate)) in &rows {
                writeln!(out, "{}\t{}\t{}", rid, samples, rate)?;
            }
            out.flush()?;
            drop(out);
            fs::rename(&tmp, &output)?;
            println!("{} reads from {} POD5 file(s)", rows.len(), files);
        }
        "scale" => {
            let lengths = PathBuf::from(one(&mut map, "lengths")?);
            let reads = PathBuf::from(one(&mut map, "reads")?);
            let bps: f64 = if map.contains_key("bases-per-second") { one(&mut map, "bases-per-second")?.parse()? } else { BASES_PER_SECOND };
            let s = scale_from(&read_lengths(&lengths)?, &read_tokens(&reads)?, bps)?;
            let text = format!("key\tvalue\nbp_per_token\t{:.6}\npooled_bp_per_token\t{:.6}\nreads\t{}\nreads_without_signal\t{}\n\
                                bases_per_second\t{:.6}\nsample_rates\t{}\nlengths\t{}\nreads_table\t{}\n",
                               s.bp_per_token, s.pooled_bp_per_token, s.reads, s.reads_without_signal, bps, s.sample_rates,
                               resolve(&lengths)?, resolve(&reads)?);
            fs::write(&tmp, text)?;
            fs::rename(&tmp, &output)?;
            println!("bp per token {:.4} (pooled {:.4}) over {} reads", s.bp_per_token, s.pooled_bp_per_token, s.reads);
        }
        "counts" => {
            let tokens = PathBuf::from(one(&mut map, "tokens")?);
            let pigz: usize = if map.contains_key("pigz") { one(&mut map, "pigz")?.parse()? } else { 0 };
            let counts = token_counts(&tokens, pigz)?;
            let mut gz = GzEncoder::new(BufWriter::new(File::create(&tmp)?), Compression::default());
            gz.write_all(b"read_id\ttoken_count\n")?;
            for (rid, n) in &counts {
                gz.write_all(rid)?;
                writeln!(gz, "\t{}", n)?;
            }
            gz.finish()?.flush()?;
            fs::rename(&tmp, &output)?;
        }
        _ => unreachable!(),
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scale_is_the_median_of_signal_bases_per_token() {
        // 5 kHz at 400 bases/s: 12.5 samples per base.
        let lengths: HashMap<String, (i64, i64)> = [("a", (12500, 5000)), ("b", (25000, 5000)), ("c", (5000, 5000))]
            .iter().map(|(r, x)| (r.to_string(), *x)).collect();
        let tokens: Vec<(String, i64)> = [("a", 100), ("b", 100), ("c", 100), ("d", 7), ("e", 0)]
            .iter().map(|(r, n)| (r.to_string(), *n)).collect();
        let s = scale_from(&lengths, &tokens, 400.0).unwrap();
        assert_eq!((s.bp_per_token, s.reads, s.reads_without_signal), (10.0, 3, 1));
        assert!((s.pooled_bp_per_token - 3400.0 / 300.0).abs() < 1e-12);
        assert_eq!(s.sample_rates, "5000");
    }

    #[test]
    fn even_count_median_is_the_mean_of_the_middle_two() {
        let lengths: HashMap<String, (i64, i64)> = [("a", (1000, 4000)), ("b", (3000, 4000))]
            .iter().map(|(r, x)| (r.to_string(), *x)).collect();
        let tokens: Vec<(String, i64)> = vec![("a".into(), 10), ("b".into(), 10)];
        assert_eq!(scale_from(&lengths, &tokens, 400.0).unwrap().bp_per_token, 20.0);
    }

    #[test]
    fn bp_to_tokens_rounds_half_to_even() {
        assert_eq!(bp_to_tokens(25, 10.0).unwrap(), 2);
        assert_eq!(bp_to_tokens(35, 10.0).unwrap(), 4);
        assert_eq!(bp_to_tokens(1, 10.0).unwrap(), 1);
        assert_eq!(bp_to_tokens(50000, 10.81).unwrap(), 4625);
        assert!(bp_to_tokens(0, 10.0).is_err());
    }
}
