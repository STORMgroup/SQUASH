//! Sparse token selection for the assembly (`sparse-repeat`): the 1-in-N subset of a token table's values (by a value
//! hash), then metaMDBG's cut of repetitive values, in one pass over the input (the sparse reads kept in memory, 1 in
//! N tokens).
//!
//! - The sparse subset: keep a token when `splitmix64(value) < 2^64 / N` (all tokens when N = 1).
//! - The repeat cut: drop the `max(int(fraction x distinct values), 1)` most frequent token values of the subset (ties:
//!   smaller value first), metaMDBG's `determineRepetitiveMinimizers`.
//!
//! Writes the uncompressed token table (header + tab-separated rows, no quoting); the pipeline gzips it. Reads left
//! without tokens are dropped and the kept reads renumbered (after each filter). Rows must be grouped by read in
//! ordinal order; a read that reappears later is refused. `--pigz T` (T > 0) decompresses the input with
//! `pigz -dc -p T`.

use crate::io::{Result, open_text};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};


pub const HEADER: &str = "read_id\tread_ordinal\tsource_length\ttoken_ordinal\ttoken_id\tstart\tend\tcoordinate_space";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Token {
    pub value: u32,
    pub start: i64,
    pub end: i64,
}

/// One read of a token table.
pub struct Read {
    pub id: String,
    pub length: i64,
    pub tokens: Vec<Token>,
}

/// SplitMix64's output function: the value hash of the sparse selection.
pub fn splitmix64(value: u64) -> u64 {
    let mut x = value.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

/// Whether a token value is in the value-hash subset of 1 in `inverse_density` (every value when it is 1).
pub fn keep_token(value: u32, inverse_density: u64) -> bool {
    inverse_density == 1 || (splitmix64(value as u64) as u128) < (1u128 << 64) / inverse_density as u128
}

/// The `max(int(fraction x distinct), 1)` most frequent values as (value, count), most frequent first, ties by value.
pub fn repetitive_tokens(counts: &HashMap<u32, u64>, fraction: f64) -> Vec<(u32, u64)> {
    let mut ranked: Vec<(u32, u64)> = counts.iter().map(|(&v, &c)| (v, c)).collect();
    ranked.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let n = ((fraction * ranked.len() as f64) as usize).max(1).min(ranked.len());
    ranked.truncate(n);
    ranked
}

/// The input as text: through `pigz -dc -p threads` when `pigz_threads > 0` (the child is returned so that its exit
/// status can be checked once the text is read), else as `open_text`.
fn open_input(path: &Path, pigz_threads: usize) -> Result<(Box<dyn BufRead>, Option<std::process::Child>)> {
    if pigz_threads == 0 || !path.extension().map_or(false, |x| x == "gz") {
        return Ok((open_text(path)?, None));
    }
    let mut child = std::process::Command::new("pigz")
        .args(["-dc", "-p", &pigz_threads.to_string()])
        .arg(path)
        .stdout(std::process::Stdio::piped())
        .spawn()
        .map_err(|e| format!("pigz: {}", e))?;
    let stdout = child.stdout.take().ok_or("pigz: no stdout")?;
    Ok((Box::new(BufReader::with_capacity(8 << 20, stdout)), Some(child)))
}

fn int(field: &str, what: &str) -> Result<i64> {
    field.parse::<i64>().map_err(|e| format!("bad {} {:?}: {}", what, field, e).into())
}

/// Calls `each` for every read of the token table read from `input`, in order, with the table's coordinate space;
/// returns that space. Checks rows grouped by read in ordinal order, contiguous token ordinals, one coordinate space
/// (q_symbol or base) and at least one read. No heap allocation per row (at most 32 columns).
pub fn for_each_read_from(mut input: Box<dyn BufRead>, mut each: impl FnMut(Read, &str) -> Result<()>) -> Result<String> {
    let mut line = String::new();
    if input.read_line(&mut line)? == 0 {
        return Err("empty token bundle".into());
    }
    let header = line.trim_end_matches(['\n', '\r']).to_string();
    let columns: HashMap<&str, usize> = header.split('\t').enumerate().map(|(i, c)| (c, i)).collect();
    let col = |name: &str| columns.get(name).copied().ok_or_else(|| format!("missing column {}", name));
    let (c_id, c_ord, c_len, c_tok) = (col("read_id")?, col("read_ordinal")?, col("source_length")?, col("token_ordinal")?);
    let (c_val, c_start, c_end, c_space) = (col("token_id")?, col("start")?, col("end")?, col("coordinate_space")?);
    let width = columns.len();
    if width > 32 {
        return Err("more than 32 columns".into());
    }
    let mut coordinate: Option<String> = None;
    let mut current: Option<Read> = None;
    let mut reads = 0usize;
    loop {
        line.clear();
        if input.read_line(&mut line)? == 0 {
            break;
        }
        let row = line.strip_suffix('\n').unwrap_or(&line);
        let row = row.strip_suffix('\r').unwrap_or(row);
        if row.contains('"') {
            return Err("quoted fields are not supported".into());
        }
        let mut f: [&str; 32] = [""; 32];
        let mut n = 0usize;
        for field in row.split('\t') {
            if n < 32 {
                f[n] = field;
            }
            n += 1;
        }
        if n != width {
            return Err(format!("row with {} fields, header has {}", n, width).into());
        }
        match &coordinate {
            None => coordinate = Some(f[c_space].to_string()),
            Some(c) if c != f[c_space] => return Err("mixed coordinate spaces".into()),
            _ => {}
        }
        let ordinal = int(f[c_ord], "read ordinal")? as usize;
        let same = current.as_ref().map_or(false, |r| r.id == f[c_id]);
        if !same {
            if ordinal != reads {
                return Err("generic token rows are not ordinal-major (or a read reappears)".into());
            }
            if let Some(read) = current.take() {
                each(read, coordinate.as_deref().unwrap())?;
            }
            current = Some(Read { id: f[c_id].to_string(), length: 0, tokens: Vec::new() });
            reads += 1;
        } else if ordinal != reads - 1 {
            return Err("generic token rows are not ordinal-major".into());
        }
        let read = current.as_mut().unwrap();
        read.length = int(f[c_len], "source length")?;
        if int(f[c_tok], "token ordinal")? as usize != read.tokens.len() {
            return Err("noncontiguous token ordinal".into());
        }
        let value = f[c_val].parse::<u32>().map_err(|e| format!("bad token id {:?}: {}", f[c_val], e))?;
        read.tokens.push(Token { value, start: int(f[c_start], "start")?, end: int(f[c_end], "end")? });
    }
    if let Some(read) = current.take() {
        each(read, coordinate.as_deref().unwrap())?;
    }
    match coordinate {
        Some(c) if reads > 0 && (c == "q_symbol" || c == "base") => Ok(c),
        _ => Err("empty/invalid generic token bundle".into()),
    }
}

/// Writes reads as token-table rows (renumbered; reads without tokens dropped).
pub struct TableWriter {
    out: BufWriter<File>,
    pub reads: u64,
    pub tokens: u64,
}

impl TableWriter {
    pub fn create(path: &Path) -> Result<TableWriter> {
        let mut out = BufWriter::with_capacity(1 << 20, File::create(path)?);
        writeln!(out, "{}", HEADER)?;
        Ok(TableWriter { out, reads: 0, tokens: 0 })
    }
    pub fn write(&mut self, id: &str, length: i64, tokens: &[Token], coordinate: &str) -> Result<()> {
        if tokens.is_empty() {
            return Ok(());
        }
        // the table has no quoting: an id with a tab, quote or line break cannot be written
        if id.contains(|c: char| c == '\t' || c == '"' || c == '\n' || c == '\r') {
            return Err(format!("read id {:?} would need quoting", id).into());
        }
        for (i, t) in tokens.iter().enumerate() {
            writeln!(self.out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", id, self.reads, length, i, t.value, t.start, t.end, coordinate)?;
        }
        self.reads += 1;
        self.tokens += tokens.len() as u64;
        Ok(())
    }
    pub fn finish(mut self) -> Result<(u64, u64)> {
        self.out.flush()?;
        Ok((self.reads, self.tokens))
    }
}

pub fn run() -> Result<()> {
    const USAGE: &str = "usage: squash __sparse sparse-repeat --tokens IN.tsv.gz --inverse-density N --fraction F --output TABLE.tsv \
[--pigz THREADS]";
    let mut argv = crate::argv::args().skip(1);
    let command = argv.next().ok_or(USAGE)?;
    if command != "sparse-repeat" {
        return Err(USAGE.into());
    }
    let mut map: HashMap<String, String> = HashMap::new();
    while let Some(key) = argv.next() {
        if !key.starts_with("--") {
            return Err(format!("unexpected argument {:?}\n{}", key, USAGE).into());
        }
        let value = argv.next().ok_or_else(|| format!("{} needs a value\n{}", key, USAGE))?;
        map.insert(key[2..].to_string(), value);
    }
    let mut need = |key: &str| map.remove(key).ok_or_else(|| format!("missing --{}\n{}", key, USAGE));
    let input = PathBuf::from(need("tokens")?);
    let output = PathBuf::from(need("output")?);
    if output.exists() {
        return Err(format!("refusing to overwrite {}", output.display()).into());
    }
    let n: u64 = need("inverse-density")?.parse()?;
    let fraction: f64 = need("fraction")?.parse()?;
    let pigz: usize = map.remove("pigz").map(|v| v.parse()).transpose()?.unwrap_or(0);
    if let Some(key) = map.keys().next() {
        return Err(format!("unknown option --{}\n{}", key, USAGE).into());
    }
    if n == 0 {
        return Err("--inverse-density must be at least 1".into());
    }
    if !(fraction > 0.0 && fraction < 1.0) {
        return Err("the repeat fraction must be in (0, 1)".into());
    }
    sparse_repeat(&input, n, fraction, pigz, &output)
}

/// `sparse-repeat`: the sparse subset of `input`'s tokens (1 in `n`), then the repeat cut (`fraction`) on it, into
/// `output`.
pub fn sparse_repeat(input: &Path, n: u64, fraction: f64, pigz: usize, output: &Path) -> Result<()> {
    let mut kept_reads: Vec<Read> = Vec::new();
    let (reader, child) = open_input(input, pigz)?;
    let coordinate = for_each_read_from(reader, |read, _| {
        let kept: Vec<Token> = read.tokens.iter().copied().filter(|t| keep_token(t.value, n)).collect();
        // the table has no quoting (as `TableWriter::write`, which writes the reads with tokens)
        if !kept.is_empty() && read.id.contains(|c: char| c == '\t' || c == '"' || c == '\n' || c == '\r') {
            return Err(format!("read id {:?} would need quoting", read.id).into());
        }
        if !kept.is_empty() {
            kept_reads.push(Read { id: read.id, length: read.length, tokens: kept });
        }
        Ok(())
    })?;
    if let Some(mut child) = child {
        let status = child.wait()?;
        if !status.success() {
            return Err(format!("pigz -dc {}: {}", input.display(), status).into());
        }
    }
    // the repeat cut on the sparse reads (kept_reads, in order)
    let mut counts: HashMap<u32, u64> = HashMap::new();
    for read in &kept_reads {
        for t in &read.tokens {
            *counts.entry(t.value).or_insert(0) += 1;
        }
    }
    let repetitive = repetitive_tokens(&counts, fraction);
    let drop: std::collections::HashSet<u32> = repetitive.iter().map(|&(v, _)| v).collect();
    let mut filtered = TableWriter::create(output)?;
    for read in &kept_reads {
        let kept: Vec<Token> = read.tokens.iter().copied().filter(|t| !drop.contains(&t.value)).collect();
        filtered.write(&read.id, read.length, &kept, &coordinate)?;
    }
    filtered.finish()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splitmix64_matches_its_reference_values() {
        // SplitMix64's first outputs for the seeds 0 and 1.
        assert_eq!(splitmix64(0), 0xE220_A839_7B1D_CDAF);
        assert_eq!(splitmix64(1), 0x910A_2DEC_8902_5CC1);
        assert!(keep_token(5, 1));
    }

    #[test]
    fn repetitive_tokens_rank_by_count_then_value_and_keep_at_least_one() {
        let counts: HashMap<u32, u64> = [(7, 4), (8, 3), (9, 3), (10, 1)].iter().copied().collect();
        assert_eq!(repetitive_tokens(&counts, 0.5), vec![(7, 4), (8, 3)]);
        assert_eq!(repetitive_tokens(&counts, 0.0001), vec![(7, 4)]);
    }

    /// `sparse-repeat`'s table equals the two filters applied one after the other to the reads (reads left without
    /// tokens dropped and the rest renumbered after each), in process and through pigz when it is installed; the reader
    /// gives the table's reads.
    #[test]
    fn sparse_repeat_equals_the_two_filters() {
        use flate2::write::GzEncoder;
        let dir = std::env::temp_dir().join(format!("select_fused_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let input = dir.join("in.tsv.gz");
        let mut text = format!("{}\n", HEADER);
        let mut reads: Vec<(String, i64, Vec<Token>)> = Vec::new();
        for r in 0..200u32 {
            if r % 17 == 3 || r % 23 == 0 {
                continue; // reads absent from the table do not exist; ordinals stay contiguous
            }
            let mut tokens = Vec::new();
            for i in 0..(r % 23) {
                let value = r.wrapping_mul(2654435761u32).wrapping_add(i * 40503) % 3000 + (i % 3) * 7;
                text.push_str(&format!("r{}\t{}\t{}\t{}\t{}\t{}\t{}\tq_symbol\n", r, reads.len(), 900 + r, i, value, 10 * i, 10 * i + 12));
                tokens.push(Token { value, start: 10 * i as i64, end: 10 * i as i64 + 12 });
            }
            reads.push((format!("r{}", r), 900 + r as i64, tokens));
        }
        let mut gz = GzEncoder::new(File::create(&input).unwrap(), flate2::Compression::default());
        gz.write_all(text.as_bytes()).unwrap();
        gz.finish().unwrap();
        // The two filters on the reads.
        let sparse: Vec<(&String, i64, Vec<Token>)> = reads.iter()
            .map(|(id, len, t)| (id, *len, t.iter().copied().filter(|t| keep_token(t.value, 3)).collect::<Vec<Token>>()))
            .filter(|r| !r.2.is_empty())
            .collect();
        let mut counts: HashMap<u32, u64> = HashMap::new();
        for r in &sparse {
            for t in &r.2 {
                *counts.entry(t.value).or_insert(0) += 1;
            }
        }
        let drop: std::collections::HashSet<u32> = repetitive_tokens(&counts, 0.01).iter().map(|&(v, _)| v).collect();
        assert!(!drop.is_empty());
        let mut want = format!("{}\n", HEADER);
        let mut ordinal = 0;
        for (id, len, tokens) in &sparse {
            let kept: Vec<&Token> = tokens.iter().filter(|t| !drop.contains(&t.value)).collect();
            if kept.is_empty() {
                continue;
            }
            for (i, t) in kept.iter().enumerate() {
                want.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\tq_symbol\n", id, ordinal, len, i, t.value, t.start, t.end));
            }
            ordinal += 1;
        }
        // Fused, in process and through pigz when it is installed.
        for pigz in [0usize, 2] {
            if pigz > 0 && std::process::Command::new("pigz").arg("--version").output().is_err() {
                continue;
            }
            let out = dir.join(format!("rf{}.tsv", pigz));
            sparse_repeat(&input, 3, 0.01, pigz, &out).unwrap();
            assert_eq!(std::fs::read_to_string(&out).unwrap(), want);
        }
        let mut got = Vec::new();
        for_each_read_from(open_text(&input).unwrap(), |r, c| {
            assert_eq!(c, "q_symbol");
            got.push((r.id, r.length, r.tokens));
            Ok(())
        }).unwrap();
        assert_eq!(got, reads);
        std::fs::remove_dir_all(&dir).ok();
    }
}
