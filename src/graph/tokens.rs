//! The token arrays: `prepare` (the q-symbols sketched into tokens), `adapt-tokens` (a tokens table as token arrays),
//! and the reader of both.

use crate::events::{self, RawParams, RawSketchToken, RawSymbols};
use crate::io::{READS_TABLE, TOKEN_SPANS, TOKEN_VALUES, open_text};
use rayon::prelude::*;
use structopt::StructOpt;
use std::collections::{BTreeMap, HashSet};
use std::convert::{TryFrom, TryInto};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};

use super::caps::stop_on_cap;
use super::files::{
    BackgroundGz, create_plain, csv_field, finish_plain, invalid_data, invalid_input, require_threads, u64_usize,
    write_u32, write_u64,
};

const QSYMBOLS_TAG: &str = crate::io::QSYMBOLS_TAG;

const QSYMBOLS_COLUMNS: &str = "#columns=read_id\tsymbols";

/// Parse a generic `tokens.tsv.gz` and write the prepared token arrays and their reads table.
#[derive(Debug, StructOpt)]
pub(super) struct AdaptTokensOpt {
    #[structopt(long, parse(from_os_str))]
    tokens: PathBuf,
    /// Writes DIR/tokens.u32le, DIR/token_spans.u64le and DIR/reads.tsv.gz (read_id, read_ordinal, source_symbols,
    /// token_offset, token_count).
    #[structopt(long, parse(from_os_str))]
    output_dir: PathBuf,
}

#[derive(Debug, StructOpt)]
pub(super) struct PrepareOpt {
    /// The q-symbols table: `#squash-qsymbols-v1`, the quantizer line, `#columns=read_id\tsymbols`, then one read per
    /// row.
    #[structopt(long, parse(from_os_str))]
    q_symbols: PathBuf,
    /// Q-symbols per token.
    #[structopt(long)]
    signal_events: u32,
    /// Minimizer window: one token per window of this many hashes.
    #[structopt(long)]
    signal_window: u32,
    /// Writes DIR/tokens.u32le, DIR/token_spans.u64le and DIR/reads.tsv.gz.
    #[structopt(long, parse(from_os_str))]
    output_dir: PathBuf,
    #[structopt(long, default_value = "1")]
    pub(super) threads: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct TokenRead {
    id: String,
    pub(super) ordinal: u64,
    source_symbols: u64,
    pub(super) tokens: Vec<RawSketchToken>,
}

/// One read's reads-table row, tokens and spans at token `offset`; returns the next offset.
fn write_prepared_read<T: Write, S: Write, R: Write>(tw: &mut T, sw: &mut S, rw: &mut R, read: &TokenRead, offset: u64) -> io::Result<u64> {
    writeln!(
        rw,
        "{}\t{}\t{}\t{}\t{}",
        read.id,
        read.ordinal,
        read.source_symbols,
        offset,
        read.tokens.len()
    )?;
    for token in &read.tokens {
        let hash = u32::try_from(token.hash)
            .map_err(|_| invalid_data("raw token exceeds u32 packing contract"))?;
        write_u32(tw, hash)?;
        write_u64(sw, token.start as u64)?;
        write_u64(sw, token.end as u64)?;
    }
    offset
        .checked_add(read.tokens.len() as u64)
        .ok_or_else(|| invalid_data("token offset overflow"))
}

pub(super) fn load_prepared(dir: &Path) -> io::Result<Vec<TokenRead>> {
    let reads_path = dir.join(READS_TABLE);
    let tokens_path = dir.join(TOKEN_VALUES);
    let spans_path = dir.join(TOKEN_SPANS);
    let mut tr = BufReader::with_capacity(1 << 20, File::open(tokens_path)?);
    let mut sr = BufReader::with_capacity(1 << 20, File::open(spans_path)?);
    let mut rr = open_text(&reads_path)?;
    let mut line = String::new();
    rr.read_line(&mut line)?;
    if line.trim_end_matches(['\r', '\n'])
        != "read_id\tread_ordinal\tsource_symbols\ttoken_offset\ttoken_count"
    {
        return Err(invalid_data("malformed read index header"));
    }
    let mut result = Vec::new();
    let mut expected_offset = 0u64;
    let mut token_bytes: Vec<u8> = Vec::new();
    let mut span_bytes: Vec<u8> = Vec::new();
    for (row_number, row) in rr.lines().enumerate() {
        let row = row?;
        let f: Vec<&str> = row.split('\t').collect();
        if f.len() != 5 {
            return Err(invalid_data(format!(
                "malformed read index row {}",
                row_number + 2
            )));
        }
        let ordinal = f[1]
            .parse::<u64>()
            .map_err(|_| invalid_data("invalid read ordinal"))?;
        let source_symbols = f[2]
            .parse::<u64>()
            .map_err(|_| invalid_data("invalid source length"))?;
        let offset = f[3]
            .parse::<u64>()
            .map_err(|_| invalid_data("invalid token offset"))?;
        let count = f[4]
            .parse::<u64>()
            .map_err(|_| invalid_data("invalid token count"))?;
        if ordinal != result.len() as u64 || offset != expected_offset || f[0].is_empty() {
            return Err(invalid_data("noncontiguous read index"));
        }
        // Bulk-read both arrays for the whole read instead of three small reads per token (this runs at every process
        // start, over tens of millions of tokens).
        let n = u64_usize(count, "token count")?;
        let mut tokens = Vec::with_capacity(n);
        token_bytes.resize(n * 4, 0);
        span_bytes.resize(n * 16, 0);
        tr.read_exact(&mut token_bytes)?;
        sr.read_exact(&mut span_bytes)?;
        let mut previous = 0usize;
        for (i, (token_chunk, span_chunk)) in token_bytes
            .chunks_exact(4)
            .zip(span_bytes.chunks_exact(16))
            .enumerate()
        {
            let hash = u32::from_le_bytes(token_chunk.try_into().unwrap()) as u64;
            let start = u64_usize(
                u64::from_le_bytes(span_chunk[..8].try_into().unwrap()),
                "token start",
            )?;
            let end = u64_usize(
                u64::from_le_bytes(span_chunk[8..].try_into().unwrap()),
                "token end",
            )?;
            if start >= end || end > source_symbols as usize || (i > 0 && start < previous) {
                return Err(invalid_data("invalid prepared token span"));
            }
            previous = start;
            tokens.push(RawSketchToken { hash, start, end });
        }
        result.push(TokenRead {
            id: f[0].to_owned(),
            ordinal,
            source_symbols,
            tokens,
        });
        expected_offset = expected_offset
            .checked_add(count)
            .ok_or_else(|| invalid_data("token offset overflow"))?;
    }
    let mut extra = [0u8; 1];
    if tr.read(&mut extra)? != 0 || sr.read(&mut extra)? != 0 {
        return Err(invalid_data("trailing bytes in prepared token arrays"));
    }
    if result.is_empty() {
        return Err(invalid_data("the prepared tokens are empty"));
    }
    Ok(result)
}

/// The sketch parameters of `graph prepare`: `events` q-symbols per token and a minimizer window of `window`, with the
/// quantizer of the q-symbols table's header line (`quant_bit`, `n_buckets` and the fine range, as bit patterns).
fn corrected_params(fingerprint: &str, events: u32, window: u32) -> io::Result<RawParams> {
    let mut fields = BTreeMap::new();
    for item in fingerprint.trim_start_matches('#').split_whitespace() {
        let (k, v) = item
            .split_once('=')
            .ok_or_else(|| invalid_data("malformed corrected-symbol fingerprint"))?;
        fields.insert(k, v);
    }
    let get = |k: &str| {
        fields
            .get(k)
            .copied()
            .ok_or_else(|| invalid_data(format!("corrected-symbol fingerprint lacks {}", k)))
    };
    let quant = get("quant_bit")?
        .parse::<u32>()
        .map_err(|_| invalid_data("bad quant_bit"))?;
    let buckets = get("n_buckets")?
        .parse::<u32>()
        .map_err(|_| invalid_data("bad n_buckets"))?;
    // The sketcher's own budgets, checked here rather than inside a worker thread.
    if events == 0 || events * quant > 58 {
        return Err(invalid_input("--signal-events exceeds the 58-bit packing budget"));
    }
    if window == 0 || window >= 256 {
        return Err(invalid_input("--signal-window must be in 1..255"));
    }
    let parse_bits = |k: &str| -> io::Result<f32> {
        let bits =
            u32::from_str_radix(get(k)?, 16).map_err(|_| invalid_data(format!("bad {}", k)))?;
        Ok(f32::from_bits(bits))
    };
    Ok(RawParams {
        w: window,
        e: events,
        quant_bit: quant,
        n_buckets: buckets,
        bucket_diff: 0,
        sig_diff: -1.0,
        fine_min: parse_bits("fine_min_bits")?,
        fine_max: parse_bits("fine_max_bits")?,
        fine_range: parse_bits("fine_range_bits")?,
        ..RawParams::default()
    })
}

/// One q-symbols row, parsed (`prepare_streamed`).
enum ParsedCorrected {
    Malformed,
    Row(String, Result<Vec<u8>, &'static str>),
}

fn parse_corrected_row(row: &str, n_buckets: u32) -> ParsedCorrected {
    match row.split_once('\t') {
        None => ParsedCorrected::Malformed,
        Some((id, field)) => {
            let mut values = Vec::with_capacity(field.len() / 2 + 1);
            let mut error = None;
            if !field.is_empty() {
                for value in field.split(',') {
                    match value.parse::<u8>() {
                        Err(_) => {
                            error = Some("invalid corrected q-symbol");
                            break;
                        }
                        Ok(symbol) if symbol as u32 >= n_buckets => {
                            error = Some("corrected q-symbol outside the quantizer's alphabet");
                            break;
                        }
                        Ok(symbol) => values.push(symbol),
                    }
                }
            }
            ParsedCorrected::Row(id.to_owned(), error.map_or(Ok(values), Err))
        }
    }
}

/// The corrected-symbol file's three header lines: the quantizer fingerprint.
fn read_corrected_header<R: BufRead>(reader: &mut R) -> io::Result<String> {
    let mut line = String::new();
    reader.read_line(&mut line)?;
    if line.trim_end_matches(['\r', '\n']) != QSYMBOLS_TAG {
        return Err(invalid_data("missing corrected-symbol magic"));
    }
    line.clear();
    reader.read_line(&mut line)?;
    let fingerprint = line.trim_end_matches(['\r', '\n']).to_string();
    line.clear();
    reader.read_line(&mut line)?;
    if line.trim_end_matches(['\r', '\n']) != QSYMBOLS_COLUMNS {
        return Err(invalid_data("malformed corrected-symbol header"));
    }
    Ok(fingerprint)
}

/// `graph prepare`: rows read per batch (at most), and text bytes per batch (at least one row).
const PREPARE_BATCH_ROWS: usize = 16384;

const PREPARE_BATCH_BYTES: usize = 256 << 20;

/// The q-symbols read, checked, sketched and written a batch of reads at a time (rows parsed and reads sketched on the
/// pool, written in file order), into the token array, the span array and the reads table. Returns the token count. A
/// bad row stops the run when its batch is read.
fn prepare_streamed(
    path: &Path,
    events: u32,
    window: u32,
    tokens_path: &Path,
    spans_path: &Path,
    reads_path: &Path,
) -> io::Result<u64> {
    let mut reader = BufReader::with_capacity(1 << 20, File::open(path)?);
    let fingerprint = read_corrected_header(&mut reader)?;
    let params = corrected_params(&fingerprint, events, window)?;
    let mut tw = create_plain(tokens_path)?;
    let mut sw = create_plain(spans_path)?;
    let mut rw = BackgroundGz::create(reads_path)?;
    writeln!(rw, "read_id\tread_ordinal\tsource_symbols\ttoken_offset\ttoken_count")?;
    let mut seen = HashSet::new();
    let (mut ordinal, mut offset) = (0u64, 0u64);
    let mut lines = reader.lines();
    loop {
        let (mut batch, mut bytes) = (Vec::new(), 0usize);
        while batch.len() < PREPARE_BATCH_ROWS && bytes < PREPARE_BATCH_BYTES {
            match lines.next() {
                None => break,
                Some(row) => {
                    let row = row?;
                    bytes += row.len();
                    batch.push(row);
                }
            }
        }
        if batch.is_empty() {
            break;
        }
        let parsed: Vec<ParsedCorrected> = batch.par_iter().map(|row| parse_corrected_row(row, params.n_buckets)).collect();
        drop(batch);
        let mut raw = Vec::with_capacity(parsed.len());
        for row in parsed {
            match row {
                ParsedCorrected::Malformed => {
                    return Err(invalid_data(format!("malformed corrected-symbol row {}", ordinal + 4)));
                }
                ParsedCorrected::Row(id, values) => {
                    if id.is_empty() || !seen.insert(id.clone()) {
                        return Err(invalid_data("empty or duplicate corrected-symbol read ID"));
                    }
                    let values = values.map_err(invalid_data)?;
                    raw.push((id, ordinal, values));
                }
            }
            ordinal += 1;
        }
        let reads: Vec<TokenRead> = raw.into_par_iter().map(|(id, o, values)| sketch_corrected_read(id, o, values, &params)).collect();
        for read in &reads {
            offset = write_prepared_read(&mut tw, &mut sw, &mut rw, read, offset)?;
        }
    }
    if ordinal == 0 {
        return Err(invalid_data("corrected-symbol input is empty"));
    }
    finish_plain(tw)?;
    finish_plain(sw)?;
    rw.finish()?;
    Ok(offset)
}

/// One read's tokens from its q-symbols (`prepare_streamed`).
fn sketch_corrected_read(id: String, ordinal: u64, values: Vec<u8>, params: &RawParams) -> TokenRead {
    let source_symbols = values.len() as u64;
    let symbols = RawSymbols::from_values(values);
    let filtered = events::filter_symbol_buckets(&symbols, params);
    let tokens = events::sketch_symbols_with_spans(&filtered, params);
    TokenRead { id, ordinal, source_symbols, tokens }
}

/// `graph prepare`: the q-symbols sketched into tokens, a batch of reads at a time (`prepare_streamed`). The two raw
/// arrays take twenty bytes per token: more than `PREPARE_MAX_SCRATCH_BYTES` of them, or of all three files, is a stop.
pub(super) fn prepare(opt: PrepareOpt) -> io::Result<()> {
    require_threads(opt.threads)?;
    std::fs::create_dir_all(&opt.output_dir)?;
    let tokens_path = opt.output_dir.join(TOKEN_VALUES);
    let spans_path = opt.output_dir.join(TOKEN_SPANS);
    let reads_path = opt.output_dir.join(READS_TABLE);
    let tokens = prepare_streamed(&opt.q_symbols, opt.signal_events, opt.signal_window, &tokens_path, &spans_path,
                                  &reads_path)?;
    let minimum_array_bytes = tokens
        .checked_mul(20)
        .ok_or_else(|| invalid_data("prepared-array byte estimate overflow"))?;
    if minimum_array_bytes > PREPARE_MAX_SCRATCH_BYTES {
        return stop_on_cap(format!("minimum_scratch_bytes={} exceeds {}", minimum_array_bytes, PREPARE_MAX_SCRATCH_BYTES));
    }
    let mut scratch = 0u64;
    for path in [&tokens_path, &spans_path, &reads_path] {
        scratch = scratch
            .checked_add(path.metadata()?.len())
            .ok_or_else(|| invalid_data("scratch-byte total overflow"))?;
    }
    if scratch > PREPARE_MAX_SCRATCH_BYTES {
        return stop_on_cap(format!("scratch_bytes={} exceeds {}", scratch, PREPARE_MAX_SCRATCH_BYTES));
    }
    Ok(())
}

/// `graph prepare`'s resource cap (100 GiB).
const PREPARE_MAX_SCRATCH_BYTES: u64 = 107374182400;

pub(super) fn adapt_tokens(opt: AdaptTokensOpt) -> io::Result<()> {
    let mut reader = open_text(&opt.tokens)?;
    let mut header = String::new();
    reader.read_line(&mut header)?;
    let columns: Vec<&str> = header.trim_end_matches(['\r', '\n']).split('\t').collect();
    let column = |name: &str| -> io::Result<usize> {
        columns
            .iter()
            .position(|c| *c == name)
            .ok_or_else(|| invalid_data(format!("generic token table lacks {}", name)))
    };
    let (c_id, c_ordinal, c_length, c_token_ordinal, c_value, c_start, c_end, c_space) = (
        column("read_id")?,
        column("read_ordinal")?,
        column("source_length")?,
        column("token_ordinal")?,
        column("token_id")?,
        column("start")?,
        column("end")?,
        column("coordinate_space")?,
    );
    let width = columns.len();
    std::fs::create_dir_all(&opt.output_dir)?;
    let mut tw = create_plain(&opt.output_dir.join(TOKEN_VALUES))?;
    let mut sw = create_plain(&opt.output_dir.join(TOKEN_SPANS))?;
    // (read_id, source_length of its last row, token count)
    let mut reads: Vec<(String, String, u64)> = Vec::new();
    let mut coordinate: Option<String> = None;
    let mut line = String::new();
    let number = |text: &str, what: &str| -> io::Result<u64> {
        text.parse::<u64>()
            .map_err(|_| invalid_data(format!("generic token table has a non-numeric {}", what)))
    };
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let row = line.trim_end_matches(['\r', '\n']);
        let f: Vec<&str> = row.split('\t').collect();
        if f.len() != width {
            return Err(invalid_data("malformed generic token row"));
        }
        let ordinal = number(f[c_ordinal], "read_ordinal")? as usize;
        let rid = f[c_id];
        if ordinal == reads.len() {
            reads.push((rid.to_owned(), String::new(), 0));
        } else if ordinal + 1 != reads.len() || reads[ordinal].0 != rid {
            // A read's rows must be contiguous: prepared tables are always read-major.
            return Err(invalid_data("generic token rows are not ordinal-major"));
        }
        let space = f[c_space];
        match &coordinate {
            None => coordinate = Some(space.to_owned()),
            Some(seen) if seen != space => return Err(invalid_data("mixed coordinate spaces")),
            _ => {}
        }
        let read = reads.last_mut().unwrap();
        // The last source_length seen per read is kept.
        read.1.clear();
        read.1.push_str(&number(f[c_length], "source_length")?.to_string());
        if number(f[c_token_ordinal], "token_ordinal")? != read.2 {
            return Err(invalid_data("noncontiguous token ordinal"));
        }
        read.2 += 1;
        let value = number(f[c_value], "token_id")?;
        let value = u32::try_from(value).map_err(|_| invalid_data("token exceeds u32"))?;
        write_u32(&mut tw, value)?;
        write_u64(&mut sw, number(f[c_start], "start")?)?;
        write_u64(&mut sw, number(f[c_end], "end")?)?;
    }
    match coordinate.as_deref() {
        Some("q_symbol" | "base") if !reads.is_empty() => {}
        _ => return Err(invalid_data("empty/invalid generic token bundle")),
    }
    finish_plain(tw)?;
    finish_plain(sw)?;
    let mut rw = BackgroundGz::create(&opt.output_dir.join(READS_TABLE))?;
    writeln!(rw, "read_id\tread_ordinal\tsource_symbols\ttoken_offset\ttoken_count")?;
    let mut offset = 0u64;
    for (ordinal, (rid, length, count)) in reads.iter().enumerate() {
        writeln!(rw, "{}\t{}\t{}\t{}\t{}", csv_field(rid), ordinal, length, offset, count)?;
        offset += count;
    }
    rw.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn n_buckets_is_explicitly_eight() {
        let fp = "#quant_bit=3 n_buckets=8 fine_min_bits=bfc00000 fine_max_bits=3fc00000 fine_range_bits=3f400000";
        let params = corrected_params(fp, 7, 2).unwrap();
        assert_eq!((params.n_buckets, params.quant_bit, params.e, params.w), (8, 3, 7, 2));
        assert_ne!(params.n_buckets, RawParams::default().n_buckets);
        // the sketcher's budgets: 7 x 3 bits fit in 58, 20 x 3 do not; the window is under 256
        assert!(corrected_params(fp, 20, 2).is_err());
        assert!(corrected_params(fp, 7, 256).is_err());
    }
}
