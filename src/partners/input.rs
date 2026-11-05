//! The search's inputs: the uncorrected tokens and the excluded reads.

use super::*;
use crate::io::{READS_TABLE, TOKEN_SPANS, TOKEN_VALUES};

fn int(field: &str) -> Result<i64> {
    field.trim().parse::<i64>().map_err(|e| format!("bad integer {:?}: {}", field, e).into())
}

/// Loads a folder of token arrays -> (read ids, lengths in q-symbols, tokens), in catalog order: reads.tsv.gz (read_id,
/// read_ordinal, source_symbols, token_offset, token_count), tokens.u32le (one value per token) and token_spans.u64le
/// (start and end per token, q-symbols). The token and span arrays are read a read at a
/// time (seeking when a read's range does not follow the last one), so the whole files are never held in memory.
pub fn load_prepared(dir: &Path) -> Result<(Vec<String>, Vec<i64>, Vec<Vec<Token>>)> {
    use std::io::{Seek, SeekFrom};
    let mut lines = open_text(&dir.join(READS_TABLE))?.lines();
    let header = lines.next().ok_or("empty prepared catalog")??;
    let columns: HashMap<String, usize> = header.split('\t').enumerate().map(|(i, c)| (c.to_string(), i)).collect();
    let col = |name: &str| columns.get(name).copied().ok_or_else(|| format!("missing column {}", name));
    let (c_id, c_ord, c_len) = (col("read_id")?, col("read_ordinal")?, col("source_symbols")?);
    let (c_off, c_count) = (col("token_offset")?, col("token_count")?);
    let (token_path, span_path) = (dir.join(TOKEN_VALUES), dir.join(TOKEN_SPANS));
    let (token_bytes, span_bytes) = (std::fs::metadata(&token_path)?.len() as usize, std::fs::metadata(&span_path)?.len() as usize);
    if token_bytes % 4 != 0 || span_bytes % 16 != 0 || token_bytes / 4 != span_bytes / 16 {
        return Err("malformed prepared arrays".into());
    }
    let n_tokens = token_bytes / 4;
    let mut token_file = BufReader::with_capacity(1 << 22, File::open(&token_path)?);
    let mut span_file = BufReader::with_capacity(1 << 24, File::open(&span_path)?);
    let (mut at, mut raw_tokens, mut raw_spans) = (0usize, Vec::new(), Vec::new());
    let (mut ids, mut lengths, mut tokens) = (Vec::new(), Vec::new(), Vec::new());
    for (expected, line) in lines.enumerate() {
        let line = line?;
        let row: Vec<&str> = line.split('\t').collect();
        if int(row[c_ord])? != expected as i64 {
            return Err("noncontiguous prepared ordinals".into());
        }
        let (offset, count) = (int(row[c_off])? as usize, int(row[c_count])? as usize);
        if offset + count > n_tokens {
            return Err("prepared token range outside arrays".into());
        }
        if offset != at {
            token_file.seek(SeekFrom::Start(4 * offset as u64))?;
            span_file.seek(SeekFrom::Start(16 * offset as u64))?;
        }
        raw_tokens.resize(4 * count, 0);
        raw_spans.resize(16 * count, 0);
        token_file.read_exact(&mut raw_tokens)?;
        span_file.read_exact(&mut raw_spans)?;
        at = offset + count;
        let value = |i: usize| u32::from_le_bytes(raw_tokens[4 * i..4 * i + 4].try_into().unwrap());
        let span = |i: usize| u64::from_le_bytes(raw_spans[8 * i..8 * i + 8].try_into().unwrap()) as i64;
        ids.push(row[c_id].to_string());
        lengths.push(int(row[c_len])?);
        tokens.push((0..count).map(|i| Token::new(value(i), span(2 * i), span(2 * i + 1))).collect());
    }
    Ok((ids, lengths, tokens))
}


/// The excluded reads (`--exclude-reads`): a table whose first column is read_id (header first), as a flag
/// per read ordinal of `ids`. Every listed read must be in the catalog, once.
pub fn load_excluded(path: &Path, ids: &[String]) -> Result<Vec<bool>> {
    let ordinal: HashMap<&str, usize> = ids.iter().enumerate().map(|(i, x)| (x.as_str(), i)).collect();
    let mut lines = open_text(path)?.lines();
    let header = lines.next().ok_or("empty excluded-read table")??;
    if header.split('\t').next() != Some("read_id") {
        return Err(format!("{}: the first column must be read_id", path.display()).into());
    }
    let mut excluded = vec![false; ids.len()];
    for line in lines {
        let line = line?;
        let id = line.split('\t').next().unwrap_or("");
        let &r = ordinal.get(id).ok_or_else(|| format!("{}: read {:?} is not in the catalog", path.display(), id))?;
        if excluded[r] {
            return Err(format!("{}: read {} listed twice", path.display(), id).into());
        }
        excluded[r] = true;
    }
    Ok(excluded)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The prepared arrays read a read at a time give the tokens of the whole-file reader, with read ranges out of
    /// order in the arrays (a seek per read) and in order.
    #[test]
    fn streamed_prepared_arrays_equal_the_whole_file_reader() {
        use flate2::write::GzEncoder;
        let dir = std::env::temp_dir().join(format!("prepared_stream_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (_, tokens) = noisy_reads(10, 8);
        for reversed in [false, true] {
            let prefix = dir.join(if reversed { "rev" } else { "fwd" });
            std::fs::create_dir_all(&prefix).unwrap();
            let order: Vec<usize> = if reversed { (0..tokens.len()).rev().collect() } else { (0..tokens.len()).collect() };
            let (mut values, mut spans, mut offset) = (Vec::new(), Vec::new(), vec![0usize; tokens.len()]);
            for &r in &order {
                offset[r] = values.len() / 4;
                for t in &tokens[r] {
                    values.extend(t.value.to_le_bytes());
                    spans.extend((t.s() as u64).to_le_bytes());
                    spans.extend((t.e() as u64).to_le_bytes());
                }
            }
            std::fs::write(prefix.join(TOKEN_VALUES), &values).unwrap();
            std::fs::write(prefix.join(TOKEN_SPANS), &spans).unwrap();
            let mut gz = GzEncoder::new(File::create(prefix.join(READS_TABLE)).unwrap(), flate2::Compression::fast());
            writeln!(gz, "read_id\tread_ordinal\tsource_symbols\ttoken_offset\ttoken_count").unwrap();
            for (r, read) in tokens.iter().enumerate() {
                writeln!(gz, "read{}\t{}\t{}\t{}\t{}", r, r, 10 * read.len() + 3, offset[r], read.len()).unwrap();
            }
            gz.finish().unwrap();
            let (ids, lengths, back) = load_prepared(&prefix).unwrap();
            assert_eq!(back, tokens);
            assert!(back.iter().all(|r| r.capacity() == r.len()));
            assert_eq!(ids, (0..tokens.len()).map(|r| format!("read{}", r)).collect::<Vec<_>>());
            assert_eq!(lengths, tokens.iter().map(|r| 10 * r.len() as i64 + 3).collect::<Vec<_>>());
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
