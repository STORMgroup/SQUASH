//! File helpers shared by the steps: text input (plain or gzip), the binary token and span arrays, header-addressed
//! TSV rows, and floats written as Python writes them. Errors name the file.

use flate2::read::MultiGzDecoder;
use std::convert::TryInto;
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// The steps' result type: any error, with a message.
pub type Result<T> = std::result::Result<T, Box<dyn Error + Send + Sync>>;

/// A set of reads as tokens is a folder of three files: the reads table (read_id, read_ordinal, source_symbols,
/// token_offset, token_count; gzip), every token's value (u32), and every token's q-symbol span (start, end; u64 pairs).
pub const READS_TABLE: &str = "reads.tsv.gz";
pub const TOKEN_VALUES: &str = "tokens.u32le";
pub const TOKEN_SPANS: &str = "token_spans.u64le";

/// The first line of a q-symbols table.
pub const QSYMBOLS_TAG: &str = "#squash-qsymbols-v1";

/// The error `e` of file `path`, with the path in its message.
pub fn file_error(path: &Path, e: impl std::fmt::Display) -> io::Error {
    io::Error::new(io::ErrorKind::Other, format!("{}: {}", path.display(), e))
}

/// A text file as lines, gunzipped when its name ends in `.gz` (1-MiB buffer).
pub fn open_text(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|e| file_error(path, e))?;
    Ok(if path.extension().map_or(false, |x| x == "gz") {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, file))
    })
}

/// A file of little-endian u32 values (the token arrays).
pub fn read_u32le(path: &Path) -> io::Result<Vec<u32>> {
    let raw = fs::read(path).map_err(|e| file_error(path, e))?;
    if raw.len() % 4 != 0 {
        return Err(file_error(path, "not whole u32 values"));
    }
    Ok(raw.chunks_exact(4).map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]])).collect())
}

/// A file of little-endian u64 pairs (the token spans: each token's start and end q-symbol).
pub fn read_u64_pairs(path: &Path) -> io::Result<Vec<[u64; 2]>> {
    let raw = fs::read(path).map_err(|e| file_error(path, e))?;
    if raw.len() % 16 != 0 {
        return Err(file_error(path, "not whole (start, end) u64 pairs"));
    }
    Ok(raw
        .chunks_exact(16)
        .map(|b| [u64::from_le_bytes(b[..8].try_into().unwrap()), u64::from_le_bytes(b[8..].try_into().unwrap())])
        .collect())
}

/// Calls `each` with every row of a TSV with a header (gzip or plain), as the fields of `columns` in that order. A
/// missing column or a short row is an error.
pub fn tsv_rows(path: &Path, columns: &[&str], mut each: impl FnMut(&[&str]) -> Result<()>) -> Result<()> {
    let mut lines = open_text(path)?.lines();
    let head = lines.next().ok_or_else(|| file_error(path, "empty"))?.map_err(|e| file_error(path, e))?;
    let head: Vec<&str> = head.split('\t').collect();
    let at: Vec<usize> = columns
        .iter()
        .map(|c| head.iter().position(|h| h == c).ok_or_else(|| file_error(path, format!("no column {}", c))))
        .collect::<io::Result<_>>()?;
    for line in lines {
        let line = line.map_err(|e| file_error(path, e))?;
        let f: Vec<&str> = line.split('\t').collect();
        let mut row: Vec<&str> = Vec::with_capacity(at.len());
        for &i in &at {
            row.push(f.get(i).copied().ok_or_else(|| file_error(path, "short row"))?);
        }
        each(&row)?;
    }
    Ok(())
}

/// `x` as Python's `str(float)` writes it: the shortest digits that read back as `x`, a `.0` on whole numbers, and
/// the exponent form (`1e-05`, `1e+16`) below 1e-4 and from 1e16 on.
pub fn py_float(x: f64) -> String {
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf" } else { "-inf" }.into();
    }
    let a = x.abs();
    if a != 0.0 && !(1e-4..1e16).contains(&a) {
        let e = format!("{:e}", x);
        let (mantissa, exp) = e.split_once('e').unwrap();
        let exp: i32 = exp.parse().unwrap();
        return format!("{}e{}{:02}", mantissa, if exp < 0 { '-' } else { '+' }, exp.abs());
    }
    let s = format!("{}", x);
    if s.contains('.') { s } else { s + ".0" }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn py_float_writes_as_python() {
        for (x, s) in [(0.0, "0.0"), (-0.0, "-0.0"), (1.0, "1.0"), (0.1, "0.1"), (1.5e-5, "1.5e-05"), (1e-4, "0.0001"),
                       (123456.0, "123456.0"), (1e16, "1e+16"), (9999999999999998.0, "9999999999999998.0"),
                       (2.5e300, "2.5e+300"), (-3.25e-7, "-3.25e-07"), (f64::NAN, "nan"), (f64::INFINITY, "inf"),
                       (f64::NEG_INFINITY, "-inf")] {
            assert_eq!(py_float(x), s, "{}", x);
        }
    }

    #[test]
    fn arrays_and_rows_read_back() {
        let dir = std::env::temp_dir().join(format!("squash_io_test_{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let u = dir.join("a.u32le");
        fs::write(&u, [1u8, 0, 0, 0, 2, 0, 0, 0]).unwrap();
        assert_eq!(read_u32le(&u).unwrap(), vec![1, 2]);
        fs::write(&u, [1u8, 0, 0]).unwrap();
        assert!(read_u32le(&u).is_err());
        let p = dir.join("b.u64le");
        let raw: Vec<u8> = [3u64, 9].iter().flat_map(|v| v.to_le_bytes()).collect();
        fs::write(&p, raw).unwrap();
        assert_eq!(read_u64_pairs(&p).unwrap(), vec![[3, 9]]);
        let t = dir.join("c.tsv.gz");
        let mut gz = flate2::write::GzEncoder::new(File::create(&t).unwrap(), flate2::Compression::default());
        std::io::Write::write_all(&mut gz, b"a\tb\tc\n1\t2\t3\n4\t5\t6\n").unwrap();
        gz.finish().unwrap();
        let mut got = Vec::new();
        tsv_rows(&t, &["c", "a"], |r| {
            got.push(format!("{}{}", r[0], r[1]));
            Ok(())
        })
        .unwrap();
        assert_eq!(got, ["31", "64"]);
        assert!(tsv_rows(&t, &["d"], |_| Ok(())).is_err());
        fs::remove_dir_all(&dir).unwrap();
    }
}
