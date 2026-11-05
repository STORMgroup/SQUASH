//! The round's token table (formatted in parallel blocks, uncompressed or through pigz) and the final round's
//! sparse table made in memory (`--sparse-last`).

use super::*;

/// The token table of the corrected reads (`tokens.tsv`).
pub(super) fn write_token_table<W: Write>(out: &mut W, ids: &[String], lengths: &[i64], corrected: &[Vec<Token>]) -> Result<()> {
    write_table_blocks(out, corrected.len(), |buf, r| format_rows(buf, &ids[r], r, lengths[r], &corrected[r]))
}

/// The header, then `row(buf, i)` for i in 0..n, in order: the rows are formatted in blocks of reads on a pool of
/// their own (`table_pool`) and written in read order, so the text does not depend on the pool's size.
fn write_table_blocks<W: Write>(out: &mut W, n: usize, row: impl Fn(&mut Vec<u8>, usize) + Sync) -> Result<()> {
    writeln!(out, "read_id\tread_ordinal\tsource_length\ttoken_ordinal\ttoken_id\tstart\tend\tcoordinate_space")?;
    const BLOCK: usize = 64; // reads per formatted block
    const BATCH: usize = 32; // blocks formatted at once (~200 MB of text with long reads)
    let pool = table_pool();
    let format_batch = |r0: usize| -> Vec<Vec<u8>> {
        let r1 = (r0 + BLOCK * BATCH).min(n);
        let starts: Vec<usize> = (r0..r1).step_by(BLOCK).collect();
        pool.install(|| starts.par_iter().map(|&a| {
            let mut buf = Vec::with_capacity(1 << 20);
            for r in a..(a + BLOCK).min(r1) {
                row(&mut buf, r);
            }
            buf
        }).collect())
    };
    // The next batch is formatted while this one is written (two batches at most in flight), so formatting and the
    // sink (pigz) overlap; the same bytes in the same order.
    std::thread::scope(|scope| -> Result<()> {
        let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<Vec<u8>>>(1);
        let format_batch = &format_batch;
        let producer = scope.spawn(move || {
            let mut r0 = 0;
            while r0 < n {
                if tx.send(format_batch(r0)).is_err() {
                    break; // the writer stopped (an error)
                }
                r0 += BLOCK * BATCH;
            }
        });
        let mut result = Ok(());
        for blocks in rx.iter() {
            if result.is_ok() {
                for b in &blocks {
                    if let Err(e) = out.write_all(b) {
                        result = Err(e);
                        break;
                    }
                }
            }
            if result.is_err() {
                break;
            }
        }
        drop(rx);
        producer.join().map_err(|_| "token table formatter panicked")?;
        result?;
        Ok(())
    })?;
    out.flush()?;
    Ok(())
}

/// `write_token_table`'s threads.
const TABLE_THREADS: usize = 16;

/// `write_token_table`'s pool of `TABLE_THREADS` threads.
fn table_pool() -> &'static rayon::ThreadPool {
    static POOL: std::sync::OnceLock<rayon::ThreadPool> = std::sync::OnceLock::new();
    POOL.get_or_init(|| {
        rayon::ThreadPoolBuilder::new().num_threads(TABLE_THREADS).thread_name(|i| format!("table-{}", i)).build()
            .expect("the token table's thread pool")
    })
}

/// One read's rows of the token table, as `writeln!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\tq_symbol", id, r, length, k,
/// value, start, end)` writes them.
fn format_rows(buf: &mut Vec<u8>, id: &str, r: usize, length: i64, read: &[Token]) {
    let mut prefix = Vec::with_capacity(id.len() + 32);
    prefix.extend_from_slice(id.as_bytes());
    prefix.push(b'\t');
    push_int(&mut prefix, r as i64);
    prefix.push(b'\t');
    push_int(&mut prefix, length);
    prefix.push(b'\t');
    for (k, t) in read.iter().enumerate() {
        buf.extend_from_slice(&prefix);
        push_int(buf, k as i64);
        buf.push(b'\t');
        push_int(buf, t.value as i64);
        buf.push(b'\t');
        push_int(buf, t.s());
        buf.push(b'\t');
        push_int(buf, t.e());
        buf.extend_from_slice(b"\tq_symbol\n");
    }
}

/// An integer in decimal, as `{}` formats it.
fn push_int(buf: &mut Vec<u8>, v: i64) {
    if v < 0 {
        buf.push(b'-');
    }
    let mut x = v.unsigned_abs();
    let mut tmp = [0u8; 20];
    let mut i = tmp.len();
    loop {
        i -= 1;
        tmp[i] = b'0' + (x % 10) as u8;
        x /= 10;
        if x == 0 {
            break;
        }
    }
    buf.extend_from_slice(&tmp[i..]);
}

/// `rounds --tokens-gzip-dir DIR`: each round's token table written straight through `pigz -c -LEVEL -n -p THREADS`
/// into DIR/roundK.tokens.tsv.gz, never as uncompressed text (a round's table can take tens of GB). The same bytes as
/// pigz on the uncompressed file.
#[derive(Clone)]
pub(super) struct TokensGzip {
    pub(super) dir: PathBuf,
    pub(super) threads: usize,
    pub(super) level: u32,
}

impl TokensGzip {
    pub(super) fn write(&self, k: usize, ids: &[String], lengths: &[i64], corrected: &[Vec<Token>]) -> Result<()> {
        through_pigz(&self.dir.join(format!("round{}.tokens.tsv.gz", k)), self.level, self.threads,
                     |out| write_token_table(out, ids, lengths, corrected))
    }
}

/// The text `write` writes, through `pigz -c -LEVEL -n -p THREADS` into PATH (as `.NAME.part` beside it, then renamed).
fn through_pigz(path: &Path, level: u32, threads: usize,
                write: impl FnOnce(&mut BufWriter<std::process::ChildStdin>) -> Result<()>) -> Result<()> {
    let name = path.file_name().ok_or("pigz: no file name")?.to_string_lossy();
    let part = path.with_file_name(format!(".{}.part", name));
    let mut child = std::process::Command::new("pigz")
        .args(["-c", &format!("-{}", level), "-n", "-p", &threads.to_string()])
        .stdin(std::process::Stdio::piped())
        .stdout(File::create(&part)?)
        .spawn()
        .map_err(|e| format!("pigz: {}", e))?;
    let written = {
        // pigz's stdin closes when `out` is dropped at the end of this block
        let mut out = BufWriter::with_capacity(1 << 20, child.stdin.take().ok_or("pigz: no stdin")?);
        write(&mut out)
    };
    let status = child.wait()?;
    written?;
    if !status.success() {
        return Err(format!("pigz for {}: {}", name, status).into());
    }
    std::fs::rename(&part, path)?;
    Ok(())
}

/// The token table of `rows` (a read's index in ids / lengths, its kept tokens) alone, renumbered 0.. in order:
/// the text sparse's TableWriter writes for them.
fn write_subset_table<W: Write>(out: &mut W, ids: &[String], lengths: &[i64], rows: &[(usize, Vec<Token>)]) -> Result<()> {
    write_table_blocks(out, rows.len(), |buf, i| format_rows(buf, &ids[rows[i].0], i, lengths[rows[i].0], &rows[i].1))
}

/// `rounds --sparse-last INV:FRACTION`: the stopping round's sparse assembly table made from its corrected reads in
/// memory, as `sparse sparse-repeat` makes it from the round's table (the same rule code: `keep_token`,
/// `repetitive_tokens`; the same text): 1 in INV tokens by value hash, then metaMDBG's repetitive values (FRACTION)
/// dropped from those; reads left without tokens dropped and the rest renumbered. Into DIR/roundK.repeat.tsv.gz,
/// through `pigz -c -9 -n -p THREADS`.
#[allow(clippy::too_many_arguments)]
pub(super) fn write_sparse_last(dir: &Path, k: usize, ids: &[String], lengths: &[i64], corrected: &[Vec<Token>], inverse_density: u64,
                     fraction: f64, threads: usize) -> Result<()> {
    use crate::sparse::{keep_token, repetitive_tokens};
    let pool = table_pool();
    // sparse's errors on such tables: no rows, or a read without rows before one with rows (its ordinals skip)
    if corrected.iter().all(Vec::is_empty) {
        return Err("empty/invalid generic token bundle".into());
    }
    if let Some(e) = corrected.iter().position(Vec::is_empty) {
        if corrected[e..].iter().any(|r| !r.is_empty()) {
            return Err("generic token rows are not ordinal-major (or a read reappears)".into());
        }
    }
    let sparse: Vec<(usize, Vec<Token>)> = pool.install(|| corrected.par_iter().enumerate().filter_map(|(r, read)| {
        let kept: Vec<Token> = read.iter().copied().filter(|t| keep_token(t.value, inverse_density)).collect();
        (!kept.is_empty()).then_some((r, kept))
    }).collect());
    if let Some(&(r, _)) = sparse.iter().find(|(r, _)| ids[*r].contains(|c: char| c == '\t' || c == '"' || c == '\n' || c == '\r')) {
        return Err(format!("read id {:?} would need quoting", ids[r]).into());   // sparse refuses it too
    }
    let counts: HashMap<u32, u64> = pool.install(|| sparse.par_iter().fold(HashMap::new, |mut m: HashMap<u32, u64>, (_, read)| {
        for t in read {
            *m.entry(t.value).or_insert(0) += 1;
        }
        m
    }).reduce(HashMap::new, |mut a, b| {
        for (v, c) in b {
            *a.entry(v).or_insert(0) += c;
        }
        a
    }));
    let repetitive = repetitive_tokens(&counts, fraction);
    let drop: HashSet<u32> = repetitive.iter().map(|&(v, _)| v).collect();
    let filtered: Vec<(usize, Vec<Token>)> = pool.install(|| sparse.par_iter().filter_map(|(r, read)| {
        let kept: Vec<Token> = read.iter().copied().filter(|t| !drop.contains(&t.value)).collect();
        (!kept.is_empty()).then_some((*r, kept))
    }).collect());
    std::mem::drop(sparse);
    through_pigz(&dir.join(format!("round{}.repeat.tsv.gz", k)), 9, threads.max(1),
                 |out| write_subset_table(out, ids, lengths, &filtered))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The token table formatted in parallel blocks equals the one-row-at-a-time `writeln!` formatting (empty reads,
    /// read counts around the block and batch sizes, extreme values).
    #[test]
    fn parallel_token_table_equals_writeln() {
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut next = |m: u64| -> u64 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % m
        };
        for n in [0usize, 1, 63, 64, 65, 2047, 2048, 2049, 5000] {
            let ids: Vec<String> = (0..n).map(|r| format!("read-{:x}-{}", next(1 << 40), r)).collect();
            let lengths: Vec<i64> = (0..n).map(|_| next(1 << 31) as i64).collect();
            let corrected: Vec<Vec<Token>> = (0..n).map(|_| {
                let len = if next(10) == 0 { 0 } else { next(30) as usize };
                (0..len).map(|_| {
                    let value = match next(5) { 0 => 0, 1 => u32::MAX, _ => next(1 << 32) as u32 };
                    let start = next(1 << 21) as i64;
                    Token::new(value, start, start + next(1 << 11) as i64)
                }).collect()
            }).collect();
            let mut got = Vec::new();
            write_token_table(&mut got, &ids, &lengths, &corrected).unwrap();
            let mut want = Vec::new();
            writeln!(want, "read_id\tread_ordinal\tsource_length\ttoken_ordinal\ttoken_id\tstart\tend\tcoordinate_space").unwrap();
            for (r, read) in corrected.iter().enumerate() {
                for (k, t) in read.iter().enumerate() {
                    writeln!(want, "{}\t{}\t{}\t{}\t{}\t{}\t{}\tq_symbol", ids[r], r, lengths[r], k, t.value, t.s(), t.e()).unwrap();
                }
            }
            assert_eq!(got, want, "{} reads", n);
        }
        for v in [0i64, 9, 10, 99, 100, -1, -10, i64::MAX, i64::MIN + 1] {
            let mut b = Vec::new();
            push_int(&mut b, v);
            assert_eq!(String::from_utf8(b).unwrap(), v.to_string());
        }
    }

    /// `write_sparse_last` writes the text that `sparse sparse-repeat` writes from the round's table (empty reads,
    /// reads emptied by either step, few values so that some are repetitive; several densities and fractions).
    #[test]
    fn sparse_last_equals_the_sparse_step() {
        if std::process::Command::new("pigz").arg("--version").output().is_err() {
            eprintln!("pigz not installed: skipped");
            return;
        }
        let mut state: u64 = 0x1234_5678_9ABC_DEF1;
        let mut next = |m: u64| -> u64 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % m
        };
        let dir = std::env::temp_dir().join(format!("sparse_last_{}", std::process::id()));
        // empty: 0 = no read without tokens, 1 = the last reads without tokens, 2 = one in the middle (both refuse)
        for (case, &(n, values, inv, fraction, empty)) in [(0usize, 50u64, 3u64, 0.1f64, 0), (1, 50, 1, 0.5, 0),
                                                            (300, 40, 3, 0.05, 1), (3000, 5000, 8, 0.0001, 0),
                                                            (500, 3, 2, 0.4, 1), (50, 40, 3, 0.05, 2)].iter().enumerate() {
            let d = dir.join(case.to_string());
            std::fs::create_dir_all(&d).unwrap();
            let ids: Vec<String> = (0..n).map(|r| format!("read-{:x}-{}", next(1 << 40), r)).collect();
            let lengths: Vec<i64> = (0..n).map(|_| next(1 << 31) as i64).collect();
            let mut corrected: Vec<Vec<Token>> = (0..n).map(|_| {
                (0..1 + next(40) as usize).map(|_| {
                    let start = next(1 << 20) as i64;
                    Token::new(next(values) as u32, start, start + next(1 << 10) as i64)
                }).collect()
            }).collect();
            match empty {
                1 => corrected.iter_mut().skip(n - 5).for_each(Vec::clear),
                2 => corrected[n / 2].clear(),
                _ => {}
            }
            let mut table = BufWriter::new(File::create(d.join("in.tsv")).unwrap());
            write_token_table(&mut table, &ids, &lengths, &corrected).unwrap();
            drop(table);
            let p = |name: &str| d.join(name);
            let want = crate::sparse::sparse_repeat(&p("in.tsv"), inv, fraction, 0, &p("want.repeat.tsv"));
            let got = write_sparse_last(&d, 7, &ids, &lengths, &corrected, inv, fraction, 3);
            if n == 0 || empty == 2 {
                assert_eq!(want.unwrap_err().to_string(), got.unwrap_err().to_string(), "case {}", case);
                continue;
            }
            want.unwrap();
            got.unwrap();
            let unzip = |name: &str| -> Vec<u8> {
                let mut text = Vec::new();
                std::io::Read::read_to_end(&mut MultiGzDecoder::new(File::open(p(name)).unwrap()), &mut text).unwrap();
                text
            };
            let read = |name: &str| std::fs::read(p(name)).unwrap();
            assert_eq!(unzip("round7.repeat.tsv.gz"), read("want.repeat.tsv"), "case {}", case);
            // the same .gz bytes as pigz -9 -n -p 1 on the text
            let gz = std::process::Command::new("pigz").args(["-c", "-9", "-n", "-p", "1"]).arg(p("want.repeat.tsv")).output().unwrap();
            assert_eq!(read("round7.repeat.tsv.gz"), gz.stdout, "case {}", case);
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
