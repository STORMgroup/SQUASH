//! The signal step: POD5 signal -> RawHash2 R10.4.1 events -> 3-bit q-symbols, and the events' POD5 sample spans for
//! the signal layout (the detector and quantizer are in `events`).
//!
//! One pass over the POD5 writes every read's uncorrected q-symbols (`--signal-q-symbols`) and its events' POD5 sample
//! spans (`--dump-event-spans`). Reads are met in a fixed order (`pod5_fast5_order`, blocks of `READ_BLOCK` reads) and
//! written in that order, whatever the thread count.

use crate::events;
use crate::events::RawParams;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use structopt::StructOpt;

// ---- events lines and the POD5 signal

/// One events-file line: (read id, event means), fields split on spaces, tabs or commas; None for an empty or `#`
/// line.
pub fn parse_events_line(line: &str) -> Option<(String, Vec<f32>)> {
    let trimmed = line.trim_end_matches(['\n', '\r']);
    if trimmed.is_empty() || trimmed.starts_with('#') {
        return None;
    }
    let mut tokens = trimmed.split([' ', '\t', ',']).filter(|s| !s.is_empty());
    // A line that survives the empty/`#` check always has at least one token.
    let id = tokens
        .next()
        .expect("events-file: no read_id on non-empty line")
        .to_string();
    let events: Vec<f32> = tokens
        .map(|s| {
            s.parse::<f32>().unwrap_or_else(|_| {
                panic!("events-file: bad f32 '{}' (read {})", s, id);
            })
        })
        .collect();
    Some((id, events))
}


/// The order in which the reads of one POD5 file are processed, given the file's
/// read ids in POD5 order: the order of a FAST5 conversion read back file by
/// file (`pod5 convert to_fast5 --file-read-count chunk`, reads sorted by id
/// within a file). 0 = one file (all reads sorted by id).
pub fn pod5_fast5_order(read_ids: &[String], chunk: usize) -> Vec<usize> {
    let n = read_ids.len();
    let chunk = if chunk == 0 { n.max(1) } else { chunk };
    // File names `<stem>.<c>_0.fast5` sort as strings: chunk 10 before chunk 2.
    let mut chunks: Vec<(String, usize)> = (0..(n + chunk - 1) / chunk)
        .map(|c| (format!("{}_0.fast5", c), c))
        .collect();
    chunks.sort();
    let mut order = Vec::with_capacity(n);
    for (_, c) in chunks {
        let mut idx: Vec<usize> = (c * chunk..((c + 1) * chunk).min(n)).collect();
        idx.sort_by(|&a, &b| read_ids[a].cmp(&read_ids[b]));
        order.extend(idx);
    }
    order
}


/// The read's signal in pA (samples outside 30-200 pA dropped when `pa_filter`), with, for every kept sample, its
/// index in the POD5 signal, and the POD5 signal's length (for `--dump-event-spans`).
pub fn pod5_signal_pa_indexed(
    read: &pod5lib::ReadRecord,
    pa_filter: bool,
    pa_min: f32,
    pa_max: f32,
) -> Result<(Vec<f32>, Vec<u32>, usize), String> {
    let adc = read
        .signal()
        .map_err(|e| format!("POD5: cannot read the signal of {}: {}", read.read_id_str(), e))?;
    let (off, scale) = (read.calibration.offset, read.calibration.scale);
    let mut out: Vec<f32> = Vec::with_capacity(adc.len());
    let mut idx: Vec<u32> = Vec::with_capacity(adc.len());
    for (i, &v) in adc.iter().enumerate() {
        let pa = (v as f32 + off) * scale;
        if !pa_filter || (pa > pa_min && pa < pa_max) {
            out.push(pa);
            idx.push(i as u32);
        }
    }
    Ok((out, idx, adc.len()))
}

// ---- q-symbols files

/// The quantizer line of a q-symbols file's header (the float parameters as their f32 bits).
fn quantizer_line(params: &RawParams) -> String {
    format!(
        "#quant_bit={} n_buckets={} fine_min_bits={:08x} fine_max_bits={:08x} fine_range_bits={:08x} positions=implicit_zero_based",
        params.quant_bit,
        params.n_buckets,
        params.fine_min.to_bits(),
        params.fine_max.to_bits(),
        params.fine_range.to_bits(),
    )
}

/// Appends `symbols` as comma-separated decimals.
fn push_symbols(out: &mut Vec<u8>, symbols: &[u8]) {
    for (i, &v) in symbols.iter().enumerate() {
        if i > 0 {
            out.push(b',');
        }
        if v >= 100 {
            out.push(b'0' + v / 100);
        }
        if v >= 10 {
            out.push(b'0' + (v / 10) % 10);
        }
        out.push(b'0' + v % 10);
    }
}

/// One events-file line (`read_id<TAB>event...`) as its q-symbols line (`read_id<TAB>s,s,...\n`), the read id and the
/// symbol count; None for a comment or empty line. The fused POD5 route (`--signal-q-symbols`) formats each read's
/// events line in memory and passes it here, so its q-symbols come from the same text as the events file's.
pub fn q_symbols_line(line: &str, params: &RawParams) -> Option<(String, Vec<u8>, u64)> {
    parse_events_line(line).map(|(id, events)| {
        let q = events::quantize_events(&events, params);
        let mut out = Vec::with_capacity(id.len() + 2 * q.values.len() + 2);
        out.extend_from_slice(id.as_bytes());
        out.push(b'\t');
        push_symbols(&mut out, &q.values);
        out.push(b'\n');
        (id, out, q.values.len() as u64)
    })
}

/// The header lines of a q-symbols file for `params`: the format tag, the quantizer, the columns.
pub fn q_symbols_header(params: &RawParams) -> String {
    format!("{}\n{}\n#columns=read_id\tsymbols\n", crate::io::QSYMBOLS_TAG, quantizer_line(params))
}






// ---- the two POD5 routes

const DETECTED_EVENTS_BATCH: usize = 256;

/// One events-file line: `read_id`, then each event with six decimals.
fn detected_events_line(read_id: &str, events: &[f32]) -> Vec<u8> {
    let mut line = Vec::with_capacity(read_id.len() + 10 * events.len() + 1);
    write!(line, "{}", read_id).expect("format detected events read_id");
    for ev in events {
        write!(line, "\t{:.6}", ev).expect("format detected event");
    }
    line.push(b'\n');
    line
}


/// The signal step: POD5 -> events -> q-symbols. Reads go in `pod5_fast5_order`; per batch, the pool converts each
/// read's signal to pA, detects and merges its events, formats its events line (`detected_events_line`) and turns the
/// line into its q-symbols line (`q_symbols_line`); reads without events are skipped, and with `--keep-reads` so are
/// reads not listed. Written in read order: the q-symbols file (`q_symbols_header`, then one line per read), the
/// and each written read's event spans (`--dump-event-spans`;
/// little-endian: name length u32, name, POD5 signal length u32, events u32, then every event's first sample u32,
/// then every event's last sample + 1 u32).
#[allow(clippy::too_many_arguments)]
fn dump_signal_q_symbols_pod5(
    pod5_path: &PathBuf,
    qsym_path: &PathBuf,
    keep_path: Option<&PathBuf>,
    rp: &events::RawParams,
    event_merge: f32,
    sort_chunk: usize,
    threads: usize,
    spans_path: &PathBuf,
) {
    let files: Vec<PathBuf> = if pod5_path.is_dir() {
        let mut v: Vec<PathBuf> = std::fs::read_dir(pod5_path)
            .unwrap_or_else(|e| panic!("Cannot read POD5 directory {}: {}", pod5_path.display(), e))
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.is_file() && p.extension().and_then(|s| s.to_str()) == Some("pod5"))
            .collect();
        v.sort();
        if v.is_empty() {
            panic!("no *.pod5 files in directory {}", pod5_path.display());
        }
        v
    } else {
        vec![pod5_path.clone()]
    };
    let keep: Option<std::collections::HashSet<String>> = keep_path.map(|p| {
        std::fs::read_to_string(p)
            .unwrap_or_else(|e| panic!("Cannot read --keep-reads {}: {}", p.display(), e))
            .lines()
            .map(|l| l.trim().to_string())
            .filter(|l| !l.is_empty())
            .collect()
    });
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .expect("build the POD5 q-symbol pool");
    let create = |p: &PathBuf, what: &str| {
        BufWriter::with_capacity(1 << 24, File::create(p).unwrap_or_else(|e| panic!("Cannot create {} {}: {}", what, p.display(), e)))
    };
    let mut qsym = create(qsym_path, "signal-q-symbols output");
    qsym.write_all(q_symbols_header(rp).as_bytes()).expect("write q-symbols header");
    let mut spans_out = create(spans_path, "dump-event-spans output");
    println!(
        "[signal-q-symbols] detector win1={} win2={} thr1={} thr2={} peak={} min_seg={} max_seg={} pa_filter=true merge_short=0 event_merge={}; {}",
        rp.window_length1,
        rp.window_length2,
        rp.threshold1,
        rp.threshold2,
        rp.peak_height,
        rp.min_seg_len,
        rp.max_seg_len,
        event_merge,
        q_symbols_header(rp).lines().nth(1).unwrap_or("")
    );
    let batch = (threads.max(1) * 16).max(DETECTED_EVENTS_BATCH);
    let mut written = 0u64;
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    let (mut seen_reads, mut skipped, mut empty_events, mut total_events, mut total_symbols, mut event_bytes) =
        (0u64, 0u64, 0u64, 0u64, 0u64, 0u64);
    for fpath in &files {
        let reader = pod5lib::Reader::open(fpath)
            .unwrap_or_else(|e| panic!("Cannot open POD5 input {}: {}", fpath.display(), e));
        let recs = reader.reads();
        let ids: Vec<String> = recs.iter().map(|r| r.read_id_str()).collect();
        let order = pod5_fast5_order(&ids, sort_chunk);
        println!(
            "[signal-q-symbols] {} read(s) from {} (pod5, sort chunk {})",
            recs.len(),
            fpath.display(),
            sort_chunk
        );
        for part in order.chunks(batch) {
            // (kept, events, events line, q-symbols line, symbols, spans record)
            #[allow(clippy::type_complexity)]
            let lines: Vec<(bool, usize, Option<Vec<u8>>, Option<Vec<u8>>, u64, Vec<u8>)> = pool.install(|| {
                part.par_iter()
                    .map(|&i| {
                        if let Some(k) = &keep {
                            if !k.contains(&ids[i]) {
                                return (false, 0, None, None, 0, Vec::new());
                            }
                        }
                        // The pA filter keeps samples strictly inside 30-200 pA; `kept` maps them back to POD5 samples.
                        let (signal, kept, raw_len) = pod5_signal_pa_indexed(&recs[i], true, 30.0, 200.0)
                            .unwrap_or_else(|e| panic!("{} ({})", e, fpath.display()));
                        let (events, spans) = events::detect_events_merged_spans(&signal, rp, event_merge);
                        if events.is_empty() {
                            return (true, 0, None, None, 0, Vec::new());
                        }
                        let name = ids[i].as_bytes();
                        let mut rec = Vec::with_capacity(12 + name.len() + 8 * spans.len());
                        rec.extend_from_slice(&(name.len() as u32).to_le_bytes());
                        rec.extend_from_slice(name);
                        rec.extend_from_slice(&(raw_len as u32).to_le_bytes());
                        rec.extend_from_slice(&(spans.len() as u32).to_le_bytes());
                        for &(s, _) in &spans {
                            rec.extend_from_slice(&kept[s as usize].to_le_bytes());
                        }
                        for &(_, e) in &spans {
                            rec.extend_from_slice(&(kept[e as usize - 1] + 1).to_le_bytes());
                        }
                        let line = detected_events_line(&ids[i], &events);
                        let text = std::str::from_utf8(&line).expect("events line is UTF-8");
                        let (_, q, n) = q_symbols_line(text, rp).expect("an events line with events");
                        (true, events.len(), Some(line), Some(q), n, rec)
                    })
                    .collect()
            });
            for (j, (kept, n_events, line, q, n_symbols, rec)) in lines.into_iter().enumerate() {
                seen_reads += 1;
                if !kept {
                    skipped += 1;
                    continue;
                }
                let (line, q) = match (line, q) {
                    (Some(line), Some(q)) => (line, q),
                    _ => {
                        empty_events += 1;
                        continue;
                    }
                };
                let id = &ids[part[j]];
                if !seen.insert(id.clone()) {
                    panic!("--signal-q-symbols: duplicate read_id '{}'", id);
                }
                event_bytes += line.len() as u64;
                qsym.write_all(&q).expect("write q-symbols");
                spans_out.write_all(&rec).expect("write event spans");
                written += 1;
                total_events += n_events as u64;
                total_symbols += n_symbols;
            }
        }
    }
    qsym.flush().expect("flush q-symbols");
    spans_out.flush().expect("flush event spans");
    if let Some(k) = &keep {
        let missing = k.iter().filter(|id| !seen.contains(*id)).count();
        println!("[signal-q-symbols] --keep-reads: {} listed, {} written, {} listed but not written (not in the POD5 or no events)",
                 k.len(), seen.len(), missing);
    }
    println!(
        "[signal-q-symbols] wrote {} read(s), {} symbol(s) to {} (pod5_reads={} not_kept={} empty_events={} events={} events_text_bytes={})",
        written,
        total_symbols,
        qsym_path.display(),
        seen_reads,
        skipped,
        empty_events,
        total_events,
        event_bytes,
    );
}

/// The signal step's options (`squash __signal FILE ...`): the event detector, the quantizer's bits, the outputs.
#[derive(Debug, StructOpt)]
#[structopt(name = "signal")]
struct SignalOpt {
    /// POD5 file, or a directory of *.pod5 files
    #[structopt(parse(from_os_str))]
    reads: PathBuf,
    /// Event detector: short t-test window (samples)
    #[structopt(long = "evdet-win1")]
    evdet_win1: u32,
    /// Event detector: long t-test window (samples)
    #[structopt(long = "evdet-win2")]
    evdet_win2: u32,
    /// Event detector: t-statistic threshold of the short window
    #[structopt(long = "evdet-thr1")]
    evdet_thr1: f32,
    /// Event detector: t-statistic threshold of the long window
    #[structopt(long = "evdet-thr2")]
    evdet_thr2: f32,
    /// Event detector: peak height
    #[structopt(long = "evdet-peak")]
    evdet_peak: f32,
    /// Merge consecutive events within this of their run's mean (z units; 0 = off)
    #[structopt(long = "evdet-event-merge", default_value = "0")]
    evdet_event_merge: f32,
    /// Bits per q-symbol (2^bits quantization levels)
    #[structopt(long = "quant-bit")]
    quant_bit: u32,
    /// Threads (the outputs do not depend on them)
    #[structopt(long = "threads", default_value = "8")]
    threads: usize,
    /// Write the q-symbols file here
    #[structopt(long = "signal-q-symbols", parse(from_os_str))]
    signal_q_symbols: PathBuf,
    /// Only the reads listed in this file (one id per line), after the read order is fixed over every read
    #[structopt(long = "keep-reads", parse(from_os_str))]
    keep_reads: Option<PathBuf>,
    /// Write each event's span of POD5 samples here (binary)
    #[structopt(long = "dump-event-spans", parse(from_os_str))]
    dump_event_spans: PathBuf,
}

/// Reads of a POD5 are taken in blocks of this many in file order (`pod5_fast5_order`).
pub const READ_BLOCK: usize = 40_000;

/// `RawParams` from the options (RawHash2's defaults for the rest).
fn raw_params(opt: &SignalOpt) -> RawParams {
    let mut rp = events::RawParams::default();
    if !(1..=7).contains(&opt.quant_bit) {
        panic!("--quant-bit must be in [1, 7]; got {}.", opt.quant_bit);
    }
    rp.quant_bit = opt.quant_bit;
    rp.n_buckets = 1u32 << opt.quant_bit;
    rp.window_length1 = opt.evdet_win1;
    rp.window_length2 = opt.evdet_win2;
    rp.threshold1 = opt.evdet_thr1;
    rp.threshold2 = opt.evdet_thr2;
    rp.peak_height = opt.evdet_peak;
    rp
}

/// `squash __signal FILE ...`: the q-symbols and the event spans of a POD5 input.
pub fn run() -> Result<(), String> {
    let opt = SignalOpt::from_iter(crate::argv::args());
    let rp = raw_params(&opt);
    dump_signal_q_symbols_pod5(
        &opt.reads,
        &opt.signal_q_symbols,
        opt.keep_reads.as_ref(),
        &rp,
        opt.evdet_event_merge,
        READ_BLOCK,
        opt.threads.max(1),
        &opt.dump_event_spans,
    );
    Ok(())
}
