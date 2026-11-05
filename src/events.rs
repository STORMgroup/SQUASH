//! Signal events and q-symbols: RawHash2's event detector (t-statistic peaks of two windows over the z-normalized
//! signal; each event is the mean of its segment), its three-zone quantizer (event means -> q-symbols) and its e-mer
//! sketcher (`e` consecutive q-symbols packed and hashed, all of them or windowed minimizers).
//!
//! Ported from RawHash2 (github.com/STORMgroup/RawHash2, `revent.c` and `rsketch.c` at commit 92366af); function
//! names follow their C counterparts. The sketch reproduces RawHash2's hash stream bit for bit (the golden-vector
//! test `rawhash2_sketch_parity`), so the floating-point widths follow the C code: f64 for the running sums
//! (mean_sum / std_dev_sum, sqrt), f32 everywhere else.

use std::cmp::Ordering;

#[derive(Clone, Debug)]
pub struct RawParams {
    // Event detector (revent.c)
    pub window_length1: u32,
    pub window_length2: u32,
    pub threshold1: f32,
    pub threshold2: f32,
    pub peak_height: f32,
    pub min_seg_len: u32,
    pub max_seg_len: u32,
    // Sketcher (rsketch.c)
    pub w: u32,         // window for minimizer selection; 0 = dense (ri_sketch_reg)
    pub e: u32,         // number of events packed per hash
    pub quant_bit: u32, // bits per quantized event; cap on the bucket value packed into the hash
    // Number of distinct quantization buckets dynamic_quantize returns. Independent
    // of `quant_bit` (RawHash2's default is 7, not a power of 2); the signal step's
    // options derive one from the other.
    pub n_buckets: u32,
    pub fine_min: f32,
    pub fine_max: f32,
    pub fine_range: f32,
    // Integer threshold for bucket-level event merging (RawHash2's `diff`): skip an
    // event if its bucket differs from the previous kept event's by at most
    // `bucket_diff`. 0 skips runs of consecutive same-bucket events; -1 disables.
    pub bucket_diff: i32,
    // Optional float threshold on |s[f_pos] - s[last_kept]| in normalized event
    // space (an earlier RawHash2 filter): an event closer than this is skipped.
    // Applied after `bucket_diff`; when both are enabled, an event must survive
    // both filters to contribute to the hash. Negative disables it. Read only by the RawHash2 reference sketch the
    // tests compare against (`filter_event_symbols`).
    #[cfg_attr(not(test), allow(dead_code))]
    pub sig_diff: f32,
}

impl Default for RawParams {
    // RawHash2's defaults (roptions.c `ri_idxopt_init`), with its R10.4.1 detector settings.
    fn default() -> Self {
        RawParams {
            window_length1: 3,
            window_length2: 9,
            threshold1: 4.0,
            threshold2: 3.5,
            peak_height: 0.4,
            min_seg_len: 0,
            max_seg_len: 500,
            w: 2,
            e: 8,
            quant_bit: 3,
            n_buckets: 7,
            fine_min: -1.6,
            fine_max: 1.6,
            fine_range: 0.75,
            bucket_diff: 0,
            sig_diff: -1.0,
        }
    }
}

// -----------------------------------------------------------------------------
// hash64 — RawHash2's integer hash of a packed e-mer (rsketch.c). The mask is
// (1<<32)-1 at every call site here.
// -----------------------------------------------------------------------------
#[inline]
pub fn hash64(key: u64, mask: u64) -> u64 {
    let mut k = key & mask;
    k = ((!k).wrapping_add(k.wrapping_shl(21))) & mask;
    k ^= k >> 24;
    k = (k
        .wrapping_add(k.wrapping_shl(3))
        .wrapping_add(k.wrapping_shl(8)))
        & mask;
    k ^= k >> 14;
    k = (k
        .wrapping_add(k.wrapping_shl(2))
        .wrapping_add(k.wrapping_shl(4)))
        & mask;
    k ^= k >> 28;
    k = k.wrapping_add(k.wrapping_shl(31)) & mask;
    k
}

// -----------------------------------------------------------------------------
// dynamic_quantize — RawHash2's quantizer (rsketch.c): a normalized event mean -> q-symbol.
//
// Three-zone adaptive quantizer: signal in [fine_min, fine_max] gets `n_fine`
// buckets, signal outside that band gets the remaining `n_buckets - n_fine`
// buckets split evenly between a low zone ([-3, fine_min)) and a high zone
// ((fine_max, 3]). Each zone is mapped linearly onto its bucket sub-range.
// -----------------------------------------------------------------------------
#[inline]
pub fn dynamic_quantize(
    signal: f32,
    fine_min: f32,
    fine_max: f32,
    fine_range: f32,
    n_buckets: u32,
) -> u32 {
    if n_buckets == 0 {
        return 0;
    }

    let nb_f = n_buckets as f32 * fine_range + 0.5;
    let n_fine_i = (nb_f as i32).clamp(1, n_buckets as i32);
    let n_fine = n_fine_i as u32;
    let n_low = (n_buckets - n_fine) / 2;
    let n_high = n_buckets - n_fine - n_low;

    let lo_in: f32 = -3.0;
    let hi_in: f32 = 3.0;

    if signal < fine_min {
        if n_low == 0 {
            return 0;
        }
        let denom = fine_min - lo_in;
        let t = if denom > 0.0 {
            ((signal - lo_in) / denom).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let b = (t * n_low as f32) as u32;
        if b >= n_low {
            n_low - 1
        } else {
            b
        }
    } else if signal <= fine_max {
        let denom = fine_max - fine_min;
        let t = if denom > 0.0 {
            ((signal - fine_min) / denom).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let b = (t * n_fine as f32) as u32;
        n_low + if b >= n_fine { n_fine - 1 } else { b }
    } else {
        if n_high == 0 {
            return n_buckets - 1;
        }
        let denom = hi_in - fine_max;
        let t = if denom > 0.0 {
            ((signal - fine_max) / denom).clamp(0.0, 1.0)
        } else {
            1.0
        };
        let b = (t * n_high as f32) as u32;
        n_low + n_fine + if b >= n_high { n_high - 1 } else { b }
    }
}

// -----------------------------------------------------------------------------
// Event detector — revent.c
// -----------------------------------------------------------------------------

/// Running accumulator shared across reads when streaming many signals through
/// a single normalizer. For a one-shot read, start at zero.
#[derive(Clone, Debug, Default)]
pub struct NormalizeState {
    pub mean_sum: f64,
    pub std_dev_sum: f64,
    // u64, not u32: one read never approaches 2^32 samples, but a state shared
    // across every read of a data set (millions of reads × ~10^4–10^5 samples)
    // does. No per-read effect: the value is only read as f64.
    pub n_events_sum: u64,
}

#[cfg(test)]
/// RawHash2's `normalize_signal` (revent.c). Z-normalizes `sig` using the
/// cumulative running statistics in `state`, drops samples outside ±3σ, returns
/// the surviving normalized values.
pub fn normalize_signal(sig: &[f32], state: &mut NormalizeState) -> Vec<f32> {
    let mut sum = state.mean_sum;
    let mut sum2 = state.std_dev_sum;
    for &v in sig {
        let vd = v as f64;
        sum += vd;
        sum2 += vd * vd;
    }
    state.n_events_sum += sig.len() as u64;
    state.mean_sum = sum;
    state.std_dev_sum = sum2;

    let n = state.n_events_sum as f64;
    let mean = sum / n;
    let std_dev = (sum2 / n - mean * mean).sqrt();

    let mut out = Vec::with_capacity(sig.len());
    for &v in sig {
        let norm = ((v as f64 - mean) / std_dev) as f32;
        if norm < 3.0 && norm > -3.0 {
            out.push(norm);
        }
    }
    out
}






fn comp_prefix_prefixsq(sig: &[f32]) -> (Vec<f32>, Vec<f32>) {
    let n = sig.len();
    let mut ps = vec![0.0f32; n + 1];
    let mut pss = vec![0.0f32; n + 1];
    for i in 0..n {
        ps[i + 1] = ps[i] + sig[i];
        pss[i + 1] = pss[i] + sig[i] * sig[i];
    }
    (ps, pss)
}

fn comp_tstat(prefix_sum: &[f32], prefix_sum_square: &[f32], s_len: u32, w_len: u32) -> Vec<f32> {
    let eta = f32::MIN_POSITIVE;
    let mut tstat = vec![0.0f32; (s_len + 1) as usize];
    if s_len < 2 * w_len || w_len < 2 {
        return tstat;
    }
    let w = w_len as usize;
    let s = s_len as usize;
    for i in w..=s - w {
        let mut sum1 = prefix_sum[i];
        let mut sumsq1 = prefix_sum_square[i];
        if i > w {
            sum1 -= prefix_sum[i - w];
            sumsq1 -= prefix_sum_square[i - w];
        }
        let sum2 = prefix_sum[i + w] - prefix_sum[i];
        let sumsq2 = prefix_sum_square[i + w] - prefix_sum_square[i];
        let wf = w_len as f32;
        let mean1 = sum1 / wf;
        let mean2 = sum2 / wf;
        let mut combined_var = (sumsq1 / wf - mean1 * mean1 + sumsq2 / wf - mean2 * mean2) / wf;
        if combined_var < eta {
            combined_var = eta;
        }
        let delta_mean = mean2 - mean1;
        // C uses sqrt (double) and fabs (promotes to double) — match that width.
        tstat[i] = ((delta_mean as f64).abs() / (combined_var as f64).sqrt()) as f32;
    }
    // fudge trailing boundary — C: memset(tstat + s_len - w_len + 1, 0, w_len * sizeof(float))
    let start = s - w + 1;
    for i in start..start + w {
        if i < tstat.len() {
            tstat[i] = 0.0;
        }
    }
    tstat
}

struct Detector<'a> {
    def_peak_pos: i32,
    def_peak_val: f32,
    sig: &'a [f32],
    threshold: f32,
    window_length: u32,
    masked_to: u32,
    peak_pos: i32,
    peak_value: f32,
    valid_peak: bool,
}

impl<'a> Detector<'a> {
    fn new(sig: &'a [f32], threshold: f32, window_length: u32) -> Self {
        Detector {
            def_peak_pos: -1,
            def_peak_val: f32::MAX,
            sig,
            threshold,
            window_length,
            masked_to: 0,
            peak_pos: -1,
            peak_value: f32::MAX,
            valid_peak: false,
        }
    }
}

fn gen_peaks(detectors: &mut [Detector<'_>], peak_height: f32) -> Vec<u32> {
    let s_len = detectors[0].sig.len();
    let w0 = detectors[0].window_length;
    let mut peaks = Vec::with_capacity(s_len);
    for i in 0..s_len {
        for k in 0..detectors.len() {
            if detectors[k].masked_to as usize >= i {
                continue;
            }
            let current_value = detectors[k].sig[i];

            if detectors[k].peak_pos == detectors[k].def_peak_pos {
                // CASE 1 — tracking the running minimum until we rise by peak_height.
                if current_value < detectors[k].peak_value {
                    detectors[k].peak_value = current_value;
                } else if current_value - detectors[k].peak_value > peak_height {
                    detectors[k].peak_value = current_value;
                    detectors[k].peak_pos = i as i32;
                }
            } else {
                // CASE 2 — refining an existing peak candidate.
                if current_value > detectors[k].peak_value {
                    detectors[k].peak_value = current_value;
                    detectors[k].peak_pos = i as i32;
                }
                // Mask out all *later* detectors if this one has crossed threshold.
                if detectors[k].peak_value > detectors[k].threshold {
                    let mask = detectors[k].peak_pos as u32 + w0;
                    for d in &mut detectors[k + 1..] {
                        d.masked_to = mask;
                        d.peak_pos = d.def_peak_pos;
                        d.peak_value = d.def_peak_val;
                        d.valid_peak = false;
                    }
                }
                if detectors[k].peak_value - current_value > peak_height
                    && detectors[k].peak_value > detectors[k].threshold
                {
                    detectors[k].valid_peak = true;
                }
                if detectors[k].valid_peak
                    && (i as i32 - detectors[k].peak_pos) as u32 > detectors[k].window_length / 2
                {
                    peaks.push(detectors[k].peak_pos as u32);
                    detectors[k].peak_pos = detectors[k].def_peak_pos;
                    detectors[k].peak_value = current_value;
                    detectors[k].valid_peak = false;
                }
            }
        }
    }
    peaks
}

/// The mean of a segment's values within [q1 - iqr, q3 + iqr] (sorts the segment).
fn calculate_mean_of_filtered_segment(segment: &mut [f32]) -> f32 {
    // qsort with compare_floats in C returns 0 for NaN comparisons; we treat
    // NaN as Equal here for determinism. Real nanopore signal normalised to
    // ±3σ should never produce NaN, so this just avoids panicking in edge tests.
    segment.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let n = segment.len();
    if n == 0 {
        return 0.0;
    }
    let q1 = segment[n / 4];
    let q3 = segment[3 * n / 4];
    let iqr = q3 - q1;
    let lower = q1 - iqr;
    let upper = q3 + iqr;
    let mut sum = 0.0f32;
    let mut count: u32 = 0;
    for &v in segment.iter() {
        if v >= lower && v <= upper {
            sum += v;
            count += 1;
        }
    }
    if count == 0 {
        0.0
    } else {
        sum / count as f32
    }
}

#[cfg(test)]
fn gen_events(sig: &mut [f32], peaks: &[u32], min_seg_len: u32, max_seg_len: u32) -> Vec<f32> {
    let s_len = sig.len() as u32;
    let mut events = Vec::with_capacity(peaks.len());
    let mut start_idx: u32 = 0;
    for &pk in peaks {
        if !(pk > 0 && pk < s_len) {
            continue;
        }
        let seg_len = pk - start_idx;
        if seg_len >= min_seg_len && seg_len < max_seg_len {
            let slice = &mut sig[start_idx as usize..pk as usize];
            events.push(calculate_mean_of_filtered_segment(slice));
        }
        start_idx = pk;
    }
    events
}

#[cfg(test)]
/// Top-level event detector. Runs normalize_signal + t-stat peak detection on
/// a raw pA signal and returns the normalized event-mean sequence that the
/// sketcher consumes.
pub fn detect_events(sig: &[f32], p: &RawParams) -> Vec<f32> {
    let mut state = NormalizeState::default();
    let mut norm = normalize_signal(sig, &mut state);
    if norm.is_empty() {
        return Vec::new();
    }
    let (ps, pss) = comp_prefix_prefixsq(&norm);
    let n_sig = norm.len() as u32;
    let tstat1 = comp_tstat(&ps, &pss, n_sig, p.window_length1);
    let tstat2 = comp_tstat(&ps, &pss, n_sig, p.window_length2);
    let short_det = Detector::new(&tstat1, p.threshold1, p.window_length1);
    let long_det = Detector::new(&tstat2, p.threshold2, p.window_length2);
    let mut detectors = [short_det, long_det];
    let peaks = gen_peaks(&mut detectors, p.peak_height);
    if peaks.is_empty() {
        return Vec::new();
    }
    gen_events(&mut norm, &peaks, p.min_seg_len, p.max_seg_len)
}

#[cfg(test)]
/// The event detector with the merge of close events (the reference the tests check the signal step's spans path
/// against).
///
/// `event_merge > 0`: consecutive events whose means stay within
/// `event_merge` (z units) of the running mean of the current run are averaged
/// into one event, before quantization; 0 gives exactly `detect_events`.
pub fn detect_events_merged(sig: &[f32], p: &RawParams, event_merge: f32) -> Vec<f32> {
    if event_merge <= 0.0 {
        return detect_events(sig, p);
    }
    let mut state = NormalizeState::default();
    let mut norm = normalize_signal(sig, &mut state);
    if norm.is_empty() {
        return Vec::new();
    }
    let (ps, pss) = comp_prefix_prefixsq(&norm);
    let n_sig = norm.len() as u32;
    let tstat1 = comp_tstat(&ps, &pss, n_sig, p.window_length1);
    let tstat2 = comp_tstat(&ps, &pss, n_sig, p.window_length2);
    let short_det = Detector::new(&tstat1, p.threshold1, p.window_length1);
    let long_det = Detector::new(&tstat2, p.threshold2, p.window_length2);
    let mut detectors = [short_det, long_det];
    let peaks = gen_peaks(&mut detectors, p.peak_height);
    if peaks.is_empty() {
        return Vec::new();
    }
    let events = gen_events(&mut norm, &peaks, p.min_seg_len, p.max_seg_len);
    merge_close_events(&events, event_merge)
}


#[cfg(test)]
/// Running-mean merge of close events: an event within `threshold` of the
/// current run's mean joins the run; each run becomes one event, its mean.
fn merge_close_events(events: &[f32], threshold: f32) -> Vec<f32> {
    if threshold <= 0.0 || events.is_empty() {
        return events.to_vec();
    }
    let mut out = Vec::with_capacity(events.len());
    let mut sum = events[0] as f64;
    let mut cnt = 1u32;
    for &v in &events[1..] {
        let mean = sum / cnt as f64;
        if (v as f64 - mean).abs() <= threshold as f64 {
            sum += v as f64;
            cnt += 1;
        } else {
            out.push((sum / cnt as f64) as f32);
            sum = v as f64;
            cnt = 1;
        }
    }
    out.push((sum / cnt as f64) as f32);
    out
}

/// `normalize_signal` that also returns the index in `sig` of every sample it keeps (same values).
fn normalize_signal_indexed(sig: &[f32], state: &mut NormalizeState) -> (Vec<f32>, Vec<u32>) {
    let mut sum = state.mean_sum;
    let mut sum2 = state.std_dev_sum;
    for &v in sig {
        let vd = v as f64;
        sum += vd;
        sum2 += vd * vd;
    }
    state.n_events_sum += sig.len() as u64;
    state.mean_sum = sum;
    state.std_dev_sum = sum2;

    let n = state.n_events_sum as f64;
    let mean = sum / n;
    let std_dev = (sum2 / n - mean * mean).sqrt();

    let mut out = Vec::with_capacity(sig.len());
    let mut idx = Vec::with_capacity(sig.len());
    for (i, &v) in sig.iter().enumerate() {
        let norm = ((v as f64 - mean) / std_dev) as f32;
        if norm < 3.0 && norm > -3.0 {
            out.push(norm);
            idx.push(i as u32);
        }
    }
    (out, idx)
}

/// `gen_events` that also returns each event's segment [start, end) in `sig`'s indices.
fn gen_events_spans(sig: &mut [f32], peaks: &[u32], min_seg_len: u32, max_seg_len: u32) -> (Vec<f32>, Vec<(u32, u32)>) {
    let s_len = sig.len() as u32;
    let mut events = Vec::with_capacity(peaks.len());
    let mut spans = Vec::with_capacity(peaks.len());
    let mut start_idx: u32 = 0;
    for &pk in peaks {
        if !(pk > 0 && pk < s_len) {
            continue;
        }
        let seg_len = pk - start_idx;
        if seg_len >= min_seg_len && seg_len < max_seg_len {
            let slice = &mut sig[start_idx as usize..pk as usize];
            events.push(calculate_mean_of_filtered_segment(slice));
            spans.push((start_idx, pk));
        }
        start_idx = pk;
    }
    (events, spans)
}


/// `merge_close_events` that also joins the merged events' spans (first start, last end).
fn merge_close_events_spans(events: &[f32], spans: &[(u32, u32)], threshold: f32) -> (Vec<f32>, Vec<(u32, u32)>) {
    if threshold <= 0.0 || events.is_empty() {
        return (events.to_vec(), spans.to_vec());
    }
    let mut out = Vec::with_capacity(events.len());
    let mut out_spans = Vec::with_capacity(events.len());
    let mut sum = events[0] as f64;
    let mut cnt = 1u32;
    let mut span = spans[0];
    for (&v, &sp) in events[1..].iter().zip(&spans[1..]) {
        let mean = sum / cnt as f64;
        if (v as f64 - mean).abs() <= threshold as f64 {
            sum += v as f64;
            cnt += 1;
            span.1 = sp.1;
        } else {
            out.push((sum / cnt as f64) as f32);
            out_spans.push(span);
            sum = v as f64;
            cnt = 1;
            span = sp;
        }
    }
    out.push((sum / cnt as f64) as f32);
    out_spans.push(span);
    (out, out_spans)
}

/// `detect_events_merged` with each event's samples (the signal step's `--dump-event-spans`): the same events, and
/// for each the span (first sample, last sample + 1) in `sig`'s indices, through the samples the normalisation drops
/// (outside +-3 sd) and the events merged by `event_merge`. Segments the detector drops (too short or too long)
/// leave gaps between spans.
pub fn detect_events_merged_spans(sig: &[f32], p: &RawParams, event_merge: f32) -> (Vec<f32>, Vec<(u32, u32)>) {
    let mut state = NormalizeState::default();
    let (mut norm, idx) = normalize_signal_indexed(sig, &mut state);
    if norm.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let (ps, pss) = comp_prefix_prefixsq(&norm);
    let n_sig = norm.len() as u32;
    let tstat1 = comp_tstat(&ps, &pss, n_sig, p.window_length1);
    let tstat2 = comp_tstat(&ps, &pss, n_sig, p.window_length2);
    let short_det = Detector::new(&tstat1, p.threshold1, p.window_length1);
    let long_det = Detector::new(&tstat2, p.threshold2, p.window_length2);
    let mut detectors = [short_det, long_det];
    let peaks = gen_peaks(&mut detectors, p.peak_height);
    if peaks.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let (events, segs) = gen_events_spans(&mut norm, &peaks, p.min_seg_len, p.max_seg_len);
    let (events, segs) = merge_close_events_spans(&events, &segs, event_merge);
    let spans = segs.iter().map(|&(s, e)| (idx[s as usize], idx[e as usize - 1] + 1)).collect();
    (events, spans)
}

// -----------------------------------------------------------------------------
// Sketcher — rsketch.c
// -----------------------------------------------------------------------------

// Both sketch variants use mask = (1<<32)-1 for hash64, so the produced hash
// word fits in the low 32 bits (returned as u64).
const HASH_MASK: u64 = (1u64 << 32) - 1;

/// Quantized event symbols and their coordinates in the source event stream.
///
/// Correction operates on `values`; `event_positions` keeps each symbol's source
/// event, so that the exact `(position, hash)` stream can be regenerated. A
/// packed e-mer is located at the first symbol it covers (RawHash2's C ring
/// buffer has this start-coordinate behavior even though some of its comments
/// call it a "lastPos").
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RawSymbols {
    pub values: Vec<u8>,
    pub event_positions: Vec<usize>,
}

/// One sketch token (an e-mer hash) of the graph with its exact source q-symbol span.
///
/// `start` and `end` are zero-based, half-open coordinates in the unfiltered
/// q-symbol stream.  In particular, `end` is the source position of the last
/// retained symbol covered by the e-mer plus one.  It therefore cannot in
/// general be reconstructed as `start + e` after bucket/signal filtering has
/// removed symbols between retained positions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RawSketchToken {
    pub start: usize,
    pub end: usize,
    pub hash: u64,
}

impl RawSymbols {
    pub fn new(values: Vec<u8>, event_positions: Vec<usize>) -> Self {
        assert_eq!(
            values.len(),
            event_positions.len(),
            "raw symbols and provenance positions must have equal lengths"
        );
        for pair in event_positions.windows(2) {
            assert!(
                pair[0] <= pair[1],
                "raw symbol provenance positions must be monotone"
            );
        }
        RawSymbols {
            values,
            event_positions,
        }
    }

    /// Construct a corrected/synthetic symbol stream whose coordinates are
    /// symbol ordinals. This deliberately performs no quantization or
    /// bucket/signal-difference filtering.
    pub fn from_values(values: Vec<u8>) -> Self {
        let event_positions = (0..values.len()).collect();
        RawSymbols::new(values, event_positions)
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
}

fn validate_symbol_params(p: &RawParams) {
    assert!(
        (1..=8).contains(&p.quant_bit),
        "quant_bit must be in [1, 8]"
    );
    assert!(
        p.n_buckets >= 1 && p.n_buckets <= (1u32 << p.quant_bit),
        "n_buckets must fit in quant_bit"
    );
}

/// Quantize every event without applying `bucket_diff` or `sig_diff`.
///
/// This is the q-symbol stream correction works on: other `e`, `w` and graph
/// `k` values can be derived from it with the same quantizer (the one its
/// file's header records). Changing the quantizer still requires the
/// normalized event means, not this discrete symbol stream.
pub fn quantize_events(s_values: &[f32], p: &RawParams) -> RawSymbols {
    validate_symbol_params(p);
    let mask_quant_bit = (1u32 << p.quant_bit) - 1;
    let values = s_values
        .iter()
        .map(|&value| {
            (dynamic_quantize(value, p.fine_min, p.fine_max, p.fine_range, p.n_buckets)
                & mask_quant_bit) as u8
        })
        .collect();
    RawSymbols::from_values(values)
}

#[cfg(test)]
/// Apply the two pre-sketch filters (`bucket_diff`, then `sig_diff`) to an
/// already-quantized event stream. The first event always survives. Both
/// comparisons are anchored at the last event that survived *both* filters; a
/// rejection advances neither anchor.
pub fn filter_event_symbols(symbols: &RawSymbols, s_values: &[f32], p: &RawParams) -> RawSymbols {
    assert_eq!(
        symbols.len(),
        s_values.len(),
        "event-value filtering requires one source value per quantized symbol"
    );
    if symbols.is_empty() {
        return RawSymbols::new(Vec::new(), Vec::new());
    }

    let mut values = Vec::with_capacity(symbols.len());
    let mut positions = Vec::with_capacity(symbols.len());
    values.push(symbols.values[0]);
    positions.push(symbols.event_positions[0]);

    let mut previous_q = symbols.values[0];
    let mut last_signal_index = 0usize;
    for i in 1..symbols.len() {
        let q = symbols.values[i];
        if p.bucket_diff >= 0 {
            let qdiff = (q as i32 - previous_q as i32).abs();
            if qdiff <= p.bucket_diff {
                continue;
            }
        }
        if p.sig_diff >= 0.0 && (s_values[i] - s_values[last_signal_index]).abs() < p.sig_diff {
            continue;
        }
        values.push(q);
        positions.push(symbols.event_positions[i]);
        previous_q = q;
        last_signal_index = i;
    }
    RawSymbols::new(values, positions)
}

/// Apply only the integer `bucket_diff` filter to a q-symbol stream.
///
/// Corrected symbols no longer have a one-to-one set of floating-point event
/// means (insertions in particular have no meaningful value to dequantize),
/// so the `sig_diff` layer cannot be replayed after correction. The assembly
/// uses this helper to apply the bucket filter after consensus.
pub fn filter_symbol_buckets(symbols: &RawSymbols, p: &RawParams) -> RawSymbols {
    if symbols.is_empty() || p.bucket_diff < 0 {
        return symbols.clone();
    }

    let mut values = Vec::with_capacity(symbols.len());
    let mut positions = Vec::with_capacity(symbols.len());
    values.push(symbols.values[0]);
    positions.push(symbols.event_positions[0]);

    let mut previous_q = symbols.values[0];
    for i in 1..symbols.len() {
        let q = symbols.values[i];
        let qdiff = (q as i32 - previous_q as i32).abs();
        if qdiff <= p.bucket_diff {
            continue;
        }
        values.push(q);
        positions.push(symbols.event_positions[i]);
        previous_q = q;
    }
    RawSymbols::new(values, positions)
}

fn visit_dense_emers(symbols: &RawSymbols, p: &RawParams, mut emit: impl FnMut(RawSketchToken)) {
    validate_symbol_params(p);
    let e = p.e as usize;
    assert!(e > 0, "signal-events/e must be positive");
    assert!(
        e * p.quant_bit as usize <= 64,
        "dense e-mer packing exceeds 64 bits"
    );
    let packed_bits = p.quant_bit * p.e;
    let mask_events = if packed_bits == 64 {
        u64::MAX
    } else {
        (1u64 << packed_bits) - 1
    };
    let mut quant_value = 0u64;
    for (i, &symbol) in symbols.values.iter().enumerate() {
        assert!(
            (symbol as u32) < p.n_buckets,
            "raw symbol {} at index {} is outside quantizer alphabet [0, {})",
            symbol,
            i,
            p.n_buckets
        );
        quant_value = ((quant_value << p.quant_bit) | symbol as u64) & mask_events;
        if i + 1 >= e {
            let end = symbols.event_positions[i]
                .checked_add(1)
                .expect("raw e-mer source endpoint overflowed usize");
            emit(RawSketchToken {
                start: symbols.event_positions[i + 1 - e],
                end,
                hash: hash64(quant_value, HASH_MASK),
            });
        }
    }
}

fn dense_output_len(symbols: &RawSymbols, e: usize) -> usize {
    if symbols.len() >= e {
        symbols.len() - e + 1
    } else {
        0
    }
}

#[cfg(test)]
/// Pack each consecutive `e`-symbol word and hash it. Output coordinates are
/// the source position of the word's first symbol.
pub fn pack_dense_emers(symbols: &RawSymbols, p: &RawParams) -> Vec<(usize, u64)> {
    let mut dense = Vec::with_capacity(dense_output_len(symbols, p.e as usize));
    visit_dense_emers(symbols, p, |token| {
        dense.push((token.start, token.hash));
    });
    dense
}

/// Pack dense e-mers while retaining their exact half-open source intervals.
pub fn pack_dense_emers_with_spans(symbols: &RawSymbols, p: &RawParams) -> Vec<RawSketchToken> {
    let mut dense = Vec::with_capacity(dense_output_len(symbols, p.e as usize));
    visit_dense_emers(symbols, p, |token| dense.push(token));
    dense
}

/// RawHash2's circular minimizer selection over a dense e-mer stream.
///
/// `<=` deliberately selects the rightmost minimum. Equal hashes are emitted
/// in RawHash2's ring order and are neither sorted nor deduplicated. The
/// selector is generic so that the graph path keeps exact e-mer endpoints while
/// `window_emers` returns RawHash2's `(start, hash)` pairs.
#[allow(clippy::needless_range_loop)]
fn window_emers_impl<T: Copy>(
    dense: &[T],
    window: usize,
    fields: impl Fn(T) -> (usize, u64),
) -> Vec<T> {
    assert!(
        window > 0 && window < 256,
        "signal window must be in (0, 256)"
    );
    if dense.is_empty() {
        return Vec::new();
    }

    #[derive(Copy, Clone)]
    struct Slot<T: Copy> {
        hash: u64,
        pos: usize,
        token: Option<T>,
    }
    let invalid = Slot {
        hash: u64::MAX,
        pos: usize::MAX,
        token: None,
    };
    let mut buf = vec![invalid; window];
    let mut min = invalid;
    let mut min_pos = 0usize;
    let mut buf_pos = 0usize;
    let mut out = Vec::new();

    for (dense_index, &token) in dense.iter().enumerate() {
        let (pos, hash) = fields(token);
        let seen = dense_index + 1;
        let info = Slot {
            hash,
            pos,
            token: Some(token),
        };
        buf[buf_pos] = info;

        if seen == window && min.token.is_some() {
            for j in (buf_pos + 1)..window {
                if min.hash == buf[j].hash && min.pos != buf[j].pos {
                    out.push(buf[j].token.unwrap());
                }
            }
            for j in 0..buf_pos {
                if min.hash == buf[j].hash && min.pos != buf[j].pos {
                    out.push(buf[j].token.unwrap());
                }
            }
        }
        if min.token.is_none() || info.hash <= min.hash {
            if seen > window && min.token.is_some() {
                out.push(min.token.unwrap());
            }
            min = info;
            min_pos = buf_pos;
        } else if buf_pos == min_pos {
            if seen >= window && min.token.is_some() {
                out.push(min.token.unwrap());
            }
            min = invalid;
            for j in (buf_pos + 1)..window {
                if buf[j].token.is_some() && (min.token.is_none() || min.hash >= buf[j].hash) {
                    min = buf[j];
                    min_pos = j;
                }
            }
            for j in 0..=buf_pos {
                if buf[j].token.is_some() && (min.token.is_none() || min.hash >= buf[j].hash) {
                    min = buf[j];
                    min_pos = j;
                }
            }
            if seen >= window && min.token.is_some() {
                for j in (buf_pos + 1)..window {
                    if buf[j].token.is_some() && min.hash == buf[j].hash && min.pos != buf[j].pos {
                        out.push(buf[j].token.unwrap());
                    }
                }
                for j in 0..=buf_pos {
                    if buf[j].token.is_some() && min.hash == buf[j].hash && min.pos != buf[j].pos {
                        out.push(buf[j].token.unwrap());
                    }
                }
            }
        }
        buf_pos += 1;
        if buf_pos == window {
            buf_pos = 0;
        }
    }
    if min.token.is_some() {
        out.push(min.token.unwrap());
    }
    out
}

#[cfg(test)]
pub fn window_emers(dense: &[(usize, u64)], window: usize) -> Vec<(usize, u64)> {
    window_emers_impl(dense, window, |token| token)
}

fn window_emers_with_spans(dense: &[RawSketchToken], window: usize) -> Vec<RawSketchToken> {
    window_emers_impl(dense, window, |token| (token.start, token.hash))
}

#[cfg(test)]
/// Sketch a symbol stream that has already undergone the desired correction-
/// space filtering. This function never quantizes or reapplies either diff
/// filter.
pub fn sketch_symbols(symbols: &RawSymbols, p: &RawParams) -> Vec<(usize, u64)> {
    let dense = pack_dense_emers(symbols, p);
    if p.w == 0 {
        dense
    } else {
        assert!(
            p.e * p.quant_bit <= 64 - 6,
            "windowed e-mer packing exceeds RawHash2's 58-bit budget"
        );
        window_emers(&dense, p.w as usize)
    }
}

/// Sketch a filtered q-symbol stream while retaining each selected e-mer's
/// exact half-open coordinates in the unfiltered source symbol stream.
pub fn sketch_symbols_with_spans(symbols: &RawSymbols, p: &RawParams) -> Vec<RawSketchToken> {
    let dense = pack_dense_emers_with_spans(symbols, p);
    if p.w == 0 {
        dense
    } else {
        assert!(
            p.e * p.quant_bit <= 64 - 6,
            "windowed e-mer packing exceeds RawHash2's 58-bit budget"
        );
        window_emers_with_spans(&dense, p.w as usize)
    }
}

#[cfg(test)]
/// RawHash2's dense sketch (`w` = 0): quantize, filter, hash every e-mer.
pub fn ri_sketch_reg(s_values: &[f32], p: &RawParams) -> Vec<(usize, u64)> {
    let symbols = quantize_events(s_values, p);
    let filtered = filter_event_symbols(&symbols, s_values, p);
    pack_dense_emers(&filtered, p)
}

#[cfg(test)]
/// RawHash2's windowed-minimizer sketch (`w` > 0).
pub fn ri_sketch_min(s_values: &[f32], p: &RawParams) -> Vec<(usize, u64)> {
    assert!(p.w > 0, "ri_sketch_min requires a positive window");
    let symbols = quantize_events(s_values, p);
    let filtered = filter_event_symbols(&symbols, s_values, p);
    sketch_symbols(&filtered, p)
}

#[cfg(test)]
/// Sketches normalized event means as RawHash2 does: (event position, hash) pairs.
pub fn sketch_signal(s_values: &[f32], p: &RawParams) -> Vec<(usize, u64)> {
    if s_values.is_empty() {
        return Vec::new();
    }
    if p.w == 0 {
        ri_sketch_reg(s_values, p)
    } else {
        ri_sketch_min(s_values, p)
    }
}

// -----------------------------------------------------------------------------
// Tests — the golden vectors (`RAW_HASH2_GOLDEN`) are RawHash2's own outputs (its C code at commit 92366af).
// -----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    #[test]
    fn merged_detector_options_are_opt_in_and_merge() {
        // A noisy staircase: 10 levels of 14 samples with a one-sample spike
        // inside some levels (a short split candidate).
        let mut sig = Vec::new();
        let levels = [80.0f32, 120.0, 95.0, 140.0, 70.0, 110.0, 90.0, 130.0, 60.0, 100.0];
        for (i, &l) in levels.iter().enumerate() {
            for j in 0..14 {
                let jitter = (((i * 7 + j * 3) % 5) as f32 - 2.0) * 0.8;
                let spike = if j == 6 && i % 2 == 0 { 25.0 } else { 0.0 };
                sig.push(l + jitter + spike);
            }
        }
        let p = RawParams { window_length1: 3, window_length2: 6, threshold1: 5.5, threshold2: 3.5,
                            peak_height: 0.2, ..RawParams::default() };
        let plain = detect_events(&sig, &p);
        assert_eq!(detect_events_merged(&sig, &p, 0.0), plain);
        let merged = detect_events_merged(&sig, &p, 10.0);
        assert_eq!(merged.len(), 1, "a huge threshold merges everything");
        let close = merge_close_events(&[0.0, 0.1, 0.15, 1.0, 1.05, -0.5], 0.2);
        assert_eq!(close.len(), 3);
        assert!((close[0] - 0.083333).abs() < 1e-4 && (close[1] - 1.025).abs() < 1e-4 && close[2] == -0.5);
    }

    use super::*;

    /// Reference outputs of RawHash2's C quantizer and sketcher (`rsketch.c` at commit 92366af, compiled with
    /// `-ffp-contract=off -fno-fast-math` and driven over the listed signal values; f32 values as their bit patterns).
    const RAW_HASH2_GOLDEN: &str = r#"# RawHash2 raw-sketch parity fixture v1
# PARAM name value
# QUANT signal_f32_bits expected_bucket
# EVENT stream index signal_f32_bits
# CASE name bucket_diff window events_per_hash output_count
# SKETCH case position hash
PARAM rawhash2_commit 92366afb09df9fa791e7705414e1372f487912ed
PARAM quant_bit 3
PARAM n_buckets 7
PARAM k_per_event 6
PARAM fine_min_bits bfcccccd
PARAM fine_max_bits 3fcccccd
PARAM fine_range_bits 3f400000
QUANT c0400001 0
QUANT c0400000 0
QUANT c03fffff 0
QUANT bfccccce 0
QUANT bfcccccd 1
QUANT bfcccccc 1
QUANT bf75c290 1
QUANT bf75c28f 2
QUANT bf75c28e 2
QUANT bea3d70b 2
QUANT bea3d70a 2
QUANT bea3d709 3
QUANT 3ea3d709 3
QUANT 3ea3d70a 4
QUANT 3ea3d70b 4
QUANT 3f75c28e 4
QUANT 3f75c28f 4
QUANT 3f75c290 4
QUANT 3fcccccc 5
QUANT 3fcccccd 5
QUANT 3fccccce 6
QUANT 403fffff 6
QUANT 40400000 6
QUANT 40400001 6
QUANT ff800000 0
QUANT 7f800000 6
EVENT main 0 c0000000
EVENT main 1 c0000000
EVENT main 2 bfa66666
EVENT main 3 bfa66666
EVENT main 4 bf333333
EVENT main 5 bf333333
EVENT main 6 bdcccccd
EVENT main 7 bdcccccd
EVENT main 8 3f000000
EVENT main 9 3f000000
EVENT main 10 3f99999a
EVENT main 11 3f99999a
EVENT main 12 40000000
EVENT main 13 40000000
EVENT ties 0 00000000
EVENT ties 1 00000000
EVENT ties 2 00000000
EVENT ties 3 00000000
EVENT ties 4 00000000
EVENT ties 5 00000000
EVENT ties 6 00000000
EVENT ties 7 00000000
CASE dense_diff_minus1 -1 0 3 12
SKETCH dense_diff_minus1 0 3079993582
SKETCH dense_diff_minus1 1 4103259927
SKETCH dense_diff_minus1 2 2484412683
SKETCH dense_diff_minus1 3 3507638712
SKETCH dense_diff_minus1 4 1892933000
SKETCH dense_diff_minus1 5 768683448
SKETCH dense_diff_minus1 6 3446239576
SKETCH dense_diff_minus1 7 2319514396
SKETCH dense_diff_minus1 8 2850539067
SKETCH dense_diff_minus1 9 3873719402
SKETCH dense_diff_minus1 10 2260549206
SKETCH dense_diff_minus1 11 3282807565
CASE dense_diff_zero 0 0 3 5
SKETCH dense_diff_zero 0 2888324037
SKETCH dense_diff_zero 2 2292723904
SKETCH dense_diff_zero 4 1702667295
SKETCH dense_diff_zero 6 3252034901
SKETCH dense_diff_zero 8 2661514514
CASE window_diff_minus1 -1 4 3 5
SKETCH window_diff_minus1 2 2484412683
SKETCH window_diff_minus1 4 1892933000
SKETCH window_diff_minus1 5 768683448
SKETCH window_diff_minus1 7 2319514396
SKETCH window_diff_minus1 10 2260549206
CASE window_diff_zero 0 4 3 1
SKETCH window_diff_zero 4 1702667295
CASE window_ties_diff_minus1 -1 4 2 6
SKETCH window_ties_diff_minus1 0 3719860965
SKETCH window_ties_diff_minus1 1 3719860965
SKETCH window_ties_diff_minus1 3 3719860965
SKETCH window_ties_diff_minus1 4 3719860965
SKETCH window_ties_diff_minus1 5 3719860965
SKETCH window_ties_diff_minus1 6 3719860965
"#;

    fn golden_sketch(case_name: &str) -> Vec<(usize, u64)> {
        RAW_HASH2_GOLDEN
            .lines()
            .filter_map(|line| {
                let mut fields = line.split_whitespace();
                if fields.next()? != "SKETCH" || fields.next()? != case_name {
                    return None;
                }
                Some((
                    fields.next().unwrap().parse().unwrap(),
                    fields.next().unwrap().parse().unwrap(),
                ))
            })
            .collect()
    }

    #[test]
    fn hash64_known_values() {
        // Distinct keys give distinct hashes within the mask (the exact values
        // are checked against RawHash2 by the golden sketches).
        let mask = (1u64 << 32) - 1;
        let h0 = hash64(0, mask);
        let h1 = hash64(1, mask);
        let h2 = hash64(2, mask);
        assert_ne!(h0, h1);
        assert_ne!(h1, h2);
        assert!(h0 <= mask);
        assert!(h1 <= mask);
    }

    #[test]
    fn dynamic_quantize_three_zones() {
        // RawParams::default() values: n_buckets=7, fine_min=-1.6, fine_max=1.6,
        // fine_range=0.75 → n_fine=5, n_low=1, n_high=1.
        // Low zone: [-3, -1.6) gets 1 bucket -> id 0.
        // Fine zone: [-1.6, 1.6] gets 5 buckets -> ids 1..=5.
        // High zone: (1.6, 3] gets 1 bucket -> id 6.
        let n_buckets = 7;
        let (fmin, fmax, frange) = (-1.6_f32, 1.6_f32, 0.75_f32);
        assert_eq!(dynamic_quantize(-3.0, fmin, fmax, frange, n_buckets), 0);
        assert_eq!(dynamic_quantize(-2.0, fmin, fmax, frange, n_buckets), 0);
        // Lower fine edge: t=0 → first bucket of the fine zone.
        assert_eq!(dynamic_quantize(-1.6, fmin, fmax, frange, n_buckets), 1);
        // Center: t=0.5 → trunc(0.5*5)=2 → n_low + 2 = 3.
        assert_eq!(dynamic_quantize(0.0, fmin, fmax, frange, n_buckets), 3);
        // Upper fine edge: clamp to last bucket of the fine zone (n_low + n_fine - 1 = 5).
        assert_eq!(dynamic_quantize(1.6, fmin, fmax, frange, n_buckets), 5);
        assert_eq!(dynamic_quantize(2.0, fmin, fmax, frange, n_buckets), 6);
        assert_eq!(dynamic_quantize(3.0, fmin, fmax, frange, n_buckets), 6);
    }

    #[test]
    fn dynamic_quantize_matches_rawhash2_boundaries() {
        assert!(RAW_HASH2_GOLDEN
            .contains("PARAM rawhash2_commit 92366afb09df9fa791e7705414e1372f487912ed"));
        for line in RAW_HASH2_GOLDEN.lines() {
            let mut fields = line.split_whitespace();
            if fields.next() != Some("QUANT") {
                continue;
            }
            let bits = u32::from_str_radix(fields.next().unwrap(), 16).unwrap();
            let expected = fields.next().unwrap().parse::<u32>().unwrap();
            assert_eq!(
                dynamic_quantize(f32::from_bits(bits), -1.6, 1.6, 0.75, 7),
                expected,
                "quantizer mismatch for f32 bits {:08x}",
                bits
            );
        }
    }

    #[test]
    fn normalize_signal_drops_outliers() {
        // 99 zeros + 1 big outlier. mean≈1, std≈9.95, so the outlier lands
        // near +9.95σ and gets dropped; the 99 near-zero samples are within ±3σ.
        let mut sig: Vec<f32> = vec![0.0; 99];
        sig.push(100.0);
        let mut st = NormalizeState::default();
        let norm = normalize_signal(&sig, &mut st);
        assert_eq!(norm.len(), sig.len() - 1);
    }

    #[test]
    fn normalize_signal_keeps_inliers() {
        // Gaussian-ish signal with no outliers — everything should be kept.
        let sig: Vec<f32> = (0..100).map(|i| ((i as f32) * 0.1).sin()).collect();
        let mut st = NormalizeState::default();
        let norm = normalize_signal(&sig, &mut st);
        assert_eq!(norm.len(), sig.len());
    }

    #[test]
    fn sketch_reg_emits_expected_count() {
        // Both merge filters disabled — every event contributes to the hash
        // stream, so the count is purely a function of `e` and input length.
        let p = RawParams {
            e: 3,
            w: 0,
            quant_bit: 3,
            bucket_diff: -1,
            sig_diff: -1.0,
            ..RawParams::default()
        };
        let events: Vec<f32> = (0..10).map(|i| i as f32 * 0.3 - 1.5).collect();
        let sk = sketch_signal(&events, &p);
        // n events = 10, e = 3 -> 8 hashes. RawHash2 locates each packed
        // word at its first event, so the coordinates are 0..=7.
        assert_eq!(sk.len(), 8);
        assert_eq!(
            sk.iter().map(|&(pos, _)| pos).collect::<Vec<_>>(),
            (0..8).collect::<Vec<_>>()
        );
    }

    #[test]
    fn quantize_filter_and_pack_are_explicit_stages() {
        let events = vec![
            -2.0, -2.0, -1.3, -1.3, -0.7, -0.7, -0.1, -0.1, 0.5, 0.5, 1.2, 1.2, 2.0, 2.0,
        ];
        let base = RawParams {
            e: 3,
            w: 0,
            quant_bit: 3,
            n_buckets: 7,
            bucket_diff: -1,
            sig_diff: -1.0,
            fine_min: -1.6,
            fine_max: 1.6,
            fine_range: 0.75,
            ..RawParams::default()
        };

        let symbols = quantize_events(&events, &base);
        assert_eq!(
            symbols.values,
            vec![0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]
        );
        assert_eq!(
            symbols.event_positions,
            (0..events.len()).collect::<Vec<_>>()
        );

        let unfiltered = filter_event_symbols(&symbols, &events, &base);
        assert_eq!(unfiltered, symbols);

        let collapsed_params = RawParams {
            bucket_diff: 0,
            ..base.clone()
        };
        let collapsed = filter_event_symbols(&symbols, &events, &collapsed_params);
        assert_eq!(collapsed.values, vec![0, 1, 2, 3, 4, 5, 6]);
        assert_eq!(collapsed.event_positions, vec![0, 2, 4, 6, 8, 10, 12]);

        let custom = RawSymbols::new(vec![1, 2, 3, 4], vec![3, 7, 10, 20]);
        assert_eq!(
            pack_dense_emers(&custom, &base),
            vec![(3, 2_292_723_904), (7, 1_702_667_295)]
        );
        assert_eq!(
            RawSymbols::from_values(vec![4, 5]).event_positions,
            vec![0, 1]
        );
    }

    #[test]
    #[should_panic(expected = "outside quantizer alphabet")]
    fn corrected_symbols_must_fit_the_quantizer_alphabet() {
        let p = RawParams {
            e: 2,
            w: 0,
            quant_bit: 3,
            n_buckets: 7,
            ..RawParams::default()
        };
        let symbols = RawSymbols::from_values(vec![0, 7]);
        let _ = pack_dense_emers(&symbols, &p);
    }

    #[test]
    fn rawhash2_sketch_parity() {
        // The golden sketches come from RawHash2's C sketcher
        // against RawHash2 at commit 92366afb09df9fa791e7705414e1372f487912ed.
        let events = vec![
            -2.0, -2.0, -1.3, -1.3, -0.7, -0.7, -0.1, -0.1, 0.5, 0.5, 1.2, 1.2, 2.0, 2.0,
        ];
        let base = RawParams {
            e: 3,
            w: 0,
            quant_bit: 3,
            n_buckets: 7,
            bucket_diff: -1,
            sig_diff: -1.0,
            fine_min: -1.6,
            fine_max: 1.6,
            fine_range: 0.75,
            ..RawParams::default()
        };

        assert_eq!(
            sketch_signal(&events, &base),
            golden_sketch("dense_diff_minus1")
        );

        let dense_diff_zero_params = RawParams {
            bucket_diff: 0,
            ..base.clone()
        };
        assert_eq!(
            sketch_signal(&events, &dense_diff_zero_params),
            golden_sketch("dense_diff_zero")
        );

        let windowed_params = RawParams {
            w: 4,
            ..base.clone()
        };
        assert_eq!(
            sketch_signal(&events, &windowed_params),
            golden_sketch("window_diff_minus1")
        );

        let windowed_diff_zero_params = RawParams {
            w: 4,
            bucket_diff: 0,
            ..base
        };
        assert_eq!(
            sketch_signal(&events, &windowed_diff_zero_params),
            golden_sketch("window_diff_zero")
        );
    }

    #[test]
    fn rawhash2_equal_hash_window_ties_are_preserved() {
        let p = RawParams {
            e: 2,
            w: 4,
            quant_bit: 3,
            n_buckets: 7,
            bucket_diff: -1,
            sig_diff: -1.0,
            fine_min: -1.6,
            fine_max: 1.6,
            fine_range: 0.75,
            ..RawParams::default()
        };
        let got = sketch_signal(&vec![0.0; 8], &p);
        assert_eq!(got, golden_sketch("window_ties_diff_minus1"));
    }



    #[test]
    fn detect_events_merged_spans_gives_the_same_events() {
        // A synthetic staircase with noise and outliers: the spans version must give detect_events_merged's values,
        // spans inside the signal, in order, one per event.
        let mut sig = Vec::new();
        let mut x: u64 = 12345;
        let mut rnd = || {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((x >> 33) as f32 / (1u64 << 31) as f32) - 0.5
        };
        for level in 0..400 {
            let base = 80.0 + 30.0 * (((level * 7919) % 13) as f32 / 13.0);
            let len = 6 + (level * 31) % 9;
            for _ in 0..len {
                sig.push(base + 3.0 * rnd());
            }
            if level % 50 == 0 {
                sig.push(900.0); // an outlier the normalisation drops
            }
        }
        let p = RawParams::default();
        for &em in &[0.0f32, 0.3] {
            let plain = detect_events_merged(&sig, &p, em);
            let (ev, spans) = detect_events_merged_spans(&sig, &p, em);
            assert_eq!(plain, ev, "event_merge {}", em);
            assert_eq!(ev.len(), spans.len());
            assert!(!ev.is_empty());
            for w in spans.windows(2) {
                assert!(w[0].0 < w[0].1 && w[0].1 <= w[1].0, "{:?}", w);
            }
            assert!(spans.last().unwrap().1 as usize <= sig.len());
        }
    }


    #[test]
    fn sketch_reg_positions_monotone() {
        let p = RawParams {
            e: 3,
            w: 0,
            quant_bit: 3,
            bucket_diff: -1,
            sig_diff: -1.0,
            ..RawParams::default()
        };
        let events: Vec<f32> = (0..30).map(|i| ((i as f32) * 0.21 - 2.0).sin()).collect();
        let sk = sketch_signal(&events, &p);
        assert!(!sk.is_empty());
        for w in sk.windows(2) {
            assert!(w[0].0 <= w[1].0);
        }
    }

    #[test]
    fn bucket_diff_skips_same_bucket_runs() {
        // 30-event signal where every value falls in the same fine-zone bucket.
        // With bucket_diff = 0, only the first event survives the filter, so
        // the e=3 ring never fills and no hash is emitted.
        let p = RawParams {
            e: 3,
            w: 0,
            quant_bit: 3,
            bucket_diff: 0,
            sig_diff: -1.0,
            ..RawParams::default()
        };
        let plateau: Vec<f32> = vec![0.0; 30];
        let sk = sketch_signal(&plateau, &p);
        assert!(
            sk.is_empty(),
            "bucket_diff=0 should collapse a same-bucket plateau"
        );

        // Toggling bucket_diff off restores dense behavior.
        let p_off = RawParams {
            bucket_diff: -1,
            ..p
        };
        let sk_off = sketch_signal(&plateau, &p_off);
        // 30 events, e=3 → 30 - (e-1) = 28 emitted hashes.
        assert_eq!(sk_off.len(), 28);
    }

    #[test]
    fn corrected_symbol_bucket_filter_preserves_fresh_coordinates() {
        let symbols = RawSymbols::from_values(vec![1, 1, 2, 2, 4, 3, 3]);
        let params = RawParams {
            bucket_diff: 0,
            sig_diff: -1.0,
            ..RawParams::default()
        };
        let filtered = filter_symbol_buckets(&symbols, &params);
        assert_eq!(filtered.values, vec![1, 2, 4, 3]);
        assert_eq!(filtered.event_positions, vec![0, 2, 4, 5]);

        let disabled = RawParams {
            bucket_diff: -1,
            ..params
        };
        assert_eq!(filter_symbol_buckets(&symbols, &disabled), symbols);
    }

    #[test]
    fn gapped_plateaus_preserve_exact_dense_and_windowed_endpoints() {
        let symbols = RawSymbols::from_values(vec![1, 1, 2, 2, 4, 4, 3]);
        let dense_params = RawParams {
            e: 3,
            w: 0,
            quant_bit: 3,
            n_buckets: 7,
            bucket_diff: 0,
            sig_diff: -1.0,
            ..RawParams::default()
        };
        let filtered = filter_symbol_buckets(&symbols, &dense_params);
        assert_eq!(filtered.event_positions, vec![0, 2, 4, 6]);

        let spanned = pack_dense_emers_with_spans(&filtered, &dense_params);
        assert_eq!(
            spanned
                .iter()
                .map(|token| (token.start, token.end))
                .collect::<Vec<_>>(),
            vec![(0, 5), (2, 7)]
        );
        assert_eq!(
            pack_dense_emers(&filtered, &dense_params),
            spanned
                .iter()
                .map(|token| (token.start, token.hash))
                .collect::<Vec<_>>()
        );
        // A two-token graph window covers [0, 7), not [0, 2 + e) = [0, 5).
        assert_eq!(spanned.first().unwrap().start, 0);
        assert_eq!(spanned.last().unwrap().end, 7);

        let windowed_params = RawParams {
            w: 2,
            ..dense_params
        };
        let exact_windowed = sketch_symbols_with_spans(&filtered, &windowed_params);
        assert_eq!(
            sketch_symbols(&filtered, &windowed_params),
            exact_windowed
                .iter()
                .map(|token| (token.start, token.hash))
                .collect::<Vec<_>>()
        );
        for token in exact_windowed {
            let dense_token = spanned
                .iter()
                .find(|candidate| candidate.start == token.start && candidate.hash == token.hash)
                .unwrap();
            assert_eq!(token.end, dense_token.end);
        }
    }

    #[test]
    fn sig_diff_layers_on_top_of_bucket_diff() {
        // Construct a signal where consecutive events alternate buckets (so
        // bucket_diff = 0 lets every event through) but the raw signal deltas
        // stay tiny (so a sig_diff above that delta drops them all after the
        // first kept event).
        //
        // Bucket layout with defaults (n_buckets=7, fine_min=-1.6, fine_max=1.6,
        // fine_range=0.75): one low bucket, then five fine buckets of width 0.64,
        // so bucket 3 is [-0.32, 0.32) and the 3/4 boundary is at signal 0.32
        // (t=0.6, trunc(0.6*5)=3). Pick values straddling it.
        let mut events = Vec::new();
        for _ in 0..20 {
            events.push(0.31_f32);
            events.push(0.33_f32);
        }
        let p_bucket_only = RawParams {
            e: 3,
            w: 0,
            quant_bit: 3,
            bucket_diff: 0,
            sig_diff: -1.0,
            ..RawParams::default()
        };
        let sk_bucket = sketch_signal(&events, &p_bucket_only);
        // bucket_diff=0 keeps the alternating pattern (buckets differ by 1).
        assert!(!sk_bucket.is_empty());

        let p_layered = RawParams {
            sig_diff: 0.05,
            ..p_bucket_only
        };
        let sk_layered = sketch_signal(&events, &p_layered);
        // |0.33 - 0.31| = 0.02 < 0.05, so every event after the first is
        // dropped by the sig_diff layer regardless of bucket transitions.
        assert!(
            sk_layered.is_empty(),
            "sig_diff should drop events the bucket filter let through"
        );
    }
}
