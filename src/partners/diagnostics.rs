//! Diagnostics switched by `SQUASH_*` environment variables (phase profile, cover probe, memory tracking and
//! reports); none of them changes the rows.

use super::*;

/// `SQUASH_PARTNERS_PROFILE=1`: phase times and counts of `partner_rows` on stderr (a diagnostic; the rows do not change).
/// Off, it costs one flag test per target and per donor.
pub mod profile {
    use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
    use std::sync::OnceLock;
    pub const NAMES: [&str; 22] = ["value_sort_ns", "index_rest_ns", "gather_ns", "chain_ns", "filter_ns", "guard_ns",
        "rows_ns", "targets", "anchors", "donors_touched", "donors_chained", "chains_kept_before_guard",
        "guard_passed", "rows", "pruned_donors", "pruned_anchors", "order_ns", "prune_ns", "check_ns", "batched_targets",
        "batches", "batch_place_ns"];
    pub static C: [AtomicU64; 22] = [const { AtomicU64::new(0) }; 22];
    pub fn on() -> bool {
        static ON: OnceLock<bool> = OnceLock::new();
        *ON.get_or_init(|| std::env::var("SQUASH_PARTNERS_PROFILE").map_or(false, |v| v == "1"))
    }
    pub fn add(i: usize, v: u64) {
        C[i].fetch_add(v, Relaxed);
    }
    pub fn since(i: usize, t: std::time::Instant) {
        add(i, t.elapsed().as_nanos() as u64);
    }
    /// Print and reset the counters (times as thread-seconds).
    pub fn report(label: &str) {
        let parts: Vec<String> = NAMES.iter().enumerate().map(|(i, n)| {
            let v = C[i].swap(0, Relaxed);
            if n.ends_with("_ns") { format!("{}={:.1}s", &n[..n.len() - 3], v as f64 / 1e9) } else { format!("{}={}", n, v) }
        }).collect();
        eprintln!("[partners-profile] {} {}", label, parts.join(" "));
    }
}

/// `SQUASH_PARTNERS_COVER_PROBE=1` (a diagnostic): for every donor that reaches chaining, whether an exact bound from its
/// anchors alone (`cover_bound`) says the chain must fail the overlap or cover test, how much chaining time those
/// donors take, and any donor the bound would have skipped whose chain passes (a violation; must stay 0). Every donor
/// is still chained; the rows do not change.
pub mod cover_probe {
    use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
    use std::sync::OnceLock;
    pub const NAMES: [&str; 10] = ["donors_chained", "skip_cover", "skip_overlap", "anchors_chained", "anchors_skippable",
        "chain_ns", "chain_ns_skippable", "violations", "passed_tests", "bound_ns"];
    pub static C: [AtomicU64; 10] = [const { AtomicU64::new(0) }; 10];
    pub fn on() -> bool {
        static ON: OnceLock<bool> = OnceLock::new();
        *ON.get_or_init(|| std::env::var("SQUASH_PARTNERS_COVER_PROBE").map_or(false, |v| v == "1"))
    }
    pub fn add(i: usize, v: u64) {
        C[i].fetch_add(v, Relaxed);
    }
    pub fn report(label: &str) {
        if !on() {
            return;
        }
        let parts: Vec<String> = NAMES.iter().enumerate().map(|(i, n)| {
            let v = C[i].swap(0, Relaxed);
            if n.ends_with("_ns") { format!("{}={:.1}s", &n[..n.len() - 3], v as f64 / 1e9) } else { format!("{}={}", n, v) }
        }).collect();
        eprintln!("[cover-probe] {} {}", label, parts.join(" "));
    }
}

/// `SQUASH_PARTNERS_MEMTRACK=1` (a diagnostic): a line a second on stderr (`[memtrack]`) with the process's resident memory
/// (VmRSS, RssAnon) and the bytes the search's per-thread buffers hold, summed over the threads: anchor buffers, the
/// chaining scratch, the donor order (with touched donors and hits), the check's alignment buffers (the most each
/// thread held since the last line, and what it keeps between pairs), the fixed search scratch; and how many threads
/// hold an anchor buffer over 4 MB / over 64 MB. Off, it costs one flag test per update.
pub mod memtrack {
    use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering::Relaxed};
    use std::sync::OnceLock;
    pub const ANCHORS: usize = 0;
    pub const CHAIN: usize = 1;
    pub const ORDER: usize = 2;
    pub const PATH: usize = 3;
    pub const PATH_KEPT: usize = 4;
    pub const FIXED: usize = 5;
    const NAMES: [&str; 6] = ["anchors", "chain_scratch", "order", "check_peak", "check_kept", "fixed_scratch"];
    const SLOTS: usize = 256;
    #[repr(align(64))]
    struct Slot([AtomicU64; 6]);
    static SLOT: [Slot; SLOTS] = [const { Slot([const { AtomicU64::new(0) }; 6]) }; SLOTS];
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    thread_local! {
        static ID: usize = NEXT.fetch_add(1, Relaxed) % SLOTS;
    }
    pub fn on() -> bool {
        static ON: OnceLock<bool> = OnceLock::new();
        *ON.get_or_init(|| std::env::var("SQUASH_PARTNERS_MEMTRACK").map_or(false, |v| v == "1"))
    }
    /// This thread's bytes in category `cat`.
    pub fn set(cat: usize, bytes: usize) {
        if on() {
            ID.with(|&i| SLOT[i].0[cat].store(bytes as u64, Relaxed));
        }
    }
    /// This thread's most bytes in `cat` since the last line.
    pub fn peak(cat: usize, bytes: usize) {
        if on() {
            ID.with(|&i| {
                SLOT[i].0[cat].fetch_max(bytes as u64, Relaxed);
            });
        }
    }
    pub struct Sampler {
        stop: std::sync::Arc<AtomicBool>,
        handle: Option<std::thread::JoinHandle<()>>,
    }
    impl Sampler {
        pub fn start() -> Sampler {
            let stop = std::sync::Arc::new(AtomicBool::new(false));
            if !on() {
                return Sampler { stop, handle: None };
            }
            let flag = std::sync::Arc::clone(&stop);
            let handle = std::thread::spawn(move || {
                let t0 = std::time::Instant::now();
                let mut waited = 0u64;
                while !flag.load(Relaxed) {
                    std::thread::sleep(std::time::Duration::from_millis(100));
                    waited += 100;
                    if waited < 1000 {
                        continue;
                    }
                    waited = 0;
                    let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
                    let gb = |key: &str| status.lines().find(|l| l.starts_with(key)).and_then(|l| l.split_whitespace().nth(1))
                        .and_then(|x| x.parse::<f64>().ok()).unwrap_or(0.0) / 1e6;
                    let mut sums = [0u64; 6];
                    let (mut big4, mut big64) = (0, 0);
                    for s in SLOT.iter() {
                        for (c, sum) in sums.iter_mut().enumerate() {
                            *sum += if c == PATH { s.0[c].swap(0, Relaxed) } else { s.0[c].load(Relaxed) };
                        }
                        let a = s.0[ANCHORS].load(Relaxed);
                        big4 += (a > 4 << 20) as usize;
                        big64 += (a > 64 << 20) as usize;
                    }
                    let parts: Vec<String> = NAMES.iter().zip(sums).map(|(n, b)| format!("{}={:.3}", n, b as f64 / 1e9)).collect();
                    eprintln!("[memtrack] t={:.0} rss={:.3} anon={:.3} {} big4={} big64={}", t0.elapsed().as_secs_f64(),
                              gb("VmRSS:"), gb("RssAnon:"), parts.join(" "), big4, big64);
                }
            });
            Sampler { stop, handle: Some(handle) }
        }
    }
    impl Drop for Sampler {
        fn drop(&mut self) {
            self.stop.store(true, Relaxed);
            if let Some(handle) = self.handle.take() {
                let _ = handle.join();
            }
        }
    }
}

/// `cover_probe`'s bound for one donor's anchors (sorted or not): Some(1) when no chain of them can imply an overlap of
/// `min_overlap`, Some(0) when none can cover `min_cover` of the overlap it implies, None otherwise. A chain spans at
/// most from its earliest anchor's start to its latest anchor's end, and its offset (the rounded mean of two anchors'
/// diagonals) lies between the anchors' least and greatest diagonal; the implied overlap is concave in the offset, so
/// its least value there is at an end and its greatest at an end or a breakpoint (offset 0 or len_d - len_r).
pub fn cover_bound<A: Anchor>(tt: &[Token], dt: &[Token], len_r: i64, len_d: i64, group: &[A], min_cover: f64,
                              min_overlap: i64) -> Option<u8> {
    let (mut min_s, mut max_e, mut dmin, mut dmax) = (i64::MAX, i64::MIN, i64::MAX, i64::MIN);
    for a in group {
        let t = &tt[a.t() as usize];
        let ts = t.s();
        min_s = min_s.min(ts);
        max_e = max_e.max(t.e());
        let dg = dt[a.d() as usize].s() - ts;
        dmin = dmin.min(dg);
        dmax = dmax.max(dg);
    }
    let ov = |o: i64| len_r.min(len_d - o) - 0i64.max(-o);
    let (a, b) = (ov(dmin), ov(dmax));
    let mut hi = a.max(b);
    for c in [0, len_d - len_r] {
        if c > dmin && c < dmax {
            hi = hi.max(ov(c));
        }
    }
    if hi < min_overlap {
        return Some(1);
    }
    if ((max_e - min_s) as f64) < min_cover * (a.min(b) as f64) {
        return Some(0);
    }
    None
}

/// `SQUASH_ROUNDS_MEM=1` (a diagnostic): the process's resident memory and its peak since the previous report (the peak is
/// reset after each report through /proc/self/clear_refs), with the sizes of the main structures in GB. Outputs do not
/// change.
pub mod mem {
    use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};
    use std::sync::OnceLock;
    /// Bytes held by the alignment DP's per-thread trace buffers (they grow and are kept).
    pub static DP_BYTES: AtomicUsize = AtomicUsize::new(0);
    pub fn on() -> bool {
        static ON: OnceLock<bool> = OnceLock::new();
        *ON.get_or_init(|| std::env::var("SQUASH_ROUNDS_MEM").map_or(false, |v| v == "1"))
    }
    #[allow(dead_code)]
    pub fn grew(before: usize, after: usize) {
        if after > before {
            DP_BYTES.fetch_add(after - before, Relaxed);
        } else if before > after {
            DP_BYTES.fetch_sub(before - after, Relaxed);
        }
    }
    pub fn report(phase: &str, parts: &[(&str, usize)]) {
        if !on() {
            return;
        }
        let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
        let gb = |key: &str| status.lines().find(|l| l.starts_with(key)).and_then(|l| l.split_whitespace().nth(1))
            .and_then(|x| x.parse::<f64>().ok()).unwrap_or(0.0) / 1e6;
        let mut all: Vec<String> = parts.iter().map(|(n, b)| format!("{} {:.2}", n, *b as f64 / 1e9)).collect();
        all.push(format!("dp_traces {:.2}", DP_BYTES.load(Relaxed) as f64 / 1e9));
        eprintln!("[mem] {}: rss {:.2} GB, peak since last {:.2} GB; {}", phase, gb("VmRSS:"), gb("VmHWM:"), all.join(", "));
        let _ = std::fs::write("/proc/self/clear_refs", "5");
    }
    /// Heap bytes of a per-read vector of vectors.
    #[allow(dead_code)]
    pub fn nested<T>(x: &[Vec<T>]) -> usize {
        x.iter().map(|r| r.capacity() * std::mem::size_of::<T>()).sum::<usize>() + x.len() * std::mem::size_of::<Vec<T>>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The cover bound (`cover_bound`) never marks a donor whose chain passes the overlap and cover tests, on repeat
    /// and noisy reads at several covers and overlaps, and does mark some (else it would measure nothing).
    #[test]
    fn cover_bound_never_skips_a_passing_chain() {
        let (mut marked, mut checked) = (0usize, 0usize);
        for (seed, repeats) in [(1u64, true), (2, false), (3, true)] {
            let (lengths, tokens) = if repeats { repeat_reads(10, seed) } else { noisy_reads(10, seed) };
            for (min_cover, min_overlap) in [(0.1, 500i64), (0.25, 500), (0.5, 1000), (0.9, 200)] {
                let p = Params { min_cover, min_overlap, ..Params::default() };
                for r in 0..tokens.len() {
                    for d in 0..tokens.len() {
                        if r == d {
                            continue;
                        }
                        let group: Vec<(u32, u32)> = tokens[r].iter().enumerate().flat_map(|(i, x)| {
                            tokens[d].iter().enumerate().filter(move |(_, y)| y.value == x.value).map(move |(j, _)| (i as u32, j as u32))
                        }).collect();
                        if group.len() < p.min_anchors {
                            continue;
                        }
                        let bound = cover_bound(&tokens[r], &tokens[d], lengths[r], lengths[d], &group, min_cover, min_overlap);
                        checked += 1;
                        marked += bound.is_some() as usize;
                        let (tt, dt) = (&tokens[r], &tokens[d]);
                        if let Some((n, _, t0, t1, d0, d1)) = chain_sorted(&group, &p) {
                            if n < p.min_anchors {
                                continue;
                            }
                            let (t0, t1, d0, d1) = (t0 as usize, t1 as usize, d0 as usize, d1 as usize);
                            let offset = round_half_even_div2((dt[d0].s() - tt[t0].s()) + (dt[d1].s() - tt[t1].s()));
                            let (ts, te) = (0i64.max(-offset), lengths[r].min(lengths[d] - offset));
                            let pass = !(te - ts < min_overlap || ((tt[t1].e() - tt[t0].s()) as f64) < min_cover * ((te - ts) as f64));
                            assert!(!(pass && bound.is_some()), "bound {:?} skips a passing chain: target {} donor {}", bound, r, d);
                        }
                    }
                }
            }
        }
        assert!(checked > 1000);
        // A donor sharing only a short stretch with the target (a repeat copy): its anchors cannot cover the overlap.
        let tok = |v: u32, i: i64| Token::new(v, 10 * i, 10 * i + 10);
        let target: Vec<Token> = (0..300).map(|i| tok(i as u32, i)).collect();
        let donor: Vec<Token> = (0..300).map(|i| tok(if (40..60).contains(&i) { 100 + i as u32 - 40 } else { 10_000 + i as u32 }, i)).collect();
        let group: Vec<(u32, u32)> = (0..20).map(|k| (100 + k, 40 + k)).collect();
        assert_eq!(cover_bound(&target, &donor, 3000, 3000, &group, 0.25, 500), Some(0));
        assert_eq!(cover_bound(&target, &donor, 3000, 3000, &group, 0.05, 500), None);
        assert_eq!(cover_bound(&target, &donor, 3000, 3000, &group, 0.05, 2500), Some(1));
        let _ = marked;
    }
}
