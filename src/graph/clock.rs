//! The wall-clock and memory profile of the steps (`SQUASH_GRAPH_TIMING`, `SQUASH_GRAPH_MEM`); it never touches an output.

use std::fs;

/// Wall-clock profile (`SQUASH_GRAPH_TIMING=1`): one `graph-timing` line per phase on stderr; it never touches an output.
/// `SQUASH_GRAPH_MEM=1` also turns it on and appends the process's resident memory at the mark and its peak since the
/// previous mark of any clock (VmHWM, reset through `/proc/self/clear_refs`; phases written in parallel with the next
/// order's compute share one peak).
pub(super) struct PhaseClock {
    label: &'static str,
    last: std::time::Instant,
    last_cpu: f64,
    on: bool,
    mem: bool,
}

/// (VmRSS, VmHWM) in GB from /proc/self/status, for the memory profile only.
fn process_rss_hwm_gb() -> (f64, f64) {
    let mut out = (0.0, 0.0);
    if let Ok(text) = fs::read_to_string("/proc/self/status") {
        for line in text.lines() {
            let kb = |v: &str| v.split_whitespace().next().and_then(|x| x.parse::<f64>().ok()).unwrap_or(0.0) / 1048576.0;
            if let Some(v) = line.strip_prefix("VmRSS:") {
                out.0 = kb(v);
            } else if let Some(v) = line.strip_prefix("VmHWM:") {
                out.1 = kb(v);
            }
        }
    }
    out
}

/// Reset the process's peak resident memory (VmHWM) to its current size.
fn reset_peak_rss() {
    let _ = fs::write("/proc/self/clear_refs", "5");
}

/// Process CPU seconds (user + system), for the timing profile only.
fn process_cpu_seconds() -> f64 {
    fs::read_to_string("/proc/self/stat")
        .ok()
        .and_then(|text| {
            let fields: Vec<&str> = text.rsplit(')').next()?.split_whitespace().collect();
            let ticks = fields.get(11)?.parse::<f64>().ok()? + fields.get(12)?.parse::<f64>().ok()?;
            Some(ticks / 100.0)
        })
        .unwrap_or(0.0)
}

impl PhaseClock {
    pub(super) fn new(label: &'static str) -> Self {
        let mem = std::env::var_os("SQUASH_GRAPH_MEM").map_or(false, |v| v == "1");
        let on = mem || std::env::var_os("SQUASH_GRAPH_TIMING").map_or(false, |v| v == "1");
        Self {
            label,
            last: std::time::Instant::now(),
            last_cpu: if on { process_cpu_seconds() } else { 0.0 },
            on,
            mem,
        }
    }
    pub(super) fn mark(&mut self, phase: &str) {
        let now = std::time::Instant::now();
        if self.on {
            let cpu = process_cpu_seconds();
            let mem = if self.mem {
                let (rss, hwm) = process_rss_hwm_gb();
                reset_peak_rss();
                format!("\trss_gb {:.2}\tpeak_gb {:.2}", rss, hwm)
            } else {
                String::new()
            };
            eprintln!(
                "graph-timing\t{}\t{}\t{:.3}\tcpu {:.1}{}",
                self.label,
                phase,
                (now - self.last).as_secs_f64(),
                cpu - self.last_cpu,
                mem
            );
            self.last_cpu = cpu;
        }
        self.last = now;
    }
}
