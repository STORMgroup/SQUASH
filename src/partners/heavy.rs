//! Heavy targets apart (setting THREADS:HITS:BUDGET_GB[:MIN_MB], fixed to `SETTING`). Some reads carry a stretch where
//! nearly every token is one value (pore stalls, satellite runs); against another such read one donor has count x count
//! anchors (hundreds of millions), placed and chained as one unit (8 bytes an anchor plus the chaining scratch: several
//! GB on one thread), and such targets, the slowest, would otherwise meet at the end of the search. So `search_targets`
//! (in a search of all targets over all donors) first estimates each target's work (its tokens' occurrences in the
//! index); the targets with at least HITS go to THREADS threads of their own, heaviest first, the others to the
//! search's pool. On those threads `target_chains` takes the bytes of its largest placement (the anchors placed at once
//! and the chaining scratch of its largest donor) from a budget of BUDGET_GB before placing anything (needs under
//! MIN_MB, default 256, take nothing), waiting while the budget is full (a need larger than the budget runs alone), and
//! gives them back once the target's buffers are trimmed. Every target is computed the same way wherever it runs and
//! its rows land in its own slot, so the table and the verdicts do not depend on the split
//! (`heavy_targets_apart_keep_the_rows`).

use std::cell::{Cell, RefCell};
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};
use std::sync::{Arc, Condvar, Mutex};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Config {
    pub threads: usize,
    pub hits: u64,
    pub budget: usize,
    pub min_bytes: usize,
}

impl Config {
    /// THREADS:HITS:BUDGET_GB[:MIN_MB] (GB = 1e9 bytes, MB = 1e6).
    pub fn parse(s: &str) -> Option<Config> {
        let f: Vec<&str> = s.split(':').collect();
        if f.len() < 3 || f.len() > 4 {
            return None;
        }
        let threads = f[0].parse::<usize>().ok().filter(|&t| t > 0)?;
        let hits = f[1].parse::<f64>().ok().filter(|&h| h >= 0.0)? as u64;
        let budget = (f[2].parse::<f64>().ok().filter(|&g| g > 0.0)? * 1e9) as usize;
        let min_bytes = match f.get(3) {
            Some(m) => (m.parse::<f64>().ok().filter(|&m| m >= 0.0)? * 1e6) as usize,
            None => 256_000_000,
        };
        Some(Config { threads, hits, budget, min_bytes })
    }
}

/// The setting of every search that covers all targets and all donors (`search_targets`: no previous rows, no allowed
/// donors): always `SETTING`; no environment variable is read.
#[cfg(not(test))]
pub fn config() -> Option<Config> {
    static C: std::sync::OnceLock<Option<Config>> = std::sync::OnceLock::new();
    *C.get_or_init(|| Config::parse(SETTING))
}

/// The heavy-target setting (THREADS:HITS:BUDGET_GB): 4 threads, targets with >= 1e8 hits, an 8-GB budget.
pub const SETTING: &str = "4:1e8:8";

#[cfg(test)]
thread_local! {
    pub static TEST: Cell<Option<Config>> = const { Cell::new(None) };
}

#[cfg(test)]
pub fn config() -> Option<Config> {
    TEST.with(|c| c.get())
}

/// The bytes the heavy threads hold (`take` / `give`), with counts for the report.
pub struct Budget {
    state: Mutex<State>,
    freed: Condvar,
    cap: usize,
    min_bytes: usize,
}

#[derive(Default)]
struct State {
    used: usize,
    takes: u64,
    waits: u64,
    alone: u64,
    wait_ns: u128,
    peak: usize,
}

impl Budget {
    pub fn new(cap: usize, min_bytes: usize) -> Budget {
        Budget { state: Mutex::new(State::default()), freed: Condvar::new(), cap, min_bytes }
    }

    /// Takes `bytes` (none under `min_bytes`), waiting while they do not fit (alone, when larger than the budget);
    /// returns what it took.
    pub fn take(&self, bytes: usize) -> usize {
        if bytes < self.min_bytes {
            return 0;
        }
        let clock = std::time::Instant::now();
        let mut st = self.state.lock().expect("the heavy budget");
        let mut waited = false;
        while !(st.used == 0 || st.used + bytes <= self.cap) {
            waited = true;
            st = self.freed.wait(st).expect("the heavy budget");
        }
        st.used += bytes;
        st.takes += 1;
        st.waits += waited as u64;
        st.alone += (bytes > self.cap) as u64;
        if waited {
            st.wait_ns += clock.elapsed().as_nanos();
        }
        st.peak = st.peak.max(st.used);
        bytes
    }

    pub fn give(&self, bytes: usize) {
        if bytes == 0 {
            return;
        }
        self.state.lock().expect("the heavy budget").used -= bytes;
        self.freed.notify_all();
    }

    #[cfg(test)]
    pub fn used(&self) -> usize {
        self.state.lock().expect("the heavy budget").used
    }

    fn report(&self) -> String {
        let st = self.state.lock().expect("the heavy budget");
        format!("{} takes, {} waited ({:.1} s), {} alone, peak {:.2} GB held", st.takes, st.waits, st.wait_ns as f64 / 1e9,
                st.alone, st.peak as f64 / 1e9)
    }
}

thread_local! {
    static CURRENT: RefCell<Option<Arc<Budget>>> = const { RefCell::new(None) };
    static HELD: Cell<usize> = const { Cell::new(0) };
}

/// Whether this thread is a heavy thread (`target_chains` then measures its need).
pub fn active() -> bool {
    CURRENT.with(|c| c.borrow().is_some())
}

/// `target_chains`, on a heavy thread: takes `bytes` before placing, held until `release`.
pub fn take(bytes: usize) {
    CURRENT.with(|c| {
        if let Some(b) = &*c.borrow() {
            let t = b.take(bytes);
            HELD.with(|h| h.set(h.get() + t));
        }
    });
}

/// `search_targets`, after the target's buffers are trimmed: gives back what the thread holds.
pub fn release() {
    CURRENT.with(|c| {
        if let Some(b) = &*c.borrow() {
            b.give(HELD.with(|h| h.replace(0)));
        }
    });
}

/// `search_targets`'s search with heavy targets apart: `body(r)` for every target, the targets whose `estimate` is at
/// least `cfg.hits` on `cfg.threads` threads of their own in decreasing estimate (ties by ordinal), the others on the
/// current rayon pool; the results in target order.
pub fn run<T: Send>(n: usize, cfg: Config, estimate: &(dyn Fn(usize) -> u64 + Sync),
                    body: &(dyn Fn(usize) -> T + Sync)) -> Vec<T> {
    use rayon::prelude::*;
    let clock = std::time::Instant::now();
    let est: Vec<u64> = (0..n).into_par_iter().map(estimate).collect();
    let mut big: Vec<usize> = (0..n).filter(|&r| est[r] >= cfg.hits).collect();
    big.sort_by(|&a, &b| est[b].cmp(&est[a]).then(a.cmp(&b)));
    let mut mark = vec![false; n];
    for &r in &big {
        mark[r] = true;
    }
    let total: u64 = est.iter().sum();
    let heavy: u64 = big.iter().map(|&r| est[r]).sum();
    eprintln!("[partners] heavy targets: {} of {} with >= {} hits ({:.3e} of {:.3e} hits; largest {:.3e}), {} threads, \
               budget {:.1} GB, take from {:.0} MB; estimate {:.1} s", big.len(), n, cfg.hits, heavy as f64, total as f64,
              big.first().map_or(0, |&r| est[r]) as f64, cfg.threads, cfg.budget as f64 / 1e9, cfg.min_bytes as f64 / 1e6,
              clock.elapsed().as_secs_f64());
    let budget = Arc::new(Budget::new(cfg.budget, cfg.min_bytes));
    let slots: Vec<Mutex<Option<T>>> = (0..n).map(|_| Mutex::new(None)).collect();
    let cursor = AtomicUsize::new(0);
    let heavy_done = Mutex::new(None::<f64>);
    std::thread::scope(|s| {
        for _ in 0..cfg.threads.min(big.len()) {
            s.spawn(|| {
                CURRENT.with(|c| *c.borrow_mut() = Some(Arc::clone(&budget)));
                loop {
                    let i = cursor.fetch_add(1, Relaxed);
                    if i >= big.len() {
                        break;
                    }
                    let r = big[i];
                    let out = body(r);
                    *slots[r].lock().expect("a target's slot") = Some(out);
                }
                CURRENT.with(|c| *c.borrow_mut() = None);
                let mut d = heavy_done.lock().expect("the heavy threads' end");
                *d = Some(d.unwrap_or(0.0).max(clock.elapsed().as_secs_f64()));
            });
        }
        (0..n).into_par_iter().filter(|&r| !mark[r]).for_each(|r| {
            let out = body(r);
            *slots[r].lock().expect("a target's slot") = Some(out);
        });
        eprintln!("[partners] heavy targets apart: the other targets done at {:.1} s", clock.elapsed().as_secs_f64());
    });
    eprintln!("[partners] heavy targets apart: heavy threads done at {:.1} s; {}",
              heavy_done.lock().expect("the heavy threads' end").unwrap_or(0.0), budget.report());
    slots.into_iter().map(|s| s.into_inner().expect("a target's slot").expect("every target searched")).collect()
}

#[cfg(test)]
mod tests {
    use super::super::*;

    /// The heavy-target settings parse as THREADS:HITS:BUDGET_GB[:MIN_MB]; anything else is refused (None).
    #[test]
    fn heavy_config_parses() {
        assert_eq!(heavy::Config::parse("4:1e8:8"),
                   Some(heavy::Config { threads: 4, hits: 100_000_000, budget: 8_000_000_000, min_bytes: 256_000_000 }));
        assert_eq!(heavy::Config::parse("2:50000000:0.5:16"),
                   Some(heavy::Config { threads: 2, hits: 50_000_000, budget: 500_000_000, min_bytes: 16_000_000 }));
        for bad in ["", "4", "4:1e8", "0:1e8:8", "4:x:8", "4:1e8:0", "4:1e8:8:-1", "4:1e8:8:1:2"] {
            assert_eq!(heavy::Config::parse(bad), None, "{:?}", bad);
        }
    }

    /// The fixed setting of heavy targets apart, 4:1e8:8 (the non-test `heavy::config()` returns it without reading the
    /// environment).
    #[test]
    fn heavy_setting_is_fixed() {
        assert_eq!(heavy::Config::parse(heavy::SETTING),
                   Some(heavy::Config { threads: 4, hits: 100_000_000, budget: 8_000_000_000, min_bytes: 256_000_000 }));
    }

    /// The heavy budget: needs under the floor take nothing; a need waits while it does not fit and runs alone when
    /// larger than the budget; every take is given back.
    #[test]
    fn heavy_budget_waits_and_returns() {
        let b = heavy::Budget::new(100, 10);
        assert_eq!(b.take(5), 0);
        assert_eq!(b.take(60), 60);
        assert_eq!(b.take(40), 40);
        assert_eq!(b.used(), 100);
        let b = std::sync::Arc::new(b);
        let waiter = {
            let b = std::sync::Arc::clone(&b);
            std::thread::spawn(move || b.take(30))   // waits until 30 fits
        };
        std::thread::sleep(std::time::Duration::from_millis(50));
        assert_eq!(b.used(), 100, "the third take must wait");
        b.give(40);
        assert_eq!(waiter.join().unwrap(), 30);
        assert_eq!(b.used(), 90);
        b.give(60);
        b.give(30);
        assert_eq!(b.take(500), 500, "larger than the budget: alone, when nothing is held");
        b.give(500);
        assert_eq!(b.used(), 0);
    }

    /// The first search with heavy targets apart gives the rows of the plain search, whichever targets are heavy (none,
    /// some, all), on one or several heavy threads, with a budget that never fits (every need alone) or always does.
    #[test]
    fn heavy_targets_apart_keep_the_rows() {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(3).build().unwrap();
        let mut checked = 0;
        for (seed, noise, repeats) in [(1, 10, true), (2, 30, false), (3, 5, true)] {
            let (lengths, tokens) = if repeats { repeat_reads(noise, seed) } else { noisy_reads(noise, seed) };
            let ids: Vec<String> = (0..tokens.len()).map(|r| format!("read{:02}", r)).collect();
            let p = Params { max_reads: 0, min_anchors: 3, min_overlap: 500, min_cover: 0.1, coverage: 2, ..Params::default() };
            heavy::TEST.with(|c| c.set(None));
            let plain = partner_rows(&ids, &lengths, &tokens, None, &p, &[], &pool);
            let plain_packed = partner_rows_packed(&ids, &lengths, &tokens, None, &p, &[], &pool, None, None).unpack();
            assert_eq!(plain_packed, plain);
            for hits in [0u64, 40, 1 << 40] {
                for threads in [1usize, 4] {
                    for (budget, min_bytes) in [(1usize, 0usize), (1 << 40, 0), (1 << 30, 1 << 20)] {
                        heavy::TEST.with(|c| c.set(Some(heavy::Config { threads, hits, budget, min_bytes })));
                        let got = partner_rows(&ids, &lengths, &tokens, None, &p, &[], &pool);
                        let got_packed = partner_rows_packed(&ids, &lengths, &tokens, None, &p, &[], &pool, None, None).unpack();
                        heavy::TEST.with(|c| c.set(None));
                        assert_eq!(got, plain, "seed {} hits {} threads {} budget {}", seed, hits, threads, budget);
                        assert_eq!(got_packed, plain, "packed: seed {} hits {} threads {} budget {}", seed, hits, threads, budget);
                        checked += 1;
                    }
                }
            }
        }
        assert_eq!(checked, 54);
    }
}
