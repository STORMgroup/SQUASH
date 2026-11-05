//! The partner search's options (`Params`, `Guard`) and their parsing.

use super::*;

/// The search's options (`take_params`, the `--max-reads` ... options).
#[derive(Clone, Copy, Debug)]
pub struct Params {
    /// Seeds are the values held by 2..=max_reads distinct reads (0 = no cap).
    pub max_reads: usize,
    /// The fewest anchors of a kept chain.
    pub min_anchors: usize,
    /// The largest step between chained anchors on either read (tokens).
    pub max_gap: i64,
    /// The largest drift |dd - dt| between chained anchors: max_drift + drift_fraction x max(dt, dd) (tokens).
    pub max_drift: i64,
    pub drift_fraction: f64,
    /// How many anchors before an anchor the chaining DP looks back.
    pub lookback: usize,
    /// The shortest implied overlap of a kept chain (q-symbols).
    pub min_overlap: i64,
    /// The least fraction of the implied overlap the chain's span on the target must cover.
    pub min_cover: f64,
    /// The coverage cap: a chain is kept while some target token of its span has fewer than this many kept chains.
    pub coverage: usize,
    /// The fraction of distinct values, held by the most reads, that are not seeds (0 = none; nanoMDBG's repeat rule).
    pub repeat_fraction: f64,
    /// With the signal check (-1 = off): the chain cover test of this search uses this instead of `min_cover`
    /// while a signal check is installed (`signal`): many more candidates, which the signal check then filters. Without
    /// a signal check it has no effect.
    pub signal_min_cover: f64,
}

impl Default for Params {
    fn default() -> Self {
        Params { max_reads: 60, min_anchors: 3, max_gap: 120, max_drift: 6, drift_fraction: 0.25, lookback: 60,
                 min_overlap: 1000, min_cover: 0.5, coverage: 20, repeat_fraction: 0.0,
                 signal_min_cover: -1.0 }
    }
}


pub(super) fn get<T: std::str::FromStr>(map: &mut HashMap<String, String>, key: &str, default: T) -> Result<T>
where
    T::Err: std::fmt::Display,
{
    match map.remove(key) {
        None => Ok(default),
        Some(v) => v.parse::<T>().map_err(|e| format!("--{} {:?}: {}", key, v, e).into()),
    }
}

/// The partner search's options (`--max-reads` ... `--repeat-fraction`, `--signal-min-cover`), removed from MAP.
pub fn take_params(map: &mut HashMap<String, String>) -> Result<Params> {
    let d = Params::default();
    let p = Params {
        max_reads: get(map, "max-reads", d.max_reads)?,
        min_anchors: get(map, "min-anchors", d.min_anchors)?,
        max_gap: get(map, "max-gap", d.max_gap)?,
        max_drift: get(map, "max-drift", d.max_drift)?,
        drift_fraction: get(map, "drift-fraction", d.drift_fraction)?,
        lookback: get(map, "lookback", d.lookback)?,
        min_overlap: get(map, "min-overlap", d.min_overlap)?,
        min_cover: get(map, "min-cover", d.min_cover)?,
        coverage: get(map, "coverage", d.coverage)?,
        repeat_fraction: get(map, "repeat-fraction", d.repeat_fraction)?,
        signal_min_cover: get(map, "signal-min-cover", d.signal_min_cover)?,
        ..d
    };
    if !(p.signal_min_cover == -1.0 || (0.0..=1.0).contains(&p.signal_min_cover)) {
        return Err("--signal-min-cover must be in [0, 1] (or -1 = off)".into());
    }
    Ok(p)
}

/// The raw-evidence guard's options (`--guard-min-anchors`, `--guard-band`, `--guard-band-fraction`), removed from
/// MAP.
pub fn take_guard(map: &mut HashMap<String, String>) -> Result<Guard> {
    let g = Guard::default();
    Ok(Guard {
        min_anchors: get(map, "guard-min-anchors", g.min_anchors)?,
        band: get(map, "guard-band", g.band)?,
        band_fraction: get(map, "guard-band-fraction", g.band_fraction)?,
    })
}
