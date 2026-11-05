//! Every tunable value of a run, with its default under each preset.
//!
//! The table below is the single source of the options of `squash`: their names, help, units and defaults (the
//! values of the SQUASH paper: `small` for genomes under 10 Mb, `large` from 10 Mb). Values are kept as strings,
//! checked against their kind, and passed to the steps as written (only bp suffixes are expanded). Values that are
//! not options (fixed parts of the method) live in the step code.

use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Preset {
    Small,
    Large,
}

impl Preset {
    pub fn name(self) -> &'static str {
        match self {
            Preset::Small => "small",
            Preset::Large => "large",
        }
    }
}

/// What a value must look like (checked before a run starts).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// A non-negative integer.
    Int,
    /// A number.
    Float,
    /// A length in bp; `k` / `m` suffixes allowed (`50k`), stored as plain bp.
    Bp,
    /// A non-negative integer, or `off`.
    IntOrOff,
    /// A fraction in [0, 1], or `off`.
    FractionOrOff,
    /// Comma-separated integers.
    IntList,
    /// Comma-separated numbers.
    FloatList,
}

pub struct Param {
    pub name: &'static str,
    pub short: Option<char>,
    pub heading: &'static str,
    pub small: &'static str,
    pub large: &'static str,
    pub kind: Kind,
    pub help: &'static str,
}

pub const SIGNAL: &str = "Signal: events and q-symbols";
pub const CORRECTION: &str = "Correction in token space";
pub const ASSEMBLY: &str = "Assembly (multi-k minimizer-space de Bruijn graph)";
pub const THREADING: &str = "Read threading and contig joins";
pub const LAYOUT: &str = "Signal layout";

macro_rules! p {
    ($name:expr, $short:expr, $heading:expr, $small:expr, $large:expr, $kind:expr, $help:expr) => {
        Param { name: $name, short: $short, heading: $heading, small: $small, large: $large, kind: $kind, help: $help }
    };
}

pub const PARAMS: &[Param] = &[
    // Signal (RawHash2's R10.4.1 event detector, q-symbols, tokens).
    p!("evdet-win", None, SIGNAL, "3,6", "3,6", Kind::IntList, "Event detector window lengths, short,long (RawHash2 window_length1,2)"),
    p!("evdet-thr", None, SIGNAL, "5.5,3.5", "5.5,3.5", Kind::FloatList, "Event detector t-statistic thresholds, short,long (RawHash2 threshold1,2)"),
    p!("evdet-peak", None, SIGNAL, "0.2", "0.2", Kind::Float, "Event detector peak height (RawHash2 peak_height)"),
    p!("event-merge", None, SIGNAL, "0.30", "0.30", Kind::Float, "Merge consecutive events whose means stay within this many z-units (0 = off)"),
    p!("quant-bits", Some('q'), SIGNAL, "3", "3", Kind::Int, "Bits per q-symbol (2^q quantization levels of event means)"),
    p!("events", Some('e'), SIGNAL, "9", "9", Kind::Int, "q-symbols per token (each token hashes e consecutive q-symbols)"),
    p!("window", Some('w'), SIGNAL, "16", "16", Kind::Int, "Minimizer window: one token per window of w hashes"),
    p!("bases-per-second", None, SIGNAL, "400", "400", Kind::Float, "Translocation speed of the chemistry; converts signal length to bp for every bp limit (R10.4.1 400 bps)"),
    // Correction (nanoMDBG-style rounds in token space).
    p!("low-complexity", None, CORRECTION, "0.2", "0.2", Kind::FractionOrOff, "Leave out of correction (as targets and donors) a read whose most frequent token is at least this fraction of its tokens"),
    p!("min-anchors", None, CORRECTION, "4", "4", Kind::Int, "Shared tokens a partner's chain needs"),
    p!("min-overlap", None, CORRECTION, "500", "500", Kind::Int, "Overlap a partner needs, in q-symbols"),
    p!("min-cover", None, CORRECTION, "0.5", "0.5", Kind::Float, "Fraction of the implied overlap a partner's chain must cover"),
    p!("coverage", None, CORRECTION, "20", "20", Kind::Int, "Coverage cap: chains are taken best first while some target token in their span has fewer than N"),
    p!("seed-repeat-fraction", None, CORRECTION, "0.0001", "0.0001", Kind::Float, "Fraction of the most frequent token values not used as seeds (nanoMDBG's repeat rule)"),
    p!("rounds", None, CORRECTION, "12", "16", Kind::Int, "Correction rounds at most"),
    p!("stop-below", None, CORRECTION, "0.0005", "0.0005", Kind::Float, "Stop early when fewer than this fraction of tokens change in a round"),
    p!("band", None, CORRECTION, "256", "256", Kind::Int, "Band of the token-level POA alignment"),
    p!("hops-from", None, CORRECTION, "off", "2", Kind::IntOrOff, "From this round, look for partners two hops away (partners of partners)"),
    p!("changed-from", None, CORRECTION, "off", "2", Kind::IntOrOff, "From this round, search new partners only for reads that changed"),
    p!("skip-from", None, CORRECTION, "off", "8", Kind::IntOrOff, "From this round, correct only the segments that changed"),
    p!("skip-pad", None, CORRECTION, "150", "150", Kind::Int, "Q-symbols of padding around a changed segment (with --skip-from)"),
    p!("qsym-check", None, CORRECTION, "off", "0.667", Kind::FractionOrOff, "Q-symbol alignment check of partners: largest weighted edit cost per alignment column"),
    // Assembly (metaMDBG-style multi-k graph on 1 in `density` tokens).
    p!("density", None, ASSEMBLY, "8", "8", Kind::Int, "Assemble with 1 in N corrected tokens, chosen by value hash"),
    p!("repeat-ppm", None, ASSEMBLY, "100", "100", Kind::Int, "Leave out the most frequent token values, in parts per million (metaMDBG's repetitive-minimizer cut)"),
    p!("k-min", None, ASSEMBLY, "4", "4", Kind::Int, "Smallest k of the multi-k graph"),
    p!("k-max", None, ASSEMBLY, "32", "32", Kind::Int, "Largest k of the multi-k graph (the final graph's k)"),
    p!("abundance-cutoffs", None, ASSEMBLY, "2,3", "2,3", Kind::IntList, "Progressive abundance cutoffs of the graph cleaning"),
    p!("tip-max", None, ASSEMBLY, "50000", "50000", Kind::Bp, "Longest tip removed, in bp"),
    p!("superbubble-max", None, ASSEMBLY, "50000", "50000", Kind::Bp, "Longest superbubble popped, in bp"),
    p!("max-path-rows", None, ASSEMBLY, "100000000", "1000000000", Kind::Int, "Resource cap of the graph construction (no effect unless reached)"),
    // Threading and joins.
    p!("min-support", None, THREADING, "3", "3", Kind::Int, "Reads a threaded pairing of unitigs needs"),
    p!("max-conflict", None, THREADING, "0.2", "0.2", Kind::Float, "Largest fraction of conflicting reads of a threaded pairing"),
    p!("anchor-nodes", None, THREADING, "200", "200", Kind::Int, "Atomic nodes a unitig needs to anchor a tangle (200 nodes ~ 20 kb)"),
    p!("anchor-edge", None, THREADING, "40", "40", Kind::Int, "Nodes at an anchor's end where reads are counted"),
    p!("tangle-dominance", None, THREADING, "3", "3", Kind::Int, "A tangle's pairing needs this many times the reads of the next best"),
    p!("join-paths", None, THREADING, "2", "2", Kind::Int, "Independent read walks a contig join needs"),
    p!("join-threshold", None, THREADING, "0.62", "0.62", Kind::Float, "Largest weighted edit cost per column of a join hop's q-symbol alignment"),
    p!("join-min-nodes", None, THREADING, "200", "200", Kind::Int, "Nodes a unitig needs for its ends to be joined, or to be landed on"),
    p!("join-head", None, THREADING, "200", "200", Kind::Int, "A join lands only in the first N tokens of the other unitig"),
    // Signal layout.
    p!("layout-max-gap", None, LAYOUT, "32", "32", Kind::Int, "Unmatched windows allowed inside a read's piece of a unitig"),
    p!("layout-own-gap", None, LAYOUT, "2000", "2000", Kind::Int, "Largest gap, in q-symbols, between a piece's tokens read from its own signal"),
];


/// Checks a value against its kind and returns it as the step will see it (bp suffixes expanded; else unchanged).
pub fn normalize(p: &Param, value: &str) -> Result<String, String> {
    let bad = || format!("--{}: '{}' is not {}", p.name, value, describe(p.kind));
    let v = value.trim();
    let int_ok = |s: &str| !s.is_empty() && s.chars().all(|c| c.is_ascii_digit());
    let float_ok = |s: &str| s.parse::<f64>().map(|x| x.is_finite()).unwrap_or(false);
    match p.kind {
        Kind::Int => int_ok(v).then(|| v.to_string()).ok_or_else(bad),
        Kind::Float => float_ok(v).then(|| v.to_string()).ok_or_else(bad),
        Kind::IntOrOff => (v == "off" || int_ok(v)).then(|| v.to_string()).ok_or_else(bad),
        Kind::FractionOrOff => {
            if v == "off" {
                return Ok(v.to_string());
            }
            match v.parse::<f64>() {
                Ok(x) if (0.0..=1.0).contains(&x) => Ok(v.to_string()),
                _ => Err(bad()),
            }
        }
        Kind::IntList => v.split(',').all(int_ok).then(|| v.to_string()).ok_or_else(bad),
        Kind::FloatList => v.split(',').all(float_ok).then(|| v.to_string()).ok_or_else(bad),
        Kind::Bp => parse_bp(v).map(|n| n.to_string()).ok_or_else(bad),
    }
}

pub fn describe(kind: Kind) -> &'static str {
    match kind {
        Kind::Int => "a non-negative integer",
        Kind::Float => "a number",
        Kind::Bp => "a length in bp (e.g. 2500, 50k, 1.5m)",
        Kind::IntOrOff => "a non-negative integer or 'off'",
        Kind::FractionOrOff => "a fraction in [0, 1] or 'off'",
        Kind::IntList => "comma-separated integers",
        Kind::FloatList => "comma-separated numbers",
    }
}

/// `2500`, `50k`, `1.5m`, `2g` → bp.
pub fn parse_bp(s: &str) -> Option<u64> {
    let s = s.trim().to_ascii_lowercase();
    let (num, mult) = match s.chars().last()? {
        'k' => (&s[..s.len() - 1], 1e3),
        'm' => (&s[..s.len() - 1], 1e6),
        'g' => (&s[..s.len() - 1], 1e9),
        _ => (&s[..], 1.0),
    };
    if mult == 1.0 {
        return num.parse::<u64>().ok();
    }
    let x = num.parse::<f64>().ok()?;
    if !(x.is_finite() && x >= 0.0) {
        return None;
    }
    Some((x * mult).round() as u64)
}

/// The resolved values of a run: each parameter's value and whether the user set it.
#[derive(Clone, Debug)]
pub struct Values {
    pub preset: Preset,
    values: BTreeMap<&'static str, (String, bool)>,
}

impl Values {
    /// The preset's values, then the user's (already normalized) values on top.
    pub fn resolve(preset: Preset, user: &[(&'static str, String)]) -> Values {
        let mut values = BTreeMap::new();
        for p in PARAMS {
            let d = match preset {
                Preset::Small => p.small,
                Preset::Large => p.large,
            };
            values.insert(p.name, (d.to_string(), false));
        }
        for (k, v) in user {
            values.insert(k, (v.clone(), true));
        }
        Values { preset, values }
    }

    pub fn get(&self, name: &str) -> &str {
        &self.values.get(name).unwrap_or_else(|| panic!("no parameter {}", name)).0
    }

    #[cfg(test)]
    /// `None` for `off`, else the integer.
    pub fn int_or_off(&self, name: &str) -> Option<u64> {
        match self.get(name) {
            "off" => None,
            v => Some(v.parse().expect("checked integer")),
        }
    }

    #[cfg(test)]
    pub fn is_off(&self, name: &str) -> bool {
        self.get(name) == "off"
    }

    pub fn iter(&self) -> impl Iterator<Item = (&'static str, &str, bool)> + '_ {
        self.values.iter().map(|(k, (v, user))| (*k, v.as_str(), *user))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_are_valid_and_differ_only_where_the_presets_do() {
        let differ: Vec<&str> = PARAMS.iter().filter(|p| p.small != p.large).map(|p| p.name).collect();
        assert_eq!(differ, ["rounds", "hops-from", "changed-from", "skip-from", "qsym-check", "max-path-rows"]);
        for p in PARAMS {
            for v in [p.small, p.large] {
                assert_eq!(normalize(p, v).as_deref(), Ok(v), "{}", p.name);
            }
        }
        let names: std::collections::HashSet<&str> = PARAMS.iter().map(|p| p.name).collect();
        assert_eq!(names.len(), PARAMS.len(), "duplicate parameter names");
    }

    #[test]
    fn bp_suffixes() {
        assert_eq!(parse_bp("2500"), Some(2500));
        assert_eq!(parse_bp("50k"), Some(50_000));
        assert_eq!(parse_bp("1.5m"), Some(1_500_000));
        assert_eq!(parse_bp("2G"), Some(2_000_000_000));
        assert_eq!(parse_bp("x"), None);
        assert_eq!(parse_bp("-1k"), None);
    }

    #[test]
    fn user_values_override_the_preset() {
        let v = Values::resolve(Preset::Large, &[("rounds", "20".to_string())]);
        assert_eq!(v.get("rounds"), "20");
        assert_eq!(v.get("qsym-check"), "0.667");
        assert_eq!(v.int_or_off("hops-from"), Some(2));
        let s = Values::resolve(Preset::Small, &[]);
        assert_eq!(s.int_or_off("hops-from"), None);
        assert!(s.is_off("qsym-check"));
    }
}

/// The options table of docs/USAGE.md (`squash __params-md` prints it).
pub fn markdown() -> String {
    let mut s = String::new();
    let mut heading = "";
    for p in PARAMS {
        if p.heading != heading {
            heading = p.heading;
            s += &format!("\n### {}\n\n| Option | small | large | What it sets |\n|---|---|---|---|\n", heading);
        }
        let name = match p.short {
            Some(c) => format!("`-{}`, `--{}`", c, p.name),
            None => format!("`--{}`", p.name),
        };
        s += &format!("| {} | {} | {} | {} |\n", name, p.small, p.large, p.help);
    }
    s
}
