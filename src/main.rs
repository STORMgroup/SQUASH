//! SQUASH: assembling raw nanopore signals with minimizer-space de Bruijn graphs.
//!
//! One binary. `squash -x small|large -o DIR READS.pod5` runs the whole pipeline (`pipeline.rs`). Each step of the
//! pipeline keeps its own option parser and runs in a child process as `squash __STEP ARGS...` (`run_step`), with the
//! fixed arguments the pipeline builds for it.
#![allow(clippy::all)]

mod argv;
mod assemble;
mod bp_scale;
mod cli;
mod correct;
mod events;
mod exec;
mod graph;
mod io;
mod joins;
mod layout;
mod partners;
mod pipeline;
mod prepare;
mod presets;
mod rounds;
mod signal;
mod sparse;
mod threading;

/// The steps a child process can run, each with a note on what it is (a test checks that every step the drivers run
/// is here and dispatched).
pub const STEPS: &[(&str, &str)] = &[
    ("signal", "POD5 -> events -> q-symbols, and the events' POD5 sample spans (signal.rs)"),
    ("low-complexity", "the reads left out of correction (prepare.rs)"),
    ("partners", "the partner search (partners.rs)"),
    ("rounds", "the correction rounds (rounds.rs)"),
    ("sparse", "the sparse tokens and the repetitive-value cut (sparse.rs)"),
    ("bp-scale", "signal lengths and bp per token (bp_scale.rs)"),
    ("graph", "the multi-k graph (graph/)"),
    ("thread", "read threading and contig joins (threading.rs, joins.rs)"),
    ("layout", "the signal layout of the final graph (layout.rs)"),
];

/// Runs one step in this process with `args` (`[STEP, ARGS...]`) and returns its exit code: 0 on success, 1 on an
/// error (printed on standard error), 2 when the graph step stops at a resource cap (`RESOURCE_CAP`), 64 for an
/// unknown step.
fn run_step(step: &str, args: Vec<String>) -> i32 {
    argv::set(args.clone());
    let report = |name: &str, e: &dyn std::fmt::Display| {
        eprintln!("{}: {}", name, e);
        1
    };
    match step {
        "signal" => match signal::run() {
            Ok(()) => 0,
            Err(e) => report("signal", &e),
        },
        "low-complexity" => match prepare::run_low_complexity() {
            Ok(()) => 0,
            Err(e) => report("low-complexity", &e),
        },
        "partners" => match partners::run() {
            Ok(()) => 0,
            Err(e) => report("partners", &e),
        },
        "rounds" => match rounds::run() {
            Ok(()) => 0,
            Err(e) => report("rounds", &e),
        },
        "sparse" => match sparse::run() {
            Ok(()) => 0,
            Err(e) => report("sparse", &e),
        },
        "bp-scale" => match bp_scale::run() {
            Ok(()) => 0,
            Err(e) => report("bp_scale", &e),
        },
        "graph" => match graph::run() {
            Ok(()) => 0,
            Err(error) => {
                if error.to_string().starts_with("RESOURCE_CAP:") {
                    eprintln!("stopped at a resource cap: {}", error.to_string().trim_start_matches("RESOURCE_CAP: "));
                    2
                } else {
                    eprintln!("integrity error: {}", error);
                    1
                }
            }
        },
        "params-md" => {
            print!("{}", presets::markdown());
            0
        }
        "layout" => match layout::run() {
            Ok(()) => 0,
            Err(e) => report("layout", &e),
        },
        "thread" => match threading::run() {
            Ok(()) => 0,
            Err(e) => report("threading", &e),
        },
        _ => {
            eprintln!("squash: unknown step '{}'", step);
            64
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if let Some(step) = args.get(1).and_then(|a| a.strip_prefix("__")) {
        let step = step.to_string();
        let mut step_args = vec![step.clone()];
        step_args.extend(args[2..].iter().cloned());
        std::process::exit(run_step(&step, step_args));
    }
    std::process::exit(cli::main(args));
}

/// `1.0.0 (commit abc1234)` when built from a git checkout, else the package version.
pub fn version() -> &'static str {
    match option_env!("SQUASH_GIT_COMMIT") {
        Some(c) if !c.is_empty() => Box::leak(format!("{} (commit {})", env!("CARGO_PKG_VERSION"), c).into_boxed_str()),
        _ => env!("CARGO_PKG_VERSION"),
    }
}

#[cfg(test)]
mod tests {
    use super::STEPS;

    /// The step names the stage drivers pass to `Runner::step` / `step_lines` (their second argument).
    fn called_steps() -> Vec<String> {
        let call = regex_lite_calls();
        let mut v = Vec::new();
        for src in [include_str!("pipeline.rs"), include_str!("prepare.rs"), include_str!("correct.rs"), include_str!("assemble.rs")] {
            v.extend(call(src));
        }
        v.sort();
        v.dedup();
        v
    }

    /// `.step(NAME, "STEP", ...)` / `.step_lines(NAME, "STEP", ...)`, NAME a literal or a `&format!(...)`.
    fn regex_lite_calls() -> impl Fn(&str) -> Vec<String> {
        |src: &str| {
            let mut out = Vec::new();
            for (i, _) in src.match_indices(".step") {
                let rest = &src[i + 5..];
                let rest = rest.strip_prefix("_lines").unwrap_or(rest);
                let Some(rest) = rest.strip_prefix('(') else { continue };
                // skip the first argument: up to the comma at depth 0
                let (mut depth, mut in_str, mut end) = (0i32, false, None);
                for (j, c) in rest.char_indices() {
                    match c {
                        '"' => in_str = !in_str,
                        '(' if !in_str => depth += 1,
                        ')' if !in_str => {
                            if depth == 0 {
                                break;
                            }
                            depth -= 1
                        }
                        ',' if !in_str && depth == 0 => {
                            end = Some(j);
                            break;
                        }
                        _ => {}
                    }
                }
                let Some(end) = end else { continue };
                let second = rest[end + 1..].trim_start();
                if let Some(s) = second.strip_prefix('"') {
                    out.push(s[..s.find('"').unwrap()].to_string());
                }
            }
            out
        }
    }

    #[test]
    fn every_step_the_drivers_run_is_dispatched() {
        let called = called_steps();
        assert!(called.len() >= 9, "{:?}", called);
        let main = include_str!("main.rs");
        for s in &called {
            assert!(STEPS.iter().any(|(n, _)| n == s), "step {} is not in STEPS", s);
            assert!(main.contains(&format!("\"{}\" =>", s)), "step {} is not dispatched by run_step", s);
        }
    }
}
