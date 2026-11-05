//! The command line of `squash`: one required choice (`-x small|large`), the input and output, and every tunable value
//! of `presets.rs` as an option whose default is the chosen preset's value.

use crate::pipeline::{self, Config, STAGES};
use crate::presets::{self, Kind, Preset, Values, PARAMS};
use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use std::path::PathBuf;

const ABOUT: &str = "SQUASH: assembling raw nanopore signals with minimizer-space de Bruijn graphs";

const LONG_ABOUT: &str = "\
SQUASH assembles nanopore reads from their raw current signal (POD5), without basecalling. Each read's signal becomes
events and 3-bit q-symbols, then minimizer tokens; the tokens are corrected in rounds against overlapping reads, as
nanoMDBG corrects minimizers, and assembled into a multi-k directed minimizer-space de Bruijn graph, as metaMDBG does;
read threading and contig joins finish the graph. The output is the final graph at k = 32 (token unitigs, no bases),
each read's path through it, and each unitig's signal layout.

Choose the preset by genome size: -x small for genomes under 10 Mb (bacteria, archaea), -x large from 10 Mb
(eukaryotes, single chromosomes). Every other value defaults to the paper's for that preset; `squash --help` lists
them all and `--dry-run` prints the resolved values. Both strands assemble on their own (raw signal has no reverse
complement), so give at least 20x per strand (40x of the genome); the paper uses 30x and 50x per strand.";

const AFTER_HELP: &str = "Examples:
  squash -x small -o ecoli ecoli.pod5
  squash -x large -t 64 --tmp-dir /scratch/$USER -o dmel dmel_pod5/
  squash -x large --rounds 20 -o dmel_r20 dmel_pod5/
Documentation: docs/USAGE.md (every option), docs/OUTPUTS.md (the files), https://github.com/STORMgroup/SQUASH";

pub fn command() -> Command {
    let mut cmd = Command::new("squash")
        .version(crate::version())
        .about(ABOUT)
        .long_about(LONG_ABOUT)
        .after_help(AFTER_HELP)
        .next_line_help(false)
        .arg(
            Arg::new("preset")
                .short('x')
                .long("preset")
                .value_name("small|large")
                .value_parser(["small", "large"])
                .required(true)
                .help("Genome size class: small (< 10 Mb) or large (>= 10 Mb); sets every default below"),
        )
        .arg(
            Arg::new("out-dir")
                .short('o')
                .long("out-dir")
                .value_name("DIR")
                .value_parser(value_parser!(PathBuf))
                .required(true)
                .help("Output folder (created if missing)"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .value_name("N")
                .value_parser(value_parser!(usize))
                .help("Threads [default: every core this process may use]; outputs do not depend on it"),
        )
        .arg(
            Arg::new("reads")
                .value_name("READS")
                .value_parser(value_parser!(PathBuf))
                .num_args(1..)
                .required(true)
                .help("POD5 files, or folders of them (a folder's *.pod5 in name order)"),
        )
        .arg(
            Arg::new("genome-size")
                .short('g')
                .long("genome-size")
                .value_name("SIZE")
                .help("Expected genome size (e.g. 4.6m), only for the coverage estimate and its warnings"),
        )
        .arg(
            Arg::new("read-ids")
                .long("read-ids")
                .value_name("FILE")
                .value_parser(value_parser!(PathBuf))
                .help("Assemble only these reads (one id per line, plain or gzipped); the read order is fixed over every read first"),
        )
        .arg(
            Arg::new("tmp-dir")
                .long("tmp-dir")
                .value_name("DIR")
                .value_parser(value_parser!(PathBuf))
                .help("Scratch folder for large temporary files [default: OUT/work/tmp]; a local disk helps"),
        )
        .arg(
            Arg::new("keep-work")
                .long("keep-work")
                .action(ArgAction::SetTrue)
                .help("Keep the intermediate files in OUT/work (deleted after a successful run otherwise)"),
        )
        .arg(
            Arg::new("resume")
                .long("resume")
                .action(ArgAction::SetTrue)
                .help("Continue a stopped run in OUT: finished stages run with the same values are kept"),
        )
        .arg(
            Arg::new("stop-after")
                .long("stop-after")
                .value_name("STAGE")
                .value_parser(STAGES.to_vec())
                .help("Stop after this stage (implies --keep-work)"),
        )
        .arg(Arg::new("dry-run").long("dry-run").action(ArgAction::SetTrue).help("Print every resolved value and stop"))
        .arg(
            Arg::new("no-threading")
                .long("no-threading")
                .action(ArgAction::SetTrue)
                .help_heading(presets::THREADING)
                .hide_short_help(true)
                .help("Output the assembly graph before read threading and joins"),
        )
        .arg(
            Arg::new("no-tangles")
                .long("no-tangles")
                .action(ArgAction::SetTrue)
                .help_heading(presets::THREADING)
                .hide_short_help(true)
                .help("Thread reads through repeats but not through tangles"),
        )
        .arg(
            Arg::new("no-joins")
                .long("no-joins")
                .action(ArgAction::SetTrue)
                .help_heading(presets::THREADING)
                .hide_short_help(true)
                .help("No contig joins"),
        );
    for p in PARAMS {
        let default = if p.small == p.large {
            format!("[default: {}]", p.small)
        } else {
            format!("[small: {}, large: {}]", p.small, p.large)
        };
        let mut arg = Arg::new(p.name)
            .long(p.name)
            .value_name(value_name(p.kind))
            .help(format!("{} {}", p.help, default))
            .help_heading(p.heading)
            .hide_short_help(true);
        if let Some(c) = p.short {
            arg = arg.short(c);
        }
        cmd = cmd.arg(arg);
    }
    cmd
}

fn value_name(kind: Kind) -> &'static str {
    match kind {
        Kind::Int => "N",
        Kind::Float => "X",
        Kind::Bp => "BP",
        Kind::IntOrOff => "N|off",
        Kind::FractionOrOff => "F|off",
        Kind::IntList => "N,N",
        Kind::FloatList => "X,X",
    }
}

/// Parses the command line, then runs the pipeline (or prints the values). Returns the exit code.
pub fn main(args: Vec<String>) -> i32 {
    let m = match command().try_get_matches_from(args) {
        Ok(m) => m,
        Err(e) => {
            let _ = e.print();
            return if e.use_stderr() { 2 } else { 0 };
        }
    };
    match config(&m) {
        Ok(cfg) => {
            if m.get_flag("dry-run") {
                print!("{}", pipeline::params_tsv(&cfg));
                return 0;
            }
            match pipeline::run(&cfg) {
                Ok(()) => 0,
                Err(e) => {
                    eprintln!("squash: {}", e);
                    1
                }
            }
        }
        Err(e) => {
            eprintln!("squash: {}", e);
            2
        }
    }
}

fn config(m: &ArgMatches) -> Result<Config, String> {
    let preset = match m.get_one::<String>("preset").map(String::as_str) {
        Some("small") => Preset::Small,
        Some("large") => Preset::Large,
        _ => unreachable!("clap checks the preset"),
    };
    let mut user = Vec::new();
    for p in PARAMS {
        if let Some(v) = m.get_one::<String>(p.name) {
            user.push((p.name, presets::normalize(p, v)?));
        }
    }
    let values = Values::resolve(preset, &user);
    let genome_size = match m.get_one::<String>("genome-size") {
        Some(g) => Some(presets::parse_bp(g).ok_or_else(|| format!("--genome-size: '{}' is not a size (e.g. 4.6m)", g))?),
        None => None,
    };
    let threads = m
        .get_one::<usize>("threads")
        .copied()
        .unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1))
        .max(1);
    let stop_after = m.get_one::<String>("stop-after").cloned();
    Ok(Config {
        reads: m.get_many::<PathBuf>("reads").unwrap().cloned().collect(),
        out_dir: m.get_one::<PathBuf>("out-dir").unwrap().clone(),
        threads,
        tmp_dir: m.get_one::<PathBuf>("tmp-dir").cloned(),
        read_ids: m.get_one::<PathBuf>("read-ids").cloned(),
        genome_size,
        keep_work: m.get_flag("keep-work") || stop_after.is_some(),
        resume: m.get_flag("resume"),
        stop_after,
        threading: !m.get_flag("no-threading"),
        tangles: !m.get_flag("no-tangles"),
        joins: !m.get_flag("no-joins"),
        values,
        command_line: std::env::args().collect::<Vec<_>>().join(" "),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(args: &[&str]) -> Result<Config, String> {
        let m = command().try_get_matches_from(args.iter().map(|s| s.to_string())).map_err(|e| e.to_string())?;
        config(&m)
    }

    #[test]
    fn the_preset_is_required() {
        assert!(parse(&["squash", "-o", "out", "r.pod5"]).is_err());
        let c = parse(&["squash", "-x", "large", "-o", "out", "r.pod5", "--rounds", "20", "--tip-max", "30k"]).unwrap();
        assert_eq!(c.values.get("rounds"), "20");
        assert_eq!(c.values.get("tip-max"), "30000");
        assert_eq!(c.values.get("qsym-check"), "0.667");
    }

    #[test]
    fn values_are_checked() {
        assert!(parse(&["squash", "-x", "small", "-o", "o", "r.pod5", "--rounds", "x"]).is_err());
        assert!(parse(&["squash", "-x", "small", "-o", "o", "r.pod5", "--qsym-check", "2"]).is_err());
        assert!(parse(&["squash", "-x", "small", "-o", "o", "r.pod5", "--hops-from", "off"]).is_ok());
    }

    #[test]
    fn help_builds() {
        command().debug_assert();
    }
}
