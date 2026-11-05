//! The assemble, thread and layout stages, from the final correction round to the output files: the sparse table's bp
//! scale; the multi-k directed minimizer-space de Bruijn graph (k = 4..32 on the sparse tokens, metaMDBG's graph
//! rules; `graph/`); read threading and contig joins (`threading.rs`); the signal layout (`layout.rs`); and the final
//! GFA with its tables. Each step runs as a `squash __STEP` child through `exec::Runner`; pruning and the output files
//! are done in process here.
//!
//! Files (ROOT = the work folder, ROUND = ROOT/correction/roundN, the final correction round):
//!
//! | step (its log name) | reads | writes |
//! |---|---|---|
//! | `bp-scale-lengths` (signal stage) | POD5 (reads table only) | ROOT/bp_scale/signal_lengths.tsv |
//! | `bp-scale-counts` | ROUND/sparse/tokens.tsv.gz | ROOT/bp_scale/token_counts.tsv.gz |
//! | `bp-scale` | signal_lengths.tsv, token_counts.tsv.gz | ROOT/bp_scale/bp_per_token.tsv |
//! | `graph-input` | ROUND/sparse/tokens.tsv.gz | ROOT/graph/input/ (reads.tsv.gz, tokens.u32le, token_spans.u64le) |
//! | `graph` | ROOT/graph/input/ | ROOT/graph/kNN.{node_tokens.u32le,unitigs.tsv.gz,unitigs.gfa} |
//! | prune (in process) | | removes ROOT/graph/kNN.* for NN != 32 |
//! | `thread` (thread stage) | ROOT/graph/k32.*, ROOT/graph/input/, ROOT/signal/qsymbols.tsv, ROOT/correction/partners.tsv.gz | ROOT/threaded/ |
//! | `layout` | ROOT/threaded/k32.*, ROOT/graph/input/, ROUND/tokens.tsv.gz, ROOT/tokens/, ROOT/signal/spans.bin | ROOT/layout/ |
//! | output files (in process) | ROOT/threaded/k32.unitigs.gfa, k32.unitigs.tsv.gz, unitig_map.tsv, ROOT/layout/ | ROOT/assembly_graph.gfa, assembly_info.tsv, read_paths.tsv.gz, layout.tsv.gz |
//!
//! With `--no-threading` the final graph is ROOT/graph itself. The assemble stage ends after the prune, the thread
//! stage after threading.

use crate::exec::Runner;
use crate::graph::{GFA, NODE_TOKENS, UNITIGS, order_file};
use crate::pipeline::work;
use crate::bp_scale::{bp_per_token, bp_to_tokens};
use flate2::write::GzEncoder;
use flate2::Compression;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

type Result<T> = std::result::Result<T, String>;

fn err<E: std::fmt::Display>(what: &str) -> impl Fn(E) -> String + '_ {
    move |e| format!("{}: {}", what, e)
}

/// The values of the assemble, thread and layout stages; `small()` / `large()` hold the presets' values.
#[derive(Clone, Debug)]
pub struct AssembleParams {
    /// signal bases per second (R10.4.1 chemistry)
    pub bases_per_second: String,
    /// the graph's first and last orders
    pub k_min: u16,
    pub k_max: u16,
    /// the progressive abundance cutoffs (`2,3`: metaMDBG's, the graph step's own)
    pub abundance_cutoffs: String,
    /// metaMDBG's limits in bp, converted with the sparse table's bp-per-token scale
    pub tip_max_bp: i64,
    pub superbubble_max_bp: i64,
    /// the graph's path-row cap (a resource cap: a stop when exceeded, no other effect); None = the graph step's own
    /// 1e8
    pub max_path_rows: Option<u64>,
    /// read threading and two-walk joins
    pub min_support: u32,
    pub max_conflict: f64,
    pub anchor_nodes: u32,
    pub anchor_edge: u32,
    pub tangle_dominance: f64,
    pub join_threshold: f64,
    pub join_min_nodes: u32,
    pub join_head: u32,
    pub join_min_paths: u32,
    /// signal layout
    pub layout_max_gap: u32,
    pub own_gap: u32,
    pub own_min: u32,
    /// tangle threading with closed single tangles (`--tangles`) and the joins (`--join`)
    pub tangles: bool,
    pub joins: bool,
}

impl AssembleParams {
    /// Preset `small` (< 10 Mb), with the graph step's own path-row cap.
    pub fn small() -> AssembleParams {
        AssembleParams {
            bases_per_second: "400".to_string(),
            k_min: 4,
            k_max: 32,
            abundance_cutoffs: "2,3".to_string(),
            tip_max_bp: 50_000,
            superbubble_max_bp: 50_000,
            max_path_rows: None,
            min_support: 3,
            max_conflict: 0.2,
            anchor_nodes: 200,
            anchor_edge: 40,
            tangle_dominance: 3.0,
            join_threshold: 0.62,
            join_min_nodes: 200,
            join_head: 200,
            join_min_paths: 2,
            layout_max_gap: 32,
            own_gap: 2_000,
            own_min: 2,
            tangles: true,
            joins: true,
        }
    }

    /// Preset `large` (>= 10 Mb): as small, with a path-row cap of 1e9.
    pub fn large() -> AssembleParams {
        AssembleParams { max_path_rows: Some(1_000_000_000), ..AssembleParams::small() }
    }
}

/// The token limits of the assembly at the sparse table's bp-per-token scale.
#[derive(Clone, Debug, PartialEq)]
pub struct TokenLimits {
    pub bp_per_token: f64,
    pub tip_max_tokens: i64,
    pub superbubble_max_tokens: i64,
}

/// Every limit as tokens: max(1, round_half_even(bp / scale)), with the scale as read back from its table (written
/// with 6 decimals).
pub fn token_limits(scale: f64, p: &AssembleParams) -> Result<TokenLimits> {
    let t = |bp: i64| bp_to_tokens(bp, scale).map_err(|e| e.to_string());
    Ok(TokenLimits {
        bp_per_token: scale,
        tip_max_tokens: t(p.tip_max_bp)?,
        superbubble_max_tokens: t(p.superbubble_max_bp)?,
    })
}

/// The run's folders and the inputs of these steps (absolute). `new` sets every input to its place in the work
/// folder; a caller may point a field elsewhere.
#[derive(Clone, Debug)]
pub struct RunPaths {
    /// the work folder
    pub root: PathBuf,
    /// the uncorrected tokens (ROOT/tokens): the layout
    pub uncorrected: PathBuf,
    /// round N's corrected tokens (ROUND/tokens.tsv.gz): the layout
    pub corrected_tokens: PathBuf,
    /// round N's sparse table (ROUND/sparse/tokens.tsv.gz): the graph's input
    pub sparse_tokens: PathBuf,
    /// the uncorrected q-symbols (ROOT/signal/qsymbols.tsv): the joins
    pub symbols: PathBuf,
    /// the first partner search's table (ROOT/correction/partners.tsv.gz): the joins
    pub partners: PathBuf,
}

impl RunPaths {
    pub fn new(root: &Path, round: &Path) -> Result<RunPaths> {
        let root = fs::canonicalize(root).map_err(err(&root.display().to_string()))?;
        let round = fs::canonicalize(round).map_err(err(&round.display().to_string()))?;
        Ok(RunPaths {
            uncorrected: root.join(work::TOKENS),
            corrected_tokens: round.join("tokens.tsv.gz"),
            sparse_tokens: round.join("sparse/tokens.tsv.gz"),
            symbols: root.join(work::QSYMBOLS),
            partners: root.join(work::PARTNERS),
            root,
        })
    }
    pub fn lengths(&self) -> PathBuf {
        self.root.join(work::SIGNAL_LENGTHS)
    }
    /// The bp scale's folder (the signal lengths, and the sparse table's token counts and bp per token).
    pub fn scale_dir(&self) -> PathBuf {
        self.root.join("bp_scale")
    }
    /// The multi-k graph: its input reads (`input/`) and its orders' files.
    pub fn graph_dir(&self) -> PathBuf {
        self.root.join("graph")
    }
    /// The threaded graph: the K-max order's files, unitig_map.tsv, joins.tsv, stats.tsv.
    pub fn threaded_dir(&self) -> PathBuf {
        self.root.join("threaded")
    }
}

fn s(p: &Path) -> String {
    p.display().to_string()
}

// ---------------------------------------------------------------- bp scale

/// The bp scale of the sparse table: its tokens per read, then its bp per token (the scale of every limit in bp);
/// returns the limits as tokens.
pub fn bp_scale_sparse(run: &Runner, paths: &RunPaths, p: &AssembleParams) -> Result<TokenLimits> {
    let dir = paths.scale_dir();
    fs::create_dir_all(&dir).map_err(err(&s(&dir)))?;
    let (counts, scale) = (dir.join("token_counts.tsv.gz"), dir.join("bp_per_token.tsv"));
    run.step("bp-scale-counts", "bp-scale",
             &["counts".to_string(), "--tokens".into(), s(&paths.sparse_tokens), "--output".into(),
               s(&counts)])?;
    run.step("bp-scale", "bp-scale",
             &["scale".to_string(), "--lengths".into(), s(&paths.lengths()), "--reads".into(), s(&counts),
               "--bases-per-second".into(), p.bases_per_second.clone(), "--output".into(), s(&scale)])?;
    let value = bp_per_token(&scale).map_err(|e| e.to_string())?;
    token_limits(value, p)
}

// ---------------------------------------------------------------- the graph

/// `assemble`'s arguments (the orders, then tips, superbubbles, other abundance cutoffs than 2,3, and the large
/// preset's path-row cap), every length in tokens.
pub fn assemble_args(graph_dir: &Path, threads: usize, p: &AssembleParams, t: &TokenLimits) -> Vec<String> {
    let mut a: Vec<String> = vec![
        "assemble".into(), "--reads-dir".into(), s(&graph_dir.join("input")), "--k-start".into(),
        p.k_min.to_string(), "--k-max".into(), p.k_max.to_string(), "--graph-dir".into(), s(graph_dir),
        "--threads".into(), threads.to_string(), "--tip-max-tokens".into(), t.tip_max_tokens.to_string(),
        "--superbubble-max-tokens".into(), t.superbubble_max_tokens.to_string(),
    ];
    if p.abundance_cutoffs != "2,3" {
        a.extend(["--progressive-cutoffs".to_string(), p.abundance_cutoffs.clone()]);
    }
    if let Some(cap) = p.max_path_rows {
        a.extend(["--max-path-rows".to_string(), cap.to_string()]);
    }
    a
}

/// The graph, orders K-min..K-max, on the sparse table: `adapt-tokens` (the table as the graph's input token arrays,
/// `input/`), then `assemble`. Returns the graph's folder.
pub fn assemble_graph(run: &Runner, paths: &RunPaths, threads: usize, p: &AssembleParams, t: &TokenLimits) -> Result<PathBuf> {
    let graph_dir = paths.graph_dir();
    if graph_dir.exists() {
        return Err(format!("refusing to overwrite {}", graph_dir.display()));
    }
    let input = graph_dir.join("input");
    fs::create_dir_all(&input).map_err(err(&s(&input)))?;
    run.step("graph-input", "graph",
             &["adapt-tokens".to_string(), "--tokens".into(), s(&paths.sparse_tokens), "--output-dir".into(), s(&input)])?;
    run.step("graph", "graph", &assemble_args(&graph_dir, threads, p, t))?;
    Ok(graph_dir)
}

// ---------------------------------------------------------------- prune

/// Removes the graph's files of every order but K (`kNN.*`); kK.* and `input/` stay. Returns the bytes freed.
pub fn prune_graph(graph_dir: &Path, keep_k: u16) -> Result<u64> {
    let order_dir = graph_dir;
    for keep in [order_file(order_dir, keep_k as usize, UNITIGS), order_file(order_dir, keep_k as usize, NODE_TOKENS)] {
        if !keep.is_file() {
            return Err(format!("refusing to prune {}: missing {}", graph_dir.display(), keep.display()));
        }
    }
    let mut files: Vec<PathBuf> = Vec::new();
    for entry in fs::read_dir(&order_dir).map_err(err(&s(&order_dir)))? {
        let entry = entry.map_err(err(&s(&order_dir)))?;
        let name = entry.file_name().to_string_lossy().to_string();
        let order = name.strip_prefix('k').and_then(|r| {
            let digits: String = r.chars().take_while(|c| c.is_ascii_digit()).collect();
            (!digits.is_empty() && r[digits.len()..].starts_with('.')).then(|| digits.parse::<u64>().ok()).flatten()
        });
        if entry.file_type().map_err(err(&name))?.is_file() && order.map_or(false, |k| k != keep_k as u64) {
            files.push(entry.path());
        }
    }
    files.sort();
    let mut freed = 0;
    for f in files {
        freed += f.metadata().map_err(err(&s(&f)))?.len();
        fs::remove_file(&f).map_err(err(&s(&f)))?;
    }
    Ok(freed)
}

// ---------------------------------------------------------------- threading and joins

/// The threading step's arguments: read threading, tangles with closed single tangles, and two-walk joins, less the
/// options `--no-tangles` / `--no-joins` leave out.
pub fn thread_args(graph_dir: &Path, out: &Path, symbols: &Path, partners: &Path, threads: usize, p: &AssembleParams) -> Vec<String> {
    let all: Vec<String> = vec![
        "--threads".into(), threads.to_string(), "--graph-dir".into(), s(graph_dir), "--k".into(), p.k_max.to_string(),
        "--min-support".into(), p.min_support.to_string(), "--max-conflict".into(), p.max_conflict.to_string(),
        "--tangles".into(), "--anchor-nodes".into(), p.anchor_nodes.to_string(), "--anchor-edge".into(),
        p.anchor_edge.to_string(), "--tangle-dominance".into(), p.tangle_dominance.to_string(),
        "--join".into(), "--join-symbols".into(), s(symbols), "--join-partners".into(), s(partners),
        "--join-threshold".into(), p.join_threshold.to_string(), "--join-min-nodes".into(), p.join_min_nodes.to_string(),
        "--join-head".into(), p.join_head.to_string(), "--join-min-paths".into(), p.join_min_paths.to_string(),
        "--output-dir".into(), s(out),
    ];
    if p.tangles && p.joins {
        return all;
    }
    // --no-tangles drops `--tangles` and `--tangle-dominance` with its value; --no-joins drops
    // `--join` and the options that only the joins read, with their values
    let tangle_opts = ["--tangles"];
    let tangle_vals = ["--tangle-dominance"];
    let join_opts = ["--join"];
    let join_vals = ["--join-symbols", "--join-partners", "--join-threshold", "--join-min-nodes", "--join-head", "--join-min-paths"];
    let mut out_args = Vec::new();
    let mut i = 0;
    while i < all.len() {
        let a = all[i].as_str();
        if (!p.tangles && tangle_opts.contains(&a)) || (!p.joins && join_opts.contains(&a)) {
            i += 1;
        } else if (!p.tangles && tangle_vals.contains(&a)) || (!p.joins && join_vals.contains(&a)) {
            i += 2;
        } else {
            out_args.push(all[i].clone());
            i += 1;
        }
    }
    out_args
}

/// Read threading, tangle threading (closed single tangles) and two-walk joins on the K-max graph. The joins read
/// the uncorrected q-symbols and the first partner search's table. Returns the threaded bundle.
pub fn thread(run: &Runner, paths: &RunPaths, graph_dir: &Path, threads: usize, p: &AssembleParams) -> Result<PathBuf> {
    let out = paths.threaded_dir();
    run.step("thread", "thread", &thread_args(graph_dir, &out, &paths.symbols, &paths.partners, threads, p))?;
    Ok(out)
}

// ---------------------------------------------------------------- layout and the final graph

/// The layout's arguments: the final graph, the graph's input reads, own tokens from the final round's full table
/// against the uncorrected tokens, POD5 sample columns from `spans` (the signal step's `--dump-event-spans`).
pub fn layout_args(threaded: &Path, graph_dir: &Path, corrected: &Path, uncorrected: &Path, spans: &Path, out: &Path,
                   threads: usize, p: &AssembleParams) -> Vec<String> {
    let mut a = vec!["--graph-dir".to_string(), s(threaded), "--reads-dir".into(), s(&graph_dir.join("input")),
                     "--spans".into(), s(spans)];
    a.extend([
        "--k".to_string(), p.k_max.to_string(), "--max-gap".into(), p.layout_max_gap.to_string(), "--own-gap".into(),
        p.own_gap.to_string(), "--own-min".into(), p.own_min.to_string(), "--corrected-tokens".into(), s(corrected),
        "--uncorrected-dir".into(), s(uncorrected), "--threads".into(), threads.to_string(), "--output-dir".into(), s(out),
    ]);
    a
}

/// The signal layout of the final graph into `out` (layout.tsv, pile.tsv.gz, pile_unverified.tsv.gz, stats.tsv).
pub fn layout(run: &Runner, paths: &RunPaths, graph_dir: &Path, threaded: &Path, spans: &Path, out: &Path, threads: usize,
              p: &AssembleParams) -> Result<PathBuf> {
    run.step("layout", "layout", &layout_args(threaded, graph_dir, &paths.corrected_tokens, &paths.uncorrected, spans, out,
                                             threads, p))?;
    Ok(out.to_path_buf())
}

/// Universal newlines: "\r\n" and "\r" become "\n".
fn universal_newlines(text: &str) -> String {
    if text.contains('\r') {
        text.replace("\r\n", "\n").replace('\r', "\n")
    } else {
        text.to_string()
    }
}

/// The ids of the order-k final unitigs flagged circular, in table order.
pub fn circular_unitigs(order_dir: &Path, k: u16) -> Result<Vec<String>> {
    let path = order_file(order_dir, k as usize, UNITIGS);
    let file = File::open(&path).map_err(err(&s(&path)))?;
    let mut text = String::new();
    BufReader::new(flate2::read::MultiGzDecoder::new(file)).read_to_string(&mut text).map_err(err(&s(&path)))?;
    let text = universal_newlines(&text);
    let mut lines = text.split_inclusive('\n');
    let head: Vec<&str> = lines.next().ok_or("empty P3 table")?.trim_end_matches('\n').split('\t').collect();
    let (iu, ic) = (head.iter().position(|h| *h == "unitig_id").ok_or("P3 table: no unitig_id")?,
                    head.iter().position(|h| *h == "circular").ok_or("P3 table: no circular")?);
    let mut out = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.strip_suffix('\n').unwrap_or(line).splitn(head.len(), '\t').collect();
        let (u, c) = (f.get(iu).ok_or("short P3 row")?, f.get(ic).ok_or("short P3 row")?);
        if *c == "1" {
            out.push(u.to_string());
        }
    }
    Ok(out)
}

/// The final graph SQUASH hands to users: the final graph's kK.unitigs.gfa followed by a closing self-link
/// `L u + u + 0M` for each circular unitig (the graph step's GFA has none).
pub fn write_final_graph(threaded: &Path, k: u16, out: &Path) -> Result<()> {
    let order_dir = threaded;
    let gfa = order_file(order_dir, k as usize, GFA);
    let text = universal_newlines(&fs::read_to_string(&gfa).map_err(err(&s(&gfa)))?);
    let mut body = text;
    for u in circular_unitigs(&order_dir, k)? {
        body.push_str(&format!("L\t{}\t+\t{}\t+\t0M\n", u, u));
    }
    let tmp = PathBuf::from(format!("{}.tmp", s(out)));
    fs::write(&tmp, body).map_err(err(&s(&tmp)))?;
    fs::rename(&tmp, out).map_err(err(&s(out)))?;
    Ok(())
}

// ---------------------------------------------------------------- all of it

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("squash_assemble_{}_{}", name, std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn presets_differ_only_in_the_path_row_cap() {
        let (a, b) = (AssembleParams::small(), AssembleParams::large());
        assert_eq!(a.max_path_rows, None);
        assert_eq!(b.max_path_rows, Some(1_000_000_000));
        assert_eq!(format!("{:?}", AssembleParams { max_path_rows: None, ..b }), format!("{:?}", a));
    }

    #[test]
    fn limits_at_the_paper_scales() {
        let p = AssembleParams::small();
        // the bp-per-token scales of two real runs
        let t = token_limits(89.879905, &p).unwrap();
        assert_eq!((t.tip_max_tokens, t.superbubble_max_tokens), (556, 556));
        let t = token_limits(94.654755, &p).unwrap();
        assert_eq!((t.tip_max_tokens, t.superbubble_max_tokens), (528, 528));
    }

    #[test]
    fn assemble_and_thread_arguments() {
        let p = AssembleParams::small();
        let t = token_limits(89.879905, &p).unwrap();
        let a = assemble_args(Path::new("/m"), 64, &p, &t).join(" ");
        assert_eq!(a, "assemble --reads-dir /m/input --k-start 4 --k-max 32 --graph-dir /m \
--threads 64 --tip-max-tokens 556 --superbubble-max-tokens 556");
        let a = assemble_args(Path::new("/m"), 64, &AssembleParams::large(), &t);
        assert_eq!(a[a.len() - 2..], ["--max-path-rows".to_string(), "1000000000".to_string()]);
        let cut = AssembleParams { abundance_cutoffs: "1,2".into(), ..AssembleParams::small() };
        let a = assemble_args(Path::new("/m"), 64, &cut, &t);
        assert_eq!(a[a.len() - 2..], ["--progressive-cutoffs".to_string(), "1,2".to_string()]);
        let th = thread_args(Path::new("/m"), Path::new("/o"), Path::new("/s"), Path::new("/q.gz"), 64, &p).join(" ");
        assert_eq!(th, "--threads 64 --graph-dir /m --k 32 --min-support 3 --max-conflict 0.2 --tangles --anchor-nodes 200 \
--anchor-edge 40 --tangle-dominance 3 --join --join-symbols /s --join-partners /q.gz --join-threshold 0.62 \
--join-min-nodes 200 --join-head 200 --join-min-paths 2 --output-dir /o");
    }

    #[test]
    fn prune_keeps_k_max_and_other_names() {
        let d = tmp("prune");
        let c = d.clone();
        fs::create_dir_all(c.join("input")).unwrap();
        for f in ["k32.unitigs.tsv.gz", "k32.node_tokens.u32le", "final.gfa", "k04.node_tokens.u32le", "k31.unitigs.gfa",
                  "kx.txt", "k5"] {
            fs::write(c.join(f), b"x").unwrap();
        }
        assert_eq!(prune_graph(&d, 32).unwrap(), 2);
        let mut left: Vec<String> = fs::read_dir(&c).unwrap().map(|e| e.unwrap().file_name().to_string_lossy().to_string()).collect();
        left.sort();
        assert_eq!(left, ["final.gfa", "input", "k32.node_tokens.u32le", "k32.unitigs.tsv.gz", "k5", "kx.txt"]);
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn final_graph_adds_circular_self_links() {
        let d = tmp("gfa");
        let c = d.clone();
        fs::write(c.join("k32.unitigs.gfa"), "H\tVN:Z:1.0\nS\t0\t*\tLN:i:40\tIA:i:3\nS\t1\t*\tLN:i:50\tIA:i:4\n").unwrap();
        let mut gz = GzEncoder::new(File::create(c.join("k32.unitigs.tsv.gz")).unwrap(), Compression::default());
        gz.write_all(b"k\tstate\tunitig_id\tcircular\tnodes\ttoken_count\tabundance2\tindegree\toutdegree\tmembers\n\
32\tP3\t0\t0\t9\t40\t6\t0\t0\t1,2\n32\tP3\t1\t1\t19\t50\t8\t0\t0\t3,4\n").unwrap();
        gz.finish().unwrap();
        write_final_graph(&d, 32, &d.join("g.gfa")).unwrap();
        assert_eq!(fs::read_to_string(d.join("g.gfa")).unwrap(),
                   "H\tVN:Z:1.0\nS\t0\t*\tLN:i:40\tIA:i:3\nS\t1\t*\tLN:i:50\tIA:i:4\nL\t1\t+\t1\t+\t0M\n");
        let _ = fs::remove_dir_all(&d);
    }
}

// ---------------------------------------------------------------- the pipeline's stages

use crate::pipeline::Run;
use crate::presets::Preset;

/// The assemble, thread and layout stages' values from a run's resolved values (the fixed ones from the preset).
pub fn params_of(r: &Run) -> Result<AssembleParams> {
    let base = match r.preset() {
        Preset::Small => AssembleParams::small(),
        Preset::Large => AssembleParams::large(),
    };
    let int = |n: &str| r.v(n).parse::<u32>().map_err(|_| format!("--{}: '{}' is not an integer", n, r.v(n)));
    let big = |n: &str| r.v(n).parse::<i64>().map_err(|_| format!("--{}: '{}' is not an integer", n, r.v(n)));
    let num = |n: &str| r.v(n).parse::<f64>().map_err(|_| format!("--{}: '{}' is not a number", n, r.v(n)));
    let k_max = int("k-max")?;
    let k_min = int("k-min")?;
    // the graph step's own path-row cap is 1e8: `--max-path-rows` is passed only for another value
    let rows: u64 = r.v("max-path-rows").parse().map_err(|_| "--max-path-rows: not an integer".to_string())?;
    Ok(AssembleParams {
        bases_per_second: r.v("bases-per-second").to_string(),
        k_min: <u16 as std::convert::TryFrom<u32>>::try_from(k_min).map_err(|_| "--k-min: too large".to_string())?,
        k_max: <u16 as std::convert::TryFrom<u32>>::try_from(k_max).map_err(|_| "--k-max: too large".to_string())?,
        abundance_cutoffs: r.v("abundance-cutoffs").to_string(),
        tip_max_bp: big("tip-max")?,
        superbubble_max_bp: big("superbubble-max")?,
        max_path_rows: if rows == 100_000_000 { None } else { Some(rows) },
        min_support: int("min-support")?,
        max_conflict: num("max-conflict")?,
        anchor_nodes: int("anchor-nodes")?,
        anchor_edge: int("anchor-edge")?,
        tangle_dominance: num("tangle-dominance")?,
        join_threshold: num("join-threshold")?,
        join_min_nodes: int("join-min-nodes")?,
        join_head: int("join-head")?,
        join_min_paths: int("join-paths")?,
        layout_max_gap: int("layout-max-gap")?,
        own_gap: int("layout-own-gap")?,
        tangles: r.cfg.tangles,
        joins: r.cfg.joins,
        ..base
    })
}

fn paths_of(r: &Run) -> Result<RunPaths> {
    // the correction's own rule: the one round with a corrected table and its sparse tables
    let round = crate::correct::final_round_dir(r)?;
    // the corrected and sparse tables of the final round, the uncorrected q-symbols and the partner table: where
    // the correction stage left them
    RunPaths::new(&r.root, &round)
}

/// The token limits of a run whose sparse scale exists.
fn limits_of(paths: &RunPaths, p: &AssembleParams) -> Result<TokenLimits> {
    let scale = paths.scale_dir().join("bp_per_token.tsv");
    let value = bp_per_token(&scale).map_err(|e| e.to_string())?;
    token_limits(value, p)
}

/// Stage `assemble`: the sparse table's bp scale, the graph and its pruning. The signal lengths were made by the signal
/// stage.
pub fn run(r: &Run) -> Result<()> {
    let p = params_of(r)?;
    let paths = paths_of(r)?;
    // a rerun of the stage starts from the correction's outputs (and the signal lengths)
    let scale = paths.scale_dir();
    for f in [scale.join("token_counts.tsv.gz"), scale.join("bp_per_token.tsv")] {
        if f.exists() {
            fs::remove_file(&f).map_err(err(&s(&f)))?;
        }
    }
    for d in [paths.graph_dir(), paths.threaded_dir()] {
        if d.exists() {
            fs::remove_dir_all(&d).map_err(err(&s(&d)))?;
        }
    }
    let limits = bp_scale_sparse(&r.runner, &paths, &p)?;
    let graph_dir = assemble_graph(&r.runner, &paths, r.cfg.threads, &p, &limits)?;
    prune_graph(&graph_dir, p.k_max)?;
    Ok(())
}

/// Stage `thread`: read threading, tangles and joins.
pub fn thread_stage(r: &Run) -> Result<()> {
    let p = params_of(r)?;
    let paths = paths_of(r)?;
    let out = paths.threaded_dir();
    if out.exists() {
        fs::remove_dir_all(&out).map_err(err(&s(&out)))?;
    }
    thread(&r.runner, &paths, &paths.graph_dir(), r.cfg.threads, &p)?;
    Ok(())
}

/// The final graph's folder: the threaded graph, or the graph with --no-threading.
fn final_bundle(r: &Run, paths: &RunPaths) -> PathBuf {
    if r.cfg.threading {
        paths.threaded_dir()
    } else {
        paths.graph_dir()
    }
}

/// Stage `layout`: the signal layout of the final graph (from the events' POD5 sample spans the signal stage wrote,
/// `signal/spans.bin`), and the output files: the final graph with its circular self-links, assembly_info.tsv,
/// read_paths.tsv.gz and layout.tsv.gz.
pub fn layout_stage(r: &Run) -> Result<()> {
    let p = params_of(r)?;
    let paths = paths_of(r)?;
    let limits = limits_of(&paths, &p)?;
    let graph = final_bundle(r, &paths);
    let spans = r.root.join(work::SPANS);
    let dir = r.root.join("layout");
    if dir.exists() {
        fs::remove_dir_all(&dir).map_err(err(&s(&dir)))?;
    }
    layout(&r.runner, &paths, &paths.graph_dir(), &graph, &spans, &dir, r.cfg.threads, &p)?;
    write_final_graph(&graph, p.k_max, &r.root.join("assembly_graph.gfa"))?;
    write_assembly_info(&graph, p.k_max, limits.bp_per_token, &r.root.join("assembly_info.tsv"))?;
    write_read_paths(&dir, &r.root.join("read_paths.tsv.gz"))?;
    gzip_file(&dir.join("layout.tsv"), &r.root.join("layout.tsv.gz"))?;
    Ok(())
}

/// assembly_info.tsv: one row per unitig of the final graph: id, tokens, an estimate of its length in bp (tokens x the
/// run's bp per token), its abundance (the GFA's IA tag), circular, and the graph's unitigs it was made from.
pub fn write_assembly_info(bundle: &Path, k: u16, bp_per_token: f64, out: &Path) -> Result<()> {
    let gfa = order_file(bundle, k as usize, GFA);
    let mut ia: HashMap<String, String> = HashMap::new();
    for line in fs::read_to_string(&gfa).map_err(err(&s(&gfa)))?.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.first() == Some(&"S") && f.len() >= 2 {
            let v = f.iter().find_map(|x| x.strip_prefix("IA:i:")).unwrap_or("");
            ia.insert(f[1].to_string(), v.to_string());
        }
    }
    let map = bundle.join("unitig_map.tsv");
    let mut text = String::from("unitig\ttokens\tapprox_bp\tabundance\tcircular\toriginal_unitigs\n");
    if map.is_file() {
        // the threaded bundle: unitig_id, tokens, circular, original_unitigs
        for line in fs::read_to_string(&map).map_err(err(&s(&map)))?.lines().skip(1) {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 4 {
                return Err(format!("{}: short row", map.display()));
            }
            let tokens: f64 = f[1].parse().map_err(err("unitig_map tokens"))?;
            text += &format!("{}\t{}\t{}\t{}\t{}\t{}\n", f[0], f[1], (tokens * bp_per_token).round() as u64,
                             ia.get(f[0]).map(String::as_str).unwrap_or(""), f[2], f[3]);
        }
    } else {
        // the graph (--no-threading): its final unitig table
        let path = order_file(bundle, k as usize, UNITIGS);
        let file = File::open(&path).map_err(err(&s(&path)))?;
        let mut body = String::new();
        BufReader::new(flate2::read::MultiGzDecoder::new(file)).read_to_string(&mut body).map_err(err(&s(&path)))?;
        let mut lines = body.lines();
        let head: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
        let col = |n: &str| head.iter().position(|h| *h == n).ok_or_else(|| format!("{}: no {}", path.display(), n));
        let (iu, it, ic) = (col("unitig_id")?, col("token_count")?, col("circular")?);
        for line in lines {
            let f: Vec<&str> = line.split('\t').collect();
            let tokens: f64 = f[it].parse().map_err(err("token_count"))?;
            text += &format!("{}\t{}\t{}\t{}\t{}\t{}\n", f[iu], f[it], (tokens * bp_per_token).round() as u64,
                             ia.get(f[iu]).map(String::as_str).unwrap_or(""), f[ic], f[iu]);
        }
    }
    fs::write(out, text).map_err(err(&s(out)))
}

/// read_paths.tsv.gz: every piece of a read on a unitig (the layout's pile: verified pieces, then unverified ones),
/// with a `verified` column.
pub fn write_read_paths(layout_dir: &Path, out: &Path) -> Result<()> {
    let mut w = GzEncoder::new(BufWriter::new(File::create(out).map_err(err(&s(out)))?), Compression::new(6));
    let mut header_done = false;
    for (name, flag) in [("pile.tsv.gz", "1"), ("pile_unverified.tsv.gz", "0")] {
        let path = layout_dir.join(name);
        if !path.is_file() {
            continue;
        }
        let file = File::open(&path).map_err(err(&s(&path)))?;
        let mut body = String::new();
        BufReader::new(flate2::read::MultiGzDecoder::new(file)).read_to_string(&mut body).map_err(err(&s(&path)))?;
        let mut lines = body.lines();
        let head = lines.next().unwrap_or("");
        if !header_done {
            writeln!(w, "{}\tverified", head).map_err(err(&s(out)))?;
            header_done = true;
        }
        for line in lines {
            writeln!(w, "{}\t{}", line, flag).map_err(err(&s(out)))?;
        }
    }
    w.finish().map_err(err(&s(out)))?.flush().map_err(err(&s(out)))
}

fn gzip_file(src: &Path, dst: &Path) -> Result<()> {
    let mut r = BufReader::new(File::open(src).map_err(err(&s(src)))?);
    let mut w = GzEncoder::new(BufWriter::new(File::create(dst).map_err(err(&s(dst)))?), Compression::new(6));
    std::io::copy(&mut r, &mut w).map_err(err(&s(dst)))?;
    w.finish().map_err(err(&s(dst)))?.flush().map_err(err(&s(dst)))
}

/// Every name `final_files` can give (a run without `--resume` removes them from OUT first).
pub const OUTPUT_NAMES: &[&str] =
    &["assembly_graph.gfa", "assembly_info.tsv", "read_paths.tsv.gz", "layout.tsv.gz", "joins.tsv", "stats.tsv"];

/// The files a user gets (source in the work folder, name in the output folder).
pub fn final_files(r: &Run) -> Result<Vec<(PathBuf, String)>> {
    let paths = paths_of(r)?;
    let bundle = final_bundle(r, &paths);
    let mut v = vec![
        (r.root.join("assembly_graph.gfa"), "assembly_graph.gfa".to_string()),
        (r.root.join("assembly_info.tsv"), "assembly_info.tsv".to_string()),
        (r.root.join("read_paths.tsv.gz"), "read_paths.tsv.gz".to_string()),
        (r.root.join("layout.tsv.gz"), "layout.tsv.gz".to_string()),
    ];
    if r.cfg.threading {
        if r.cfg.joins {
            v.push((bundle.join("joins.tsv"), "joins.tsv".to_string()));
        }
        v.push((bundle.join("stats.tsv"), "stats.tsv".to_string()));
    }
    Ok(v)
}
