//! Signal layout of the final graph (`squash __layout`): each unitig of the final K-max graph is described by the
//! stretches of raw signal that make it up, as miniasm lists the reads of a unitig (docs/OUTPUTS.md, `layout.tsv.gz`):
//!
//! - Paths: every read of the graph's input (`--reads-dir`: its assembly tokens and each token's q-symbol span)
//!   is matched to the graph by its K-token windows (a window counts only when it equals a node's tokens; it is looked
//!   up at the first node of its hash only).
//! - Pile: a run of a read's matched windows that follows one unitig collinearly (one offset node - window, modulo N
//!   on a circular unitig; at most --max-gap unmatched windows between consecutive matches) is a piece; a piece made
//!   only of nodes that occur at several places of the graph is dropped.
//! - Own signal (--own-gap Q > 0): a corrected token that equals the read's own uncorrected token with the same
//!   q-symbol span (`--uncorrected-dir`, the correction's input) is an own token (the full-density final round's table,
//!   `--corrected-tokens`). Each piece is split where consecutive own tokens are more than Q q-symbols apart, trimmed
//!   to its first and last own token and kept with at least --own-min own tokens and one whole K-min-mer (verified).
//! - Tiling: from the unitig's first node, the piece that contains the current node and reaches furthest is taken
//!   (verified pieces first, unverified ones only where no verified piece reaches past the current node), until the
//!   unitig's end (a circle: back to its first node); nodes no piece covers are gaps.
//! - Coordinates: each piece's contribution in the unitig (nodes, tokens) and in its read (tokens, q-symbols = events
//!   and POD5 signal samples from --spans; end exclusive).
//!
//! Output (--output-dir, written into a staging folder and renamed once complete): layout.tsv, pile.tsv.gz (with
//! --own-gap the verified pieces) and pile_unverified.tsv.gz (with --own-gap). The spans file is streamed once for the
//! samples the rows need (a read listed twice takes its last record).

use crate::graph::{NODE_TOKENS, UNITIGS, order_file};
use crate::io::{READS_TABLE, Result, TOKEN_SPANS, TOKEN_VALUES, read_u32le, read_u64_pairs, tsv_rows};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::convert::TryInto;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};


const LAYOUT_COLUMNS: [&str; 18] = [
    "unitig", "circular", "nodes", "rank", "read_id", "read_ordinal", "node_from", "node_to", "token_from", "token_to",
    "read_token_from", "read_token_to", "q_from", "q_to", "sample_from", "sample_to", "piece_node_from",
    "piece_node_to",
];
const PILE_HEADER: &str = "unitig\tread_id\tnode_from\tnode_to\toffset\n";

// ---------------------------------------------------------------- graph and reads

/// The final graph's unitigs as node lists (in the unitig table's order; a unitig id listed twice keeps its first
/// place and its last row) and where each node occurs: (unitig index, node index).
pub struct Unitigs {
    pub k: usize,
    pub ids: Vec<i64>,
    pub members: Vec<Vec<u32>>,
    pub circular: Vec<bool>,
    places: HashMap<u32, Vec<(usize, i64)>>,
}

impl Unitigs {
    pub fn new(k: usize, rows: Vec<(i64, Vec<u32>, bool)>) -> Unitigs {
        let (mut ids, mut members, mut circular) = (Vec::new(), Vec::new(), Vec::new());
        let mut at: HashMap<i64, usize> = HashMap::new();
        for (u, m, c) in rows {
            match at.get(&u) {
                Some(&i) => {
                    members[i] = m;
                    circular[i] = c;
                }
                None => {
                    at.insert(u, ids.len());
                    ids.push(u);
                    members.push(m);
                    circular.push(c);
                }
            }
        }
        let mut places: HashMap<u32, Vec<(usize, i64)>> = HashMap::new();
        for (x, m) in members.iter().enumerate() {
            for (i, &node) in m.iter().enumerate() {
                places.entry(node).or_default().push((x, i as i64));
            }
        }
        Unitigs { k, ids, members, circular, places }
    }

    fn unique(&self, node: u32) -> bool {
        self.places.get(&node).map_or(false, |p| p.len() == 1)
    }

    fn n(&self, x: usize) -> i64 {
        self.members[x].len() as i64
    }

    /// Tokens the layout covers: N + k - 1 for a linear unitig, N for a circular one (each token once).
    pub fn tokens(&self, x: usize) -> i64 {
        let n = self.n(x);
        if self.circular[x] {
            n
        } else {
            n + self.k as i64 - 1
        }
    }
}

/// The order-k graph of a graph folder and its node tokens.
pub fn load_unitigs(order_dir: &Path, k: usize) -> Result<(Unitigs, Vec<u32>)> {
    let node_tokens = read_u32le(&order_file(order_dir, k, NODE_TOKENS))?;
    let mut rows = Vec::new();
    tsv_rows(&order_file(order_dir, k, UNITIGS), &["unitig_id", "members", "token_count", "circular"], |r| {
        let m: Vec<u32> = r[1].split(',').map(|x| x.parse()).collect::<std::result::Result<_, _>>()?;
        if m.len() + k - 1 != r[2].parse::<usize>()? {
            return Err(format!("unitig {}: members and token_count disagree", r[0]).into());
        }
        rows.push((r[0].parse::<i64>()?, m, r[3] == "1"));
        Ok(())
    })?;
    Ok((Unitigs::new(k, rows), node_tokens))
}

/// One read of the graph's input: id, ordinal, and the range of its tokens / spans in the arrays.
pub struct GraphRead {
    pub id: String,
    pub ordinal: i64,
    pub range: (usize, usize),
}

/// The graph's input reads (a folder of token arrays): the reads, the token values and their q-symbol spans.
pub fn load_reads(reads_dir: &Path) -> Result<(Vec<GraphRead>, Vec<u32>, Vec<[u64; 2]>)> {
    let tokens = read_u32le(&reads_dir.join(TOKEN_VALUES))?;
    let spans = read_u64_pairs(&reads_dir.join(TOKEN_SPANS))?;
    if spans.len() != tokens.len() {
        return Err("token and span arrays differ in length".into());
    }
    let mut reads = Vec::new();
    let n = tokens.len();
    tsv_rows(&reads_dir.join(READS_TABLE), &["read_id", "read_ordinal", "token_offset", "token_count"], |r| {
        let (o, c): (usize, usize) = (r[2].parse()?, r[3].parse()?);
        // a range past the end of the arrays is cut at their end
        reads.push(GraphRead { id: r[0].to_string(), ordinal: r[1].parse()?, range: (o.min(n), (o + c).min(n).max(o.min(n))) });
        Ok(())
    })?;
    Ok((reads, tokens, spans))
}

// ---------------------------------------------------------------- pile

/// The graph's nodes keyed by their K-token window hash, sorted by (hash, node).
fn node_index(graph: &Unitigs, node_tokens: &[u32]) -> Result<Vec<(u64, u32)>> {
    let k = graph.k;
    let rows = node_tokens.len() / k;
    if node_tokens.len() % k != 0 {
        return Err("the atomic tokens are not whole K-token nodes".into());
    }
    let mut nodes: Vec<u32> = graph.places.keys().copied().collect();
    nodes.sort_unstable();
    if let Some(&last) = nodes.last() {
        if last as usize >= rows {
            return Err(format!("atomic node {} out of range", last).into());
        }
    }
    let mut keyed: Vec<(u64, u32)> = nodes
        .par_iter()
        .map(|&x| (crate::threading::window_hash(&node_tokens[x as usize * k..(x as usize + 1) * k]), x))
        .collect();
    keyed.par_sort_unstable();
    Ok(keyed)
}

/// (window index, node) of every K-token window of `tok` equal to a graph node, in window order. A window is
/// compared with the first node of its hash only.
fn read_matches(keyed: &[(u64, u32)], node_tokens: &[u32], k: usize, tok: &[u32]) -> Vec<(i64, u32)> {
    let mut out = Vec::new();
    if tok.len() < k || keyed.is_empty() {
        return out;
    }
    let top = (1..k).fold(1u64, |p, _| p.wrapping_mul(1_099_511_628_211));
    let mut h = crate::threading::window_hash(&tok[..k]);
    for w in 0..tok.len() - k + 1 {
        if w > 0 {
            h = h
                .wrapping_sub((tok[w - 1] as u64 + 1).wrapping_mul(top))
                .wrapping_mul(1_099_511_628_211)
                .wrapping_add(tok[w + k - 1] as u64 + 1);
        }
        let j = keyed.partition_point(|x| x.0 < h);
        if j == keyed.len() || keyed[j].0 != h {
            continue;
        }
        let x = keyed[j].1 as usize;
        if node_tokens[x * k..(x + 1) * k] == tok[w..w + k] {
            out.push((w as i64, keyed[j].1));
        }
    }
    out
}

/// The read's pieces (unitig index, a, b, w0): collinear runs of matches (node index = window + one offset, modulo N
/// on a circular unitig) with at most `max_gap` unmatched windows between consecutive matches and at least one node
/// that occurs once in the graph. a..b is the run's node range (b inclusive and unwrapped on a circle), w0 the read
/// window at node a. Output order: a run as a later match of the same offset is too far from it (that match starts a
/// new run in its place), then the open runs in the order their offsets first appeared.
fn pieces_of_read(graph: &Unitigs, matches: &[(i64, u32)], max_gap: i64) -> Vec<(usize, i64, i64, i64)> {
    struct Run {
        a: i64,
        b: i64,
        w_first: i64,
        w_last: i64,
        uniq: bool,
    }
    let mut runs: Vec<((usize, i64), Run)> = Vec::new();
    let mut at: HashMap<(usize, i64), usize> = HashMap::new();
    let mut out = Vec::new();
    for &(w, node) in matches {
        let uniq = graph.unique(node);
        let Some(places) = graph.places.get(&node) else { continue };
        for &(u, i) in places {
            let n = graph.n(u);
            let d = if graph.circular[u] { (i - w).rem_euclid(n) } else { i - w };
            let key = (u, d);
            if let Some(&slot) = at.get(&key) {
                let run = &mut runs[slot].1;
                if w - run.w_last <= max_gap + 1 {
                    run.b = run.a + (w - run.w_first);
                    run.w_last = w;
                    run.uniq = run.uniq || uniq;
                    continue;
                }
                if run.uniq {
                    out.push((u, run.a, run.b, run.w_first));
                }
                runs[slot].1 = Run { a: i, b: i, w_first: w, w_last: w, uniq };
            } else {
                at.insert(key, runs.len());
                runs.push((key, Run { a: i, b: i, w_first: w, w_last: w, uniq }));
            }
        }
    }
    for ((u, _), run) in &runs {
        if run.uniq {
            out.push((*u, run.a, run.b, run.w_first));
        }
    }
    out
}

/// Binary searches on a non-decreasing slice: `search_left` = the first index with v[i] >= x, `search_right` = the
/// first with v[i] > x.
fn search_left(v: &[i64], x: i64) -> usize {
    v.partition_point(|&y| y < x)
}

fn search_right(v: &[i64], x: i64) -> usize {
    v.partition_point(|&y| y <= x)
}

/// The verified sub-pieces (a', b', w0') of piece (a, b, w0): the read's own token starts inside the piece's read
/// tokens, split where consecutive ones are more than `max_gap` q-symbols apart; each group with >= `min_own` own
/// tokens becomes the read tokens whose starts lie between its first and last own token, kept when they hold a whole
/// K-token window. `starts` / `ends` = the read's token q-symbol starts / ends.
fn verified_pieces(a: i64, b: i64, w0: i64, starts: &[i64], ends: &[i64], own: &[i64], k: i64, max_gap: i64,
                   min_own: usize) -> Vec<(i64, i64, i64)> {
    let (t_lo, t_hi) = (w0, w0 + (b - a) + k);
    let (q_lo, q_hi) = (starts[t_lo as usize], ends[(t_hi - 1) as usize]);
    // own is sorted: the starts inside [q_lo, q_hi) are one slice of it
    let x = &own[search_left(own, q_lo)..search_left(own, q_hi).max(search_left(own, q_lo))];
    let mut out = Vec::new();
    if x.is_empty() {
        return out;
    }
    let mut lo = 0usize;
    for hi in 1..=x.len() {
        if hi < x.len() && x[hi] - x[hi - 1] <= max_gap {
            continue;
        }
        let g = &x[lo..hi];
        lo = hi;
        if g.len() < min_own {
            continue;
        }
        let i0 = t_lo.max(search_left(starts, g[0]) as i64);
        let i1 = t_hi.min(search_right(starts, g[g.len() - 1]) as i64);
        if i1 - i0 < k {
            continue;
        }
        out.push((a + (i0 - w0), a + (i1 - k - w0), i0));
    }
    out
}

// ---------------------------------------------------------------- tiling

/// Greedy tiling of nodes 0..n-1 (and back to node 0, i.e. n, for a circle) by pieces (a, b, id), b inclusive and
/// unwrapped; `fallback` pieces only where no piece of `pieces` reaches past the current node. Returns (piece id,
/// takeover node, a, b) in order, a and b as the tiling used them (a circle's pieces also shifted by -n / +n).
pub fn tile(pieces: &[(i64, i64, usize)], n: i64, circular: bool, fallback: &[(i64, i64, usize)]) -> Vec<(usize, i64, i64, i64)> {
    let shifted = |ps: &[(i64, i64, usize)]| -> Vec<(i64, i64, usize)> {
        let mut cand: Vec<(i64, i64, usize)> = ps.to_vec();
        if circular {
            cand.extend(ps.iter().filter(|&&(_, b, _)| b >= n).map(|&(a, b, p)| (a - n, b - n, p)));
            cand.extend(ps.iter().map(|&(a, b, p)| (a + n, b + n, p)));
        }
        cand
    };
    let mut tiers = vec![shifted(pieces)];
    if !fallback.is_empty() {
        tiers.push(shifted(fallback));
    }
    let hand_back = |tiers: &Vec<Vec<(i64, i64, usize)>>, x: i64, b: i64| -> i64 {
        tiers[0].iter().filter(|&&(a, _, _)| x < a && a <= b).map(|&(a, _, _)| a).min().unwrap_or(b)
    };
    let end = if circular { n } else { n - 1 };
    let (mut chosen, mut cur): (Vec<(usize, i64, i64, i64)>, i64) = (Vec::new(), 0);
    loop {
        let mut best: Option<(i64, i64, usize)> = None;
        let mut tier = 0usize;
        for (t, cand) in tiers.iter().enumerate() {
            tier = t;
            for &(a, b, p) in cand {
                if a <= cur && cur <= b && best.map_or(true, |x| b > x.1) {
                    best = Some((a, b, p));
                }
            }
            if let Some(x) = best {
                if chosen.is_empty() || x.1 > cur {
                    break;
                }
            }
        }
        match best {
            Some((a, b, p)) if chosen.is_empty() || b > cur => {
                chosen.push((p, cur, a, b));
                cur = if tier == 0 { b } else { hand_back(&tiers, cur, b) };
            }
            _ => {
                // the earliest start after cur: (a, tier, -b, p) smallest
                let mut next: Option<(i64, usize, i64, usize)> = None;
                for (t, cand) in tiers.iter().enumerate() {
                    for &(a, b, p) in cand {
                        if a > cur {
                            let key = (a, t, -b, p);
                            if next.map_or(true, |x| key < x) {
                                next = Some(key);
                            }
                        }
                    }
                }
                let Some((a, t, nb, p)) = next else { break };
                chosen.push((p, a, a, -nb));
                cur = if t == 0 { -nb } else { hand_back(&tiers, a, -nb) };
            }
        }
        if cur >= end {
            break;
        }
    }
    chosen
}

/// Token ranges of the tiling: (index in `chosen` or None for a gap, token from, token to), clipped to [0, total).
pub fn contributions(chosen: &[(usize, i64, i64, i64)], k: i64, total: i64) -> Vec<(Option<usize>, i64, i64)> {
    let (mut out, mut prev) = (Vec::new(), 0i64);
    for (j, &(_, x, _, b)) in chosen.iter().enumerate() {
        let start = x.max(prev);
        if start > prev {
            out.push((None, prev, start.min(total)));
        }
        let next = chosen.get(j + 1).map(|c| c.1);
        let stop = match next {
            Some(nx) if nx <= b => nx,
            Some(nx) => (b + k).min(nx),
            None => (b + k).min(total),
        };
        let stop = stop.min(total);
        if stop > start {
            out.push((Some(j), start, stop));
            prev = stop;
        }
        if prev >= total {
            break;
        }
    }
    if prev < total {
        out.push((None, prev, total));
    }
    out
}

// ---------------------------------------------------------------- own signal

/// Reads until `buf` is full or the stream ends; the bytes read.
fn read_full<R: Read>(r: &mut R, buf: &mut [u8]) -> std::io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match r.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(n)
}

/// Galloping lower bound: the first index with v[i] >= x, searched from `from` when v[from - 1] < x (the same index
/// as a binary search over the whole slice, which v must allow: non-decreasing).
fn lower_bound_from(v: &[u64], x: u64, from: usize) -> usize {
    if from == 0 || from > v.len() || v[from - 1] >= x {
        return v.partition_point(|&y| y < x);
    }
    let (mut lo, mut step) = (from, 1usize);
    let mut hi = from;
    while hi < v.len() && v[hi] < x {
        lo = hi + 1;
        hi = hi.saturating_add(step).min(v.len());
        step = step.saturating_mul(2);
    }
    lo + v[lo..hi.min(v.len())].partition_point(|&y| y < x)
}

/// read id -> sorted q-symbol starts of its own tokens: the corrected tokens (`corrected`: the final round's full table,
/// columns 2, 5, 6, 7 = read ordinal, value, start, end) equal to an uncorrected token of the same read ordinal, start,
/// end and value (the folder `uncorrected`). An uncorrected token's ordinal is the number of distinct token offsets of the reads
/// after the first that are <= its index: an offset counts once, so after a read without tokens every later read's
/// tokens are keyed one ordinal below the read's own. A corrected token is compared with the first uncorrected token
/// of its (ordinal, start) key in key order.
pub fn own_token_starts(corrected: &Path, uncorrected: &Path) -> Result<HashMap<String, Vec<i64>>> {
    let (mut ids, mut offsets) = (Vec::new(), Vec::new());
    tsv_rows(&uncorrected.join(READS_TABLE), &["read_id", "read_ordinal", "token_offset"], |r| {
        if r[1].parse::<usize>()? != ids.len() {
            return Err("prepared reads are not in ordinal order".into());
        }
        ids.push(r[0].to_string());
        offsets.push(r[2].parse::<i64>()?);
        Ok(())
    })?;
    let val = read_u32le(&uncorrected.join(TOKEN_VALUES))?;
    let n = val.len();
    let spans_path = uncorrected.join(TOKEN_SPANS);
    let spans_file = File::open(&spans_path).map_err(|e| format!("{}: {}", spans_path.display(), e))?;
    if spans_file.metadata()?.len() != 16 * n as u64 {
        return Err("prepared token and span arrays differ in length".into());
    }
    // a token's ordinal = the distinct offsets of the reads after the first that are <= its index (an offset listed
    // twice counts once)
    let mut bounds: Vec<usize> = Vec::new();
    for &o in offsets.iter().skip(1) {
        if o < 0 {
            return Err("negative token offset".into());
        }
        if (o as usize) < n {
            bounds.push(o as usize);
        }
    }
    bounds.sort_unstable();
    bounds.dedup();
    // key = (ordinal << 32) | start (the spans streamed), ends alongside; stably sorted by key (already sorted when
    // every read's starts are non-decreasing)
    let (mut key, mut end): (Vec<u64>, Vec<i64>) = (Vec::with_capacity(n), Vec::with_capacity(n));
    let mut rd = BufReader::with_capacity(1 << 22, spans_file);
    let (mut ord, mut next) = (0i64, 0usize);
    let mut pair = [0u8; 16];
    for i in 0..n {
        while next < bounds.len() && bounds[next] <= i {
            ord += 1;
            next += 1;
        }
        rd.read_exact(&mut pair)?;
        let (start, stop) = (u64::from_le_bytes(pair[..8].try_into().unwrap()), u64::from_le_bytes(pair[8..].try_into().unwrap()));
        if start >= 1u64 << 32 {
            return Err("q-symbol start beyond 2^32".into());
        }
        key.push(((ord << 32) | start as i64) as u64);
        end.push(stop as i64);
    }
    if !key.windows(2).all(|w| w[0] <= w[1]) {
        let mut order: Vec<usize> = (0..n).collect();
        order.par_sort_by_key(|&i| (key[i], i));
        key = order.iter().map(|&i| key[i]).collect();
        end = order.iter().map(|&i| end[i]).collect();
        let v: Vec<u32> = order.iter().map(|&i| val[i]).collect();
        return own_matches(corrected, &ids, &key, &end, &v);
    }
    own_matches(corrected, &ids, &key, &end, &val)
}

/// The own tokens of the corrected table against the uncorrected tokens' sorted keys (`own_token_starts`).
fn own_matches(corrected: &Path, ids: &[String], key: &[u64], end: &[i64], val: &[u32]) -> Result<HashMap<String, Vec<i64>>> {
    let n = key.len();
    if n == 0 {
        return Err("no uncorrected tokens".into());
    }
    // the corrected table, streamed (columns 2, 5, 6, 7 of every line after the header): one thread decompresses
    // blocks of whole lines, the pool parses and matches them (found is sorted below: its order does not matter)
    let file = File::open(corrected).map_err(|e| format!("{}: {}", corrected.display(), e))?;
    let (tx, rx) = std::sync::mpsc::sync_channel::<std::io::Result<Vec<u8>>>(16);
    let reader = std::thread::spawn(move || {
        let mut gz = MultiGzDecoder::new(file);
        let mut rest: Vec<u8> = Vec::new();
        let mut header = true;
        loop {
            let mut block = std::mem::take(&mut rest);
            let have = block.len();
            block.resize(have + (8 << 20), 0);
            let got = match read_full(&mut gz, &mut block[have..]) {
                Ok(n) => n,
                Err(e) => {
                    let _ = tx.send(Err(e));
                    return;
                }
            };
            block.truncate(have + got);
            if header {
                match block.iter().position(|&c| c == b'\n') {
                    Some(i) => {
                        block.drain(..=i);
                        header = false;
                    }
                    None if got == 0 => return,
                    None => {
                        rest = block;
                        continue;
                    }
                }
            }
            if got == 0 {
                if !block.is_empty() && tx.send(Ok(block)).is_err() {
                    return;
                }
                return;
            }
            match block.iter().rposition(|&c| c == b'\n') {
                Some(i) => {
                    rest = block.split_off(i + 1);
                    if tx.send(Ok(block)).is_err() {
                        return;
                    }
                }
                None => rest = block,
            }
        }
    });
    let parse_block = |block: &[u8]| -> Result<Vec<(i64, i64)>> {
        let mut out = Vec::new();
        let mut last = 0usize;
        let parse = |b: &[u8]| -> Result<i64> {
            std::str::from_utf8(b)?.trim().parse::<i64>().map_err(|e| format!("{}: {}", corrected.display(), e).into())
        };
        for text in block.split(|&c| c == b'\n') {
            if text.iter().all(|c| c.is_ascii_whitespace()) {
                continue;
            }
            let mut f = text.split(|&c| c == b'\t');
            let _ = f.next();
            let o = parse(f.next().ok_or("short corrected row")?)?;
            let _ = f.next();
            let _ = f.next();
            let v = parse(f.next().ok_or("short corrected row")?)?;
            let s = parse(f.next().ok_or("short corrected row")?)?;
            let e = parse(f.next().ok_or("short corrected row")?)?;
            let k = (o << 32) | s;
            if k < 0 {
                continue; // below every key: the first uncorrected token differs
            }
            let j = lower_bound_from(key, k as u64, last).min(n - 1);
            last = j;
            if key[j] == k as u64 && end[j] == e && val[j] as i64 == v {
                out.push((o, s));
            }
        }
        Ok(out)
    };
    let mut found: Vec<(i64, i64)> = Vec::new();
    let batch = 2 * rayon::current_num_threads().max(1);
    let mut blocks: Vec<Vec<u8>> = Vec::with_capacity(batch);
    let mut done = false;
    while !done {
        blocks.clear();
        while blocks.len() < batch {
            match rx.recv() {
                Ok(Ok(b)) => blocks.push(b),
                Ok(Err(e)) => return Err(format!("{}: {}", corrected.display(), e).into()),
                Err(_) => {
                    done = true;
                    break;
                }
            }
        }
        let parts: Vec<Result<Vec<(i64, i64)>>> = blocks.par_iter().map(|b| parse_block(b)).collect();
        for part in parts {
            found.extend(part?);
        }
    }
    reader.join().map_err(|_| "the corrected table's reader failed")?;
    if !found.windows(2).all(|w| w[0] <= w[1]) {
        found.par_sort_unstable();
    }
    let mut out: HashMap<String, Vec<i64>> = HashMap::new();
    let mut i = 0;
    while i < found.len() {
        let o = found[i].0;
        let mut j = i;
        while j < found.len() && found[j].0 == o {
            j += 1;
        }
        let id = ids.get(o as usize).ok_or("own token of an unknown read ordinal")?;
        out.insert(id.clone(), found[i..j].iter().map(|x| x.1).collect());
        i = j;
    }
    Ok(out)
}

// ---------------------------------------------------------------- layout

/// One pile entry of a unitig: (a, b, w0, read index, verified).
type Piece = (i64, i64, i64, usize, bool);

/// A layout row before the POD5 samples: the 18 columns (sample columns filled later) and, with own tokens,
/// (verified, own tokens).
pub struct Row {
    pub unitig: i64,
    pub circular: bool,
    pub nodes: i64,
    pub rank: usize,
    pub read: Option<usize>,
    pub x: i64,
    pub y: i64,
    pub t0: i64,
    pub t1: i64,
    pub r0: i64,
    pub r1: i64,
    pub q0: i64,
    pub q1: i64,
    pub s0: i64,
    pub s1: i64,
    pub a: i64,
    pub b: i64,
    pub check: Option<(bool, usize)>,
}

pub struct Layout {
    pub rows: Vec<Row>,
    /// (unitig, read index, a, b, w0): the pile (with own tokens: the verified pieces) and the unverified pieces
    pub pile: Vec<(i64, usize, i64, i64, i64)>,
    pub unverified: Vec<(i64, usize, i64, i64, i64)>,
}

/// The layout rows and pile of every unitig, in unitig id order.
pub fn layout(graph: &Unitigs, node_tokens: &[u32], reads: &[GraphRead], tokens: &[u32], spans: &[[u64; 2]], max_gap: i64,
              own: Option<&HashMap<String, Vec<i64>>>, own_gap: i64, own_min: usize) -> Result<Layout> {
    let k = graph.k;
    let keyed = node_index(graph, node_tokens)?;
    let check = own.is_some() && own_gap > 0;
    let empty: Vec<i64> = Vec::new();
    // paths, pile and own signal per read (in parallel; collected in read order)
    let per_read: Vec<Vec<(usize, Piece)>> = reads
        .par_iter()
        .enumerate()
        .map(|(r, read)| {
            let tok = &tokens[read.range.0..read.range.1];
            let qs = &spans[read.range.0..read.range.1];
            let matches = read_matches(&keyed, node_tokens, k, tok);
            let mut out = Vec::new();
            let (starts, ends): (Vec<i64>, Vec<i64>) = if check {
                qs.iter().map(|s| (s[0] as i64, s[1] as i64)).unzip()
            } else {
                (Vec::new(), Vec::new())
            };
            // the binary searches need sorted token starts (the graph step refuses reads whose starts decrease)
            if !starts.windows(2).all(|w| w[0] <= w[1]) {
                return Err(format!("read {}: token starts decrease", read.id));
            }
            for (u, a, b, w0) in pieces_of_read(graph, &matches, max_gap) {
                if !check {
                    out.push((u, (a, b, w0, r, true)));
                    continue;
                }
                let n = graph.n(u);
                let o = own.unwrap().get(&read.id).unwrap_or(&empty);
                let sub = verified_pieces(a, b, w0, &starts, &ends, o, k as i64, own_gap, own_min);
                for &(a2, b2, w2) in &sub {
                    let (a2, b2) = if graph.circular[u] && a2 >= n { (a2 - n, b2 - n) } else { (a2, b2) };
                    out.push((u, (a2, b2, w2, r, true)));
                }
                if sub.is_empty() || sub.iter().map(|&(x, y, _)| y - x + 1).sum::<i64>() < b - a + 1 {
                    out.push((u, (a, b, w0, r, false)));
                }
            }
            Ok(out)
        })
        .collect::<std::result::Result<Vec<_>, String>>()?;
    let mut pile: Vec<Vec<Piece>> = vec![Vec::new(); graph.ids.len()];
    for list in per_read {
        for (u, p) in list {
            pile[u].push(p);
        }
    }
    let mut order: Vec<usize> = (0..graph.ids.len()).collect();
    order.sort_by_key(|&x| graph.ids[x]);
    // tiling and coordinates per unitig (in parallel; collected in id order)
    let per_unitig: Vec<Result<Vec<Row>>> = order
        .par_iter()
        .map(|&x| {
            let (n, circ, u) = (graph.n(x), graph.circular[x], graph.ids[x]);
            let ps = &pile[x];
            let chosen = if check {
                let v: Vec<(i64, i64, usize)> = ps.iter().enumerate().filter(|(_, p)| p.4).map(|(j, p)| (p.0, p.1, j)).collect();
                let f: Vec<(i64, i64, usize)> = ps.iter().enumerate().filter(|(_, p)| !p.4).map(|(j, p)| (p.0, p.1, j)).collect();
                tile(&v, n, circ, &f)
            } else {
                let v: Vec<(i64, i64, usize)> = ps.iter().enumerate().map(|(j, p)| (p.0, p.1, j)).collect();
                tile(&v, n, circ, &[])
            };
            let total = graph.tokens(x);
            let mut rows = Vec::new();
            for (rank, (j, t0, t1)) in contributions(&chosen, k as i64, total).into_iter().enumerate() {
                let (xn, yn) = (t0.min(n), t1.min(n));
                let Some(j) = j else {
                    rows.push(Row { unitig: u, circular: circ, nodes: n, rank, read: None, x: xn, y: yn, t0, t1, r0: -1,
                                    r1: -1, q0: -1, q1: -1, s0: -1, s1: -1, a: -1, b: -1,
                                    check: if check { Some((false, 0)) } else { None } });
                    continue;
                };
                let (p, _, a, b) = chosen[j];
                let (_, _, w0, r, verified) = ps[p];
                let read = &reads[r];
                let len = (read.range.1 - read.range.0) as i64;
                let (r0, r1) = (t0 - a + w0, t1 - a + w0);
                if !(0 <= r0 && r0 < r1 && r1 <= len) {
                    return Err(format!("unitig {}: read {} tokens [{}, {}) outside the read ({})", u, read.id, r0, r1, len).into());
                }
                let qs = &spans[read.range.0..read.range.1];
                let (q0, q1) = (qs[r0 as usize][0] as i64, qs[(r1 - 1) as usize][1] as i64);
                let check_cols = if check {
                    let o = own.unwrap().get(&read.id).unwrap_or(&empty);
                    Some((verified, o.iter().filter(|&&v| v >= q0 && v < q1).count()))
                } else {
                    None
                };
                rows.push(Row { unitig: u, circular: circ, nodes: n, rank, read: Some(r), x: xn, y: yn, t0, t1, r0, r1, q0,
                                q1, s0: -1, s1: -1, a, b, check: check_cols });
            }
            Ok(rows)
        })
        .collect();
    let mut rows = Vec::new();
    for r in per_unitig {
        rows.extend(r?);
    }
    let (mut verified_pile, mut unverified) = (Vec::new(), Vec::new());
    for &x in &order {
        for &(a, b, w0, r, v) in &pile[x] {
            if !check || v {
                verified_pile.push((graph.ids[x], r, a, b, w0));
            } else {
                unverified.push((graph.ids[x], r, a, b, w0));
            }
        }
    }
    Ok(Layout { rows, pile: verified_pile, unverified })
}

/// Fills the rows' POD5 sample columns from the signal step's --dump-event-spans file (per read: name length u32,
/// name, a u32, the event count n u32, n event first samples u32, n event end samples u32), streamed once: s0 = the
/// first sample of event q0, s1 = the end of event q1 - 1, when 0 <= q0 < q1 <= n (a read listed twice: its last
/// record).
/// A read the file does not hold keeps -1.
pub fn fill_samples(path: &Path, rows: &mut [Row], reads: &[GraphRead]) -> Result<()> {
    let mut want: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, row) in rows.iter().enumerate() {
        if let Some(r) = row.read {
            want.entry(reads[r].id.as_str()).or_default().push(i);
        }
    }
    let all: HashSet<&str> = reads.iter().map(|r| r.id.as_str()).collect();
    let file = File::open(path).map_err(|e| format!("{}: {}", path.display(), e))?;
    let size = file.metadata()?.len();
    let mut f = BufReader::with_capacity(1 << 22, file);
    let mut pos = 0u64;
    let mut u32buf = [0u8; 4];
    let mut name = Vec::new();
    let mut values: Vec<(i64, i64)> = Vec::new();
    while pos < size {
        f.read_exact(&mut u32buf)?;
        let nl = u32::from_le_bytes(u32buf) as usize;
        name.resize(nl, 0);
        f.read_exact(&mut name)?;
        let mut two = [0u8; 8];
        f.read_exact(&mut two)?;
        let n = u32::from_le_bytes(two[4..8].try_into().unwrap()) as usize;
        pos += 4 + nl as u64 + 8;
        let id = std::str::from_utf8(&name).map_err(|_| "spans file: a read name is not UTF-8")?;
        let rows_of = if all.contains(id) { want.get(id) } else { None };
        match rows_of {
            Some(list) => {
                let mut st = vec![0u8; 4 * n];
                let mut en = vec![0u8; 4 * n];
                f.read_exact(&mut st)?;
                f.read_exact(&mut en)?;
                let at = |b: &[u8], i: usize| u32::from_le_bytes(b[4 * i..4 * i + 4].try_into().unwrap()) as i64;
                values.clear();
                for &i in list {
                    let (q0, q1) = (rows[i].q0, rows[i].q1);
                    values.push(if 0 <= q0 && q0 < q1 && q1 as usize <= n {
                        (at(&st, q0 as usize), at(&en, (q1 - 1) as usize))
                    } else {
                        (-1, -1)
                    });
                }
                for (&i, &(s0, s1)) in list.iter().zip(values.iter()) {
                    rows[i].s0 = s0;
                    rows[i].s1 = s1;
                }
            }
            None => {
                f.seek_relative(8 * n as i64)?;
            }
        }
        pos += 8 * n as u64;
    }
    if pos != size {
        return Err(format!("{}: truncated spans record", path.display()).into());
    }
    Ok(())
}

// ---------------------------------------------------------------- output

fn write_gz(path: &Path, header: &str, rows: &[(i64, usize, i64, i64, i64)], reads: &[GraphRead]) -> Result<()> {
    let mut gz = GzEncoder::new(BufWriter::with_capacity(1 << 20, File::create(path)?), Compression::default());
    gz.write_all(header.as_bytes())?;
    for &(u, r, a, b, w0) in rows {
        writeln!(gz, "{}\t{}\t{}\t{}\t{}", u, reads[r].id, a, b, w0)?;
    }
    gz.finish()?.flush()?;
    Ok(())
}

/// layout.tsv, pile.tsv.gz and, with own tokens (`check`), pile_unverified.tsv.gz into `dir`.
pub fn write(dir: &Path, reads: &[GraphRead], lay: &Layout, check: bool) -> Result<()> {
    let mut out = BufWriter::with_capacity(1 << 20, File::create(dir.join("layout.tsv"))?);
    let mut head: Vec<&str> = LAYOUT_COLUMNS.to_vec();
    if check {
        head.extend(["verified", "own_tokens"]);
    }
    writeln!(out, "{}", head.join("\t"))?;
    for r in &lay.rows {
        let (id, ord) = match r.read {
            Some(i) => (reads[i].id.as_str(), reads[i].ordinal),
            None => ("-", -1),
        };
        write!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", r.unitig, r.circular as i32,
               r.nodes, r.rank, id, ord, r.x, r.y, r.t0, r.t1, r.r0, r.r1, r.q0, r.q1, r.s0, r.s1, r.a, r.b)?;
        match r.check {
            Some((v, o)) if r.read.is_some() => writeln!(out, "\t{}\t{}", v as i32, o)?,
            Some(_) => writeln!(out, "\t-1\t-1")?,
            None => writeln!(out)?,
        }
    }
    out.flush()?;
    write_gz(&dir.join("pile.tsv.gz"), PILE_HEADER, &lay.pile, reads)?;
    if check {
        write_gz(&dir.join("pile_unverified.tsv.gz"), PILE_HEADER, &lay.unverified, reads)?;
    }
    Ok(())
}

// ---------------------------------------------------------------- main

const USAGE: &str = "usage: squash __layout --graph-dir GRAPH_BUNDLE --reads-dir READS_BUNDLE --spans SPANS.bin \
--corrected-tokens TOKENS.tsv.gz --uncorrected-dir TOKENS_DIR --output-dir OUT [--k 32] [--max-gap 32] [--own-gap 0] \
[--own-min 2] [--threads N]";

pub fn run() -> Result<()> {
    let mut map: HashMap<String, String> = HashMap::new();
    let mut argv = crate::argv::args().skip(1);
    while let Some(key) = argv.next() {
        let Some(name) = key.strip_prefix("--") else { return Err(format!("unexpected argument {:?}\n{}", key, USAGE).into()) };
        if !["graph-dir", "reads-dir", "spans", "k", "max-gap", "own-gap", "own-min", "corrected-tokens", "uncorrected-dir",
             "output-dir", "threads"].contains(&name) {
            return Err(format!("unknown option {}\n{}", key, USAGE).into());
        }
        let value = argv.next().ok_or_else(|| format!("{} needs a value\n{}", key, USAGE))?;
        map.insert(name.to_string(), value);
    }
    let path = |k: &str| map.get(k).map(PathBuf::from);
    let int = |k: &str, d: i64| -> Result<i64> { Ok(match map.get(k) { Some(v) => v.parse()?, None => d }) };
    let graph_dir = path("graph-dir").ok_or_else(|| format!("missing --graph-dir\n{}", USAGE))?;
    let output = path("output-dir").ok_or_else(|| format!("missing --output-dir\n{}", USAGE))?;
    let (k, max_gap, own_gap, own_min) = (int("k", 32)?, int("max-gap", 32)?, int("own-gap", 0)?, int("own-min", 2)?);
    if k < 1 || own_min < 0 {
        return Err("--k must be at least 1 and --own-min at least 0".into());
    }
    let threads = int("threads", 0)?;
    if threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(threads as usize).build_global()?;
    }
    if output.exists() {
        return Err(format!("refusing to overwrite {}", output.display()).into());
    }
    let (graph, node_tokens) = load_unitigs(&graph_dir, k as usize)?;
    let src = path("reads-dir").ok_or_else(|| format!("missing --reads-dir\n{}", USAGE))?;
    let (reads, tokens, spans) = load_reads(&src)?;
    let required = |k: &str| path(k).ok_or_else(|| format!("missing --{}\n{}", k, USAGE));
    let (spans_file, corrected, uncorrected) = (required("spans")?, required("corrected-tokens")?, required("uncorrected-dir")?);
    let own = if own_gap > 0 { Some(own_token_starts(&corrected, &uncorrected)?) } else { None };
    let mut lay = layout(&graph, &node_tokens, &reads, &tokens, &spans, max_gap, own.as_ref(), own_gap, own_min as usize)?;
    fill_samples(&spans_file, &mut lay.rows, &reads)?;
    let parent = match output.parent() {
        Some(p) if !p.as_os_str().is_empty() => p.to_path_buf(),
        _ => PathBuf::from("."),
    };
    fs::create_dir_all(&parent)?;
    let name = output.file_name().ok_or("--output-dir has no name")?.to_string_lossy().to_string();
    let nanos = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)?.subsec_nanos();
    let stage = parent.join(format!(".{}.tmp.{}{:08x}", name, std::process::id(), nanos));
    fs::create_dir(&stage)?;
    match write(&stage, &reads, &lay, own.is_some()) {
        Ok(()) => {
            fs::rename(&stage, &output)?;
            let gaps = lay.rows.iter().filter(|r| r.read.is_none()).count();
            println!("{} unitigs, {} layout pieces from {} pile pieces, {} gaps -> {}", graph.ids.len(),
                     lay.rows.len() - gaps, lay.pile.len(), gaps, output.display());
            Ok(())
        }
        Err(e) => {
            let _ = fs::remove_dir_all(&stage);
            Err(e)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(rows: Vec<(i64, Vec<u32>, bool)>, k: usize) -> Unitigs {
        Unitigs::new(k, rows)
    }

    #[test]
    fn duplicate_unitig_ids_keep_the_first_place_and_the_last_row() {
        let g = graph(vec![(5, vec![0, 1], false), (2, vec![2], true), (5, vec![3], false)], 2);
        assert_eq!(g.ids, vec![5, 2]);
        assert_eq!(g.members[0], vec![3]);
        assert!(g.places.get(&0).is_none() && g.unique(3));
        assert_eq!((g.tokens(0), g.tokens(1)), (2, 1));
    }

    #[test]
    fn matches_need_equal_tokens_at_the_first_node_of_a_hash() {
        // nodes (k = 2): 0 = [1, 2], 1 = [2, 3], 2 = [9, 9]
        let node_tokens = vec![1, 2, 2, 3, 9, 9];
        let g = graph(vec![(0, vec![0, 1], false), (1, vec![2], false)], 2);
        let keyed = node_index(&g, &node_tokens).unwrap();
        assert_eq!(read_matches(&keyed, &node_tokens, 2, &[7, 1, 2, 3, 9, 9]), vec![(1, 0), (2, 1), (4, 2)]);
        assert!(read_matches(&keyed, &node_tokens, 2, &[1]).is_empty());
        // a collision: two nodes of one hash, only the first (lowest id) is compared
        let fake = vec![(5u64, 0u32), (5, 1)];
        assert!(read_matches(&fake, &node_tokens, 2, &[2, 3]).is_empty());
    }

    #[test]
    fn pieces_follow_offsets_and_close_in_python_dict_order() {
        // unitig 0: nodes 10..15 linear; node 12 also in unitig 1 (a repeat copy)
        let g = graph(vec![(0, vec![10, 11, 12, 13, 14, 15], false), (1, vec![12, 20], false)], 3);
        // windows 0, 1, 2 match nodes 10, 11, 12 (offset 0); window 40 matches node 13 (offset 3 - 40 = -37, a new
        // run); the run of unitig 1 holds only the repeated node 12 and is dropped; open runs close in first-seen order
        let m = vec![(0, 10), (1, 11), (2, 12), (40, 13)];
        assert_eq!(pieces_of_read(&g, &m, 32), vec![(0, 0, 2, 0), (0, 3, 3, 40)]);
        // the same offset after more than max_gap unmatched windows: the old run is closed then, the new one at the end
        let g3 = graph(vec![(0, (0..60).collect(), false)], 3);
        assert_eq!(pieces_of_read(&g3, &[(0, 0), (1, 1), (50, 50)], 32), vec![(0, 0, 1, 0), (0, 50, 50, 50)]);
        // a gap of exactly max_gap unmatched windows keeps the run
        assert_eq!(pieces_of_read(&g3, &[(0, 0), (33, 33)], 32), vec![(0, 0, 33, 0)]);
    }

    #[test]
    fn circular_offsets_wrap() {
        let g = graph(vec![(0, vec![0, 1, 2, 3], true)], 2);
        // windows 0..5 matching nodes 2, 3, 0, 1, 2: one run, b unwrapped past n
        let m = vec![(0, 2), (1, 3), (2, 0), (3, 1), (4, 2)];
        assert_eq!(pieces_of_read(&g, &m, 32), vec![(0, 2, 6, 0)]);
    }

    #[test]
    fn verified_pieces_split_trim_and_need_a_window() {
        // 10 read tokens starting every 10 q-symbols, each 5 long
        let starts: Vec<i64> = (0..10).map(|i| i * 10).collect();
        let ends: Vec<i64> = starts.iter().map(|s| s + 5).collect();
        // piece nodes 0..7 (k = 3) = read tokens 0..9; own starts 10, 20, 30 then 90 (a gap of 60 > 50)
        let own = vec![10, 20, 30, 90];
        assert_eq!(verified_pieces(0, 7, 0, &starts, &ends, &own, 3, 50, 2), vec![(1, 1, 1)]);
        assert_eq!(verified_pieces(0, 7, 0, &starts, &ends, &own, 3, 100, 2), vec![(1, 7, 1)]);
        assert!(verified_pieces(0, 7, 0, &starts, &ends, &[], 3, 50, 2).is_empty());
        assert!(verified_pieces(0, 7, 0, &starts, &ends, &[10, 20], 3, 50, 2).is_empty()); // 2 tokens < k
    }

    #[test]
    fn tiling_takes_the_furthest_piece_then_fallback_then_gaps() {
        // linear, n = 10: verified (0..4), (3..6), fallback (5..9)
        let chosen = tile(&[(0, 4, 0), (3, 6, 1)], 10, false, &[(5, 9, 2)]);
        assert_eq!(chosen, vec![(0, 0, 0, 4), (1, 4, 3, 6), (2, 6, 5, 9)]);
        // a gap: nothing contains node 5
        let chosen = tile(&[(0, 4, 0), (6, 9, 1)], 10, false, &[]);
        assert_eq!(chosen, vec![(0, 0, 0, 4), (1, 6, 6, 9)]);
        assert_eq!(contributions(&chosen, 3, 12), vec![(Some(0), 0, 6), (Some(1), 6, 12)]);
        // a fallback piece hands back at the first verified start on it
        let chosen = tile(&[(0, 2, 0), (6, 9, 1)], 10, false, &[(1, 8, 2)]);
        assert_eq!(chosen, vec![(0, 0, 0, 2), (2, 2, 1, 8), (1, 6, 6, 9)]);
        // a circle of 4: one piece over the origin
        let chosen = tile(&[(2, 5, 0)], 4, true, &[]);
        assert_eq!(chosen, vec![(0, 0, -2, 1), (0, 2, 2, 5)]);
        assert_eq!(contributions(&chosen, 2, 4), vec![(Some(0), 0, 2), (Some(1), 2, 4)]);
        assert!(tile(&[], 3, false, &[]).is_empty());
        assert_eq!(contributions(&[], 2, 5), vec![(None, 0, 5)]);
    }

    #[test]
    fn contributions_keep_python_overlapping_gaps() {
        // a piece whose stop does not pass its start leaves prev behind: the next gap starts at prev again
        let chosen = vec![(0, 3, 3, 3), (1, 3, 3, 3), (2, 6, 6, 8)];
        assert_eq!(contributions(&chosen, 1, 10),
                   vec![(None, 0, 3), (None, 0, 3), (Some(1), 3, 4), (None, 4, 6), (Some(2), 6, 9), (None, 9, 10)]);
    }

    #[test]
    fn galloping_lower_bound_is_a_binary_search() {
        let v: Vec<u64> = vec![1, 3, 3, 3, 7, 9, 9, 20];
        for x in 0..25u64 {
            let want = v.partition_point(|&y| y < x);
            for from in 0..=v.len() {
                assert_eq!(lower_bound_from(&v, x, from), want, "x {} from {}", x, from);
            }
        }
    }

    #[test]
    fn own_tokens_keep_numpys_buffered_ordinals() {
        // reads 0..2 of the uncorrected input, read 1 without tokens: offsets 0, 2, 2. Offset 2 counts once, so read
        // 2's tokens are keyed with ordinal 1.
        let d = std::env::temp_dir().join(format!("squash_layout_own_{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        let gz = |path: &Path, text: &str| {
            let mut e = GzEncoder::new(File::create(path).unwrap(), Compression::default());
            e.write_all(text.as_bytes()).unwrap();
            e.finish().unwrap();
        };
        fs::create_dir_all(d.join("p")).unwrap();
        gz(&d.join("p").join(READS_TABLE), "read_id\tread_ordinal\tsource_symbols\ttoken_offset\ttoken_count\n\
r0\t0\t9\t0\t2\nr1\t1\t9\t2\t0\nr2\t2\t9\t2\t1\n");
        let vals: Vec<u8> = [5u32, 6, 7].iter().flat_map(|v| v.to_le_bytes()).collect();
        fs::write(d.join("p").join(TOKEN_VALUES), vals).unwrap();
        let spans: Vec<u8> = [(0u64, 3u64), (4, 6), (1, 2)].iter().flat_map(|&(a, b)| [a.to_le_bytes(), b.to_le_bytes()].concat()).collect();
        fs::write(d.join("p").join(TOKEN_SPANS), spans).unwrap();
        gz(&d.join("c.tsv.gz"), "read_id\tread_ordinal\tsource_length\ttoken_ordinal\ttoken_id\tstart\tend\tcoordinate_space\n\
r0\t0\t9\t0\t5\t0\t3\tq\nr0\t0\t9\t1\t6\t4\t7\tq\nr1\t1\t9\t0\t7\t1\t2\tq\nr2\t2\t9\t0\t7\t1\t2\tq\n");
        let own = own_token_starts(&d.join("c.tsv.gz"), &d.join("p")).unwrap();
        let mut got: Vec<(String, Vec<i64>)> = own.into_iter().collect();
        got.sort();
        // r0's second token differs in its end; r2's own token is found under ordinal 1 (r1)
        assert_eq!(got, vec![("r0".to_string(), vec![0]), ("r1".to_string(), vec![1])]);
        let _ = fs::remove_dir_all(&d);
    }
}
