//! Contig joins after read threading (`threading --join`): a dead end A (no successor) of a long unitig is joined
//! to the start of another long unitig B when the reads say so, by two rules applied in order.
//!
//! - D, signal-guided read hopping: the reads hanging off A's end share >= 3 six-token words on one diagonal with A's
//!   last 100 tokens; of the 10 that run furthest past it, the one most others align with on their uncorrected
//!   q-symbols (3,000 ending 200 before A's end) starts. Each hop takes, among the current read's same-strand partners
//!   of the first partner search that run >= 500 q-symbols past its end (the 15 that run furthest), the one reaching
//!   furthest whose q-symbols align with the current read's last 3,000 (200 before its end left out): a weighted edit
//!   distance (|a - b| per substitution, gap 2) in a band of 400 that follows the path from the window's middle row,
//!   placed there by a subsequence alignment within +-1,000 of the partners' diagonal, at a cost <= `threshold` per
//!   alignment column. The walk stops when a read's tokens past the previous read hold >= 3 words unique to another
//!   unitig: a join when that unitig has >= `min_nodes` nodes and the first word lies in its first `head` tokens. A
//!   start reached from two ends, or one that has a predecessor, is joined to neither.
//!   With `--join-min-paths N` (N >= 2) the walk is repeated from the next start reads, each walk barred from the
//!   reads of the walks before it, and a join needs N walks landing on the same unitig (so a single fold-back read
//!   cannot make a join); the first of them is the join's path. N = 1 is the single walk.
//! - A, word bridges, on the ends D left: >= 3 reads that touch A's last 400 tokens and then B's first 400 tokens by
//!   >= 3 words each, words unique to one end region; each the other's best partner by >= 3x the next.
//!
//! Both need A and B of >= `min_nodes` nodes. A join becomes one unitig: A up to the last word it shares with the first
//! read, each read's tokens between the word it shares with the previous piece and the word it shares with the next
//! (shared words on the pair's dominant diagonal; two reads sharing no word are spliced at the q-symbol where the
//! hop's verified window ends), B from the first word it shares with the last read; the k-token windows across the
//! junctions become new nodes, appended to the node tokens. D's joins are applied first, then A's; a join
//! that would close a circle is skipped. The rules read no reference.

use crate::io::{READS_TABLE, Result, TOKEN_SPANS, tsv_rows};
use crate::threading::Graph;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::convert::TryInto;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::os::unix::fs::FileExt;
use std::path::Path;

const KW: usize = 6; // word length (tokens)
const MAXOCC: usize = 200; // words in more places are not indexed
const TAIL: usize = 100; // recruitment tail (tokens)
const MINHIT: usize = 3; // hits needed: words on one diagonal, words per end and reads per pair (rule A)
const DIAG: i64 = 8; // tokens a shared word may lie off the median diagonal
const MINOV: usize = 8; // tokens a recruited read must run past the tail
const START_READS: usize = 10; // reads hanging off a dead end that may start a walk
const WIN: i64 = 3000; // q-symbols of a hop's alignment window
const MARGIN: i64 = 200; // q-symbols before a read's end left out of the window
const MAXCAND: usize = 15; // hop candidates aligned (furthest reaching first)
const MAXHOPS: usize = 150; // hops before a walk stalls
const MIN_EXT: i64 = 500; // q-symbols a hop candidate must run past the current read's end
const BRIDGE_W: usize = 400; // end regions of rule A and of the splice (tokens)

pub struct Params {
    /// q-symbol alignment cost per column at or below which two reads align
    pub threshold: f64,
    /// nodes both unitigs of a join need
    pub min_nodes: usize,
    /// a walk must land within this many tokens of the unitig's start
    pub head: usize,
    /// independent walks (disjoint reads) that must land on the same unitig; 1 = the single walk
    pub min_paths: usize,
}

/// The reads (the graph's input: corrected tokens), their q-symbol positions, uncorrected q-symbols and partners (the
/// first partner search's).
pub struct Reads {
    tok: Vec<u32>,
    ranges: Vec<(usize, usize)>,
    q: Vec<u32>,
    pub ids: Vec<String>,
    qlen: Vec<i64>,
    sym: SymFile,
    partners: Vec<Vec<(u32, i64)>>,
}

impl Reads {
    fn read(&self, i: usize) -> &[u32] {
        &self.tok[self.ranges[i].0..self.ranges[i].1]
    }

    fn len(&self, i: usize) -> usize {
        self.ranges[i].1 - self.ranges[i].0
    }

    /// q-symbol start of token t of read i (clamped into the read).
    fn qat(&self, i: usize, t: i64) -> i64 {
        let n = self.len(i) as i64;
        self.q[self.ranges[i].0 + t.clamp(0, (n - 1).max(0)) as usize] as i64
    }

    /// The first token of read i at or after q-symbol q.
    fn tok_at(&self, i: usize, q: i64) -> usize {
        self.q[self.ranges[i].0..self.ranges[i].1].partition_point(|&x| (x as i64) < q)
    }
}

/// tokens + ranges as threading::load_reads gives them; the rest from the graph's input reads, the symbols table
/// and the partner table.
pub fn load(reads_dir: &Path, tokens: Vec<u32>, ranges: Vec<(usize, usize)>, symbols: &Path, partners: &Path) -> Result<Reads> {
    let raw = std::fs::read(reads_dir.join(TOKEN_SPANS))?;
    if raw.len() != tokens.len() * 16 {
        return Err("the input token spans do not match the tokens".into());
    }
    let q: Vec<u32> = raw.chunks_exact(16).map(|b| u64::from_le_bytes(b[..8].try_into().unwrap()) as u32).collect();
    let (mut ids, mut qlen) = (Vec::new(), Vec::new());
    tsv_rows(&reads_dir.join(READS_TABLE), &["read_id", "source_symbols"], |r| {
        ids.push(r[0].to_string());
        qlen.push(r[1].parse()?);
        Ok(())
    })?;
    if ids.len() != ranges.len() {
        return Err("reads table and token ranges differ".into());
    }
    let idx: HashMap<&str, usize> = ids.iter().enumerate().map(|(i, r)| (r.as_str(), i)).collect();
    let sym = SymFile::open(symbols, &idx, ids.len())?;
    let mut part = vec![Vec::new(); ids.len()];
    tsv_rows(partners, &["target_id", "donor_id", "orientation", "target_start", "donor_start"], |r| {
        let (Some(&t), Some(&d)) = (idx.get(r[0]), idx.get(r[1])) else { return Ok(()) };
        if r[2] != "+" {
            return Ok(());
        }
        let o = r[3].parse::<i64>()? - r[4].parse::<i64>()?;
        part[t].push((d as u32, o));
        part[d].push((t as u32, -o));
        Ok(())
    })?;
    Ok(Reads { tok: tokens, ranges, q, ids, qlen, sym, partners: part })
}

/// The reads' uncorrected q-symbols, left in the symbols table and read when a walk needs them (parsing every read
/// into memory costs gigabytes on a large genome). One pass records where each read's symbols start and how many
/// there are; a line of single-digit values ("5,4,6,...") is read in place (value k at byte 2k), any other line is
/// parsed whole into `other` at that pass.
pub struct SymFile {
    file: File,
    /// per read: (byte offset of its first symbol, symbol count); count 0 = not in the table
    at: Vec<(u64, usize)>,
    other: HashMap<usize, Vec<u8>>,
}

impl SymFile {
    fn open(path: &Path, idx: &HashMap<&str, usize>, n: usize) -> Result<SymFile> {
        let file = File::open(path).map_err(|e| format!("{}: {}", path.display(), e))?;
        let mut rd = BufReader::with_capacity(1 << 24, File::open(path)?);
        let (mut at, mut other) = (vec![(0u64, 0usize); n], HashMap::new());
        let (mut line, mut off) = (Vec::new(), 0u64);
        loop {
            line.clear();
            let got = rd.read_until(b'\n', &mut line)?;
            if got == 0 {
                break;
            }
            let start = off;
            off += got as u64;
            let body = if line.last() == Some(&b'\n') { &line[..line.len() - 1] } else { &line[..] };
            if body.first() == Some(&b'#') {
                continue;
            }
            let tab = body.iter().position(|&b| b == b'\t').ok_or("symbols: no tab")?;
            let rid = std::str::from_utf8(&body[..tab]).map_err(|_| "symbols: read id not UTF-8")?;
            let Some(&i) = idx.get(rid) else { continue };
            let rest = &body[tab + 1..];
            let fixed = rest.len() % 2 == 1
                && rest.iter().enumerate().all(|(k, &b)| if k % 2 == 0 { b.is_ascii_digit() } else { b == b',' });
            if fixed {
                at[i] = (start + tab as u64 + 1, (rest.len() + 1) / 2);
            } else {
                let text = std::str::from_utf8(rest).map_err(|_| "symbols: not UTF-8")?;
                let v: Vec<u8> = if text.is_empty() { Vec::new() } else {
                    text.split(',').map(|v| v.trim().parse::<u8>()).collect::<std::result::Result<_, _>>()?
                };
                at[i] = (0, v.len());
                other.insert(i, v);
            }
        }
        Ok(SymFile { file, at, other })
    }

    /// Read i's number of q-symbols.
    fn len(&self, i: usize) -> usize {
        self.at[i].1
    }

    /// Read i's q-symbols lo..hi as f32.
    fn values(&self, i: usize, lo: usize, hi: usize) -> Vec<f32> {
        if let Some(v) = self.other.get(&i) {
            return v[lo..hi].iter().map(|&x| x as f32).collect();
        }
        let mut buf = vec![0u8; 2 * (hi - lo) - 1];
        self.file.read_exact_at(&mut buf, self.at[i].0 + 2 * lo as u64).expect("reading the symbols table");
        buf.iter().step_by(2).map(|&b| (b - b'0') as f32).collect()
    }
}

fn word_hashes(t: &[u32]) -> Vec<u64> {
    if t.len() < KW {
        return Vec::new();
    }
    (0..=t.len() - KW).map(|i| crate::threading::window_hash(&t[i..i + KW])).collect()
}

/// The median of integers; of an even count, the mean of the middle two truncated toward zero.
fn py_median(v: &mut [i64]) -> i64 {
    v.sort_unstable();
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        ((v[n / 2 - 1] + v[n / 2]) as f64 / 2.0).trunc() as i64
    }
}

/// Every 6-token word of every read: (hash, read, position) sorted by hash; words in more than MAXOCC places dropped.
struct WordIndex {
    h: Vec<u64>,
    r: Vec<u32>,
    p: Vec<u32>,
}

impl WordIndex {
    fn new(reads: &Reads) -> Self {
        let per: Vec<Vec<(u64, u32, u32)>> = (0..reads.ranges.len()).into_par_iter()
            .map(|i| word_hashes(reads.read(i)).into_iter().enumerate().map(|(p, h)| (h, i as u32, p as u32)).collect())
            .collect();
        let mut all: Vec<(u64, u32, u32)> = per.into_iter().flatten().collect();
        all.par_sort_by_key(|x| x.0);
        let (mut h, mut r, mut p) = (Vec::new(), Vec::new(), Vec::new());
        let mut i = 0;
        while i < all.len() {
            let mut j = i;
            while j < all.len() && all[j].0 == all[i].0 {
                j += 1;
            }
            if j - i <= MAXOCC {
                for x in &all[i..j] {
                    h.push(x.0);
                    r.push(x.1);
                    p.push(x.2);
                }
            }
            i = j;
        }
        WordIndex { h, r, p }
    }

    /// Reads aligned to the tail and running past its end: [(read, start of the overhang in the read)], by read.
    fn recruit(&self, reads: &Reads, tail: &[u32]) -> Vec<(usize, usize)> {
        let mut hits: Vec<(u32, i64)> = Vec::new();
        for (j, h) in word_hashes(tail).into_iter().enumerate() {
            let lo = self.h.partition_point(|&x| x < h);
            let hi = self.h.partition_point(|&x| x <= h);
            for k in lo..hi {
                hits.push((self.r[k], self.p[k] as i64 - j as i64));
            }
        }
        hits.sort_by_key(|x| x.0);
        let mut out = Vec::new();
        let mut a = 0;
        while a < hits.len() {
            let mut b = a;
            while b < hits.len() && hits[b].0 == hits[a].0 {
                b += 1;
            }
            if b - a >= MINHIT {
                let mut offs: Vec<i64> = hits[a..b].iter().map(|x| x.1).collect();
                let med = py_median(&mut offs);
                if offs.iter().filter(|&&o| (o - med).abs() <= DIAG).count() >= MINHIT {
                    let start = med + tail.len() as i64;
                    let rid = hits[a].0 as usize;
                    if start > 0 && reads.len(rid) as i64 - start >= MINOV as i64 {
                        out.push((rid, start as usize));
                    }
                }
            }
            a = b;
        }
        out
    }
}

/// q-symbol alignment: a weighted edit distance (substitution |a - b|, gap 2) in a band of 400 whose centre follows the
/// previous row's best cell (at most 50 off the diagonal step per row), seeded by a subsequence alignment.
mod qalign {
    const GAP: f32 = 2.0;

    #[inline]
    fn cst(a: f32, b: f32) -> f32 {
        (a - b).abs()
    }

    #[inline]
    fn step(up: f32, dg: f32, left: Option<f32>, c: f32) -> (f32, u8) {
        let (mut b, mut d) = (dg + c, 1u8);
        if up + GAP < b {
            b = up + GAP;
            d = 0;
        }
        if let Some(l) = left {
            if l + GAP < b {
                b = l + GAP;
                d = 2;
            }
        }
        (b, d)
    }

    #[allow(clippy::too_many_arguments)]
    fn tracked(x: &[f32], y: &[f32], c0: i64, anchor: Option<i64>, band: i64, shift: i64, dir: &mut Vec<u8>,
               lefts: &mut Vec<i64>, path: &mut Vec<(usize, i64)>) -> Option<f32> {
        const INF: f32 = f32::INFINITY;
        let n = x.len();
        let m = y.len() as i64;
        let width = (2 * band + 1) as usize;
        dir.clear();
        dir.resize(n * width, 255);
        lefts.clear();
        let mut pc = vec![INF; width];
        let mut cc = vec![INF; width];
        let mut centre = anchor.unwrap_or(c0);
        let mut pl = centre - band;
        lefts.push(pl);
        let (mut best_k, mut best_v) = (band as usize, INF);
        for k in 0..width {
            let j = pl + k as i64;
            if j >= 0 && j < m && anchor.map_or(true, |a| a == j) {
                pc[k] = cst(x[0], y[j as usize]);
                dir[k] = 3;
                if pc[k] < best_v {
                    best_v = pc[k];
                    best_k = k;
                }
            }
        }
        if best_v.is_infinite() {
            return None;
        }
        for i in 1..n {
            let want = pl + best_k as i64 + 1;
            centre = want.clamp(centre + 1 - shift, centre + 1 + shift);
            let l = centre - band;
            lefts.push(l);
            let xi = x[i];
            best_v = INF;
            for k in 0..width {
                let j = l + k as i64;
                if j < 0 || j >= m {
                    cc[k] = INF;
                    continue;
                }
                let prev = |jj: i64| -> f32 {
                    let kk = jj - pl;
                    if kk >= 0 && (kk as usize) < width { pc[kk as usize] } else { INF }
                };
                let (b, dd) = step(prev(j), prev(j - 1), if k > 0 { Some(cc[k - 1]) } else { None }, cst(xi, y[j as usize]));
                if b.is_infinite() {
                    cc[k] = INF;
                    continue;
                }
                cc[k] = b;
                dir[i * width + k] = dd;
                if cc[k] < best_v {
                    best_v = cc[k];
                    best_k = k;
                }
            }
            if best_v.is_infinite() {
                return None;
            }
            std::mem::swap(&mut pc, &mut cc);
            pl = l;
        }
        path.clear();
        let (mut i, mut j) = (n - 1, pl + best_k as i64);
        loop {
            path.push((i, j));
            let k = (j - lefts[i]) as usize;
            match dir[i * width + k] {
                3 => break,
                0 => i -= 1,
                1 => {
                    i -= 1;
                    j -= 1;
                }
                2 => j -= 1,
                _ => panic!("bad traceback"),
            }
        }
        Some(best_v)
    }

    /// Semi-global subsequence alignment (free start and end in y): (cost, start column of the best end's path).
    fn subseq(x: &[f32], y: &[f32]) -> (f32, usize) {
        let m = y.len();
        let mut pc: Vec<f32> = y.iter().map(|&v| cst(x[0], v)).collect();
        let mut ps: Vec<u32> = (0..m as u32).collect();
        let mut cc = vec![0f32; m];
        let mut cs = vec![0u32; m];
        for (r, &xi) in x[1..].iter().enumerate() {
            let fresh = (r + 1) as f32 * GAP + cst(xi, y[0]);
            if pc[0] + GAP <= fresh {
                cc[0] = pc[0] + GAP;
                cs[0] = ps[0];
            } else {
                cc[0] = fresh;
                cs[0] = 0;
            }
            for j in 1..m {
                let (b, d) = step(pc[j], pc[j - 1], Some(cc[j - 1]), cst(xi, y[j]));
                cc[j] = b;
                cs[j] = match d {
                    0 => ps[j],
                    1 => ps[j - 1],
                    _ => cs[j - 1],
                };
            }
            std::mem::swap(&mut pc, &mut cc);
            std::mem::swap(&mut ps, &mut cs);
        }
        let (mut bj, mut bv) = (0usize, f32::INFINITY);
        for (j, &v) in pc.iter().enumerate() {
            if v < bv {
                bv = v;
                bj = j;
            }
        }
        (bv, ps[bj] as usize)
    }

    #[allow(clippy::too_many_arguments)]
    fn from_middle(x: &[f32], y: &[f32], line_c: &[i64], mid: usize, band: i64, shift: i64, seed: i64,
                   dir: &mut Vec<u8>, lefts: &mut Vec<i64>, path: &mut Vec<(usize, i64)>) -> Option<f32> {
        let n = x.len();
        let anchor = if seed > 0 {
            let w = 1000.min(n - mid);
            let m = y.len() as i64;
            let (l, h) = ((line_c[mid] - seed).max(0), (line_c[mid] + w as i64 + seed).min(m));
            if h - l < w as i64 / 2 {
                None
            } else {
                let (_, st) = subseq(&x[mid..mid + w], &y[l as usize..h as usize]);
                Some(l + st as i64)
            }
        } else {
            None
        };
        let total_f = tracked(&x[mid..], y, line_c[mid], anchor, band, shift, dir, lefts, path)?;
        let mut cells: Vec<(usize, i64)> = path.iter().map(|&(i, j)| (i + mid, j)).collect();
        let j_mid = cells.iter().filter(|c| c.0 == mid).map(|c| c.1).min().unwrap();
        let mut total = total_f;
        if mid > 0 {
            let xr: Vec<f32> = x[..=mid].iter().rev().copied().collect();
            let yr: Vec<f32> = y.iter().rev().copied().collect();
            let m = y.len() as i64;
            let a = m - 1 - j_mid;
            let tb = tracked(&xr, &yr, a, Some(a), band, shift, dir, lefts, path)?;
            total += tb - cst(x[mid], y[j_mid as usize]);
            for &(i, j) in path.iter() {
                if i > 0 {
                    cells.push((mid - i, m - 1 - j));
                }
            }
        }
        path.clear();
        path.extend(cells);
        Some(total)
    }

    /// Cost per alignment column of x (positions lo..lo + x.len() of its read) against the whole of y, the band started
    /// along `line` (two points (x position, y position); slope clamped to 0.9..1.1) and run from x's middle both ways:
    /// the path's cost, rounded to 3 decimals, over its cells. None when no path.
    pub fn cost(x: &[f32], y: &[f32], lo: i64, line: (i64, i64, i64, i64)) -> Option<f64> {
        let (band, shift, seed) = (400i64, 50i64, 1000i64);
        let n = x.len();
        let m = y.len() as i64;
        if n == 0 || m == 0 {
            return None;
        }
        let (a0, b0, a1, b1) = line;
        let slope = if a1 > a0 { ((b1 - b0) as f64 / (a1 - a0) as f64).clamp(0.9, 1.1) } else { 1.0 };
        let pf = |i: f64| -> f64 { b0 as f64 + (lo as f64 + i - a0 as f64) * slope };
        let line_c: Vec<i64> = (0..n).map(|i| (pf(i as f64).round() as i64).clamp(0, m - 1)).collect();
        let (mut dir, mut lefts, mut path) = (Vec::new(), Vec::new(), Vec::new());
        let mid = (((a0 + a1) / 2 - lo).clamp(0, n as i64 - 1)) as usize;
        from_middle(x, y, &line_c, mid, band, shift, seed, &mut dir, &mut lefts, &mut path)?;
        let mut rowcost = vec![0f64; n];
        let mut fw: Vec<(usize, i64)> = path.iter().copied().filter(|c| c.0 >= mid).collect();
        fw.sort_unstable();
        let mut bw: Vec<(usize, i64)> = path.iter().copied().filter(|c| c.0 <= mid).collect();
        bw.sort_unstable_by(|a, b| b.cmp(a));
        for (part, first) in [(&fw, true), (&bw, false)] {
            let mut last: Option<(usize, i64)> = None;
            for &(i, j) in part.iter() {
                let c = match last {
                    Some((pi, pj)) if (i as i64 - pi as i64).abs() == 1 && (j - pj).abs() == 1 => cst(x[i], y[j as usize]) as f64,
                    Some(_) => GAP as f64,
                    None => if first { cst(x[i], y[j as usize]) as f64 } else { 0.0 },
                };
                rowcost[i] += c;
                last = Some((i, j));
            }
        }
        let total: f64 = rowcost.iter().sum();
        let cells = path.len();
        if cells == 0 {
            return None;
        }
        let printed: f64 = format!("{:.3}", total as f32).parse().ok()?;
        Some(printed / cells as f64)
    }
}

fn spell(g: &Graph, node_tokens: &[u32], k: usize, u: i64) -> Vec<u32> {
    let m = &g.units[&u].members;
    let mut t: Vec<u32> = node_tokens[m[0] as usize * k..(m[0] as usize + 1) * k].to_vec();
    for &x in &m[1..] {
        t.push(node_tokens[(x as usize + 1) * k - 1]);
    }
    t
}

/// One proposed join: dead end `a` to the start of `b` through `reads` (in order).
#[derive(Clone, Debug)]
pub struct Join {
    pub rule: &'static str,
    pub a: i64,
    pub b: i64,
    pub reads: Vec<usize>,
    /// per hop: (q-symbol where the current read's verified window ends, offset q_current - q_next)
    pub hop_q: Vec<(i64, i64)>,
    pub hops: usize,
    pub note: String,
}

/// One applied join, for joins.tsv.
#[derive(Clone, Debug)]
pub struct Applied {
    pub join: Join,
    pub a_last_node: u32,
    pub b_first_node: u32,
    pub merged: i64,
    pub dropped_a: usize,
    pub dropped_b: usize,
    pub gap_tokens: usize,
    pub new_nodes: usize,
    pub status: String,
}

struct Ctx<'a> {
    g: &'a Graph,
    reads: &'a Reads,
    index: &'a WordIndex,
    body: HashMap<u64, (i64, usize)>, // words unique to one unitig: hash -> (unitig, position)
    toks: HashMap<i64, Vec<u32>>,     // each unitig's tokens
    p: &'a Params,
}

impl<'a> Ctx<'a> {
    /// Cost per column of read i's 3,000 q-symbols ending MARGIN before q-position e against read j (q_j = q_i - off).
    fn window_cost(&self, i: usize, e: i64, j: usize, off: i64) -> Option<f64> {
        let (lo, hi) = (e - MARGIN - WIN, e - MARGIN);
        if lo < 0 || lo - off < 0 {
            return None;
        }
        let (xn, yn) = (self.reads.sym.len(i), self.reads.sym.len(j));
        let hi = hi.min(xn as i64);
        if xn == 0 || yn == 0 || hi <= lo {
            return None;
        }
        let x = self.reads.sym.values(i, lo as usize, hi as usize);
        let y = self.reads.sym.values(j, 0, yn);
        qalign::cost(&x, &y, lo, (lo + 1000, lo + 1000 - off, lo + 2000, lo + 2000 - off))
    }

    /// q-symbol offset of read j against read i (q_i - q_j) from their shared corrected tokens in the 30 before si / sj
    /// (the median over >= 3 shared values).
    fn offset(&self, i: usize, si: usize, j: usize, sj: usize) -> Option<i64> {
        let r = self.reads;
        let dict = |x: usize, s: usize| -> HashMap<u32, i64> {
            let t = r.read(x);
            let from = s.saturating_sub(30);
            (from..s.min(t.len())).map(|p| (t[p], r.qat(x, p as i64))).collect()
        };
        let (ti, tj) = (dict(i, si), dict(j, sj));
        let mut d: Vec<i64> = ti.iter().filter_map(|(v, qi)| tj.get(v).map(|qj| qi - qj)).collect();
        if d.len() >= 3 { Some(py_median(&mut d)) } else { None }
    }

    /// >= 3 words of the read's tokens from `from` on unique to one unitig other than a: (unitig, its first word's
    /// position in that unitig).
    fn landing(&self, i: usize, from: usize, a: i64) -> Option<(i64, usize)> {
        let t = self.reads.read(i);
        let from = from.min(t.len());
        let mut cnt: HashMap<i64, (usize, usize)> = HashMap::new();
        for h in word_hashes(&t[from..]) {
            if let Some(&(u, pos)) = self.body.get(&h) {
                if u != a {
                    let e = cnt.entry(u).or_insert((0, pos));
                    e.0 += 1;
                    if e.0 >= 3 {
                        return Some((u, e.1));
                    }
                }
            }
        }
        None
    }

    /// The walk's landing unitig when it may be joined: landed within the head of a unitig of >= min_nodes nodes.
    fn accept(&self, w: &Walk) -> Option<i64> {
        let (u, pos) = w.end?;
        if self.g.units[&u].members.len() < self.p.min_nodes || pos > self.p.head {
            return None;
        }
        Some(u)
    }

    /// Rule D from dead end a.
    fn hop(&self, a: i64) -> Option<Join> {
        let r = self.reads;
        let ta = &self.toks[&a];
        let mut ov = self.index.recruit(r, &ta[ta.len().saturating_sub(TAIL)..]);
        ov.sort_by_key(|&(i, s)| std::cmp::Reverse(r.len(i) as i64 - s as i64));
        ov.truncate(START_READS);
        if ov.is_empty() {
            return None;
        }
        let mut start = ov[0];
        let mut order: Vec<usize> = (0..ov.len()).collect();
        if ov.len() > 1 {
            let mut support = vec![0usize; ov.len()];
            for x in 0..ov.len() {
                for y in 0..ov.len() {
                    if x == y {
                        continue;
                    }
                    let ((i, si), (j, sj)) = (ov[x], ov[y]);
                    if let Some(off) = self.offset(i, si, j, sj) {
                        if let Some(c) = self.window_cost(i, r.qat(i, si as i64), j, off) {
                            if c <= self.p.threshold {
                                support[x] += 1;
                            }
                        }
                    }
                }
            }
            let mut best = 0;
            for x in 1..ov.len() {
                let key = |z: usize| (support[z], r.len(ov[z].0) as i64 - ov[z].1 as i64);
                if key(x) > key(best) {
                    best = x;
                }
            }
            start = ov[best];
            // the other start reads for further walks: by the same key, ties by position (best first)
            order.sort_by_key(|&z| std::cmp::Reverse((support[z], r.len(ov[z].0) as i64 - ov[z].1 as i64)));
            order.retain(|&z| z != best);
            order.insert(0, best);
        }
        if self.p.min_paths <= 1 {
            let w = self.walk(a, start, &HashSet::new());
            let u = self.accept(&w)?;
            return Some(Join { rule: "D hop", a, b: u, reads: w.path, hop_q: w.hop_q, hops: w.hops, note: String::new() });
        }
        // independent walks from the next start reads, each barred from the reads of the walks before it
        let mut barred: HashSet<usize> = HashSet::new();
        let mut walks: Vec<(Walk, Option<i64>)> = Vec::new();
        let mut votes: HashMap<i64, usize> = HashMap::new();
        for &z in &order {
            if barred.contains(&ov[z].0) {
                continue;
            }
            let w = self.walk(a, ov[z], &barred);
            barred.extend(w.path.iter().copied());
            let ok = self.accept(&w);
            walks.push((w, ok));
            if let Some(u) = ok {
                let n = votes.entry(u).or_insert(0);
                *n += 1;
                if *n >= self.p.min_paths {
                    let (w, _) = walks.into_iter().find(|(_, x)| *x == Some(u)).unwrap();
                    let note = format!("{} paths", self.p.min_paths);
                    return Some(Join { rule: "D hop", a, b: u, reads: w.path, hop_q: w.hop_q, hops: w.hops, note });
                }
            }
        }
        None
    }

    /// One walk from `start`, never through a barred read: the landing (unitig, word position), or None when it stalls.
    fn walk(&self, a: i64, start: (usize, usize), barred: &HashSet<usize>) -> Walk {
        let r = self.reads;
        let (mut cur, mut land_from) = start;
        let mut visited: HashSet<usize> = std::iter::once(cur).collect();
        let mut path = vec![cur];
        let mut hop_q = Vec::new();
        let mut hops = 0;
        let end = loop {
            if hops >= MAXHOPS {
                break None;
            }
            if let Some(landed) = self.landing(cur, land_from, a) {
                break Some(landed);
            }
            let e_cur = r.qat(cur, r.len(cur) as i64 - 1);
            let mut cands: Vec<(usize, usize, i64)> = r.partners[cur].iter()
                .filter(|&&(c, off)| !visited.contains(&(c as usize)) && !barred.contains(&(c as usize))
                        && r.qlen[c as usize] - (e_cur - off) >= MIN_EXT)
                .map(|&(c, off)| (c as usize, r.tok_at(c as usize, e_cur - off), off)).collect();
            cands.sort_by_key(|&(c, _, off)| std::cmp::Reverse(r.qlen[c] - (e_cur - off)));
            cands.truncate(MAXCAND);
            if cands.is_empty() {
                break None;
            }
            let mut chosen: Option<(usize, usize, i64)> = None;
            for &(c, sc, off) in &cands {
                if let Some(cost) = self.window_cost(cur, e_cur, c, off) {
                    if cost <= self.p.threshold
                        && chosen.map_or(true, |(b, _, bo)| r.qlen[c] - (e_cur - off) > r.qlen[b] - (e_cur - bo)) {
                        chosen = Some((c, sc, off));
                    }
                }
            }
            let Some((c, sc, off)) = chosen else {
                break None;
            };
            hop_q.push((e_cur - MARGIN, off));
            visited.insert(c);
            path.push(c);
            hops += 1;
            cur = c;
            land_from = sc;
        };
        Walk { path, hop_q, hops, end }
    }
}

/// One walk of rule D: its reads, hop points, hop count and landing (unitig, word position; None when it stalled).
struct Walk {
    path: Vec<usize>,
    hop_q: Vec<(i64, i64)>,
    hops: usize,
    end: Option<(i64, usize)>,
}

/// Rule A on the given tails / heads: (a, b, best bridging read) for every accepted pair.
fn bridges(ctx: &Ctx, tails: &[i64], heads: &[i64]) -> Vec<(i64, i64, usize, usize)> {
    let mut owner: [HashMap<u64, Option<i64>>; 2] = [HashMap::new(), HashMap::new()];
    for (side, units) in [(0usize, tails), (1usize, heads)] {
        for &u in units {
            let t = &ctx.toks[&u];
            let region = if side == 0 { &t[t.len().saturating_sub(BRIDGE_W)..] } else { &t[..t.len().min(BRIDGE_W)] };
            let set: HashSet<u64> = word_hashes(region).into_iter().collect();
            for h in set {
                owner[side].entry(h).and_modify(|o| *o = None).or_insert(Some(u));
            }
        }
    }
    let n = ctx.reads.ranges.len();
    // per read: the (a, b) pairs it bridges, with its hit count
    let per: Vec<Vec<(i64, i64, usize)>> = (0..n).into_par_iter().map(|i| {
        let hs = word_hashes(ctx.reads.read(i));
        if hs.len() < 2 * MINHIT {
            return Vec::new();
        }
        let (mut tn, mut tl, mut hn, mut hf): (HashMap<i64, usize>, HashMap<i64, usize>, HashMap<i64, usize>, HashMap<i64, usize>) =
            Default::default();
        for (p, h) in hs.iter().enumerate() {
            if let Some(Some(u)) = owner[0].get(h) {
                *tn.entry(*u).or_insert(0) += 1;
                tl.insert(*u, p);
            }
            if let Some(Some(u)) = owner[1].get(h) {
                *hn.entry(*u).or_insert(0) += 1;
                hf.entry(*u).or_insert(p);
            }
        }
        let mut out = Vec::new();
        for (&a, &na) in &tn {
            if na < MINHIT {
                continue;
            }
            for (&b, &nb) in &hn {
                if b != a && nb >= MINHIT && hf[&b] + BRIDGE_W > tl[&a] {
                    out.push((a, b, na + nb));
                }
            }
        }
        out
    }).collect();
    let mut pairs: HashMap<(i64, i64), Vec<(usize, usize)>> = HashMap::new();
    for (i, v) in per.into_iter().enumerate() {
        for (a, b, hits) in v {
            pairs.entry((a, b)).or_default().push((i, hits));
        }
    }
    let (mut by_a, mut by_b): (HashMap<i64, Vec<(i64, usize)>>, HashMap<i64, Vec<(i64, usize)>>) = Default::default();
    for (&(a, b), rs) in &pairs {
        by_a.entry(a).or_default().push((b, rs.len()));
        by_b.entry(b).or_default().push((a, rs.len()));
    }
    let top2 = |v: &Vec<(i64, usize)>| -> (i64, usize, usize) {
        let mut w = v.clone();
        w.sort_by_key(|&(x, n)| (std::cmp::Reverse(n), x));
        (w[0].0, w[0].1, w.get(1).map_or(0, |x| x.1))
    };
    let mut out = Vec::new();
    let mut a_ids: Vec<i64> = by_a.keys().copied().collect();
    a_ids.sort_unstable();
    for a in a_ids {
        let (b, n, second) = top2(&by_a[&a]);
        if n < MINHIT || n < 3 * second {
            continue;
        }
        let (back, _, back_second) = top2(&by_b[&b]);
        if back != a || n < 3 * back_second {
            continue;
        }
        if ctx.g.units[&a].members.len() < ctx.p.min_nodes || ctx.g.units[&b].members.len() < ctx.p.min_nodes {
            continue;
        }
        let rs = &pairs[&(a, b)];
        let &(read, _) = rs.iter().max_by_key(|&&(i, h)| (h, std::cmp::Reverse(i))).unwrap();
        out.push((a, b, read, n));
    }
    out
}

/// The splice of a join: (tokens of X kept, the tokens between, Y's entry token), or a reason it fails. `fallback[k]`
/// (consecutive reads k, k + 1): the token positions of one q-symbol in both, used when they share no word (their
/// corrected tokens can differ where they were verified on the q-symbols).
fn splice(x: &[u32], y: &[u32], reads: &[&[u32]], fallback: &[(usize, usize)]) -> std::result::Result<(usize, Vec<u32>, usize), String> {
    // pieces: X (its last 400 tokens may be left), the reads, Y (entered within its first 400 tokens). Piece k keeps
    // tokens [start[k], exit[k] + KW): the path holds every token before start[k], and its next word starts at >= lo[k].
    let mut pieces: Vec<&[u32]> = vec![x];
    pieces.extend_from_slice(reads);
    pieces.push(y);
    let last = pieces.len() - 1;
    let mut start = vec![0usize; pieces.len()];
    let mut lo = vec![0usize; pieces.len()];
    let mut exit = vec![0usize; pieces.len()];
    lo[0] = x.len().saturating_sub(BRIDGE_W);
    for k in 0..last {
        let (p, q) = (pieces[k], pieces[k + 1]);
        if p.len() < KW || q.len() < KW {
            return Err(format!("piece {} or {} shorter than a word", k, k + 1));
        }
        let q_hi = if k + 1 == last { q.len().min(BRIDGE_W) } else { q.len() };
        let mut at: HashMap<&[u32], Vec<usize>> = HashMap::new();
        for i in lo[k]..=p.len() - KW {
            at.entry(&p[i..i + KW]).or_default().push(i);
        }
        let mut m: Vec<(usize, usize)> = Vec::new();
        for j in 0..q_hi.saturating_sub(KW - 1) {
            if let Some(v) = at.get(&q[j..j + KW]) {
                m.extend(v.iter().map(|&i| (i, j)));
            }
        }
        if m.is_empty() {
            // two reads of a hop: splice at the q-symbol where the current read's verified window ends
            if k >= 1 && k + 1 < last {
                let (tp, tq) = fallback[k - 1];
                if tp >= KW && tp - KW >= lo[k] && tp <= p.len() && tq + KW <= q.len() {
                    exit[k] = tp - KW;
                    start[k + 1] = tq;
                    lo[k + 1] = tq;
                    continue;
                }
            }
            return Err(format!("no shared word between pieces {} and {}", k, k + 1));
        }
        let mut d: Vec<i64> = m.iter().map(|&(i, j)| j as i64 - i as i64).collect();
        let med = py_median(&mut d);
        m.retain(|&(i, j)| (j as i64 - i as i64 - med).abs() <= DIAG);
        if m.is_empty() {
            // the median of an even count can fall between two diagonals that are both far from it
            return Err(format!("no shared word near the median diagonal between pieces {} and {}", k, k + 1));
        }
        let &(i, j) = if k + 1 == last {
            m.iter().min_by_key(|&&(i, j)| (j, std::cmp::Reverse(i))).unwrap()
        } else {
            m.iter().max_by_key(|&&(i, j)| (i, std::cmp::Reverse(j))).unwrap()
        };
        exit[k] = i;
        start[k + 1] = j + KW;
        lo[k + 1] = j;
    }
    let mut mids = Vec::new();
    for k in 1..last {
        if exit[k] + KW < start[k] {
            return Err(format!("read {} leaves before it is entered", k));
        }
        mids.extend_from_slice(&pieces[k][start[k]..exit[k] + KW]);
    }
    Ok((exit[0] + KW, mids, start[last] - KW))
}

/// Rules D then A on the threaded graph; applies the joins (`node_tokens` grows by the junction nodes). Returns (applied
/// joins, the number of dead ends D walked from).
pub fn run(g: &mut Graph, node_tokens: &mut Vec<u32>, k: usize, reads: &Reads, p: &Params) -> Result<(Vec<Applied>, usize)> {
    let index = WordIndex::new(reads);
    let mut body: HashMap<u64, Option<(i64, usize)>> = HashMap::new();
    let mut toks: HashMap<i64, Vec<u32>> = HashMap::new();
    for &u in g.units.keys() {
        let t = spell(g, node_tokens, k, u);
        for (pos, h) in word_hashes(&t).into_iter().enumerate() {
            match body.get(&h) {
                Some(Some((v, _))) if *v != u => {
                    body.insert(h, None);
                }
                None => {
                    body.insert(h, Some((u, pos)));
                }
                _ => {}
            }
        }
        toks.insert(u, t);
    }
    let body: HashMap<u64, (i64, usize)> = body.into_iter().filter_map(|(h, v)| v.map(|v| (h, v))).collect();
    let ctx = Ctx { g, reads, index: &index, body, toks, p };
    let dead_out: Vec<i64> = g.units.iter()
        .filter(|(&u, x)| !x.circular && g.succ(u).is_empty() && x.members.len() >= p.min_nodes).map(|(&u, _)| u).collect();
    let walked: Vec<Option<Join>> = dead_out.par_iter().map(|&a| ctx.hop(a)).collect();
    let mut d_joins: Vec<Join> = walked.into_iter().flatten().collect();
    let mut reach: HashMap<i64, usize> = HashMap::new();
    for j in &d_joins {
        *reach.entry(j.b).or_insert(0) += 1;
    }
    let mut skipped: Vec<Join> = Vec::new();
    d_joins.retain(|j| {
        if reach[&j.b] > 1 || !g.pred(j.b).is_empty() {
            skipped.push(j.clone());
            false
        } else {
            true
        }
    });
    // rule A on the ends D left
    let used_a: HashSet<i64> = d_joins.iter().map(|j| j.a).collect();
    let used_b: HashSet<i64> = d_joins.iter().map(|j| j.b).collect();
    let tails: Vec<i64> = g.units.iter().filter(|(&u, x)| !x.circular && g.succ(u).is_empty() && !used_a.contains(&u))
        .map(|(&u, _)| u).collect();
    let heads: Vec<i64> = g.units.iter().filter(|(&u, x)| !x.circular && g.pred(u).is_empty() && !used_b.contains(&u))
        .map(|(&u, _)| u).collect();
    let a_joins: Vec<Join> = bridges(&ctx, &tails, &heads).into_iter()
        .map(|(a, b, read, n)| Join { rule: "A bridge", a, b, reads: vec![read], hop_q: Vec::new(), hops: 0,
                                      note: format!("{} reads", n) })
        .collect();
    drop(ctx);
    let first_node: HashMap<i64, u32> = g.units.iter().map(|(&u, x)| (u, x.members[0])).collect();
    let last_node: HashMap<i64, u32> = g.units.iter().map(|(&u, x)| (u, *x.members.last().unwrap())).collect();
    // apply, D first
    let mut cur: HashMap<i64, i64> = g.units.keys().map(|&u| (u, u)).collect();
    let find = |cur: &HashMap<i64, i64>, mut u: i64| -> i64 {
        while cur[&u] != u {
            u = cur[&u];
        }
        u
    };
    let mut nxt = g.units.keys().next_back().map_or(0, |m| m + 1);
    let mut applied = Vec::new();
    for j in skipped {
        applied.push(Applied { a_last_node: last_node[&j.a], b_first_node: first_node[&j.b], join: j, merged: -1, dropped_a: 0, dropped_b: 0, gap_tokens: 0, new_nodes: 0,
                               status: "skipped (a start reached from two ends)".into() });
    }
    for j in d_joins.into_iter().chain(a_joins) {
        let (x, y) = (find(&cur, j.a), find(&cur, j.b));
        let mut rec = Applied { join: j.clone(), a_last_node: last_node[&j.a], b_first_node: first_node[&j.b], merged: -1, dropped_a: 0, dropped_b: 0,
                                gap_tokens: 0, new_nodes: 0, status: String::new() };
        if x == y {
            rec.status = "skipped (would close a circle)".into();
            applied.push(rec);
            continue;
        }
        let (xt, yt) = (spell(g, node_tokens, k, x), spell(g, node_tokens, k, y));
        let rd: Vec<&[u32]> = j.reads.iter().map(|&i| reads.read(i)).collect();
        let fallback: Vec<(usize, usize)> = j.hop_q.iter().enumerate()
            .map(|(h, &(q, off))| (reads.tok_at(j.reads[h], q), reads.tok_at(j.reads[h + 1], q - off))).collect();
        match splice(&xt, &yt, &rd, &fallback) {
            Err(why) => {
                rec.status = format!("skipped ({})", why);
            }
            Ok((keep_x, mids, entry_y)) => {
                let kx = keep_x as i64 - k as i64 + 1;
                let (xm, ym) = (g.units[&x].members.clone(), g.units[&y].members.clone());
                if kx < 1 || entry_y >= ym.len() {
                    rec.status = "skipped (the splice leaves no node of A or B)".into();
                    applied.push(rec);
                    continue;
                }
                let kx = kx as usize;
                let mut s: Vec<u32> = xt[..keep_x].to_vec();
                s.extend_from_slice(&mids);
                s.extend_from_slice(&yt[entry_y + KW..]);
                let sy = s.len() - (yt.len() - entry_y);
                let mut members: Vec<u32> = xm[..kx].to_vec();
                let first_new = (node_tokens.len() / k) as u32;
                for pidx in kx..sy {
                    members.push(first_new + (pidx - kx) as u32);
                    node_tokens.extend_from_slice(&s[pidx..pidx + k]);
                }
                members.extend_from_slice(&ym[entry_y..]);
                let (ux, uy) = (g.units[&x].clone(), g.units[&y].clone());
                let new_n = sy - kx;
                let (wx, wy, wn) = (kx as f64, (ym.len() - entry_y) as f64, new_n as f64);
                let mean = |vx: f64, vy: f64| (vx * wx + vy * wy + vx.min(vy) * wn) / (wx + wy + wn);
                let mut origin = ux.origin.clone();
                origin.extend_from_slice(&uy.origin);
                g.add(nxt, members, mean(ux.abundance, uy.abundance), false, Some(origin),
                      Some(mean(ux.abundance2, uy.abundance2)));
                for v in g.pred(x).iter().copied().collect::<Vec<_>>() {
                    g.link(v, nxt);
                }
                for v in g.succ(y).iter().copied().collect::<Vec<_>>() {
                    g.link(nxt, v);
                }
                g.remove(x);
                g.remove(y);
                cur.insert(x, nxt);
                cur.insert(y, nxt);
                cur.insert(nxt, nxt);
                rec.merged = nxt;
                rec.dropped_a = xt.len() - keep_x;
                rec.dropped_b = entry_y;
                rec.gap_tokens = mids.len();
                rec.new_nodes = new_n;
                rec.status = "joined".into();
                nxt += 1;
            }
        }
        applied.push(rec);
    }
    // a unitig made by one join can be merged again by a later one: report the final unitig
    for rec in &mut applied {
        if rec.merged >= 0 {
            rec.merged = find(&cur, rec.merged);
        }
    }
    Ok((applied, dead_out.len()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_median_truncates_the_middle_mean() {
        assert_eq!(py_median(&mut [3, 1, 2]), 2);
        assert_eq!(py_median(&mut [1, 2]), 1);
        assert_eq!(py_median(&mut [-1, -2]), -1);
        assert_eq!(py_median(&mut [5, 8, 1, 2]), 3);
    }

    #[test]
    fn identical_q_symbols_cost_nothing_and_unrelated_ones_cost() {
        let s: Vec<f32> = (0..6000u32).map(|i| ((i * 2654435761u32) >> 29) as f32).collect();
        let x = &s[1000..4000];
        let c = qalign::cost(x, &s, 1000, (2000, 2000, 3000, 3000)).unwrap();
        assert!(c < 1e-9, "{}", c);
        let u: Vec<f32> = (0..6000u32).map(|i| ((i.wrapping_mul(40503) ^ 0x5bd1) >> 13 & 7) as f32).collect();
        assert!(qalign::cost(x, &u, 1000, (2000, 2000, 3000, 3000)).unwrap() > 0.5);
    }

    #[test]
    fn splice_keeps_x_the_reads_between_shared_words_and_y() {
        // X = 0..100, read = 90..200 (with 150 changed in the read), Y = 190..300: one path 0..300
        let x: Vec<u32> = (0..100).collect();
        let mut r: Vec<u32> = (90..200).collect();
        r[60] = 9999;
        let y: Vec<u32> = (190..300).collect();
        let (keep, mids, entry) = splice(&x, &y, &[&r], &[]).unwrap();
        let mut s = x[..keep].to_vec();
        s.extend_from_slice(&mids);
        s.extend_from_slice(&y[entry + KW..]);
        let mut want: Vec<u32> = (0..300).collect();
        want[150] = 9999;
        assert_eq!(s, want);
        assert_eq!((keep, entry), (100, 0));
    }

    #[test]
    fn splice_refuses_words_only_far_from_the_median_diagonal() {
        // two of X's words recur in the read on diagonals -60 and 0 (median -30): neither is within DIAG of it
        let x: Vec<u32> = (0..100).collect();
        let mut r: Vec<u32> = (1000..1200).collect();
        r[..KW].copy_from_slice(&x[60..60 + KW]);
        r[150 - 60..150 - 60 + KW].copy_from_slice(&x[90..90 + KW]);
        let y: Vec<u32> = (2000..2100).collect();
        assert!(splice(&x, &y, &[&r], &[]).is_err());
    }

    #[test]
    fn splice_refuses_pieces_without_a_shared_word() {
        let x: Vec<u32> = (0..100).collect();
        let r: Vec<u32> = (1000..1100).collect();
        let y: Vec<u32> = (2000..2100).collect();
        assert!(splice(&x, &y, &[&r], &[]).is_err());
    }
}

#[cfg(test)]
mod fallback_tests {
    use super::*;

    #[test]
    fn reads_without_a_shared_word_splice_at_the_hop_point() {
        // X 0..100; read 1 = 90..200; read 2 = 5000.. (its own tokens) then 300..400 shared with Y = 300..500
        let x: Vec<u32> = (0..100).collect();
        let r1: Vec<u32> = (90..200).collect();
        let mut r2: Vec<u32> = (5000..5050).collect();
        r2.extend(300..400);
        let y: Vec<u32> = (300..500).collect();
        // hop point: read 1 token 80 (value 170) = read 2 token 20 (value 5020)
        let (keep, mids, entry) = splice(&x, &y, &[&r1, &r2], &[(80, 20)]).unwrap();
        assert_eq!(keep, 100);
        // the exit word into Y (300..306) closes the middle; Y follows from its entry word's end
        let mut want: Vec<u32> = (100..170).collect();
        want.extend(5020..5050);
        want.extend(300..306);
        assert_eq!(mids, want);
        assert_eq!(entry, 0);
    }
}
