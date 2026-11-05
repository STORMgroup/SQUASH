//! Read threading on the final graph (the K-max unitigs): repeats and tangles resolved by the reads that cross
//! them, then, with `--join`, contig joins (`joins`).
//!
//! The reads (the graph's input tokens) are mapped to the graph by exact k-token windows; each read becomes its ordered
//! unitig visits. A repeat unitig R (>= 2 predecessors or successors, as many of each, no self link) is split into one
//! copy per predecessor when the reads entering from a predecessor A and leaving to a successor B pair the predecessors
//! and successors one to one, every pair has >= `--min-support` reads, and the reads of the other pairs of an A (or
//! a B) stay <= `--max-conflict` times its pair's. With `--tangles`, the non-circular unitigs of >= `--anchor-nodes`
//! nodes are anchors, and a connected component of the other unitigs that >= 2 anchors enter or leave is a tangle:
//! an entering anchor A and a leaving anchor B are paired when they are each other's best by the reads crossing from
//! A's end to B's start (within `--anchor-edge` nodes), with >= `--min-support` reads and >= `--tangle-dominance`
//! times the next best; the commonest read route between them that is a path of the graph is copied for the pair
//! (a closed tangle that one anchor enters and another leaves is taken too, `closed_single`).
//! Tangles are threaded first, then the repeats; last, a unitig whose only successor has it as only predecessor is
//! merged with it.
//!
//! The output folder holds the threaded K-max graph (the unitig table and the GFA, in the graph's formats, and the node
//! tokens), unitig_map.tsv (each unitig's original unitigs), joins.tsv (`--join`) and stats.tsv. The reads are mapped
//! on `--threads` threads (default: all); the visits, and so the output, do not depend on the thread count.

use crate::graph::{GFA, NODE_TOKENS, UNITIGS, order_file};
use crate::io::{READS_TABLE, Result, TOKEN_VALUES, py_float, read_u32le, tsv_rows};
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};


const HEAD: &str = "k\tstate\tunitig_id\tcircular\tnodes\ttoken_count\tabundance2\tindegree\toutdegree\tmembers";
/// The base of the rolling window hash (a polynomial hash of each token + 1).
const BASE: u64 = 1_099_511_628_211;
static EMPTY: BTreeSet<i64> = BTreeSet::new();

/// A read's visit to a unitig: (unitig, first node index, last node index) of a run of matched windows.
pub type Visit = (i64, i64, i64);

// ---------------------------------------------------------------- graph

#[derive(Clone, Debug)]
pub struct Unitig {
    pub members: Vec<u32>, // node ids
    pub abundance: f64,    // the GFA's IA
    pub abundance2: f64,   // the unitig table's abundance2
    pub circular: bool,
    pub origin: Vec<i64>, // the input unitigs it is made of (a copy repeats its source's)
}

/// Unitigs as lists of node ids (consecutive nodes overlap k - 1 tokens) with directed links.
pub struct Graph {
    pub k: usize,
    pub units: BTreeMap<i64, Unitig>,
    succ: HashMap<i64, BTreeSet<i64>>,
    pred: HashMap<i64, BTreeSet<i64>>,
}

impl Graph {
    pub fn new(k: usize) -> Self {
        Graph { k, units: BTreeMap::new(), succ: HashMap::new(), pred: HashMap::new() }
    }

    pub fn add(&mut self, uid: i64, members: Vec<u32>, abundance: f64, circular: bool, origin: Option<Vec<i64>>,
               abundance2: Option<f64>) {
        let abundance2 = abundance2.unwrap_or(2.0 * abundance);
        let origin = origin.unwrap_or_else(|| vec![uid]);
        self.units.insert(uid, Unitig { members, abundance, abundance2, circular, origin });
    }

    pub fn link(&mut self, a: i64, b: i64) {
        self.succ.entry(a).or_default().insert(b);
        self.pred.entry(b).or_default().insert(a);
    }

    pub fn unlink(&mut self, a: i64, b: i64) {
        if let Some(s) = self.succ.get_mut(&a) {
            s.remove(&b);
        }
        if let Some(p) = self.pred.get_mut(&b) {
            p.remove(&a);
        }
    }

    pub fn remove(&mut self, u: i64) {
        let outs: Vec<i64> = self.succ(u).iter().copied().collect();
        for v in outs {
            self.unlink(u, v);
        }
        let ins: Vec<i64> = self.pred(u).iter().copied().collect();
        for v in ins {
            self.unlink(v, u);
        }
        self.units.remove(&u);
        self.succ.remove(&u);
        self.pred.remove(&u);
    }

    pub fn succ(&self, u: i64) -> &BTreeSet<i64> {
        self.succ.get(&u).unwrap_or(&EMPTY)
    }

    pub fn pred(&self, u: i64) -> &BTreeSet<i64> {
        self.pred.get(&u).unwrap_or(&EMPTY)
    }

    fn len_of(&self, u: i64) -> i64 {
        self.units[&u].members.len() as i64
    }

    pub fn tokens(&self, u: i64) -> usize {
        self.k + self.units[&u].members.len() - 1
    }

    fn links(&self) -> usize {
        self.succ.values().map(|s| s.len()).sum()
    }

    fn circulars(&self) -> usize {
        self.units.values().filter(|x| x.circular).count()
    }
}

/// One repeat candidate's decision (`resolve_candidates`).
#[derive(Clone, Debug)]
pub struct Decision {
    pub r: i64,
    pub preds: Vec<i64>,
    pub succs: Vec<i64>,
    pub counts: BTreeMap<(i64, i64), u64>,
    pub matching: Option<BTreeMap<i64, i64>>,
}

/// The entry of `keys` with the largest count, ties to the first (the smallest key, as `keys` are sorted).
fn best_by(keys: &[i64], count: impl Fn(i64) -> u64) -> i64 {
    let mut best = keys[0];
    for &x in &keys[1..] {
        if count(x) > count(best) {
            best = x;
        }
    }
    best
}

/// Every repeat unitig's decision, in unitig order.
pub fn resolve_candidates(g: &Graph, triples: &HashMap<(i64, i64, i64), u64>, min_support: i64, max_conflict: f64)
                          -> Vec<Decision> {
    let mut out = Vec::new();
    for &r in g.units.keys() {
        let preds: Vec<i64> = g.pred(r).iter().copied().collect();
        let succs: Vec<i64> = g.succ(r).iter().copied().collect();
        if preds.len() < 2 && succs.len() < 2 {
            continue;
        }
        let mut counts = BTreeMap::new();
        for &a in &preds {
            for &b in &succs {
                counts.insert((a, b), triples.get(&(a, r, b)).copied().unwrap_or(0));
            }
        }
        let decision = |matching| Decision { r, preds: preds.clone(), succs: succs.clone(), counts: counts.clone(), matching };
        if preds.contains(&r) || succs.contains(&r) || preds.len() != succs.len() {
            out.push(decision(None));
            continue;
        }
        let c = |a: i64, b: i64| counts[&(a, b)];
        let matching: BTreeMap<i64, i64> = preds.iter().map(|&a| (a, best_by(&succs, |b| c(a, b)))).collect();
        // one to one, every pair with >= min_support reads, and no conflict above max_conflict on either side
        let mut resolved = matching.values().collect::<HashSet<_>>().len() == preds.len()
            && (preds.iter().map(|&a| c(a, matching[&a])).min().unwrap() as i64) >= min_support;
        if resolved {
            let back: HashMap<i64, i64> = matching.iter().map(|(&a, &b)| (b, a)).collect();
            for &a in &preds {
                let other: u64 = succs.iter().filter(|&&b| b != matching[&a]).map(|&b| c(a, b)).sum();
                if other as f64 > max_conflict * c(a, matching[&a]) as f64 {
                    resolved = false;
                }
            }
            for &b in &succs {
                let other: u64 = preds.iter().filter(|&&a| a != back[&b]).map(|&a| c(a, b)).sum();
                if other as f64 > max_conflict * c(back[&b], b) as f64 {
                    resolved = false;
                }
            }
        }
        out.push(decision(if resolved { Some(matching) } else { None }));
    }
    out
}

/// Split the resolved repeats in unitig order (one pass); a repeat whose neighbourhood (R, its predecessors and
/// successors) meets an earlier split repeat's or `taken` (unitigs already changed by tangle threading) is blocked.
/// Returns (split R with their copy ids, blocked R).
pub fn apply(g: &mut Graph, decisions: &[Decision], taken: &BTreeSet<i64>) -> (Vec<(i64, Vec<i64>)>, Vec<i64>) {
    let mut taken = taken.clone();
    let (mut chosen, mut blocked) = (Vec::new(), Vec::new());
    for d in decisions {
        if d.matching.is_none() {
            continue;
        }
        let hood: BTreeSet<i64> = d.preds.iter().chain(d.succs.iter()).copied().chain(std::iter::once(d.r)).collect();
        if hood.iter().any(|u| taken.contains(u)) {
            blocked.push(d.r);
            continue;
        }
        chosen.push(d);
        taken.extend(hood);
    }
    let mut done = Vec::new();
    let mut nxt = g.units.keys().next_back().map_or(0, |m| m + 1);
    for d in chosen {
        let matching = d.matching.as_ref().unwrap();
        let mut total: u64 = matching.iter().map(|(&a, &b)| d.counts[&(a, b)]).sum();
        if total == 0 {
            total = 1;
        }
        let r = g.units[&d.r].clone();
        let mut copies = Vec::new();
        for (&a, &b) in matching {
            let share = d.counts[&(a, b)] as f64 / total as f64;
            g.add(nxt, r.members.clone(), r.abundance * share, false, Some(r.origin.clone()), Some(r.abundance2 * share));
            g.link(a, nxt);
            g.link(nxt, b);
            copies.push(nxt);
            nxt += 1;
        }
        g.remove(d.r);
        done.push((d.r, copies));
    }
    (done, blocked)
}

/// Merge u -> v when u's only successor is v and v's only predecessor is u (u != v); a chain that closes on itself
/// becomes one circular unitig.
pub fn compact(g: &mut Graph) {
    let mut merged = true;
    while merged {
        merged = false;
        let ids: Vec<i64> = g.units.keys().copied().collect();
        for u in ids {
            if !g.units.contains_key(&u) || g.succ(u).len() != 1 {
                continue;
            }
            let v = *g.succ(u).iter().next().unwrap();
            if v == u {
                if g.pred(u).len() == 1 {
                    g.unlink(u, u);
                    g.units.get_mut(&u).unwrap().circular = true;
                    merged = true;
                }
                continue;
            }
            if g.pred(v).len() != 1 {
                continue;
            }
            let uv = g.units[&v].clone();
            let x = g.units.get_mut(&u).unwrap();
            let (nu, nv) = (x.members.len() as f64, uv.members.len() as f64);
            x.abundance = (x.abundance * nu + uv.abundance * nv) / (nu + nv);
            x.abundance2 = (x.abundance2 * nu + uv.abundance2 * nv) / (nu + nv);
            x.members.extend_from_slice(&uv.members);
            x.origin.extend_from_slice(&uv.origin);
            let outs: Vec<i64> = g.succ(v).iter().copied().collect();
            g.remove(v);
            for w in outs {
                g.link(u, if w == v { u } else { w });
            }
            merged = true;
        }
    }
}

// ---------------------------------------------------------------- tangles (`--tangles`)

pub type Tangles = (BTreeMap<i64, i64>, BTreeMap<i64, BTreeSet<i64>>, BTreeMap<i64, BTreeSet<i64>>);

/// Connected components of the non-anchor unitigs (links in either direction) that >= 2 anchors enter or leave:
/// (unitig -> tangle id = its smallest unitig, tangle -> entering anchors, tangle -> leaving anchors), with a closed tangle
/// that one anchor enters and another leaves (`closed_single`).
pub fn tangles(g: &Graph, anchors: &HashSet<i64>) -> Tangles {
    let mut comp: BTreeMap<i64, i64> = BTreeMap::new();
    for &start in g.units.keys() {
        if anchors.contains(&start) || comp.contains_key(&start) {
            continue;
        }
        let (mut stack, mut seen) = (vec![start], BTreeSet::from([start]));
        while let Some(u) = stack.pop() {
            for &v in g.succ(u).iter().chain(g.pred(u).iter()) {
                if !anchors.contains(&v) && seen.insert(v) {
                    stack.push(v);
                }
            }
        }
        let low = *seen.iter().next().unwrap();
        for u in seen {
            comp.insert(u, low);
        }
    }
    let mut entries: BTreeMap<i64, BTreeSet<i64>> = BTreeMap::new();
    let mut exits: BTreeMap<i64, BTreeSet<i64>> = BTreeMap::new();
    for &a in anchors {
        for v in g.succ(a) {
            if let Some(&c) = comp.get(v) {
                entries.entry(c).or_default().insert(a);
            }
        }
        for v in g.pred(a) {
            if let Some(&c) = comp.get(v) {
                exits.entry(c).or_default().insert(a);
            }
        }
    }
    let ids: BTreeSet<i64> = comp.values().copied().collect();
    let keep: BTreeSet<i64> = ids.into_iter().filter(|c| {
        let (ins, outs) = (entries.get(c).unwrap_or(&EMPTY), exits.get(c).unwrap_or(&EMPTY));
        ins.len() >= 2 || outs.len() >= 2 || closed_single(g, &comp, *c, ins, outs)
    }).collect();
    (comp.into_iter().filter(|(_, c)| keep.contains(c)).collect(),
     keep.iter().map(|&c| (c, entries.get(&c).cloned().unwrap_or_default())).collect(),
     keep.iter().map(|&c| (c, exits.get(&c).cloned().unwrap_or_default())).collect())
}

/// One anchor A enters tangle c, another anchor B leaves it, A -> B is a link, and nothing else leaves A or enters B.
pub fn closed_single(g: &Graph, comp: &BTreeMap<i64, i64>, c: i64, ins: &BTreeSet<i64>, outs: &BTreeSet<i64>) -> bool {
    if ins.len() != 1 || outs.len() != 1 || ins == outs {
        return false;
    }
    let (a, b) = (*ins.iter().next().unwrap(), *outs.iter().next().unwrap());
    g.succ(a).contains(&b) && g.succ(a).iter().all(|v| comp.get(v) == Some(&c) || *v == b)
        && g.pred(b).iter().all(|v| comp.get(v) == Some(&c) || *v == a)
}

pub type PairCounts = HashMap<(i64, i64, i64), u64>;
pub type Routes = HashMap<(i64, i64, i64), HashMap<Vec<i64>, u64>>;

/// Reads crossing a tangle from an anchor's end to an anchor's start: pair counts {(tangle, A, B): n} and routes
/// {(tangle, A, B): {route: n}}. `direct` (closed single tangles): {(A, B): tangle}; a read going from A's end
/// straight into B's start counts for it with the empty route.
#[allow(clippy::too_many_arguments)]
pub fn count_anchor_pairs(g: &Graph, visits: &[Vec<Visit>], anchors: &HashSet<i64>, comp: &BTreeMap<i64, i64>,
                          entries: &BTreeMap<i64, BTreeSet<i64>>, exits: &BTreeMap<i64, BTreeSet<i64>>, edge: i64,
                          direct: &HashMap<(i64, i64), i64>) -> (PairCounts, Routes) {
    let (mut pairs, mut routes): (PairCounts, Routes) = (HashMap::new(), HashMap::new());
    for seq in visits {
        let mut last: Option<(i64, i64)> = None;
        let mut middle: Vec<i64> = Vec::new();
        for &(u, first, last_pos) in seq {
            if anchors.contains(&u) {
                if let Some((a, a_pos)) = last {
                    if !direct.is_empty() && middle.is_empty() {
                        if let Some(&c) = direct.get(&(a, u)) {
                            if a_pos >= g.len_of(a) - 1 - edge && first <= edge {
                                *pairs.entry((c, a, u)).or_insert(0) += 1;
                                *routes.entry((c, a, u)).or_default().entry(Vec::new()).or_insert(0) += 1;
                            }
                        }
                    }
                    if !middle.is_empty() {
                        let cs: HashSet<Option<i64>> = middle.iter().map(|x| comp.get(x).copied()).collect();
                        if cs.len() == 1 && !cs.contains(&None) && a_pos >= g.len_of(a) - 1 - edge && first <= edge {
                            let c = cs.into_iter().next().unwrap().unwrap();
                            if entries[&c].contains(&a) && exits[&c].contains(&u) {
                                *pairs.entry((c, a, u)).or_insert(0) += 1;
                                *routes.entry((c, a, u)).or_default().entry(middle.clone()).or_insert(0) += 1;
                            }
                        }
                    }
                }
                last = Some((u, last_pos));
                middle.clear();
            } else if last.is_some() {
                middle.push(u);
            }
        }
    }
    (pairs, routes)
}

fn is_path(g: &Graph, a: i64, route: &[i64], b: i64) -> bool {
    let chain: Vec<i64> = std::iter::once(a).chain(route.iter().copied()).chain(std::iter::once(b)).collect();
    chain.windows(2).all(|w| g.succ(w[0]).contains(&w[1]))
}

/// One tangle's decision (`resolve_tangles`): its entering and leaving anchors and the accepted pairs (A, B, reads,
/// route), revised by `apply_tangles`.
#[derive(Clone, Debug)]
pub struct TangleDecision {
    pub c: i64,
    pub ins: Vec<i64>,
    pub outs: Vec<i64>,
    pub accepted: Vec<(i64, i64, u64, Vec<i64>)>,
}

pub fn resolve_tangles(g: &Graph, entries: &BTreeMap<i64, BTreeSet<i64>>, exits: &BTreeMap<i64, BTreeSet<i64>>,
                       pairs: &PairCounts, routes: &Routes, min_support: i64, dominance: f64) -> Vec<TangleDecision> {
    let mut out = Vec::new();
    let ids: BTreeSet<i64> = entries.keys().chain(exits.keys()).copied().collect();
    for c in ids {
        let ins: Vec<i64> = entries.get(&c).map_or(Vec::new(), |s| s.iter().copied().collect());
        let outs: Vec<i64> = exits.get(&c).map_or(Vec::new(), |s| s.iter().copied().collect());
        let mut counts = BTreeMap::new();
        for &a in &ins {
            for &b in &outs {
                counts.insert((a, b), pairs.get(&(c, a, b)).copied().unwrap_or(0));
            }
        }
        let n_of = |a: i64, b: i64| counts[&(a, b)];
        let mut accepted = Vec::new();
        for &a in &ins {
            if outs.is_empty() {
                continue;
            }
            let b = best_by(&outs, |x| n_of(a, x));
            let n = n_of(a, b);
            let second_out = outs.iter().filter(|&&x| x != b).map(|&x| n_of(a, x)).max().unwrap_or(0);
            let second_in = ins.iter().filter(|&&y| y != a).map(|&y| n_of(y, b)).max().unwrap_or(0);
            // >= min_support reads, mutual best, and >= dominance x the next best on both sides
            if (n as i64) < min_support
                || best_by(&ins, |y| n_of(y, b)) != a
                || (n as f64) < dominance * second_out as f64
                || (n as f64) < dominance * second_in as f64
            {
                continue;
            }
            let best = routes.get(&(c, a, b)).and_then(|rs| {
                rs.iter().filter(|(r, _)| is_path(g, a, r, b))
                    .max_by_key(|(r, &m)| (m, r.iter().map(|v| -v).collect::<Vec<i64>>()))
                    .map(|(r, _)| r.clone())
            });
            if let Some(route) = best {
                accepted.push((a, b, n, route));
            }
        }
        out.push(TangleDecision { c, ins, outs, accepted });
    }
    out
}

/// Thread every accepted pair; remove the tangles whose anchors are all paired: (changed unitigs, tangles removed).
pub fn apply_tangles(g: &mut Graph, comp: &BTreeMap<i64, i64>, decisions: &mut [TangleDecision]) -> (BTreeSet<i64>, usize) {
    let mut members_of: BTreeMap<i64, Vec<i64>> = BTreeMap::new();
    for (&u, &c) in comp {
        members_of.entry(c).or_default().push(u);
    }
    let mut nxt = g.units.keys().next_back().map_or(0, |m| m + 1);
    let (mut touched, mut removed) = (BTreeSet::new(), 0);
    let (mut paired_a, mut paired_b) = (HashSet::new(), HashSet::new());
    for d in decisions.iter_mut() {
        // An anchor next to two tangles is paired in the first only.
        d.accepted.retain(|x| !paired_a.contains(&x.0) && !paired_b.contains(&x.1));
        paired_a.extend(d.accepted.iter().map(|x| x.0));
        paired_b.extend(d.accepted.iter().map(|x| x.1));
        if d.accepted.is_empty() {
            continue;
        }
        let mut used: HashMap<i64, u64> = HashMap::new();
        for (_, _, n, route) in &d.accepted {
            for &u in route {
                *used.entry(u).or_insert(0) += n;
            }
        }
        for (a, b, n, route) in d.accepted.clone() {
            for v in g.succ(a).iter().copied().collect::<Vec<_>>() {
                g.unlink(a, v);
            }
            for v in g.pred(b).iter().copied().collect::<Vec<_>>() {
                g.unlink(v, b);
            }
            let mut prev = a;
            for &u in &route {
                let share = n as f64 / used[&u] as f64;
                let x = g.units[&u].clone();
                g.add(nxt, x.members, x.abundance * share, false, Some(x.origin), Some(x.abundance2 * share));
                g.link(prev, nxt);
                prev = nxt;
                nxt += 1;
            }
            g.link(prev, b);
            touched.insert(a);
            touched.insert(b);
        }
        let members = members_of.get(&d.c).cloned().unwrap_or_default();
        touched.extend(members.iter().copied());
        let got_a: BTreeSet<i64> = d.accepted.iter().map(|x| x.0).collect();
        let got_b: BTreeSet<i64> = d.accepted.iter().map(|x| x.1).collect();
        if got_a == d.ins.iter().copied().collect() && got_b == d.outs.iter().copied().collect() {
            for u in members {
                g.remove(u);
            }
            removed += 1;
        }
    }
    (touched, removed)
}

// ---------------------------------------------------------------- reads -> unitig visits

pub fn window_hash(tokens: &[u32]) -> u64 {
    tokens.iter().fold(0u64, |h, &t| h.wrapping_mul(BASE).wrapping_add(t as u64 + 1))
}

/// The nodes that belong to a unitig, keyed by their window hash: (hash, node), sorted. A read's window is
/// looked up at the first entry of its hash only.
pub fn node_index(node_tokens: &[u32], k: usize, node_unitig: &[i64]) -> Vec<(u64, u32)> {
    let mut keyed: Vec<(u64, u32)> = (0..node_unitig.len()).into_par_iter().filter(|&i| node_unitig[i] >= 0)
        .map(|i| (window_hash(&node_tokens[i * k..(i + 1) * k]), i as u32)).collect();
    keyed.par_sort_unstable();
    keyed
}

/// One read's ordered unitig visits (runs of one unitig collapsed; unmatched windows skipped; a window counts only
/// when it equals a node's tokens): (visits, matched windows, windows).
pub fn read_visit(tok: &[u32], k: usize, keyed: &[(u64, u32)], node_tokens: &[u32], node_unitig: &[i64], node_pos: &[i64])
                  -> (Vec<Visit>, u64, u64) {
    if tok.len() < k {
        return (Vec::new(), 0, 0);
    }
    let n = tok.len() - k + 1;
    if keyed.is_empty() {
        return (Vec::new(), 0, n as u64);
    }
    let top = (1..k).fold(1u64, |p, _| p.wrapping_mul(BASE));
    let mut h = window_hash(&tok[..k]);
    let (mut seq, mut matched): (Vec<Visit>, u64) = (Vec::new(), 0);
    for i in 0..n {
        if i > 0 {
            h = h.wrapping_sub((tok[i - 1] as u64 + 1).wrapping_mul(top)).wrapping_mul(BASE)
                .wrapping_add(tok[i + k - 1] as u64 + 1);
        }
        let j = keyed.partition_point(|x| x.0 < h);
        if j == keyed.len() || keyed[j].0 != h {
            continue;
        }
        let x = keyed[j].1 as usize;
        if node_tokens[x * k..(x + 1) * k] != tok[i..i + k] {
            continue;
        }
        matched += 1;
        let (u, p) = (node_unitig[x], node_pos[x]);
        match seq.last_mut() {
            Some(last) if last.0 == u => {
                last.1 = last.1.min(p);
                last.2 = last.2.max(p);
            }
            _ => seq.push((u, p, p)),
        }
    }
    (seq, matched, n as u64)
}

/// For each visit to a repeat unitig R: the read's entry A (a predecessor of R) and exit B (a successor of R), counted
/// as (A, R, B) only when both are decided outside R's k - 1-node junction zone: visits to successors that end within
/// their first k - 1 nodes, and to predecessors that start within their last k - 1 nodes, are passed over.
pub fn count_triples(g: &Graph, visits: &[Vec<Visit>]) -> HashMap<(i64, i64, i64), u64> {
    let mut c = HashMap::new();
    let margin = |u: i64| (g.k as i64 - 1).min(g.len_of(u) - 1);
    for seq in visits {
        for (i, &(r, _, _)) in seq.iter().enumerate() {
            if !g.units.contains_key(&r) || (g.pred(r).len() < 2 && g.succ(r).len() < 2) {
                continue;
            }
            let mut j = i + 1;
            while j < seq.len() && g.succ(r).contains(&seq[j].0) && seq[j].2 < margin(seq[j].0) {
                j += 1;
            }
            let mut h = i as i64 - 1;
            while h >= 0 && g.pred(r).contains(&seq[h as usize].0)
                && seq[h as usize].1 > g.len_of(seq[h as usize].0) - 1 - margin(seq[h as usize].0) {
                h -= 1;
            }
            if j >= seq.len() || h < 0 {
                continue;
            }
            let (a, b) = (seq[h as usize].0, seq[j].0);
            if g.pred(r).contains(&a) && g.succ(r).contains(&b) {
                *c.entry((a, r, b)).or_insert(0) += 1;
            }
        }
    }
    c
}

// ---------------------------------------------------------------- I/O

/// The K-max graph of a graph folder (unitigs in table order) and its node tokens.
pub fn load_graph(order_dir: &Path, k: usize) -> Result<(Graph, Vec<u32>, Vec<i64>)> {
    let node_tokens = read_u32le(&order_file(order_dir, k, NODE_TOKENS))?;
    let mut g = Graph::new(k);
    let mut order = Vec::new();
    tsv_rows(&order_file(order_dir, k, UNITIGS),
            &["unitig_id", "members", "token_count", "circular", "abundance2"], |r| {
        let members: Vec<u32> = r[1].split(',').map(|x| x.parse()).collect::<std::result::Result<_, _>>()?;
        if members.len() + k - 1 != r[2].parse::<usize>()? {
            return Err(format!("unitig {}: members and token_count disagree", r[0]).into());
        }
        let u: i64 = r[0].parse()?;
        g.add(u, members, 0.0, r[3] == "1", None, Some(r[4].parse()?));
        order.push(u);
        Ok(())
    })?;
    let gfa = order_file(order_dir, k, GFA);
    for line in fs::read_to_string(&gfa).map_err(|e| format!("{}: {}", gfa.display(), e))?.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f[0] == "S" {
            let ia = f.get(3..).unwrap_or(&[]).iter().find(|x| x.starts_with("IA:i:"))
                .ok_or_else(|| format!("{}: segment {} without IA", gfa.display(), f[1]))?;
            let value: f64 = ia[5..].parse()?;
            if let Some(x) = g.units.get_mut(&f[1].parse::<i64>()?) {
                x.abundance = value;
            }
        }
        if f[0] == "L" {
            if f[2] != "+" || f[4] != "+" {
                return Err("the graph is directional: every link is + to +".into());
            }
            g.link(f[1].parse()?, f[3].parse()?);
        }
    }
    Ok((g, node_tokens, order))
}

/// The tokens of the graph's input reads (a folder of token arrays), one slice range per read in table order.
fn load_reads(reads_dir: &Path) -> Result<(Vec<u32>, Vec<(usize, usize)>)> {
    let tokens = read_u32le(&reads_dir.join(TOKEN_VALUES))?;
    let mut spans = Vec::new();
    tsv_rows(&reads_dir.join(READS_TABLE), &["token_offset", "token_count"], |r| {
        let (o, n): (usize, usize) = (r[0].parse()?, r[1].parse()?);
        // a range past the end of the tokens is clipped
        spans.push((o.min(tokens.len()), (o + n).min(tokens.len())));
        Ok(())
    })?;
    Ok((tokens, spans))
}

fn round_even(x: f64) -> i64 {
    x.round_ties_even() as i64
}

fn join<T: ToString>(xs: impl IntoIterator<Item = T>, sep: &str) -> String {
    xs.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join(sep)
}

/// Renumber 0..n-1 (by smallest original id, then first member) and write the unitig table and GFA in the graph's
/// formats.
pub fn write_graph(g: &Graph, order_dir: &Path, k: usize) -> Result<HashMap<i64, usize>> {
    let mut order: Vec<i64> = g.units.keys().copied().collect();
    order.sort_by_key(|u| (*g.units[u].origin.iter().min().unwrap(), g.units[u].members[0], *u));
    let new: HashMap<i64, usize> = order.iter().enumerate().map(|(i, &u)| (u, i)).collect();
    let file = File::create(order_file(order_dir, k, UNITIGS))?;
    let mut gz = GzEncoder::new(BufWriter::new(file), Compression::default());
    writeln!(gz, "{}", HEAD)?;
    for &u in &order {
        let x = &g.units[&u];
        writeln!(gz, "{}\tP3\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", k, new[&u], x.circular as i32, x.members.len(), g.tokens(u),
                 round_even(x.abundance2).max(1), g.pred(u).len(), g.succ(u).len(), join(&x.members, ","))?;
    }
    gz.finish()?.flush()?;
    let mut gfa = BufWriter::new(File::create(order_file(order_dir, k, GFA))?);
    writeln!(gfa, "H\tVN:Z:1.0")?;
    for &u in &order {
        writeln!(gfa, "S\t{}\t*\tLN:i:{}\tIA:i:{}", new[&u], g.tokens(u), round_even(g.units[&u].abundance).max(1))?;
    }
    for &u in &order {
        let mut outs: Vec<usize> = g.succ(u).iter().map(|v| new[v]).collect();
        outs.sort_unstable();
        for v in outs {
            writeln!(gfa, "L\t{}\t+\t{}\t+\t{}M", new[&u], v, k - 1)?;
        }
    }
    gfa.flush()?;
    Ok(new)
}

// ---------------------------------------------------------------- main

struct Args {
    graph_dir: PathBuf,
    k: usize,
    min_support: i64,
    max_conflict: f64,
    tangles: bool,
    anchor_nodes: i64,
    anchor_edge: i64,
    tangle_dominance: f64,
    output_dir: PathBuf,
    threads: usize,
    join: bool,
    join_symbols: Option<PathBuf>,
    join_partners: Option<PathBuf>,
    join_threshold: f64,
    join_min_nodes: usize,
    join_head: usize,
    join_min_paths: usize,
}

const USAGE: &str = "usage: squash __thread --graph-dir GRAPH --output-dir OUT [--k 32] [--min-support 3] [--max-conflict 0.2] \
[--tangles [--anchor-nodes 200] [--anchor-edge 40] [--tangle-dominance 3.0]] [--threads N] \
[--join --join-symbols ROOT/signal/qsymbols.tsv --join-partners ROOT/correction/partners.tsv.gz \
[--join-threshold 0.62] [--join-min-nodes 200] [--join-head 200] [--join-min-paths 1]]";

fn parse_args() -> Result<Args> {
    let mut a = Args { graph_dir: PathBuf::new(), k: 32, min_support: 3, max_conflict: 0.2, tangles: false, anchor_nodes: 200,
                       anchor_edge: 40, tangle_dominance: 3.0, output_dir: PathBuf::new(), threads: 0,
                       join: false, join_symbols: None, join_partners: None, join_threshold: 0.62, join_min_nodes: 200,
                       join_head: 200, join_min_paths: 1 };
    let (mut graph_dir, mut out) = (None, None);
    let mut argv = crate::argv::args().skip(1);
    while let Some(key) = argv.next() {
        match key.as_str() {
            "--tangles" => a.tangles = true,
            "--join" => a.join = true,
            "-h" | "--help" => return Err(USAGE.into()),
            _ => {
                let value = argv.next().ok_or_else(|| format!("{} needs a value\n{}", key, USAGE))?;
                match key.as_str() {
                    "--graph-dir" => graph_dir = Some(PathBuf::from(value)),
                    "--output-dir" => out = Some(PathBuf::from(value)),
                    "--k" => a.k = value.parse()?,
                    "--min-support" => a.min_support = value.parse()?,
                    "--max-conflict" => a.max_conflict = value.parse()?,
                    "--anchor-nodes" => a.anchor_nodes = value.parse()?,
                    "--anchor-edge" => a.anchor_edge = value.parse()?,
                    "--tangle-dominance" => a.tangle_dominance = value.parse()?,
                    "--threads" => a.threads = value.parse()?,
                    "--join-symbols" => a.join_symbols = Some(PathBuf::from(value)),
                    "--join-partners" => a.join_partners = Some(PathBuf::from(value)),
                    "--join-threshold" => a.join_threshold = value.parse()?,
                    "--join-min-nodes" => a.join_min_nodes = value.parse()?,
                    "--join-head" => a.join_head = value.parse()?,
                    "--join-min-paths" => a.join_min_paths = value.parse()?,
                    _ => return Err(format!("unknown option {}\n{}", key, USAGE).into()),
                }
            }
        }
    }
    a.graph_dir = graph_dir.ok_or_else(|| format!("missing --graph-dir\n{}", USAGE))?;
    a.output_dir = out.ok_or_else(|| format!("missing --output-dir\n{}", USAGE))?;
    if a.k == 0 {
        return Err("--k must be at least 1".into());
    }
    if a.join && (a.join_symbols.is_none() || a.join_partners.is_none()) {
        return Err("--join needs --join-symbols and --join-partners".into());
    }
    Ok(a)
}

pub fn run() -> Result<()> {
    let a = parse_args()?;
    if a.output_dir.exists() {
        return Err(format!("refusing to overwrite {}", a.output_dir.display()).into());
    }
    if a.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(a.threads).build_global()?;
    }
    let k = a.k;
    let src = a.graph_dir.clone();
    let (mut g, mut node_tokens, order) = load_graph(&src, k)?;
    let before = (g.units.len(), g.links(), g.circulars());
    let n_nodes = node_tokens.len() / k;
    let (mut node_unitig, mut node_pos) = (vec![-1i64; n_nodes], vec![-1i64; n_nodes]);
    for &u in &order {
        for (i, &m) in g.units[&u].members.iter().enumerate() {
            let m = m as usize;
            if m >= n_nodes {
                return Err(format!("atomic node {} out of range", m).into());
            }
            if node_unitig[m] >= 0 {
                return Err(format!("atomic node {} in two unitigs", m).into());
            }
            node_unitig[m] = u;
            node_pos[m] = i as i64;
        }
    }
    let reads_dir = a.graph_dir.join("input");
    let (tokens, spans) = load_reads(&reads_dir)?;
    let keyed = node_index(&node_tokens, k, &node_unitig);
    let per_read: Vec<(Vec<Visit>, u64, u64)> = spans.par_iter()
        .map(|&(s, e)| read_visit(&tokens[s..e], k, &keyed, &node_tokens, &node_unitig, &node_pos)).collect();
    let matched: u64 = per_read.iter().map(|x| x.1).sum();
    let windows: u64 = per_read.iter().map(|x| x.2).sum();
    let visits: Vec<Vec<Visit>> = per_read.into_iter().map(|x| x.0).collect();
    drop(keyed);
    let triples = count_triples(&g, &visits);
    let decisions = resolve_candidates(&g, &triples, a.min_support, a.max_conflict);
    let (mut tangle_decisions, mut touched, mut tangles_removed) = (Vec::new(), BTreeSet::new(), 0);
    if a.tangles {
        let anchors: HashSet<i64> = g.units.iter()
            .filter(|(_, x)| x.members.len() as i64 >= a.anchor_nodes && !x.circular).map(|(&u, _)| u).collect();
        let (comp, entries, exits) = tangles(&g, &anchors);
        let mut direct = HashMap::new();
        for (&c, ins) in &entries {
            if closed_single(&g, &comp, c, ins, &exits[&c]) {
                let (x, y) = (*ins.iter().next().unwrap(), *exits[&c].iter().next().unwrap());
                if g.succ(x).contains(&y) {
                    direct.insert((x, y), c);
                }
            }
        }
        let (pairs, routes) = count_anchor_pairs(&g, &visits, &anchors, &comp, &entries, &exits, a.anchor_edge, &direct);
        tangle_decisions = resolve_tangles(&g, &entries, &exits, &pairs, &routes, a.min_support, a.tangle_dominance);
        let (t, r) = apply_tangles(&mut g, &comp, &mut tangle_decisions);
        touched = t;
        tangles_removed = r;
    }
    let (done, blocked) = apply(&mut g, &decisions, &touched);
    compact(&mut g);
    // `--join`: contig joins (joins) after threading; the node tokens grow by the junction nodes
    let n_node_tokens = node_tokens.len();
    let (join_applied, join_walked, join_ids) = if a.join {
        let reads = crate::joins::load(&reads_dir, tokens, spans, a.join_symbols.as_ref().unwrap(),
                                             a.join_partners.as_ref().unwrap())?;
        let params = crate::joins::Params { threshold: a.join_threshold, min_nodes: a.join_min_nodes, head: a.join_head,
                                                  min_paths: a.join_min_paths };
        let (applied, walked) = crate::joins::run(&mut g, &mut node_tokens, k, &reads, &params)?;
        (applied, walked, reads.ids)
    } else {
        drop((tokens, spans));
        (Vec::new(), 0, Vec::new())
    };

    let parent = match a.output_dir.parent() {
        Some(p) if !p.as_os_str().is_empty() => p.to_path_buf(),
        _ => PathBuf::from("."),
    };
    fs::create_dir_all(&parent)?;
    let name = a.output_dir.file_name().ok_or("--output-dir has no name")?.to_string_lossy().to_string();
    let nanos = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)?.subsec_nanos();
    let stage = parent.join(format!(".{}.tmp.{}{:08x}", name, std::process::id(), nanos));
    fs::create_dir(&stage)?;
    let written = (|| -> Result<()> {
        let order_dir = stage.clone();
        let (from, to) = (order_file(&src, k, NODE_TOKENS), order_file(&order_dir, k, NODE_TOKENS));
        if node_tokens.len() > n_node_tokens {
            // joins added junction nodes: the graph's nodes, then the new ones
            let mut bytes = Vec::with_capacity(node_tokens.len() * 4);
            for v in &node_tokens {
                bytes.extend_from_slice(&v.to_le_bytes());
            }
            fs::write(&to, bytes)?;
        } else if fs::hard_link(&from, &to).is_err() {
            fs::copy(&from, &to)?;
        }
        let new = write_graph(&g, &order_dir, k)?;
        let mut by_new: Vec<(usize, i64)> = new.iter().map(|(&u, &i)| (i, u)).collect();
        by_new.sort_unstable();
        let mut fh = BufWriter::new(File::create(stage.join("unitig_map.tsv"))?);
        writeln!(fh, "unitig_id\ttokens\tcircular\toriginal_unitigs")?;
        for (i, u) in by_new {
            writeln!(fh, "{}\t{}\t{}\t{}", i, g.tokens(u), g.units[&u].circular as i32, join(&g.units[&u].origin, ","))?;
        }
        fh.flush()?;
        if a.join {
            let mut fh = BufWriter::new(File::create(stage.join("joins.tsv"))?);
            writeln!(fh, "rule\tfrom_last_node\tto_first_node\tunitig_id\treads\thops\tnote\tstatus\tdropped_end_tokens\t\
                          dropped_start_tokens\tgap_tokens\tnew_nodes")?;
            for x in &join_applied {
                let merged = if x.merged >= 0 { new.get(&x.merged).map_or("-".to_string(), |v| v.to_string()) } else { "-".into() };
                writeln!(fh, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", x.join.rule, x.a_last_node, x.b_first_node, merged,
                         join(x.join.reads.iter().map(|&i| join_ids[i].clone()), ","), x.join.hops,
                         if x.join.note.is_empty() { "-" } else { &x.join.note }, x.status, x.dropped_a, x.dropped_b,
                         x.gap_tokens, x.new_nodes)?;
            }
            fh.flush()?;
        }
        let after = (g.units.len(), g.links(), g.circulars());
        let threaded: usize = tangle_decisions.iter().map(|d| d.accepted.len()).sum();
        let mut stats: Vec<(&str, String)> = vec![
            ("reads", visits.len().to_string()), ("read_windows", windows.to_string()),
            ("matched_windows", matched.to_string()),
            ("reads_visiting_2plus_unitigs", visits.iter().filter(|v| v.len() >= 2).count().to_string()),
            ("repeat_candidates", decisions.len().to_string()), ("resolved", done.len().to_string()),
            ("blocked", blocked.len().to_string()), ("unitigs_before", before.0.to_string()),
            ("links_before", before.1.to_string()), ("circular_before", before.2.to_string()),
            ("unitigs_after", after.0.to_string()), ("links_after", after.1.to_string()),
            ("circular_after", after.2.to_string()), ("min_support", a.min_support.to_string()),
            ("max_conflict", py_float(a.max_conflict)), ("k", k.to_string())];
        if a.tangles {
            stats.extend([("tangles", tangle_decisions.len().to_string()), ("tangle_pairs_threaded", threaded.to_string()),
                          ("tangles_removed", tangles_removed.to_string()), ("anchor_nodes", a.anchor_nodes.to_string()),
                          ("anchor_edge", a.anchor_edge.to_string()), ("tangle_dominance", py_float(a.tangle_dominance))]);
            stats.push(("tangle_single", "1".to_string()));
        }
        let joined = |rule: &str| join_applied.iter().filter(|x| x.join.rule == rule && x.status == "joined").count();
        if a.join {
            stats.extend([("join_ends_walked", join_walked.to_string()), ("joins_hop", joined("D hop").to_string()),
                          ("joins_bridge", joined("A bridge").to_string()),
                          ("joins_skipped", join_applied.iter().filter(|x| x.status != "joined").count().to_string()),
                          ("join_new_nodes", (node_tokens.len() - n_node_tokens).checked_div(k).unwrap_or(0).to_string()),
                          ("join_threshold", py_float(a.join_threshold)), ("join_min_nodes", a.join_min_nodes.to_string()),
                          ("join_head", a.join_head.to_string())]);
            if a.join_min_paths > 1 {
                stats.push(("join_min_paths", a.join_min_paths.to_string()));
            }
        }
        fs::write(stage.join("stats.tsv"),
                  format!("metric\tvalue\n{}", stats.iter().map(|(m, v)| format!("{}\t{}\n", m, v)).collect::<String>()))?;
        Ok(())
    })();
    match written {
        Ok(()) => {
            fs::rename(&stage, &a.output_dir)?;
            println!("{} repeat candidates, {} resolved, {} blocked; graph {} -> {} unitigs -> {}", decisions.len(),
                     done.len(), blocked.len(), before.0, g.units.len(), a.output_dir.display());
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

    const K: usize = 4;

    /// A circular genome A -> R -> B -> R -> A whose repeat R the graph cannot cross.
    fn ring_with_repeat() -> Graph {
        let mut g = Graph::new(K);
        g.add(0, vec![0, 1, 2], 30.0, false, None, None);
        g.add(1, vec![3, 4], 60.0, false, None, None);
        g.add(2, vec![5, 6, 7], 30.0, false, None, None);
        for (a, b) in [(0, 1), (1, 2), (2, 1), (1, 0)] {
            g.link(a, b);
        }
        g
    }

    fn triples(xs: &[((i64, i64, i64), u64)]) -> HashMap<(i64, i64, i64), u64> {
        xs.iter().copied().collect()
    }

    #[test]
    fn two_copy_repeat_closes_the_ring() {
        let mut g = ring_with_repeat();
        let d = resolve_candidates(&g, &triples(&[((0, 1, 2), 5), ((2, 1, 0), 4)]), 3, 0.2);
        assert_eq!(d.iter().map(|x| (x.r, x.matching.is_some())).collect::<Vec<_>>(), vec![(1, true)]);
        let (done, blocked) = apply(&mut g, &d, &BTreeSet::new());
        assert_eq!((done.len(), blocked.len()), (1, 0));
        compact(&mut g);
        assert_eq!(g.units.len(), 1);
        let (&u, x) = g.units.iter().next().unwrap();
        assert!(x.circular && g.succ(u).is_empty());
        let m = &x.members;
        let i = m.iter().position(|&v| v == 0).unwrap();
        let rotated: Vec<u32> = m[i..].iter().chain(m[..i].iter()).copied().collect();
        assert_eq!(rotated, vec![0, 1, 2, 3, 4, 5, 6, 7, 3, 4]);
        let mut origin = x.origin.clone();
        origin.sort();
        assert_eq!(origin, vec![0, 1, 1, 2]);
    }

    #[test]
    fn unresolved_cases_leave_the_graph_unchanged() {
        // support < 3, conflicting reads, not one-to-one
        for t in [triples(&[((0, 1, 2), 2), ((2, 1, 0), 4)]),
                  triples(&[((0, 1, 2), 5), ((2, 1, 0), 4), ((0, 1, 0), 2)]),
                  triples(&[((0, 1, 0), 5), ((2, 1, 0), 4)])] {
            let mut g = ring_with_repeat();
            let d = resolve_candidates(&g, &t, 3, 0.2);
            assert!(d[0].matching.is_none());
            let (done, _) = apply(&mut g, &d, &BTreeSet::new());
            compact(&mut g);
            assert!(done.is_empty());
            assert_eq!(g.units.keys().copied().collect::<Vec<_>>(), vec![0, 1, 2]);
        }
    }

    #[test]
    fn self_link_and_unequal_degrees_are_not_split() {
        let mut g = Graph::new(K);
        for u in 0..4 {
            g.add(u, vec![u as u32], 10.0, false, None, None);
        }
        for (a, b) in [(0, 1), (1, 1), (1, 2), (3, 1)] {
            g.link(a, b);
        }
        let d = resolve_candidates(&g, &triples(&[((0, 1, 2), 9)]), 3, 0.2);
        assert_eq!(d.iter().map(|x| (x.r, x.matching.is_some())).collect::<Vec<_>>(), vec![(1, false)]);
        let mut g2 = Graph::new(K);
        for u in 0..4 {
            g2.add(u, vec![u as u32], 10.0, false, None, None);
        }
        for (a, b) in [(0, 1), (3, 1), (1, 2)] {
            g2.link(a, b);
        }
        assert!(resolve_candidates(&g2, &HashMap::new(), 3, 0.2)[0].matching.is_none());
    }

    #[test]
    fn reads_map_to_unitigs_by_exact_windows() {
        let seq: Vec<u32> = (100..112).collect();
        let mut node_tokens: Vec<u32> = (0..=seq.len() - K).flat_map(|i| seq[i..i + K].to_vec()).collect();
        node_tokens.extend([7, 8, 9, 10]);
        let node_unitig: Vec<i64> = [0i64; 9].iter().copied().chain([1]).collect();
        let node_pos: Vec<i64> = (0..9).chain([0]).collect();
        let keyed = node_index(&node_tokens, K, &node_unitig);
        let mut read_a: Vec<u32> = seq[2..9].to_vec();
        read_a.extend([7, 8, 9, 10, 99]);
        let a = read_visit(&read_a, K, &keyed, &node_tokens, &node_unitig, &node_pos);
        assert_eq!(a, (vec![(0, 2, 5), (1, 0, 0)], 5, (12 - K + 1) as u64));
        let b = read_visit(&[1, 2, 3, 4, 5], K, &keyed, &node_tokens, &node_unitig, &node_pos);
        assert_eq!(b, (vec![], 0, 2));
    }

    #[test]
    fn rolling_hash_equals_the_direct_hash() {
        let tok: Vec<u32> = (0..50).map(|i| (i * 7919 % 101) as u32).collect();
        for k in 1..8 {
            let node_tokens: Vec<u32> = (0..=tok.len() - k).flat_map(|i| tok[i..i + k].to_vec()).collect();
            let n = node_tokens.len() / k;
            let node_unitig = vec![0i64; n];
            let node_pos: Vec<i64> = (0..n as i64).collect();
            let keyed = node_index(&node_tokens, k, &node_unitig);
            let (seq, matched, windows) = read_visit(&tok, k, &keyed, &node_tokens, &node_unitig, &node_pos);
            assert_eq!((matched, windows), (n as u64, n as u64));
            assert_eq!(seq.len(), 1);
        }
    }

    #[test]
    fn closed_single_tangle_needs_the_direct_link() {
        // A -> B directly and through x <-> y, z a dead end hanging from x.
        let mut g = Graph::new(K);
        g.add(0, vec![0, 1, 2], 30.0, false, None, None);
        g.add(1, vec![3], 300.0, false, None, None);
        g.add(2, vec![5], 300.0, false, None, None);
        g.add(3, vec![6], 5.0, false, None, None);
        g.add(4, vec![7, 8, 9], 30.0, false, None, None);
        for (a, b) in [(0, 4), (0, 1), (1, 2), (2, 1), (1, 3), (1, 4)] {
            g.link(a, b);
        }
        let anchors: HashSet<i64> = [0i64, 4].iter().copied().collect();
        assert_eq!(tangles(&g, &anchors).1.keys().copied().collect::<Vec<_>>(), vec![1]);
        g.unlink(0, 4);
        assert!(tangles(&g, &anchors).1.is_empty());
    }

    #[test]
    fn rounds_ties_to_even() {
        assert_eq!(round_even(2.5), 2);
        assert_eq!(round_even(3.5), 4);
    }
}
