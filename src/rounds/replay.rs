//! One read's replay: the donors' votes into the target's layers (alleles and transitions), the heaviest path
//! (`content_path`) and the corrected tokens on the target's axis.

use super::*;

/// `(inserted tokens before this layer, token at this layer)`; `None` is DEL.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Allele {
    pub insertion: Vec<u32>,
    pub site: Option<u32>,
}

impl Allele {
    fn emitted(&self) -> impl DoubleEndedIterator<Item = u32> + '_ {
        self.insertion.iter().copied().chain(self.site)
    }
    fn emitted_len(&self) -> usize {
        self.insertion.len() + self.site.is_some() as usize
    }
}

/// One donor row: the overlap's q-symbol bounds on the target and on the donor.
#[derive(Clone, Debug)]
pub struct DonorRow {
    pub donor: usize,
    pub target_start: i64,
    pub target_end: i64,
    pub donor_start: i64,
    pub donor_end: i64,
}

/// Per-layer votes. `alleles` (with their support) in order of first sight: the
/// target's own allele first, then donor alleles.
#[derive(Clone, Debug, Default)]
pub struct Layer {
    pub alleles: Vec<(Allele, i64)>,
    /// `(allele id at the previous layer, allele id here) -> support`.
    pub edges: Edges,
}

/// A layer's edge supports as a sorted vector `(left, right, support)`: a layer holds a few edges, so this is smaller
/// and faster than a hash map per layer. Looked up by key only.
#[derive(Clone, Debug, Default)]
pub struct Edges(Vec<(u32, u32, u32)>);

impl Edges {
    fn find(&self, left: u32, right: u32) -> std::result::Result<usize, usize> {
        self.0.binary_search_by(|&(l, r, _)| (l, r).cmp(&(left, right)))
    }
    /// One more read through the edge.
    pub fn add(&mut self, left: u32, right: u32) {
        match self.find(left, right) {
            Ok(k) => self.0[k].2 += 1,
            Err(k) => self.0.insert(k, (left, right, 1)),
        }
    }
    #[cfg(test)]
    pub fn get(&self, left: u32, right: u32) -> i64 {
        self.find(left, right).map_or(0, |k| self.0[k].2 as i64)
    }
}

impl Layer {
    fn intern(&mut self, allele: Allele) -> u32 {
        if let Some(k) = self.alleles.iter().position(|(a, _)| *a == allele) {
            self.alleles[k].1 += 1;
            return k as u32;
        }
        self.alleles.push((allele, 1));
        (self.alleles.len() - 1) as u32
    }
}

#[derive(Clone, Debug)]
pub struct Replay {
    pub layers: Vec<Layer>,
    /// Chosen allele id per layer (id 0 is always the target's own allele).
    pub chosen: Vec<u32>,
}

impl Replay {
    #[cfg(test)]
    pub fn output_values(&self) -> Vec<u32> {
        self.chosen_alleles().flat_map(|a| a.emitted()).collect()
    }
    fn chosen_alleles(&self) -> impl Iterator<Item = &Allele> + '_ {
        self.layers.iter().zip(&self.chosen).map(|(l, &c)| &l.alleles[c as usize].0)
    }
    pub fn changed_layers(&self) -> usize {
        self.chosen.iter().filter(|&&c| c != 0).count()
    }
}

/// First and last index of tokens overlapping `[start, end)`.
fn overlapping(tokens: &[Token], start: i64, end: i64) -> Option<(usize, usize)> {
    let hit = |t: &Token| t.e() > start && t.s() < end;
    Some((tokens.iter().position(hit)?, tokens.iter().rposition(hit)?))
}

/// ceil(a * b / c) on ints, through a correctly rounded f64 division.
pub(super) fn ceil_ratio(a: i64, b: i64, c: i64) -> i64 {
    ((a * b) as f64 / c as f64).ceil() as i64
}

/// One read's replay: every donor row aligned to the target and voted into its layers, then the heaviest path
/// (`content_path`); see the module documentation.
pub fn replay_read(
    target: &[Token],
    target_length: i64,
    target_symbols: i64,
    rows: &[DonorRow],
    tokens: &[Vec<Token>],
    lengths: &[i64],
    symbols: &[i64],
    pad: usize,
    band: usize,
) -> Replay {
    let mut layers: Vec<Layer> = target
        .iter()
        .map(|t| Layer { alleles: vec![(Allele { insertion: vec![], site: Some(t.value) }, 1)], edges: Edges::default() })
        .collect();
    for layer in layers.iter_mut().skip(1) {
        layer.edges.add(0, 0);
    }
    let margin = Some(2 * pad);
    let target_values: Vec<u32> = target.iter().map(|t| t.value).collect();
    for row in rows {
        let donor = &tokens[row.donor];
        let ts = row.target_start * target_length / target_symbols;
        let te = ceil_ratio(row.target_end, target_length, target_symbols);
        let ds = row.donor_start * lengths[row.donor] / symbols[row.donor];
        let de = ceil_ratio(row.donor_end, lengths[row.donor], symbols[row.donor]);
        let (ti, di) = match (overlapping(target, ts, te), overlapping(donor, ds, de)) {
            (Some(ti), Some(di)) => (ti, di),
            _ => continue,
        };
        let t0 = ti.0.saturating_sub(pad);
        let t1 = (ti.1 + 1 + pad).min(target.len());
        let d0 = di.0.saturating_sub(pad);
        let d1 = (di.1 + 1 + pad).min(donor.len());
        let donor_values: Vec<u32> = donor[d0..d1].iter().map(|t| t.value).collect();
        let full = if band == 0 {
            align_overlap(&target_values[t0..t1], &donor_values, margin)
        } else {
            align_overlap_banded(&target_values[t0..t1], &donor_values, margin, band)
        };
        let core = aligned_core(&full);
        if core.is_empty() {
            continue;
        }
        // Absolute target layer -> allele id; the core covers a contiguous target range.
        let mut previous: Option<(usize, u32)> = None;
        let mut pending: Vec<u32> = Vec::new();
        for &(a, b) in core {
            let a = match a {
                None => {
                    pending.push(donor_values[b.unwrap()]);
                    continue;
                }
                Some(a) => a,
            };
            let absolute = t0 + a;
            let site = b.map(|b| donor_values[b]);
            let id = layers[absolute].intern(Allele { insertion: std::mem::take(&mut pending), site });
            if let Some((left, left_id)) = previous {
                if absolute == left + 1 {
                    layers[absolute].edges.add(left_id, id);
                }
            }
            previous = Some((absolute, id));
        }
    }
    let chosen = content_path(&layers);
    Replay { layers, chosen }
}

/// A DP path as a persistent list: one node per (layer, chosen allele).
struct PathNode {
    parent: Option<usize>,
    layer: usize,
    allele: u32,
    /// Emitted tokens on the path up to and including this node.
    length: usize,
}

struct Paths<'a> {
    layers: &'a [Layer],
    nodes: Vec<PathNode>,
}

impl<'a> Paths<'a> {
    fn allele(&self, node: usize) -> &'a Allele {
        let n = &self.nodes[node];
        &self.layers[n.layer].alleles[n.allele as usize].0
    }
    fn length(&self, node: Option<usize>) -> usize {
        node.map_or(0, |n| self.nodes[n].length)
    }
    /// Move a reverse cursor past exhausted nodes; `remaining` counts tokens of `node` not yet read.
    fn settle(&self, cursor: &mut (Option<usize>, usize)) {
        while let (Some(node), 0) = *cursor {
            let parent = self.nodes[node].parent;
            *cursor = (parent, parent.map_or(0, |p| self.allele(p).emitted_len()));
        }
    }
    fn token(&self, cursor: &mut (Option<usize>, usize)) -> u32 {
        let node = cursor.0.unwrap();
        cursor.1 -= 1;
        self.allele(node).emitted().nth(cursor.1).unwrap()
    }
    /// `suffix_greater`: compare emitted tokens backwards; a longer path wins a full match.
    fn suffix_greater(&self, x: Option<usize>, y: Option<usize>) -> bool {
        let start = |node: Option<usize>| (node, node.map_or(0, |n| self.allele(n).emitted_len()));
        let (mut cx, mut cy) = (start(x), start(y));
        loop {
            self.settle(&mut cx);
            self.settle(&mut cy);
            if cx.0.is_none() || cy.0.is_none() {
                break;
            }
            if cx == cy {
                // Same remaining sequence on both: no difference left, and equal lengths.
                return false;
            }
            let (a, b) = (self.token(&mut cx), self.token(&mut cy));
            if a != b {
                return a > b;
            }
        }
        self.length(x) > self.length(y)
    }
}

/// `content_path`: heaviest target-layer path (node plus transition support),
/// ties resolved by emitted token content compared backwards.
pub fn content_path(layers: &[Layer]) -> Vec<u32> {
    let mut paths = Paths { layers, nodes: Vec::new() };
    // (score, path node) per allele at the previous layer, in allele order.
    let mut previous: Vec<(i64, Option<usize>)> = vec![(0, None)];
    // The layer's edge supports as a (previous allele x allele) table: the values of `Layer::edge` without a search
    // per pair.
    let mut table: Vec<i64> = Vec::new();
    for (position, layer) in layers.iter().enumerate() {
        let width = layer.alleles.len();
        table.clear();
        table.resize(previous.len() * width, 0);
        if position > 0 {
            for &(left, right, support) in &layer.edges.0 {
                table[left as usize * width + right as usize] = support as i64;
            }
        }
        let mut current = Vec::with_capacity(layer.alleles.len());
        for (id, (allele, support)) in layer.alleles.iter().enumerate() {
            let mut best: Option<(i64, Option<usize>)> = None;
            for (pid, &(before, node)) in previous.iter().enumerate() {
                let edge = table[pid * width + id];
                let score = before + support + edge;
                let wins = match best {
                    None => true,
                    Some((b, bnode)) => score > b || (score == b && paths.suffix_greater(node, bnode)),
                };
                if wins {
                    best = Some((score, node));
                }
            }
            let (score, parent) = best.unwrap();
            let length = paths.length(parent) + allele.emitted_len();
            paths.nodes.push(PathNode { parent, layer: position, allele: id as u32, length });
            current.push((score, Some(paths.nodes.len() - 1)));
        }
        previous = current;
    }
    let top = previous.iter().map(|x| x.0).max().unwrap_or(0);
    let mut best: Option<Option<usize>> = None;
    for &(score, node) in &previous {
        if score == top && best.map_or(true, |b| paths.suffix_greater(node, b)) {
            best = Some(node);
        }
    }
    let mut chosen = vec![0u32; layers.len()];
    let mut node = best.flatten();
    while let Some(n) = node {
        chosen[paths.nodes[n].layer] = paths.nodes[n].allele;
        node = paths.nodes[n].parent;
    }
    chosen
}

/// The corrected tokens on the target's axis: a layer's token keeps the target token's span; inserted tokens are
/// spread between the previous token's start and the layer's.
pub fn spanned_tokens(target: &[Token], replay: &Replay, source_length: i64) -> Vec<Token> {
    spanned_from_alleles(target, replay.chosen_alleles(), source_length)
}

/// `spanned_tokens` for any chosen allele per target layer.
pub(super) fn spanned_from_alleles<'a>(target: &[Token], alleles: impl Iterator<Item = &'a Allele>, source_length: i64) -> Vec<Token> {
    let mut result: Vec<Token> = Vec::new();
    for (layer, allele) in target.iter().zip(alleles) {
        let width = (layer.e() - layer.s()).max(1);
        let low = result.last().map_or((layer.s() - width).max(0), |t| t.s());
        let count = allele.insertion.len() as i64;
        for (index, &value) in allele.insertion.iter().enumerate() {
            let start = (source_length - 1).min(low + ((layer.s() - low) * (index as i64 + 1)).div_euclid(count + 1));
            result.push(Token::new(value, start, source_length.min((start + 1).max(start + width))));
        }
        if let Some(site) = allele.site {
            result.push(layer.with_value(site));
        }
    }
    // no spare capacity: the round's tokens stay in memory as the next round's input
    result.shrink_to_fit();
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn edges_match_a_hash_map() {
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut next = |m: u64| -> u64 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % m
        };
        for _ in 0..200 {
            let (mut edges, mut map) = (Edges::default(), HashMap::new());
            let ids = 1 + next(6);
            for _ in 0..next(60) {
                let (l, r) = (next(ids) as u32, next(ids) as u32);
                edges.add(l, r);
                *map.entry((l, r)).or_insert(0i64) += 1;
            }
            for l in 0..ids as u32 + 1 {
                for r in 0..ids as u32 + 1 {
                    assert_eq!(edges.get(l, r), *map.get(&(l, r)).unwrap_or(&0));
                }
            }
        }
    }

    /// A test cohort: unit spans, every read uses all others as donors over whole reads.
    fn outputs(reads: &[&[u32]], pad: usize) -> Vec<Vec<u32>> {
        let tokens: Vec<Vec<Token>> = reads.iter().map(|r| unit(r)).collect();
        let lengths: Vec<i64> = reads.iter().map(|r| r.len() as i64).collect();
        (0..reads.len())
            .map(|r| {
                let rows: Vec<DonorRow> = (0..reads.len())
                    .filter(|&d| d != r)
                    .map(|d| DonorRow { donor: d, target_start: 0, target_end: lengths[r], donor_start: 0, donor_end: lengths[d] })
                    .collect();
                replay_read(&tokens[r], lengths[r], lengths[r], &rows, &tokens, &lengths, &lengths, pad, 0).output_values()
            })
            .collect()
    }

    #[test]
    fn donor_end_is_not_a_deletion() {
        let tokens = vec![unit(&[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]), unit(&[1, 2, 3, 4, 5, 6])];
        let rows = [DonorRow { donor: 1, target_start: 0, target_end: 10, donor_start: 0, donor_end: 6 }];
        let replay = replay_read(&tokens[0], 10, 10, &rows, &tokens, &[10, 6], &[10, 6], 2, 0);
        assert!(replay.layers[9].alleles.len() == 1);
        assert_eq!(replay.layers[5].alleles[0].1, 2);
    }

    #[test]
    fn content_tie_break_is_frame_independent() {
        let reads: [&[u32]; 4] = [&[1, 2, 3, 4, 5, 6, 7], &[1, 2, 3, 4, 5, 6, 7], &[1, 2, 9, 5, 6, 7], &[1, 2, 9, 5, 6, 7]];
        let out = outputs(&reads, 8);
        assert!(out.iter().all(|x| *x == out[0]));
    }

    #[test]
    fn suffix_comparison_uses_the_most_recent_difference() {
        // Three single-token layers per path: (0, 3, 5) vs (9, 2, 5), and (7, 3) vs (3,).
        let layer = |values: &[u32]| Layer {
            alleles: values.iter().map(|&v| (Allele { insertion: vec![], site: Some(v) }, 1)).collect(),
            edges: Edges::default(),
        };
        let layers = vec![layer(&[0, 9, 7]), layer(&[3, 2]), layer(&[5])];
        let mut paths = Paths { layers: &layers, nodes: Vec::new() };
        let mut push = |parent: Option<usize>, layer: usize, allele: u32, len: usize| {
            paths.nodes.push(PathNode { parent, layer, allele, length: len });
            paths.nodes.len() - 1
        };
        let a0 = push(None, 0, 0, 1);
        let a1 = push(Some(a0), 1, 0, 2);
        let a2 = push(Some(a1), 2, 0, 3);
        let b0 = push(None, 0, 1, 1);
        let b1 = push(Some(b0), 1, 1, 2);
        let b2 = push(Some(b1), 2, 0, 3);
        let c0 = push(None, 0, 2, 1);
        let c1 = push(Some(c0), 1, 0, 2);
        let d = push(None, 1, 0, 1);
        assert!(paths.suffix_greater(Some(a2), Some(b2)));
        assert!(!paths.suffix_greater(Some(b2), Some(a2)));
        assert!(paths.suffix_greater(Some(c1), Some(d)));
        assert!(!paths.suffix_greater(Some(d), Some(c1)));
        assert!(!paths.suffix_greater(Some(a2), Some(a2)));
        assert!(paths.suffix_greater(Some(d), None));
    }

    #[test]
    fn balanced_site_goes_to_the_larger_recent_token() {
        // Balanced 1-vs-1 site with equal content lengths: the backwards content comparison decides
        // (larger token wins), so both reads pick the same allele whatever their own token is.
        let reads: [&[u32]; 2] = [&[1, 2, 3], &[1, 7, 3]];
        let out = outputs(&reads, 8);
        assert_eq!(out, vec![vec![1, 7, 3], vec![1, 7, 3]]);
    }

    #[test]
    fn spanned_tokens_keep_real_coordinates() {
        let target = vec![Token::new(10, 0, 5), Token::new(20, 10, 15), Token::new(30, 20, 25)];
        let layer = |alleles: Vec<Allele>| Layer { alleles: alleles.into_iter().map(|a| (a, 1)).collect(), edges: Edges::default() };
        let replay = Replay {
            layers: vec![
                layer(vec![Allele { insertion: vec![], site: Some(10) }]),
                layer(vec![Allele { insertion: vec![], site: Some(20) }, Allele { insertion: vec![99], site: Some(21) }]),
                layer(vec![Allele { insertion: vec![], site: Some(30) }, Allele { insertion: vec![], site: None }]),
            ],
            chosen: vec![0, 1, 1],
        };
        let tokens = spanned_tokens(&target, &replay, 30);
        assert_eq!(tokens.iter().map(|t| t.value).collect::<Vec<_>>(), vec![10, 99, 21]);
        assert_eq!((tokens[0].s(), tokens[0].e()), (0, 5));
        assert_eq!((tokens[2].s(), tokens[2].e()), (10, 15));
        assert!(0 <= tokens[1].s() && tokens[1].s() < 10);
        assert_eq!((tokens[1].s(), tokens[1].e()), (5, 10));
    }

    #[test]
    fn inserted_spans_floor_like_python_when_starts_go_backwards() {
        // The previous emitted token starts after this layer: the division floors toward -inf.
        let target = vec![Token::new(1, 10, 12), Token::new(2, 7, 9)];
        let layer = |a: Allele| Layer { alleles: vec![(a, 1)], edges: Edges::default() };
        let replay = Replay {
            layers: vec![layer(Allele { insertion: vec![], site: Some(1) }), layer(Allele { insertion: vec![5, 6], site: Some(2) })],
            chosen: vec![0, 0],
        };
        let tokens = spanned_tokens(&target, &replay, 30);
        // low = 10; starts: 10 + floor(-3 * 1 / 3) = 9, 10 + floor(-3 * 2 / 3) = 8.
        assert_eq!(tokens.iter().map(|t| (t.s(), t.e())).collect::<Vec<_>>(), vec![(10, 12), (9, 11), (8, 10), (7, 9)]);
    }

    #[test]
    fn window_bounds_match_python_rounding() {
        assert_eq!(ceil_ratio(7, 3, 2), 11);
        assert_eq!(ceil_ratio(4, 3, 2), 6);
        assert_eq!(ceil_ratio(0, 3, 2), 0);
    }
}
