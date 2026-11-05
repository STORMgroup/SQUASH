//! An order's graph: its nodes, edges, read paths, the previous order's contigs as paths, and its contigs.

use std::convert::TryFrom;

/// A node's identity: its k tokens at the seed order, the pair of previous-order nodes it joins above it.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub(super) enum NodeKey {
    Seed(Vec<u32>),
    Parents(u32, u32),
}

/// A node of an order's graph (a k-min-mer) with its abundance and where it was seen.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Node {
    pub(super) key: NodeKey,
    pub(super) abundance_in: u32,
    pub(super) abundance_out: u32,
    pub(super) direct_read_count: u32,
    pub(super) read_occurrences: u64,
    pub(super) contig_count: u32,
    pub(super) contig_occurrences: u64,
    pub(super) source_class: u8, // 1 reads, 2 contig paths, 3 both
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Edge {
    pub(super) from: u32,
    pub(super) to: u32,
    pub(super) read_count: u32,
    pub(super) read_occurrences: u64,
    pub(super) contig_paths: u32,
    pub(super) contig_occurrences: u64,
    pub(super) inferred_only: bool,
    pub(super) multiplicity: u64,
}

/// One step of a read's path through an order's nodes, packed in 16 bytes (read-path steps dominate an order's
/// state in memory): the node as u32 with `NO_NODE` for none, and the ordinal and q span as u32 (the accessors give
/// u64). A value outside u32 panics (`PathStep::new`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct PathStep {
    node_id: u32,
    ordinal32: u32,
    q_start32: u32,
    q_end32: u32,
}

const NO_NODE: u32 = u32::MAX;

impl PathStep {
    pub(super) fn new(node: Option<u32>, path_ordinal: u64, q_start: u64, q_end: u64) -> PathStep {
        let narrow = |x: u64, what: &str| -> u32 {
            u32::try_from(x).unwrap_or_else(|_| panic!("path step {} {} outside the compact layout (u32)", what, x))
        };
        let node_id = match node {
            Some(id) if id == NO_NODE => panic!("path step node id {} is the compact layout's empty marker", id),
            Some(id) => id,
            None => NO_NODE,
        };
        PathStep {
            node_id,
            ordinal32: narrow(path_ordinal, "ordinal"),
            q_start32: narrow(q_start, "q start"),
            q_end32: narrow(q_end, "q end"),
        }
    }
    pub(super) fn node(&self) -> Option<u32> {
        if self.node_id == NO_NODE { None } else { Some(self.node_id) }
    }
    pub(super) fn path_ordinal(&self) -> u64 {
        self.ordinal32 as u64
    }
    pub(super) fn q_start(&self) -> u64 {
        self.q_start32 as u64
    }
    pub(super) fn q_end(&self) -> u64 {
        self.q_end32 as u64
    }
}

/// A read's path through an order's nodes.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ReadPath {
    pub(super) read_ordinal: u64,
    pub(super) steps: Vec<PathStep>,
}

/// One step of a contig path: the node (if any) of the contig's tokens `token_start..token_end`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ContigStep {
    pub(super) node: Option<u32>,
    pub(super) token_start: u64,
    pub(super) token_end: u64,
}

/// A contig of the previous order as a path through this order's nodes (metaMDBG threads the contigs of each order
/// into the next one like reads).
#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ContigPath {
    pub(super) source_contig_id: u64,
    pub(super) circular: bool,
    /// Number of unique token starts in the circular source. Zero for a
    /// linear source. Circular contig path coordinates are unrolled against this
    /// period so the closing window remains representable without inventing
    /// a read.
    pub(super) token_period: u64,
    pub(super) steps: Vec<ContigStep>,
}

/// A contig of an order: a unitig selected from the cleaned masks (`select_contigs_bounded`), passed to the next order.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Contig {
    pub(super) id: u64,
    pub(super) nodes: Vec<u32>,
    pub(super) abundance2: u32,
    pub(super) circular: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct State {
    pub(super) k: u16,
    pub(super) nodes: Vec<Node>,
    pub(super) edges: Vec<Edge>,
    pub(super) read_paths: Vec<ReadPath>,
    pub(super) contig_paths: Vec<ContigPath>,
    pub(super) active: Vec<bool>,
    pub(super) contigs: Vec<Contig>,
    pub(super) previous_prefix: String,
    /// `assemble --tip-max-tokens`: the longest tip cut at every order, in tokens of the unitig (nodes + K - 1);
    /// None = 8 x K nodes.
    pub(super) tip_max_tokens: Option<u32>,
    /// `assemble --progressive-cutoffs`: the first and second abundance cutoffs instead of 2 and 3.
    pub(super) progressive_cutoffs: Option<(u32, u32)>,
    /// `assemble --superbubble-max-tokens`: the longest superbubble path popped at every order, in tokens (nodes +
    /// K - 1; metaMDBG pops up to 50 kb at every order); None = 8 x K nodes.
    pub(super) superbubble_max_tokens: Option<u32>,
}

impl State {
    /// Longest superbubble path, in nodes, popped at this order.
    pub(super) fn superbubble_max(&self) -> usize {
        match self.superbubble_max_tokens {
            Some(t) => (t as usize).saturating_sub(self.k as usize - 1).max(1),
            None => 8 * self.k as usize,
        }
    }
    /// The progressive filter's cutoffs (first, second): metaMDBG's 2 and 3 unless `progressive_cutoffs` is set.
    pub(super) fn cutoffs(&self) -> (u32, u32) {
        self.progressive_cutoffs.unwrap_or((2, 3))
    }
    /// Longest tip, in nodes, cut at this order.
    pub(super) fn tip_max(&self) -> usize {
        match self.tip_max_tokens {
            Some(t) => (t as usize).saturating_sub(self.k as usize - 1).max(1),
            None => 8 * self.k as usize,
        }
    }
}

#[derive(Clone, Debug)]
pub(super) struct Unitig {
    pub(super) id: u32,
    pub(super) nodes: Vec<u32>,
    pub(super) circular: bool,
    pub(super) abundance2: u32,
    pub(super) indegree: u32,
    pub(super) outdegree: u32,
}

pub(super) fn node_overlap_ok(left: &Node, right: &Node) -> bool {
    match (&left.key, &right.key) {
        (NodeKey::Seed(a), NodeKey::Seed(b)) => {
            a.len() == b.len() && a.len() >= 2 && a[1..] == b[..b.len() - 1]
        }
        (NodeKey::Parents(_, suffix), NodeKey::Parents(prefix, _)) => suffix == prefix,
        _ => false,
    }
}

pub(super) fn node_adjacencies(nodes: &[u32], circular: bool) -> Vec<(usize, u32, u32)> {
    let mut out: Vec<_> = nodes
        .windows(2)
        .enumerate()
        .map(|(rank, pair)| (rank, pair[0], pair[1]))
        .collect();
    if circular && !nodes.is_empty() {
        out.push((nodes.len() - 1, *nodes.last().unwrap(), nodes[0]));
    }
    out
}

pub(super) fn contig_path_adjacencies(path: &ContigPath) -> Vec<(usize, ContigStep, ContigStep)> {
    (0..contig_path_adjacency_count(path))
        .map(|rank| {
            let (left, right) = contig_path_adjacency_at(path, rank);
            (rank, left, right)
        })
        .collect()
}

/// Number of adjacencies `contig_path_adjacencies` would emit for `path`.
fn contig_path_adjacency_count(path: &ContigPath) -> usize {
    let windows = path.steps.len().saturating_sub(1);
    if path.circular && !path.steps.is_empty() {
        windows + 1
    } else {
        windows
    }
}

/// The adjacency at `rank`, with the circular closing pair projected onto the unrolled coordinate axis exactly as
/// `contig_path_adjacencies` does. An inverted index can then store adjacency ranks and rebuild only the pairs a query
/// touches, instead of materialising every step of every path.
fn contig_path_adjacency_at(path: &ContigPath, rank: usize) -> (ContigStep, ContigStep) {
    let windows = path.steps.len().saturating_sub(1);
    if rank < windows {
        return (path.steps[rank].clone(), path.steps[rank + 1].clone());
    }
    debug_assert!(path.circular && !path.steps.is_empty() && rank == path.steps.len() - 1);
    let left = path.steps[path.steps.len() - 1].clone();
    let mut right = path.steps[0].clone();
    right.token_start = right.token_start.saturating_add(path.token_period);
    right.token_end = right.token_end.saturating_add(path.token_period);
    (left, right)
}

/// A state is a seed state exactly when its nodes carry token-vector keys rather than parent pairs. An empty node set
/// falls back to "has no predecessor", the only other way a seed state is distinguishable.
pub(super) fn state_is_seed(state: &State) -> bool {
    state
        .nodes
        .first()
        .map(|node| matches!(node.key, NodeKey::Seed(_)))
        .unwrap_or(state.previous_prefix.is_empty())
}
