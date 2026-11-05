//! Fixtures the graph's tests share.

use super::caps::ResourceCaps;
use super::state::{Edge, Node, NodeKey, State};

pub(super) fn generous_caps() -> ResourceCaps {
    ResourceCaps {
        max_rss_bytes: u64::MAX,
        max_nodes: u64::MAX,
        max_links: u64::MAX,
        max_path_rows: u64::MAX,
        max_contigs: u64::MAX,
        max_contig_token_ratio: u64::MAX,
        max_scratch_bytes: u64::MAX,
    }
}

pub(super) fn node(id: u32, abundance: u32) -> Node {
    let mut v = vec![0; 20];
    v[19] = id;
    Node {
        key: NodeKey::Seed(v),
        abundance_in: abundance,
        abundance_out: abundance,
        direct_read_count: abundance,
        read_occurrences: abundance as u64,
        contig_count: 0,
        contig_occurrences: 0,
        source_class: 1,
    }
}

pub(super) fn overlap_node(start: u32, abundance: u32) -> Node {
    Node {
        key: NodeKey::Seed((start..start + 20).collect()),
        abundance_in: abundance,
        abundance_out: abundance,
        direct_read_count: abundance,
        read_occurrences: abundance as u64,
        contig_count: 0,
        contig_occurrences: 0,
        source_class: 1,
    }
}

pub(super) fn edge(a: u32, b: u32) -> Edge {
    Edge {
        from: a,
        to: b,
        read_count: 1,
        read_occurrences: 1,
        contig_paths: 0,
        contig_occurrences: 0,
        inferred_only: false,
        multiplicity: 1,
    }
}

pub(super) fn state(nodes: Vec<Node>, edges: Vec<Edge>) -> State {
    State {
        k: 20,
        active: vec![true; nodes.len()],
        nodes,
        edges,
        read_paths: Vec::new(),
        contig_paths: Vec::new(),
        contigs: Vec::new(),
        previous_prefix: String::new(),
        tip_max_tokens: None,
        progressive_cutoffs: None,
        superbubble_max_tokens: None,
    }
}

/// Two variant sites closer than K: 0 -> {1, 2} -> {3, 4} -> 5, every branch linked to both next branches.
pub(super) fn crossed_bubble() -> (Vec<Node>, Vec<Edge>) {
    let nodes = vec![node(0, 30), node(1, 30), node(2, 8), node(3, 10), node(4, 30), node(5, 30)];
    let edges = vec![edge(0, 1), edge(0, 2), edge(1, 3), edge(1, 4), edge(2, 3), edge(2, 4), edge(3, 5), edge(4, 5)];
    (nodes, edges)
}
