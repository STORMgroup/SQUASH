//! The graph of an order: the seed order from the reads' token windows, each next order from the previous one's nodes,
//! read paths and contigs, and the edges both infer.

use rayon::prelude::*;
use std::cmp::min;
use std::io;
use std::path::Path;

use super::caps::{ResourceCaps, cap_error, enforce_rss_cap};
use super::clock::PhaseClock;
use super::files::{invalid_data, invalid_input, usize_u32};
use super::state::{
    Contig, ContigPath, ContigStep, Edge, Node, NodeKey, PathStep, ReadPath, State, contig_path_adjacencies,
    node_adjacencies, node_overlap_ok,
};
use super::tokens::TokenRead;

fn inferred_edges(
    nodes: &[Node],
    read_paths: &[ReadPath],
    contig_paths: &[ContigPath],
) -> io::Result<Vec<Edge>> {
    inferred_edges_bounded(nodes, read_paths, contig_paths, None)
}

/// Exact overlap graph over the nodes, annotated with observed support.
///
/// Node IDs ascend with `NodeKey`, so the nodes are sorted by their `(k-1)`-prefix: each node's successors are a
/// contiguous range found by binary search, and walking `from` in ascending order emits the candidate pairs in
/// ascending `(from, to)`. Observed support is attached by sorting flat triples instead of keeping a per-edge set of
/// read IDs.
fn inferred_edges_bounded(
    nodes: &[Node],
    read_paths: &[ReadPath],
    contig_paths: &[ContigPath],
    caps: Option<&ResourceCaps>,
) -> io::Result<Vec<Edge>> {
    let mut read_obs: Vec<(u32, u32, u64)> = Vec::new();
    for path in read_paths {
        for pair in path.steps.windows(2) {
            if pair[0].path_ordinal() + 1 != pair[1].path_ordinal() {
                continue;
            }
            if let (Some(a), Some(b)) = (pair[0].node(), pair[1].node()) {
                read_obs.push((a, b, path.read_ordinal));
            }
        }
    }
    let mut contig_obs: Vec<(u32, u32, u64)> = Vec::new();
    for path in contig_paths {
        for (_, left, right) in contig_path_adjacencies(path) {
            if let (Some(a), Some(b)) = (left.node, right.node) {
                contig_obs.push((a, b, path.source_contig_id));
            }
        }
    }
    read_obs.par_sort_unstable();
    contig_obs.par_sort_unstable();

    let is_k20 = nodes
        .first()
        .map(|n| matches!(n.key, NodeKey::Seed(_)))
        .unwrap_or(false);
    // Successor range for `from`: the nodes whose (k-1)-prefix equals its
    // (k-1)-suffix.
    let successors = |from: usize| -> (usize, usize) {
        if is_k20 {
            let key = match &nodes[from].key {
                NodeKey::Seed(v) => v,
                _ => return (0, 0),
            };
            let suffix = &key[1..];
            let lo = nodes.partition_point(|node| match &node.key {
                NodeKey::Seed(v) => v[..v.len() - 1] < *suffix,
                _ => true,
            });
            let hi = nodes.partition_point(|node| match &node.key {
                NodeKey::Seed(v) => v[..v.len() - 1] <= *suffix,
                _ => true,
            });
            (lo, hi)
        } else {
            let suffix = match nodes[from].key {
                NodeKey::Parents(_, suffix) => suffix,
                _ => return (0, 0),
            };
            let lo = nodes.partition_point(|node| match node.key {
                NodeKey::Parents(prefix, _) => prefix < suffix,
                _ => true,
            });
            let hi = nodes.partition_point(|node| match node.key {
                NodeKey::Parents(prefix, _) => prefix <= suffix,
                _ => true,
            });
            (lo, hi)
        }
    };

    // Both observed lists and the emitted candidates are in ascending
    // (from, to), so support is attached by advancing cursors in lockstep.
    let mut read_cursor = 0usize;
    let mut contig_cursor = 0usize;
    let mut edges: Vec<Edge> = Vec::new();
    let mut emitted = 0u64;
    for from in 0..nodes.len() {
        let (lo, hi) = successors(from);
        for to in lo..hi {
            let (from_id, to_id) = (from as u32, to as u32);
            emitted += 1;
            if caps.map(|c| emitted > c.max_links).unwrap_or(false) {
                return Err(cap_error(format!(
                    "links>{} during exact edge expansion",
                    caps.unwrap().max_links
                )));
            }
            if edges.len() & 0xffff == 0 {
                enforce_rss_cap(caps)?;
            }
            if !node_overlap_ok(&nodes[from], &nodes[to]) {
                return Err(invalid_data(
                    "compact parent-pair overlap consistency failure",
                ));
            }
            while read_cursor < read_obs.len()
                && (read_obs[read_cursor].0, read_obs[read_cursor].1) < (from_id, to_id)
            {
                return Err(invalid_data(format!(
                    "observed path edge {:?} is absent from exact overlap graph",
                    (read_obs[read_cursor].0, read_obs[read_cursor].1)
                )));
            }
            while contig_cursor < contig_obs.len()
                && (contig_obs[contig_cursor].0, contig_obs[contig_cursor].1) < (from_id, to_id)
            {
                return Err(invalid_data(format!(
                    "observed path edge {:?} is absent from exact overlap graph",
                    (contig_obs[contig_cursor].0, contig_obs[contig_cursor].1)
                )));
            }
            let mut read_occurrences = 0u64;
            let mut read_count = 0u32;
            let mut last_read: Option<u64> = None;
            while read_cursor < read_obs.len()
                && (read_obs[read_cursor].0, read_obs[read_cursor].1) == (from_id, to_id)
            {
                read_occurrences = read_occurrences
                    .checked_add(1)
                    .ok_or_else(|| invalid_data("biological edge occurrence overflow"))?;
                if last_read != Some(read_obs[read_cursor].2) {
                    read_count += 1;
                    last_read = Some(read_obs[read_cursor].2);
                }
                read_cursor += 1;
            }
            let mut contig_occurrences = 0u64;
            let mut contig_path_count = 0u32;
            let mut last_contig: Option<u64> = None;
            while contig_cursor < contig_obs.len()
                && (contig_obs[contig_cursor].0, contig_obs[contig_cursor].1) == (from_id, to_id)
            {
                contig_occurrences = contig_occurrences
                    .checked_add(1)
                    .ok_or_else(|| invalid_data("contig path edge occurrence overflow"))?;
                if last_contig != Some(contig_obs[contig_cursor].2) {
                    contig_path_count += 1;
                    last_contig = Some(contig_obs[contig_cursor].2);
                }
                contig_cursor += 1;
            }
            let inferred = read_count == 0 && contig_path_count == 0;
            let multiplicity = if inferred {
                1
            } else {
                read_occurrences
                    .checked_add(contig_occurrences)
                    .ok_or_else(|| invalid_data("edge multiplicity overflow"))?
            };
            edges.push(Edge {
                from: from_id,
                to: to_id,
                read_count,
                read_occurrences,
                contig_paths: contig_path_count,
                contig_occurrences,
                inferred_only: inferred,
                multiplicity,
            });
        }
    }
    if read_cursor < read_obs.len() {
        return Err(invalid_data(format!(
            "observed path edge {:?} is absent from exact overlap graph",
            (read_obs[read_cursor].0, read_obs[read_cursor].1)
        )));
    }
    if contig_cursor < contig_obs.len() {
        return Err(invalid_data(format!(
            "observed path edge {:?} is absent from exact overlap graph",
            (contig_obs[contig_cursor].0, contig_obs[contig_cursor].1)
        )));
    }
    Ok(edges)
}

/// The seed graph: every k-token window of the reads seen in at least two distinct reads, numbered in ascending
/// token-vector order, with each read's complete window path.
pub(super) fn seed_state(reads: &[TokenRead], k: u16) -> io::Result<State> {
    let width = k as usize;
    if width < 2 {
        return Err(invalid_input("seed order must be at least 2"));
    }
    let total_tokens: usize = reads.iter().map(|read| read.tokens.len()).sum();
    let mut flat: Vec<u32> = Vec::with_capacity(total_tokens);
    let mut token_base: Vec<u32> = Vec::with_capacity(reads.len());
    for read in reads {
        token_base.push(usize_u32(flat.len(), "prepared token offset")?);
        flat.extend(read.tokens.iter().map(|token| token.hash as u32));
    }

    // One entry per k-mer window: where its tokens start, and which read and
    // path ordinal it came from.
    let mut window_base: Vec<u32> = Vec::with_capacity(reads.len());
    let mut total_windows = 0usize;
    for read in reads {
        window_base.push(usize_u32(total_windows, "window index")?);
        total_windows += read.tokens.len().saturating_sub(width - 1);
    }
    let mut entries: Vec<(u32, u32, u32)> = Vec::with_capacity(total_windows);
    for (read_index, read) in reads.iter().enumerate() {
        let read_index = usize_u32(read_index, "read index")?;
        let count = read.tokens.len().saturating_sub(width - 1);
        for offset in 0..count {
            entries.push((
                token_base[read_index as usize] + usize_u32(offset, "window offset")?,
                read_index,
                usize_u32(offset, "path ordinal")?,
            ));
        }
    }

    let window = |start: u32| -> &[u32] { &flat[start as usize..start as usize + width] };
    entries.par_sort_unstable_by(|a, b| {
        window(a.0)
            .cmp(window(b.0))
            .then_with(|| (a.1, a.2).cmp(&(b.1, b.2)))
    });

    // Group equal windows. Groups appear in ascending token-vector order, so assigning IDs in group order numbers the
    // nodes in ascending key order.
    let mut nodes: Vec<Node> = Vec::new();
    let mut node_of_window: Vec<u32> = vec![u32::MAX; total_windows];
    let mut group_start = 0usize;
    while group_start < entries.len() {
        let mut group_end = group_start + 1;
        while group_end < entries.len()
            && window(entries[group_end].0) == window(entries[group_start].0)
        {
            group_end += 1;
        }
        let occurrences = (group_end - group_start) as u64;
        // Entries in a group are ordered by (read, ordinal), so distinct
        // reads are counted by transitions instead of a per-node set.
        let mut distinct = 1u32;
        for pair in entries[group_start..group_end].windows(2) {
            if pair[0].1 != pair[1].1 {
                distinct += 1;
            }
        }
        if distinct >= 2 {
            let id = usize_u32(nodes.len(), "atomic node id")?;
            for &(_, read_index, ordinal) in &entries[group_start..group_end] {
                let index = window_base[read_index as usize] as usize + ordinal as usize;
                node_of_window[index] = id;
            }
            nodes.push(Node {
                key: NodeKey::Seed(window(entries[group_start].0).to_vec()),
                abundance_in: distinct,
                abundance_out: distinct,
                direct_read_count: distinct,
                read_occurrences: occurrences,
                contig_count: 0,
                contig_occurrences: 0,
                source_class: 1,
            });
        }
        group_start = group_end;
    }
    drop(entries);

    let mut read_paths = Vec::with_capacity(reads.len());
    for (read_index, read) in reads.iter().enumerate() {
        let count = read.tokens.len().saturating_sub(width - 1);
        let base = window_base[read_index] as usize;
        let mut steps = Vec::with_capacity(count);
        for offset in 0..count {
            let id = node_of_window[base + offset];
            steps.push(PathStep::new(if id == u32::MAX { None } else { Some(id) }, offset as u64, read.tokens[offset].start as u64, read.tokens[offset + width - 1].end as u64));
        }
        read_paths.push(ReadPath {
            read_ordinal: read.ordinal,
            steps,
        });
    }
    let contig_paths = Vec::new();
    let edges = inferred_edges(&nodes, &read_paths, &contig_paths)?;
    Ok(State {
        k,
        nodes,
        edges,
        read_paths,
        contig_paths,
        active: Vec::new(),
        contigs: Vec::new(),
        previous_prefix: String::new(),
        tip_max_tokens: None,
        progressive_cutoffs: None,
        superbubble_max_tokens: None,
    })
}

#[cfg(test)]
fn next_state(previous: &State, previous_prefix: &Path) -> io::Result<State> {
    next_state_with(previous, previous_prefix, None)
}

pub(super) fn next_state_with(previous: &State, previous_prefix: &Path, caps: Option<&ResourceCaps>) -> io::Result<State> {
    // The schedule bound (k_max) is the caller's; this is the representational limit of the order counter itself.
    if previous.k == u16::MAX {
        return Err(invalid_input("cannot advance beyond the maximum order"));
    }
    enforce_rss_cap(caps)?;
    let mut clock = PhaseClock::new("next_state");
    let k = previous.k + 1;
    let next_path_rows: u64 = previous
        .read_paths
        .iter()
        .map(|path| path.steps.len().saturating_sub(1) as u64)
        .sum();
    if let Some(caps) = caps {
        if next_path_rows > caps.max_path_rows {
            return Err(cap_error(format!(
                "path_rows={} exceeds {}",
                next_path_rows, caps.max_path_rows
            )));
        }
    }
    // Each path emits flat (prefix, suffix, source) triples; sorting them reduces in parallel and yields the ascending
    // (prefix, suffix) grouping.
    let admissible = |a: u32, b: u32| -> bool {
        previous.nodes[a as usize].abundance_out >= 2
            && previous.nodes[b as usize].abundance_out >= 2
            && node_overlap_ok(&previous.nodes[a as usize], &previous.nodes[b as usize])
    };
    let partials: Vec<Vec<(u32, u32, u64)>> = {
        previous
            .read_paths
            .par_iter()
            .map(|path| {
                let mut local = Vec::new();
                for pair in path.steps.windows(2) {
                    if pair[0].path_ordinal() + 1 != pair[1].path_ordinal() {
                        continue;
                    }
                    if let (Some(a), Some(b)) = (pair[0].node(), pair[1].node()) {
                        if admissible(a, b) {
                            local.push((a, b, path.read_ordinal));
                        }
                    }
                }
                local
            })
            .collect()
    };
    let mut read_obs: Vec<(u32, u32, u64)> = Vec::new();
    // The resident-memory cap is checked every 4,096 paths and once after (one /proc/self/status read per path is slow
    // over millions of paths).
    for (index, partial) in partials.into_iter().enumerate() {
        read_obs.extend_from_slice(&partial);
        if index & 0xfff == 0 {
            enforce_rss_cap(caps)?;
        }
    }
    enforce_rss_cap(caps)?;
    clock.mark("observations");
    read_obs.par_sort_unstable();
    clock.mark("observations sorted");

    // Whether a previous contig's path keeps its closing link: a circular contig of one node (a period-1
    // cycle) is passed on without it, as metaMDBG's ToMinspace wraps a circular contig only when it has more than k
    // minimizers (with the link the loop would grow by one token per order with no read holding the run).
    let wraps = |contig: &Contig| contig.circular && contig.nodes.len() >= 2;
    let mut contig_obs: Vec<(u32, u32, u64)> = Vec::new();
    for contig in &previous.contigs {
        for (_, a, b) in node_adjacencies(&contig.nodes, wraps(contig)) {
            if admissible(a, b) {
                contig_obs.push((a, b, contig.id));
            }
        }
    }
    contig_obs.par_sort_unstable();

    /// One grouped candidate: the parent pair, its occurrence count and the
    /// number of distinct sources that produced it.
    fn group_pairs(observed: &[(u32, u32, u64)]) -> Vec<((u32, u32), u64, u32)> {
        let mut out: Vec<((u32, u32), u64, u32)> = Vec::new();
        let mut index = 0usize;
        while index < observed.len() {
            let key = (observed[index].0, observed[index].1);
            let mut occurrences = 0u64;
            let mut distinct = 0u32;
            let mut last: Option<u64> = None;
            while index < observed.len() && (observed[index].0, observed[index].1) == key {
                occurrences += 1;
                if last != Some(observed[index].2) {
                    distinct += 1;
                    last = Some(observed[index].2);
                }
                index += 1;
            }
            out.push((key, occurrences, distinct));
        }
        out
    }
    let read_groups = group_pairs(&read_obs);
    let contig_groups = group_pairs(&contig_obs);
    drop(read_obs);
    drop(contig_obs);

    // Merge the two ascending key lists into their ascending union.
    let mut nodes: Vec<Node> = Vec::with_capacity(read_groups.len() + contig_groups.len());
    let mut lookup: Vec<(u64, u32)> = Vec::with_capacity(nodes.capacity());
    let (mut bi, mut ci) = (0usize, 0usize);
    while bi < read_groups.len() || ci < contig_groups.len() {
        let key = match (read_groups.get(bi), contig_groups.get(ci)) {
            (Some(b), Some(c)) => min(b.0, c.0),
            (Some(b), None) => b.0,
            (None, Some(c)) => c.0,
            (None, None) => unreachable!("loop guard guarantees a remaining group"),
        };
        let (read_occ, read_distinct) = match read_groups.get(bi) {
            Some(b) if b.0 == key => {
                bi += 1;
                (b.1, b.2)
            }
            _ => (0, 0),
        };
        let (contig_occ, contig_distinct) = match contig_groups.get(ci) {
            Some(c) if c.0 == key => {
                ci += 1;
                (c.1, c.2)
            }
            _ => (0, 0),
        };
        let (a, b) = key;
        let id = usize_u32(nodes.len(), "atomic node id")?;
        if let Some(caps) = caps {
            if nodes.len() as u64 + 1 > caps.max_nodes {
                return Err(cap_error(format!(
                    "nodes>{} during candidate reduction",
                    caps.max_nodes
                )));
            }
        }
        let abundance = min(
            previous.nodes[a as usize].abundance_out,
            previous.nodes[b as usize].abundance_out,
        );
        let source = (if read_occ > 0 { 1 } else { 0 }) | (if contig_occ > 0 { 2 } else { 0 });
        lookup.push(((a as u64) << 32 | b as u64, id));
        nodes.push(Node {
            key: NodeKey::Parents(a, b),
            abundance_in: abundance,
            abundance_out: abundance,
            direct_read_count: read_distinct,
            read_occurrences: read_occ,
            contig_count: contig_distinct,
            contig_occurrences: contig_occ,
            source_class: source,
        });
    }
    enforce_rss_cap(caps)?;
    clock.mark("groups and nodes");
    // A flat sorted table beats a tree here: every retained path step resolves
    // its successor node through this, tens of millions of times per order.
    let keys = |a: u32, b: u32| -> Option<u32> {
        let probe = (a as u64) << 32 | b as u64;
        lookup
            .binary_search_by(|entry| entry.0.cmp(&probe))
            .ok()
            .map(|index| lookup[index].1)
    };
    let read_paths: Vec<ReadPath> = {
        previous
            .read_paths
            .par_iter()
            .map(|path| {
                let mut steps = Vec::with_capacity(path.steps.len().saturating_sub(1));
                for pair in path.steps.windows(2) {
                    let node = if pair[0].path_ordinal() + 1 == pair[1].path_ordinal() {
                        match (pair[0].node(), pair[1].node()) {
                            (Some(a), Some(b)) => keys(a, b),
                            _ => None,
                        }
                    } else {
                        None
                    };
                    steps.push(PathStep::new(node, pair[0].path_ordinal(), pair[0].q_start(), pair[1].q_end()));
                }
                ReadPath {
                    read_ordinal: path.read_ordinal,
                    steps,
                }
            })
            .collect()
    };
    let contig_paths: Vec<ContigPath> = previous
        .contigs
        .par_iter()
        .map(|contig| {
            let steps = node_adjacencies(&contig.nodes, wraps(contig))
                .into_iter()
                .map(|(i, a, b)| ContigStep {
                    node: keys(a, b),
                    token_start: i as u64,
                    token_end: (i + k as usize) as u64,
                })
                .collect();
            ContigPath {
                source_contig_id: contig.id,
                circular: wraps(contig),
                token_period: if wraps(contig) { contig.nodes.len() as u64 } else { 0 },
                steps,
            }
        })
        .collect();
    enforce_rss_cap(caps)?;
    clock.mark("paths");
    let edges = inferred_edges_bounded(&nodes, &read_paths, &contig_paths, caps)?;
    clock.mark("edges");
    Ok(State {
        k,
        nodes,
        edges,
        read_paths,
        contig_paths,
        active: Vec::new(),
        contigs: Vec::new(),
        previous_prefix: previous_prefix.display().to_string(),
        tip_max_tokens: previous.tip_max_tokens,
        progressive_cutoffs: previous.progressive_cutoffs,
        superbubble_max_tokens: previous.superbubble_max_tokens,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::test_support::*;
    use super::super::caps::cap_reason;

    #[test]
    fn prefix_suffix_recurrence_uses_min_and_distinct_rc() {
        let mut prev = state(vec![overlap_node(0, 7), overlap_node(1, 3)],
            Vec::new(),
        );
        let steps = vec![
            PathStep::new(Some(0), 0, 0, 20),
            PathStep::new(Some(1), 1, 1, 21),
        ];
        prev.read_paths = vec![
            ReadPath {
                read_ordinal: 9,
                steps: steps.clone(),
            },
            ReadPath {
                read_ordinal: 9,
                steps: steps.clone(),
            },
            ReadPath {
                read_ordinal: 10,
                steps,
            },
        ];
        let next = next_state(&prev, Path::new("prev")).unwrap();
        assert_eq!(next.nodes.len(), 1);
        assert_eq!(next.nodes[0].abundance_in, 3);
        assert_eq!(next.nodes[0].direct_read_count, 2);
        assert_eq!(next.nodes[0].read_occurrences, 3);
        assert_eq!(next.nodes[0].contig_occurrences, 0);
    }

    #[test]
    fn a_contig_only_node_has_zero_read_abundance_and_parent_min() {
        let mut prev = state(vec![overlap_node(0, 8), overlap_node(1, 4)],
            Vec::new(),
        );
        prev.contigs.push(Contig {
            id: 4,
            nodes: vec![0, 1],
            abundance2: 8,
            circular: false,
        });
        let next = next_state(&prev, Path::new("prev")).unwrap();
        assert_eq!(next.nodes.len(), 1);
        let n = &next.nodes[0];
        assert_eq!(
            (
                n.abundance_in,
                n.direct_read_count,
                n.read_occurrences,
                n.contig_count,
                n.source_class
            ),
            (4, 0, 0, 1, 2)
        );
    }

    #[test]
    fn abundance_one_or_missing_parent_is_excluded() {
        let mut prev = state(vec![overlap_node(0, 5), overlap_node(1, 1)],
            Vec::new(),
        );
        prev.read_paths.push(ReadPath {
            read_ordinal: 0,
            steps: vec![
                PathStep::new(Some(0), 0, 0, 20),
                PathStep::new(Some(1), 1, 1, 21),
                PathStep::new(None, 2, 2, 22),
            ],
        });
        assert!(next_state(&prev, Path::new("prev"))
            .unwrap()
            .nodes
            .is_empty());
    }

    #[test]
    fn observed_edge_multiplicity_counts_occurrences() {
        let nodes = vec![overlap_node(0, 3), overlap_node(1, 3)];
        let steps = vec![
            PathStep::new(Some(0), 0, 0, 20),
            PathStep::new(Some(1), 1, 1, 21),
        ];
        let edges = inferred_edges(
            &nodes,
            &[
                ReadPath {
                    read_ordinal: 0,
                    steps: steps.clone(),
                },
                ReadPath {
                    read_ordinal: 0,
                    steps,
                },
            ],
            &[],
        )
        .unwrap();
        assert_eq!(
            edges
                .iter()
                .find(|x| x.from == 0 && x.to == 1)
                .unwrap()
                .multiplicity,
            2
        );
    }

    fn seed_node(tokens: Vec<u32>) -> Node {
        Node {
            key: NodeKey::Seed(tokens),
            abundance_in: 5,
            abundance_out: 5,
            direct_read_count: 5,
            read_occurrences: 5,
            contig_count: 0,
            contig_occurrences: 0,
            source_class: 1,
        }
    }

    fn circular_contig(nodes: Vec<u32>) -> Contig {
        Contig {
            id: 17,
            nodes,
            abundance2: 10,
            circular: true,
        }
    }

    #[test]
    fn a_one_node_cycle_reaches_the_next_order_without_its_closing_link() {
        // metaMDBG's ToMinspace wraps a circular contig only with more than k minimizers: a one-node
        // cycle (one token repeated) is not grown into the next order from the contig alone.
        let mut previous = state(vec![seed_node(vec![7; 20])], vec![edge(0, 0)]);
        previous.contigs.push(circular_contig(vec![0]));
        let next = next_state(&previous, Path::new("k20")).unwrap();
        assert_eq!(next.nodes.len(), 0);
        assert_eq!(next.edges.len(), 0);
        assert_eq!(next.contig_paths.len(), 1);
        assert!(!next.contig_paths[0].circular);
        assert_eq!(next.contig_paths[0].token_period, 0);
        assert!(next.contig_paths[0].steps.is_empty());
    }

    #[test]
    fn a_longer_cycle_keeps_its_closing_link() {
        let a: Vec<u32> = (0..20).map(|i| 1 + (i % 2)).collect(); // 1,2,1,2,...
        let b: Vec<u32> = (0..20).map(|i| 2 - (i % 2)).collect(); // 2,1,2,1,...
        let mut previous = state(vec![seed_node(a), seed_node(b)], vec![edge(0, 1), edge(1, 0)]);
        previous.contigs.push(circular_contig(vec![0, 1]));
        let new = next_state(&previous, Path::new("k20")).unwrap();
        assert_eq!(new.contig_paths.len(), 1);
        assert!(new.contig_paths[0].circular);
        assert_eq!(new.contig_paths[0].token_period, 2);
        assert_eq!(new.nodes.len(), 2);
    }

    #[test]
    fn recurrence_stops_before_path_or_node_cap_allocation() {
        let mut previous = state(vec![overlap_node(0, 3), overlap_node(1, 3)],
            vec![edge(0, 1)],
        );
        previous.read_paths = (0..3u64)
            .map(|read_ordinal| ReadPath {
                read_ordinal,
                steps: vec![
                    PathStep::new(Some(0), 0, 0, 20),
                    PathStep::new(Some(1), 1, 1, 21),
                ],
            })
            .collect();
        let mut caps = generous_caps();
        caps.max_path_rows = 2;
        let error =
            next_state_with(&previous, Path::new("k20"), Some(&caps)).unwrap_err();
        assert!(cap_reason(&error).unwrap().starts_with("path_rows="));
    }

    #[test]
    fn parallel_recurrence_is_identical_at_1_4_16_threads() {
        let mut previous = state(vec![overlap_node(0, 6), overlap_node(1, 5), overlap_node(2, 4)],
            vec![edge(0, 1), edge(1, 2)],
        );
        for read_ordinal in 0..8u64 {
            previous.read_paths.push(ReadPath {
                read_ordinal,
                steps: (0..3)
                    .map(|path_ordinal| PathStep::new(Some(path_ordinal as u32), path_ordinal, path_ordinal, path_ordinal + 20))
                    .collect(),
            });
        }
        previous.contigs.push(Contig {
            id: 5,
            nodes: vec![0, 1, 2],
            abundance2: 8,
            circular: false,
        });
        let states = [1usize, 4, 16]
            .iter()
            .map(|&threads| {
                // `run` installs one pool around the whole dispatch, so the
                // recurrence inherits whatever pool is ambient.  Install the
                // width under test here to keep exercising it.
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap()
                    .install(|| {
                        next_state_with(&previous, Path::new("k20"), None).unwrap()
                    })
            })
            .collect::<Vec<_>>();
        assert_eq!(states[0], states[1]);
        assert_eq!(states[0], states[2]);
    }
}
