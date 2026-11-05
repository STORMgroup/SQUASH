//! Unitigs and contigs: the compacted unitigs of a mask, the contigs an order passes to the next one, and the final
//! graph's unitig table and GFA.

use flate2::Compression;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, Write};
use std::path::Path;

use super::caps::{ResourceCaps, cap_error, enforce_rss_cap};
use super::files::{BackgroundGz, create_plain, finish_plain};
use super::state::{Contig, Edge, Node, State, Unitig};

const UNITIG_HEADER: &str =
    "k\tstate\tunitig_id\tcircular\tnodes\ttoken_count\tabundance2\tindegree\toutdegree\tmembers";

pub(super) fn abundance2(nodes: &[u32], atomic: &[Node]) -> u32 {
    let mut values: Vec<u32> = nodes
        .iter()
        .map(|&id| atomic[id as usize].abundance_in)
        .collect();
    let mid = values.len() / 2;
    if values.len() % 2 == 1 {
        let (_, median, _) = values.select_nth_unstable(mid);
        median.saturating_mul(2)
    } else {
        let (low, high, _) = values.select_nth_unstable(mid);
        let upper = *high;
        // `low` holds exactly the mid smallest values in arbitrary order.
        let lower = *low
            .iter()
            .max()
            .expect("even-length unitig has a lower half");
        lower.saturating_add(upper)
    }
}

/// Lossless node-partition compaction. Branch/merge vertices remain singleton
/// unitigs, while maximal 1-in/1-out interiors are compacted. Pure cycles are
/// emitted once, rotated to their lowest node ID.
pub(super) fn compact_unitigs(nodes: &[Node], edges: &[Edge], active: &[bool]) -> Vec<Unitig> {
    let n = nodes.len();
    const NONE: u32 = u32::MAX;
    // A compactable edge u->v is unique out of u and unique into v. Only that predicate and the unique neighbour are
    // needed, so the distinct degree is tracked with a first-seen value plus a "saw something else" flag instead of an
    // adjacency list per node (this runs several times per order over millions of nodes).
    let mut first_out = vec![NONE; n];
    let mut first_in = vec![NONE; n];
    let mut multi_out = vec![false; n];
    let mut multi_in = vec![false; n];
    for edge in edges {
        if !active[edge.from as usize] || !active[edge.to as usize] {
            continue;
        }
        let (from, to) = (edge.from as usize, edge.to as usize);
        if first_out[from] == NONE {
            first_out[from] = edge.to;
        } else if first_out[from] != edge.to {
            multi_out[from] = true;
        }
        if first_in[to] == NONE {
            first_in[to] = edge.from;
        } else if first_in[to] != edge.from {
            multi_in[to] = true;
        }
    }
    let mut next = vec![NONE; n];
    let mut prev = vec![NONE; n];
    for u in 0..n {
        if !active[u] || multi_out[u] || first_out[u] == NONE {
            continue;
        }
        let v = first_out[u] as usize;
        if !multi_in[v] {
            next[u] = first_out[u];
            prev[v] = u as u32;
        }
    }
    let mut visited = vec![false; n];
    let mut raw: Vec<(Vec<u32>, bool)> = Vec::new();
    for start in 0..n {
        if !active[start] || visited[start] || prev[start] != NONE {
            continue;
        }
        let mut path = Vec::new();
        let mut cur = start;
        loop {
            if visited[cur] {
                break;
            }
            visited[cur] = true;
            path.push(cur as u32);
            if next[cur] == NONE {
                break;
            }
            cur = next[cur] as usize;
        }
        raw.push((path, false));
    }
    // Remaining compactable components are pure cycles. Node IDs are
    // assigned by sorted exact keys, so the lowest key is also the exact
    // token-vector minimal rotation anchor -- which makes the lowest member ID
    // the anchor too, since ID order and key order coincide.
    debug_assert!(nodes.windows(2).all(|w| w[0].key < w[1].key));
    for start in 0..n {
        if !active[start] || visited[start] {
            continue;
        }
        let mut cycle = Vec::new();
        let mut cur = start;
        loop {
            if visited[cur] {
                break;
            }
            visited[cur] = true;
            cycle.push(cur as u32);
            debug_assert!(next[cur] != NONE, "cycle component must be closed");
            cur = next[cur] as usize;
        }
        if !cycle.is_empty() {
            let p = cycle
                .iter()
                .enumerate()
                .min_by_key(|(_, node)| **node)
                .map(|(i, _)| i)
                .unwrap();
            cycle.rotate_left(p);
            raw.push((cycle, true));
        }
    }
    // Unitigs partition the active nodes, so their first members are distinct
    // and already determine the lexicographic order of the full member vectors.
    raw.sort_by_key(|(members, _)| members[0]);
    let mut member_to = vec![u32::MAX; n];
    for (id, (members, _)) in raw.iter().enumerate() {
        for &node in members {
            member_to[node as usize] = id as u32;
        }
    }
    // Only the distinct cross-component neighbour counts are needed, from one flat sort/dedup.
    let mut crossings: Vec<(u32, u32)> = Vec::new();
    for edge in edges {
        if !active[edge.from as usize] || !active[edge.to as usize] {
            continue;
        }
        let from = member_to[edge.from as usize];
        let to = member_to[edge.to as usize];
        if from != to {
            crossings.push((from, to));
        }
    }
    crossings.sort_unstable();
    crossings.dedup();
    let mut indegree = vec![0u32; raw.len()];
    let mut outdegree = vec![0u32; raw.len()];
    for &(from, to) in &crossings {
        outdegree[from as usize] += 1;
        indegree[to as usize] += 1;
    }
    raw.into_iter()
        .enumerate()
        .map(|(id, (members, circular))| Unitig {
            id: id as u32,
            abundance2: abundance2(&members, nodes),
            nodes: members,
            circular,
            indegree: indegree[id],
            outdegree: outdegree[id],
        })
        .collect()
}

#[cfg(test)]
fn select_contigs(nodes: &[Node], edges: &[Edge], cleaned: &[bool], first: &[bool], second: &[bool]) -> Vec<Contig> {
    select_contigs_bounded(nodes, edges, cleaned, first, second, None).expect("uncapped contig selection cannot fail")
}

/// The contigs of an order: the unitigs after the second cutoff of median abundance at least 6, then those after the
/// first cutoff of at least 4, then those of the cleaned unfiltered graph (metaMDBG dumps each cutoff's unitigs only
/// after simplification), each mask as `select_contigs_in_turn`.
pub(super) fn select_contigs_bounded(
    nodes: &[Node],
    edges: &[Edge],
    cleaned: &[bool],
    first_cleaned: &[bool],
    second_cleaned: &[bool],
    caps: Option<&ResourceCaps>,
) -> io::Result<Vec<Contig>> {
    let phases: [(f64, &[bool]); 3] = [(6.0, second_cleaned), (4.0, first_cleaned), (1.0, cleaned)];
    select_contigs_in_turn(nodes, edges, &phases, caps)
}

/// The contigs of `phases` (minimum median abundance, mask), in order: each mask's unitigs, longest first,
/// are taken when their median (abundance2 / 2) is at least the minimum and they share no node with a taken one.
fn select_contigs_in_turn(
    nodes: &[Node],
    edges: &[Edge],
    phases: &[(f64, &[bool])],
    caps: Option<&ResourceCaps>,
) -> io::Result<Vec<Contig>> {
    let mut claimed = BTreeSet::new();
    let mut out = Vec::new();
    for &(minimum, mask) in phases {
        let mut rows = compact_unitigs(nodes, edges, mask);
        // Node IDs are assigned in ascending NodeKey order by both state constructors (`seed_state` and
        // `next_state_with` group sorted keys), so `nodes[i].key < nodes[j].key` exactly when `i < j`: comparing the
        // member ID vectors is the same as comparing their key vectors.
        debug_assert!(nodes.windows(2).all(|w| w[0].key < w[1].key));
        rows.sort_by(|a, b| {
            b.nodes
                .len()
                .cmp(&a.nodes.len())
                .then_with(|| b.abundance2.cmp(&a.abundance2))
                .then_with(|| a.circular.cmp(&b.circular))
                .then_with(|| a.nodes.cmp(&b.nodes))
                .then_with(|| a.id.cmp(&b.id))
        });
        for u in rows {
            let enough = u.abundance2 as f64 >= minimum * 2.0;
            let overlap = u.nodes.iter().any(|id| claimed.contains(id));
            if enough && !overlap {
                if let Some(caps) = caps {
                    if out.len() as u64 >= caps.max_contigs {
                        return Err(cap_error(format!(
                            "contigs={} exceeds {}",
                            out.len() as u64 + 1,
                            caps.max_contigs
                        )));
                    }
                }
                for &id in &u.nodes {
                    claimed.insert(id);
                }
                let id = out.len() as u64;
                out.push(Contig {
                    id,
                    nodes: u.nodes,
                    abundance2: u.abundance2,
                    circular: u.circular,
                });
            }
        }
        enforce_rss_cap(caps)?;
    }
    Ok(out)
}

pub(super) fn mask_from_unitigs(n: usize, unitigs: &[Unitig]) -> Vec<bool> {
    let mut mask = vec![false; n];
    for u in unitigs {
        for &id in &u.nodes {
            mask[id as usize] = true;
        }
    }
    mask
}

/// The unitig table of one mask (`UNITIG_HEADER`): one row per unitig, its nodes in `members`.
pub(super) fn write_unitig_table(path: &Path, k: u16, state_name: &str, unitigs: &[Unitig], level: Compression) -> io::Result<()> {
    let mut uw = BackgroundGz::create_level(path, level)?;
    writeln!(uw, "{}", UNITIG_HEADER)?;
    for u in unitigs {
        let tokens = k as usize + u.nodes.len().saturating_sub(1);
        writeln!(
            uw,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            k,
            state_name,
            u.id,
            usize::from(u.circular),
            u.nodes.len(),
            tokens,
            u.abundance2,
            u.indegree,
            u.outdegree,
            u.nodes
                .iter()
                .map(u32::to_string)
                .collect::<Vec<_>>()
                .join(",")
        )?;
    }
    uw.finish()
}

/// The unitig graph of one mask in GFA: a segment per unitig (LN = tokens, IA = rounded-up median abundance) and a
/// link per pair of unitigs joined by an edge of the mask, in (from, to) order.
pub(super) fn write_unitig_graph(path: &Path, state: &State, mask: &[bool], unitigs: &[Unitig]) -> io::Result<()> {
    let mut uw = create_plain(path)?;
    writeln!(uw, "H\tVN:Z:1.0")?;
    let mut member_to = BTreeMap::new();
    for u in unitigs {
        for &id in &u.nodes {
            member_to.insert(id, u.id);
        }
        writeln!(
            uw,
            "S\t{}\t*\tLN:i:{}\tIA:i:{}",
            u.id,
            state.k as usize + u.nodes.len().saturating_sub(1),
            (u.abundance2 + 1) / 2
        )?;
    }
    let mut links = BTreeSet::new();
    for e in &state.edges {
        if mask[e.from as usize] && mask[e.to as usize] {
            let a = member_to[&e.from];
            let b = member_to[&e.to];
            if a != b {
                links.insert((a, b));
            }
        }
    }
    for (a, b) in links {
        writeln!(uw, "L\t{}\t+\t{}\t+\t{}M", a, b, state.k - 1)?;
    }
    finish_plain(uw)
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::test_support::*;

    #[test]
    fn compacts_complete_linear_chain() {
        let nodes = (0..4).map(|i| node(i, 5)).collect::<Vec<_>>();
        let u = compact_unitigs(
            &nodes,
            &[edge(0, 1), edge(1, 2), edge(2, 3)],
            &vec![true; 4],
        );
        assert_eq!(u.len(), 1);
        assert_eq!(u[0].nodes, vec![0, 1, 2, 3]);
    }

    #[test]
    fn fork_and_merge_are_node_partitioned() {
        let nodes = (0..6).map(|i| node(i, 5)).collect::<Vec<_>>();
        let u = compact_unitigs(
            &nodes,
            &[
                edge(0, 1),
                edge(1, 2),
                edge(1, 3),
                edge(2, 4),
                edge(3, 4),
                edge(4, 5),
            ],
            &vec![true; 6],
        );
        let mut all = u.iter().flat_map(|x| x.nodes.clone()).collect::<Vec<_>>();
        all.sort_unstable();
        assert_eq!(all, (0..6).collect::<Vec<_>>());
        assert_eq!(u.iter().map(|x| x.nodes.len()).sum::<usize>(), 6);
    }

    #[test]
    fn self_loop_and_pure_cycle_emit_once() {
        let one = vec![node(0, 5)];
        let u = compact_unitigs(&one, &[edge(0, 0)], &[true]);
        assert_eq!(u.len(), 1);
        assert!(u[0].circular);
        let nodes = (0..3).map(|i| node(i, 5)).collect::<Vec<_>>();
        let u = compact_unitigs(
            &nodes,
            &[edge(0, 1), edge(1, 2), edge(2, 0)],
            &vec![true; 3],
        );
        assert_eq!(u.len(), 1);
        assert!(u[0].circular);
        assert_eq!(u[0].nodes, vec![0, 1, 2]);
    }

    #[test]
    fn contigs_follow_the_mask_order_and_take_whole_unitigs() {
        let nodes = vec![node(0, 7), node(1, 7)];
        let c = select_contigs(&nodes, &[edge(0, 1)], &[true, true], &[true, true], &[true, true]);
        assert_eq!(c.len(), 1);
        assert_eq!(c[0].nodes, vec![0, 1]);
    }
}
