//! Graph cleaning: tips, bubbles and superbubbles, and the progressive abundance filters.

use std::collections::{BTreeMap, BTreeSet};

use super::clock::PhaseClock;
use super::contigs::compact_unitigs;
use super::state::{Edge, Node, State};

/// metaMDBG's progressive abundance filter on one mask: unitigs whose median abundance is at most `threshold` are removed,
/// and the graph recompacted, until none is (metaMDBG: the cutoffs 2 and 3; `assemble --progressive-cutoffs` sets
/// others).
fn progressive_filter(nodes: &[Node], edges: &[Edge], start: &[bool], threshold: u32) -> Vec<bool> {
    let mut active = start.to_vec();
    let mut round = 0u32;
    loop {
        let unitigs = compact_unitigs(nodes, edges, &active);
        let mut remove = Vec::new();
        for u in &unitigs {
            if u.abundance2 <= threshold.saturating_mul(2) {
                remove.extend_from_slice(&u.nodes);
            }
        }
        if remove.is_empty() {
            break;
        }
        for id in remove {
            active[id as usize] = false;
        }
        round += 1;
        if round > nodes.len() as u32 {
            return active;
        }
    }
    active
}

/// Graph cleaning on one mask, as metaMDBG v1.2 simplifies its graph: tips, simple bubbles and superbubbles are removed
/// and the graph recompacted until nothing changes.
///
/// * Tip (metaMDBG's `TipRemover`): every non-circular unitig of at most `tip_max` nodes with no successor (and
///   a predecessor), or no predecessor (and a successor), in order of node count, abundance2 and first member, as long
///   as each neighbour it hangs from keeps another link this round or has it as its only successor (predecessor) from
///   the start of the round; otherwise it is retried after recompaction.
/// * Bubble: two or more simple unitigs (one predecessor `a`, one successor `c`, at most `bubble_max` nodes) that share
///   `a` and `c`; all but the strongest branch are removed, where strength orders by (abundance2, node count, lowest
///   first member).
/// * Superbubble: see `find_superbubble` (paths up to `superbubble_max` nodes). The unitigs off its widest path
///   (`heaviest_path`) are removed, unless the source links straight to the exit or the bubble touches a unitig already
///   removed this round (it is then retried on the recompacted graph).
///
/// Removed unitigs stay in the unfiltered graph, so a later contig selection step can still take them if nothing else
/// claims their nodes. `cutoff` names the mask (0 = the unfiltered graph, where a superbubble source may have at most
/// five successors; 2 = after the first abundance cutoff; 3 = after the second).
fn clean_graph(
    nodes: &[Node],
    edges: &[Edge],
    start: &[bool],
    tip_max: usize,
    bubble_max: usize,
    superbubble_max: usize,
    cutoff: u32,
) -> Vec<bool> {
    let mut active = start.to_vec();
    let mut round = 0u32;
    loop {
        let unitigs = compact_unitigs(nodes, edges, &active);
        let mut member_to = vec![u32::MAX; nodes.len()];
        for u in &unitigs {
            for &id in &u.nodes {
                member_to[id as usize] = u.id;
            }
        }
        let mut succ: Vec<BTreeSet<u32>> = vec![BTreeSet::new(); unitigs.len()];
        let mut pred: Vec<BTreeSet<u32>> = vec![BTreeSet::new(); unitigs.len()];
        for edge in edges {
            if !active[edge.from as usize] || !active[edge.to as usize] {
                continue;
            }
            let (a, b) = (member_to[edge.from as usize], member_to[edge.to as usize]);
            if a != b {
                succ[a as usize].insert(b);
                pred[b as usize].insert(a);
            }
        }
        let strength = |id: u32| {
            let u = &unitigs[id as usize];
            (u.abundance2, u.nodes.len(), std::cmp::Reverse(u.nodes[0]))
        };
        let mut removed: BTreeSet<u32> = BTreeSet::new();
        let mut out_live: Vec<usize> = succ.iter().map(|x| x.len()).collect();
        let mut in_live: Vec<usize> = pred.iter().map(|x| x.len()).collect();
        let mut tips: Vec<u32> = unitigs
            .iter()
            .filter(|u| {
                let i = u.id as usize;
                !u.circular
                    && u.nodes.len() <= tip_max
                    && ((succ[i].is_empty() && !pred[i].is_empty()) || (pred[i].is_empty() && !succ[i].is_empty()))
            })
            .map(|u| u.id)
            .collect();
        tips.sort_by_key(|&id| {
            let u = &unitigs[id as usize];
            (u.nodes.len(), u.abundance2, u.nodes[0])
        });
        for id in tips {
            let i = id as usize;
            let dead_end = succ[i].is_empty();
            // A dead end that is its predecessors' only successor from the start of the round is cut too (a dead end
            // reached from several unitigs); a sibling left alone by an earlier cut this round still waits, as metaMDBG
            // first recompacts it with the predecessor.
            let hangs_free = if dead_end {
                pred[i].iter().all(|&p| out_live[p as usize] >= 2 || succ[p as usize].len() == 1)
            } else {
                succ[i].iter().all(|&s| in_live[s as usize] >= 2 || pred[s as usize].len() == 1)
            };
            if !hangs_free {
                continue;
            }
            removed.insert(id);
            if dead_end {
                for &p in &pred[i] {
                    out_live[p as usize] -= 1;
                }
            } else {
                for &s in &succ[i] {
                    in_live[s as usize] -= 1;
                }
            }
        }
        for a in 0..unitigs.len() {
            let mut groups: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
            for &b in &succ[a] {
                let bi = b as usize;
                if b as usize == a || unitigs[bi].circular || unitigs[bi].nodes.len() > bubble_max {
                    continue;
                }
                if pred[bi].len() == 1 && succ[bi].len() == 1 {
                    let c = *succ[bi].iter().next().unwrap();
                    if c as usize != a && c != b {
                        groups.entry(c).or_default().push(b);
                    }
                }
            }
            for branches in groups.values() {
                if branches.len() < 2 {
                    continue;
                }
                let keep = *branches.iter().max_by_key(|&&b| strength(b)).unwrap();
                for &b in branches {
                    if b != keep {
                        removed.insert(b);
                    }
                }
            }
        }
        let lengths: Vec<usize> = unitigs.iter().map(|u| u.nodes.len()).collect();
        for source in 0..unitigs.len() as u32 {
            if succ[source as usize].len() < 2 || removed.contains(&source) {
                continue;
            }
            let (exit, inner) = match find_superbubble(source, &succ, &pred, &lengths, superbubble_max, cutoff == 0) {
                Some(found) => found,
                None => continue,
            };
            if succ[source as usize].contains(&exit)
                || removed.contains(&exit)
                || inner.iter().any(|x| removed.contains(x))
            {
                continue;
            }
            let keep = heaviest_path(
                source,
                exit,
                &inner,
                &pred,
                |id| strength(id),
                |id| unitigs[id as usize].abundance2 as u64,
                |id| unitigs[id as usize].nodes.len() as u64,
            );
            for &b in &inner {
                if !keep.contains(&b) {
                    removed.insert(b);
                }
            }
        }
        if removed.is_empty() {
            break;
        }
        for &id in &removed {
            for &node in &unitigs[id as usize].nodes {
                active[node as usize] = false;
            }
        }
        round += 1;
        if round > nodes.len() as u32 {
            break;
        }
    }
    active
}

/// Superbubble search from `source` on the unitig graph, as metaMDBG v1.2's
/// `isSuperbubble` (after Onodera et al. 2013): a depth-first walk that enters
/// a unitig only once all its predecessors are visited. It fails on a dead end,
/// a cycle, an edge back to the source, a path longer than `limit` nodes,
/// or, on the unfiltered graph (`unfiltered`), a unitig with more than five
/// successors. Returns the exit and the unitigs strictly inside.
fn find_superbubble(
    source: u32,
    succ: &[BTreeSet<u32>],
    pred: &[BTreeSet<u32>],
    lengths: &[usize],
    limit: usize,
    unfiltered: bool,
) -> Option<(u32, Vec<u32>)> {
    let mut visited = BTreeSet::new();
    let mut seen = BTreeSet::new();
    let mut length: BTreeMap<u32, usize> = BTreeMap::new();
    length.insert(source, 0);
    let mut stack = vec![source];
    while let Some(&v) = stack.last() {
        if unfiltered && succ[v as usize].len() > 5 {
            return None;
        }
        stack.pop();
        if length[&v] > limit {
            return None;
        }
        visited.insert(v);
        seen.remove(&v);
        if succ[v as usize].is_empty() {
            return None;
        }
        for &u in &succ[v as usize] {
            if u == source || visited.contains(&u) {
                return None;
            }
            seen.insert(u);
            let base = if v == source { 0 } else { length[&v] };
            length.insert(u, base + lengths[u as usize]);
        }
        for &u in &succ[v as usize] {
            if pred[u as usize].iter().all(|p| visited.contains(p)) {
                stack.push(u);
            }
            if stack.len() == 1 && seen.len() == 1 && seen.contains(&stack[0]) {
                let exit = stack[0];
                if succ[exit as usize].contains(&source) {
                    return None;
                }
                return Some((exit, visited.into_iter().filter(|&x| x != source).collect()));
            }
        }
    }
    None
}

/// The widest source-to-exit path through a superbubble: the path that maximizes (weakest inner unitig weight, sum of
/// weight x length over inner unitigs), found by a topological pass (the bubble is acyclic; metaMDBG takes the largest
/// sum of unitig abundances instead, which rewards a path split into more unitigs, while the weakest unitig bounds the
/// reads that can follow the whole path). The pair is exact as a DP key because after any unitig two partial paths
/// differ only through it. Equal keys go to the stronger predecessor. Returns the inner unitigs on that path.
fn heaviest_path<F: Fn(u32) -> S, S: Ord>(
    source: u32,
    exit: u32,
    inner: &[u32],
    pred: &[BTreeSet<u32>],
    strength: F,
    weight: impl Fn(u32) -> u64,
    length: impl Fn(u32) -> u64,
) -> BTreeSet<u32> {
    let members: BTreeSet<u32> = inner.iter().copied().chain([exit]).collect();
    // (key, predecessor); key = (weakest, length-weighted sum).
    let mut best: BTreeMap<u32, ((u64, u64), u32)> = BTreeMap::new();
    let mut pending: BTreeSet<u32> = members.clone();
    while !pending.is_empty() {
        let ready: Vec<u32> = pending
            .iter()
            .copied()
            .filter(|&v| pred[v as usize].iter().all(|p| *p == source || best.contains_key(p)))
            .collect();
        assert!(!ready.is_empty(), "superbubble is acyclic");
        for v in ready {
            let score = |p: u32| if p == source { (u64::MAX, 0) } else { best[&p].0 };
            let from = *pred[v as usize]
                .iter()
                .max_by(|&&a, &&b| (score(a), strength(a)).cmp(&(score(b), strength(b))))
                .unwrap();
            let (x, y) = score(from);
            let key = if v == exit { (x, y) } else { (x.min(weight(v)), y + weight(v) * length(v)) };
            best.insert(v, (key, from));
            pending.remove(&v);
        }
    }
    let mut keep = BTreeSet::new();
    let mut v = best[&exit].1;
    while v != source {
        keep.insert(v);
        v = best[&v].1;
    }
    keep
}

/// One cleaning of `mask` (named by `cutoff`, see `clean_graph`) with the state's limits: tips up to `tip_max`, bubble
/// branches up to 8 x K nodes, superbubbles up to `superbubble_max` (`cleaned_masks`).
fn clean_mask(state: &State, mask: &[bool], cutoff: u32) -> Vec<bool> {
    clean_graph(&state.nodes, &state.edges, mask, state.tip_max(), 8 * state.k as usize, state.superbubble_max(), cutoff)
}

/// The three cleaned masks: all nodes cleaned; that filtered by the first abundance cutoff and cleaned; that filtered
/// by the second cutoff and cleaned (the order's graph).
pub(super) fn cleaned_masks(state: &State) -> (Vec<bool>, Vec<bool>, Vec<bool>) {
    let all = vec![true; state.nodes.len()];
    let (c2, c3) = state.cutoffs();
    let mut clock = PhaseClock::new("cleaned_masks");
    let cleaned = clean_mask(state, &all, 0);
    clock.mark("clean");
    let first_cut = progressive_filter(&state.nodes, &state.edges, &cleaned, c2);
    clock.mark("first cut");
    let first_cleaned = clean_mask(state, &first_cut, 2);
    clock.mark("clean first cut");
    let second_cut = progressive_filter(&state.nodes, &state.edges, &first_cleaned, c3);
    clock.mark("second cut");
    let second_cleaned = clean_mask(state, &second_cut, 3);
    clock.mark("clean second cut");
    (cleaned, first_cleaned, second_cleaned)
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::test_support::*;
    use super::super::orders::process_state;

    #[test]
    fn progressive_is_simultaneous_global_cutoff() {
        let nodes = vec![node(0, 2), node(1, 2), node(2, 5)];
        let mask = progressive_filter(&nodes, &[edge(0, 1)], &vec![true; 3], 2);
        assert_eq!(mask, vec![false, false, true]);
    }

    #[test]
    fn progressive_threshold_is_separate_from_the_mask_label() {
        // Abundance-2 unitigs go at metaMDBG's cutoff 2 but stay at threshold 1.
        let nodes = vec![node(0, 2), node(1, 2), node(2, 5)];
        let mask = progressive_filter(&nodes, &[edge(0, 1)], &vec![true; 3], 1);
        assert_eq!(mask, vec![true, true, true]);
        let mask = progressive_filter(&nodes, &[edge(0, 1)], &vec![true; 3], 5);
        assert_eq!(mask, vec![false, false, false]);
    }

    #[test]
    fn a_wider_superbubble_limit_pops_a_bubble_the_default_keeps() {
        let (nodes, edges) = crossed_bubble();
        let mut st = state(nodes.clone(), edges.clone());
        st.k = 2;
        // Default 8 x 2 = 16 nodes pops it; a 2-token limit (1 node) does not (tips off: limit 0).
        let all = vec![true; 6];
        assert_eq!(clean_graph(&nodes, &edges, &all, 0, 10, st.superbubble_max(), 3), vec![true, true, false, false, true, true]);
        st.superbubble_max_tokens = Some(2);
        assert_eq!(clean_graph(&nodes, &edges, &all, 0, 10, st.superbubble_max(), 3), all);
    }

    #[test]
    fn clean_graph_pops_the_weaker_bubble_branch() {
        let nodes = vec![node(0, 5), node(1, 5), node(2, 2), node(3, 5)];
        let edges = vec![edge(0, 1), edge(0, 2), edge(1, 3), edge(2, 3)];
        let mask = clean_graph(&nodes, &edges, &vec![true; 4], 0, 10, 10, 3);
        assert_eq!(mask, vec![true, true, false, true]);
    }

    #[test]
    fn clean_graph_cuts_a_weak_tip_but_keeps_the_main_path() {
        // 0 -> 1 -> 2 is the main path; 1 -> 3 is a weak dead end. Tip limit 1: the 2-node unitig 0-1 is no tip, and
        // after 3 goes 2 is its only successor.
        let nodes = vec![node(0, 5), node(1, 5), node(2, 5), node(3, 2)];
        let edges = vec![edge(0, 1), edge(1, 2), edge(1, 3)];
        let mask = clean_graph(&nodes, &edges, &vec![true; 4], 1, 10, 10, 2);
        assert_eq!(mask, vec![true, true, true, false]);
    }

    #[test]
    fn tip_max_tokens_cuts_a_tip_longer_than_the_default_limit() {
        // Main path 0 -> 1 -> 2 -> ... -> 61 (strong, 62 nodes); a weak dead end 62..=81 (20 nodes) hangs off 40. At
        // K = 2 tips up to 8K = 16 nodes are cut, so it stays; with --tip-max-tokens 40 (39 nodes) it goes, and its
        // 21-node sibling 41..=61 then stays as the only successor of 0..=40 (41 nodes, no tip).
        let mut nodes: Vec<Node> = (0..62).map(|i| node(i, 20)).collect();
        let mut edges: Vec<Edge> = (1..62).map(|i| edge(i - 1, i)).collect();
        edges.push(edge(40, 62));
        for i in 62..82 {
            nodes.push(node(i, 6));
            if i > 62 {
                edges.push(edge(i - 1, i));
            }
        }
        let mut plain = state(nodes, edges);
        plain.k = 2;
        let mut long = plain.clone();
        long.tip_max_tokens = Some(40);
        assert_eq!((plain.tip_max(), long.tip_max()), (16, 39));
        let a = process_state(plain);
        let b = process_state(long);
        assert!(a.active.iter().all(|x| *x));
        assert_eq!(b.active, (0..82).map(|i| i < 62).collect::<Vec<_>>());
        // The limit never drops below one node, whatever K.
        let mut tiny = state(vec![node(0, 1)], Vec::new());
        tiny.k = 16;
        tiny.tip_max_tokens = Some(3);
        assert_eq!(tiny.tip_max(), 1);
    }

    #[test]
    fn clean_graph_pops_a_crossed_superbubble() {
        let (nodes, edges) = crossed_bubble();
        // the simple-bubble rule does not see it (each branch has two successors); tips off (limit 0)
        assert_eq!(clean_graph(&nodes, &edges, &vec![true; 6], 0, 10, 0, 3), vec![true; 6]);
        assert_eq!(clean_graph(&nodes, &edges, &vec![true; 6], 0, 10, 10, 3), vec![true, true, false, false, true, true]);
    }

    #[test]
    fn superbubble_search_finds_the_exit_and_the_inside() {
        let (nodes, edges) = crossed_bubble();
        let unitigs = compact_unitigs(&nodes, &edges, &vec![true; 6]);
        assert_eq!(unitigs.len(), 6);
        let (mut succ, mut pred) = (vec![BTreeSet::new(); 6], vec![BTreeSet::new(); 6]);
        for e in &edges {
            succ[e.from as usize].insert(e.to);
            pred[e.to as usize].insert(e.from);
        }
        let lengths = vec![1; 6];
        assert_eq!(find_superbubble(0, &succ, &pred, &lengths, 10, false), Some((5, vec![1, 2, 3, 4])));
        // Too long for the limit: the path 0 -> 1 -> 3 has two nodes.
        assert_eq!(find_superbubble(0, &succ, &pred, &lengths, 1, false), None);
        // 1 and 2 are not sources: each reaches 3 and 4, which do not meet before 5 only from 1.
        assert_eq!(find_superbubble(1, &succ, &pred, &lengths, 10, false), None);
        let keep = heaviest_path(0, 5, &[1, 2, 3, 4], &pred, |id| id, |id| [30, 30, 8, 10, 30, 30][id as usize], |_| 1);
        assert_eq!(keep, [1u32, 4].iter().copied().collect::<BTreeSet<u32>>());
    }

    #[test]
    fn superbubbles_leave_a_bubble_with_a_direct_edge_and_a_cycle_alone() {
        // 0 -> 1 -> 2 plus the direct edge 0 -> 2: metaMDBG skips it (the source reaches the exit directly).
        let nodes = vec![node(0, 30), node(1, 8), node(2, 30)];
        let edges = vec![edge(0, 1), edge(1, 2), edge(0, 2)];
        let mask = clean_graph(&nodes, &edges, &vec![true; 3], 0, 0, 10, 3);
        assert_eq!(mask, vec![true; 3]);
        // A cycle 1 -> 3 -> 1 inside the would-be bubble: no superbubble.
        let nodes = vec![node(0, 30), node(1, 30), node(2, 8), node(3, 30), node(4, 30)];
        let edges = vec![edge(0, 1), edge(0, 2), edge(1, 3), edge(3, 1), edge(1, 4), edge(2, 4)];
        let mask = clean_graph(&nodes, &edges, &vec![true; 5], 0, 0, 10, 3);
        assert_eq!(mask, vec![true; 5]);
    }

    #[test]
    fn superbubbles_keep_the_widest_path_not_the_most_fragmented() {
        // Path a = 1 -> 2 -> 3 (abundance 10 each), path b = 4 (abundance 20), source 0, exit 5: a has the larger sum of
        // abundances, b the larger weakest unitig.
        let mut pred = vec![BTreeSet::new(); 6];
        for (from, to) in [(0u32, 1usize), (1, 2), (2, 3), (0, 4), (3, 5), (4, 5)] {
            pred[to].insert(from);
        }
        let weight = |id: u32| [0, 10, 10, 10, 20, 0][id as usize];
        let widest = heaviest_path(0, 5, &[1, 2, 3, 4], &pred, |id| id, weight, |_| 1);
        assert_eq!(widest, [4u32].iter().copied().collect::<BTreeSet<u32>>());
        // Equal weakest unitig: the longer, better-covered branch wins.
        let mut pred = vec![BTreeSet::new(); 4];
        for (from, to) in [(0u32, 1usize), (0, 2), (1, 3), (2, 3)] {
            pred[to].insert(from);
        }
        let keep = heaviest_path(0, 3, &[1, 2], &pred, |id| id, |_| 10, |id| [0, 5, 2, 0][id as usize]);
        assert_eq!(keep, [1u32].iter().copied().collect::<BTreeSet<u32>>());
    }

    #[test]
    fn tips_cut_a_short_dead_end_whose_sibling_is_weaker() {
        // 5 -> 6 -> 7 -> 0 -> {1 -> 2 -> 3 (weak, 3 nodes), 4 (a 1-node dead end, stronger)}; tip limit 2: only 4 goes.
        let nodes = vec![node(0, 30), node(1, 8), node(2, 8), node(3, 30), node(4, 12), node(5, 30), node(6, 30), node(7, 30)];
        let edges = vec![edge(5, 6), edge(6, 7), edge(7, 0), edge(0, 1), edge(1, 2), edge(2, 3), edge(0, 4)];
        let mask = clean_graph(&nodes, &edges, &vec![true; 8], 2, 0, 0, 3);
        assert_eq!(mask, vec![true, true, true, true, false, true, true, true]);
        // Two short dead ends from one node: only the shorter goes (the other then compacts).
        let nodes = vec![node(0, 30), node(1, 30), node(2, 12), node(3, 12), node(4, 30), node(5, 30), node(6, 30)];
        let edges = vec![edge(4, 5), edge(5, 6), edge(6, 0), edge(0, 1), edge(0, 2), edge(2, 3)];
        let mask = clean_graph(&nodes, &edges, &vec![true; 7], 2, 0, 0, 3);
        assert_eq!(mask, vec![true, false, true, true, true, true, true]);
    }

    #[test]
    fn tips_unblock_a_superbubble_held_open_by_a_dead_end() {
        // 8 -> 9 -> 10 -> 0 -> {1, 2}; 1 -> 3; 2 -> {7 -> 3, 6 (dead end, stronger than its sibling 7)}; 3 -> 4 -> 11 ->
        // 12. Tip limit 2: the 4-node entry and exit unitigs are not tips.
        let ab = [30, 30, 8, 30, 30, 30, 12, 8, 30, 30, 30, 30, 30];
        let nodes: Vec<Node> = ab.iter().enumerate().map(|(i, &a)| node(i as u32, a)).collect();
        let edges = vec![edge(8, 9), edge(9, 10), edge(10, 0), edge(0, 1), edge(1, 3), edge(0, 2), edge(2, 7), edge(7, 3),
                         edge(2, 6), edge(3, 4), edge(4, 11), edge(11, 12)];
        let mut start = vec![true; 13];
        start[5] = false;
        let mask = clean_graph(&nodes, &edges, &start, 2, 16, 16, 3);
        let gone: Vec<usize> = (0..13).filter(|&i| start[i] && !mask[i]).collect();
        assert_eq!(gone, vec![2, 6, 7]);
    }

    #[test]
    fn tips_cut_a_dead_end_reached_from_two_branches_and_pop_the_bubble() {
        // 10..12 -> 0 -> {1 (strong) -> {2 -> 4, 3 -> 7}, 5 (weak) -> {6 -> 4, 8 -> 7}} and 4 -> 9 -> 13..15; 7 is a
        // dead end that is the only successor of 3 and 8. Tip limit 2 nodes, so the 4-node entry and 5-node exit
        // unitigs are not tips.
        let ab = [46, 46, 46, 8, 46, 18, 18, 8, 8, 46, 46, 46, 46, 46, 46, 46];
        let nodes: Vec<Node> = ab.iter().enumerate().map(|(i, &a)| node(i as u32, a)).collect();
        let edges = vec![edge(10, 11), edge(11, 12), edge(12, 0), edge(0, 1), edge(1, 2), edge(1, 3), edge(2, 4),
                         edge(3, 7), edge(0, 5), edge(5, 6), edge(5, 8), edge(6, 4), edge(8, 7), edge(4, 9), edge(9, 13),
                         edge(13, 14), edge(14, 15)];
        let all = vec![true; 16];
        let mm = clean_graph(&nodes, &edges, &all, 2, 16, 16, 3);
        let gone: Vec<usize> = (0..16).filter(|&i| !mm[i]).collect();
        assert_eq!(gone, vec![3, 5, 6, 7, 8]);                      // dead end, its two feeders, the weak arm
        // Two dead ends from one node: only the shorter goes, then the other compacts with the 5-node entry (metaMDBG
        // recompacts the predecessor first). Tip limit 4.
        let nodes: Vec<Node> = [30, 30, 12, 12, 30, 30, 30, 30].iter().enumerate().map(|(i, &a)| node(i as u32, a)).collect();
        let edges = vec![edge(4, 5), edge(5, 6), edge(6, 7), edge(7, 0), edge(0, 1), edge(0, 2), edge(2, 3)];
        let mask = clean_graph(&nodes, &edges, &vec![true; 8], 4, 0, 0, 3);
        assert_eq!(mask, vec![true, false, true, true, true, true, true, true]);
    }
}
