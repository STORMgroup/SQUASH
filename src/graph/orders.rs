//! The order loop of `assemble`: each order processed (cleaning, contigs) and written while the next one is
//! computed.

use flate2::Compression;
use structopt::StructOpt;
use std::cmp::max;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use super::build::{next_state_with, seed_state};
use super::caps::{ResourceCaps, check_caps, enforce_rss_cap, stop_on_cap, stop_or_error};
use super::clean::cleaned_masks;
use super::clock::PhaseClock;
use super::contigs::{
    compact_unitigs, mask_from_unitigs, select_contigs_bounded, write_unitig_graph, write_unitig_table,
};
use super::files::{GFA, NODE_TOKENS, UNITIGS, create_plain, finish_plain, invalid_data, invalid_input, require_threads,
                   suffix_path};
use super::state::{NodeKey, State, state_is_seed};
use super::tokens::{TokenRead, load_prepared};

/// One process for the whole multi-k graph: the seed order `k_start` built from the prepared tokens, every order up to
/// `k_max` (each computed from the one before, while that one is written).
#[derive(Debug, StructOpt)]
pub(super) struct AssembleOpt {
    /// The input reads' folder: reads.tsv.gz, tokens.u32le, token_spans.u64le.
    #[structopt(long, parse(from_os_str))]
    reads_dir: PathBuf,
    /// The seed order (tokens per node).
    #[structopt(long)]
    k_start: u16,
    /// The last order.
    #[structopt(long)]
    k_max: u16,
    /// Orders are written to DIR/kNN.
    #[structopt(long, parse(from_os_str))]
    graph_dir: PathBuf,
    #[structopt(long, default_value = "1")]
    pub(super) threads: usize,
    /// Resource cap: read-path rows of an order (a stop when exceeded, no other effect).
    #[structopt(long, default_value = "100000000")]
    max_path_rows: u64,
    /// Cut tips up to this many tokens (unitig length) at every order (metaMDBG cuts tips up to 50 kb). Absent = 8 x K
    /// nodes.
    #[structopt(long)]
    tip_max_tokens: Option<u32>,
    /// The two progressive abundance cutoffs "A,B" with 0 <= A <= B, instead of metaMDBG's 2 and 3: the first filter
    /// drops the unitigs of median abundance <= A and the second those <= B (abundance2 <= 2A, 2B). The contigs'
    /// minimum abundances stay.
    #[structopt(long)]
    progressive_cutoffs: Option<String>,
    /// The longest superbubble path popped, in tokens of the path (nodes + K - 1; metaMDBG pops superbubbles up
    /// to 50 kb at every order). Absent = 8 x K nodes.
    #[structopt(long)]
    superbubble_max_tokens: Option<u32>,
}

#[cfg(test)]
pub(super) fn process_state(state: State) -> State {
    process_state_bounded(state, None).expect("uncapped state processing cannot fail")
}

/// An order's processing: the three cleaned masks (`cleaned_masks`), the contigs passed to the next order
/// (`select_contigs_bounded`), the last mask as the active nodes, and each node's `abundance_out` (its contig's median
/// abundance, at least 2 for a contig longer than k; otherwise `abundance_in`).
fn process_state_bounded(
    mut state: State,
    caps: Option<&ResourceCaps>,
) -> io::Result<State> {
    enforce_rss_cap(caps)?;
    let mut clock = PhaseClock::new("process_state");
    let (cleaned, first_cleaned, second_cleaned) = cleaned_masks(&state);
    clock.mark("cleaned masks");
    let contigs = select_contigs_bounded(&state.nodes, &state.edges, &cleaned, &first_cleaned, &second_cleaned, caps)?;
    state.active = second_cleaned;
    clock.mark("contigs");
    // Abundance smoothing is a property of a selected contig, not of merely
    // surviving a cleaning. Unclaimed nodes retain a_in exactly.
    for node in &mut state.nodes {
        node.abundance_out = node.abundance_in;
    }
    for contig in &contigs {
        let mut value = (contig.abundance2.saturating_add(1)) / 2;
        if contig.nodes.len() > state.k as usize {
            value = max(value, 2);
        }
        for &id in &contig.nodes {
            state.nodes[id as usize].abundance_out = value;
        }
    }
    state.contigs = contigs;
    Ok(state)
}

fn write_atomic_token_file(
    path: &Path,
    state: &State,
    previous_prefix: Option<&Path>,
    previous_header: Option<(u16, usize)>,
) -> io::Result<()> {
    let mut writer = create_plain(path)?;
    // A state is a seed state exactly when its nodes carry token-vector keys (`state_is_seed`), so any seed order
    // works.
    match state_is_seed(state) {
        true => {
            let width = state.k as usize;
            let mut record = vec![0u8; width * 4];
            for node in &state.nodes {
                if let NodeKey::Seed(v) = &node.key {
                    if v.len() != width {
                        return Err(invalid_data("seed atomic vector length mismatch"));
                    }
                    for (slot, &x) in record.chunks_exact_mut(4).zip(v) {
                        slot.copy_from_slice(&x.to_le_bytes());
                    }
                    writer.write_all(&record)?;
                } else {
                    return Err(invalid_data("seed state contains parent-pair node"));
                }
            }
        }
        false => {
            let k = state.k;
            let (prefix, (prev_k, prev_count)) = match (previous_prefix, previous_header) {
                (Some(prefix), Some(header)) => (prefix, header),
                _ => return Err(invalid_input("higher-k token publication requires the previous order")),
            };
            if prev_k + 1 != k {
                return Err(invalid_data("previous token order mismatch"));
            }
            // The previous order's token array is read whole (a few hundred MB
            // at most) instead of with two positioned reads per node.
            let previous = fs::read(suffix_path(prefix, NODE_TOKENS)?)?;
            let width = prev_k as usize * 4;
            if Some(previous.len()) != prev_count.checked_mul(width) {
                return Err(invalid_data("atomic token file length mismatch"));
            }
            let vector = |id: u32| -> io::Result<&[u8]> {
                if id as usize >= prev_count {
                    return Err(invalid_data("atomic token ID out of range"));
                }
                Ok(&previous[id as usize * width..(id as usize + 1) * width])
            };
            for node in &state.nodes {
                if let NodeKey::Parents(a, b) = node.key {
                    let left = vector(a)?;
                    let right = vector(b)?;
                    if left[4..] != right[..width - 4] {
                        return Err(invalid_data("parent token vectors violate exact overlap"));
                    }
                    writer.write_all(left)?;
                    writer.write_all(&right[width - 4..])?;
                } else {
                    return Err(invalid_data("higher-k state contains a seed node key"));
                }
            }
        }
    }
    finish_plain(writer)
}

/// Writes order `state.k` under `prefix`: its node tokens (`.node_tokens.u32le`, each node's k tokens, built
/// from the previous order's: `previous` is that order's prefix, order and node count) and, at K-max (`final_order`),
/// the final graph: its unitig table (`.unitigs.tsv.gz`) and GFA (`.unitigs.gfa`).
fn write_order(
    prefix: &Path,
    previous: Option<&(PathBuf, (u16, usize))>,
    state: &State,
    caps: &ResourceCaps,
    final_order: bool,
) -> io::Result<()> {
    let mut clock = PhaseClock::new("write_order");
    let atomic_token_bytes = (state.nodes.len() as u64)
        .checked_mul(state.k as u64)
        .and_then(|value| value.checked_mul(4))
        .ok_or_else(|| invalid_data("iteration atomic-token byte estimate overflow"))?;
    if atomic_token_bytes > caps.max_scratch_bytes {
        return stop_on_cap(format!(
            "minimum_scratch_bytes={} exceeds {}",
            atomic_token_bytes, caps.max_scratch_bytes
        ));
    }
    if let Some(parent) = prefix.parent() {
        fs::create_dir_all(parent)?;
    }
    let tokens = suffix_path(prefix, NODE_TOKENS)?;
    let (previous_prefix, previous_header) = (previous.map(|p| p.0.as_path()), previous.map(|p| p.1));
    write_atomic_token_file(&tokens, state, previous_prefix, previous_header)?;
    clock.mark("tokens");
    if final_order {
        let unitigs = compact_unitigs(&state.nodes, &state.edges, &state.active);
        write_unitig_table(&suffix_path(prefix, UNITIGS)?, state.k, "P3", &unitigs, Compression::fast())?;
        let mask = mask_from_unitigs(state.nodes.len(), &unitigs);
        write_unitig_graph(&suffix_path(prefix, GFA)?, state, &mask, &unitigs)?;
        clock.mark("final_graph");
    }
    Ok(())
}

/// A computed order's state: the order's graph processed (`process_state_bounded`: cleaning, masks, contigs).
fn compute_order(raw: State, caps: &ResourceCaps) -> io::Result<State> {
    process_state_bounded(raw, Some(caps))
}

/// The order loop of `assemble`. Order k is written once order k+1 is computed (while k+2 is computed from it).
fn order_loop(
    opt: &AssembleOpt,
    caps: &ResourceCaps,
    reads: &[TokenRead],
    first: State,
    clock: &mut PhaseClock,
) -> io::Result<()> {
    let order_prefix = |k: u16| opt.graph_dir.join(format!("k{:02}", k));
    let mut current = first;
    let mut k = opt.k_start;
    let mut pending: Option<(u16, State)> = None;
    // The last order written: its prefix, order and node count.
    let mut previous: Option<(PathBuf, (u16, usize))> = None;
    loop {
        if let Err(reason) = check_caps(&current, caps, reads) {
            return stop_on_cap(reason);
        }
        let prefix = order_prefix(k);
        let next_k = k + 1;
        // `SQUASH_GRAPH_SERIAL=1` (profiling only): the pending order is written first and the next order computed after it,
        // so each phase's process CPU belongs to that phase alone. The same files.
        let serial = std::env::var_os("SQUASH_GRAPH_SERIAL").map_or(false, |v| v == "1");
        let write_pending = || -> io::Result<()> {
            let Some((pk, order)) = &pending else { return Ok(()) };
            write_order(&order_prefix(*pk), previous.as_ref(), order, caps, false)
        };
        let compute_next = || {
            if next_k > opt.k_max {
                return None;
            }
            Some(next_state_with(&current, &prefix, Some(caps)).and_then(|raw| compute_order(raw, caps)))
        };
        let (written, next) = if serial {
            let written = write_pending();
            (written, compute_next())
        } else {
            rayon::join(write_pending, compute_next)
        };
        written?;
        if let Some((pk, order)) = pending.take() {
            clock.mark(&format!("k{:02} write (+k{:02} compute)", pk, next_k));
            previous = Some((order_prefix(pk), (pk, order.nodes.len())));
        }
        // The written order's memory back to the system (glibc keeps freed memory per thread).
        unsafe {
            libc::malloc_trim(0);
        }
        pending = Some((k, current));
        match next {
            None => break,
            Some(result) => {
                current = result.map_err(stop_or_error)?;
                k = next_k;
            }
        }
    }
    // K-max: the final graph's files.
    let (pk, order) = pending.take().expect("the loop leaves the last order pending");
    write_order(&order_prefix(pk), previous.as_ref(), &order, caps, true)?;
    clock.mark(&format!("k{:02} write", pk));
    drop(order);
    unsafe {
        libc::malloc_trim(0);
    }
    Ok(())
}

pub(super) fn assemble(opt: AssembleOpt) -> io::Result<()> {
    let mut clock = PhaseClock::new("assemble");
    require_threads(opt.threads)?;
    if opt.max_path_rows == 0 {
        return Err(invalid_input("--max-path-rows must be positive"));
    }
    let caps = ResourceCaps::with_path_rows(opt.max_path_rows);
    if opt.k_start < 2 || opt.k_max < opt.k_start {
        return Err(invalid_input("assemble needs 2 <= --k-start <= --k-max"));
    }
    if opt.tip_max_tokens == Some(0) {
        return Err(invalid_input("--tip-max-tokens needs a positive value"));
    }
    let progressive_cutoffs = match &opt.progressive_cutoffs {
        None => None,
        Some(text) => {
            let parts: Vec<&str> = text.split(',').collect();
            let parsed: Option<Vec<u32>> = parts.iter().map(|x| x.trim().parse().ok()).collect();
            match parsed.as_deref() {
                Some([a, b]) if a <= b => Some((*a, *b)),
                _ => return Err(invalid_input("--progressive-cutoffs needs A,B with A <= B")),
            }
        }
    };
    if opt.superbubble_max_tokens == Some(0) {
        return Err(invalid_input("--superbubble-max-tokens needs a positive value"));
    }
    let reads = load_prepared(&opt.reads_dir)?;
    clock.mark("reads");
    let windows: u64 = reads
        .iter()
        .map(|r| r.tokens.len().saturating_sub(opt.k_start as usize - 1) as u64)
        .sum();
    if windows > caps.max_path_rows {
        // The seed graph stops on the first path row over the cap.
        return stop_on_cap(format!("path_rows={} exceeds {}", caps.max_path_rows + 1, caps.max_path_rows));
    }
    // The seed order: every k_start-token window of the reads seen in at least two distinct reads.
    let mut seed = seed_state(&reads, opt.k_start)?;
    seed.active = vec![true; seed.nodes.len()];
    seed.tip_max_tokens = opt.tip_max_tokens;
    seed.progressive_cutoffs = progressive_cutoffs;
    seed.superbubble_max_tokens = opt.superbubble_max_tokens;
    let first = compute_order(seed, &caps).map_err(stop_or_error)?;
    clock.mark(&format!("k{:02} compute", opt.k_start));
    order_loop(&opt, &caps, &reads, first, &mut clock)
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::test_support::*;

    #[test]
    fn even_median_uses_rational_half_and_ceil_smoothing() {
        let raw = state(vec![node(0, 2), node(1, 3)],
            vec![edge(0, 1)],
        );
        let out = process_state(raw);
        assert_eq!(out.contigs[0].abundance2, 5);
        assert_eq!(out.nodes[0].abundance_out, 3);
        assert_eq!(out.nodes[1].abundance_out, 3);
    }

    #[test]
    fn a_long_contig_lifts_abundance_one_to_two() {
        let nodes = (0..21).map(|i| node(i, 1)).collect::<Vec<_>>();
        let edges = (0..20).map(|i| edge(i, i + 1)).collect();
        let out = process_state(state(nodes, edges));
        assert_eq!(out.contigs[0].nodes.len(), 21);
        assert!(out.nodes.iter().all(|n| n.abundance_out == 2));
    }

    #[test]
    fn contigs_come_from_the_cleaned_graph_so_popped_branches_are_not_contigs() {
        // The crossed bubble behind a 3-node entry and a 3-node exit, with a tip limit of one node: only the
        // superbubble rule changes it.
        let (mut nodes, mut edges) = crossed_bubble();
        for i in 6..12u32 {
            nodes.push(node(i, 30));
        }
        edges.extend([edge(6, 7), edge(7, 8), edge(8, 0), edge(5, 9), edge(9, 10), edge(10, 11)]);
        let mut st = state(nodes, edges);
        st.k = 2;
        st.tip_max_tokens = Some(1);
        let out = process_state(st);
        let has = |id: u32| out.contigs.iter().any(|c| c.nodes.contains(&id));
        assert!(!out.active[2] && !out.active[3]);
        assert!(!has(2) && !has(3), "popped branches are not contigs of the cleaned graph");
    }
}
