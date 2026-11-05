//! Resource caps: reaching one stops the step with a `RESOURCE_CAP` error; no other effect.

use std::fs;
use std::io;

use super::files::invalid_data;
use super::state::State;
use super::tokens::TokenRead;

/// The graph's resource caps: reaching one stops the step (`stop_on_cap`); no other effect.
#[derive(Clone, Debug)]
pub(super) struct ResourceCaps {
    pub(super) max_rss_bytes: u64,
    pub(super) max_nodes: u64,
    pub(super) max_links: u64,
    pub(super) max_path_rows: u64,
    pub(super) max_contigs: u64,
    pub(super) max_contig_token_ratio: u64,
    pub(super) max_scratch_bytes: u64,
}

impl ResourceCaps {
    /// The caps with `max_path_rows` read-path rows (64 GiB resident, 2e7 nodes, 4e7 links, 1e6 contigs passed on, 4 x
    /// the reads' tokens in those contigs, 100 GiB of order files).
    pub(super) fn with_path_rows(max_path_rows: u64) -> Self {
        Self {
            max_rss_bytes: 68719476736,
            max_nodes: 20000000,
            max_links: 40000000,
            max_path_rows,
            max_contigs: 1000000,
            max_contig_token_ratio: 4,
            max_scratch_bytes: 107374182400,
        }
    }
}

const CAP_ERROR_PREFIX: &str = "CAP: ";

pub(super) fn cap_error(message: impl Into<String>) -> io::Error {
    io::Error::new(
        io::ErrorKind::Other,
        format!("{}{}", CAP_ERROR_PREFIX, message.into()),
    )
}

pub(super) fn cap_reason(error: &io::Error) -> Option<String> {
    error
        .to_string()
        .strip_prefix(CAP_ERROR_PREFIX)
        .map(str::to_owned)
}

pub(super) fn enforce_rss_cap(caps: Option<&ResourceCaps>) -> io::Result<()> {
    if let Some(caps) = caps {
        if let Ok(rss) = current_rss_bytes() {
            if rss > caps.max_rss_bytes {
                return Err(cap_error(format!(
                    "rss_bytes={} exceeds {}",
                    rss, caps.max_rss_bytes
                )));
            }
        }
    }
    Ok(())
}

fn current_rss_bytes() -> io::Result<u64> {
    let text = fs::read_to_string("/proc/self/status")?;
    for line in text.lines() {
        if let Some(v) = line.strip_prefix("VmRSS:") {
            let kb = v
                .split_whitespace()
                .next()
                .ok_or_else(|| invalid_data("malformed VmRSS"))?
                .parse::<u64>()
                .map_err(|_| invalid_data("malformed VmRSS value"))?;
            return kb
                .checked_mul(1024)
                .ok_or_else(|| invalid_data("VmRSS overflow"));
        }
    }
    Err(invalid_data("/proc/self/status lacks VmRSS"))
}

pub(super) fn check_caps(state: &State, caps: &ResourceCaps, reads: &[TokenRead]) -> Result<(), String> {
    let path_rows: usize = state.read_paths.iter().map(|p| p.steps.len()).sum();
    let contig_tokens: u64 = state
        .contigs
        .iter()
        .map(|c| (c.nodes.len() + state.k as usize - 1) as u64)
        .sum();
    let read_tokens: u64 = reads.iter().map(|r| r.tokens.len() as u64).sum();
    let contig_token_cap = read_tokens.saturating_mul(caps.max_contig_token_ratio);
    let checks = [
        ("nodes", state.nodes.len() as u64, caps.max_nodes),
        ("links", state.edges.len() as u64, caps.max_links),
        ("path_rows", path_rows as u64, caps.max_path_rows),
        (
            "contigs",
            state.contigs.len() as u64,
            caps.max_contigs,
        ),
        ("contig_tokens", contig_tokens, contig_token_cap),
    ];
    for (name, value, limit) in checks {
        if value > limit {
            return Err(format!("{}={} exceeds {}", name, value, limit));
        }
    }
    if let Ok(rss) = current_rss_bytes() {
        if rss > caps.max_rss_bytes {
            return Err(format!("rss_bytes={} exceeds {}", rss, caps.max_rss_bytes));
        }
    }
    Ok(())
}

/// A resource cap reached: the step ends with the message `RESOURCE_CAP: REASON` (exit status 2, `main.rs`).
pub(super) fn stop_on_cap<T>(reason: String) -> io::Result<T> {
    Err(io::Error::new(io::ErrorKind::Other, format!("RESOURCE_CAP: {}", reason)))
}

/// A resource cap reached while computing an order stops the step (`stop_on_cap`).
pub(super) fn stop_or_error(error: io::Error) -> io::Error {
    match cap_reason(&error) {
        Some(reason) => stop_on_cap::<()>(reason).unwrap_err(),
        None => error,
    }
}
