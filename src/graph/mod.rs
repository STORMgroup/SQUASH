//! Multi-k minimizer-space de Bruijn graph assembly over signal tokens.
//!
//! Input: the prepared token stream (each read's corrected q-symbols sketched into minimizer tokens). The graph is
//! built at the seed order k_start and the contigs of each order are threaded into the next one as paths, up to k_max;
//! at every order the graph is cleaned (tips, bubbles, superbubbles, progressive abundance filters) and its contigs are
//! selected. Outputs: each order's node tokens (read by the next order) and the K-max order's unitig table and GFA (the
//! final graph).
//! The multi-k contig paths, cleaning and contig rules follow metaMDBG (github.com/GaetanBenoitDev/metaMDBG, v1.2); the
//! minimizer-space de Bruijn graph follows rust-mdbg (github.com/ekimb/rust-mdbg).
//!
//! The module never reads alignments or a reference. Seed-order nodes keep their exact token vectors; a node of a
//! higher order k is identified by the compact `(prefix_(k-1), suffix_(k-1))` pair of previous-order nodes, which
//! avoids storing an O(k) vector for every node while keeping the identity collision free.
//!
//! Files: `tokens` (the prepared tokens: the `prepare` and `adapt-tokens` steps), `state` (an order's graph), `build`
//! (the seed order and each next order), `clean` (tips, bubbles, superbubbles, abundance filters), `contigs` (unitigs,
//! the contigs passed to the next order, the unitig table and GFA), `orders` (the order loop, the `assemble` step),
//! `caps` (resource caps), `clock` (the profile), `files` (file helpers).

mod tokens;
mod state;
mod build;
mod clean;
mod contigs;
mod orders;
mod caps;
mod clock;
mod files;
#[cfg(test)]
mod test_support;

use structopt::StructOpt;
use std::io;

use files::{invalid_input, require_threads};
pub(crate) use files::{GFA, NODE_TOKENS, UNITIGS, order_file};
use orders::{AssembleOpt, assemble};
use tokens::{AdaptTokensOpt, PrepareOpt, adapt_tokens, prepare};

#[derive(Debug, StructOpt)]
#[structopt(name = "graph", about = "the multi-k minimizer-space de Bruijn graph")]
enum Command {
    Prepare(PrepareOpt),
    #[structopt(name = "adapt-tokens")]
    AdaptTokens(AdaptTokensOpt),
    #[structopt(name = "assemble")]
    Assemble(AssembleOpt),
}

/// The thread budget the selected subcommand was given.
fn command_threads(command: &Command) -> usize {
    match command {
        Command::Prepare(opt) => opt.threads,
        Command::AdaptTokens(_) => 1,
        Command::Assemble(opt) => opt.threads,
    }
}

pub fn run() -> io::Result<()> {
    let command = Command::from_iter(crate::argv::args());
    let threads = command_threads(&command);
    require_threads(threads)?;
    // One pool for the whole process, installed around the entire dispatch, so every rayon call site (bare
    // `par_sort_unstable` calls included, which would otherwise use rayon's global pool sized to the machine) honours
    // `--threads`.
    //
    // Thread count does not affect output. Every parallel site is either an order-preserving `map`/`collect` or an
    // unstable sort over a totally ordered key, and `parallel_recurrence_is_identical_at_1_4_16_threads` pins
    // that. Workers get 64 MiB stacks (reserved, touched only as used): a worker waiting inside one contig's nested
    // parallel work runs other queued work on top of its stack, and with many small tasks queued the default 2 MiB can
    // overflow.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(64 << 20)
        .build()
        .map_err(|e| invalid_input(format!("cannot create the graph's thread pool: {}", e)))?;
    pool.install(|| match command {
        Command::Prepare(opt) => prepare(opt),
        Command::AdaptTokens(opt) => adapt_tokens(opt),
        Command::Assemble(opt) => assemble(opt),
    })
}
