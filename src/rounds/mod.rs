//! The correction rounds: each read's tokens are replaced by a token-level partial-order consensus of the read (the
//! target) and its partners (the donors), as nanoMDBG corrects reads in minimizer space.
//!
//! One replay (`replay_read`): every donor is aligned to the target over their overlap (the donor row's q-symbol
//! bounds as tokens, widened by `pad` tokens; `align_overlap`: +2 match, -2 mismatch, -3 gap, free end gaps within
//! 2 x `pad` tokens of the window's ends). Between its first and last paired columns, an
//! alignment votes at every target layer for an allele (the donor tokens inserted before the layer, then the donor's
//! token there or a deletion) and for the transition from the previous layer's allele; the target votes for its own
//! tokens. The corrected read is the heaviest path through the layers (allele plus transition support), ties broken
//! by the emitted tokens compared backwards (`content_path`); its tokens keep the target's coordinates, inserted
//! tokens spread before their layer (`spanned_tokens`).
//!
//! `rounds rounds` runs every round in one process (`run_rounds`) and writes each round's token table.
//!
//! Reads are replayed independently and written in catalog order: the output does not depend on the thread count.

use crate::io::{Result, open_text};
#[cfg(test)]
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};


/// The partner search's token (value and q span), shared so that the rounds search partners on the corrected tokens in
/// memory.
pub use crate::partners::Token;

mod align;
mod replay;
mod skip;
mod input;
mod tables;
mod step;
#[cfg(test)]
mod testing;

use self::align::*;
use self::replay::*;
use self::skip::*;
use self::input::*;
use self::tables::*;
#[cfg(test)]
use self::testing::*;

pub use self::step::run;
