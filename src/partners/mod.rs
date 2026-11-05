//! Partner search of the correction (step `partners`; `rounds` calls it every round). For every read (the target), the
//! reads that share token values with it (donors) are found through an index of the values, their shared tokens
//! (anchors) are chained, and a chain is kept when the overlap it implies and the part of that overlap it covers are
//! long enough; the kept donors are ranked by chain score and taken best first under a coverage cap. These are
//! nanoMDBG's correction rules (repeat values are not seeds, overlap, cover, coverage cap) on signal tokens instead of
//! base minimizers. Optional filters: the raw-evidence guard of the rounds (the uncorrected tokens of the pair
//! must still agree) and the q-symbol alignment check (`signal`: the pair's q-symbols aligned over the implied
//! overlap).
//!
//! Output: the selected-donor table (`COLUMNS`; header, then one row per kept donor, targets in catalog order),
//! written uncompressed; the pipeline gzips it. Each target is searched on its own and its chains are ranked by fixed
//! rules (score, then donor id), so the table does not depend on the thread count (the check's verdicts are per
//! directed pair).

use crate::io::{Result, open_text};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::convert::TryInto;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};


mod diagnostics;
mod tokens;
mod guard;
mod params;
mod index;
mod chain;
pub mod heavy;
mod search;
mod donors;
mod free;
mod rows;
mod input;
mod targets;
mod step;
pub mod signal;
#[cfg(test)]
mod testing;

use self::diagnostics::*;
use self::guard::*;
use self::params::*;
use self::index::*;
use self::chain::*;
use self::search::*;
use self::donors::*;
use self::free::*;
use self::rows::*;
#[cfg(test)]
use self::targets::*;
#[cfg(test)]
use self::testing::*;

// the items the correction rounds use
pub use self::diagnostics::mem;
pub use self::guard::{load_guard_view, guard_tokens, spill_guard_view, GuardRead};
pub use self::index::kept_values;
pub use self::input::{load_excluded, load_prepared};
pub use self::params::{take_guard, take_params};
pub use self::rows::{write_packed_rows, PackedRows, COLUMNS};
pub use self::step::{run, set_mmap_threshold};
pub use self::targets::partner_rows_packed;
pub use self::tokens::Token;
