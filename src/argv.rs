//! The argument list of a step. Each step of the pipeline keeps its own option parser; `squash` runs it as
//! `squash __STEP ARGS...` in a child process, and the step reads `[STEP, ARGS...]` here instead of
//! `std::env::args()`, so its parser sees the step name as the program name followed by the step's own arguments.

use std::sync::OnceLock;

static ARGS: OnceLock<Vec<String>> = OnceLock::new();

/// Set once, by `main`, before a step runs.
pub fn set(args: Vec<String>) {
    ARGS.set(args).expect("step arguments set twice");
}

/// The step's argument list (program name first), or the process's own when no step is running.
pub fn args() -> std::vec::IntoIter<String> {
    match ARGS.get() {
        Some(a) => a.clone().into_iter(),
        None => std::env::args().collect::<Vec<_>>().into_iter(),
    }
}
