//! Running one step: `squash __STEP ARGS...` in a child process, its standard output and error in `LOGS/NAME.log` and
//! its command line in `LOGS/NAME.command.txt`. A step that fails stops the run with an error naming its log.

use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

pub struct Runner {
    exe: PathBuf,
    logs: PathBuf,
}

impl Runner {
    /// A runner of this executable's steps that writes its records to `logs`.
    pub fn new(logs: &Path) -> io::Result<Runner> {
        fs::create_dir_all(logs)?;
        Ok(Runner { exe: std::env::current_exe()?, logs: logs.to_path_buf() })
    }


    /// The child's command (`squash __STEP ARGS...`, its standard error to the log), with the command record written,
    /// and the log. The child gets `LC_ALL=C`; the diagnostics' `SQUASH_*` variables (docs/USAGE.md) reach it, and no
    /// variable changes a step's output.
    fn command(&self, name: &str, step: &str, args: &[String]) -> Result<(Command, File, PathBuf), String> {
        let log_path = self.logs.join(format!("{}.log", name));
        let mut cmd_file = File::create(self.logs.join(format!("{}.command.txt", name)))
            .map_err(|e| format!("{}: cannot write its command record: {}", name, e))?;
        let quoted: Vec<String> = std::iter::once(format!("squash __{}", step)).chain(args.iter().map(|a| shell_quote(a))).collect();
        writeln!(cmd_file, "{}", quoted.join(" ")).map_err(|e| format!("{}: {}", name, e))?;
        let log = File::create(&log_path).map_err(|e| format!("{}: cannot create {}: {}", name, log_path.display(), e))?;
        let err = log.try_clone().map_err(|e| format!("{}: {}", name, e))?;
        let mut cmd = Command::new(&self.exe);
        cmd.arg(format!("__{}", step)).args(args).stdin(Stdio::null()).stderr(err);
        cmd.env("LC_ALL", "C");
        Ok((cmd, log, log_path))
    }

    /// Runs `squash __STEP ARGS...` and waits for it (standard output and error in the log).
    pub fn step(&self, name: &str, step: &str, args: &[String]) -> Result<(), String> {
        let (mut cmd, log, log_path) = self.command(name, step, args)?;
        cmd.stdout(log);
        let status = cmd.status().map_err(|e| format!("{}: cannot start: {}", name, e))?;
        if status.success() {
            Ok(())
        } else {
            Err(format!("step {} failed ({}); see {}", name, status, log_path.display()))
        }
    }

    /// As `step`, for a step whose standard output reports progress (the correction rounds' `ROUND` / `SPARSE` lines):
    /// each line of it goes to the log, then to `on_line` while the step runs. An error from `on_line` stops the
    /// child and is returned.
    pub fn step_lines(&self, name: &str, step: &str, args: &[String],
                      mut on_line: impl FnMut(&str) -> Result<(), String>) -> Result<(), String> {
        let (mut cmd, mut log, log_path) = self.command(name, step, args)?;
        cmd.stdout(Stdio::piped());
        let mut child = cmd.spawn().map_err(|e| format!("{}: cannot start: {}", name, e))?;
        let stdout = child.stdout.take().ok_or_else(|| format!("{}: no standard output", name))?;
        let mut failed: Option<String> = None;
        for line in BufReader::new(stdout).lines() {
            let line = match line {
                Ok(line) => line,
                Err(e) => {
                    failed = Some(format!("{}: reading its output: {}", name, e));
                    break;
                }
            };
            if let Err(e) = writeln!(log, "{}", line) {
                failed = Some(format!("{}: {}", name, e));
                break;
            }
            if let Err(e) = on_line(&line) {
                failed = Some(format!("step {}: {}; see {}", name, e, log_path.display()));
                break;
            }
        }
        if failed.is_some() {
            let _ = child.kill();
        }
        let status = child.wait().map_err(|e| format!("{}: {}", name, e))?;
        match failed {
            Some(e) => Err(e),
            None if status.success() => Ok(()),
            None => Err(format!("step {} failed ({}); see {}", name, status, log_path.display())),
        }
    }
}

/// A shell word for the command record (the record is for reading; steps never go through a shell).
pub fn shell_quote(a: &str) -> String {
    if !a.is_empty() && a.chars().all(|c| c.is_ascii_alphanumeric() || "-_./:=,+%@".contains(c)) {
        a.to_string()
    } else {
        format!("'{}'", a.replace('\'', "'\\''"))
    }
}
