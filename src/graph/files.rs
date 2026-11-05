//! Helpers for the engine's files: errors, paths, plain and gzip writers (gzip deflated on a worker thread), binary
//! arrays.

use flate2::write::GzEncoder;
use flate2::{Compression, GzBuilder};
use std::convert::TryFrom;
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

/// The files of order k in a graph folder, `kNN` + suffix: the node tokens (k tokens per node, in node order), the
/// unitig table and the GFA.
pub(crate) const NODE_TOKENS: &str = ".node_tokens.u32le";
pub(crate) const UNITIGS: &str = ".unitigs.tsv.gz";
pub(crate) const GFA: &str = ".unitigs.gfa";

/// The file `suffix` of order k in graph folder `dir`.
pub(crate) fn order_file(dir: &Path, k: usize, suffix: &str) -> PathBuf {
    dir.join(format!("k{:02}{}", k, suffix))
}

pub(super) fn invalid_data(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

pub(super) fn invalid_input(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message.into())
}

pub(super) fn suffix_path(prefix: &Path, suffix: &str) -> io::Result<PathBuf> {
    let text = prefix
        .to_str()
        .ok_or_else(|| invalid_input("output prefix is not UTF-8"))?;
    Ok(PathBuf::from(format!("{}{}", text, suffix)))
}

pub(super) fn require_threads(threads: usize) -> io::Result<()> {
    if threads == 0 {
        Err(invalid_input("--threads must be at least one"))
    } else {
        Ok(())
    }
}

pub(super) fn create_plain(path: &Path) -> io::Result<BufWriter<File>> {
    let file = OpenOptions::new().write(true).create_new(true).open(path)?;
    Ok(BufWriter::with_capacity(1 << 20, file))
}

/// A gzip sink with a buffer in front of the deflate stage, so the encoder receives megabyte blocks instead of every
/// `writeln!` fragment. Deflate output depends on the byte sequence and the flush points, not on how the writes are
/// chunked, so the compressed bytes are those of an unbuffered encoder.
type GzWriter = BufWriter<GzEncoder<BufWriter<File>>>;

#[cfg(test)]
fn create_gzip(path: &Path) -> io::Result<GzWriter> {
    create_gzip_level(path, Compression::default())
}

fn create_gzip_level(path: &Path, level: Compression) -> io::Result<GzWriter> {
    let file = OpenOptions::new().write(true).create_new(true).open(path)?;
    Ok(BufWriter::with_capacity(
        1 << 20,
        GzBuilder::new()
            .mtime(0)
            .write(BufWriter::with_capacity(1 << 20, file), level),
    ))
}

pub(super) fn finish_plain(mut writer: BufWriter<File>) -> io::Result<()> {
    writer.flush()?;
    writer.get_ref().sync_all()
}

fn finish_gzip(writer: GzWriter) -> io::Result<()> {
    let encoder = writer.into_inner().map_err(|error| error.into_error())?;
    let mut inner = encoder.finish()?;
    inner.flush()?;
    inner.get_ref().sync_all()
}

/// A gzip artifact whose deflate stage runs on its own thread.
///
/// The producer only formats rows and copies them into a block buffer; a worker thread owns the encoder, so
/// independent artifacts compress concurrently and overlap the producer's own work. The output bytes are those of a
/// single-threaded writer: the same encoder, at the same level, consumes the same byte sequence with no intermediate
/// flushes.
pub(super) struct BackgroundGz {
    sender: Option<std::sync::mpsc::SyncSender<Vec<u8>>>,
    worker: Option<std::thread::JoinHandle<io::Result<()>>>,
    pending: Vec<u8>,
}

impl BackgroundGz {
    const BLOCK: usize = 1 << 20;
    /// Blocks a sink may queue ahead of its deflate thread (16 MiB).
    const QUEUE_BLOCKS: usize = 16;

    pub(super) fn create(path: &Path) -> io::Result<Self> {
        Self::create_level(path, Compression::default())
    }

    pub(super) fn create_level(path: &Path, level: Compression) -> io::Result<Self> {
        let mut writer = create_gzip_level(path, level)?;
        let (sender, receiver) = std::sync::mpsc::sync_channel::<Vec<u8>>(Self::QUEUE_BLOCKS);
        let worker = std::thread::Builder::new()
            .name("graph-gzip".to_owned())
            .spawn(move || -> io::Result<()> {
                for block in receiver {
                    writer.write_all(&block)?;
                }
                finish_gzip(writer)
            })?;
        Ok(Self {
            sender: Some(sender),
            worker: Some(worker),
            pending: Vec::with_capacity(Self::BLOCK),
        })
    }

    fn send_pending(&mut self) -> io::Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }
        let block = std::mem::replace(&mut self.pending, Vec::with_capacity(Self::BLOCK));
        match self.sender.as_ref() {
            Some(sender) => sender.send(block).map_err(|_| {
                io::Error::new(io::ErrorKind::Other, "gzip worker stopped accepting output")
            }),
            None => Err(invalid_data("gzip sink already finished")),
        }
    }

    /// Flushes the tail block and propagates any worker error.  Must be called
    /// for every sink; dropping without finishing abandons buffered output.
    pub(super) fn finish(mut self) -> io::Result<()> {
        self.send_pending()?;
        self.sender = None;
        match self.worker.take() {
            Some(handle) => handle
                .join()
                .map_err(|_| io::Error::new(io::ErrorKind::Other, "gzip worker panicked"))?,
            None => Ok(()),
        }
    }
}

impl Write for BackgroundGz {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.pending.extend_from_slice(buf);
        if self.pending.len() >= Self::BLOCK {
            self.send_pending()?;
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl Drop for BackgroundGz {
    fn drop(&mut self) {
        // A sink is normally consumed by `finish`.  On an error path the
        // artifact is discarded anyway, so just let the worker wind down.
        self.sender = None;
        if let Some(handle) = self.worker.take() {
            let _ = handle.join();
        }
    }
}

pub(super) fn write_u32<W: Write>(w: &mut W, v: u32) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

pub(super) fn write_u64<W: Write>(w: &mut W, v: u64) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

pub(super) fn usize_u32(v: usize, label: &str) -> io::Result<u32> {
    u32::try_from(v).map_err(|_| invalid_data(format!("{} exceeds u32", label)))
}

pub(super) fn u64_usize(v: u64, label: &str) -> io::Result<usize> {
    usize::try_from(v).map_err(|_| invalid_data(format!("{} exceeds usize", label)))
}

/// A TSV field with minimal quoting: quoted, with doubled quotes, only when it holds a tab, a quote or a line break.
pub(super) fn csv_field(x: &str) -> String {
    if x.contains(['\t', '"', '\n', '\r']) {
        format!("\"{}\"", x.replace('"', "\"\""))
    } else {
        x.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn deterministic_gzip_has_zero_mtime() {
        let dir = std::env::temp_dir().join(format!("squash_graph_gzip_test_{}", std::process::id()));
        if dir.exists() {
            fs::remove_dir_all(&dir).unwrap();
        }
        fs::create_dir_all(&dir).unwrap();
        let gz = dir.join("x.gz");
        let mut w = create_gzip(&gz).unwrap();
        writeln!(w, "x").unwrap();
        finish_gzip(w).unwrap();
        let bytes = fs::read(&gz).unwrap();
        assert_eq!(&bytes[4..8], &[0, 0, 0, 0]);
        fs::remove_dir_all(dir).unwrap();
    }
}
