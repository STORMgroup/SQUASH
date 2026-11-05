"""Shared readers of the evaluation's Python scripts: SQUASH's signal layout (layout.tsv.gz), the graph's input reads
(work/graph/input/: each read's tokens and their q-symbol spans), the signal step's event spans (work/signal/spans.bin),
FASTA files. Evaluation only: nothing here is read by SQUASH."""
from __future__ import annotations

import gzip
import struct
from pathlib import Path

import numpy as np


def open_text(path: Path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def layout_rows(path: Path) -> list[dict]:
    """The rows of SQUASH's layout.tsv.gz (or a plain layout.tsv) as dicts of strings, in file order."""
    with open_text(path) as fh:
        lines = fh.read().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, x.split("\t"))) for x in lines[1:]]


def layout_reads(rows: list[dict]) -> set[str]:
    return {r["read_id"] for r in rows if r["read_id"] != "-"}


def graph_reads(folder: Path, want: set[str] | None = None) -> dict[str, np.ndarray]:
    """read id -> its tokens' q-symbol spans (array [n, 2]: first q-symbol, end) from a folder of token arrays
    (reads.tsv.gz: read_id, read_ordinal, token_offset, token_count; token_spans.u64le: two u64 per token), for the
    reads in `want` (every read when None); a read listed twice keeps its last row."""
    folder = Path(folder)
    spans = np.frombuffer((folder / "token_spans.u64le").read_bytes(), dtype="<u8").reshape(-1, 2)
    n_tokens = (folder / "tokens.u32le").stat().st_size // 4
    if n_tokens != len(spans):
        raise ValueError(f"{folder}: token and span arrays differ in length")
    out = {}
    with gzip.open(folder / "reads.tsv.gz", "rt") as fh:
        head = next(fh).rstrip("\n").split("\t")
        i, o, c = head.index("read_id"), head.index("token_offset"), head.index("token_count")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if want is None or f[i] in want:
                a, n = int(f[o]), int(f[c])
                out[f[i]] = spans[a:a + n]
    return out


def event_spans(path: Path, want: set[str] | None = None) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """read id -> (first sample, end sample) of each of its events, from the signal step's spans file (per read: u32
    name length, name, u32 samples, u32 events n, n u32 first samples, n u32 end samples), streamed: only the reads in
    `want` (every read when None) are kept; a read listed twice keeps its last record."""
    out = {}
    with open(path, "rb", buffering=1 << 24) as fh:
        while True:
            h = fh.read(4)
            if not h:
                break
            nl = struct.unpack("<I", h)[0]
            name = fh.read(nl).decode()
            n = struct.unpack("<II", fh.read(8))[1]
            if want is None or name in want:
                body = fh.read(8 * n)
                out[name] = (np.frombuffer(body, dtype="<u4", count=n), np.frombuffer(body, dtype="<u4", count=n, offset=4 * n))
            else:
                fh.seek(8 * n, 1)
    return out


def cut_spans(src: Path, want: set[str], out: Path) -> int:
    """Copies the records of the reads in `want` from a spans file to `out` (same format, source order); returns how
    many were written."""
    n_out = 0
    with open(src, "rb", buffering=1 << 24) as fh, open(out, "wb") as o:
        while True:
            h = fh.read(4)
            if not h:
                break
            nl = struct.unpack("<I", h)[0]
            name = fh.read(nl)
            hd = fh.read(8)
            n = struct.unpack("<II", hd)[1]
            if name.decode() in want:
                o.write(h + name + hd + fh.read(8 * n))
                n_out += 1
            else:
                fh.seek(8 * n, 1)
    return n_out


def fasta_lengths(path: Path) -> dict[str, int]:
    """record name (first word) -> sequence length, in file order."""
    out, name = {}, None
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = 0
            elif name is not None:
                out[name] += len(line.strip())
    return out


def write_fasta_record(fh, header: str, seq: str, width: int = 80) -> None:
    fh.write(f">{header}\n")
    fh.writelines(seq[i:i + width] + "\n" for i in range(0, len(seq), width))
