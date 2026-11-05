#!/usr/bin/env python3
"""Tests of select_reads.py on synthetic POD5s written here (no data file; pod5 0.3.39).
Run: python3 test_select.py  (or python3 -m unittest test_select)."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

import gzip
import hashlib
import statistics
import subprocess
import tempfile
import unittest
import uuid
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pod5

SELECT = Path(HERE) / "select_reads.py"
OFFSET, SCALE, MEDIAN_BEFORE = 10.0, 0.2, 200.0


def run_info(acq, rate):
    t = datetime(2024, 1, 1, tzinfo=timezone.utc)
    return pod5.RunInfo(acquisition_id=acq, acquisition_start_time=t, adc_max=4095, adc_min=-4096,
                        context_tags={}, experiment_name="test", flow_cell_id="FC", flow_cell_product_code="P",
                        protocol_name="p", protocol_run_id="r", protocol_start_time=t, sample_id="s",
                        sample_rate=rate, sequencing_kit="k", sequencer_position="x", sequencer_position_type="y",
                        software="t", system_name="n", system_type="z", tracking_id={})


def signal(kind, n, rng):
    """int16 samples: 'clean' = levels of ~10 samples with a little noise, 'noisy' = white noise, 'open' = clean with
    an open-pore-like stretch (>= 0.85 x median_before) of 60 samples in the middle"""
    if kind == "noisy":
        x = rng.normal(500, 100, n)
    else:
        x = np.repeat(rng.normal(500, 60, n // 10 + 1), 10)[:n] + rng.normal(0, 5, n)
        if kind == "open":
            x[n // 2:n // 2 + 60] = 220 / SCALE - OFFSET
    return np.clip(np.round(x), -4000, 4000).astype(np.int16)


def make_reads(spec, rng, ri):
    """spec: (kind, samples) per read"""
    reads = []
    for i, (kind, n) in enumerate(spec):
        reads.append(pod5.Read(read_id=uuid.UUID(int=int(rng.integers(0, 2**62)) * 7919 + i),
                               signal=signal(kind, n, rng),
                               pore=pod5.Pore(channel=1 + i % 50, well=1, pore_type="pore"),
                               calibration=pod5.Calibration(offset=OFFSET, scale=SCALE),
                               end_reason=pod5.EndReason(pod5.EndReasonEnum.SIGNAL_POSITIVE, False),
                               read_number=i, run_info=ri, median_before=MEDIAN_BEFORE, start_sample=i * 1000))
    return reads


def write(path, reads):
    with pod5.Writer(path) as w:
        for r in reads:
            w.add_read(r)


def md5(rid):
    return hashlib.md5(rid.encode()).hexdigest()


def ref_length(n, rate, bps=400, digits=0):
    if digits:
        return round(float(f"{n / rate:.{digits}g}") * bps)
    return round(n * bps / rate)


def ref_features(r):
    """noise and open-pore run as select_reads.py's docstring defines them, written independently"""
    off, scale = np.float32(r.calibration.offset), np.float32(r.calibration.scale)
    x = ((r.signal.astype(np.float32) + off) * scale).astype(np.float64)
    q = float(f"{np.median(np.abs(np.diff(x))) / x.std():.6g}")
    core = x[2000:-2000] if len(x) > 4010 else x
    best = run = 0
    for v in core >= 0.85 * r.median_before:
        run = run + 1 if v else 0
        best = max(best, run)
    return q, best


def ref_select(reads, G, cov, min_len=1000, pool_cov=None, noise_max=1.15, open_pore=0, digits=0):
    """the rules of select_reads.py's docstring: (kept ids in md5 order, pool ids, statuses)"""
    length = {str(r.read_id): ref_length(len(r.signal), r.run_info.sample_rate, digits=digits) for r in reads}
    byid = {str(r.read_id): r for r in reads}
    long_ = sorted((i for i, n in length.items() if n >= min_len), key=md5)
    pool, acc = [], 0
    for i in long_:
        if pool_cov is not None and acc >= pool_cov * G:
            break
        pool.append(i)
        acc += length[i]
    feat = {i: ref_features(byid[i]) for i in pool}
    cut = noise_max * statistics.median(feat[i][0] for i in pool)
    kept, acc, status = [], 0, {}
    for i in pool:
        if acc >= cov * G:
            status[i] = "unused"
        elif feat[i][0] > cut:
            status[i] = "noise"
        elif open_pore and feat[i][1] > open_pore:
            status[i] = "open_pore"
        else:
            status[i] = "kept"
            kept.append(i)
            acc += length[i]
    return kept, pool, status, length


def py(*args, check=True):
    return subprocess.run([sys.executable, *map(str, args)], capture_output=True, text=True, check=check)


def read_gz(path):
    return gzip.decompress(Path(path).read_bytes()).decode().split()


def table(path):
    rows = [line.split("\t") for line in Path(path).read_text().splitlines()]
    return {r[0]: dict(zip(rows[0], r)) for r in rows[1:]}


class Select(unittest.TestCase):
    G = 20000

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        d = cls.d = Path(cls.tmp.name)
        rng = np.random.default_rng(11)
        kinds = ["clean"] * 7 + ["noisy", "open"]
        spec_a = [(kinds[i % 9], int(rng.integers(5000, 40000))) for i in range(160)]
        spec_b = [(kinds[i % 9], int(rng.integers(5000, 30000))) for i in range(60)]
        spec_b += [("clean", 9998), ("clean", 9999)]  # 999.8 and 999.9 bases: 1,000 when rounded
        spec_b += [("clean", 40005), ("clean", 40035), ("clean", 40045)]  # 6-digit durations round these otherwise
        cls.a = make_reads(spec_a, rng, run_info("acq-a", 5000))
        cls.b = make_reads(spec_b, rng, run_info("acq-b", 4000))
        cls.kind = {str(r.read_id): k for r, (k, _) in zip(cls.a + cls.b, spec_a + spec_b)}
        (d / "src" / "x").mkdir(parents=True)
        write(d / "src" / "x" / "a.pod5", cls.a)
        write(d / "src" / "b.pod5", cls.b)
        cls.reads = cls.a + cls.b

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_select(self, name, *opts, src=None, check=True):
        out = self.d / name
        res = py(SELECT, *(src or ["--pod5", self.d / "src"]), "-G", self.G, "--min-length", 1000, "--out", out,
                 "--threads", 2, *opts, check=check)
        return out, res

    def test_1_rules(self):
        out, res = self.run_select("k5", "--coverage", 5, "--pool-coverage", 8)
        kept, pool, status, length = ref_select(self.reads, self.G, 5, pool_cov=8)
        self.assertEqual(read_gz(out / "reads.txt.gz"), kept)
        t = table(out / "selection.tsv")
        self.assertEqual(list(t), pool)
        self.assertEqual({i: r["status"] for i, r in t.items()}, status)
        self.assertEqual({i: int(r["length"]) for i, r in t.items()}, {i: length[i] for i in pool})
        self.assertTrue(all(status[i] in ("noise", "unused") for i in pool if self.kind[i] == "noisy"))
        self.assertIn("unused", status.values())
        # the target is reached and the read reaching it is kept
        acc = sum(length[i] for i in kept)
        self.assertGreaterEqual(acc, 5 * self.G)
        self.assertLess(acc - length[kept[-1]], 5 * self.G)
        params = dict(line.split("\t") for line in (out / "params.tsv").read_text().splitlines()[1:])
        self.assertEqual(int(params["kept"]), len(kept))
        self.assertEqual(int(params["kept_bases"]), acc)
        self.assertEqual(int(params["pool_reads"]), len(pool))
        self.assertEqual(res.stdout.split("\t")[0], str(len(kept)))
        # gzip -9 -n: no name, no time stamp
        raw = (out / "reads.txt.gz").read_bytes()
        self.assertEqual(raw[4:8], b"\0\0\0\0")
        self.assertEqual(raw[3] & 0x08, 0)

    def test_2_open_pore(self):
        out, _ = self.run_select("b5", "--coverage", 5, "--open-pore", 4)
        kept, pool, status, _ = ref_select(self.reads, self.G, 5, open_pore=4)
        self.assertEqual(read_gz(out / "reads.txt.gz"), kept)
        self.assertIn("open_pore", status.values())
        t = table(out / "selection.tsv")
        self.assertTrue(all(int(t[i]["open_pore_run"]) >= 60 for i in pool if self.kind[i] == "open"))
        self.assertTrue(all(status[i] != "kept" for i in pool if self.kind[i] == "open"))

    def test_3_table_reselect(self):
        # a 3x set from the 8x pool's table equals the 3x set read from the POD5s
        big, _ = self.run_select("t8", "--coverage", 8, "--pool-coverage", 8)
        out, _ = self.run_select("t3", "--coverage", 3, "--pool-coverage", 8, src=["--table", big / "selection.tsv"])
        ref, _ = self.run_select("p3", "--coverage", 3, "--pool-coverage", 8)
        self.assertEqual(read_gz(out / "reads.txt.gz"), read_gz(ref / "reads.txt.gz"))
        self.assertEqual((out / "selection.tsv").read_text(), (ref / "selection.tsv").read_text())

    def test_4_pool_list_and_missing(self):
        ids = sorted(str(r.read_id) for r in self.reads)[::2]
        lst = self.d / "pool.txt.gz"
        lst.write_bytes(gzip.compress(("\n".join(ids) + "\n").encode()))
        out, _ = self.run_select("l3", "--coverage", 3, "--pool-ids", lst)
        sub = [r for r in self.reads if str(r.read_id) in set(ids)]
        self.assertEqual(read_gz(out / "reads.txt.gz"), ref_select(sub, self.G, 3)[0])
        bad = self.d / "bad.txt"
        bad.write_text("\n".join(ids + [str(uuid.uuid4())]) + "\n")
        _, res = self.run_select("never", "--coverage", 3, "--pool-ids", bad, check=False)
        self.assertNotEqual(res.returncode, 0)
        self.assertFalse((self.d / "never" / "reads.txt.gz").exists())
        out, res = self.run_select("l3m", "--coverage", 3, "--pool-ids", bad, "--missing-ok")
        self.assertEqual(res.returncode, 0)
        self.assertEqual(read_gz(out / "reads.txt.gz"), read_gz(self.d / "l3" / "reads.txt.gz"))

    def test_5_summary_and_duration_digits(self):
        # a sequencing summary: only reads with duration x 400 >= 4000 bases
        summ = self.d / "sequencing_summary.txt"
        rows = ["filename\tread_id\tduration"]
        rows += [f"f\t{r.read_id}\t{len(r.signal) / r.run_info.sample_rate:.4f}" for r in self.reads]
        summ.write_text("\n".join(rows) + "\n")
        out, _ = self.run_select("s3", "--coverage", 3, "--summary", summ, "--summary-min-length", 4000)
        sub = [r for r in self.reads if len(r.signal) / r.run_info.sample_rate * 400 >= 4000]
        self.assertEqual(read_gz(out / "reads.txt.gz"), ref_select(sub, self.G, 3)[0])
        # lengths from the duration at 6 significant digits differ from the exact ones at 4 kHz
        out, _ = self.run_select("d3", "--coverage", 3, "--duration-digits", 6)
        t = table(out / "selection.tsv")
        want = {str(r.read_id): ref_length(len(r.signal), r.run_info.sample_rate, digits=6) for r in self.reads}
        self.assertEqual({i: int(r["length"]) for i, r in t.items()}, {i: want[i] for i in t})
        exact = {str(r.read_id): ref_length(len(r.signal), r.run_info.sample_rate) for r in self.reads}
        self.assertNotEqual(want, exact)
        self.assertEqual(read_gz(out / "reads.txt.gz"), ref_select(self.reads, self.G, 3, digits=6)[0])

    def test_6_pool_out(self):
        lst = self.d / "long.txt.gz"
        res = py(SELECT, "--pod5", self.d / "src" / "x" / "a.pod5", self.d / "src" / "b.pod5", "--min-length", 2000,
                 "--pool-out", lst)
        self.assertEqual(res.returncode, 0)
        want = sorted(str(r.read_id) for r in self.reads if ref_length(len(r.signal), r.run_info.sample_rate) >= 2000)
        self.assertEqual(read_gz(lst), want)


if __name__ == "__main__":
    unittest.main(verbosity=2)
