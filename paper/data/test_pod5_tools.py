#!/usr/bin/env python3
"""Tests of pod5_cut.py and pod5_digest.py on synthetic POD5s written here (no data file; pod5 0.3.39).
Run: python3 test_pod5_tools.py  (or python3 -m unittest test_pod5_tools)."""
import gzip
import hashlib
import random
import struct
import subprocess
import sys
import tempfile
import unittest
import uuid
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pod5

HERE = Path(__file__).resolve().parent
CUT, DIGEST = HERE / "pod5_cut.py", HERE / "pod5_digest.py"


def run_info(acq, rate):
    t = datetime(2024, 1, 1, tzinfo=timezone.utc)
    return pod5.RunInfo(acquisition_id=acq, acquisition_start_time=t, adc_max=4095, adc_min=-4096,
                        context_tags={}, experiment_name="test", flow_cell_id="FC", flow_cell_product_code="P",
                        protocol_name="p", protocol_run_id="r", protocol_start_time=t, sample_id="s",
                        sample_rate=rate, sequencing_kit="k", sequencer_position="x", sequencer_position_type="y",
                        software="t", system_name="n", system_type="z", tracking_id={})


def make_reads(n, rng, acqs):
    reads = []
    for i in range(n):
        ri = acqs[i % len(acqs)]
        sig = rng.integers(-2000, 2000, size=int(rng.integers(500, 250000)), dtype=np.int16)
        reads.append(pod5.Read(read_id=uuid.UUID(int=int(rng.integers(0, 2**62)) * 7919 + i), signal=sig,
                               pore=pod5.Pore(channel=1 + i % 50, well=1, pore_type="pore"),
                               calibration=pod5.Calibration(offset=float(rng.uniform(-300, 0)), scale=float(rng.uniform(0.1, 0.3))),
                               end_reason=pod5.EndReason(pod5.EndReasonEnum.SIGNAL_POSITIVE, False),
                               read_number=i, run_info=ri, median_before=200.0, start_sample=i * 1000))
    return reads


def write(path, reads):
    with pod5.Writer(path) as w:
        for r in reads:
            w.add_read(r)


def ref_digest(reads):
    """the digest as pod5_digest.py's docstring defines it, written independently"""
    h = hashlib.sha256(b"squash-pod5-digest-v1\n")
    for r in reads:
        h.update(str(r.read_id).encode())
        h.update(struct.pack("<f", r.calibration.offset) + struct.pack("<f", r.calibration.scale))
        h.update(struct.pack("<H", r.run_info.sample_rate) + struct.pack("<Q", len(r.signal)))
        h.update(r.signal.astype("<i2").tobytes())
    return h.hexdigest()


def py(*args, check=True):
    return subprocess.run([sys.executable, *map(str, args)], capture_output=True, text=True, check=check)


class Pod5Tools(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        d = cls.d = Path(cls.tmp.name)
        rng = np.random.default_rng(7)
        acqs = [run_info("acq-a", 5000), run_info("acq-b", 4000)]
        cls.a = make_reads(1300, rng, acqs[:1])          # > 1000 reads: several read batches
        cls.b = make_reads(700, rng, acqs)               # two run infos, one of them new
        (d / "src" / "x").mkdir(parents=True)
        write(d / "src" / "x" / "a.pod5", cls.a)
        write(d / "src" / "b.pod5", cls.b)
        byid = {str(r.read_id): r for r in cls.a + cls.b}
        ids = list(byid)
        random.Random(3).shuffle(ids)
        cls.want = ids[:900]                              # a list in an order no source has
        cls.want_reads = [byid[x] for x in cls.want]
        (d / "list.txt.gz").write_bytes(gzip.compress(("\n".join(cls.want) + "\n").encode()))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_1_digest_matches_reference(self):
        out = py(DIGEST, self.d / "src" / "b.pod5").stdout.split()
        self.assertEqual(out[0], "700")
        self.assertEqual(out[2], ref_digest(self.b))

    def test_2_cut_order_and_digest(self):
        dst = self.d / "cut.pod5"
        py(CUT, "cut", self.d / "list.txt.gz", dst, self.d / "src", "--procs", "2")
        with pod5.Reader(dst) as r:
            self.assertEqual(r.read_ids, self.want)
        res = py(DIGEST, "--ids", self.d / "list.txt.gz", "--expect", ref_digest(self.want_reads), dst)
        self.assertEqual(res.returncode, 0)
        self.assertFalse((self.d / "cut.pod5.parts").exists())

    def test_3_part_and_join(self):
        parts = []
        for i, src in enumerate([self.d / "src" / "x" / "a.pod5", self.d / "src" / "b.pod5"]):
            parts.append(self.d / f"p{i}.pod5")
            py(CUT, "part", self.d / "list.txt.gz", src, parts[-1])
        dst = self.d / "joined.pod5"
        py(CUT, "join", self.d / "list.txt.gz", dst, *parts[::-1])
        self.assertEqual(py(DIGEST, dst).stdout.split()[2], ref_digest(self.want_reads))

    def test_4_checks_fail(self):
        dst = self.d / "cut.pod5"
        if not dst.exists():
            py(CUT, "cut", self.d / "list.txt.gz", dst, self.d / "src")
        bad = self.d / "bad.txt"
        bad.write_text("\n".join([self.want[1], self.want[0]] + self.want[2:]) + "\n")
        res = py(DIGEST, "--ids", bad, dst, check=False)
        self.assertEqual(res.returncode, 2)
        self.assertIn("at read 1", res.stderr)
        res = py(DIGEST, "--expect", "0" * 64, dst, check=False)
        self.assertEqual(res.returncode, 1)
        missing = self.d / "missing.txt"
        missing.write_text(self.want[0] + "\n" + str(uuid.uuid4()) + "\n")
        res = py(CUT, "cut", missing, self.d / "never.pod5", self.d / "src", check=False)
        self.assertNotEqual(res.returncode, 0)
        self.assertFalse((self.d / "never.pod5").exists())

    def test_5_select(self):
        # the list's reads straight from the sources, in source order (a, then b, each in file order)
        order = [str(r.read_id) for r in self.a + self.b if str(r.read_id) in set(self.want)]
        lst = self.d / "src_order.txt"
        lst.write_text("\n".join(order) + "\n")
        byid = {str(r.read_id): r for r in self.a + self.b}
        res = py(DIGEST, "--select", lst, self.d / "src" / "x" / "a.pod5", self.d / "src" / "b.pod5")
        self.assertEqual(res.stdout.split()[2], ref_digest([byid[x] for x in order]))
        res = py(DIGEST, "--select", self.d / "list.txt.gz", self.d / "src" / "x" / "a.pod5", self.d / "src" / "b.pod5",
                 check=False)
        self.assertEqual(res.returncode, 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
