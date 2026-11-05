#!/usr/bin/env python3
"""Tests of cpu_subset.py on a synthetic POD5 written here (no data file; pod5 0.3.39).
Run: python3 test_cpu_subset.py  (or python3 -m unittest test_cpu_subset)."""
import hashlib
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
SUBSET = HERE / "cpu_subset.py"


def write_pod5(path, n, seed=1):
    rng = np.random.default_rng(seed)
    t = datetime(2024, 1, 1, tzinfo=timezone.utc)
    ri = pod5.RunInfo(acquisition_id="acq", acquisition_start_time=t, adc_max=4095, adc_min=-4096, context_tags={},
                      experiment_name="test", flow_cell_id="FC", flow_cell_product_code="P", protocol_name="p",
                      protocol_run_id="r", protocol_start_time=t, sample_id="s", sample_rate=5000,
                      sequencing_kit="k", sequencer_position="x", sequencer_position_type="y", software="t",
                      system_name="n", system_type="z", tracking_id={})
    lengths = {}
    with pod5.Writer(path) as w:
        for i in range(n):
            rid = uuid.UUID(int=int(rng.integers(0, 2**62)) * 7919 + i)
            sig = rng.integers(-2000, 2000, size=int(rng.integers(100, 5000)), dtype=np.int16)
            w.add_read(pod5.Read(read_id=rid, signal=sig, pore=pod5.Pore(channel=1, well=1, pore_type="pore"),
                                 calibration=pod5.Calibration(offset=0.0, scale=0.2),
                                 end_reason=pod5.EndReason(pod5.EndReasonEnum.SIGNAL_POSITIVE, False),
                                 read_number=i, run_info=ri, median_before=200.0, start_sample=0))
            lengths[str(rid)] = len(sig)
    return lengths


def run(key, target, pod, out):
    subprocess.run([sys.executable, str(SUBSET), key, target, str(pod), str(out)], check=True, capture_output=True)
    ids = (out / f"{key}.cpu_subset.ids.txt").read_text().split()
    rows = [l.split("\t") for l in (out / f"{key}.cpu_subset.tsv").read_text().splitlines()]
    return ids, dict(zip(rows[0], rows[1]))


class CpuSubsetTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)
        self.pod = self.d / "set.pod5"
        self.lengths = write_pod5(self.pod, 60)
        self.order = sorted(self.lengths, key=lambda x: hashlib.md5(x.encode()).hexdigest())

    def tearDown(self):
        self.tmp.cleanup()

    def test_all(self):
        ids, row = run("k", "all", self.pod, self.d)
        self.assertEqual(ids, self.order)
        self.assertEqual(int(row["subset_samples"]), sum(self.lengths.values()))
        self.assertEqual(int(row["set_reads"]), 60)

    def test_target_keeps_the_crossing_read(self):
        target = sum(self.lengths[r] for r in self.order[:10]) - 1
        ids, row = run("k", str(target), self.pod, self.d)
        self.assertEqual(ids, self.order[:10])
        self.assertEqual(int(row["subset_samples"]), target + 1)
        self.assertEqual(row["sample_rate"], "5000")

    def test_target_on_a_boundary(self):
        target = sum(self.lengths[r] for r in self.order[:10])
        ids, _ = run("k", str(target), self.pod, self.d)
        self.assertEqual(ids, self.order[:10])


if __name__ == "__main__":
    unittest.main()
