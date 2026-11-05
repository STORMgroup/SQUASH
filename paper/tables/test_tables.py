#!/usr/bin/env python3
"""paper/tables/test_tables.py: tests of lib.py's readers, formats and ranking on a small OUT_DIR made here (no data
needed). Run: python3 tables/test_tables.py"""
import os
import sys
import tempfile
import unittest
from pathlib import Path

TMP = tempfile.TemporaryDirectory()
os.environ["OUT_DIR"] = TMP.name
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib  # noqa: E402

TIME = ("\tUser time (seconds): {u}\n\tSystem time (seconds): {s}\n\tPercent of CPU this job got: {p}%\n"
        "\tElapsed (wall clock) time (h:mm:ss or m:ss): {w}\n\tMaximum resident set size (kbytes): {r}\n"
        "\tExit status: {e}\n")


def put(rel, text):
    p = Path(TMP.name) / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    return p


class Readers(unittest.TestCase):
    def test_time(self):
        p = put("a.time.txt", TIME.format(u=10.5, s=1.5, p=400, w="1:02:03.50", r=2048, e=0))
        self.assertEqual(lib.time_v(p), (3723.5, 12.0, 2048))
        self.assertEqual(lib.cpu_share(p), 4.0)
        put("b.time.txt", TIME.format(u=1, s=1, p=100, w="0:10", r=1, e=1))
        self.assertIsNone(lib.time_v(Path(TMP.name) / "b.time.txt"))          # a failed step is not a time
        self.assertIsNone(lib.time_v(Path(TMP.name) / "none.time.txt"))

    def test_total_and_routes(self):
        for n, (w, r) in {"rawsamble_index": ("10:00", 100), "rawsamble_overlap": ("20:00", 300), "miniasm": ("0:30", 50)}.items():
            put(f"assembly/k1/rawsamble/logs/{n}.time.txt", TIME.format(u=60, s=0, p=100, w=w, r=r, e=0))
        self.assertEqual(lib.tool_time("k1", "rawsamble"), (1830.0, 180.0, 300))
        self.assertTrue(lib.ran("k1", "rawsamble"))
        self.assertFalse(lib.ran("k1", "flye"))
        put("timing/k1/basecall_hac_cpu.time.txt", TIME.format(u=1, s=0, p=6400, w="1:00:00", r=5000, e=0))
        put("timing/k1/basecall_hac_cpu.tsv", "key\tsamples_timed\tdorado_s\tset_samples\tset_est_s\nk1\t10\t3600\t40\t14400\n")
        self.assertEqual(lib.hac_cpu("k1"), (14400.0, 14400.0 * 64, 5000))

    def test_quast(self):
        put("eval/k1/squash/quast/report.tsv",
            "Assembly\tSQUASH (spelled; forward-strand unitigs)\tSQUASH (spelled; reverse-strand unitigs)\n"
            "NGA50\t5128000\t-\n# misassemblies\t1234\t0\nGenome fraction (%)\t99.95\t0.2\n")
        rep = lib.quast(Path(TMP.name) / "eval/k1/squash/quast/report.tsv")
        fw, rv = lib.squash_column(rep, "spelled", "forward"), lib.squash_column(rep, "spelled", "reverse")
        self.assertEqual(lib.strand_cells(rep, [fw, rv]), ["5.13\\,/\\,--", "1,234\\,/\\,0", "100.0\\,/\\,0.2"])
        self.assertEqual(lib.quast_cell(rep, "Flye", "NGA50"), lib.NA)


class Formats(unittest.TestCase):
    def test_formats(self):
        self.assertEqual(lib.hms(370.5), "0:06:10")       # rounded half to even, as the paper's tables
        self.assertEqual(lib.hms(90061), "25:01:01")
        self.assertEqual(lib.gb(3113012), "3.1")
        self.assertEqual(lib.cpu_h(3877920), "1,077.2")
        self.assertEqual(lib.mb(5151986), "5.15")

    def test_rank(self):
        self.assertEqual(lib.rank(["204", "208"], low=True), [0])
        self.assertEqual(lib.rank(["1.2", "1.2"], low=True), [0, 1])
        self.assertEqual(lib.rank(["3.47", lib.NA], low=False), [])
        self.assertEqual(lib.rank(["0610", "2526", "13749"], low=True), [0])

    def test_contiguity_best(self):
        head = "key\tassembler\tsegments\tlargest_bp\tNG50\tlargest_component_bp\n"
        put("eval/k2/squash/tier1.tsv", head + "k2\tSQUASH\t204\t5153678\t5151986\t5163717\n")
        put("eval/k2/compare/tier1.tsv", head + "k2\tRawsamble + miniasm\t208\t281722\t0\t9834394\n"
            "k2\tnanoMDBG\t509\t5132227\t5132227\t-\n")
        (Path(TMP.name) / "assembly/k2/rawsamble").mkdir(parents=True)
        rows = lib.contiguity_cells("k2", lib.Values("test"))
        self.assertEqual(rows[0], ["\\best{204}", "\\best{5.15}", "\\best{5.15}", "5.16"])
        self.assertEqual(rows[1], ["208", "0.28", "--", "\\best{9.83}"])     # NG50 0 = undefined, and loses
        self.assertEqual(rows[2], ["509", "5.13", "5.13", "5.13"])           # no links: the longest contig
        self.assertEqual(rows[3], [lib.NA] * 4)


if __name__ == "__main__":
    try:
        unittest.main()
    finally:
        TMP.cleanup()
