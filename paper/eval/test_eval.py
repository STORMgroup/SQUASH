#!/usr/bin/env python3
"""Tests of the evaluation's scripts on small made-up inputs (numpy; pysam and Remora are not needed).
Usage: test_eval.py (exit 0 when every test passes)"""
import gzip
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common  # noqa: E402
import mask_arrays  # noqa: E402
import mask_dmel  # noqa: E402
import project  # noqa: E402
import spell  # noqa: E402
import tier1  # noqa: E402


def spans_record(name, starts, ends, samples=0):
    n = len(starts)
    return (struct.pack("<I", len(name)) + name.encode() + struct.pack("<II", samples, n)
            + np.asarray(starts, "<u4").tobytes() + np.asarray(ends, "<u4").tobytes())


class Files(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_spans_stream_and_cut(self):
        p = self.d / "spans.bin"
        p.write_bytes(spans_record("a", [0, 10, 20], [10, 20, 30]) + spans_record("b", [5], [9])
                      + spans_record("a", [1], [2]))
        sp = common.event_spans(p, {"a"})
        self.assertEqual(list(sp), ["a"])
        self.assertEqual(sp["a"][0].tolist(), [1])        # a read listed twice: its last record
        self.assertEqual(common.cut_spans(p, {"b"}, self.d / "cut.bin"), 1)
        self.assertEqual(common.event_spans(self.d / "cut.bin")["b"][1].tolist(), [9])

    def test_graph_reads(self):
        g = self.d / "input"
        g.mkdir()
        np.asarray([1, 2, 3, 4, 5], "<u4").tofile(g / "tokens.u32le")
        np.asarray([[0, 9], [3, 12], [6, 15], [0, 9], [2, 11]], "<u8").tofile(g / "token_spans.u64le")
        with gzip.open(g / "reads.tsv.gz", "wt") as f:
            f.write("read_id\tread_ordinal\tq_symbols\ttoken_offset\ttoken_count\nr0\t0\t20\t0\t3\nr1\t1\t15\t3\t2\n")
        qs = common.graph_reads(g, {"r1"})
        self.assertEqual(list(qs), ["r1"])
        self.assertEqual(qs["r1"].tolist(), [[0, 9], [2, 11]])


class Spell(unittest.TestCase):
    def test_moves_blocks(self):
        m = spell.Moves("ACGT", 5, np.array([1, 0, 1, 1, 0, 1], np.int8), ts=10, sp=0)
        self.assertEqual((m.lo, m.hi), (10, 40))
        self.assertEqual(m.bases(10, 40), "ACGT")
        self.assertEqual(m.bases(20, 30), "CG")          # blocks 2 and 3 start bases 2 and 3
        self.assertEqual(spell.read_bases([m], 0, 15), "N" * 1 + "A")   # 10 uncovered samples at 4/30 bases per sample

    def test_piece_ends_where_the_next_takes_over(self):
        rows = [{"unitig": "0", "rank": "0", "circular": "0", "read_id": "a", "read_token_to": "2", "token_from": "0",
                 "token_to": "4", "sample_from": "100", "sample_to": "400"},
                {"unitig": "0", "rank": "1", "circular": "0", "read_id": "b", "read_token_to": "9", "token_from": "4",
                 "token_to": "8", "sample_from": "50", "sample_to": "250"}]
        qs = {"a": np.array([[0, 9], [3, 12], [7, 16]]), "b": np.array([[0, 9]])}
        st = {"a": np.arange(0, 2000, 10), "b": np.arange(0, 2000, 10)}
        # row 0 ends at the first sample of read a's token 2 (q-symbol 7 -> sample 70, but not before its start 100)
        self.assertEqual(spell.piece_samples(rows, 0, False, qs.__getitem__, st.__getitem__), (100, 100))
        qs["a"] = np.array([[0, 9], [3, 12], [30, 39]])
        self.assertEqual(spell.piece_samples(rows, 0, False, qs.__getitem__, st.__getitem__), (100, 300))
        # the last piece of a linear unitig: to the end of its last token
        self.assertEqual(spell.piece_samples(rows, 1, False, qs.__getitem__, st.__getitem__), (50, 250))

    def test_gap_sized_by_median_bases_per_token(self):
        rows = [{"unitig": "3", "rank": "0", "circular": "0", "read_id": "a", "read_token_to": "9", "token_from": "0",
                 "token_to": "2", "sample_from": "0", "sample_to": "20", "verified": "1"},
                {"unitig": "3", "rank": "1", "circular": "0", "read_id": "-", "read_token_to": "-1", "token_from": "2",
                 "token_to": "5", "sample_from": "-1", "sample_to": "-1"}]
        m = spell.Moves("ACGT", 5, np.array([1, 1, 1, 1], np.int8), ts=0, sp=0)
        out, table, bpt = spell.spell(rows, {"a": np.zeros((1, 2))}.__getitem__, lambda r: np.zeros(1), {"a": [m]})
        self.assertEqual(bpt, 2.0)
        self.assertEqual(out[3], ("ACGT" + "N" * 6, False, 1, 1))


class Project(unittest.TestCase):
    def test_blocks_and_bridged_ns(self):
        pr = project.Projector(500)
        for p in (100, 101, 102):
            pr.add_event("c", "+", p)
        pr.add_n(5, "unplaced")
        pr.add_event("c", "+", 110)                       # continues the block across the Ns
        pr.add_event("d", "+", 7)                         # another contig: a new block
        pr.drop_bridged_ns()
        self.assertEqual([it[0] for it in pr.items], ["B", "B"])
        ref = {"c": "A" * 100 + "ACGTACGTACGT", "d": "T" * 8}
        self.assertEqual(pr.sequence(ref), "ACGTACGTACG" + "T")

    def test_reverse_strand(self):
        pr = project.Projector(500)
        pr.add_event("c", "-", 2)
        pr.add_event("c", "-", 0)
        self.assertEqual(pr.sequence({"c": "AACG"}), "GTT")


class Tier1(unittest.TestCase):
    def test_metrics(self):
        row = tier1.metrics("x", "bp", {"a": 60.0, "b": 30.0, "c": 10.0}, 200, {("a", "+", "b", "+")}, {"c"})
        self.assertEqual((row["N50"], row["NG50"], row["largest_bp"], row["total_bp"]), (60, 10, 60, 100))
        self.assertEqual(row["auN"], round((3600 + 900 + 100) / 100))
        self.assertEqual(row["auNG"], round((60 * 60 + 30 * 30 + 10 * 10) / 200))
        self.assertEqual((row["components"], row["largest_component_bp"], row["largest_component_share"]), (2, 90, "0.9000"))

    def test_gfa_links_canonical(self):
        with tempfile.NamedTemporaryFile("w", suffix=".gfa", delete=False) as f:
            f.write("S\t1\t*\tLN:i:5\nS\t2\t*\tLN:i:7\nL\t1\t+\t2\t+\t31M\nL\t2\t-\t1\t-\t31M\nL\t2\t+\t2\t+\t0M\n")
        names, links, circ = tier1.read_gfa_graph(Path(f.name))
        Path(f.name).unlink()
        self.assertEqual((names, len(links), circ), (["1", "2"], 2, {"2"}))

    def test_rawsamble_units(self):
        self.assertEqual(tier1.piece_units(10, 50, "+", 30), (10, 40))
        self.assertEqual(tier1.piece_units(10, 50, "-", 30), (20, 50))


class Masks(unittest.TestCase):
    def test_arrays_found_and_kept_apart(self):
        rng = np.random.default_rng(1)
        unit = "".join(rng.choice(list("ACGT"), 180))
        flank = lambda n: "".join(rng.choice(list("ACGT"), n))  # noqa: E731
        seq = flank(300_000) + unit * 1200 + flank(300_000)      # a 216-kb array of a 180-bp unit
        iv = mask_arrays.arrays(seq)
        self.assertEqual(len(iv), 1)
        a, b = iv[0]
        self.assertTrue(290_000 <= a <= 310_000 and 500_000 <= b <= 530_000, iv)

    def test_dmel_regions(self):
        length = {"chr2L": 1000, "chr4": 100, "chrY": 50}
        band = {"chr2L": (100, 800)}
        rr = {"chr2L": [(30000, 30010)], "chr4": [(40, 45)]}
        reg = mask_dmel.regions(length, band, rr)
        self.assertEqual(reg["chr2L"], [(0, 100, "tel"), (100, 800, "eu"), (800, 1000, "het")])
        self.assertEqual(reg["chr4"], [(0, 100, "rdna")])  # within 20 kb of an rRNA element
        self.assertEqual(reg["chrY"], [(0, 50, "chrY")])


class Scripts(unittest.TestCase):
    def test_strands_and_closed(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            (d / "asm").mkdir()
            (d / "strands").mkdir()
            (d / "asm/squash_spelled.fa").write_text(">unitig_0 LN:i:100 circular=1\n" + "A" * 100 + "\n"
                                                       ">unitig_1 LN:i:50 circular=0\n" + "C" * 50 + "\n"
                                                       ">unitig_2 LN:i:10 circular=0\n" + "G" * 10 + "\n")
            (d / "asm/squash_placed.fa").write_text(">unitig_0\nA\n>unitig_1\nC\n>unitig_2\nG\n")
            paf = d / "strands/spelled.paf"
            paf.write_text("unitig_0\t100\t0\t100\t-\tchr\t100\t0\t100\t99\t100\t60\n"
                           "unitig_1\t50\t0\t50\t+\tchr\t100\t0\t50\t49\t50\t60\n")
            subprocess.run([sys.executable, "-B", str(HERE / "strands.py"), str(paf), str(d / "asm"), str(d / "strands")], check=True)
            self.assertEqual((d / "strands/squash_placed_rev.fa").read_text(), ">unitig_0\nA\n")
            self.assertEqual((d / "strands/squash_spelled_fwd.fa").read_text().count(">"), 1)
            summary = (d / "strands/summary.tsv").read_text().splitlines()
            self.assertEqual(summary[1:], ["all\t3\t160", "forward\t1\t50", "reverse\t1\t100", "neither\t1\t10"])
            out = subprocess.run([sys.executable, "-B", str(HERE / "closed.py"), "k", str(d)], check=True, capture_output=True,
                                 text=True).stdout.splitlines()
            self.assertEqual(out[1], "k\tchr\t100\tno\tyes\tunitig_0:100")


if __name__ == "__main__":
    unittest.main()
