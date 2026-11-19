"""Small scientific sanity checks, runnable without network access."""

import unittest

from Bio import Align
from Bio.Seq import Seq

from src.align_analysis import (annotated_cds, gc_percent, identity_percent,
                                protein_aligner, scan_orfs, sliding_gc)


class SequenceTests(unittest.TestCase):
    def test_gc_percentage(self):
        self.assertEqual(gc_percent("GGAT"), 50.0)

    def test_gc_ignores_ambiguous(self):
        self.assertEqual(gc_percent("GGNNAT"), 50.0)

    def test_empty_sequence_rejected(self):
        with self.assertRaises(ValueError):
            gc_percent("")

    def test_sliding_windows(self):
        points = sliding_gc("GGGGAAAA", window=4, step=4)
        self.assertEqual([p["gc_percent"] for p in points], [100.0, 0.0])
        self.assertEqual([p["position_nt"] for p in points], [2.5, 6.5])

    def test_invalid_window(self):
        with self.assertRaises(ValueError):
            sliding_gc("ATCG", window=0)

    def test_orf_scanning(self):
        orfs = scan_orfs("CCCATGAAACCCGGGTAA", min_aa=4)
        self.assertEqual(len(orfs), 1)
        self.assertEqual(orfs[0]["start_nt"], 4)
        self.assertEqual(orfs[0]["end_nt"], 18)
        self.assertEqual(orfs[0]["length_aa_excluding_stop"], 4)

    def test_orf_no_stop_not_counted(self):
        self.assertEqual(scan_orfs("ATGAAACCC", min_aa=1), [])


class AlignmentTests(unittest.TestCase):
    def test_identical_seqs_are_100_percent(self):
        alignment = protein_aligner().align("MAKPG", "MAKPG")[0]
        self.assertEqual(identity_percent(alignment), 100.0)

    def test_mismatched_seqs_have_less_than_100_percent(self):
        alignment = protein_aligner().align("MAKPG", "MAQPG")[0]
        self.assertLess(identity_percent(alignment), 100.0)



    def test_identity_with_and_without_gaps(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 2
        aligner.mismatch_score = -5
        aligner.open_gap_score = -1

        alignment = aligner.align("ACG", "ACTG")[0]

        self.assertEqual(identity_percent(alignment), 75.0)
        self.assertEqual(
            identity_percent(alignment, include_gaps=False),
            100.0
        )
    if __name__ == "__main__":
        unittest.main()