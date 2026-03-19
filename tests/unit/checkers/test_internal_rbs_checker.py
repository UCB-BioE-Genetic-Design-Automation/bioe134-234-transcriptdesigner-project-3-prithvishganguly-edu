"""
test_internal_rbs_checker.py

Unit tests for InternalRBSChecker.

Run with:
    pytest tests/unit/checkers/test_internal_rbs_checker.py -v
"""

import pytest
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker


@pytest.fixture(scope="module")
def checker():
    c = InternalRBSChecker()
    c.initiate()
    return c


# ---------------------------------------------------------------------------
# Sequences that SHOULD pass
# ---------------------------------------------------------------------------

def test_normal_cds_passes(checker):
    """A typical CDS with no internal SD+ATG should pass."""
    seq = "ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTCTTGCATTTATCGCAGCATAA"
    ok, msg = checker.run(seq)
    assert ok
    assert msg is None

def test_sd_at_start_not_flagged(checker):
    """
    The legitimate SD sequence upstream is not part of the CDS.
    An ATG at the very start should not be flagged as internal.
    """
    # Starts with ATG (position 0) — this is the intended start, not internal
    seq = "ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTCTTGCATAA"
    ok, _ = checker.run(seq)
    assert ok

def test_agg_far_from_atg_passes(checker):
    """An AGG motif with no ATG within spacer range should pass."""
    # AGG at position 15, next ATG is 30+ nt away
    seq = "ATGCTTGACAAGTCGAGGCCTATGGAACGTAAGCTTTGACCTAAGCTTGACTAA"
    ok, _ = checker.run(seq)
    assert ok

def test_empty_sequence_passes(checker):
    """Empty sequence has no internal RBS."""
    ok, _ = checker.run("")
    assert ok

def test_mkkkkkkk_passes(checker):
    """
    MKKKKKKK uses AAA for K — lots of A's but no SD+ATG combination.
    This is the professor's example edge case peptide.
    """
    # Best codon for K is AAA — sequence is ATGAAAAAAAAAAAAAAAAAAAAATAA
    seq = "ATGAAAAAAAAAAAAAAAAAAAAATAA"
    ok, _ = checker.run(seq)
    assert ok


# ---------------------------------------------------------------------------
# Sequences that SHOULD fail
# ---------------------------------------------------------------------------

def test_perfect_sd_with_atg_fails(checker):
    """Perfect AGGAGG followed by ATG within spacer should be flagged."""
    # AGGAGG at position 20, ATG 8 nt later
    seq = "ATGCTTGACAAGTCGCCTATGAGGAGGTTTCCCATGCTTGACTAA"
    ok, msg = checker.run(seq)
    assert not ok
    assert "Internal RBS" in msg

def test_partial_sd_aaggag_fails(checker):
    """AAGGAG variant followed by ATG should be flagged."""
    seq = "ATGCTTGACAAGTCGCCTATGAAGGAGTTTCCCATGCTTGACTAA"
    ok, msg = checker.run(seq)
    assert not ok

def test_sd_with_atg_at_min_spacer_fails(checker):
    """SD motif with ATG exactly at minimum spacer distance."""
    # AGGAGG then 5 nt then ATG
    seq = "ATGCTTGACAAGTCGCCTATGAGGAGGAAAAATGCTTGACTAA"
    ok, msg = checker.run(seq)
    assert not ok

def test_message_contains_position_info(checker):
    """Failure message should include position information."""
    seq = "ATGCTTGACAAGTCGCCTATGAGGAGGTTTCCCATGCTTGACTAA"
    ok, msg = checker.run(seq)
    assert not ok
    assert "position" in msg.lower()


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

def test_short_sequence_passes(checker):
    """Sequence shorter than search window cannot have internal RBS."""
    seq = "ATGCTTGAC"
    ok, _ = checker.run(seq)
    assert ok

def test_agg_motif_no_downstream_atg_passes(checker):
    """AGG present but no ATG within short spacer — should pass."""
    # AGG at position 12, no ATG within 8 nt
    seq = "ATGCTTGACAAGTCGAGGCCCCCCCCCCCCCCCCCCCCCTAA"
    ok, _ = checker.run(seq)
    assert ok
