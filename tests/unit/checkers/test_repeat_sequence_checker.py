"""
test_repeat_sequence_checker.py

Unit tests for RepeatSequenceChecker.

Run with:
    pytest tests/unit/checkers/test_repeat_sequence_checker.py -v
"""

import pytest
from genedesign.checkers.repeat_sequence_checker import RepeatSequenceChecker


@pytest.fixture(scope="module")
def checker():
    c = RepeatSequenceChecker()
    c.initiate()
    return c


# ---------------------------------------------------------------------------
# Sequences that SHOULD pass (no tandem repeats)
# ---------------------------------------------------------------------------

def test_normal_cds_passes(checker):
    """A realistic codon-optimized CDS with no repeats should pass."""
    seq = "ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTCTTGCATTTATCGCAGCATAA"
    ok, msg = checker.run(seq)
    assert ok
    assert msg is None

def test_short_sequence_passes(checker):
    """A sequence shorter than 2x MIN_UNIT_LEN cannot have a tandem repeat."""
    seq = "ATGCTTGA"
    ok, _ = checker.run(seq)
    assert ok

def test_empty_sequence_passes(checker):
    """Empty string has no repeats."""
    ok, _ = checker.run("")
    assert ok

def test_no_repeat_diverse_sequence(checker):
    """A sequence with high codon diversity and no repeated blocks should pass."""
    seq = "ATGCTTGACAAGTCGCCTATGGAACGTAAGCTTTGACCTAAGCTTGACTAA"
    ok, _ = checker.run(seq)
    assert ok


# ---------------------------------------------------------------------------
# Sequences that SHOULD fail (tandem repeats present)
# ---------------------------------------------------------------------------

def test_exact_tandem_repeat_fails(checker):
    """Two identical 10-bp units back-to-back should be flagged."""
    unit = "ATGCTTGACA"  # 10 bp
    seq = "ATG" + unit + unit + "TAA"
    ok, msg = checker.run(seq)
    assert not ok
    assert "Tandem repeat" in msg

def test_tandem_repeat_message_contains_unit(checker):
    """The failure message should include the repeated unit."""
    unit = "ATGCTTGACAAG"  # 12 bp
    seq = unit + unit + "CGTTAA"
    ok, msg = checker.run(seq)
    assert not ok
    assert unit in msg

def test_long_tandem_repeat_fails(checker):
    """A longer repeated unit (18 bp) should also be caught."""
    unit = "ATGCTTGACAAGTCGCCT"  # 18 bp
    seq = "ATG" + unit + unit + "TAA"
    ok, msg = checker.run(seq)
    assert not ok

def test_all_same_base_fails(checker):
    """40 identical bases trivially contain tandem repeats."""
    seq = "A" * 40
    ok, _ = checker.run(seq)
    assert not ok


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

def test_repeat_unit_just_below_threshold_passes(checker):
    """
    A repeated unit of length MIN_UNIT_LEN - 1 (7 bp) should pass
    because we only flag repeats of length >= MIN_UNIT_LEN (8).
    """
    unit = "ATGCTTG"  # 7 bp — just below threshold
    seq = "ATG" + unit + unit + "CGTTAA"
    ok, _ = checker.run(seq)
    assert ok

def test_repeat_at_end_of_sequence_fails(checker):
    """Tandem repeat occurring at the very end of the sequence should be caught."""
    unit = "ATGCTTGACAAG"  # 12 bp
    seq = "GCGATCGTACCG" + unit + unit
    ok, _ = checker.run(seq)
    assert not ok
