"""
test_transcript_designer.py

Edge case tests for TranscriptDesigner.

As described in the BioE134/234 lecture, TranscriptDesigner must handle
specific challenging inputs that are prone to causing incorrect behavior.
The MKKKKKKK peptide is a reliable challenge because Lysine (K) is encoded
by AAA/AAG — repeated K residues tend to produce poly(A) runs which are
forbidden sequences.
"""

import pytest
from genedesign.transcript_designer import TranscriptDesigner
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker
from genedesign.seq_utils.translate import Translate


@pytest.fixture(scope="module")
def designer():
    d = TranscriptDesigner()
    d.initiate()
    return d

@pytest.fixture(scope="module")
def forbidden_checker():
    c = ForbiddenSequenceChecker()
    c.initiate()
    return c

@pytest.fixture(scope="module")
def promoter_checker():
    c = PromoterChecker()
    c.initiate()
    return c

@pytest.fixture(scope="module")
def rbs_checker():
    c = InternalRBSChecker()
    c.initiate()
    return c

@pytest.fixture(scope="module")
def translator():
    t = Translate()
    t.initiate()
    return t


def get_cds(transcript):
    return ''.join(transcript.codons)


# ---------------------------------------------------------------------------
# MKKKKKKK — the professor's canonical edge case
# ---------------------------------------------------------------------------

def test_mkkkkkkk_no_polya(designer, forbidden_checker):
    """
    MKKKKKKK is prone to poly(A) because AAA is the best codon for K.
    The designer must avoid producing 8+ consecutive A's.
    """
    t = designer.run("MKKKKKKK", set())
    cds = get_cds(t)
    ok, site = forbidden_checker.run(cds)
    assert ok, f"Forbidden sequence found: {site}"

def test_mkkkkkkk_correct_translation(designer, translator):
    """CDS must translate back to the original protein."""
    peptide = "MKKKKKKK"
    t = designer.run(peptide, set())
    cds = get_cds(t)
    translated = translator.run(cds)
    assert translated == peptide

def test_mkkkkkkk_starts_with_atg(designer):
    """CDS must start with ATG."""
    t = designer.run("MKKKKKKK", set())
    assert get_cds(t).startswith("ATG")

def test_mkkkkkkk_ends_with_stop(designer):
    """CDS must end with a valid stop codon."""
    t = designer.run("MKKKKKKK", set())
    assert get_cds(t).endswith(("TAA", "TGA", "TAG"))

def test_mkkkkkkk_has_rbs(designer):
    """Transcript must have an RBS assigned."""
    t = designer.run("MKKKKKKK", set())
    assert t.rbs is not None
    assert t.rbs.utr is not None

def test_mkkkkkkk_no_internal_promoter(designer, promoter_checker):
    """CDS should not contain an accidental internal promoter."""
    t = designer.run("MKKKKKKK", set())
    ok, _ = promoter_checker.run(t.rbs.utr.upper() + get_cds(t))
    assert ok

def test_mkkkkkkk_no_internal_rbs(designer, rbs_checker):
    """CDS should not contain an accidental internal RBS."""
    t = designer.run("MKKKKKKK", set())
    ok, _ = rbs_checker.run(get_cds(t))
    assert ok

# ---------------------------------------------------------------------------
# Additional edge cases
# ---------------------------------------------------------------------------

def test_single_amino_acid(designer, translator):
    """Single amino acid protein — minimal valid CDS."""
    t = designer.run("M", set())
    cds = get_cds(t)
    assert cds.startswith("ATG")
    assert cds.endswith(("TAA", "TGA", "TAG"))

def test_ignores_respected(designer):
    """When an RBS is in the ignore set, a different one should be chosen."""
    t1 = designer.run("MYPFIRTARMTV", set())
    ignores = {t1.rbs}
    t2 = designer.run("MYPFIRTARMTV", ignores)
    assert t2.rbs not in ignores

def test_standard_protein_correct_translation(designer, translator):
    """Standard protein must translate back correctly."""
    peptide = "MYPFIRTARMTVCAKKHVHLTRDAAEQLLADID"
    t = designer.run(peptide, set())
    assert translator.run(get_cds(t)) == peptide

def test_no_forbidden_sequences_standard(designer, forbidden_checker):
    """Standard protein CDS should have no forbidden sequences."""
    t = designer.run("MYPFIRTARMTVCAKKHVHLTRDAAEQLLADID", set())
    ok, site = forbidden_checker.run(t.rbs.utr.upper() + get_cds(t))
    assert ok, f"Forbidden sequence found: {site}"
