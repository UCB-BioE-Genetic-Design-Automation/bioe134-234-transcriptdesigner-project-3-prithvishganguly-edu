"""
transcript_designer.py

Improved TranscriptDesigner using a fast sliding window algorithm.

Speed strategy
--------------
The sliding window only runs CHEAP checkers (forbidden sequences and
internal RBS) during window enumeration. These are simple string searches
and run in microseconds.

The EXPENSIVE checkers (hairpin, promoter, codon) are only run once on
the final complete sequence, with a small number of weighted random
fallback retries if needed.

This makes the sliding window practical on a full proteome.

Algorithm
---------
1. Left to right, 1 codon at a time, 3aa lookahead window
2. Only top 2 codons per AA considered (keeps combinations at 2^3 = 8 max)
3. Score each combo with fast checkers only
4. Lock in best codon for current position, advance
5. Final full check with all 5 checkers
6. Weighted random fallback retries if final check fails

Checkers
--------
Window scoring (fast):
  - ForbiddenSequenceChecker
  - InternalRBSChecker

Final check (all 5):
  - CodonChecker
  - ForbiddenSequenceChecker
  - hairpin_checker
  - PromoterChecker
  - InternalRBSChecker
"""

import random
from itertools import product

from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker


class TranscriptDesigner:

    WINDOW_SIZE = 3
    TOP_N_CODONS = 2
    FALLBACK_ATTEMPTS = 20

    def __init__(self):
        self.rbsChooser = None
        self.codonChecker = None
        self.forbiddenChecker = None
        self.promoterChecker = None
        self.internalRBSChecker = None
        self.bestCodon = {}
        self.allCodons = {}
        self.topCodons = {}
        self.codonWeights = {}
        self.altStartCodons = {'V': 'GTG', 'L': 'TTG'}

    def initiate(self) -> None:

        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()

        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.forbiddenChecker.initiate()

        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()

        self.internalRBSChecker = InternalRBSChecker()
        self.internalRBSChecker.initiate()

        self.bestCodon = {
            'A': 'GCG', 'C': 'TGC', 'D': 'GAT', 'E': 'GAA', 'F': 'TTC',
            'G': 'GGT', 'H': 'CAC', 'I': 'ATC', 'K': 'AAA', 'L': 'CTG',
            'M': 'ATG', 'N': 'AAC', 'P': 'CCG', 'Q': 'CAG', 'R': 'CGT',
            'S': 'TCT', 'T': 'ACC', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAC'
        }

        codon_freq = {
            'A': [('GCG', 0.33), ('GCC', 0.26), ('GCA', 0.23), ('GCT', 0.18)],
            'C': [('TGC', 0.54), ('TGT', 0.46)],
            'D': [('GAT', 0.63), ('GAC', 0.37)],
            'E': [('GAA', 0.68), ('GAG', 0.32)],
            'F': [('TTC', 0.58), ('TTT', 0.42)],
            'G': [('GGT', 0.38), ('GGC', 0.37), ('GGG', 0.15), ('GGA', 0.09)],
            'H': [('CAC', 0.57), ('CAT', 0.43)],
            'I': [('ATC', 0.42), ('ATT', 0.28), ('ATA', 0.08)],
            'K': [('AAA', 0.74), ('AAG', 0.26)],
            'L': [('CTG', 0.54), ('TTA', 0.11), ('TTG', 0.11),
                  ('CTT', 0.10), ('CTC', 0.10), ('CTA', 0.04)],
            'M': [('ATG', 1.0)],
            'N': [('AAC', 0.51), ('AAT', 0.49)],
            'P': [('CCG', 0.55), ('CCA', 0.20), ('CCT', 0.16), ('CCC', 0.10)],
            'Q': [('CAG', 0.69), ('CAA', 0.31)],
            'R': [('CGT', 0.42), ('CGC', 0.37), ('CGA', 0.06),
                  ('CGG', 0.07), ('AGA', 0.04), ('AGG', 0.04)],
            'S': [('AGC', 0.24), ('TCT', 0.17), ('TCC', 0.15),
                  ('AGT', 0.16), ('TCA', 0.14), ('TCG', 0.15)],
            'T': [('ACC', 0.40), ('ACG', 0.28), ('ACT', 0.17), ('ACA', 0.15)],
            'V': [('GTG', 0.35), ('GTT', 0.28), ('GTC', 0.20), ('GTA', 0.17)],
            'W': [('TGG', 1.0)],
            'Y': [('TAC', 0.57), ('TAT', 0.43)],
        }

        self.allCodons = {aa: [c for c, _ in pairs] for aa, pairs in codon_freq.items()}
        self.codonWeights = {aa: [w for _, w in pairs] for aa, pairs in codon_freq.items()}
        self.topCodons = {
            aa: [c for c, _ in pairs[:self.TOP_N_CODONS]]
            for aa, pairs in codon_freq.items()
        }

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _get_start_codon(self, first_aa: str) -> str:
        if first_aa == 'M':
            return 'ATG'
        elif first_aa in self.altStartCodons:
            return self.altStartCodons[first_aa]
        else:
            return self.bestCodon[first_aa]

    def _weighted_codon(self, aa: str) -> str:
        return random.choices(
            self.allCodons[aa],
            weights=self.codonWeights[aa],
            k=1
        )[0]

    def _fast_score(self, codons: list, utr: str) -> int:
        """
        Fast scoring using only cheap string-search checkers.
        Used during sliding window enumeration.
        Returns 0-2.
        """
        cds = ''.join(codons)
        transcript_dna = utr + cds
        forbidden_ok, _ = self.forbiddenChecker.run(transcript_dna)
        rbs_ok, _       = self.internalRBSChecker.run(cds)
        return sum([forbidden_ok, rbs_ok])

    def _full_score(self, codons: list, utr: str) -> int:
        """
        Full scoring using all 5 checkers.
        Used only on final sequence and fallback retries.
        Returns 0-5.
        """
        cds = ''.join(codons)
        transcript_dna = utr + cds
        codon_ok, _, _, _ = self.codonChecker.run(codons)
        forbidden_ok, _   = self.forbiddenChecker.run(transcript_dna)
        hairpin_ok, _     = hairpin_checker(transcript_dna)
        promoter_ok, _    = self.promoterChecker.run(transcript_dna)
        rbs_ok, _         = self.internalRBSChecker.run(cds)
        return sum([codon_ok, forbidden_ok, hairpin_ok, promoter_ok, rbs_ok])

    # ------------------------------------------------------------------
    # Sliding window
    # ------------------------------------------------------------------

    def _sliding_window_design(self, peptide: str, utr: str) -> list:
        """
        Fast sliding window: only cheap checkers used during enumeration.
        """
        chosen = []

        for i in range(len(peptide)):
            aa = peptide[i]

            if i == 0:
                chosen.append(self._get_start_codon(aa))
                continue

            window_aas = peptide[i:i + self.WINDOW_SIZE]
            options_per_pos = [self.topCodons[a] for a in window_aas]

            best_score = -1
            best_codon_for_i = self.topCodons[aa][0]

            for combo in product(*options_per_pos):
                candidate = chosen + list(combo)
                score = self._fast_score(candidate, utr)
                if score > best_score:
                    best_score = score
                    best_codon_for_i = combo[0]

            chosen.append(best_codon_for_i)

        return chosen

    # ------------------------------------------------------------------
    # Main run
    # ------------------------------------------------------------------

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Designs a Transcript for the given peptide sequence.

        Parameters
        ----------
        peptide : str
            Amino acid sequence (single-letter codes, no stop symbol).
        ignores : set
            RBS options to exclude.

        Returns
        -------
        Transcript
        """
        peptide = peptide.strip().rstrip('*').upper()

        # Get UTR for context scoring
        placeholder_cds = self._get_start_codon(peptide[0]) + \
                          ''.join(self.bestCodon[aa] for aa in peptide[1:]) + 'TAA'
        candidate_rbs = self.rbsChooser.run(placeholder_cds, ignores)
        utr = candidate_rbs.utr.upper()

        # Phase 1: fast sliding window
        codons = self._sliding_window_design(peptide, utr)
        codons.append('TAA')

        # Phase 2: full score on complete sequence
        cds = ''.join(codons)
        final_rbs = self.rbsChooser.run(cds, ignores)
        final_utr = final_rbs.utr.upper()
        best_score = self._full_score(codons[:-1], final_utr)
        best_codons = codons

        # Phase 3: weighted random fallback if full check not perfect
        if best_score < 5:
            for _ in range(self.FALLBACK_ATTEMPTS):
                fallback = [self._get_start_codon(peptide[0])]
                fallback += [self._weighted_codon(aa) for aa in peptide[1:]]
                fallback.append('TAA')

                fb_cds = ''.join(fallback)
                fb_rbs = self.rbsChooser.run(fb_cds, ignores)
                fb_score = self._full_score(fallback[:-1], fb_rbs.utr.upper())

                if fb_score > best_score:
                    best_score = fb_score
                    best_codons = fallback

                if best_score == 5:
                    break

        final_cds = ''.join(best_codons)
        selected_rbs = self.rbsChooser.run(final_cds, ignores)
        return Transcript(selected_rbs, peptide, best_codons)


if __name__ == "__main__":
    designer = TranscriptDesigner()
    designer.initiate()
    ignores = set()

    t1 = designer.run("MYPFIRTARMTV", ignores)
    print("Standard:", t1.codons)

    t2 = designer.run("MKKKKKKK", ignores)
    print("MKKKKKKK:", t2.codons)
    print("Has poly-A?", "AAAAAAAA" in ''.join(t2.codons))
