import random
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker

class TranscriptDesigner:
    """
    Improved TranscriptDesigner that reverse translates a protein sequence into
    a codon-optimized DNA sequence satisfying four biological quality checks:

        1. CodonChecker             - high CAI, low rare codons, good diversity
        2. ForbiddenSequenceChecker - no restriction sites or repeat motifs
        3. hairpin_checker          - no mRNA hairpins blocking translation
        4. PromoterChecker          - no accidental internal sigma70 promoters

    Strategy: simulated annealing retry.
        - Attempt 0: greedy (best codon per amino acid).
        - Attempts 1-N: increasingly random codon sampling.
        - Each attempt checks RBS UTR + CDS together, matching the benchmarker.
        - Returns the first design passing all 4 checks, or best seen if exhausted.

    Start codon handling:
        - If protein starts with M: use ATG (standard).
        - If protein starts with anything else: use the best codon for that AA,
          which the benchmarker accepts as long as it translates back correctly.
          The CDS will still start with a valid codon for that amino acid.
          Note: GTG and TTG are alternative start codons in E. coli but encode
          Val and Leu respectively when translated — we handle this by keeping
          the first amino acid codon faithful to the input sequence.
    """

    MAX_ATTEMPTS = 50

    def __init__(self):
        self.rbsChooser = None
        self.codonChecker = None
        self.forbiddenChecker = None
        self.promoterChecker = None
        self.bestCodon = {}
        self.allCodons = {}

        # Valid start codons in E. coli mapped to the amino acid they encode
        # when used as a start codon (vs internal codon)
        # ATG → M, GTG → V (but acts as start), TTG → L (but acts as start)
        self.altStartCodons = {
            'V': 'GTG',  # Valine — valid alternative start in E. coli
            'L': 'TTG',  # Leucine — valid alternative start in E. coli
        }

    def initiate(self) -> None:
        """
        Initializes all sub-components and codon lookup tables.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()

        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.forbiddenChecker.initiate()

        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()

        # Best single codon per amino acid (highest CAI for E. coli).
        # Used as the greedy baseline on attempt 0.
        self.bestCodon = {
            'A': 'GCG', 'C': 'TGC', 'D': 'GAT', 'E': 'GAA', 'F': 'TTC',
            'G': 'GGT', 'H': 'CAC', 'I': 'ATC', 'K': 'AAA', 'L': 'CTG',
            'M': 'ATG', 'N': 'AAC', 'P': 'CCG', 'Q': 'CAG', 'R': 'CGT',
            'S': 'TCT', 'T': 'ACC', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAC'
        }

        # All synonymous codons per amino acid.
        # Every codon encodes the correct amino acid — verified against standard
        # genetic code. Used for randomization during retry attempts.
        self.allCodons = {
            'A': ['GCT', 'GCC', 'GCA', 'GCG'],
            'C': ['TGT', 'TGC'],
            'D': ['GAT', 'GAC'],
            'E': ['GAA', 'GAG'],
            'F': ['TTT', 'TTC'],
            'G': ['GGT', 'GGC', 'GGA', 'GGG'],
            'H': ['CAT', 'CAC'],
            'I': ['ATT', 'ATC', 'ATA'],
            'K': ['AAA', 'AAG'],
            'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
            'M': ['ATG'],
            'N': ['AAT', 'AAC'],
            'P': ['CCT', 'CCC', 'CCA', 'CCG'],
            'Q': ['CAA', 'CAG'],
            'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
            'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
            'T': ['ACT', 'ACC', 'ACA', 'ACG'],
            'V': ['GTT', 'GTC', 'GTA', 'GTG'],
            'W': ['TGG'],
            'Y': ['TAT', 'TAC'],
        }

    def _get_start_codon(self, first_aa: str) -> str:
        """
        Returns the correct start codon for the first amino acid.

        E. coli recognises ATG, GTG, and TTG as start codons.
        - M always uses ATG.
        - V can use GTG (preferred alternative start).
        - L can use TTG (preferred alternative start).
        - All other amino acids: use their standard best codon.
          The benchmarker checks translation correctness, so we must
          use a codon that actually encodes the right amino acid.
        """
        if first_aa == 'M':
            return 'ATG'
        elif first_aa in self.altStartCodons:
            return self.altStartCodons[first_aa]
        else:
            return self.bestCodon[first_aa]

    def _get_start_codon_random(self, first_aa: str) -> str:
        """
        Returns a random valid codon for the first amino acid during retries.
        Same logic as _get_start_codon but randomized from allCodons.
        """
        if first_aa == 'M':
            return 'ATG'
        elif first_aa in self.altStartCodons:
            return self.altStartCodons[first_aa]
        else:
            return random.choice(self.allCodons[first_aa])

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Designs a Transcript for the given peptide sequence.

        Parameters:
            peptide (str): Amino acid sequence (single-letter codes, no stop symbol).
            ignores (set): RBS options to exclude (passed to RBSChooser).

        Returns:
            Transcript: The best designed transcript found within MAX_ATTEMPTS.
        """
        # Clean the input
        peptide = peptide.strip().rstrip('*').upper()

        first_aa = peptide[0]
        rest = peptide[1:]

        best_codons = None
        best_pass_count = -1

        for attempt in range(self.MAX_ATTEMPTS):

            # ------------------------------------------------------------------
            # Step 1: Generate a candidate codon list
            # ------------------------------------------------------------------
            if attempt == 0:
                # Greedy: best codon for every position
                # First AA gets special start codon handling
                start_codon = self._get_start_codon(first_aa)
                codons = [start_codon] + [self.bestCodon[aa] for aa in rest]
            else:
                # Randomised retry with increasing randomness
                randomness = attempt / self.MAX_ATTEMPTS
                start_codon = self._get_start_codon(first_aa) if random.random() > randomness \
                              else self._get_start_codon_random(first_aa)
                codons = [start_codon]
                for aa in rest:
                    if random.random() < randomness:
                        codons.append(random.choice(self.allCodons[aa]))
                    else:
                        codons.append(self.bestCodon[aa])

            # Always append a stop codon (TAA preferred in E. coli)
            codons.append('TAA')

            # ------------------------------------------------------------------
            # Step 2: Build transcript DNA including RBS UTR
            # The benchmarker validates rbs.utr + cds, so we check the same.
            # ------------------------------------------------------------------
            cds = ''.join(codons)
            candidate_rbs = self.rbsChooser.run(cds, ignores)
            transcript_dna = candidate_rbs.utr.upper() + cds

            # ------------------------------------------------------------------
            # Step 3: Run all 4 checkers with correct input types
            #   - CodonChecker  → takes codons LIST
            #   - Others        → take transcript DNA STRING
            # ------------------------------------------------------------------
            codon_ok, _, _, _ = self.codonChecker.run(codons)
            forbidden_ok, _   = self.forbiddenChecker.run(transcript_dna)
            hairpin_ok, _     = hairpin_checker(transcript_dna)
            promoter_ok, _    = self.promoterChecker.run(transcript_dna)

            pass_count = sum([codon_ok, forbidden_ok, hairpin_ok, promoter_ok])

            # ------------------------------------------------------------------
            # Step 4: Track best result seen
            # ------------------------------------------------------------------
            if pass_count > best_pass_count:
                best_pass_count = pass_count
                best_codons = codons

            # ------------------------------------------------------------------
            # Step 5: Stop if all 4 pass
            # ------------------------------------------------------------------
            if pass_count == 4:
                break

        # ----------------------------------------------------------------------
        # Step 6: Return the best Transcript found
        # ----------------------------------------------------------------------
        final_cds = ''.join(best_codons)
        selected_rbs = self.rbsChooser.run(final_cds, ignores)
        return Transcript(selected_rbs, peptide, best_codons)


if __name__ == "__main__":
    # Test with a standard M-start protein
    peptide = "MYPFIRTARMTV"
    designer = TranscriptDesigner()
    designer.initiate()
    ignores = set()
    transcript = designer.run(peptide, ignores)
    print(transcript)

    # Test with a non-M start protein (like many in the proteome)
    peptide2 = "GYNVTMRDIK"
    transcript2 = designer.run(peptide2, ignores)
    print(transcript2)
