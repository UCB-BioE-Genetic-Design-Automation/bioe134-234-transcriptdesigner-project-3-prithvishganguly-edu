"""
internal_rbs_checker.py

Checks for internal ribosome binding sites (Shine-Dalgarno sequences)
within a coding sequence.

Biological motivation
---------------------
In E. coli, translation initiation requires a Shine-Dalgarno (SD) sequence
upstream of the start codon. The canonical SD sequence is AGGAGG, but
partial matches like AGG, AAGG, GAGG also work.

If your designed CDS accidentally contains an SD-like sequence followed
by an ATG within ~10 bp, the ribosome may initiate translation internally,
producing a truncated, wrong protein from within your gene. This is a real
failure mode in synthetic biology constructs.

This checker scans the CDS for any SD-like sequence (partial matches to
AGGAGG) that is followed by an ATG within a 5-15 bp window, flagging
potential internal RBS sites.
"""


class InternalRBSChecker:
    """
    Checks a DNA coding sequence for internal ribosome binding sites.

    An internal RBS is defined as:
        - An SD-like sequence (matches one of the known motifs)
        - Followed by ATG within 5-15 nucleotides downstream

    Usage
    -----
    checker = InternalRBSChecker()
    checker.initiate()
    passed, message = checker.run("ATGCGT...")
    """

    # Shine-Dalgarno motifs — ranked from strongest to weakest match
    # The ribosome recognises these as translation start signals
    SD_MOTIFS = [
        "AGGAGG",   # perfect SD
        "AAGGAG",   # common variant
        "AGGAG",    # 5-mer core
        "GAGG",     # partial
        "AAGG",     # partial
        "AGGA",     # partial
        "AGG",      # minimal — only flag if followed closely by ATG
    ]

    # Window after SD motif in which ATG is considered "in range"
    MIN_SPACER = 5
    MAX_SPACER = 15

    # For short motifs (3-4 bp), require a closer ATG to reduce false positives
    SHORT_MOTIF_MAX_SPACER = 8

    def initiate(self) -> None:
        """No external data needed."""
        pass

    def run(self, cds: str) -> tuple[bool, str | None]:
        """
        Scans the CDS for internal ribosome binding sites.

        Skips the first 10 nucleotides (the legitimate start region)
        to avoid flagging the intended start codon.

        Parameters
        ----------
        cds : str
            DNA coding sequence (ATG...stop).

        Returns
        -------
        (True, None)      if no internal RBS found.
        (False, message)  describing the problematic site.
        """
        seq = cds.upper()

        # Skip the first 10 nt — that's where the legitimate start codon lives
        search_start = 10

        for motif in self.SD_MOTIFS:
            motif_len = len(motif)
            max_spacer = self.MAX_SPACER if motif_len >= 5 else self.SHORT_MOTIF_MAX_SPACER

            pos = search_start
            while True:
                idx = seq.find(motif, pos)
                if idx == -1:
                    break

                # Look for ATG within spacer window after the motif
                atg_search_start = idx + motif_len + self.MIN_SPACER
                atg_search_end = idx + motif_len + max_spacer + 3

                region = seq[atg_search_start:atg_search_end]
                atg_pos = region.find("ATG")

                if atg_pos != -1:
                    absolute_atg = atg_search_start + atg_pos
                    spacer = atg_search_start + atg_pos - (idx + motif_len)
                    return (
                        False,
                        f"Internal RBS detected: SD-like motif '{motif}' at position "
                        f"{idx} followed by ATG at position {absolute_atg} "
                        f"(spacer: {spacer} nt)"
                    )

                pos = idx + 1

        return True, None


if __name__ == "__main__":
    checker = InternalRBSChecker()
    checker.initiate()

    # Should fail — has AGGAGG followed by ATG within spacer
    bad = "ATGCTTGACAAGTCGCCTATGGAACGTAGGAGGTTTCCCATGCTTGAC"
    ok, msg = checker.run(bad)
    print(f"Bad sequence: {ok}, {msg}")

    # Should pass — no internal SD+ATG combination
    good = "ATGCTTGACAAGTCGCCTATGGAACGTAAGCTTTGACCTAAGCTTGAC"
    ok, msg = checker.run(good)
    print(f"Good sequence: {ok}, {msg}")
