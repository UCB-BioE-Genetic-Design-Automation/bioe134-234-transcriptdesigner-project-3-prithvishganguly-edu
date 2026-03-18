"""
repeat_sequence_checker.py

Checks a DNA coding sequence for tandem direct repeats.

Biological motivation
---------------------
Tandem repeats — where the same sequence appears back-to-back — are hotspots
for homologous recombination in E. coli. The cell's recombination machinery
can loop out and delete everything between two identical sequences, destroying
your construct. This is a real, common failure mode in synthetic biology.

This checker scans for any repeated unit of length >= MIN_UNIT_LEN that
appears consecutively (e.g. ATGCATGC where ATGC is the unit).
"""

class RepeatSequenceChecker:
    """
    Checks a DNA sequence for tandem direct repeats.

    Usage:
        checker = RepeatSequenceChecker()
        checker.initiate()
        passed, message = checker.run("ATGCGT...")
    """

    # Minimum repeat unit length to flag (shorter repeats are less dangerous)
    MIN_UNIT_LEN = 8

    def initiate(self) -> None:
        """No external data needed."""
        pass

    def run(self, cds: str) -> tuple[bool, str | None]:
        """
        Scans the CDS for tandem direct repeats.

        Parameters:
            cds (str): DNA coding sequence.

        Returns:
            (True, None)    if no problematic tandem repeats found.
            (False, msg)    where msg names the repeated unit found.
        """
        seq = cds.upper()
        n = len(seq)

        # Slide a window of increasing size across the sequence.
        # For each position i and unit length k, check whether the k-mer
        # starting at i is identical to the k-mer immediately after it.
        for k in range(self.MIN_UNIT_LEN, n // 2 + 1):
            for i in range(n - 2 * k + 1):
                unit = seq[i:i + k]
                next_unit = seq[i + k:i + 2 * k]
                if unit == next_unit:
                    return (
                        False,
                        f"Tandem repeat found: '{unit}' repeated at positions {i} and {i + k}"
                    )

        return True, None
