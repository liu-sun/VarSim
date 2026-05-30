"""Internal module: shared codon-change classification logic.

Used by both ``missense()`` and ``codon_sub()`` to determine the HGVS
c. representation when a codon is replaced by another codon.
"""


def _classify_codon_change(orig_codon: str, new_codon: str, codon_start: int) -> str:
    """Return the HGVS c. suffix for a codon-to-codon substitution.

    Parameters
    ----------
    orig_codon : str
        The original 3-base codon (e.g. "ATG").
    new_codon : str
        The replacement 3-base codon (e.g. "CTG").
    codon_start : int
        0-indexed position of the codon start within the CDS.

    Returns
    -------
    str
        HGVS c. suffix like ``c.1A>G`` or ``c.2_3delinsTG``.
        The caller must prepend the reference accession (e.g. ``NM_000207.3:``).
    """
    c1, c2, c3 = codon_start + 1, codon_start + 2, codon_start + 3

    # Case 1: all three bases differ
    if orig_codon[0] != new_codon[0] and orig_codon[1] != new_codon[1] and orig_codon[2] != new_codon[2]:
        return f"c.{c1}_{c3}delins{new_codon}"

    # Case 2: base 1 same, bases 2+3 differ
    if orig_codon[0] == new_codon[0] and orig_codon[1] != new_codon[1] and orig_codon[2] != new_codon[2]:
        return f"c.{c2}_{c3}delins{new_codon[1:3]}"

    # Case 3: base 1 differs, base 2 same, base 3 differs
    if orig_codon[0] != new_codon[0] and orig_codon[1] == new_codon[1] and orig_codon[2] != new_codon[2]:
        return f"c.{c1}_{c3}delins{new_codon}"

    # Case 4: bases 1+2 differ, base 3 same
    if orig_codon[0] != new_codon[0] and orig_codon[1] != new_codon[1] and orig_codon[2] == new_codon[2]:
        return f"c.{c1}_{c2}delins{new_codon[0:2]}"

    # Case 5: bases 1+2 same, base 3 differs (single substitution at pos 3)
    if orig_codon[0] == new_codon[0] and orig_codon[1] == new_codon[1] and orig_codon[2] != new_codon[2]:
        return f"c.{c3}{orig_codon[2]}>{new_codon[2]}"

    # Case 6: bases 1+3 same, base 2 differs (single substitution at pos 2)
    if orig_codon[0] == new_codon[0] and orig_codon[1] != new_codon[1] and orig_codon[2] == new_codon[2]:
        return f"c.{c2}{orig_codon[1]}>{new_codon[1]}"

    # Case 7: base 1 differs, bases 2+3 same (single substitution at pos 1)
    # (orig_codon[0] != new_codon[0] and orig_codon[1] == new_codon[1] and orig_codon[2] == new_codon[2])
    return f"c.{c1}{orig_codon[0]}>{new_codon[0]}"
