"""Protein-to-nucleotide backtranslation.

Given a protein HGVS description, compute all possible single-nucleotide
variants (SNVs) that could produce the observed amino acid change.  Supports
both standalone protein HGVS and gene-anchored backtranslation that validates
against the actual MANE CDS sequence via NCBI.

Examples
--------
>>> backtranslate_protein("NP_000198.1:p.(V42G)")
['c.125T>G']

>>> backtranslate_protein("NP_000198.1:p.(V42=)")
['c.126A>C', 'c.126A>G', 'c.126A>T', 'c.126C>A', 'c.126C>G', 'c.126C>T', 'c.126G>A', 'c.126G>C', 'c.126G>T', 'c.126T>A', 'c.126T>C', 'c.126T>G']

>>> backtranslate_protein("NP_000198.1:p.(M1?)")
[]

>>> backtranslate_protein("NP_000198.1:p.(V42Gfs*17)")
[]
"""

import re

from Bio.Data.CodonTable import standard_dna_table

from . import _fetch
from ._logging import get_logger
from ._utils import genetic_code
from .parser import parse

logger = get_logger(__name__)

# Regex to extract (ref_aa)(pos)(alt_aa) from protein HGVS parentheses
_PROTEIN_SUB_RE = re.compile(
    r"\((?P<ref>[A-Z*])(?P<pos>\d+)(?P<alt>[A-Z*=?])\)"
)

# ---------------------------------------------------------------------------
# Amino acid → codon reverse mapping (includes stop codons)
# ---------------------------------------------------------------------------

_AA_TO_CODONS: dict[str, list[str]] = {}
for _codon in genetic_code:
    if _codon in standard_dna_table.forward_table:
        _aa = standard_dna_table.forward_table[_codon]
    else:
        _aa = "*"  # stop codon
    _AA_TO_CODONS.setdefault(_aa, []).append(_codon)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _hamming(c1: str, c2: str) -> int:
    """Return number of differing bases between two codons of equal length."""
    return sum(1 for a, b in zip(c1, c2) if a != b)


def _diff_position(ref_codon: str, alt_codon: str) -> tuple | None:
    """If *ref_codon* and *alt_codon* differ by exactly 1 base, return
    ``(0-based_offset, ref_base, alt_base)``; otherwise return ``None``.
    """
    diffs = [
        (i, ref_codon[i], alt_codon[i])
        for i in range(3)
        if ref_codon[i] != alt_codon[i]
    ]
    if len(diffs) != 1:
        return None
    return diffs[0]


def _codon_snv_to_c_hgvs(
    ref_codon: str, alt_codon: str, codon_start_0: int
) -> str | None:
    """Return a ``c.NNNX>Y`` suffix for a single-base codon change, or None."""
    diff = _diff_position(ref_codon, alt_codon)
    if diff is None:
        return None
    pos_in_codon, ref_base, alt_base = diff
    return f"c.{codon_start_0 + pos_in_codon + 1}{ref_base}>{alt_base}"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def _extract_protein_sub(hgvs_str: str) -> tuple[str, int, str] | None:
    """Extract ``(ref_aa, position, alt_aa)`` from a protein HGVS string.

    Returns ``None`` if the string does not contain a recognisable
    protein substitution notation like ``(V42G)`` or ``(V42=)``.
    """
    m = _PROTEIN_SUB_RE.search(hgvs_str)
    if not m:
        return None
    ref_aa = m.group("ref")
    position = int(m.group("pos"))
    alt_aa = m.group("alt")
    return ref_aa, position, alt_aa


def backtranslate_protein(hgvs_str: str) -> list[str]:
    """Compute all possible c.HGVS SNVs that could produce a protein change.

    Given a protein HGVS description (e.g. ``p.(V42G)``), enumerates every
    pair of reference and alternate codons that differ by exactly one
    nucleotide and returns the corresponding c.HGVS strings.

    Parameters
    ----------
    hgvs_str : str
        A protein HGVS description.  May include an accession prefix
        (``NP_000198.1:p.(V42G)``) or omit it (``p.(V42G)``).

    Returns
    -------
    list of str
        Deduplicated list of c.HGVS strings (e.g. ``["c.125T>G"]``).
        Returns an empty list for frameshifts, start-codon / uncertain
        variants, and unrecognised amino acids.

    Special cases
    -------------
    * **Silent** (``p.(V42=)``) — returns all synonymous SNVs at that
      position (changes that do not alter the amino acid).
    * **Start codon** (``p.(M1?)``) — returns ``[]`` (unpredictable effect).
    * **Frameshift** (``p.(V42Gfs*17)``) — returns ``[]`` (indels cause
      frameshifts that cannot be enumerated as simple SNVs).

    Examples
    --------
    >>> backtranslate_protein("NP_000198.1:p.(V42G)")
    ['c.125T>G']

    >>> backtranslate_protein("NP_000198.1:p.(V42=)")
    ['c.126A>C', 'c.126A>G', 'c.126A>T', 'c.126C>A', 'c.126C>G', 'c.126C>T', 'c.126G>A', 'c.126G>C', 'c.126G>T', 'c.126T>A', 'c.126T>C', 'c.126T>G']

    >>> backtranslate_protein("NP_000198.1:p.(M1?)")
    []

    >>> backtranslate_protein("NP_000198.1:p.(V42Gfs*17)")
    []

    >>> # Also accepts bare p.HGVS (no accession)
    >>> backtranslate_protein("p.(V42G)")
    ['c.125T>G']
    """
    logger.debug("Backtranslating: %s", hgvs_str)
    # --- Extract protein substitution notation ---
    sub = _extract_protein_sub(hgvs_str)
    if sub is None:
        # Fall back to parser for frameshift / uncertain classification
        if ":" not in hgvs_str:
            hgvs_str = f"NP_000000.0:{hgvs_str}"
        tag = parse(hgvs_str)
        if tag.variant_type == "frameshift" or tag.is_uncertain:
            return []
        return []

    ref_aa, position, alt_aa = sub

    # --- Reject unsupported cases ---
    if alt_aa == "?":
        return []  # uncertain / start-codon effect

    # Check for frameshift notation in the original string
    if "fs" in hgvs_str:
        return []

    # Silent variant: target AA = reference AA
    if alt_aa == "=":
        alt_aa = ref_aa

    ref_codons = _AA_TO_CODONS.get(ref_aa, [])
    alt_codons = _AA_TO_CODONS.get(alt_aa, [])

    if not ref_codons or not alt_codons:
        return []

    codon_start_0: int = (position - 1) * 3   # 0-based CDS offset
    seen: set[str] = set()

    for rc in ref_codons:
        for ac in alt_codons:
            if rc == ac:
                continue
            hgvs_suffix = _codon_snv_to_c_hgvs(rc, ac, codon_start_0)
            if hgvs_suffix and hgvs_suffix not in seen:
                seen.add(hgvs_suffix)

    result = sorted(seen)
    logger.debug("Found %d candidate SNVs", len(result))
    return result


def backtranslate(gene: str, p_hgvs: str) -> list[str]:
    """Backtranslate a protein HGVS for a specific gene using the real CDS.

    Fetches the MANE Select nucleotide record for *gene* via NCBI Entrez,
    extracts the CDS, validates that the reference amino acid at the given
    position matches the expected value, then generates all single-nucleotide
    variants at that codon that would produce the target amino acid.

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. ``"G6PD"``, ``"INS"``).
    p_hgvs : str
        Protein HGVS description, with or without accession prefix
        (e.g. ``"p.(V42G)"`` or ``"NP_000198.1:p.(V42G)"``).

    Returns
    -------
    list of str
        Deduplicated list of full c.HGVS strings with the NM accession
        (e.g. ``["NM_001360016.2:c.125T>G"]``).

    Raises
    ------
    ValueError
        If *p_hgvs* is a frameshift, uncertain, or start-codon variant,
        or if the reference amino acid in the HGVS does not match the
        actual CDS sequence at that position, or if the position exceeds
        the CDS length.

    Examples
    --------
    >>> backtranslate("G6PD", "p.(V42G)")  # doctest: +SKIP
    ['NM_001360016.2:c.125T>G']

    >>> backtranslate("G6PD", "p.(M1?)")  # doctest: +SKIP
    Traceback (most recent call last):
        ...
    ValueError: Uncertain protein variants cannot be backtranslated: ...

    >>> backtranslate("G6PD", "p.(V42Gfs*17)")  # doctest: +SKIP
    Traceback (most recent call last):
        ...
    ValueError: Frameshift variants cannot be backtranslated: ...
    """
    logger.info("Backtranslating %s for %s ...", p_hgvs, gene)
    # --- Extract protein substitution notation ---
    sub = _extract_protein_sub(p_hgvs)
    if sub is None:
        # Check for frameshift / uncertain via parser
        if ":" not in p_hgvs:
            p_hgvs = f"NP_000000.0:{p_hgvs}"
        tag = parse(p_hgvs)
        if tag.variant_type == "frameshift":
            raise ValueError(
                f"Frameshift variants cannot be backtranslated: {tag.original}"
            )
        if tag.is_uncertain:
            raise ValueError(
                f"Uncertain protein variants cannot be backtranslated: {tag.original}"
            )
        raise ValueError(
            f"Cannot extract protein substitution from: {p_hgvs}"
        )

    ref_aa, position, alt_aa = sub

    # --- Reject unsupported cases ---
    if alt_aa == "?":
        raise ValueError(
            f"Uncertain protein variants cannot be backtranslated: {p_hgvs}"
        )

    if "fs" in p_hgvs:
        raise ValueError(
            f"Frameshift variants cannot be backtranslated: {p_hgvs}"
        )

    # Silent variant: target AA = reference AA
    if alt_aa == "=":
        alt_aa = ref_aa

    # --- Fetch the actual CDS ---
    seqrecord = _fetch.nm(gene)

    cds_seq = None
    for feature in seqrecord.features:
        if feature.type == "CDS":
            cds_seq = feature.extract(seqrecord).seq
            break

    if cds_seq is None:
        raise ValueError(f"No CDS feature found in {seqrecord.id}")

    codon_start_0: int = (position - 1) * 3

    # Bounds check
    if codon_start_0 + 3 > len(cds_seq):
        raise ValueError(
            f"Position {position} exceeds CDS length "
            f"({len(cds_seq)} nt / {len(cds_seq) // 3} aa)"
        )

    actual_codon = str(cds_seq[codon_start_0 : codon_start_0 + 3])

    # Translate actual codon to amino acid
    if actual_codon in standard_dna_table.forward_table:
        actual_aa = standard_dna_table.forward_table[actual_codon]
    else:
        actual_aa = "*"

    # Validate reference amino acid
    if actual_aa != ref_aa:
        raise ValueError(
            f"Reference amino acid mismatch at position {position}: "
            f"expected {ref_aa}, found {actual_aa} "
            f"(codon {actual_codon} in {seqrecord.id})"
        )

    alt_codons = _AA_TO_CODONS.get(alt_aa, [])
    if not alt_codons:
        return []

    seen: set[str] = set()
    for ac in alt_codons:
        if ac == actual_codon:
            continue
        diff = _diff_position(actual_codon, ac)
        if diff is not None:
            pos_in_codon, ref_base, alt_base = diff
            c_pos = codon_start_0 + pos_in_codon + 1
            hgvs_full = f"{seqrecord.id}:c.{c_pos}{ref_base}>{alt_base}"
            if hgvs_full not in seen:
                seen.add(hgvs_full)

    return sorted(seen)
