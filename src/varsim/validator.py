"""HGVS variant description validator.

Provides syntactic and semantic validation for HGVS variant description
strings. Syntactic checks verify correct form (accession, prefix, coordinate,
allele symbols). Semantic checks verify consistency against a reference
sequence when provided.
"""

import re
from typing import Optional

from ._logging import get_logger
from .parser import parse, HGVSTag

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_VALID_ACC_PREFIXES = frozenset({'NM_', 'NP_', 'NC_', 'NG_', 'NR_'})
_VALID_PREFIXES = frozenset({'c.', 'p.', 'g.', 'n.', 'm.'})
_VALID_NUCLEOTIDES = set('ACGTRYSWKMBDHVNUacgtryswkmbdhvnu')
_VALID_AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY*')

# Accession format: two letters, underscore, digits, dot, digits
_ACC_FORMAT_RE = re.compile(r'^[A-Z]{2}_\d+\.\d+$')

# For extracting protein ref/alt when the parser doesn't (e.g. p.(M1V))
_PROTEIN_NOTATION_RE = re.compile(
    r'p\.\('
    r'(?P<ref>[A-Z*])'
    r'(?P<pos>\d+)'
    r'(?P<alt>[A-Z*=?])'
    r'(?:(?:fs|ext)[*]?\d*)?'
    r'\)$'
)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def is_valid_syntax(hgvs_str: str) -> bool:
    """Check if the HGVS string can be parsed without raising an error.

    Syntactic validation only — does not verify reference sequence
    consistency.

    Parameters
    ----------
    hgvs_str : str
        An HGVS variant description string (e.g. ``"NM_000207.3:c.1A>G"``).

    Returns
    -------
    bool
        ``True`` if the string parses successfully.

    Examples
    --------
    >>> is_valid_syntax("NM_000207.3:c.1A>G")
    True
    >>> is_valid_syntax("not a valid hgvs")
    False
    >>> is_valid_syntax("NM_000207.3:x.1A>G")
    False
    >>> is_valid_syntax("NP_000198.1:p.(M1?)")
    True
    """
    try:
        parse(hgvs_str)
        return True
    except (ValueError, IndexError):
        return False


def validate(hgvs_str: str) -> list[dict]:
    """Parse and validate an HGVS variant description string.

    Returns a list of issues found. An empty list means the variant is
    syntactically valid. Each issue is a dict with ``"severity"``
    (``"error"`` or ``"warning"``) and ``"message"`` (str).

    Checks performed:
    - Accession format (NP_, NM_, NC_, NG_, NR_ with version number)
    - Valid HGVS prefix (c., p., g., n., m.)
    - Coordinates are positive integers (negative allowed for 5' UTR)
    - Ref/alt alleles contain valid nucleotide or amino acid symbols
    - Variant type consistency (e.g. substitution has both ref *and* alt)
    - Protein variants use valid amino acid symbols (20 standard + ``*``)
    - ``?`` notation triggers a warning (ambiguous/uncertain)
    - Range positions are ordered (start <= end)

    Parameters
    ----------
    hgvs_str : str
        An HGVS variant description string.

    Returns
    -------
    list[dict]
        List of validation issues (empty when valid).

    Examples
    --------
    >>> validate("NM_000207.3:c.1A>G")
    []
    >>> issues = validate("NM_000207.3:x.1A>G")
    >>> len(issues) > 0
    True
    >>> issues[0]['severity']
    'error'
    >>> issues = validate("NM_000207.3:c.1?>G")
    >>> any(i['severity'] == 'warning' for i in issues)
    True
    >>> issues = validate("NM_000207.3:c.100_50del")
    >>> any('greater' in i['message'] for i in issues)
    True
    """
    issues: list[dict] = []
    logger.debug("Validating: %s", hgvs_str)

    # -- Parse (raise nothing) --
    try:
        tag = parse(hgvs_str)
    except ValueError as e:
        return [{"severity": "error", "message": str(e)}]
    except IndexError:
        return [{"severity": "error",
                 "message": f"Failed to parse: {hgvs_str}"}]

    # -- Accession checks --
    _check_accession(tag, issues)

    # -- Prefix check --
    if tag.prefix not in _VALID_PREFIXES:
        issues.append({
            "severity": "error",
            "message": (
                f"Invalid HGVS prefix: {tag.prefix!r}. "
                f"Expected one of {sorted(_VALID_PREFIXES)}"
            ),
        })

    # -- Coordinate checks --
    _check_coordinates(tag, issues)

    # -- Allele symbol checks --
    is_protein = tag.prefix.startswith('p')
    if is_protein:
        _check_protein_symbols(tag, issues, hgvs_str)
    else:
        _check_nucleotide_symbols(tag, issues)

    # -- Variant-type consistency --
    _check_variant_type_consistency(tag, issues, hgvs_str)

    # -- Uncertainty warnings --
    if tag.is_uncertain:
        issues.append({
            "severity": "warning",
            "message": "Variant uses '?' notation (ambiguous/uncertain)",
        })
    elif '?' in tag.original:
        # '?' may appear without the parser flagging is_uncertain
        issues.append({
            "severity": "warning",
            "message": "Variant contains '?' (ambiguous/uncertain)",
        })

    logger.debug("%d issues found", len(issues))
    return issues


def validate_semantic(hgvs_str: str,
                      ref_seq: Optional[str] = None) -> list[dict]:
    """Validate HGVS variant syntax and, optionally, reference consistency.

    When *ref_seq* is provided, additional checks are applied:

    - **c./g. variants**: the ref allele must match the reference sequence
      at the given coordinate.
    - **p. variants**: the original amino acid must match the reference at
      the given position.
    - Coordinates must fall within the reference sequence length.

    Parameters
    ----------
    hgvs_str : str
        An HGVS variant description string.
    ref_seq : str or None, optional
        Reference sequence (DNA for c./g. variants, protein for p. variants).
        When ``None`` only syntactic validation is performed.

    Returns
    -------
    list[dict]
        List of validation issues (empty when fully valid).

    Examples
    --------
    >>> validate_semantic("NM_000207.3:c.1A>G", "ATCG")
    []
    >>> issues = validate_semantic("NM_000207.3:c.1A>G", "TTCG")
    >>> any('mismatch' in i['message'] for i in issues)
    True
    >>> issues = validate_semantic("NM_000207.3:c.100A>G", "ATCG")
    >>> any('exceeds' in i['message'] for i in issues)
    True
    >>> # Without ref_seq, falls back to syntax-only:
    >>> validate_semantic("NM_000207.3:c.1A>G")
    []
    """
    issues = validate(hgvs_str)
    logger.debug("Validating semantically: %s", hgvs_str)

    if ref_seq is None:
        return issues

    # Attempt parse even if syntax errors exist — best-effort semantic checks
    try:
        tag = parse(hgvs_str)
    except (ValueError, IndexError):
        return issues

    seq_len = len(ref_seq)
    is_protein = tag.prefix.startswith('p')

    # -- Bounds checks --
    _check_bounds(tag, seq_len, is_protein, issues)

    # -- Reference allele match --
    ref = tag.ref
    if ref is None and is_protein:
        ref, _ = _extract_protein_alleles(hgvs_str)
    if ref is not None and tag.start_pos > 0:
        _check_ref_match(tag, ref, ref_seq, seq_len, is_protein, issues)

    return issues


def is_valid(hgvs_str: str, ref_seq: Optional[str] = None) -> bool:
    """Check whether an HGVS variant is fully valid.

    Returns ``True`` when both syntactic and (if *ref_seq* provided)
    semantic validation pass — i.e. no ``"error"`` severity issues
    remain.  Warnings alone do not cause ``False``.

    Parameters
    ----------
    hgvs_str : str
        An HGVS variant description string.
    ref_seq : str or None, optional
        Reference sequence for semantic validation.

    Returns
    -------
    bool
        ``True`` if all validation checks pass.

    Examples
    --------
    >>> is_valid("NM_000207.3:c.1A>G", "ATCG")
    True
    >>> is_valid("NM_000207.3:c.1A>G", "TTCG")
    False
    >>> is_valid("invalid hgvs")
    False
    >>> # Warning only — still valid:
    >>> is_valid("NM_000207.3:c.1?>G")
    True
    """
    if ref_seq is not None:
        issues = validate_semantic(hgvs_str, ref_seq)
    else:
        issues = validate(hgvs_str)

    return not any(issue['severity'] == 'error' for issue in issues)


# ---------------------------------------------------------------------------
# Internal helpers — syntactic
# ---------------------------------------------------------------------------

def _check_accession(tag: HGVSTag, issues: list[dict]) -> None:
    """Validate accession number format."""
    for label, acc in [("accession", tag.acc),
                       ("genomic accession", tag.genomic_acc)]:
        if acc is None:
            continue
        if not _ACC_FORMAT_RE.match(acc):
            issues.append({
                "severity": "error",
                "message": (
                    f"Invalid {label} format: {acc!r}. "
                    f"Expected XX_####.# (e.g. NM_000207.3)"
                ),
            })
            continue
        acc_prefix = acc[:3]
        if acc_prefix not in _VALID_ACC_PREFIXES:
            issues.append({
                "severity": "error",
                "message": (
                    f"Invalid {label} prefix: {acc_prefix!r}. "
                    f"Expected one of {sorted(_VALID_ACC_PREFIXES)}"
                ),
            })


def _check_coordinates(tag: HGVSTag, issues: list[dict]) -> None:
    """Validate coordinate fields."""
    # Range ordering
    if tag.end_pos is not None and tag.start_pos > tag.end_pos:
        issues.append({
            "severity": "error",
            "message": (
                f"Start position ({tag.start_pos}) is greater than "
                f"end position ({tag.end_pos})"
            ),
        })

    # Zero check — HGVS is 1-based; zero is invalid except for
    # special cases.  The parser will not normally produce 0 unless
    # the string explicitly used 0 (which is non-standard).
    if tag.start_pos == 0:
        issues.append({
            "severity": "error",
            "message": "Position 0 is not valid HGVS (1-based numbering)",
        })

    # Intronic offset of 0 is redundant but not an error
    if tag.start_offset is not None and tag.start_offset == 0:
        issues.append({
            "severity": "warning",
            "message": "Intronic offset of 0 is redundant",
        })
    if tag.end_offset is not None and tag.end_offset == 0:
        issues.append({
            "severity": "warning",
            "message": "Intronic offset of 0 is redundant",
        })


def _check_nucleotide_symbols(tag: HGVSTag, issues: list[dict]) -> None:
    """Validate nucleotide allele symbols (IUPAC extended allowed)."""
    for label, allele in [("ref", tag.ref), ("alt", tag.alt)]:
        if allele is None:
            continue
        # Filter out '?' — handled by uncertainty warning
        invalid = set(allele) - _VALID_NUCLEOTIDES - {'?'}
        if invalid:
            issues.append({
                "severity": "error",
                "message": (
                    f"Invalid nucleotide symbol(s) in {label} allele "
                    f"'{allele}': {sorted(invalid)}"
                ),
            })


def _extract_protein_alleles(hgvs_str: str) -> tuple[Optional[str], Optional[str]]:
    """Extract ref/alt amino acids from protein parenthetical notation.

    The parser sometimes misses ref/alt for formats like ``p.(M1V)``
    because coordinate extraction consumes part of the protein notation.
    This helper recovers them directly from the original string.
    """
    m = _PROTEIN_NOTATION_RE.search(hgvs_str)
    if m:
        return m.group('ref'), m.group('alt')
    return None, None


def _check_protein_symbols(tag: HGVSTag, issues: list[dict],
                           hgvs_str: str = '') -> None:
    """Validate amino acid symbols in ref/alt, including = and ? for alt."""
    ref = tag.ref
    alt = tag.alt

    # If the parser missed ref/alt, try to extract from original string
    if ref is None and alt is None:
        ref, alt = _extract_protein_alleles(hgvs_str)

    for label, allele in [("ref", ref), ("alt", alt)]:
        if allele is None:
            continue
        # alt may legitimately contain '=' (silent) or '?' (uncertain)
        allowed_extra = {'=', '?'} if label == 'alt' else set()
        invalid = set(allele) - _VALID_AMINO_ACIDS - allowed_extra
        if invalid:
            issues.append({
                "severity": "error",
                "message": (
                    f"Invalid amino acid symbol(s) in {label} "
                    f"'{allele}': {sorted(invalid)}"
                ),
            })


def _check_variant_type_consistency(tag: HGVSTag, issues: list[dict],
                                    hgvs_str: str = '') -> None:
    """Verify ref/alt presence matches the declared variant type."""
    vt = tag.variant_type
    ref = tag.ref
    alt = tag.alt
    is_protein = tag.prefix.startswith('p')

    # Attempt recovery for protein variants where the parser missed alleles
    if is_protein and ref is None and alt is None:
        ref, alt = _extract_protein_alleles(hgvs_str)

    if vt == 'substitution':
        if ref is None:
            issues.append({
                "severity": "error",
                "message": "Substitution variant missing ref allele",
            })
        if alt is None:
            issues.append({
                "severity": "error",
                "message": "Substitution variant missing alt allele",
            })
        if (ref is not None and alt is not None
                and ref.upper() == alt.upper()):
            issues.append({
                "severity": "warning",
                "message": (
                    f"Ref and alt alleles are identical: {ref!r}"
                ),
            })

    elif vt == 'insertion':
        if alt is None:
            issues.append({
                "severity": "error",
                "message": "Insertion variant missing inserted bases",
            })

    elif vt == 'delins':
        if alt is None:
            issues.append({
                "severity": "error",
                "message": "Deletion-insertion variant missing inserted bases",
            })

    elif vt == 'frameshift':
        if ref is None:
            issues.append({
                "severity": "error",
                "message": "Frameshift variant missing reference amino acid",
            })


# ---------------------------------------------------------------------------
# Internal helpers — semantic
# ---------------------------------------------------------------------------

def _check_bounds(tag: HGVSTag, seq_len: int, is_protein: bool,
                  issues: list[dict]) -> None:
    """Check that coordinates fall within the reference sequence length."""
    # 5' UTR negative positions — check absolute value
    if tag.start_pos < 0 and not is_protein:
        if abs(tag.start_pos) > seq_len:
            issues.append({
                "severity": "error",
                "message": (
                    f"5' UTR position {tag.start_pos} exceeds "
                    f"sequence length ({seq_len})"
                ),
            })
    elif tag.start_pos > seq_len:
        issues.append({
            "severity": "error",
            "message": (
                f"Start position {tag.start_pos} exceeds "
                f"sequence length ({seq_len})"
            ),
        })

    if tag.end_pos is not None and tag.end_pos > seq_len:
        issues.append({
            "severity": "error",
            "message": (
                f"End position {tag.end_pos} exceeds "
                f"sequence length ({seq_len})"
            ),
        })


def _check_ref_match(tag: HGVSTag, ref: str, ref_seq: str, seq_len: int,
                     is_protein: bool, issues: list[dict]) -> None:
    """Verify that the ref allele matches the reference at the given position."""
    if tag.start_pos > seq_len or tag.start_pos < 1:
        return  # bounds error already reported

    if is_protein:
        _check_protein_ref_match(tag, ref, ref_seq, issues)
    else:
        _check_nucleotide_ref_match(tag, ref, ref_seq, seq_len, issues)


def _check_protein_ref_match(tag: HGVSTag, ref: str, ref_seq: str,
                             issues: list[dict]) -> None:
    """Compare ref amino acid against protein sequence."""
    ref_expected = ref_seq[tag.start_pos - 1]
    if ref_expected.upper() != ref.upper():
        issues.append({
            "severity": "error",
            "message": (
                f"Reference amino acid mismatch at position "
                f"{tag.start_pos}: expected {ref_expected!r}, "
                f"found {ref!r}"
            ),
        })


def _check_nucleotide_ref_match(tag: HGVSTag, ref: str, ref_seq: str,
                                seq_len: int, issues: list[dict]) -> None:
    """Compare ref allele against nucleotide reference sequence."""
    if tag.end_pos is not None and tag.end_pos <= seq_len:
        # Range variant
        ref_expected = ref_seq[tag.start_pos - 1:tag.end_pos]
        if ref_expected.upper() != ref.upper():
            issues.append({
                "severity": "error",
                "message": (
                    f"Reference allele mismatch at "
                    f"{tag.start_pos}_{tag.end_pos}: "
                    f"expected {ref_expected!r}, found {ref!r}"
                ),
            })
    else:
        # Single-position variant — match exact length of ref allele
        ref_len = len(ref)
        start_idx = tag.start_pos - 1
        ref_expected = ref_seq[start_idx:start_idx + ref_len]
        if ref_expected.upper() != ref.upper():
            issues.append({
                "severity": "error",
                "message": (
                    f"Reference allele mismatch at position "
                    f"{tag.start_pos}: expected {ref_expected!r}, "
                    f"found {ref!r}"
                ),
            })
