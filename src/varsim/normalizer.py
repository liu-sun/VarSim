"""HGVS variant normalizer.

Implements HGVS nomenclature normalization rules:
- 3-prime shifting (deletions and duplications only; per the HGVS
  recommendations the 3' rule does NOT apply to substitutions)
- Insertion-to-duplication conversion (ins → dup)
- Allele representation minimization
- Range normalization (start ≤ end, single-position collapse)

All algorithms follow HGVS recommendations (v20.05 and later).
"""

from __future__ import annotations

from typing import Optional

from ._logging import get_logger
from .parser import HGVSTag, parse

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Sentinels for "no value" vs "use default"
# ---------------------------------------------------------------------------

_UNSET = object()  # sentinel: field was not explicitly provided


# ---------------------------------------------------------------------------
# HGVS string reconstruction helpers
# ---------------------------------------------------------------------------

def _reconstruct_hgvs(
    tag: HGVSTag,
    start_pos: object = _UNSET,
    end_pos: object = _UNSET,
    ref: object = _UNSET,
    alt: object = _UNSET,
    variant_type: object = _UNSET,
) -> str:
    """Reconstruct an HGVS string from a tag with optional field overrides.

    Parameters
    ----------
    tag : HGVSTag
        Original parsed tag providing defaults.
    start_pos : int or None, optional
        Override start position (1-based).  ``None`` means *no* start (rare).
    end_pos : int or None, optional
        Override end position.  ``None`` explicitly collapses to a single
        position (no range).
    ref : str or None, optional
        Override reference allele.
    alt : str or None, optional
        Override alternate allele.
    variant_type : str or None, optional
        Override variant type.

    Returns
    -------
    str
        Reconstructed HGVS string.
    """
    # --- Accession + prefix ---
    if tag.genomic_acc:
        prefix_part = f"{tag.genomic_acc}({tag.acc}):{tag.prefix}"
    else:
        prefix_part = f"{tag.acc}:{tag.prefix}"

    # --- Coordinates ---
    sp: int = start_pos if start_pos is not _UNSET else tag.start_pos  # type: ignore[assignment]
    ep: Optional[int] = end_pos if end_pos is not _UNSET else tag.end_pos  # type: ignore[assignment]

    # Build coordinate string
    if ep is not None and ep != sp:
        coord = f"{sp}_{ep}"
    else:
        coord = str(sp)

    # Append intronic offsets if present (only from original tag for simplicity)
    so = tag.start_offset
    eo = tag.end_offset
    if so is not None and eo is not None:
        coord = f"{sp}{so:+d}_{ep}{eo:+d}" if (ep is not None and ep != sp) else f"{sp}{so:+d}"
    elif so is not None:
        coord = f"{sp}{so:+d}"
    elif eo is not None:
        coord = f"{sp}_{ep}{eo:+d}" if (ep is not None and ep != sp) else f"{sp}{eo:+d}"

    # --- Variant part ---
    vt: str = variant_type if variant_type is not _UNSET else tag.variant_type  # type: ignore[assignment]
    r: Optional[str] = ref if ref is not _UNSET else tag.ref  # type: ignore[assignment]
    a: Optional[str] = alt if alt is not _UNSET else tag.alt  # type: ignore[assignment]

    if vt == "substitution":
        var_part = f"{r}>{a}"
    elif vt == "deletion":
        if r:
            var_part = f"del{r}"
        else:
            var_part = "del"
    elif vt == "insertion":
        var_part = f"ins{a}"
    elif vt == "delins":
        var_part = f"delins{a}"
    elif vt == "duplication":
        if a:
            var_part = f"dup{a}"
        else:
            var_part = "dup"
    elif vt == "inversion":
        if a:
            var_part = f"inv{a}"
        else:
            var_part = "inv"
    else:
        var_part = _fallback_variant_part(tag)

    return f"{prefix_part}{coord}{var_part}"


def _fallback_variant_part(tag: HGVSTag) -> str:
    """Extract the variant part from the original string as a last resort."""
    import re as _re
    m = _re.search(r'[cgmn]\.\d+.*$', tag.original)
    if m:
        coord_and_var = m.group(0)
        return _re.sub(r'^\d+([_]\d+)?([+-]\d+)?([_]\d+[+-]\d+)?', '', coord_and_var)
    return ""


def _make_single_pos_variant(
    tag: HGVSTag,
    new_pos: int,
    ref: object = _UNSET,
    alt: object = _UNSET,
    variant_type: object = _UNSET,
) -> str:
    """Convenience: reconstruct with single position (no range)."""
    return _reconstruct_hgvs(
        tag,
        start_pos=new_pos,
        end_pos=None,
        ref=ref,
        alt=alt,
        variant_type=variant_type,
    )


def _make_range_variant(
    tag: HGVSTag,
    start_pos: int,
    end_pos: int,
    ref: object = _UNSET,
    alt: object = _UNSET,
    variant_type: object = _UNSET,
) -> str:
    """Convenience: reconstruct with a range."""
    return _reconstruct_hgvs(
        tag,
        start_pos=start_pos,
        end_pos=end_pos,
        ref=ref,
        alt=alt,
        variant_type=variant_type,
    )


# ---------------------------------------------------------------------------
# 3-prime shift
# ---------------------------------------------------------------------------

def normalize_3prime_shift(hgvs_str: str, ref_seq: str) -> str:
    """Shift a variant as far 3' as possible while maintaining equivalence.

    Per the HGVS recommendations, the 3' rule applies to deletions,
    duplications, and insertions that are rewritten as duplications; it
    does NOT apply to substitutions, which are left unchanged.

    For deletions: shift right while the deleted sequence repeats
    immediately 3' of the deletion.
    For duplications: shift right while the duplicated sequence repeats
    immediately 3'.

    Parameters
    ----------
    hgvs_str : str
        HGVS variant string (e.g. ``"NM_000207.3:c.1A>G"``).
    ref_seq : str
        Reference sequence string used to verify equivalence.

    Returns
    -------
    str
        Normalized HGVS string with 3-prime shifted position.

    Raises
    ------
    ValueError
        If *hgvs_str* cannot be parsed or *ref_seq* is empty.

    Examples
    --------
    >>> # Substitutions are NOT shifted (3' rule does not apply)
    >>> normalize_3prime_shift("NM_000207.3:c.1A>G", "AAGC")
    'NM_000207.3:c.1A>G'

    >>> # Deletion shift
    >>> normalize_3prime_shift("NM_000207.3:c.2_3delAG", "CAGAGAG")
    'NM_000207.3:c.6_7delAG'

    >>> # Empty ref_seq raises
    >>> normalize_3prime_shift("NM_000207.3:c.1A>G", "")
    Traceback (most recent call last):
        ...
    ValueError: ref_seq must be non-empty for 3-prime shift
    """
    if not ref_seq:
        raise ValueError("ref_seq must be non-empty for 3-prime shift")

    logger.debug("3-prime shifting: %s", hgvs_str)
    tag = parse(hgvs_str)

    # Only applicable to nucleotide-level variants
    if tag.prefix not in ("c.", "g.", "n.", "m."):
        return hgvs_str

    vt = tag.variant_type

    if vt == "substitution":
        # The HGVS 3' rule does not apply to substitutions.
        return hgvs_str
    elif vt == "deletion":
        return _shift_deletion_3prime(tag, ref_seq)
    elif vt == "duplication":
        return _shift_duplication_3prime(tag, ref_seq)
    else:
        # Other variant types not shifted
        return hgvs_str


def _shift_deletion_3prime(tag: HGVSTag, ref_seq: str) -> str:
    """Shift a deletion right (3') as far as possible.

    Algorithm: shift one base at a time.  A deletion of length L at [p, q]
    can be shifted to [p+1, q+1] iff the base immediately after the deletion
    (at q+1) equals the first deleted base (at p).  Repeat.
    """
    sp = tag.start_pos  # 1-based
    ep = tag.end_pos if tag.end_pos is not None else sp  # 1-based

    # Determine deleted sequence
    if tag.ref:
        deleted = tag.ref
    else:
        deleted_len = ep - sp + 1
        deleted = ref_seq[sp - 1 : sp - 1 + deleted_len]

    if not deleted:
        return tag.original

    # Shift one base at a time: while the base just past the deletion
    # equals the first base of the deleted sequence
    while ep < len(ref_seq) and ref_seq[sp - 1] == ref_seq[ep]:
        sp += 1
        ep += 1
        # Rotate the deleted string left by 1 to maintain correct allele
        deleted = deleted[1:] + deleted[:1]

    if ep == sp:
        return _make_single_pos_variant(tag, new_pos=sp, ref=deleted, variant_type="deletion")
    return _make_range_variant(tag, start_pos=sp, end_pos=ep, ref=deleted, variant_type="deletion")


def _shift_duplication_3prime(tag: HGVSTag, ref_seq: str) -> str:
    """Shift a duplication right (3') as far as possible.

    Algorithm: shift one base at a time.  A duplication of length L at [p, q]
    can be shifted to [p+1, q+1] iff the base immediately after the duplicated
    region (at q+1) equals the first base of the duplicated sequence (at p).
    """
    sp = tag.start_pos  # 1-based
    ep = tag.end_pos if tag.end_pos is not None else sp  # 1-based

    # Determine duplicated sequence
    if tag.alt:
        dup_seq = tag.alt
    else:
        dup_len = ep - sp + 1
        dup_seq = ref_seq[sp - 1 : sp - 1 + dup_len]

    if not dup_seq:
        return tag.original

    # Shift one base at a time
    while ep < len(ref_seq) and ref_seq[sp - 1] == ref_seq[ep]:
        sp += 1
        ep += 1
        dup_seq = dup_seq[1:] + dup_seq[:1]

    if ep == sp:
        return _make_single_pos_variant(tag, new_pos=sp, alt=dup_seq if dup_seq else None, variant_type="duplication")
    return _make_range_variant(tag, start_pos=sp, end_pos=ep, alt=dup_seq if dup_seq else None, variant_type="duplication")


# ---------------------------------------------------------------------------
# Insertion to duplication conversion
# ---------------------------------------------------------------------------

def ins_to_dup(hgvs_str: str, ref_seq: str) -> str:
    """Convert an insertion to a duplication if the inserted bases match
    the sequence immediately 5' of the insertion site.

    Per HGVS recommendations, if an insertion duplicates the sequence
    preceding it, the variant should be described as a duplication (dup).

    Parameters
    ----------
    hgvs_str : str
        HGVS variant string with an insertion.
    ref_seq : str
        Reference sequence.

    Returns
    -------
    str
        Converted HGVS string (dup) if applicable, otherwise the original.

    Raises
    ------
    ValueError
        If *ref_seq* is empty.

    Examples
    --------
    >>> # Insertion duplicates the sequence immediately 5' of the site
    >>> ins_to_dup("NM_000207.3:c.4_5insA", "TAAA")
    'NM_000207.3:c.4dupA'

    >>> # Insertion duplicates the sequence immediately 3' (3'-most placed)
    >>> ins_to_dup("NM_000207.3:c.1_2insC", "CCT")
    'NM_000207.3:c.2dupC'

    >>> # Inserted bases match neither flanking sequence
    >>> ins_to_dup("NM_000207.3:c.5_6insAC", "TGACAC")
    'NM_000207.3:c.5_6insAC'

    >>> # Inserted bases do not match preceding sequence
    >>> ins_to_dup("NM_000207.3:c.4_5insG", "TAAC")
    'NM_000207.3:c.4_5insG'

    >>> # Not an insertion
    >>> ins_to_dup("NM_000207.3:c.1A>G", "ACGT")
    'NM_000207.3:c.1A>G'

    >>> # Empty ref_seq raises
    >>> ins_to_dup("NM_000207.3:c.4_5insA", "")
    Traceback (most recent call last):
        ...
    ValueError: ref_seq must be non-empty for ins-to-dup conversion
    """
    if not ref_seq:
        raise ValueError("ref_seq must be non-empty for ins-to-dup conversion")

    tag = parse(hgvs_str)

    # Only applicable to insertions on nucleotide sequences
    if tag.variant_type != "insertion":
        return hgvs_str
    if tag.prefix not in ("c.", "g.", "n.", "m."):
        return hgvs_str

    inserted = tag.alt
    if not inserted:
        return hgvs_str

    ins_len = len(inserted)
    sp = tag.start_pos

    # The insertion is between start_pos and end_pos; the sequence 5'
    # (upstream) of the insertion therefore ENDS at start_pos, i.e. it
    # occupies positions [sp - ins_len + 1, sp].
    upstream_start = sp - ins_len + 1
    up_match = False
    if upstream_start >= 1:
        upstream_seq = ref_seq[upstream_start - 1 : sp]
        up_match = upstream_seq.upper() == inserted.upper()

    # HGVS also permits ins→dup when the inserted sequence duplicates the
    # sequence immediately 3' (downstream) of the insertion site; per the
    # 3' rule the duplication is then placed at the 3'-most position.
    downstream_seq = ref_seq[sp : sp + ins_len] if sp + ins_len <= len(ref_seq) else ''
    down_match = downstream_seq.upper() == inserted.upper()

    if up_match or down_match:
        # Convert to duplication (clear alt — bases inferred from reference)
        if down_match and not up_match:
            new_start = sp + 1
        else:
            new_start = upstream_start
        if ins_len == 1:
            result = _make_single_pos_variant(
                tag, new_pos=new_start, alt=None, variant_type="duplication"
            )
        else:
            result = _make_range_variant(
                tag,
                start_pos=new_start,
                end_pos=new_start + ins_len - 1,
                alt=None,
                variant_type="duplication",
            )
        # Apply the 3' rule to the resulting duplication.
        return normalize_3prime_shift(result, ref_seq)

    return hgvs_str


# ---------------------------------------------------------------------------
# Allele minimization
# ---------------------------------------------------------------------------

def _minimize_allele(hgvs_str: str, ref_seq: str) -> str:
    """Minimize allele representation by trimming identical flanking bases.

    For deletions: if the first/last deleted base matches the flanking base,
    the variant can be described with a shorter allele representation.
    This is essentially the same as 3-prime shifting for deletions,
    so we delegate to normalize_3prime_shift for deletions.

    For insertions/duplications: similar logic applies.

    Parameters
    ----------
    hgvs_str : str
        HGVS variant string.
    ref_seq : str
        Reference sequence.

    Returns
    -------
    str
        Minimized HGVS string.
    """
    tag = parse(hgvs_str)

    if tag.prefix not in ("c.", "g.", "n.", "m."):
        return hgvs_str

    if tag.variant_type == "deletion":
        # For deletions, minimize by trimming identical flanking bases
        sp = tag.start_pos
        ep = tag.end_pos if tag.end_pos is not None else sp

        if tag.ref:
            deleted = tag.ref
        else:
            deleted_len = ep - sp + 1
            deleted = ref_seq[sp - 1 : sp - 1 + deleted_len]

        if not deleted:
            return hgvs_str

        # Trim from left: while first deleted base == base immediately 3' of deletion
        while deleted and ep < len(ref_seq) and deleted[0].upper() == ref_seq[ep].upper():
            deleted = deleted[1:]
            sp += 1
            if not deleted:
                break

        # Trim from right: while last deleted base == base immediately 5' of deletion
        while deleted and sp > 1 and deleted[-1].upper() == ref_seq[sp - 2].upper():
            deleted = deleted[:-1]
            ep -= 1
            if not deleted:
                break

        if not deleted:
            return hgvs_str  # Can't describe without deleted bases

        if sp == ep:
            return _make_single_pos_variant(tag, new_pos=sp, ref=deleted, variant_type="deletion")
        return _make_range_variant(tag, start_pos=sp, end_pos=ep, ref=deleted, variant_type="deletion")

    return hgvs_str


# ---------------------------------------------------------------------------
# Range normalization
# ---------------------------------------------------------------------------

def _normalize_range(hgvs_str: str) -> str:
    """Ensure start ≤ end and collapse single-position ranges.

    Parameters
    ----------
    hgvs_str : str
        HGVS variant string.

    Returns
    -------
    str
        Range-normalized HGVS string.
    """
    tag = parse(hgvs_str)

    sp = tag.start_pos
    ep = tag.end_pos

    # No range, nothing to do
    if ep is None:
        return hgvs_str

    # Ensure start ≤ end
    if sp > ep:
        sp, ep = ep, sp

    # Collapse single-position range (e.g. 4_4 → 4)
    if sp == ep:
        return _make_single_pos_variant(
            tag,
            new_pos=sp,
            variant_type=tag.variant_type,
        )

    # Only reconstruct if positions changed
    if sp != tag.start_pos or ep != tag.end_pos:
        return _make_range_variant(
            tag,
            start_pos=sp,
            end_pos=ep,
            variant_type=tag.variant_type,
        )

    return hgvs_str


# ---------------------------------------------------------------------------
# Full normalization pipeline
# ---------------------------------------------------------------------------

def normalize(hgvs_str: str, ref_seq: Optional[str] = None) -> str:
    """Normalize an HGVS variant description to canonical form.

    Applies the following steps in order:

    1. **3-prime shift**: Shift deletions and duplications as far 3' as
       possible while maintaining equivalence. Per the HGVS
       recommendations, substitutions are NOT shifted.
    2. **ins → dup conversion**: If an insertion duplicates the preceding
       sequence, rewrite it as a duplication.
    3. **Minimize allele representation**: Trim identical flanking bases
       to use the minimal allele length.
    4. **Range normalization**: Ensure start ≤ end, and collapse
       single-position ranges to simple position format.

    Parameters
    ----------
    hgvs_str : str
        HGVS variant description string (e.g. ``"NM_000207.3:c.1A>G"``).
    ref_seq : str or None
        Reference sequence for the transcript. Required for 3-prime shift
        and ins→dup conversion. If *None*, only range normalization is
        performed.

    Returns
    -------
    str
        Canonical HGVS variant string.

    Raises
    ------
    ValueError
        If the HGVS string cannot be parsed.

    Examples
    --------
    >>> # Substitutions are not 3'-shifted
    >>> normalize("NM_000207.3:c.1A>G", "AAGC")
    'NM_000207.3:c.1A>G'

    >>> normalize("NM_000207.3:c.4_5insA", "TAAA")
    'NM_000207.3:c.4dupA'

    >>> # No ref_seq: only range normalization
    >>> normalize("NM_000207.3:c.5_4del")
    'NM_000207.3:c.4_5del'

    >>> normalize("NM_000207.3:c.3_3A>G")
    'NM_000207.3:c.3A>G'

    >>> # Protein variants are not normalized at nucleotide level
    >>> normalize("NP_000198.1:p.M1V")
    'NP_000198.1:p.M1V'

    >>> # Invalid input
    >>> normalize("not a variant")
    Traceback (most recent call last):
        ...
    ValueError: Cannot parse accession/prefix from: not a variant
    """
    logger.debug("Normalizing: %s", hgvs_str)
    # Parse (validates input)
    tag = parse(hgvs_str)

    # Protein variants: return as-is
    if tag.prefix == "p.":
        return hgvs_str

    result = hgvs_str

    if ref_seq:
        # Step 1: 3-prime shift
        result = normalize_3prime_shift(result, ref_seq)

        # Step 2: ins → dup
        result = ins_to_dup(result, ref_seq)

        # Step 3: Minimize allele
        result = _minimize_allele(result, ref_seq)

    # Step 4: Range normalization (always applied)
    result = _normalize_range(result)

    logger.debug("Normalized → %s", result)
    return result
