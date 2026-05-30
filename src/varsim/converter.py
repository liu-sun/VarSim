"""HGVS format conversion functions.

Converts between HGVS variant descriptions and other common formats:
VCF, SPDI, and cross-system translations (coding ↔ protein).

>>> from varsim.converter import hgvs_to_vcf, vcf_to_hgvs
>>> from varsim.converter import hgvs_to_spdi, spdi_to_hgvs

Conversions
-----------
- ``hgvs_to_vcf`` — HGVS g./c. string → VCF dict
- ``vcf_to_hgvs`` — VCF record → HGVS string
- ``hgvs_to_spdi`` — HGVS string → SPDI string
- ``spdi_to_hgvs`` — SPDI string → HGVS string
- ``c_to_p`` — coding HGVS → protein HGVS
- ``p_to_c_range`` — protein HGVS → all possible coding HGVS strings
"""

from ._logging import get_logger
from .parser import parse
from . import backtranslate

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# HGVS → VCF
# ---------------------------------------------------------------------------

def hgvs_to_vcf(hgvs_str: str, chrom: str = None) -> dict:
    """Convert an HGVS g. or c. variant to VCF format.

    Parameters
    ----------
    hgvs_str : str
        HGVS variant string (e.g. ``"NC_000006.12:g.123456A>G"``
        or ``"NM_000207.3:c.1A>G"``).
    chrom : str, optional
        Chromosome name for VCF CHROM column. Required for c. variants;
        for g. variants defaults to the NC_ accession.

    Returns
    -------
    dict
        VCF-style record with keys ``CHROM``, ``POS``, ``ID``, ``REF``, ``ALT``.

    Raises
    ------
    ValueError
        If the variant type cannot be converted to VCF.

    Examples
    --------
    Substitution (g.):
    >>> hgvs_to_vcf("NC_000006.12:g.123456A>G")
    {'CHROM': 'NC_000006.12', 'POS': 123456, 'ID': '.', 'REF': 'A', 'ALT': 'G'}

    Substitution (c.):
    >>> hgvs_to_vcf("NM_000207.3:c.1A>G", chrom="chr11")
    {'CHROM': 'chr11', 'POS': 1, 'ID': '.', 'REF': 'A', 'ALT': 'G'}

    Deletion with bases:
    >>> hgvs_to_vcf("NC_000006.12:g.123456delA")
    {'CHROM': 'NC_000006.12', 'POS': 123455, 'ID': '.', 'REF': 'NA', 'ALT': 'N'}

    Insertion:
    >>> hgvs_to_vcf("NC_000006.12:g.123456_123457insACGT")
    {'CHROM': 'NC_000006.12', 'POS': 123456, 'ID': '.', 'REF': 'N', 'ALT': 'NACGT'}

    Delins:
    >>> hgvs_to_vcf("NC_000006.12:g.123456_123457delinsAC")
    {'CHROM': 'NC_000006.12', 'POS': 123456, 'ID': '.', 'REF': 'N', 'ALT': 'AC'}
    """
    logger.debug("Converting HGVS→VCF: %s", hgvs_str)
    tag = parse(hgvs_str)

    if tag.prefix not in ('g.', 'c.', 'n.', 'm.'):
        raise ValueError(
            f"Only g./c./n./m. prefixes are supported for VCF conversion, "
            f"got {tag.prefix}"
        )

    # Determine CHROM
    vcf_chrom = chrom or tag.genomic_acc or tag.acc

    # Handle intronic offsets — not directly expressible in VCF
    if tag.start_offset is not None:
        raise ValueError(
            f"Intronic offset variants (like {hgvs_str}) "
            f"cannot be directly converted to VCF without genomic context"
        )

    vt = tag.variant_type
    start = tag.start_pos
    end = tag.end_pos
    ref = tag.ref
    alt = tag.alt

    # --- Substitution ---
    if vt == 'substitution':
        if ref is None or alt is None:
            raise ValueError(f"Substitution missing ref/alt: {hgvs_str}")
        return {
            "CHROM": vcf_chrom,
            "POS": start,
            "ID": ".",
            "REF": ref.upper(),
            "ALT": alt.upper(),
        }

    # --- Deletion ---
    if vt == 'deletion':
        # VCF left-pads: POS = start-1, REF = upstream_base + deleted, ALT = upstream_base
        # Without ref_seq, use 'N' as placeholder upstream base
        pos = start - 1
        if pos < 1:
            pos = 1  # Clamp to positive
        deleted = ref if ref else ''
        return {
            "CHROM": vcf_chrom,
            "POS": pos,
            "ID": ".",
            "REF": "N" + deleted.upper(),
            "ALT": "N",
        }

    # --- Insertion ---
    if vt == 'insertion':
        if alt is None:
            raise ValueError(f"Insertion missing alt bases: {hgvs_str}")
        return {
            "CHROM": vcf_chrom,
            "POS": start,
            "ID": ".",
            "REF": "N",
            "ALT": "N" + alt.upper(),
        }

    # --- Delins ---
    if vt == 'delins':
        deleted = ref if ref else ''
        inserted = alt if alt else ''
        return {
            "CHROM": vcf_chrom,
            "POS": start,
            "ID": ".",
            "REF": deleted.upper() if deleted else 'N',
            "ALT": inserted.upper() if inserted else '.',
        }

    # --- Duplication ---
    if vt == 'duplication':
        dup_bases = alt if alt else ''
        return {
            "CHROM": vcf_chrom,
            "POS": start,
            "ID": ".",
            "REF": "N",
            "ALT": "N" + dup_bases.upper() if dup_bases else "<DUP>",
        }

    # --- Inversion ---
    if vt == 'inversion':
        return {
            "CHROM": vcf_chrom,
            "POS": start,
            "ID": ".",
            "REF": "N",
            "ALT": "<INV>",
        }

    # --- Unsupported ---
    raise ValueError(
        f"Variant type '{vt}' cannot be converted to VCF: {hgvs_str}"
    )


# ---------------------------------------------------------------------------
# VCF → HGVS
# ---------------------------------------------------------------------------

def vcf_to_hgvs(chrom: str, pos: int, ref: str, alt: str,
                acc: str = None, prefix: str = "g.") -> str:
    """Convert a VCF record to an HGVS variant description string.

    Strips the common prefix from *ref* and *alt*, then determines the
    variant type from the trimmed allele lengths.

    Parameters
    ----------
    chrom : str
        Chromosome or contig name (used as accession if *acc* not given).
    pos : int
        1-based VCF position.
    ref : str
        Reference allele string.
    alt : str
        Alternate allele string. Use ``"."`` for a deletion with no
        alternate bases.
    acc : str, optional
        HGVS accession to use instead of *chrom*.
    prefix : str
        HGVS prefix (default ``"g."``).

    Returns
    -------
    str
        HGVS variant description.

    Examples
    --------
    Substitution:
    >>> vcf_to_hgvs("chr11", 5227002, "A", "G", acc="NC_000011.10")
    'NC_000011.10:g.5227002A>G'

    Insertion:
    >>> vcf_to_hgvs("chr11", 5227002, "A", "ATCG", acc="NC_000011.10")
    'NC_000011.10:g.5227002_5227003insTCG'

    Deletion:
    >>> vcf_to_hgvs("chr11", 5227001, "ATCG", "A", acc="NC_000011.10")
    'NC_000011.10:g.5227002_5227004delTCG'

    Delins:
    >>> vcf_to_hgvs("chr11", 5227002, "AT", "GC", acc="NC_000011.10")
    'NC_000011.10:g.5227002_5227003delinsGC'
    """
    logger.debug("Converting VCF→HGVS: %s:%d %s>%s", chrom, pos, ref, alt)
    # Use provided accession or fall back to chrom
    hgvs_acc = acc if acc else chrom

    # Normalize ALT
    if alt == '.':
        alt = ''

    # Strip common prefix from ref and alt
    ref_stripped, alt_stripped = _strip_common_prefix(ref.upper(), alt.upper())

    # Adjust position: VCF position is at the start of the padded allele;
    # after stripping the common prefix, the position shifts forward
    prefix_len = len(ref) - len(ref_stripped)
    vcf_pos = pos + prefix_len

    ref_len = len(ref_stripped)
    alt_len = len(alt_stripped)

    # Determine variant type and compute HGVS coordinates
    if ref_len > 0 and alt_len > 0:
        # Delins or substitution
        if ref_len == 1 and alt_len == 1:
            hgvs = f"{hgvs_acc}:{prefix}{vcf_pos}{ref_stripped}>{alt_stripped}"
        else:
            end_pos = vcf_pos + ref_len - 1
            hgvs = f"{hgvs_acc}:{prefix}{vcf_pos}_{end_pos}delins{alt_stripped}"

    elif ref_len > 0 and alt_len == 0:
        # Deletion
        if ref_len == 1:
            hgvs = f"{hgvs_acc}:{prefix}{vcf_pos}del{ref_stripped}"
        else:
            end_pos = vcf_pos + ref_len - 1
            hgvs = f"{hgvs_acc}:{prefix}{vcf_pos}_{end_pos}del{ref_stripped}"

    elif ref_len == 0 and alt_len > 0:
        # Insertion
        ins_pos = vcf_pos - 1
        hgvs = f"{hgvs_acc}:{prefix}{ins_pos}_{ins_pos + 1}ins{alt_stripped}"

    else:
        raise ValueError(
            f"Cannot determine variant type from ref={ref!r} alt={alt!r}"
        )

    return hgvs


def _strip_common_prefix(ref: str, alt: str) -> tuple:
    """Strip the longest common prefix from *ref* and *alt*.

    >>> _strip_common_prefix("ATCG", "AT")
    ('CG', '')
    >>> _strip_common_prefix("A", "ATCG")
    ('', 'TCG')
    >>> _strip_common_prefix("ATCG", "ATGG")
    ('CG', 'GG')
    >>> _strip_common_prefix("A", "G")
    ('A', 'G')
    """
    i = 0
    min_len = min(len(ref), len(alt))
    while i < min_len and ref[i] == alt[i]:
        i += 1
    return ref[i:], alt[i:]


# ---------------------------------------------------------------------------
# HGVS ↔ SPDI
# ---------------------------------------------------------------------------

def hgvs_to_spdi(hgvs_str: str) -> str:
    """Convert an HGVS variant string to SPDI format.

    SPDI (Sequence Position Deletion Insertion) format is:
    ``ACCESSION:POSITION:REFERENCE:ALTERNATE``
    where POSITION is 0-based.

    Parameters
    ----------
    hgvs_str : str
        HGVS variant string.

    Returns
    -------
    str
        SPDI-formatted string.

    Raises
    ------
    ValueError
        If the variant type has no explicit ref/alt (e.g. length-based
        deletions without specified bases).

    Examples
    --------
    Substitution:
    >>> hgvs_to_spdi("NC_000006.12:g.123456A>G")
    'NC_000006.12:123455:A:G'

    Deletion:
    >>> hgvs_to_spdi("NC_000006.12:g.123456delA")
    'NC_000006.12:123455:A:'

    Insertion:
    >>> hgvs_to_spdi("NC_000006.12:g.123456_123457insACGT")
    'NC_000006.12:123456::ACGT'

    Delins:
    >>> hgvs_to_spdi("NC_000006.12:g.123456_123457delinsAC")
    'NC_000006.12:123455::AC'
    """
    tag = parse(hgvs_str)

    if tag.prefix not in ('g.', 'c.', 'n.', 'm.'):
        raise ValueError(
            f"SPDI conversion requires nucleotide prefix (g./c./n./m.), "
            f"got {tag.prefix}"
        )

    acc = tag.genomic_acc or tag.acc
    vt = tag.variant_type

    if vt == 'substitution':
        if tag.ref is None or tag.alt is None:
            raise ValueError(f"Substitution missing ref/alt: {hgvs_str}")
        return f"{acc}:{tag.start_pos - 1}:{tag.ref}:{tag.alt}"

    if vt == 'deletion':
        deleted = tag.ref if tag.ref else ''
        if not deleted and tag.end_pos:
            # Range deletion without explicit bases — cannot determine REF
            raise ValueError(
                f"Deletion without explicit bases cannot be converted to SPDI. "
                f"Use a bases-specified HGVS like 'delACGT' instead of 'del': "
                f"{hgvs_str}"
            )
        return f"{acc}:{tag.start_pos - 1}:{deleted}:"

    if vt == 'insertion':
        if tag.alt is None:
            raise ValueError(f"Insertion missing alt: {hgvs_str}")
        # Insertion between start_pos and end_pos (1-based).
        # 0-based interbase position is start_pos (the position after which
        # the insertion occurs).
        return f"{acc}:{tag.start_pos}::{tag.alt}"

    if vt == 'delins':
        deleted = tag.ref if tag.ref else ''
        inserted = tag.alt if tag.alt else ''
        return f"{acc}:{tag.start_pos - 1}:{deleted}:{inserted}"

    if vt == 'duplication':
        dup = tag.alt if tag.alt else ''
        return f"{acc}:{tag.start_pos - 1}:{dup}:{dup}{dup}" if dup else \
               f"{acc}:{tag.start_pos - 1}::"

    raise ValueError(
        f"Variant type '{vt}' cannot be converted to SPDI: {hgvs_str}"
    )


def spdi_to_hgvs(spdi_str: str, prefix: str = "g.") -> str:
    """Convert an SPDI string to an HGVS variant description.

    SPDI format: ``ACC:POS:REF:ALT`` where POS is 0-based.

    Parameters
    ----------
    spdi_str : str
        SPDI-formatted string (e.g. ``"NC_000006.12:123455:A:G"``).
    prefix : str
        HGVS prefix to use (default ``"g."``).

    Returns
    -------
    str
        HGVS variant description.

    Raises
    ------
    ValueError
        If the SPDI string cannot be parsed.

    Examples
    --------
    Substitution:
    >>> spdi_to_hgvs("NC_000006.12:123455:A:G")
    'NC_000006.12:g.123456A>G'

    Deletion:
    >>> spdi_to_hgvs("NC_000006.12:123455:A:")
    'NC_000006.12:g.123456delA'

    Insertion:
    >>> spdi_to_hgvs("NC_000006.12:123455::ACGT")
    'NC_000006.12:g.123455_123456insACGT'

    Delins:
    >>> spdi_to_hgvs("NC_000006.12:123455:AT:GC")
    'NC_000006.12:g.123456_123457delinsGC'
    """
    parts = spdi_str.strip().split(':')
    if len(parts) != 4:
        raise ValueError(
            f"SPDI string must have 4 colon-separated fields, "
            f"got {len(parts)}: {spdi_str!r}"
        )

    acc, pos_str, ref, alt = parts
    pos = int(pos_str)

    # SPDI position is 0-based; HGVS is 1-based
    hgvs_pos = pos + 1

    if ref and alt:
        # Substitution or delins
        if len(ref) == 1 and len(alt) == 1:
            return f"{acc}:{prefix}{hgvs_pos}{ref}>{alt}"
        else:
            end_pos = hgvs_pos + len(ref) - 1
            return f"{acc}:{prefix}{hgvs_pos}_{end_pos}delins{alt}"

    elif ref and not alt:
        # Deletion
        if len(ref) == 1:
            return f"{acc}:{prefix}{hgvs_pos}del{ref}"
        else:
            end_pos = hgvs_pos + len(ref) - 1
            return f"{acc}:{prefix}{hgvs_pos}_{end_pos}del{ref}"

    elif not ref and alt:
        # Insertion
        return f"{acc}:{prefix}{pos}_{hgvs_pos}ins{alt}"

    else:
        raise ValueError(
            f"Both REF and ALT are empty in SPDI string. "
            f"SPDI requires at least one of REF or ALT to be non-empty: "
            f"{spdi_str!r}"
        )


# ---------------------------------------------------------------------------
# Coding → Protein
# ---------------------------------------------------------------------------

def c_to_p(hgvs_str: str, gene: str) -> str:
    """Convert a coding HGVS to a protein HGVS string.

    Fetches the MANE nucleotide and protein sequences for *gene*,
    applies the coding variant to the CDS, translates the mutated
    sequence, and returns the p.HGVS description.

    Parameters
    ----------
    hgvs_str : str
        Coding HGVS string (e.g. ``"NM_000207.3:c.1A>G"``).
    gene : str
        Gene symbol (e.g. ``"INS"``).

    Returns
    -------
    str
        Protein HGVS string in one-letter notation
        (e.g. ``"NP_000198.1:p.(M1?)"``).

    Raises
    ------
    ValueError
        If the variant is not a coding variant or cannot be applied.

    Examples
    --------
    >>> result = c_to_p("NM_000207.3:c.1A>G", "INS")  # doctest: +SKIP
    >>> ":p.(" in result  # doctest: +SKIP
    True
    """
    logger.info("Converting c.→p.: %s for %s", hgvs_str, gene)
    tag = parse(hgvs_str)

    if tag.prefix not in ('c.', 'n.'):
        raise ValueError(
            f"c_to_p requires a c. or n. prefix, got {tag.prefix}"
        )

    if tag.start_offset is not None:
        raise ValueError(
            f"Intronic variants cannot be directly translated: {hgvs_str}"
        )

    # Fetch sequences
    from . import _fetch
    nm_record = _fetch.nm(gene)
    np_record = _fetch.np(gene)

    # Extract CDS
    cds_seq = None
    for feature in nm_record.features:
        if feature.type == "CDS":
            cds_seq = feature.extract(nm_record).seq
            break
    if cds_seq is None:
        raise ValueError(f"No CDS feature found for {gene}")

    # Get the reference protein sequence
    ref_protein = str(np_record.seq)
    protein_id = np_record.id

    # Apply variant to CDS and translate
    mutated_cds = _apply_coding_variant(str(cds_seq), tag)
    mutated_protein = str(mutated_cds.translate(to_stop=True))

    # Find the first differing amino acid position
    aa_pos, ref_aa, alt_aa = _find_protein_change(ref_protein, mutated_protein)

    if aa_pos is None:
        # Silent variant — no protein change
        codon_num = (tag.start_pos - 1) // 3 + 1
        orig_aa = ref_protein[codon_num - 1] if codon_num <= len(ref_protein) else '?'
        return f"{protein_id}:p.({orig_aa}{codon_num}=)"

    # Build protein HGVS
    codon_num = aa_pos  # 1-based amino acid position

    if codon_num == 1:
        return f"{protein_id}:p.(M1?)"

    if alt_aa == '*':
        # Nonsense / stop-gain
        return f"{protein_id}:p.({ref_aa}{codon_num}*)"

    if ref_aa == '*':
        # Stop-loss / extension
        remaining = mutated_protein[aa_pos - 1:]
        return f"{protein_id}:p.({ref_aa}{codon_num}{alt_aa}ext*{len(remaining)})"

    return f"{protein_id}:p.({ref_aa}{codon_num}{alt_aa})"


def _apply_coding_variant(cds: str, tag) -> str:
    """Apply a coding variant parsed from an HGVSTag to a CDS string.

    Returns the mutated CDS as a string.
    """
    start = tag.start_pos - 1  # Convert to 0-based
    end = tag.end_pos - 1 if tag.end_pos else start

    vt = tag.variant_type

    if vt == 'substitution':
        if tag.ref is None or tag.alt is None:
            raise ValueError("Substitution missing ref/alt")
        # Verify ref matches
        ref_len = len(tag.ref)
        if cds[start:start + ref_len].upper() != tag.ref.upper():
            raise ValueError(
                f"Reference mismatch at position {tag.start_pos}: "
                f"expected {tag.ref}, found {cds[start:start + ref_len]}"
            )
        return cds[:start] + tag.alt + cds[start + ref_len:]

    if vt == 'deletion':
        del_len = (end - start + 1) if tag.end_pos else 1
        if tag.ref:
            del_len = len(tag.ref)
            if cds[start:start + del_len].upper() != tag.ref.upper():
                raise ValueError(
                    f"Reference mismatch in deletion at {tag.start_pos}"
                )
        return cds[:start] + cds[start + del_len:]

    if vt == 'insertion':
        if tag.alt is None:
            raise ValueError("Insertion missing alt bases")
        return cds[:start + 1] + tag.alt + cds[start + 1:]

    if vt == 'delins':
        del_len = (end - start + 1) if tag.end_pos else len(tag.ref) if tag.ref else 0
        ins_bases = tag.alt if tag.alt else ''
        return cds[:start] + ins_bases + cds[start + del_len:]

    if vt == 'duplication':
        dup_len = (end - start + 1) if tag.end_pos else 1
        dup_seq = cds[start:start + dup_len]
        return cds[:start + dup_len] + dup_seq + cds[start + dup_len:]

    raise ValueError(f"Cannot apply variant type '{vt}' to CDS")


def _find_protein_change(ref_protein: str, alt_protein: str) -> tuple:
    """Compare two protein sequences and return (1-based_pos, ref_aa, alt_aa).

    Returns ``(None, None, None)`` if identical (silent variant).
    """
    min_len = min(len(ref_protein), len(alt_protein))
    for i in range(min_len):
        if ref_protein[i] != alt_protein[i]:
            return (i + 1, ref_protein[i], alt_protein[i])

    # If one is longer, the change is at the first extra position
    if len(ref_protein) != len(alt_protein):
        pos = min_len + 1
        ref_aa = ref_protein[min_len] if min_len < len(ref_protein) else '*'
        alt_aa = alt_protein[min_len] if min_len < len(alt_protein) else '*'
        return (pos, ref_aa, alt_aa)

    # Identical
    return (None, None, None)


# ---------------------------------------------------------------------------
# Protein → Coding
# ---------------------------------------------------------------------------

def p_to_c_range(p_hgvs: str, gene: str) -> list:
    """Return all possible c.HGVS strings for a given protein HGVS.

    Delegates to the :mod:`backtranslate` module to enumerate every
    coding variant that could produce the observed amino acid change.

    Parameters
    ----------
    p_hgvs : str
        Protein HGVS string (e.g. ``"NP_000198.1:p.(V42G)"``).
    gene : str
        Gene symbol (e.g. ``"INS"``).

    Returns
    -------
    list of str
        All possible c.HGVS strings.

    Examples
    --------
    >>> result = p_to_c_range("NP_001347945.1:p.(V42G)", "G6PD")  # doctest: +SKIP
    >>> isinstance(result, list)  # doctest: +SKIP
    True
    >>> len(result) > 0  # doctest: +SKIP
    True
    >>> all(":c." in v for v in result)  # doctest: +SKIP
    True
    """
    return backtranslate.backtranslate(gene, p_hgvs)
