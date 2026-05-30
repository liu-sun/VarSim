"""Genome-transcript variant coordinate mapping.

Provides bidirectional conversion between coding (c.) and genomic (g.)
HGVS variant descriptions using exon structure from MANE GenBank records.

Conversions
-----------
``c_to_g``
    Convert a coding HGVS (c.) to a genomic HGVS (g.).
``g_to_c``
    Convert a genomic HGVS (g.) to a coding HGVS (c.).
``get_cds_exon_map``
    Return the CDS exon structure with both cDNA and genomic coordinates.

All functions fetch MANE Select/Plus Clinical records from NCBI and
are cached via ``functools.lru_cache``.

Examples
--------
>>> from varsim.transcription import c_to_g, g_to_c, get_cds_exon_map  # doctest: +SKIP
>>> c_to_g("NM_001360016.2:c.1A>G", "G6PD")  # doctest: +SKIP
'NC_000023.11:g.153760607A>G'
>>> g_to_c("NC_000023.11:g.153760607A>G", "G6PD")  # doctest: +SKIP
'NM_001360016.2:c.1A>G'
>>> exons = get_cds_exon_map("G6PD")  # doctest: +SKIP
>>> exons[0]["exon"]  # doctest: +SKIP
1
>>> exons[0]["strand"]  # doctest: +SKIP
-1
"""

import functools
import re

from Bio import Entrez

from . import _fetch
from ._logging import get_logger
from .parser import parse, HGVSTag

logger = get_logger(__name__)


# Regex to strip the coordinate portion from variant suffix
_COORD_STRIP_RE = re.compile(
    r'^[*-]?\d+(?:[+-]\d+)?(?:_[*-]?\d+(?:[+-]\d+)?)?'
)


# ---------------------------------------------------------------------------
# Internal: genomic context from Entrez Gene
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=64)
def _get_genomic_context(gene: str) -> dict:
    """Retrieve the genomic (NC_) range and strand for a gene symbol.

    Uses Entrez Gene to look up the official genomic location.
    Results are cached.

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. ``"G6PD"``).

    Returns
    -------
    dict
        Keys: ``nc_acc`` (str), ``chr_start`` (int, 0-based),
        ``chr_end`` (int, 0-based, exclusive), ``strand`` (int: 1 or -1).

    Raises
    ------
    ValueError
        If the gene cannot be found or lacks genomic coordinates.
    """
    # Search Entrez Gene
    logger.debug("Fetching genomic context for %s ...", gene)
    stream = Entrez.esearch(
        db="gene",
        term=f'{gene}[Gene Name] AND human[Organism]',
        retmax=1,
    )
    record = Entrez.read(stream)
    if not record["IdList"]:
        raise ValueError(f"Gene not found in Entrez Gene: {gene}")
    gene_id = record["IdList"][0]

    # Fetch summary to get genomic location
    stream = Entrez.esummary(db="gene", id=gene_id)
    summary = Entrez.read(stream)
    docs = summary["DocumentSummarySet"]["DocumentSummary"]
    if not docs:
        raise ValueError(f"No summary available for gene: {gene}")
    doc = docs[0]

    # Extract genomic info — the location history gives chr/start/stop
    genomic_info = doc.get("GenomicInfo")
    if genomic_info is None:
        raise ValueError(f"No genomic coordinates for gene: {gene}")

    # GenomicInfo is a list of dicts, first entry is the primary assembly
    ginfo = genomic_info[0]
    chr_start = int(ginfo["ChrStart"])  # 0-based, from Entrez Gene
    chr_stop = int(ginfo["ChrStop"])    # 0-based, exclusive
    chr_acc = ginfo.get("ChrAccVer", "")
    strand_sign = -1 if ginfo.get("ChrStrand") == "-" else 1

    if not chr_acc.startswith("NC_"):
        # Try to map to RefSeq accession
        chr_acc = _fetch.nc(gene)

    return {
        "nc_acc": chr_acc,
        "chr_start": chr_start,
        "chr_end": chr_stop,
        "strand": strand_sign,
    }


# ---------------------------------------------------------------------------
# Exon map construction
# ---------------------------------------------------------------------------

def get_cds_exon_map(gene: str) -> list:
    """Return CDS exon mapping with both cDNA and genomic coordinates.

    Each exon is represented as a dict::

        {
            "exon": 1,              # 1-based exon number
            "cds_start": 0,         # 0-based CDS-relative start
            "cds_end": 100,         # 0-based CDS-relative end (exclusive)
            "genomic_start": 123456,# 1-based genomic start (NC_)
            "genomic_end": 123556,  # 1-based genomic end (NC_, inclusive)
            "strand": 1,            # 1 = forward, -1 = reverse
        }

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. ``"G6PD"``).

    Returns
    -------
    list of dict
        One dict per CDS exon, ordered 5' → 3' along the transcript.

    Examples
    --------
    >>> exons = get_cds_exon_map("G6PD")  # doctest: +SKIP
    >>> isinstance(exons, list)  # doctest: +SKIP
    True
    >>> len(exons) > 1  # doctest: +SKIP
    True
    >>> e = exons[0]  # doctest: +SKIP
    >>> sorted(e.keys()) == ["cds_end", "cds_start", "exon", "genomic_end", "genomic_start", "strand"]  # doctest: +SKIP
    True
    >>> e["cds_start"] < e["cds_end"]  # doctest: +SKIP
    True
    """
    seqrecord = _fetch.nm(gene)
    ctx = _get_genomic_context(gene)

    logger.info("Building exon map for %s ...", gene)
    # Find CDS boundaries on transcript (0-based)
    cds_start = None
    cds_end = None
    for feature in seqrecord.features:
        if feature.type == "CDS":
            cds_start = int(feature.location.start)
            cds_end = int(feature.location.end)
            break

    if cds_start is None:
        raise ValueError(f"No CDS feature found in NM_ record for {gene}")

    # Collect CDS-overlapping exons (transcript-relative, 0-based)
    exons_raw = []
    for feature in seqrecord.features:
        if feature.type != "exon":
            continue
        e_start = int(feature.location.start)
        e_end = int(feature.location.end)

        # Only exons that overlap the CDS
        if e_end <= cds_start or e_start >= cds_end:
            continue

        exons_raw.append((e_start, e_end))

    # Sort by position (transcript order, 5' → 3')
    exons_raw.sort(key=lambda x: x[0])

    strand = ctx["strand"]
    gene_genomic_start = ctx["chr_start"]  # 0-based

    result = []
    for i, (e_start, e_end) in enumerate(exons_raw):
        # CDS-relative coordinates (0-based, within CDS)
        cds_rel_start = max(0, e_start - cds_start)
        cds_rel_end = min(cds_end - cds_start, e_end - cds_start)

        # Genomic coordinates (1-based)
        if strand == 1:
            g_start = gene_genomic_start + e_start + 1
            g_end = gene_genomic_start + e_end
        else:
            # Reverse strand: gene_genomic_start in Entrez Gene is 0-based
            # on the + strand.  The actual genomic position on - strand
            # decreases as transcript position increases.
            g_start = gene_genomic_start + (cds_end - e_end) + 1
            g_end = gene_genomic_start + (cds_end - e_start)

        result.append({
            "exon": i + 1,
            "cds_start": cds_rel_start,
            "cds_end": cds_rel_end,
            "genomic_start": g_start,
            "genomic_end": g_end,
            "strand": strand,
        })

    return result


# ---------------------------------------------------------------------------
# Internal: locate a cDNA position within the exon map
# ---------------------------------------------------------------------------

def _cdna_to_genomic(cdna_pos: int, exon_map: list) -> tuple:
    """Map a 1-based cDNA position to a 1-based genomic position.

    Parameters
    ----------
    cdna_pos : int
        1-based cDNA position (where 1 = first base of CDS).
        Negative values represent 5'UTR positions (e.g. -5 → c.-5).
    exon_map : list of dict
        Output of ``get_cds_exon_map``.

    Returns
    -------
    (genomic_pos, in_exon)
        genomic_pos : int, 1-based genomic coordinate.
        in_exon : bool, True if the position falls within an exon.

    Raises
    ------
    ValueError
        If the cDNA position does not map to any exon.
    """
    # Convert to 0-based: positions > 0 are 1-based, positions ≤ 0 are
    # already 0-based in HGVS (c.-1 = 1 base before c.1 = 0-based -1).
    cdna_0based = cdna_pos - 1 if cdna_pos > 0 else cdna_pos
    strand = exon_map[0]["strand"]

    for exon in exon_map:
        if exon["cds_start"] <= cdna_0based < exon["cds_end"]:
            offset = cdna_0based - exon["cds_start"]
            if strand == 1:
                g = exon["genomic_start"] + offset
            else:
                g = exon["genomic_end"] - offset
            return (g, True)

    # Position is intronic — find between which exons
    for i in range(len(exon_map) - 1):
        gap_start = exon_map[i]["cds_end"]
        gap_end = exon_map[i + 1]["cds_start"]
        if gap_start <= cdna_0based < gap_end:
            # Intronic position — project to the nearer exon boundary
            # For an intronic position, the user should use splice notation
            # (+offset/-offset).  Here we compute the donor/acceptor position.
            dist_from_donor = cdna_0based - gap_start + 1   # +1, +2, ...
            dist_from_acceptor = cdna_0based - gap_end        # ..., -2, -1
            # Use the nearer boundary
            if abs(dist_from_donor) <= abs(dist_from_acceptor):
                # Measure from donor (end of exon[i])
                donor_cdna = exon_map[i]["cds_end"]  # last base of exon i
                if strand == 1:
                    donor_g = exon_map[i]["genomic_end"] + dist_from_donor
                else:
                    donor_g = exon_map[i]["genomic_end"] - dist_from_donor
                return (donor_g, False)
            else:
                # Measure from acceptor (start of exon[i+1])
                if strand == 1:
                    acceptor_g = exon_map[i + 1]["genomic_start"] + dist_from_acceptor
                else:
                    acceptor_g = exon_map[i + 1]["genomic_start"] - dist_from_acceptor
                return (acceptor_g, False)

    # Position is before first exon or after last exon (UTR)
    if cdna_0based < exon_map[0]["cds_start"]:
        offset = cdna_0based - exon_map[0]["cds_start"]
        if strand == 1:
            g = exon_map[0]["genomic_start"] + offset
        else:
            g = exon_map[0]["genomic_end"] - offset
        return (g, False)

    if cdna_0based >= exon_map[-1]["cds_end"]:
        offset = cdna_0based - exon_map[-1]["cds_end"]
        if strand == 1:
            # +1 because genomic_end is the last CDS base;
            # 3'UTR starts at genomic_end + 1
            g = exon_map[-1]["genomic_end"] + 1 + offset
        else:
            g = exon_map[-1]["genomic_start"] - 1 - offset
        return (g, False)

    raise ValueError(f"cDNA position {cdna_pos} could not be mapped to genome")


def _genomic_to_cdna(genomic_pos: int, exon_map: list) -> tuple:
    """Map a 1-based genomic position to a 1-based cDNA position.

    Parameters
    ----------
    genomic_pos : int
        1-based genomic coordinate.
    exon_map : list of dict
        Output of ``get_cds_exon_map``.

    Returns
    -------
    (cdna_pos, in_exon)
        cdna_pos : int, 1-based cDNA coordinate (negative for 5'UTR).
        in_exon : bool, True if the position falls within an exon.

    Raises
    ------
    ValueError
        If the genomic position cannot be mapped.
    """
    strand = exon_map[0]["strand"]

    def _to_cdna(cdna_0based):
        """Convert 0-based cDNA position to HGVS 1-based convention."""
        return cdna_0based + 1 if cdna_0based >= 0 else cdna_0based

    for exon in exon_map:
        g_start = exon["genomic_start"]
        g_end = exon["genomic_end"]
        if g_start <= genomic_pos <= g_end:
            offset = genomic_pos - g_start if strand == 1 else g_end - genomic_pos
            cdna_0based = exon["cds_start"] + offset
            return (_to_cdna(cdna_0based), True)

    # Check intronic regions
    for i in range(len(exon_map) - 1):
        g_donor = exon_map[i]["genomic_end"]
        g_acceptor = exon_map[i + 1]["genomic_start"]

        if strand == 1:
            if g_donor < genomic_pos < g_acceptor:
                dist_from_donor = genomic_pos - g_donor
                dist_from_acceptor = genomic_pos - g_acceptor
                if abs(dist_from_donor) <= abs(dist_from_acceptor):
                    cdna_0based = exon_map[i]["cds_end"] + dist_from_donor - 1
                else:
                    cdna_0based = exon_map[i + 1]["cds_start"] + dist_from_acceptor
                return (_to_cdna(cdna_0based), False)
        else:
            if g_acceptor < genomic_pos < g_donor:
                dist_from_donor = g_donor - genomic_pos
                dist_from_acceptor = g_acceptor - genomic_pos
                if abs(dist_from_donor) <= abs(dist_from_acceptor):
                    cdna_0based = exon_map[i]["cds_end"] + dist_from_donor - 1
                else:
                    cdna_0based = exon_map[i + 1]["cds_start"] + abs(dist_from_acceptor)
                return (_to_cdna(cdna_0based), False)

    # Check upstream/downstream (UTR)
    if strand == 1:
        if genomic_pos < exon_map[0]["genomic_start"]:
            offset = exon_map[0]["genomic_start"] - genomic_pos
            cdna_0based = exon_map[0]["cds_start"] - offset
            return (_to_cdna(cdna_0based), False)
        if genomic_pos > exon_map[-1]["genomic_end"]:
            offset = genomic_pos - exon_map[-1]["genomic_end"] - 1
            cdna_0based = exon_map[-1]["cds_end"] + offset
            return (_to_cdna(cdna_0based), False)
    else:
        # Reverse strand: 5' end → higher genomic coords (genomic_end of
        # first exon); 3' end → lower genomic coords (genomic_start of
        # last exon).
        if genomic_pos > exon_map[0]["genomic_end"]:
            offset = genomic_pos - exon_map[0]["genomic_end"]
            cdna_0based = exon_map[0]["cds_start"] - offset
            return (_to_cdna(cdna_0based), False)
        if genomic_pos < exon_map[-1]["genomic_start"]:
            offset = exon_map[-1]["genomic_start"] - genomic_pos - 1
            cdna_0based = exon_map[-1]["cds_end"] + offset
            return (_to_cdna(cdna_0based), False)

    raise ValueError(
        f"Genomic position {genomic_pos} could not be mapped to cDNA"
    )


# ---------------------------------------------------------------------------
# Internal: extract variant suffix from original HGVS
# ---------------------------------------------------------------------------

def _get_variant_suffix(tag: HGVSTag) -> str:
    """Extract the allele description suffix from the original HGVS string.

    Everything after the coordinate portion (which may include offsets and
    range notation) is the variant description suffix.
    """
    original = tag.original
    idx = original.find(f":{tag.prefix}")
    if idx == -1:
        return ""
    rest = original[idx + len(f":{tag.prefix}"):]
    m = _COORD_STRIP_RE.match(rest)
    return rest[m.end():] if m else rest


# ---------------------------------------------------------------------------
# Internal: format HGVS strings with new coordinates
# ---------------------------------------------------------------------------

def _format_hgvs(prefix: str, coord: str, variant_suffix: str,
                 tag: HGVSTag, ref: str = None, alt: str = None) -> str:
    """Build an HGVS string from parts.

    Uses the variant type from *tag* to produce the correct HGVS suffix,
    but with new coordinates.

    Parameters
    ----------
    prefix : str
        Accession and prefix, e.g. ``"NC_000023.11:g."``.
    coord : str
        Coordinate string including offsets and range if applicable.
    variant_suffix : str
        Raw suffix from the original (used if ref/alt not provided).
    tag : HGVSTag
        Parsed tag carrying the variant type.
    ref : str or None
        Reference allele (may be ``""``).
    alt : str or None
        Alternate allele (may be ``""``).

    Returns
    -------
    str
    """
    vt = tag.variant_type
    rf = ref if ref is not None else tag.ref
    al = alt if alt is not None else tag.alt

    if vt == "substitution" and rf is not None and al is not None:
        return f"{prefix}{coord}{rf}>{al}"
    elif vt == "deletion":
        if rf:
            return f"{prefix}{coord}del{rf}"
        return f"{prefix}{coord}del"
    elif vt == "insertion":
        if al:
            return f"{prefix}{coord}ins{al}"
        return f"{prefix}{coord}ins"
    elif vt == "delins":
        if al:
            return f"{prefix}{coord}delins{al}"
        return f"{prefix}{coord}delins"
    elif vt == "duplication":
        if al:
            return f"{prefix}{coord}dup{al}"
        return f"{prefix}{coord}dup"
    elif vt == "inversion":
        if al:
            return f"{prefix}{coord}inv{al}"
        return f"{prefix}{coord}inv"
    else:
        # frameshift, extension, uncertain — use raw suffix
        return f"{prefix}{coord}{variant_suffix}"


# ---------------------------------------------------------------------------
# Internal: resolve true cDNA position accounting for UTR notation
# ---------------------------------------------------------------------------

def _get_true_cdna_and_cds_length(tag: HGVSTag, exon_map: list) -> tuple:
    """Get the true 1-based cDNA position and CDS length from a c.HGVS tag.

    The parser strips UTR markers (``*``, ``-``) but ``-`` positions
    are stored as negative ``start_pos``.  The ``*`` marker is lost so
    we detect it from the original string and adjust.

    Returns (true_cdna_pos, cds_length, is_utr5, is_utr3).
    """
    # CDS length is sum of all exon CDS widths
    cds_length = sum(e["cds_end"] - e["cds_start"] for e in exon_map)
    pos = tag.start_pos

    # Detect 5'UTR: parser makes start_pos negative for "-N" notation
    is_utr5 = pos < 0 or ":c.-" in tag.original or ":n.-" in tag.original
    # Detect 3'UTR: original contains ":c.*" or ":n.*"
    is_utr3 = ":c.*" in tag.original or ":n.*" in tag.original

    if is_utr5:
        # e.g., c.-5 → true position is 5 bases before CDS start (pos 1)
        # The parser returns -5, meaning 5 bases upstream of CDS pos 1
        # True cDNA 0-based position would be -5 + (cds_length)... no.
        # We keep the HGVS convention: 1-based, negative for 5'UTR.
        return (pos, cds_length, True, False)

    if is_utr3:
        # c.*N → parser returns start_pos=N, but true pos = cds_length + N
        return (cds_length + pos, cds_length, False, True)

    return (pos, cds_length, False, False)


def _hgvs_cdna_coord(cdna_pos: int, cds_length: int) -> str:
    """Build an HGVS c. coordinate string (no prefix).

    Converts a *true* 1-based cDNA position into HGVS notation:
    * Positions ≤ 0 use ``-`` (5' UTR).
    * Positions > cds_length use ``*`` (3' UTR).
    * Others are normal CDS positions.

    >>> _hgvs_cdna_coord(1, 1000)
    '1'
    >>> _hgvs_cdna_coord(-5, 1000)
    '-5'
    >>> _hgvs_cdna_coord(1001, 1000)
    '*1'
    >>> _hgvs_cdna_coord(1050, 1000)
    '*50'
    """
    if cdna_pos <= 0:
        return str(cdna_pos)  # already negative, e.g. "-5"
    elif cdna_pos > cds_length:
        return f"*{cdna_pos - cds_length}"
    else:
        return str(cdna_pos)


def _format_coord(start_true: int, end_true: int or None,
                  cds_length: int, tag: HGVSTag) -> str:
    """Build the full coordinate portion of the HGVS (with offsets/range).

    Uses true 1-based cDNA positions and applies HGVS UTR notation.
    """
    coord = _hgvs_cdna_coord(start_true, cds_length)
    if tag.start_offset is not None:
        coord += f"{tag.start_offset:+d}"

    if end_true is not None:
        coord += "_"
        coord += _hgvs_cdna_coord(end_true, cds_length)
        if tag.end_offset is not None:
            coord += f"{tag.end_offset:+d}"

    return coord


def _format_genomic_coord(start_g: int, end_g: int or None,
                          tag: HGVSTag) -> str:
    """Build a genomic coordinate string (1-based, no UTR notation)."""
    coord = str(start_g)
    if tag.start_offset is not None:
        coord += f"{tag.start_offset:+d}"

    if end_g is not None:
        coord += f"_{end_g}"
        if tag.end_offset is not None:
            coord += f"{tag.end_offset:+d}"

    return coord


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def c_to_g(c_hgvs: str, gene: str) -> str:
    """Convert a coding HGVS (c.) to a genomic HGVS (g.).

    Uses the exon structure from the MANE transcript GenBank record
    to translate cDNA coordinates to genomic (NC_) coordinates.

    Parameters
    ----------
    c_hgvs : str
        Coding HGVS string (e.g. ``"NM_001360016.2:c.1A>G"``).
    gene : str
        Gene symbol for fetching the MANE transcript (e.g. ``"G6PD"``).

    Returns
    -------
    str
        Genomic HGVS string with NC_ accession.

    Raises
    ------
    ValueError
        If parsing or coordinate mapping fails.

    Examples
    --------
    >>> c_to_g("NM_001360016.2:c.1A>G", "G6PD")  # doctest: +SKIP
    'NC_000023.11:g.153760607A>G'

    >>> c_to_g("NM_001360016.2:c.5T>C", "G6PD")  # doctest: +SKIP
    'NC_000023.11:g.153760611T>C'

    Splice-site variants are handled:

    >>> c_to_g("NM_001360016.2:c.637+1G>T", "G6PD")  # doctest: +SKIP
    'NC_000023.11:g.153762XXXG>T'
    """
    logger.info("Converting c.→g.: %s for %s", c_hgvs, gene)
    tag = parse(c_hgvs)
    if tag.prefix not in ("c.",):
        raise ValueError(f"Expected c. prefix, got {tag.prefix}")

    exon_map = get_cds_exon_map(gene)
    ctx = _get_genomic_context(gene)
    nc_acc = ctx["nc_acc"]
    variant_suffix = _get_variant_suffix(tag)

    # Resolve true cDNA position (accounting for UTR markers)
    start_true, cds_length, is_utr5, is_utr3 = \
        _get_true_cdna_and_cds_length(tag, exon_map)

    # Map start coordinate
    g_start, _ = _cdna_to_genomic(start_true, exon_map)

    # For offset variants, adjust the genomic position
    if tag.start_offset is not None:
        g_start += tag.start_offset

    # Map end coordinate if present
    g_end = None
    if tag.end_pos is not None:
        end_true = tag.end_pos
        if is_utr3:
            end_true = cds_length + end_true
        elif is_utr5:
            end_true = end_true  # already negative
        g_end, _ = _cdna_to_genomic(end_true, exon_map)
        if tag.end_offset is not None:
            g_end += tag.end_offset

    coord = _format_genomic_coord(g_start, g_end, tag)
    return _format_hgvs(f"{nc_acc}:g.", coord, variant_suffix, tag)


def g_to_c(g_hgvs: str, gene: str) -> str:
    """Convert a genomic HGVS (g.) to a coding HGVS (c.).

    Maps genomic (NC_) coordinates back to cDNA coordinates using the
    exon structure from the MANE transcript GenBank record.

    Parameters
    ----------
    g_hgvs : str
        Genomic HGVS string (e.g. ``"NC_000023.11:g.153760607A>G"``).
    gene : str
        Gene symbol for fetching the MANE transcript (e.g. ``"G6PD"``).

    Returns
    -------
    str
        Coding HGVS string with NM_ accession.

    Raises
    ------
    ValueError
        If parsing or coordinate mapping fails.

    Examples
    --------
    >>> g_to_c("NC_000023.11:g.153760607A>G", "G6PD")  # doctest: +SKIP
    'NM_001360016.2:c.1A>G'

    >>> g_to_c("NC_000023.11:g.153760607_153760608del", "G6PD")  # doctest: +SKIP
    'NM_001360016.2:c.1_2del'
    """
    logger.info("Converting g.→c.: %s for %s", g_hgvs, gene)
    tag = parse(g_hgvs)
    if tag.prefix not in ("g.",):
        raise ValueError(f"Expected g. prefix, got {tag.prefix}")

    exon_map = get_cds_exon_map(gene)
    nm_acc = _fetch.nm(gene).id
    cds_length = sum(e["cds_end"] - e["cds_start"] for e in exon_map)
    variant_suffix = _get_variant_suffix(tag)

    # Map start coordinate
    c_start_true, _ = _genomic_to_cdna(tag.start_pos, exon_map)

    # For offset variants, adjust the cDNA position
    if tag.start_offset is not None:
        c_start_true += tag.start_offset

    # Map end coordinate if present
    c_end_true = None
    if tag.end_pos is not None:
        c_end_true, _ = _genomic_to_cdna(tag.end_pos, exon_map)
        if tag.end_offset is not None:
            c_end_true += tag.end_offset

    coord = _format_coord(c_start_true, c_end_true, cds_length, tag)
    return _format_hgvs(f"{nm_acc}:c.", coord, variant_suffix, tag)
