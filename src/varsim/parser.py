"""HGVS variant description parser.

Parses HGVS variant description strings into structured NamedTuple objects.
Supports all standard HGVS prefixes: c. (coding), p. (protein), g. (genomic),
n. (non-coding), m. (mitochondrial).

Handles variant types: substitution (>), deletion (del), insertion (ins),
duplication (dup), deletion-insertion (delins), inversion (inv), frameshift (fs),
repeat, and extension (*).
"""

import re
from typing import NamedTuple, Optional, Tuple

from ._logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Parsed variant representation
# ---------------------------------------------------------------------------

class HGVSTag(NamedTuple):
    """Structured representation of a parsed HGVS variant description.

    Attributes
    ----------
    acc : str
        Primary reference accession (e.g. ``"NM_000207.3"``).
    genomic_acc : str or None
        Genomic accession if present in NC_ACC(NM_ACC):c.X format, else None.
    prefix : str
        HGVS prefix: ``"c."``, ``"p."``, ``"g."``, ``"n."``, ``"m."``.
    start_pos : int
        Start coordinate (1-based). For ranges, the lower bound.
    end_pos : int or None
        End coordinate for range variants, else None.
    start_offset : int or None
        Intronic offset for splice-site variants (e.g. +1, -2).
    end_offset : int or None
        Intronic offset for the end of splice-site ranges.
    ref : str or None
        Reference allele string (e.g. ``"A"`` for ``A>G``).
    alt : str or None
        Alternate allele string (e.g. ``"G"`` for ``A>G``).
    variant_type : str
        One of: ``"substitution"``, ``"deletion"``, ``"insertion"``,
        ``"delins"``, ``"duplication"``, ``"inversion"``, ``"frameshift"``,
        ``"repeat"``, ``"extension"``, ``"uncertain"``.
    is_uncertain : bool
        True if the variant uses ``?`` notation (e.g. ``p.(M1?)``).
    fs_length : int or None
        For frameshift variants, the number of amino acids until stop
        (e.g. ``17`` from ``fs*17``).
    original : str
        The original unparsed HGVS string.
    """
    acc: str
    genomic_acc: Optional[str]
    prefix: str
    start_pos: int
    end_pos: Optional[int]
    start_offset: Optional[int]
    end_offset: Optional[int]
    ref: Optional[str]
    alt: Optional[str]
    variant_type: str
    is_uncertain: bool
    fs_length: Optional[int]
    original: str


# ---------------------------------------------------------------------------
# Regex patterns
# ---------------------------------------------------------------------------

# Matches HGVS prefix: "NM_000207.3:c." or "NC_000023.11(NM_000207.3):c."
_ACC_RE = re.compile(
    r'^'
    r'(?:'
    r'(?P<genomic_acc>NC_\d+\.\d+)\((?P<acc_wrapped>[A-Z]{2}_\d+\.\d+)\)'
    r'|'
    r'(?P<acc_simple>[A-Z]{2}_\d+\.\d+)'
    r')'
    r':(?P<prefix>[cngmpr]\.)'
)

# Matches the coordinate part: e.g. "123", "123+1", "123-2", "*1", "-59"
_COORD_RE = re.compile(
    r'(?P<type>[*]|[-])?'           # UTR marker (* or -)
    r'(?P<start>\d+)'               # Start coordinate
    r'(?:'                          # Optional intronic offset
    r'(?P<start_offset>[+-]\d+)'
    r')?'
    r'(?:_'                          # Optional range end
    r'(?P<end>\d+)'
    r'(?:(?P<end_offset>[+-]\d+))?'
    r')?'
)

# Variant-type specific patterns
_SUBST_RE = re.compile(r'(?P<ref>[ACGTUacgtu]+)>(?P<alt>[ACGTUacgtu]+)$')
_DEL_RE = re.compile(r'del(?P<length>\d*)(?P<bases>[ACGTUacgtu]*)$')
_INS_RE = re.compile(r'ins(?P<bases>[ACGTUacgtu]+)$')
_DELINS_RE = re.compile(r'delins(?P<bases>[ACGTUacgtu]+)$')
_DUP_RE = re.compile(r'dup(?P<bases>[ACGTUacgtu]*)$')
_INV_RE = re.compile(r'inv(?P<bases>[ACGTUacgtu]*)$')
_FS_RE = re.compile(r'(?P<ref_aa>[A-Z*])(?P<pos>\d+)(?P<alt_aa>[A-Z*]?)fs[*]?(?P<fs_len>\d*)$')
_PROTEIN_SUB_RE = re.compile(r'\((?P<ref_aa>[A-Z*])(?P<pos>\d+)(?P<alt_aa>[A-Z*=?])\)$')
_PROTEIN_EXT_RE = re.compile(r'\((?P<ref_aa>[A-Z*])(?P<pos>\d+)(?P<alt_aa>[A-Z*]?)ext[*]?(?P<ext_len>\d*)\)$')


def parse(hgvs_str: str) -> HGVSTag:
    """Parse an HGVS variant description string into a structured tag.

    Parameters
    ----------
    hgvs_str : str
        An HGVS variant description (e.g. ``"NM_000207.3:c.1A>G"``).

    Returns
    -------
    HGVSTag
        Structured representation of the parsed variant.

    Raises
    ------
    ValueError
        If the string cannot be parsed.

    Examples
    --------
    >>> tag = parse("NM_000207.3:c.1A>G")
    >>> tag.acc
    'NM_000207.3'
    >>> tag.prefix
    'c.'
    >>> tag.variant_type
    'substitution'
    >>> tag.start_pos
    1
    >>> tag.ref
    'A'
    >>> tag.alt
    'G'

    >>> tag = parse("NC_000023.11(NM_001360016.2):c.1A>G")
    >>> tag.genomic_acc
    'NC_000023.11'
    >>> tag.acc
    'NM_001360016.2'

    >>> tag = parse("NP_000198.1:p.(M1?)")
    >>> tag.prefix
    'p.'
    >>> tag.is_uncertain
    True

    >>> tag = parse("NM_000207.3:c.123_456del")
    >>> tag.variant_type
    'deletion'
    >>> tag.start_pos
    123
    >>> tag.end_pos
    456

    >>> tag = parse("NM_000207.3:c.123_124insACGT")
    >>> tag.variant_type
    'insertion'
    >>> tag.alt
    'ACGT'

    >>> tag = parse("NM_000207.3:c.123dup")
    >>> tag.variant_type
    'duplication'
    """
    original = hgvs_str.strip()
    logger.debug("Parsing HGVS: %s", hgvs_str)

    # --- Step 1: Extract accession and prefix ---
    m = _ACC_RE.match(original)
    if not m:
        raise ValueError(f"Cannot parse accession/prefix from: {original}")
    genomic_acc = m.group('genomic_acc') or None
    acc = m.group('acc_wrapped') or m.group('acc_simple')
    prefix = m.group('prefix')
    prefix_char = prefix[0]  # c, p, g, n, m
    rest = original[m.end():]

    # --- Step 3: Parse coordinate ---
    coord_m = _COORD_RE.search(rest)
    if not coord_m:
        raise ValueError(f"Cannot parse coordinates from: {rest}")
    coord_type = coord_m.group('type')  # *, -, or None
    start_raw = coord_m.group('start')
    start_offset_raw = coord_m.group('start_offset')
    end_raw = coord_m.group('end')
    end_offset_raw = coord_m.group('end_offset')

    # Compute actual start position
    start_pos = int(start_raw)
    if coord_type == '*':
        pass  # 3' UTR, keep as-is (positive, marked by prefix)
    elif coord_type == '-':
        start_pos = -start_pos  # 5' UTR

    start_offset = int(start_offset_raw) if start_offset_raw else None
    end_pos = int(end_raw) if end_raw else None
    end_offset = int(end_offset_raw) if end_offset_raw else None

    # Remove parsed coordinate from rest
    variant_part = rest[coord_m.end():]

    # --- Step 4: Determine variant type and parse ref/alt ---
    variant_type = 'substitution'  # default
    ref = None
    alt = None
    is_uncertain = False
    fs_length = None

    # Protein-specific parsing
    if prefix_char == 'p':
        is_uncertain = '?' in variant_part or '?' in original
        fs_m = _FS_RE.search(variant_part)
        ext_m = _PROTEIN_EXT_RE.search(variant_part)
        sub_m = _PROTEIN_SUB_RE.search(variant_part)
        if fs_m:
            variant_type = 'frameshift'
            ref = fs_m.group('ref_aa') or None
            alt = fs_m.group('alt_aa') or None
            fs_len_str = fs_m.group('fs_len')
            fs_length = int(fs_len_str) if fs_len_str else None
            # start_pos from coord_m is the amino acid position
            # But protein positions are typically embedded in the variant part
            # Handle case like p.(V42Gfs*17) where 42 is in p-sub
        elif ext_m:
            variant_type = 'extension'
            ref = ext_m.group('ref_aa') or None
            alt = ext_m.group('alt_aa') or None
        elif sub_m:
            ref = sub_m.group('ref_aa') or None
            alt = sub_m.group('alt_aa') or None
            if alt == '=':
                variant_type = 'substitution'  # silent
            elif alt == '?':
                is_uncertain = True
        elif is_uncertain:
            variant_type = 'uncertain'
        tag = HGVSTag(
            acc=acc, genomic_acc=genomic_acc, prefix=prefix,
            start_pos=start_pos, end_pos=end_pos,
            start_offset=start_offset, end_offset=end_offset,
            ref=ref, alt=alt, variant_type=variant_type,
            is_uncertain=is_uncertain, fs_length=fs_length, original=original,
        )
        logger.debug("Parsed → variant_type=%s pos=%d", tag.variant_type, tag.start_pos)
        return tag

    # Nucleotide-level parsing
    # Check variant types in priority order
    if _DELINS_RE.search(variant_part):
        variant_type = 'delins'
        m = _DELINS_RE.search(variant_part)
        alt = m.group('bases') if m else None
    elif _INV_RE.match(variant_part):
        variant_type = 'inversion'
        m = _INV_RE.match(variant_part)
        alt = m.group('bases') if m else None
    elif _DUP_RE.match(variant_part):
        variant_type = 'duplication'
        m = _DUP_RE.match(variant_part)
        alt = m.group('bases') if m else None
    elif _DEL_RE.match(variant_part):
        variant_type = 'deletion'
        m = _DEL_RE.match(variant_part)
        if m:
            bases = m.group('bases')
            if bases:
                ref = bases
            length = m.group('length')
            if length:
                deleted_len = int(length)
                # ref stays None for length-based del
    elif _INS_RE.search(variant_part):
        variant_type = 'insertion'
        m = _INS_RE.search(variant_part)
        alt = m.group('bases') if m else None
    elif _SUBST_RE.search(variant_part):
        variant_type = 'substitution'
        m = _SUBST_RE.search(variant_part)
        if m:
            ref = m.group('ref')
            alt = m.group('alt')
    else:
        # Unknown type but not protein
        if '?' in variant_part:
            is_uncertain = True
            variant_type = 'uncertain'

    tag = HGVSTag(
        acc=acc, genomic_acc=genomic_acc, prefix=prefix,
        start_pos=start_pos, end_pos=end_pos,
        start_offset=start_offset, end_offset=end_offset,
        ref=ref, alt=alt, variant_type=variant_type,
        is_uncertain=is_uncertain, fs_length=fs_length, original=original,
    )
    logger.debug("Parsed → variant_type=%s pos=%d", tag.variant_type, tag.start_pos)
    return tag
