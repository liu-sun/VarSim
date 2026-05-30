"""HGVS variant extraction from reference and observed sequences.

Given a reference sequence and an observed/alternate sequence, this module
computes the minimal HGVS (Human Genome Variation Society) description of
the variant by aligning the sequences and identifying the difference.
"""

from typing import Optional

from ._logging import get_logger

logger = get_logger(__name__)


def extract(
    ref_seq: str,
    obs_seq: str,
    acc: str = "NM_000207.3",
    prefix: str = "c.",
) -> str:
    """Align *ref_seq* and *obs_seq* and generate the minimal HGVS description.

    Convenience wrapper around :func:`extract_with_positions` with
    ``start_pos=1`` and accession included.

    Parameters
    ----------
    ref_seq : str
        Reference sequence (DNA).
    obs_seq : str
        Observed / alternate sequence (DNA).
    acc : str, optional
        Reference sequence accession.  Defaults to ``"NM_000207.3"``.
    prefix : str, optional
        HGVS reference-type prefix (e.g. ``"c."``, ``"g."``, ``"n."``).

    Returns
    -------
    str
        HGVS variant description string.

    Examples
    --------
    >>> extract("ATGC", "ATTC", prefix="c.")
    'NM_000207.3:c.3G>T'

    >>> extract("ATGC", "ATC", prefix="c.")
    'NM_000207.3:c.3del'

    >>> extract("ATGC", "ATGGC", prefix="c.")
    'NM_000207.3:c.3_4insG'

    >>> extract("ATGC", "ACGT", prefix="c.")
    'NM_000207.3:c.2_4delinsCGT'
    """
    return extract_with_positions(ref_seq, obs_seq, start_pos=1, acc=acc, prefix=prefix)


def extract_with_positions(
    ref_seq: str,
    obs_seq: str,
    start_pos: int = 1,
    acc: Optional[str] = None,
    prefix: str = "c.",
) -> str:
    """Align *ref_seq* and *obs_seq* relative to *start_pos* and return the
    minimal HGVS description.

    The algorithm:

    1. Finds the first and last positions where the two sequences differ.
    2. Determines the variant type (substitution, deletion, insertion, delins).
    3. For pure deletions, trims identical flanking bases from both ends to
       produce the minimal coordinate range.
    4. Formats the HGVS string with the given *acc* (if any) and *prefix*.

    Parameters
    ----------
    ref_seq : str
        Reference sequence (DNA).  Position 0 of this string is treated as
        coordinate *start_pos* in 1-based HGVS notation.
    obs_seq : str
        Observed / alternate sequence (DNA).
    start_pos : int, optional
        1-based coordinate of the first nucleotide in *ref_seq*.
        Defaults to 1.
    acc : str or None, optional
        Reference sequence accession (e.g. ``"NM_000207.3"``).  When
        provided, it is prepended to the variant description
        (``ACC:PREFIX...``).  Defaults to *None*.
    prefix : str, optional
        HGVS reference-type prefix (e.g. ``"c."``, ``"g."``, ``"n."``).

    Returns
    -------
    str
        HGVS variant description string.

    Raises
    ------
    ValueError
        If *ref_seq* and *obs_seq* are identical (no variant to describe).

    Examples
    --------
    >>> extract_with_positions("ATGC", "ATTC", prefix="c.")
    'c.3G>T'

    >>> extract_with_positions("ATGC", "ATC", prefix="c.")
    'c.3del'

    >>> extract_with_positions("ATGC", "ATGGC", prefix="c.")
    'c.3_4insG'

    >>> extract_with_positions("ATGC", "ACGT", prefix="c.")
    'c.2_4delinsCGT'

    >>> # Using a custom start position
    >>> extract_with_positions("ATGC", "ATTC", start_pos=100, prefix="g.")
    'g.102G>T'

    >>> # With an accession
    >>> extract_with_positions("A", "G", acc="NM_000207.3", prefix="c.")
    'NM_000207.3:c.1A>G'

    >>> # Longer sequences
    >>> extract_with_positions("ATGCGTACG", "ATGCGTACGT", prefix="c.")
    'c.9_10insT'
    """
    if ref_seq == obs_seq:
        raise ValueError("ref_seq and obs_seq are identical; no variant to extract")

    logger.debug("Extracting HGVS from %d bp seqs", len(ref_seq))

    # ------------------------------------------------------------------
    # Step 1: find first and last differing positions
    # ------------------------------------------------------------------
    first_diff = 0
    min_len = min(len(ref_seq), len(obs_seq))
    while first_diff < min_len and ref_seq[first_diff] == obs_seq[first_diff]:
        first_diff += 1

    # Walk backwards from the ends to find the last differing position
    last_ref = len(ref_seq) - 1
    last_obs = len(obs_seq) - 1
    while (
        last_ref >= first_diff
        and last_obs >= first_diff
        and ref_seq[last_ref] == obs_seq[last_obs]
    ):
        last_ref -= 1
        last_obs -= 1

    # Extract the differing sub-sequences
    ref_diff = ref_seq[first_diff : last_ref + 1]
    obs_diff = obs_seq[first_diff : last_obs + 1]

    # Compute 1-based coordinates
    pos_start = start_pos + first_diff      # first affected ref position (1-based)
    pos_end = start_pos + last_ref           # last affected ref position (1-based, inclusive)

    # ------------------------------------------------------------------
    # Step 2: determine variant type
    # ------------------------------------------------------------------
    if len(ref_diff) == 0:
        # --------------------------------------------------------------
        # Pure insertion — nothing removed from ref, bases added to obs
        # --------------------------------------------------------------
        # Trim identical flanking bases from obs_diff to minimise
        _obs = obs_diff
        _ref = ref_diff
        # strip common prefix
        trim = 0
        while trim < min(len(_ref), len(_obs)) and _ref[trim] == _obs[trim]:
            trim += 1
        if trim:
            _ref = _ref[trim:]
            _obs = _obs[trim:]
            pos_start += trim
            # pos_end doesn't change for insertion (no ref bases affected)
        # strip common suffix
        trim_suf = 0
        while (
            trim_suf < min(len(_ref), len(_obs))
            and _ref[-1 - trim_suf] == _obs[-1 - trim_suf]
        ):
            trim_suf += 1
        if trim_suf:
            _ref = _ref[:-trim_suf] if trim_suf < len(_ref) else ""
            _obs = _obs[:-trim_suf] if trim_suf < len(_obs) else ""
        ins_before = pos_start - 1  # insertion sits between these two coordinates
        ins_after = pos_start
        hgvs = f"{prefix}{ins_before}_{ins_after}ins{_obs}"
    elif len(obs_diff) == 0:
        # --------------------------------------------------------------
        # Pure deletion — bases removed from ref, nothing added
        # --------------------------------------------------------------
        # Step 3: trim identical flanking bases to minimise
        _ref = ref_diff
        _obs = obs_diff
        # strip common prefix
        trim = 0
        while trim < min(len(_ref), len(_obs)) and _ref[trim] == _obs[trim]:
            trim += 1
        if trim:
            _ref = _ref[trim:]
            _obs = _obs[trim:]
            pos_start += trim
        # strip common suffix
        trim_suf = 0
        while (
            trim_suf < min(len(_ref), len(_obs))
            and _ref[-1 - trim_suf] == _obs[-1 - trim_suf]
        ):
            trim_suf += 1
        if trim_suf:
            _ref = _ref[:-trim_suf] if trim_suf < len(_ref) else ""
            _obs = _obs[:-trim_suf] if trim_suf < len(_obs) else ""
        pos_end = pos_start + len(_ref) - 1
        if len(_ref) == 1:
            hgvs = f"{prefix}{pos_start}del"
        else:
            hgvs = f"{prefix}{pos_start}_{pos_end}del"
    elif len(ref_diff) == 1 and len(obs_diff) == 1:
        # --------------------------------------------------------------
        # Single-nucleotide substitution
        # --------------------------------------------------------------
        hgvs = f"{prefix}{pos_start}{ref_diff}>{obs_diff}"
    else:
        # --------------------------------------------------------------
        # Deletion-insertion (delins) — the catch-all for complex changes
        # --------------------------------------------------------------
        pos_end = start_pos + last_ref
        if pos_start == pos_end:
            hgvs = f"{prefix}{pos_start}delins{obs_diff}"
        else:
            hgvs = f"{prefix}{pos_start}_{pos_end}delins{obs_diff}"

    # ------------------------------------------------------------------
    # Step 4: prepend accession if supplied
    # ------------------------------------------------------------------
    if acc:
        hgvs = f"{acc}:{hgvs}"

    logger.debug("Extracted → %s", hgvs)
    return hgvs
