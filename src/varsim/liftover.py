"""HGVS variant liftover between genome assemblies.

Lifts genomic HGVS (g.) variants between assemblies (e.g. GRCh37 ↔ GRCh38)
using the NCBI Remap API. Also provides convenience functions for
transcript-based liftover via MANE Select/Plus Clinical transcripts.

API Reference
-------------
* NCBI Remap: https://api.ncbi.nlm.nih.gov/variation/v0/remap/
"""

import json
import urllib.error
import urllib.request
from typing import Optional

from ._logging import get_logger
from .parser import HGVSTag, parse

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# SPDI ↔ HGVS helpers
# ---------------------------------------------------------------------------


def _build_spdi(acc: str, pos_1based: int, ref: Optional[str], alt: Optional[str]) -> str:
    """Convert genomic coordinates + alleles to SPDI format.

    SPDI: ``seq_id:position:deleted_sequence:inserted_sequence`` where
    *position* is 0-based interbase.

    Parameters
    ----------
    acc : str
        NC_ accession.
    pos_1based : int
        HGVS 1-based start coordinate.
    ref : str or None
        Reference allele (deleted sequence in SPDI).
    alt : str or None
        Alternate allele (inserted sequence in SPDI).

    Returns
    -------
    str
        SPDI-formatted string for the NCBI Remap API.
    """
    pos_0based = pos_1based - 1
    deleted = ref or ""
    inserted = alt or ""
    return f"{acc}:{pos_0based}:{deleted}:{inserted}"


def _hgvs_from_spdi(spdi: dict) -> str:
    """Convert a remap SPDI result dict to HGVS g. notation.

    Parameters
    ----------
    spdi : dict
        SPDI dict with keys ``seq_id``, ``position``, ``deleted_sequence``,
        ``inserted_sequence``.

    Returns
    -------
    str
        HGVS genomic variant string (e.g. ``"NC_000001.11:g.123A>G"``).
    """
    acc = spdi["seq_id"]
    pos = spdi["position"]  # 0-based interbase
    deleted = spdi.get("deleted_sequence", "")
    inserted = spdi.get("inserted_sequence", "")

    p1 = pos + 1  # 1-based

    if deleted and inserted:
        if len(deleted) == 1 and len(inserted) == 1:
            return f"{acc}:g.{p1}{deleted}>{inserted}"
        end = pos + len(deleted)
        return f"{acc}:g.{p1}_{end}delins{inserted}"
    elif deleted:
        if len(deleted) == 1:
            return f"{acc}:g.{p1}del{deleted}"
        end = pos + len(deleted)
        return f"{acc}:g.{p1}_{end}del{deleted}"
    elif inserted:
        return f"{acc}:g.{p1}_{p1 + 1}ins{inserted}"
    else:
        return f"{acc}:g.{p1}="


def _tag_to_spdi(tag: HGVSTag) -> str:
    """Convert a parsed genomic HGVSTag to SPDI format.

    Raises
    ------
    ValueError
        If the accession is not an NC_ genomic accession.
    """
    acc = tag.genomic_acc or tag.acc
    if not acc.startswith("NC_"):
        raise ValueError(f"Genomic NC_ accession required, got: {acc}")

    return _build_spdi(acc, tag.start_pos, tag.ref, tag.alt)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def liftover_g_to_assembly(
    hgvs_str: str, target_assembly: str = "GRCh38"
) -> str:
    """Lift a genomic HGVS (g.) variant to a different genome assembly.

    Calls the `NCBI Remap API
    <https://api.ncbi.nlm.nih.gov/variation/v0/remap/>`_ to convert
    coordinates between assemblies.

    Parameters
    ----------
    hgvs_str : str
        Genomic HGVS string (e.g. ``"NC_000001.10:g.12345A>G"``).
    target_assembly : str
        Target assembly name. Default ``"GRCh38"``.
        Also accepts ``"GRCh37"``.

    Returns
    -------
    str
        The remapped HGVS string.  If remapping fails the original
        string is returned and a warning is logged.

    Examples
    --------
    >>> liftover_g_to_assembly(
    ...     "NC_000001.10:g.12345A>G", "GRCh38"
    ... )  # doctest: +SKIP
    'NC_000001.11:g.12345A>G'
    """
    tag = parse(hgvs_str)

    if tag.prefix != "g.":
        logger.warning(
            "Expected genomic (g.) HGVS, got %s. Attempting anyway.",
            tag.prefix,
        )

    logger.info("Lifting over to %s: %s", target_assembly, hgvs_str)
    # Build SPDI
    try:
        spdi = _tag_to_spdi(tag)
    except ValueError:
        logger.warning(
            "Cannot convert %s to SPDI – returning original.", hgvs_str
        )
        return hgvs_str

    url = (
        "https://api.ncbi.nlm.nih.gov/variation/v0/remap/"
        + spdi
        + "?target_assembly="
        + target_assembly
    )

    # Call NCBI Remap
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read().decode())
    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError) as exc:
        logger.warning(
            "NCBI Remap API call failed for %s: %s. Returning original.",
            hgvs_str, exc,
        )
        return hgvs_str

    # --- parse response ---
    remap_list = data.get("remap", [])
    if isinstance(remap_list, dict):
        remap_list = [remap_list]

    if not remap_list:
        logger.warning("No remapping result for %s. Returning original.", hgvs_str)
        return hgvs_str

    result = remap_list[0]

    # Prefer the server-supplied HGVS
    hgvs_result = result.get("hgvs")
    if hgvs_result:
        return hgvs_result

    # Fall back to building HGVS from returned SPDI
    spdi_result = result.get("spdi")
    if spdi_result:
        return _hgvs_from_spdi(spdi_result)

    logger.warning("Unexpected remap response format. Returning original.")
    return hgvs_str


def liftover_grch37_to_grch38(hgvs_str: str) -> str:
    """Lift a genomic HGVS from GRCh37 to GRCh38.

    Convenience wrapper around :func:`liftover_g_to_assembly`.

    Parameters
    ----------
    hgvs_str : str
        Genomic HGVS on GRCh37 (e.g. ``"NC_000001.10:g.12345A>G"``).

    Returns
    -------
    str
        Remapped HGVS on GRCh38.

    Examples
    --------
    >>> liftover_grch37_to_grch38(
    ...     "NC_000001.10:g.12345A>G"
    ... )  # doctest: +SKIP
    'NC_000001.11:g.12345A>G'
    """
    return liftover_g_to_assembly(hgvs_str, target_assembly="GRCh38")


def liftover_grch38_to_grch37(hgvs_str: str) -> str:
    """Lift a genomic HGVS from GRCh38 to GRCh37.

    Convenience wrapper around :func:`liftover_g_to_assembly`.

    Parameters
    ----------
    hgvs_str : str
        Genomic HGVS on GRCh38 (e.g. ``"NC_000001.11:g.12345A>G"``).

    Returns
    -------
    str
        Remapped HGVS on GRCh37.

    Examples
    --------
    >>> liftover_grch38_to_grch37(
    ...     "NC_000001.11:g.12345A>G"
    ... )  # doctest: +SKIP
    'NC_000001.10:g.12345A>G'
    """
    return liftover_g_to_assembly(hgvs_str, target_assembly="GRCh37")


def liftover_transcript(
    gene: str, c_hgvs: str, target_assembly: str = "GRCh38"
) -> str:
    """Lift a coding HGVS to a target assembly via the MANE transcript.

    Fetches the MANE Select/Plus Clinical transcript for *gene*,
    converts the c.HGVS to genomic coordinates using the CDS start
    position, then lifts to *target_assembly* with
    :func:`liftover_g_to_assembly`.

    .. note::

       Simple exonic variants (substitutions, deletions, insertions
       without intronic offsets) are best supported.  Intronic /
       splice-site variants may not map correctly because this
       function uses a linear CDS‐start offset rather than the full
       exon structure.

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. ``"G6PD"``).
    c_hgvs : str
        Coding HGVS string (e.g. ``"c.1A>G"``).  The accession prefix
        is optional – ``"NM_000207.3:c.1A>G"`` and ``"c.1A>G"`` are
        both accepted.
    target_assembly : str
        Target assembly name.  Default ``"GRCh38"``.

    Returns
    -------
    str
        The remapped genomic HGVS string.

    Raises
    ------
    ValueError
        If no CDS feature is found on the MANE transcript.

    Examples
    --------
    >>> liftover_transcript("G6PD", "c.1A>G", "GRCh38")  # doctest: +SKIP
    'NC_000023.11:g.154535278T>C'
    """
    from . import _fetch  # deferred – avoids circular import at module level

    logger.info("Lifting over transcript %s for %s", c_hgvs, gene)
    tag = parse(c_hgvs)
    if tag.prefix != "c.":
        logger.warning(
            "Expected coding (c.) HGVS, got %s. Proceeding anyway.",
            tag.prefix,
        )

    if tag.start_offset is not None or tag.end_offset is not None:
        logger.warning(
            "Intronic offsets in c.HGVS may not map correctly "
            "with linear CDS offset."
        )

    # ---- fetch MANE transcript & genomic accession ----
    seqrecord = _fetch.nm(gene)
    genomic_acc = _fetch.nc(gene)

    # ---- extract CDS start (0-based in Biopython) ----
    cds_features = [f for f in seqrecord.features if f.type == "CDS"]
    if not cds_features:
        raise ValueError(f"No CDS feature found on MANE transcript for {gene}")
    cds = cds_features[0]
    cds_start = int(cds.location.start)  # 0-based

    # ---- c. (1-based) → g. (1-based) ----
    g_start = cds_start + tag.start_pos

    # ---- build genomic HGVS ----
    if tag.variant_type == "substitution" and tag.ref and tag.alt:
        g_hgvs = f"{genomic_acc}:g.{g_start}{tag.ref}>{tag.alt}"
    elif tag.variant_type == "deletion":
        if tag.ref:
            g_hgvs = f"{genomic_acc}:g.{g_start}del{tag.ref}"
        elif tag.end_pos is not None:
            g_end = cds_start + tag.end_pos
            g_hgvs = f"{genomic_acc}:g.{g_start}_{g_end}del"
        else:
            g_hgvs = f"{genomic_acc}:g.{g_start}del"
    elif tag.variant_type == "insertion" and tag.alt:
        g_hgvs = f"{genomic_acc}:g.{g_start}_{g_start + 1}ins{tag.alt}"
    elif tag.variant_type == "delins" and tag.alt:
        if tag.end_pos is not None:
            g_end = cds_start + tag.end_pos
            g_hgvs = f"{genomic_acc}:g.{g_start}_{g_end}delins{tag.alt}"
        elif tag.ref:
            g_hgvs = f"{genomic_acc}:g.{g_start}del{tag.ref}ins{tag.alt}"
        else:
            g_hgvs = f"{genomic_acc}:g.{g_start}delins{tag.alt}"
    elif tag.variant_type == "duplication":
        if tag.end_pos is not None:
            g_end = cds_start + tag.end_pos
            g_hgvs = f"{genomic_acc}:g.{g_start}_{g_end}dup"
        elif tag.ref:
            g_hgvs = f"{genomic_acc}:g.{g_start}dup{tag.ref}"
        else:
            g_hgvs = f"{genomic_acc}:g.{g_start}dup"
    else:
        logger.warning(
            "Cannot convert c.HGVS to g.HGVS for variant type %r",
            tag.variant_type,
        )
        return c_hgvs

    return liftover_g_to_assembly(g_hgvs, target_assembly=target_assembly)
