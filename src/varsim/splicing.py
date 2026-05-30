"""Canonical splice site SNV simulation.

Generates all possible single-nucleotide variants at the canonical GT-AG
splice dinucleotides of each intron.
"""

from Bio.Data.IUPACData import unambiguous_dna_letters

from . import _fetch
from ._logging import get_logger

logger = get_logger(__name__)


def splice_site(gene: str) -> list:
    """Generate all possible splice site SNVs for a gene.

    Simulates substitutions at the canonical +1G, +2T (donor) and
    -2A, -1G (acceptor) positions of each intron within the CDS.

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. "G6PD").

    Returns
    -------
    list of str
        HGVS c. strings in the format
        ``NC_ACCESSION(NM_ACCESSION):c.COORDINATE±OFFSETREF>ALT``.

    Examples
    --------
    >>> result = splice_site("G6PD")
    >>> isinstance(result, list)
    True
    >>> len(result) > 0
    True
    >>> all(isinstance(v, str) for v in result)
    True
    >>> all("NC_" in v and "(NM_" in v for v in result)
    True
    >>> all(":c." in v for v in result)
    True
    >>> # Splice site variants use ± notation
    >>> any("+" in v or "-" in v for v in result)
    True
    """
    logger.info("Simulating splice site SNVs for %s ...", gene)
    variants = []
    seqrecord = _fetch.nm(gene)
    acc = _fetch.nc(gene)
    splicing = []
    for feature in seqrecord.features:
        if feature.type == "CDS":
            start = feature.location.start
            end = feature.location.end
    for feature in seqrecord.features:
        if feature.type == "exon":
            if feature.location.end < start or feature.location.start > end:
                continue
            else:
                splicing.extend(
                    (feature.location.start - start, feature.location.end - start)
                )

    for coordinate in range(1, len(splicing) - 1, 2):
        site = splicing[coordinate], splicing[coordinate] + 1
        for base in unambiguous_dna_letters:
            if base != "G":
                variants.append(f"{acc}({seqrecord.id}):c.{site[0]}+1G>{base}")
            if base != "T":
                variants.append(f"{acc}({seqrecord.id}):c.{site[0]}+2T>{base}")
            if base != "A":
                variants.append(f"{acc}({seqrecord.id}):c.{site[1]}-2A>{base}")
            if base != "G":
                variants.append(f"{acc}({seqrecord.id}):c.{site[1]}-1G>{base}")
    logger.info("splice_site(%s) → %d variants", gene, len(variants))
    return variants
