"""Codon substitution simulation.

Generates all possible codon-to-codon nucleotide changes without protein
effect annotation. Uses the shared ``_codon_change`` helper for HGVS
c. string generation.
"""

from . import _fetch
from ._codon_change import _classify_codon_change
from ._logging import get_logger
from ._utils import genetic_code

logger = get_logger(__name__)


def codon_sub(gene: str) -> list:
    """Generate all possible codon substitution HGVS c. strings.

    For each codon, every alternative codon from the genetic code is tested.
    Only the nucleotide-level HGVS is returned (no protein annotation).

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. "INS").

    Returns
    -------
    list of str
        HGVS c. strings for each codon-level substitution.

    Examples
    --------
    >>> result = codon_sub("G6PD")
    >>> isinstance(result, list)
    True
    >>> len(result) > 1000
    True
    >>> isinstance(result[0], str)
    True
    >>> all(":c." in v for v in result[:10])
    True
    >>> # Should contain both > and delins variants
    >>> any(">" in v for v in result)
    True
    >>> any("delins" in v for v in result)
    True
    """
    logger.info("Simulating codon substitutions for %s ...", gene)
    variants = []
    seqrecord = _fetch.nm(gene)
    for feature in seqrecord.features:
        if feature.type == "CDS":
            cds_seq = feature.extract(seqrecord).seq
    for codon in range(0, len(cds_seq) - 3, 3):
        for new_codon in genetic_code:
            if new_codon != cds_seq[codon : codon + 3]:
                c_hgvs_suffix = _classify_codon_change(
                    str(cds_seq[codon : codon + 3]), new_codon, codon
                )
                variants.append(f"{seqrecord.id}:{c_hgvs_suffix}")
    logger.info("codon_sub(%s) → %d variants", gene, len(variants))
    return variants
