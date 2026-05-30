"""Amino acid substitution simulation.

Generates all possible single amino acid changes for the MANE protein,
expressed in both one-letter and three-letter HGVS protein notation.
"""

from Bio.Data.IUPACData import protein_letters
from Bio.SeqUtils import seq3

from . import _fetch
from ._logging import get_logger

logger = get_logger(__name__)


def aa_sub(gene: str) -> list:
    """Generate all possible amino acid substitutions for a gene's protein.

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. "INS").

    Returns
    -------
    list of tuple
        Each tuple is ``(p_hgvs_1letter, p_hgvs_3letter)``. The first
        position (Met1) uses ``p.(M1?)`` / ``p.(Met1?)`` since start
        codon mutations have unpredictable effects.

    Examples
    --------
    >>> result = aa_sub("G6PD")
    >>> isinstance(result, list)
    True
    >>> len(result) > 1000
    True
    >>> isinstance(result[0], tuple)
    True
    >>> len(result[0]) == 2
    True
    >>> p_1, p_3 = result[0]
    >>> ":p.(" in p_1
    True
    >>> ":p.(" in p_3
    True
    >>> # Start codon uses M1? convention
    >>> any("M1?" in v[0] or "Met1?" in v[1] for v in result[:20])
    True
    """
    logger.info("Simulating amino acid substitutions for %s ...", gene)
    variants = []
    seqrecord = _fetch.np(gene)
    for index, residue in enumerate(seqrecord.seq, 1):
        for aa in protein_letters:
            if aa != residue:
                if index != 1:
                    variants.append(
                        (
                            f"{seqrecord.id}:p.({residue}{index}{aa})",
                            f"{seqrecord.id}:p.({seq3(residue)}{index}{seq3(aa)})",
                        )
                    )
                else:
                    variants.append(
                        (
                            f"{seqrecord.id}:p.({residue}{index}?)",
                            f"{seqrecord.id}:p.({seq3(residue)}{index}?)",
                        )
                    )
    logger.info("aa_sub(%s) → %d variants", gene, len(variants))
    return variants
