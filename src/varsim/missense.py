"""Missense variant simulation.

Generates all possible coding-sequence codon-level changes and annotates
their protein effect (missense, silent, or start-loss). Uses the shared
``_codon_change`` helper to classify codon substitutions.
"""

from Bio.Seq import Seq
from Bio.SeqUtils import seq3

from . import _fetch
from ._codon_change import _classify_codon_change
from ._logging import get_logger
from ._utils import genetic_code

logger = get_logger(__name__)


def missense(gene: str) -> list:
    """Generate all possible codon-level variants with protein effect.

    For each codon, every alternative codon from the genetic code is tested.
    Silent variants (no amino acid change) use ``p.(XxxN=)``; missense
    variants use ``p.(XxxNYyy)``; the start codon always uses ``p.(M1?)``.

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. "INS").

    Returns
    -------
    list of tuple
        Each tuple is ``(c_hgvs, p_hgvs_1letter, p_hgvs_3letter)``.

    Examples
    --------
    >>> result = missense("G6PD")
    >>> isinstance(result, list)
    True
    >>> len(result) > 500
    True
    >>> isinstance(result[0], tuple)
    True
    >>> len(result[0]) == 3
    True
    >>> c_hgvs, p_1, p_3 = result[0]
    >>> "NC_" in c_hgvs and "(NM_" in c_hgvs
    True
    >>> ":p.(" in p_1
    True
    >>> ":p.(" in p_3
    True
    >>> # Should have both silent (=) and missense variants
    >>> any("=" in v[1] for v in result)
    True
    >>> any("=" not in v[1] for v in result)
    True
    """
    logger.info("Simulating missense variants for %s ...", gene)
    variants = []
    seqrecord = _fetch.nm(gene)
    acc = _fetch.nc(gene)
    protein_seqrecord = _fetch.np(gene)
    protein = str(protein_seqrecord.seq)
    protein_id = protein_seqrecord.id
    for feature in seqrecord.features:
        if feature.type == "CDS":
            cds_seq = feature.extract(seqrecord).seq
    for index, codon in enumerate(range(0, len(cds_seq) - 3, 3)):
        for new_codon in genetic_code:
            if new_codon != cds_seq[codon : codon + 3]:
                seq = Seq(new_codon)
                c_hgvs_suffix = _classify_codon_change(
                    str(cds_seq[codon : codon + 3]), new_codon, codon
                )
                if index == 0:
                    variants.append(
                        (
                            f"{acc}({seqrecord.id}):{c_hgvs_suffix}",
                            f"{protein_id}:p.(M1?)",
                            f"{protein_id}:p.(Met1?)",
                        )
                    )
                else:
                    if protein[index] != seq.translate():
                        variants.append(
                            (
                                f"{acc}({seqrecord.id}):{c_hgvs_suffix}",
                                f"{protein_id}:p.({protein[index]}{index + 1}{seq.translate()})",
                                f"{protein_id}:p.({seq3(protein[index])}{index + 1}{seq3(seq.translate())})",
                            )
                        )
                    else:
                        variants.append(
                            (
                                f"{acc}({seqrecord.id}):{c_hgvs_suffix}",
                                f"{protein_id}:p.({protein[index]}{index + 1}=)",
                                f"{protein_id}:p.({seq3(protein[index])}{index + 1}=)",
                            )
                        )
    logger.info("missense(%s) → %d variants", gene, len(variants))
    return variants
