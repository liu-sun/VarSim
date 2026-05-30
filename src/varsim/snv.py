"""SNV simulation for CDS, 5'UTR, and 3'UTR regions.

Generates all possible single-nucleotide variants (SNVs) for the coding
sequence, 5' untranslated region, and 3' untranslated region of a gene's
MANE transcript.
"""

from Bio.Data.IUPACData import unambiguous_dna_letters
from Bio.Seq import Seq
from Bio.SeqFeature import SimpleLocation
from Bio.SeqUtils import seq3

from . import _fetch
from ._logging import get_logger

logger = get_logger(__name__)


def cds(gene: str) -> list:
    """Generate all possible CDS single-nucleotide variants.

    For each nucleotide position in the coding sequence, every possible
    alternative base is simulated. Each variant includes the nucleotide
    HGVS string and both one-letter and three-letter protein HGVS strings.

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. "INS", "G6PD").

    Returns
    -------
    list of tuple
        Each tuple is ``(c_hgvs, p_hgvs_1letter, p_hgvs_3letter)``.

    >>> result = cds("G6PD")
    >>> isinstance(result, list)
    True
    >>> len(result) > 1000
    True
    >>> c_hgvs, p_1, p_3 = result[0]
    >>> "NC_" in c_hgvs and "(NM_" in c_hgvs
    True
    >>> ":c." in c_hgvs
    True
    >>> ":p.(" in p_1
    True
    >>> ":p.(" in p_3
    True
    """
    logger.info("Simulating CDS SNVs for %s ...", gene)
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
        for base in unambiguous_dna_letters:
            if index == 0:
                if base != "A":
                    variants.append(
                        (
                            f"{acc}({seqrecord.id}):c.1A>{base}",
                            f"{protein_id}:p.(M1?)",
                            f"{protein_id}:p.(Met1?)",
                        )
                    )
            else:
                if base != cds_seq[codon]:
                    seq = Seq(base) + cds_seq[codon + 1 : codon + 3]
                    if protein[index] != seq.translate():
                        variants.append(
                            (
                                f"{acc}({seqrecord.id}):c.{codon + 1}{cds_seq[codon]}>{base}",
                                f"{protein_id}:p.({protein[index]}{index + 1}{seq.translate()})",
                                f"{protein_id}:p.({seq3(protein[index])}{index + 1}{seq3(seq.translate())})",
                            )
                        )
                    else:
                        variants.append(
                            (
                                f"{acc}({seqrecord.id}):c.{codon + 1}{cds_seq[codon]}>{base}",
                                f"{protein_id}:p.({protein[index]}{index + 1}=)",
                                f"{protein_id}:p.({seq3(protein[index])}{index + 1}=)",
                            )
                        )
            if index == 0:
                if base != "T":
                    variants.append(
                        (
                            f"{acc}({seqrecord.id}):c.2T>{base}",
                            f"{protein_id}:p.(M1?)",
                            f"{protein_id}:p.(Met1?)",
                        )
                    )
            else:
                if base != cds_seq[codon + 1]:
                    seq = cds_seq[codon] + Seq(base) + cds_seq[codon + 2]
                    if protein[index] != seq.translate():
                        variants.append(
                            (
                                f"{acc}({seqrecord.id}):c.{codon + 2}{cds_seq[codon + 1]}>{base}",
                                f"{protein_id}:p.({protein[index]}{index + 1}{seq.translate()})",
                                f"{protein_id}:p.({seq3(protein[index])}{index + 1}{seq3(seq.translate())})",
                            )
                        )
                    else:
                        variants.append(
                            (
                                f"{acc}({seqrecord.id}):c.{codon + 2}{cds_seq[codon + 1]}>{base}",
                                f"{protein_id}:p.({protein[index]}{index + 1}=)",
                                f"{protein_id}:p.({seq3(protein[index])}{index + 1}=)",
                            )
                        )
            if index == 0:
                if base != "G":
                    variants.append(
                        (
                            f"{acc}({seqrecord.id}):c.3G>{base}",
                            f"{protein_id}:p.(M1?)",
                            f"{protein_id}:p.(Met1?)",
                        )
                    )
            else:
                if base != cds_seq[codon + 2]:
                    seq = cds_seq[codon : codon + 2] + Seq(base)
                    if protein[index] != seq.translate():
                        variants.append(
                            (
                                f"{acc}({seqrecord.id}):c.{codon + 3}{cds_seq[codon + 2]}>{base}",
                                f"{protein_id}:p.({protein[index]}{index + 1}{seq.translate()})",
                                f"{protein_id}:p.({seq3(protein[index])}{index + 1}{seq3(seq.translate())})",
                            )
                        )
                    else:
                        variants.append(
                            (
                                f"{acc}({seqrecord.id}):c.{codon + 3}{cds_seq[codon + 2]}>{base}",
                                f"{protein_id}:p.({protein[index]}{index + 1}=)",
                                f"{protein_id}:p.({seq3(protein[index])}{index + 1}=)",
                            )
                        )
    logger.info("cds(%s) → %d variants", gene, len(variants))
    return variants


def utr5(gene: str) -> list:
    """Generate all possible 5'UTR single-nucleotide variants.

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. "INS").

    Returns
    -------
    list of str
        HGVS c. strings for each 5'UTR SNV.

    >>> result = utr5("G6PD")
    >>> isinstance(result, list)
    True
    >>> len(result) > 0
    True
    >>> isinstance(result[0], str)
    True
    >>> "NC_" in result[0] and "(NM_" in result[0]
    True
    >>> ":c.-" in result[0]
    True
    """
    logger.info("Simulating 5'UTR SNVs for %s ...", gene)
    variants = []
    seqrecord = _fetch.nm(gene)
    acc = _fetch.nc(gene)
    for feature in seqrecord.features:
        if feature.type == "CDS":
            utr5_seq = SimpleLocation(0, feature.location.start).extract(seqrecord).seq
    for index in range(len(utr5_seq)):
        for base in unambiguous_dna_letters:
            if base != utr5_seq[index]:
                variants.append(
                    f"{acc}({seqrecord.id}):c.{index - len(utr5_seq)}{utr5_seq[index]}>{base}"
                )
    logger.info("utr5(%s) → %d variants", gene, len(variants))
    return variants


def utr3(gene: str) -> list:
    """Generate all possible 3'UTR single-nucleotide variants.

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. "INS").

    Returns
    -------
    list of str
        HGVS c. strings for each 3'UTR SNV.

    >>> result = utr3("G6PD")
    >>> isinstance(result, list)
    True
    >>> len(result) > 0
    True
    >>> isinstance(result[0], str)
    True
    >>> "NC_" in result[0] and "(NM_" in result[0]
    True
    >>> ":c.*" in result[0]
    True
    """
    logger.info("Simulating 3'UTR SNVs for %s ...", gene)
    variants = []
    seqrecord = _fetch.nm(gene)
    acc = _fetch.nc(gene)
    for feature in seqrecord.features:
        if feature.type == "CDS":
            utr3_seq = (
                SimpleLocation(feature.location.end, len(seqrecord))
                .extract(seqrecord)
                .seq
            )
    for index in range(len(utr3_seq)):
        for base in unambiguous_dna_letters:
            if base != utr3_seq[index]:
                variants.append(f"{acc}({seqrecord.id}):c.*{index + 1}{utr3_seq[index]}>{base}")
    logger.info("utr3(%s) → %d variants", gene, len(variants))
    return variants
