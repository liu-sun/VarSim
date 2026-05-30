"""Frameshift indel simulation.

Generates all possible frameshift-causing insertions and deletions across
the coding sequence. A frameshift is caused by an indel whose length is
not a multiple of 3.
"""

from Bio.Data.IUPACData import unambiguous_dna_letters
from Bio.Seq import Seq
from Bio.SeqUtils import seq3

from . import _fetch
from ._logging import get_logger

logger = get_logger(__name__)


def frameshift(gene: str) -> list:
    """Generate all possible frameshift indel variants.

    For each nucleotide position in the CDS:
    - 1bp deletion (causes frameshift)
    - 1bp insertion of each of the 4 nucleotides (A, T, G, C)

    For each variant, the frameshifted protein sequence is translated
    from the mutation point until the first in-frame stop codon.

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
    >>> result = frameshift("G6PD")  # doctest: +SKIP
    >>> isinstance(result, list)  # doctest: +SKIP
    True
    >>> len(result) > 1000  # doctest: +SKIP
    True
    >>> isinstance(result[0], tuple)  # doctest: +SKIP
    True
    >>> len(result[0]) == 3  # doctest: +SKIP
    True
    >>> c_hgvs, p_1, p_3 = result[0]  # doctest: +SKIP
    >>> ":c." in c_hgvs  # doctest: +SKIP
    True
    >>> "):c." in c_hgvs  # doctest: +SKIP
    True
    >>> ":p.(" in p_1  # doctest: +SKIP
    True
    >>> ":p.(" in p_3  # doctest: +SKIP
    True
        >>> # Should have both del and ins variants
    >>> any("del" in v[0] for v in result)  # doctest: +SKIP
    True
    >>> any("ins" in v[0] for v in result)  # doctest: +SKIP
    True
        >>> # Protein effect should include fs notation
    >>> any("fs" in v[1] for v in result)  # doctest: +SKIP
    True
    >>> any("fs" in v[2] for v in result)  # doctest: +SKIP
    True
    """
    logger.info("Simulating frameshift variants for %s ...", gene)
    variants = []
    seqrecord = _fetch.nm(gene)
    protein_seqrecord = _fetch.np(gene)
    protein_id = protein_seqrecord.id
    acc = _fetch.nc(gene)
    for feature in seqrecord.features:
        if feature.type == "CDS":
            cds_seq = feature.extract(seqrecord).seq

    # 1bp deletions — frameshift at every position
    for pos in range(len(cds_seq)):
        # Build the frameshifted CDS: remove base at pos
        mutated_cds = cds_seq[:pos] + cds_seq[pos + 1 :]
        # Translate from the codon containing the deletion
        codon_index = (pos // 3) * 3
        fs_seq = mutated_cds[codon_index:]
        # Translate to first stop
        aa_seq = fs_seq.translate(to_stop=True)
        orig_aa = str(cds_seq.translate())[pos // 3]

        # HGVS c. notation
        c_hgvs = f"{acc}({seqrecord.id}):c.{pos + 1}del"

        # HGVS p. notation
        codon_num = (pos // 3) + 1
        if codon_num == 1:
            p_1 = f"{protein_id}:p.(M1?)"
            p_3 = f"{protein_id}:p.(Met1?)"
        elif len(aa_seq) == 0:
            p_1 = f"{protein_id}:p.({orig_aa}{codon_num}?)"
            p_3 = f"{protein_id}:p.({seq3(orig_aa)}{codon_num}?)"
        else:
            p_1 = f"{protein_id}:p.({orig_aa}{codon_num}{aa_seq[0]}fs*{len(aa_seq)})"
            p_3 = f"{protein_id}:p.({seq3(orig_aa)}{codon_num}{seq3(aa_seq[0])}fs*{len(aa_seq)})"

        variants.append((c_hgvs, p_1, p_3))

    # 1bp insertions — frameshift after every position
    for pos in range(len(cds_seq)):
        for base in unambiguous_dna_letters:
            # Insert base after position pos
            mutated_cds = cds_seq[: pos + 1] + Seq(base) + cds_seq[pos + 1 :]
            codon_index = (pos // 3) * 3
            fs_seq = mutated_cds[codon_index:]
            aa_seq = fs_seq.translate(to_stop=True)
            orig_aa = str(cds_seq.translate())[pos // 3]

            # HGVS c. notation: insertion between pos and pos+1
            c_hgvs = f"{acc}({seqrecord.id}):c.{pos}_{pos + 1}ins{base}"

            codon_num = (pos // 3) + 1
            if codon_num == 1:
                p_1 = f"{protein_id}:p.(M1?)"
                p_3 = f"{protein_id}:p.(Met1?)"
            elif len(aa_seq) == 0:
                p_1 = f"{protein_id}:p.({orig_aa}{codon_num}?)"
                p_3 = f"{protein_id}:p.({seq3(orig_aa)}{codon_num}?)"
            else:
                p_1 = f"{protein_id}:p.({orig_aa}{codon_num}{aa_seq[0]}fs*{len(aa_seq)})"
                p_3 = f"{protein_id}:p.({seq3(orig_aa)}{codon_num}{seq3(aa_seq[0])}fs*{len(aa_seq)})"

            variants.append((c_hgvs, p_1, p_3))

    logger.info("frameshift(%s) → %d variants", gene, len(variants))
    return variants
