import itertools
import os

from Bio import Entrez, SeqIO
from Bio.Data.IUPACData import unambiguous_dna_letters, protein_letters
from Bio.Seq import Seq
from Bio.SeqFeature import SimpleLocation
from Bio.SeqUtils import seq3

Entrez.email = os.environ["EMAIL"]
Entrez.api_key = os.environ["API_KEY"]
genetic_code = tuple(
    i + j + k for i, j, k in itertools.product(unambiguous_dna_letters, repeat=3)
)


def nm(gene: str) -> str:
    stream = Entrez.esearch(
        db="nucleotide",
        term=f'{gene}[Gene Name] AND ("MANE Select"[Keyword] OR "MANE Plus Clinical"[keyword])',
    )
    record = Entrez.read(stream)
    stream = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="gb", retmode="text"
    )
    seqrecord = SeqIO.read(stream, "genbank")
    return seqrecord


def np(gene: str) -> str:
    stream = Entrez.esearch(
        db="protein",
        term=f'{gene}[Gene Name] AND ("MANE Select"[Keyword] OR "MANE Plus Clinical"[keyword])',
    )
    record = Entrez.read(stream)
    stream = Entrez.efetch(
        db="protein", id=record["IdList"], rettype="fasta", retmode="text"
    )
    seqrecord = SeqIO.read(stream, "fasta")
    return seqrecord


def nc(gene: str) -> str:
    stream = Entrez.esearch(
        db="nucleotide",
        term=f'{gene}[Gene Name] AND "Primary Assembly"[Title] AND human[Organism]',
    )
    record = Entrez.read(stream)
    stream = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="acc", retmode="text"
    )
    acc = stream.read().strip()
    return acc


def cds(gene: str) -> list:
    variants = []
    seqrecord = nm(gene)
    protein_seqrecord = np(gene)
    protein = str(protein_seqrecord.seq)
    protein_id = protein_seqrecord.id
    for feature in seqrecord.features:
        if feature.type == "CDS":
            cds = feature.extract(seqrecord).seq
    for index, codon in enumerate(range(0, len(cds) - 3, 3)):
        for base in unambiguous_dna_letters:
            if index == 0:
                if base != "A":
                    variants.append(
                        (
                            f"{seqrecord.id}:c.1A>{base}",
                            f"{protein_id}:p.(M1?)",
                            f"{protein_id}:p.(Met1?)",
                        )
                    )
            else:
                if base != cds[codon]:
                    seq = Seq(base) + cds[codon + 1 : codon + 3]
                    if protein[index] != seq.translate():
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 1}{cds[codon]}>{base}",
                                f"{protein_id}:p.({protein[index]}{index + 1}{seq.translate()})",
                                f"{protein_id}:p.({seq3(protein[index])}{index + 1}{seq3(seq.translate())})",
                            )
                        )
                    else:
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 1}{cds[codon]}>{base}",
                                f"{protein_id}:p.({protein[index]}{index + 1}=)",
                                f"{protein_id}:p.({seq3(protein[index])}{index + 1}=)",
                            )
                        )
            if index == 0:
                if base != "T":
                    variants.append(
                        (
                            f"{seqrecord.id}:c.2T>{base}",
                            f"{protein_id}:p.(M1?)",
                            f"{protein_id}:p.(Met1?)",
                        )
                    )
            else:
                if base != cds[codon + 1]:
                    seq = cds[codon] + Seq(base) + cds[codon + 2]
                    if protein[index] != seq.translate():
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 2}{cds[codon + 1]}>{base}",
                                f"{protein_id}:p.({protein[index]}{index + 1}{seq.translate()})",
                                f"{protein_id}:p.({seq3(protein[index])}{index + 1}{seq3(seq.translate())})",
                            )
                        )
                    else:
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 2}{cds[codon + 1]}>{base}",
                                f"{protein_id}:p.({protein[index]}{index + 1}=)",
                                f"{protein_id}:p.({seq3(protein[index])}{index + 1}=)",
                            )
                        )
            if index == 0:
                if base != "G":
                    variants.append(
                        (
                            f"{seqrecord.id}:c.3G>{base}",
                            f"{protein_id}:p.(M1?)",
                            f"{protein_id}:p.(Met1?)",
                        )
                    )
            else:
                if base != cds[codon + 2]:
                    seq = cds[codon : codon + 2] + Seq(base)
                    if protein[index] != seq.translate():
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 3}{cds[codon + 2]}>{base}",
                                f"{protein_id}:p.({protein[index]}{index + 1}{seq.translate()})",
                                f"{protein_id}:p.({seq3(protein[index])}{index + 1}{seq3(seq.translate())})",
                            )
                        )
                    else:
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 3}{cds[codon + 2]}>{base}",
                                f"{protein_id}:p.({protein[index]}{index + 1}=)",
                                f"{protein_id}:p.({seq3(protein[index])}{index + 1}=)",
                            )
                        )
    return variants


def utr5(gene: str) -> list:
    variants = []
    seqrecord = nm(gene)
    for feature in seqrecord.features:
        if feature.type == "CDS":
            utr5 = SimpleLocation(0, feature.location.start).extract(seqrecord).seq
    for index in range(len(utr5)):
        for base in unambiguous_dna_letters:
            if base != utr5[index]:
                variants.append(
                    f"{seqrecord.id}:c.{index - len(utr5)}{utr5[index]}>{base}"
                )
    return variants


def utr3(gene: str) -> list:
    variants = []
    seqrecord = nm(gene)
    for feature in seqrecord.features:
        if feature.type == "CDS":
            utr3 = (
                SimpleLocation(feature.location.end, len(seqrecord))
                .extract(seqrecord)
                .seq
            )
    for index in range(len(utr3)):
        for base in unambiguous_dna_letters:
            if base != utr3[index]:
                variants.append(f"{seqrecord.id}:c.*{index + 1}{utr3[index]}>{base}")
    return variants


def splice_site(gene: str) -> list:
    variants = []
    seqrecord = nm(gene)
    acc = nc(gene)
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
    return variants


def aa_sub(gene: str) -> list:
    variants = []
    seqrecord = np(gene)
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
    return variants


def codon_sub(gene: str) -> list:
    variants = []
    seqrecord = nm(gene)
    for feature in seqrecord.features:
        if feature.type == "CDS":
            cds = feature.location.extract(seqrecord).seq
    for index, codon in enumerate(range(0, len(cds) - 3, 3)):
        for base in genetic_code:
            if base != cds[codon : codon + 3]:
                seq = Seq(base)
                if all(
                    (
                        base[0] != cds[codon],
                        base[1] != cds[codon + 1],
                        base[2] != cds[codon + 2],
                    )
                ):
                    variants.append(
                        f"{seqrecord.id}:c.{codon + 1}_{codon + 3}delins{base}"
                    )
                elif all(
                    (
                        base[0] == cds[codon],
                        base[1] != cds[codon + 1],
                        base[2] != cds[codon + 2],
                    )
                ):
                    variants.append(
                        f"{seqrecord.id}:c.{codon + 2}_{codon + 3}delins{base[1:3]}"
                    )
                elif all(
                    (
                        base[0] != cds[codon],
                        base[1] == cds[codon + 1],
                        base[2] != cds[codon + 2],
                    )
                ):
                    variants.append(
                        f"{seqrecord.id}:c.{codon + 1}_{codon + 3}delins{base}"
                    )
                elif all(
                    (
                        base[0] != cds[codon],
                        base[1] != cds[codon + 1],
                        base[2] == cds[codon + 2],
                    )
                ):
                    variants.append(
                        f"{seqrecord.id}:c.{codon + 1}_{codon + 2}delins{base[0:2]}"
                    )
                elif all(
                    (
                        base[0] == cds[codon],
                        base[1] == cds[codon + 1],
                        base[2] != cds[codon + 2],
                    )
                ):
                    variants.append(
                        f"{seqrecord.id}:c.{codon + 3}{cds[codon + 2]}>{base[2]}"
                    )
                elif all(
                    (
                        base[0] == cds[codon],
                        base[1] != cds[codon + 1],
                        base[2] == cds[codon + 2],
                    )
                ):
                    variants.append(
                        f"{seqrecord.id}:c.{codon + 2}{cds[codon + 1]}>{base[1]}"
                    )
                elif all(
                    (
                        base[0] != cds[codon],
                        base[1] == cds[codon + 1],
                        base[2] == cds[codon + 2],
                    )
                ):
                    variants.append(
                        f"{seqrecord.id}:c.{codon + 1}{cds[codon]}>{base[0]}"
                    )
    return variants


def missense(gene: str) -> list:
    variants = []
    seqrecord = nm(gene)
    protein_seqrecord = np(gene)
    protein = str(protein_seqrecord.seq)
    protein_id = protein_seqrecord.id
    for feature in seqrecord.features:
        if feature.type == "CDS":
            cds = feature.location.extract(seqrecord).seq
    for index, codon in enumerate(range(0, len(cds) - 3, 3)):
        for base in genetic_code:
            if base != cds[codon : codon + 3]:
                seq = Seq(base)
                if index == 0:
                    if all(
                        (
                            base[0] != cds[codon],
                            base[1] != cds[codon + 1],
                            base[2] != cds[codon + 2],
                        )
                    ):
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 1}_{codon + 3}delins{base}",
                                f"{protein_id}:p.(M1?)",
                                f"{protein_id}:p.(Met1?)",
                            )
                        )
                    elif all(
                        (
                            base[0] == cds[codon],
                            base[1] != cds[codon + 1],
                            base[2] != cds[codon + 2],
                        )
                    ):
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 2}_{codon + 3}delins{base[1:3]}",
                                f"{protein_id}:p.(M1?)",
                                f"{protein_id}:p.(Met1?)",
                            )
                        )
                    elif all(
                        (
                            base[0] != cds[codon],
                            base[1] == cds[codon + 1],
                            base[2] != cds[codon + 2],
                        )
                    ):
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 1}_{codon + 3}delins{base}",
                                f"{protein_id}:p.(M1?)",
                                f"{protein_id}:p.(Met1?)",
                            )
                        )
                    elif all(
                        (
                            base[0] != cds[codon],
                            base[1] != cds[codon + 1],
                            base[2] == cds[codon + 2],
                        )
                    ):
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 1}_{codon + 2}delins{base[0:2]}",
                                f"{protein_id}:p.(M1?)",
                                f"{protein_id}:p.(Met1?)",
                            )
                        )
                    elif all(
                        (
                            base[0] == cds[codon],
                            base[1] == cds[codon + 1],
                            base[2] != cds[codon + 2],
                        )
                    ):
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 3}{cds[codon + 2]}>{base[2]}",
                                f"{protein_id}:p.(M1?)",
                                f"{protein_id}:p.(Met1?)",
                            )
                        )
                    elif all(
                        (
                            base[0] == cds[codon],
                            base[1] != cds[codon + 1],
                            base[2] == cds[codon + 2],
                        )
                    ):
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 2}{cds[codon + 1]}>{base[1]}",
                                f"{protein_id}:p.(M1?)",
                                f"{protein_id}:p.(Met1?)",
                            )
                        )
                    elif all(
                        (
                            base[0] != cds[codon],
                            base[1] == cds[codon + 1],
                            base[2] == cds[codon + 2],
                        )
                    ):
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 1}{cds[codon]}>{base[0]}",
                                f"{protein_id}:p.(M1?)",
                                f"{protein_id}:p.(Met1?)",
                            )
                        )
                else:
                    if all(
                        (
                            base[0] != cds[codon],
                            base[1] != cds[codon + 1],
                            base[2] != cds[codon + 2],
                        )
                    ):
                        if protein[index] != seq.translate():
                            variants.append(
                                (
                                    f"{seqrecord.id}:c.{codon + 1}_{codon + 3}delins{base}",
                                    f"{protein_id}:p.({protein[index]}{index + 1}{seq.translate()})",
                                    f"{protein_id}:p.({seq3(protein[index])}{index + 1}{seq3(seq.translate())})",
                                )
                            )
                        else:
                            variants.append(
                                (
                                    f"{seqrecord.id}:c.{codon + 1}_{codon + 3}delins{base}",
                                    f"{protein_id}:p.({protein[index]}{index + 1}=)",
                                    f"{protein_id}:p.({seq3(protein[index])}{index + 1}=)",
                                )
                            )
                    elif all(
                        (
                            base[0] == cds[codon],
                            base[1] != cds[codon + 1],
                            base[2] != cds[codon + 2],
                        )
                    ):
                        if protein[index] != seq.translate():
                            variants.append(
                                (
                                    f"{seqrecord.id}:c.{codon + 2}_{codon + 3}delins{base[1:3]}",
                                    f"{protein_id}:p.({protein[index]}{index + 1}{seq.translate()})",
                                    f"{protein_id}:p.({seq3(protein[index])}{index + 1}{seq3(seq.translate())})",
                                )
                            )
                        else:
                            variants.append(
                                (
                                    f"{seqrecord.id}:c.{codon + 2}_{codon + 3}delins{base[1:3]}",
                                    f"{protein_id}:p.({protein[index]}{index + 1}=)",
                                    f"{protein_id}:p.({seq3(protein[index])}{index + 1}=)",
                                )
                            )
                    elif all(
                        (
                            base[0] != cds[codon],
                            base[1] == cds[codon + 1],
                            base[2] != cds[codon + 2],
                        )
                    ):
                        if protein[index] != seq.translate():
                            variants.append(
                                (
                                    f"{seqrecord.id}:c.{codon + 1}_{codon + 3}delins{base}",
                                    f"{protein_id}:p.({protein[index]}{index + 1}{seq.translate()})",
                                    f"{protein_id}:p.({seq3(protein[index])}{index + 1}{seq3(seq.translate())})",
                                )
                            )
                        else:
                            variants.append(
                                (
                                    f"{seqrecord.id}:c.{codon + 1}_{codon + 3}delins{base}",
                                    f"{protein_id}:p.({protein[index]}{index + 1}=)",
                                    f"{protein_id}:p.({seq3(protein[index])}{index + 1}=)",
                                )
                            )
                    elif all(
                        (
                            base[0] != cds[codon],
                            base[1] != cds[codon + 1],
                            base[2] == cds[codon + 2],
                        )
                    ):
                        if protein[index] != seq.translate():
                            variants.append(
                                (
                                    f"{seqrecord.id}:c.{codon + 1}_{codon + 2}delins{base[0:2]}",
                                    f"{protein_id}:p.({protein[index]}{index + 1}{seq.translate()})",
                                    f"{protein_id}:p.({seq3(protein[index])}{index + 1}{seq3(seq.translate())})",
                                )
                            )
                        else:
                            variants.append(
                                (
                                    f"{seqrecord.id}:c.{codon + 1}_{codon + 2}delins{base[0:2]}",
                                    f"{protein_id}:p.({protein[index]}{index + 1}=)",
                                    f"{protein_id}:p.({seq3(protein[index])}{index + 1}=)",
                                )
                            )
                    elif all(
                        (
                            base[0] == cds[codon],
                            base[1] == cds[codon + 1],
                            base[2] != cds[codon + 2],
                        )
                    ):
                        if protein[index] != seq.translate():
                            variants.append(
                                (
                                    f"{seqrecord.id}:c.{codon + 3}{cds[codon + 2]}>{base[2]}",
                                    f"{protein_id}:p.({protein[index]}{index + 1}{seq.translate()})",
                                    f"{protein_id}:p.({seq3(protein[index])}{index + 1}{seq3(seq.translate())})",
                                )
                            )
                        else:
                            variants.append(
                                (
                                    f"{seqrecord.id}:c.{codon + 3}{cds[codon + 2]}>{base[2]}",
                                    f"{protein_id}:p.({protein[index]}{index + 1}=)",
                                    f"{protein_id}:p.({seq3(protein[index])}{index + 1}=)",
                                )
                            )
                    elif all(
                        (
                            base[0] == cds[codon],
                            base[1] != cds[codon + 1],
                            base[2] == cds[codon + 2],
                        )
                    ):
                        if protein[index] != seq.translate():
                            variants.append(
                                (
                                    f"{seqrecord.id}:c.{codon + 2}{cds[codon + 1]}>{base[1]}",
                                    f"{protein_id}:p.({protein[index]}{index + 1}{seq.translate()})",
                                    f"{protein_id}:p.({seq3(protein[index])}{index + 1}{seq3(seq.translate())})",
                                )
                            )
                        else:
                            variants.append(
                                (
                                    f"{seqrecord.id}:c.{codon + 2}{cds[codon + 1]}>{base[1]}",
                                    f"{protein_id}:p.({protein[index]}{index + 1}=)",
                                    f"{protein_id}:p.({seq3(protein[index])}{index + 1}=)",
                                )
                            )
                    elif all(
                        (
                            base[0] != cds[codon],
                            base[1] == cds[codon + 1],
                            base[2] == cds[codon + 2],
                        )
                    ):
                        if protein[index] != seq.translate():
                            variants.append(
                                (
                                    f"{seqrecord.id}:c.{codon + 1}{cds[codon]}>{base[0]}",
                                    f"{protein_id}:p.({protein[index]}{index + 1}{seq.translate()})",
                                    f"{protein_id}:p.({seq3(protein[index])}{index + 1}{seq3(seq.translate())})",
                                )
                            )
                        else:
                            variants.append(
                                (
                                    f"{seqrecord.id}:c.{codon + 1}{cds[codon]}>{base[0]}",
                                    f"{protein_id}:p.({protein[index]}{index + 1}=)",
                                    f"{protein_id}:p.({seq3(protein[index])}{index + 1}=)",
                                )
                            )
    return variants


if __name__ == "__main__":
    import pprint

    # pprint.pprint(cds("INS"))
    # pprint.pprint(missense("INS"))
    # pprint.pprint(aa_sub("INS"))
    # pprint.pprint(utr5("INS"))
    # pprint.pprint(utr3("INS"))
    # pprint.pprint(splice_site("INS"))
    # pprint.pprint(codon_sub("INS"))
    # pprint.pprint(nc("INS"))
    # pprint.pprint(np("INS"))
    # pprint.pprint(nm("INS"))
