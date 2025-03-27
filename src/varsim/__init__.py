import os

from Bio import Entrez, SeqIO
from Bio.Data.IUPACData import (
    unambiguous_dna_letters,
    protein_letters,
    protein_letters_1to3,
    protein_letters_3to1,
)
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio.SeqFeature import SimpleLocation
from Bio.Data.CodonTable import standard_dna_table

Entrez.email = os.environ["EMAIL"]
Entrez.api_key = os.environ["API_KEY"]
codons = standard_dna_table.forward_table.keys()


def cds(gene: str) -> list:
    variants = []
    stream = Entrez.esearch(
        db="nucleotide",
        term=f'{gene}[Gene Name] "mane select"[Keyword]',
    )
    record = Entrez.read(stream)
    stream = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="gb", retmode="text"
    )
    seqrecord = SeqIO.read(stream, "genbank")
    for feature in seqrecord.features:
        if feature.type == "CDS":
            protein = "".join(feature.qualifiers.get("translation"))
            protein_id = "".join(feature.qualifiers.get("protein_id"))
            cds = feature.extract(seqrecord).seq
    for index, codon in enumerate(range(0, len(cds) - 3, 3)):
        for base in unambiguous_dna_letters:
            if base != cds[codon]:
                seq = Seq(base) + cds[codon + 1 : codon + 3]
                if protein[index] != seq.translate():
                    variants.append(
                        (
                            f"{seqrecord.id}:c.{codon + 1}{cds[codon]}>{base}",
                            f"{protein_id}:p.{protein[index]}{index + 1}{seq.translate()}",
                            f"{protein_id}:p.{seq3(protein[index])}{index + 1}{seq3(seq.translate())}",
                        )
                    )
                else:
                    variants.append(
                        (
                            f"{seqrecord.id}:c.{codon + 1}{cds[codon]}>{base}",
                            f"{protein_id}:p.{protein[index]}{index + 1}=",
                            f"{protein_id}:p.{seq3(protein[index])}{index + 1}=",
                        )
                    )
            if base != cds[codon + 1]:
                seq = cds[codon] + Seq(base) + cds[codon + 2]
                if protein[index] != seq.translate():
                    variants.append(
                        (
                            f"{seqrecord.id}:c.{codon + 2}{cds[codon + 1]}>{base}",
                            f"{protein_id}:p.{protein[index]}{index + 1}{seq.translate()}",
                            f"{protein_id}:p.{seq3(protein[index])}{index + 1}{seq3(seq.translate())}",
                        )
                    )
                else:
                    variants.append(
                        (
                            f"{seqrecord.id}:c.{codon + 2}{cds[codon + 1]}>{base}",
                            f"{protein_id}:p.{protein[index]}{index + 1}=",
                            f"{protein_id}:p.{seq3(protein[index])}{index + 1}=",
                        )
                    )
            if base != cds[codon + 2]:
                seq = cds[codon : codon + 2] + Seq(base)
                if protein[index] != seq.translate():
                    variants.append(
                        (
                            f"{seqrecord.id}:c.{codon + 3}{cds[codon + 2]}>{base}",
                            f"{protein_id}:p.{protein[index]}{index + 1}{seq.translate()}",
                            f"{protein_id}:p.{seq3(protein[index])}{index + 1}{seq3(seq.translate())}",
                        )
                    )
                else:
                    variants.append(
                        (
                            f"{seqrecord.id}:c.{codon + 3}{cds[codon + 2]}>{base}",
                            f"{protein_id}:p.{protein[index]}{index + 1}=",
                            f"{protein_id}:p.{seq3(protein[index])}{index + 1}=",
                        )
                    )
    return variants


def utr5(gene: str) -> list:
    variants = []
    stream = Entrez.esearch(
        db="nucleotide",
        term=f'{gene}[Gene Name] "mane select"[Keyword]',
    )
    record = Entrez.read(stream)
    stream = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="gb", retmode="text"
    )
    seqrecord = SeqIO.read(stream, "genbank")
    for feature in seqrecord.features:
        if feature.type == "CDS":
            utr5 = SimpleLocation(0, feature.location.start).extract(seqrecord).seq
    for index in range(len(utr5)):
        for base in unambiguous_dna_letters:
            if base != utr5[index]:
                variants.append(
                    (
                        f"{seqrecord.id}:c.{index - len(utr5)}{utr5[index]}>{base}",
                        "",
                        "",
                    )
                )
    return variants


def utr3(gene: str) -> list:
    variants = []
    stream = Entrez.esearch(
        db="nucleotide",
        term=f'{gene}[Gene Name] "mane select"[Keyword]',
    )
    record = Entrez.read(stream)
    stream = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="gb", retmode="text"
    )
    seqrecord = SeqIO.read(stream, "genbank")
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
                variants.append(
                    (
                        f"{seqrecord.id}:c.*{index + 1}{utr3[index]}>{base}",
                        "",
                        "",
                    )
                )
    return variants


def splicing(gene: str) -> list:
    variants = []
    exon = []
    stream = Entrez.esearch(
        db="nucleotide", term=f'{gene}[Gene Name] "mane select"[Keyword]'
    )
    record = Entrez.read(stream)

    stream = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="gb", retmode="text"
    )
    seqrecord = SeqIO.read(stream, "genbank")
    splicing = []
    variants = []
    start = 0
    end = 0
    for feature in seqrecord.features:
        if feature.type == "CDS":
            start = feature.location.start
            end = feature.location.end
    for feature in seqrecord.features:
        if feature.type == "exon":
            if feature.location.start < start and feature.location.end < start:
                splicing.extend(
                    (
                        feature.location.start - start - 1,
                        feature.location.end - start - 1,
                    )
                )
            elif feature.location.start < start and feature.location.end > start:
                splicing.extend(
                    (feature.location.start - start - 1, feature.location.end - start)
                )
            else:
                splicing.extend(
                    (feature.location.start - start, feature.location.end - start)
                )

    for coordinate in range(1, len(splicing) - 1, 2):
        site = splicing[coordinate], splicing[coordinate] + 1
        for base in unambiguous_dna_letters:
            if base != "G":
                variants.append((f"{seqrecord.id}:c.{site[0]}+1G>{base}"))
            if base != "T":
                variants.append((f"{seqrecord.id}:c.{site[0]}+2T>{base}"))
            if base != "A":
                variants.append((f"{seqrecord.id}:c.{site[1]}-2A>{base}"))
            if base != "G":
                variants.append((f"{seqrecord.id}:c.{site[1]}-1G>{base}"))
    return variants


def aa_sub(gene: str) -> list:
    variants = []
    term = f'{gene}[Gene Name] AND "mane select"[keyword]'
    stream = Entrez.esearch(db="protein", term=term)
    record = Entrez.read(stream)

    stream = Entrez.efetch(
        db="protein", rettype="gp", retmode="text", id=record["IdList"]
    )
    seqrecord = SeqIO.read(stream, "genbank")
    for index, residue in enumerate(seqrecord.seq, 1):
        for aa in protein_letters:
            if aa != residue:
                variants.append(
                    (
                        f"{seqrecord.id}:p.{residue}{index}{aa}",
                        f"{seqrecord.id}:p.{protein_letters_1to3[residue]}{
                            index}{protein_letters_1to3[aa]}",
                    )
                )
    return variants


def missense(gene: str) -> list:
    variants = []
    term = f'{gene}[Gene Name] "mane select"[keyword]'
    stream = Entrez.esearch(db="nucleotide", term=term)
    record = Entrez.read(stream)
    stream = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="gb", retmode="text"
    )
    seqrecord = SeqIO.read(stream, "genbank")
    for feature in seqrecord.features:
        if feature.type == "CDS":
            protein = "".join(feature.qualifiers.get("translation"))
            protein_id = "".join(feature.qualifiers.get("protein_id"))
            cds = feature.location.extract(seqrecord).seq
    for index, codon in enumerate(range(0, len(cds) - 3, 3)):
        for base in codons:
            if base != cds[codon : codon + 3]:
                seq = Seq(base)
                if protein[index] != seq.translate():
                    if (
                        base[0] == cds[codon]
                        and base[1] == cds[codon + 1]
                        and base[2] != cds[codon + 2]
                    ):
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{
                                    codon + 3}{cds[codon + 2]}>{base[2]}",
                                f"{protein_id}:p.{protein[index]}{
                                    index + 1}{seq.translate()}",
                                f"{protein_id}:p.{seq3(protein[index])}{
                                    index + 1}{seq3(seq.translate())}",
                            )
                        )
                    elif (
                        base[0] == cds[codon]
                        and base[1] != cds[codon + 1]
                        and base[2] == cds[codon + 2]
                    ):
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{
                                    codon + 2}{cds[codon + 1]}>{base[1]}",
                                f"{protein_id}:p.{protein[index]}{
                                    index + 1}{seq.translate()}",
                                f"{protein_id}:p.{seq3(protein[index])}{
                                    index + 1}{seq3(seq.translate())}",
                            )
                        )
                    elif (
                        base[0] != cds[codon]
                        and base[1] == cds[codon + 1]
                        and base[2] == cds[codon + 2]
                    ):
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{
                                    codon + 1}{cds[codon]}>{base[0]}",
                                f"{protein_id}:p.{protein[index]}{
                                    index + 1}{seq.translate()}",
                                f"{protein_id}:p.{seq3(protein[index])}{
                                    index + 1}{seq3(seq.translate())}",
                            )
                        )
                    else:
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 1}_{codon +
                                                                3}{cds[codon:codon + 3]}>{base}",
                                f"{protein_id}:p.{protein[index]}{
                                    index + 1}{seq.translate()}",
                                f"{protein_id}:p.{seq3(protein[index])}{
                                    index + 1}{seq3(seq.translate())}",
                            )
                        )
                else:
                    if (
                        base[0] == cds[codon]
                        and base[1] == cds[codon + 1]
                        and base[2] != cds[codon + 2]
                    ):
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{
                                    codon + 3}{cds[codon + 2]}>{base[2]}",
                                f"{protein_id}:p.{
                                    protein[index]}{index + 1}=",
                                f"{protein_id}:p.{seq3(protein[index])}{
                                    index + 1}=",
                            )
                        )
                    elif (
                        base[0] == cds[codon]
                        and base[1] != cds[codon + 1]
                        and base[2] == cds[codon + 2]
                    ):
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{
                                    codon + 2}{cds[codon + 1]}>{base[1]}",
                                f"{protein_id}:p.{
                                    protein[index]}{index + 1}=",
                                f"{protein_id}:p.{seq3(protein[index])}{
                                    index + 1}=",
                            )
                        )
                    elif (
                        base[0] != cds[codon]
                        and base[1] == cds[codon + 1]
                        and base[2] == cds[codon + 2]
                    ):
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{
                                    codon + 1}{cds[codon]}>{base[0]}",
                                f"{protein_id}:p.{
                                    protein[index]}{index + 1}=",
                                f"{protein_id}:p.{seq3(protein[index])}{
                                    index + 1}=",
                            )
                        )
                    else:
                        variants.append(
                            (
                                f"{seqrecord.id}:c.{codon + 1}_{codon +
                                                                3}{cds[codon:codon + 3]}>{base}",
                                f"{protein_id}:p.{
                                    protein[index]}{index + 1}=",
                                f"{protein_id}:p.{seq3(protein[index])}{
                                    index + 1}=",
                            )
                        )
    return variants


def inframe_del(gene: str) -> list:
    variants = []
    term = f'{gene}[Gene Name] "mane select"[keyword]'
    stream = Entrez.esearch(db="nucleotide", term=term)
    record = Entrez.read(stream)
    stream = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="gb", retmode="text"
    )
    seqrecord = SeqIO.read(stream, "genbank")
    for feature in seqrecord.features:
        if feature.type == "CDS":
            protein = "".join(feature.qualifiers.get("translation"))
            protein_id = "".join(feature.qualifiers.get("protein_id"))
            cds = feature.location.extract(seqrecord).seq
    for index, codon in enumerate(range(0, len(cds) - 3, 3)):
        variants.append(
            (
                f"{seqrecord.id}:c.{codon + 1}_{codon + 3}{cds[codon:codon + 3]}del",
                f"{protein_id}:p.{protein[index]}{index + 1}del",
                f"{protein_id}:p.{seq3(protein[index])}{index + 1}del",
            )
        )
    return variants


def inframe_dup(gene: str) -> list:
    variants = []
    term = f'{gene}[Gene Name] "mane select"[keyword]'
    stream = Entrez.esearch(db="nucleotide", term=term)
    record = Entrez.read(stream)
    stream = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="gb", retmode="text"
    )
    seqrecord = SeqIO.read(stream, "genbank")
    for feature in seqrecord.features:
        if feature.type == "CDS":
            protein = "".join(feature.qualifiers.get("translation"))
            protein_id = "".join(feature.qualifiers.get("protein_id"))
            cds = feature.location.extract(seqrecord).seq
    for index, codon in enumerate(range(0, len(cds) - 3, 3)):
        variants.append(
            (
                f"{seqrecord.id}:c.{codon + 1}_{codon + 3}{cds[codon:codon + 3]}dup",
                f"{protein_id}:p.{protein[index]}{index + 1}dup",
                f"{protein_id}:p.{seq3(protein[index])}{index + 1}dup",
            )
        )
    return variants


if __name__ == "__main__":
    print(splicing("TMPRSS6"))
