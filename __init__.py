import os

from Bio import Entrez, SeqIO
from Bio.Data.IUPACData import unambiguous_dna_letters, protein_letters, protein_letters_1to3
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio.SeqFeature import SimpleLocation

Entrez.email = os.environ["EMAIL"]
Entrez.api_key = os.environ["API_KEY"]
codons = (
    'AAA', 'AAT', 'AAG', 'AAC', 'ATA', 'ATT', 'ATG', 'ATC', 'AGA', 'AGT', 'AGG', 'AGC', 'ACA', 'ACT', 'ACG', 'ACC',
    'TAA', 'TAT', 'TAG', 'TAC', 'TTA', 'TTT', 'TTG', 'TTC', 'TGA', 'TGT', 'TGG', 'TGC', 'TCA', 'TCT', 'TCG', 'TCC',
    'GAA', 'GAT', 'GAG', 'GAC', 'GTA', 'GTT', 'GTG', 'GTC', 'GGA', 'GGT', 'GGG', 'GGC', 'GCA', 'GCT', 'GCG', 'GCC',
    'CAA', 'CAT', 'CAG', 'CAC', 'CTA', 'CTT', 'CTG', 'CTC', 'CGA', 'CGT', 'CGG', 'CGC', 'CCA', 'CCT', 'CCG', 'CCC')


def snv(gene: str) -> list:
    variants = []
    exon = []
    handle = Entrez.esearch(db="nucleotide", term=f'{
                            gene}[gene] "mane select"[keyword]')
    record = Entrez.read(handle)
    handle = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="gb", retmode="text")
    seqrecord = SeqIO.read(handle, "genbank")
    for feature in seqrecord.features:
        if feature.type == "CDS":
            protein = "".join(feature.qualifiers.get("translation"))
            protein_id = "".join(feature.qualifiers.get("protein_id"))
            global slice
            slice = feature.location.start
            cds = feature.location.extract(seqrecord).seq
            utr5 = SimpleLocation(
                0, feature.location.start).extract(seqrecord).seq
            utr3 = SimpleLocation(feature.location.end, len(
                seqrecord)).extract(seqrecord).seq
# coding region
    for index, codon in enumerate(range(0, len(cds) - 3, 3)):
        for base in unambiguous_dna_letters:
            if base != cds[codon]:
                seq = Seq(base) + cds[codon + 1: codon + 3]
                if protein[index] != seq.translate():
                    variants.append((f"{seqrecord.id}:c.{codon + 1}{cds[codon]}>{base}",
                                     f"{protein_id}:p.{protein[index]}{
                                         index + 1}{seq.translate()}",
                                     f"{protein_id}:p.{seq3(protein[index])}{index + 1}{seq3(seq.translate())}",))
                else:
                    variants.append((f"{seqrecord.id}:c.{codon + 1}{cds[codon]}>{base}",
                                     f"{protein_id}:p.{
                                         protein[index]}{index + 1}=",
                                     f"{protein_id}:p.{seq3(protein[index])}{index + 1}=",))
            if base != cds[codon + 1]:
                seq = cds[codon] + Seq(base) + cds[codon + 2]
                if protein[index] != seq.translate():
                    variants.append((f"{seqrecord.id}:c.{codon + 2}{cds[codon + 1]}>{base}",
                                     f"{protein_id}:p.{protein[index]}{
                                         index + 1}{seq.translate()}",
                                     f"{protein_id}:p.{seq3(protein[index])}{index + 1}{seq3(seq.translate())}",))
                else:
                    variants.append((f"{seqrecord.id}:c.{codon + 2}{cds[codon + 1]}>{base}",
                                     f"{protein_id}:p.{
                                         protein[index]}{index + 1}=",
                                     f"{protein_id}:p.{seq3(protein[index])}{index + 1}=",))
            if base != cds[codon + 2]:
                seq = cds[codon: codon + 2] + Seq(base)
                if protein[index] != seq.translate():
                    variants.append((f"{seqrecord.id}:c.{codon + 3}{cds[codon + 2]}>{base}",
                                     f"{protein_id}:p.{protein[index]}{
                                         index + 1}{seq.translate()}",
                                     f"{protein_id}:p.{seq3(protein[index])}{index + 1}{seq3(seq.translate())}",))
                else:
                    variants.append((f"{seqrecord.id}:c.{codon + 3}{cds[codon + 2]}>{base}",
                                     f"{protein_id}:p.{
                                         protein[index]}{index + 1}=",
                                     f"{protein_id}:p.{seq3(protein[index])}{index + 1}=",))
    # five prime untranslated region
    for index in range(len(utr5)):
        for base in unambiguous_dna_letters:
            if base != utr5[index]:
                variants.append((f"{seqrecord.id}:c.{index-len(utr5)}{utr5[index]}>{base}",
                                 "",
                                 "",))
    # three prime untranslated region
    for index in range(len(utr3)):
        for base in unambiguous_dna_letters:
            if base != utr3[index]:
                variants.append((f"{seqrecord.id}:c.*{index+1}{utr3[index]}>{base}",
                                 "",
                                 "",))
    # GT-AG rule
    for feature in seqrecord.features:
        if feature.type == "exon":
            exon.extend((feature.location.start, feature.location.end))
        for cord in range(1, len(exon)-1, 2):
            junction = (exon[cord], exon[cord+1])
            for base in unambiguous_dna_letters:
                if base != "G":
                    variants.append((f"{seqrecord.id}:c.{junction[0]-slice}+1{"G"}>{base}",
                                     "",
                                     "",))
                if base != "T":
                    variants.append((f"{seqrecord.id}:c.{junction[0]-slice}+2{"T"}>{base}",
                                     "",
                                     "",))
                if base != "A":
                    variants.append((f"{seqrecord.id}:c.{junction[1]+1-slice}-2{"A"}>{base}",
                                     "",
                                     "",))
                if base != "G":
                    variants.append((f"{seqrecord.id}:c.{junction[1]+1-slice}-1{"G"}>{base}",
                                     "",
                                     "",))

    return variants


def missense(gene: str) -> list:
    variants = []
    handle = Entrez.esearch(db="nucleotide", term=f'{
                            gene}[gene] "mane select"[keyword]')
    record = Entrez.read(handle)
    handle = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="gb", retmode="text")
    seqrecord = SeqIO.read(handle, "genbank")
    for feature in seqrecord.features:
        if feature.type == "CDS":
            protein = "".join(feature.qualifiers.get("translation"))
            protein_id = "".join(feature.qualifiers.get("protein_id"))
    for index, residue in enumerate(protein, 1):
        for aa in protein_letters:
            if aa != residue:
                variants.append((f"{protein_id}:p.{residue}{index}{aa}",
                                 f"{protein_id}:p.{protein_letters_1to3[residue]}{index}{protein_letters_1to3[aa]}"))
    return variants


def mnv(gene: str) -> list:
    variants = []
    handle = Entrez.esearch(db="nucleotide", term=f'{
                            gene}[gene] "mane select"[keyword]')
    record = Entrez.read(handle)
    handle = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="gb", retmode="text")
    seqrecord = SeqIO.read(handle, "genbank")
    for feature in seqrecord.features:
        if feature.type == "CDS":
            protein = "".join(feature.qualifiers.get("translation"))
            protein_id = "".join(feature.qualifiers.get("protein_id"))
            cds = feature.location.extract(seqrecord).seq
    for index, codon in enumerate(range(0, len(cds) - 3, 3)):
        for base in codons:
            if base != cds[codon:codon + 3]:
                seq = Seq(base)
                if protein[index] != seq.translate():
                    if base[0] == cds[codon] and base[1] == cds[codon + 1] and base[2] != cds[codon + 2]:
                        variants.append((f"{seqrecord.id}:c.{codon + 3}{cds[codon + 2]}>{base[2]}",
                                         f"{protein_id}:p.{protein[index]}{
                                             index + 1}{seq.translate()}",
                                         f"{protein_id}:p.{seq3(protein[index])}{index + 1}{seq3(seq.translate())}",))
                    elif base[0] == cds[codon] and base[1] != cds[codon + 1] and base[2] == cds[codon + 2]:
                        variants.append((f"{seqrecord.id}:c.{codon + 2}{cds[codon + 1]}>{base[1]}",
                                         f"{protein_id}:p.{protein[index]}{
                                             index + 1}{seq.translate()}",
                                         f"{protein_id}:p.{seq3(protein[index])}{index + 1}{seq3(seq.translate())}",))
                    elif base[0] != cds[codon] and base[1] == cds[codon + 1] and base[2] == cds[codon + 2]:
                        variants.append((f"{seqrecord.id}:c.{codon + 1}{cds[codon]}>{base[0]}",
                                         f"{protein_id}:p.{protein[index]}{
                                             index + 1}{seq.translate()}",
                                         f"{protein_id}:p.{seq3(protein[index])}{index + 1}{seq3(seq.translate())}",))
                    else:
                        variants.append((f"{seqrecord.id}:c.{codon + 1}_{codon + 3}{cds[codon:codon + 3]}>{base}",
                                         f"{protein_id}:p.{protein[index]}{
                                             index + 1}{seq.translate()}",
                                         f"{protein_id}:p.{seq3(protein[index])}{index + 1}{seq3(seq.translate())}",))
                else:
                    if base[0] == cds[codon] and base[1] == cds[codon + 1] and base[2] != cds[codon + 2]:
                        variants.append((f"{seqrecord.id}:c.{codon + 3}{cds[codon + 2]}>{base[2]}",
                                         f"{protein_id}:p.{
                                             protein[index]}{index + 1}=",
                                         f"{protein_id}:p.{seq3(protein[index])}{index + 1}=",))
                    elif base[0] == cds[codon] and base[1] != cds[codon + 1] and base[2] == cds[codon + 2]:
                        variants.append((f"{seqrecord.id}:c.{codon + 2}{cds[codon + 1]}>{base[1]}",
                                         f"{protein_id}:p.{
                                             protein[index]}{index + 1}=",
                                         f"{protein_id}:p.{seq3(protein[index])}{index + 1}=",))
                    elif base[0] != cds[codon] and base[1] == cds[codon + 1] and base[2] == cds[codon + 2]:
                        variants.append((f"{seqrecord.id}:c.{codon + 1}{cds[codon]}>{base[0]}",
                                         f"{protein_id}:p.{
                                             protein[index]}{index + 1}=",
                                         f"{protein_id}:p.{seq3(protein[index])}{index + 1}=",))
                    else:
                        variants.append((f"{seqrecord.id}:c.{codon + 1}_{codon + 3}{cds[codon:codon + 3]}>{base}",
                                         f"{protein_id}:p.{
                                             protein[index]}{index + 1}=",
                                         f"{protein_id}:p.{seq3(protein[index])}{index + 1}=",))
    return variants


if __name__ == "__main__":
    print(snv("SLC22A5"))
