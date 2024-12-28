import os

from Bio import Entrez, SeqIO
from Bio.Data.IUPACData import unambiguous_dna_letters

Entrez.email = os.environ.get("EMAIL")
Entrez.api_key = os.environ.get("API_KEY")


def splicing(gene: str) -> list:
    term = f'{gene}[GENE] AND "MANE Select"[kEYWORD]'
    stream = Entrez.esearch(db="nucleotide", term=term)
    record = Entrez.read(stream)
    stream = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="gb", retmode="text"
    )
    seqrecord = SeqIO.read(stream, "genbank")
    splicing = []
    variants = []
    start = 0
    for feature in seqrecord.features:
        if feature.type == "CDS":
            start = feature.location.start
    for feature in seqrecord.features:
        if feature.type == "exon":
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


if __name__ == "__main__":
    import pprint

    pprint.pprint(splicing("SLC22A5"))
