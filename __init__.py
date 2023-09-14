import configparser

from Bio import Entrez, SeqIO
from Bio.Seq import Seq

config = configparser.ConfigParser()
config.read("config.ini")
Entrez.email = config["DEFAULT"]["email"]
Entrez.api_key = config["DEFAULT"]["api_key"]


def simulate(gene):
    handle = Entrez.esearch(
        db='nucleotide', term=f'{gene}[gene] "mane select"[keyword]')
    record = Entrez.read(handle)
    idlist = record["IdList"]
    handle = Entrez.efetch(db='nucleotide', id=idlist,
                           rettype='gb', retmode='text')
    seqrecord = SeqIO.read(handle, 'genbank')
    for feature in seqrecord.features:
        if feature.type == "CDS":
            protein = "".join(feature.qualifiers.get("translation"))
            protein_id = "".join(feature.qualifiers.get('protein_id'))
            cds = feature.location.extract(seqrecord).seq
    for index, codon in enumerate(range(0, len(cds)-3, 3)):
        for base in "ATCG":
            if base != cds[codon]:
                seq = Seq(base)+cds[codon+1:codon+3]
                print(f"{seqrecord.id}:c.{codon+1}{cds[codon]}>{base}",
                      f"{protein_id}:p.{protein[index]}{index+1}{seq.translate()}")
            if base != cds[codon+1]:
                seq = cds[codon]+Seq(base)+cds[codon+2]
                print(f"{seqrecord.id}:c.{codon+2}{cds[codon+1]}>{base}",
                      f"{protein_id}:p.{protein[index]}{index+1}{seq.translate()}")
            if base != cds[codon+2]:
                seq = cds[codon:codon+2]+Seq(base)
                print(f"{seqrecord.id}:c.{codon+3}{cds[codon+2]}>{base}",
                      f"{protein_id}:p.{protein[index]}{index+1}{seq.translate()}")


if __name__ == "__main__":
    simulate("SLC22A5")
