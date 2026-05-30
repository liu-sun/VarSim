"""Internal module: shared constants and utilities."""

from Bio.Data.CodonTable import standard_dna_table

# All 64 codons: 61 sense codons + 3 stop codons
genetic_code = list(standard_dna_table.forward_table.keys()) + standard_dna_table.stop_codons
