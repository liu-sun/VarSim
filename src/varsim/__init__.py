"""VarSim — Sequence variant simulator for MANE transcripts.

Generates simulations for all possible single nucleotide variants (SNVs)
and frameshift indels for MANE Select/Plus Clinical transcripts,
expressed in HGVS nomenclature.

Also provides HGVS tooling: parsing, validation, normalization,
backtranslation, conversion, extraction, liftover, transcription,
and translation.
"""

from .snv import cds, utr5, utr3
from .splicing import splice_site
from .protein import aa_sub
from .codon import codon_sub
from .missense import missense
from .frameshift import frameshift
from .parser import parse
from .validator import validate, is_valid, is_valid_syntax
from .normalizer import normalize
from .backtranslate import backtranslate, backtranslate_protein
from .converter import hgvs_to_vcf, vcf_to_hgvs, hgvs_to_spdi, spdi_to_hgvs, c_to_p
from .extractor import extract
from .liftover import liftover_g_to_assembly, liftover_transcript
from .transcription import c_to_g, g_to_c
from .translation import translate_variant, get_protein_effect

__all__ = [
    # Simulation
    "cds",
    "utr5",
    "utr3",
    "splice_site",
    "aa_sub",
    "codon_sub",
    "missense",
    "frameshift",
    # HGVS tooling
    "parse",
    "validate",
    "is_valid",
    "is_valid_syntax",
    "normalize",
    "backtranslate",
    "backtranslate_protein",
    "hgvs_to_vcf",
    "vcf_to_hgvs",
    "hgvs_to_spdi",
    "spdi_to_hgvs",
    "c_to_p",
    "extract",
    "liftover_g_to_assembly",
    "liftover_transcript",
    "c_to_g",
    "g_to_c",
    "translate_variant",
    "get_protein_effect",
]