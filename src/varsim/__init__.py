"""VarSim — HGVS toolkit.

Provides HGVS tooling: parsing, validation, normalization,
backtranslation, conversion, extraction, liftover, transcription,
and translation.
"""

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