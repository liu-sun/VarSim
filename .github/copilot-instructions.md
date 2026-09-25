# Copilot instructions for VarSim

## Build and test
- Tests live in `tests/` and use `pytest`. Entrez-dependent tests are skipped unless `EMAIL` and `API_KEY` are set (see `tests/conftest.py`).
- Build: `python -m build` (standard setuptools via pyproject.toml). No CI/CD is configured.
- Install in development mode: `pip install -e .` from the repo root.
- No linting, formatting, or type-checking configuration exists. If adding tooling, prefer `ruff` for lint/format and `mypy` for type checking.

## High-level architecture
- **Multi-module package** under `src/varsim/`: `parser.py`, `validator.py`, `normalizer.py`, `backtranslate.py`, `converter.py`, `extractor.py`, `liftover.py`, `transcription.py`, `translation.py`, plus internal helpers (`_fetch.py`, `_utils.py`, `_logging.py`). `__init__.py` re-exports the public API.
- **Public API** (HGVS tooling only — no variant simulation):
  - Parsing & validation: `parse`, `validate`, `is_valid`, `is_valid_syntax`
  - Normalization: `normalize`
  - Backtranslation: `backtranslate`, `backtranslate_protein`
  - Conversion: `hgvs_to_vcf`, `vcf_to_hgvs`, `hgvs_to_spdi`, `spdi_to_hgvs`, `c_to_p`
  - Extraction: `extract`
  - Liftover: `liftover_g_to_assembly`, `liftover_transcript`
  - Transcription: `c_to_g`, `g_to_c`
  - Translation: `translate_variant`, `get_protein_effect`
- **Internal helpers**: `nm(gene)`, `np(gene)`, `nc(gene)` in `_fetch.py` fetch MANE Select/Plus Clinical records live from NCBI Entrez and are memoized with `functools.lru_cache`.
- **Dependency**: Only `biopython`. External data comes exclusively from NCBI Entrez (no local data files).

## Key conventions
- **Environment variables required at runtime**: `EMAIL` and `API_KEY` must be set before importing `varsim`. `_fetch.py` reads them at import time (lines 20-21). Any test or script using varsim must set these first.
- **Module-level global**: `genetic_code` in `_utils.py` is a list of all 64 codons (61 sense + 3 stop). It is used by `backtranslate` to enumerate codons.
