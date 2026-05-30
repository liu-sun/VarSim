# Copilot instructions for VarSim

## Build and test
- This project has no test suite. If asked to add tests, use `pytest` and place them in a `tests/` directory at the repo root.
- Build: `python -m build` (standard setuptools via pyproject.toml). No CI/CD is configured.
- Install in development mode: `pip install -e .` from the repo root.
- No linting, formatting, or type-checking configuration exists. If adding tooling, prefer `ruff` for lint/format and `mypy` for type checking.

## High-level architecture
- **Single-file package**: The entire library lives in `src/varsim/__init__.py` (~610 lines). There are no submodules.
- **Public API** (all take a `gene: str` argument, typically a gene symbol like `"INS"` or `"G6PD"`):
  - `cds(gene)` — list of (c.HGVS, p.one_letter, p.three_letter) tuples for all CDS SNVs
  - `missense(gene)` — list of (c.HGVS, p.one_letter, p.three_letter) tuples for codon-level variants with protein effect
  - `utr5(gene)` / `utr3(gene)` — list of c.HGVS strings for UTR SNVs
  - `splice_site(gene)` — list of c.HGVS strings for canonical splice site SNVs (GT-AG dinucleotides)
  - `aa_sub(gene)` — list of (p.one_letter, p.three_letter) tuples for all possible amino acid substitutions
  - `codon_sub(gene)` — list of c.HGVS strings for all possible codon substitutions
- **Internal helpers**: `nm(gene)`, `np(gene)`, `nc(gene)` fetch MANE Select/Plus Clinical records live from NCBI Entrez. They are not intended as public API but are called by every public function on every invocation (no caching).
- **Dependency**: Only `biopython`. External data comes exclusively from NCBI Entrez (no local data files).

## Key conventions
- **Environment variables required at runtime**: `EMAIL` and `API_KEY` must be set before importing `varsim`. The module reads them at import time (lines 10-11 of `__init__.py`). Any test or script using varsim must set these first.
- **Start codon special case**: The first codon (ATG/Met) receives `p.(M1?)` / `p.(Met1?)` in all protein-level outputs because a start-codon mutation's effect is unpredictable.
- **Splice site naming**: Uses the format `NC_ACCESSION(NM_ACCESSION):c.COORDINATE±OFFSET` — requires both the genomic reference (`nc()`) and the transcript record (`nm()`).
- **Code duplication**: `missense()` and `codon_sub()` share nearly identical codon-substitution classification logic (7 cases based on which bases of the 3-base codon change). Changes to one likely need mirroring in the other.
- **Return format consistency**: Functions that return protein effect (`cds`, `missense`, `aa_sub`) return triples or pairs of (one-letter, three-letter) annotations. UTR and splice functions return plain strings.
- **Module-level global**: `genetic_code` (line 12) is a list of all 64 codons (61 sense + 3 stop). It is used by `codon_sub` and `missense` to iterate over alternative codons rather than individual bases.
