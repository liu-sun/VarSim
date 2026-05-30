# VarSim — HGVS Variant Simulator & Toolkit

> **Note on Naming:** This tool, VarSim, is a sequence variant *simulator* for generating HGVS nomenclature. It is not affiliated with, and should not be confused with, the "VarSim" read simulator (PMID: 25524895).

VarSim is a comprehensive toolkit for HGVS variant nomenclature — simulating all possible SNV and frameshift variants for MANE transcripts, plus parsing, validation, normalization, backtranslation, format conversion, extraction, liftover, transcription, and translation. All powered by NCBI Entrez.

## Installation

```powershell
pip install varsim
```
## Configuration
Set two environment variables to query NCBI Entrez:

- **`EMAIL`** — (Required) A valid email so NCBI can contact you about query issues.
- **`API_KEY`** — (Recommended) An NCBI API key for higher query rates. Obtain one from your [NCBI account settings](https://www.ncbi.nlm.nih.gov/account/).

**Linux/macOS:**
```bash
export EMAIL="your.email@example.com"
export API_KEY="your_api_key_here"
```
**Windows (PowerShell):**
```powershell
$env:EMAIL="your.email@example.com"
$env:API_KEY="your_api_key_here"
```
## Usage

```python
import varsim
```
### 1. Variant Simulation
Generate all possible single-nucleotide variants for MANE transcripts. Results include nucleotide and protein HGVS where applicable.

| Function | Description |
|---|---|
| `cds(gene)` | All CDS SNVs → (c.HGVS, p.HGVS¹, p.HGVS³) |
| `utr5(gene)` / `utr3(gene)` | All 5′UTR / 3′UTR SNVs → list of c.HGVS |
| `splice_site(gene)` | Canonical splice site SNVs (±1, ±2) |
| `aa_sub(gene)` | All amino acid substitutions → (p¹, p³) |
| `codon_sub(gene)` | All codon-level substitutions → list of c.HGVS |
| `missense(gene)` | Codon variants with protein effect (missense / silent) |
| `frameshift(gene)` | All 1-bp deletion & insertion frameshift variants |

```python
>>> varsim.cds("INS")
[('NM_000207.3:c.1A>G', 'NP_000198.1:p.(M1?)', 'NP_000198.1:p.(Met1?)'), ...]

>>> varsim.splice_site("INS")
['NC_000011.10(NM_000207.3):c.187+1G>A', 'NC_000011.10(NM_000207.3):c.187+1G>T', ...]
```
### 2. HGVS Parsing & Validation
Parse HGVS strings into structured objects or validate syntax and semantics.

| Function | Description |
|---|---|
| `parse(hgvs)` | Parse HGVS → `HGVSTag` (`.acc`, `.prefix`, `.variant_type`, `.ref`, `.alt`, `.start_pos`, …) |
| `is_valid(hgvs, ref_seq=None)` | Return `True` if syntax (and optionally semantics) passes |
| `validate(hgvs)` | Detailed validation → list of `{"severity", "message"}` dicts |

```python
>>> tag = varsim.parse("NM_000207.3:c.1A>G")
>>> tag.variant_type, tag.ref, tag.alt
('substitution', 'A', 'G')
>>> varsim.validate("NM_000207.3:c.1A>G")
[]
>>> varsim.is_valid("NM_000207.3:c.1A>G", ref_seq="ATGCGTACG...")
True
```
### 3. HGVS Normalization
Normalize variants to canonical form per HGVS recommendations.

| Function | Description |
|---|---|
| `normalize(hgvs, ref_seq=None)` | Full pipeline: 3′ shift, ins→dup, allele minimization, range normalization |
| `normalize_3prime_shift(hgvs, ref_seq)` | Shift variant as far 3′ as possible |
| `ins_to_dup(hgvs, ref_seq)` | Convert insertion to duplication when applicable |

```python
>>> varsim.normalize("c.4A>G", ref_seq="AAGC")
'c.2A>G'

>>> varsim.ins_to_dup("NM_000207.3:c.4_5insA", ref_seq="TAAA")
'NM_000207.3:c.3dup'
```
### 4. Backtranslation
Determine which nucleotide changes could produce a given protein variant.

| Function | Description |
|---|---|
| `backtranslate(gene, p_hgvs)` | Protein → nucleotide backtranslation using the real MANE CDS |
| `backtranslate_protein(p_hgvs)` | Pure codon-table backtranslation (no gene fetch) |

```python
>>> varsim.backtranslate_protein("p.(V42G)")
['c.125T>G']
>>> varsim.backtranslate("G6PD", "p.(V42G)")  # validates against real CDS
['NM_001360016.2:c.125T>G']
```
### 5. Format Conversion
Convert between HGVS, VCF, and SPDI formats.

| Function | Description |
|---|---|
| `hgvs_to_vcf(hgvs, chrom=None)` | HGVS g./c. → VCF dict `{CHROM, POS, REF, ALT}` |
| `vcf_to_hgvs(chrom, pos, ref, alt, acc=None)` | VCF record → HGVS string |
| `hgvs_to_spdi(hgvs)` | HGVS → SPDI string |
| `spdi_to_hgvs(spdi, prefix="g.")` | SPDI → HGVS string |
| `c_to_p(c_hgvs, gene)` | Coding HGVS → protein HGVS |

```python
>>> varsim.hgvs_to_vcf("NC_000023.11:g.123456A>G")
{'CHROM': 'NC_000023.11', 'POS': 123456, 'ID': '.', 'REF': 'A', 'ALT': 'G'}
>>> varsim.vcf_to_hgvs("X", 123456, "A", "G", acc="NC_000023.11")
'NC_000023.11:g.123456A>G'
```
### 6. Variant Extraction
Diff two sequences and produce the minimal HGVS description.

| Function | Description |
|---|---|
| `extract(ref_seq, obs_seq, acc="NM_000207.3", prefix="c.")` | Align & diff → HGVS string |

```python
>>> varsim.extract("ATGC", "ATTC", prefix="c.")
'NM_000207.3:c.3G>T'
>>> varsim.extract("ATGC", "ATC", prefix="c.")
'NM_000207.3:c.3del'
```
### 7. Liftover
Remap genomic variants between assemblies via the NCBI Remap API.

| Function | Description |
|---|---|
| `liftover_g_to_assembly(hgvs, target_assembly="GRCh38")` | Lift g.HGVS between assemblies |
| `liftover_transcript(gene, c_hgvs, target_assembly="GRCh38")` | Transcript → genomic → liftover pipeline |

```python
>>> varsim.liftover_g_to_assembly("NC_000001.10:g.12345A>G", "GRCh38")
'NC_000001.11:g.12345A>G'
>>> varsim.liftover_transcript("G6PD", "c.1A>G", "GRCh38")
'NC_000023.11:g.153760607T>C'
```
### 8. Transcription
Convert between coding and genomic coordinate systems using exon structure.

| Function | Description |
|---|---|
| `c_to_g(c_hgvs, gene)` | Coding (c.) → genomic (g.) coordinates |
| `g_to_c(g_hgvs, gene)` | Genomic (g.) → coding (c.) coordinates |
| `get_cds_exon_map(gene)` | Exon structure mapping (cDNA + genomic coordinates) |

```python
>>> varsim.c_to_g("NM_001360016.2:c.1A>G", "G6PD")
'NC_000023.11:g.153760607A>G'
>>> varsim.get_cds_exon_map("G6PD")
[{'exon': 1, 'cds_start': 0, 'cds_end': 138, 'genomic_start': ..., 'strand': -1}, ...]
```
### 9. Translation
Translate coding variants to their protein consequences.

| Function | Description |
|---|---|
| `translate_variant(c_hgvs, gene)` | Coding → protein HGVS string |
| `translate_variants(c_hgvs_list, gene)` | Batch translation for multiple c.HGVS strings |
| `get_protein_effect(c_hgvs, gene)` | Effect dict: `effect_type`, `position`, `ref_aa`, `alt_aa`, 1-letter + 3-letter p.HGVS |

```python
>>> varsim.translate_variant("NM_000207.3:c.1A>G", "INS")
'NP_000198.1:p.(M1?)'
>>> eff = varsim.get_protein_effect("NM_000207.3:c.4A>G", "INS")
>>> eff["effect_type"]
'missense'
```

## API Reference

| Category | Function | Brief |
|---|---|---|
| **Simulation** | `cds(gene)` / `utr5(gene)` / `utr3(gene)` | SNVs for CDS, 5′UTR, 3′UTR → c.HGVS + p.HGVS |
| | `splice_site(gene)` / `aa_sub(gene)` / `codon_sub(gene)` | Splice-site SNVs / amino acid substitutions / codon substitutions |
| | `missense(gene)` / `frameshift(gene)` | Codon variants with protein effect / frameshift indels |
| **Parsing** | `parse(hgvs)` / `validate(hgvs)` / `is_valid(hgvs, ref_seq?)` | Parse → HGVSTag / detailed issues / bool check |
| **Normalization** | `normalize(hgvs, ref_seq?)` / `normalize_3prime_shift(...)` / `ins_to_dup(...)` | Full normalization / 3′-shift / ins→dup |
| **Backtranslation** | `backtranslate(gene, p_hgvs)` / `backtranslate_protein(p_hgvs)` | Protein → nucleotide via CDS / codon table |
| **Conversion** | `hgvs_to_vcf(...)` / `vcf_to_hgvs(...)` / `hgvs_to_spdi(...)` / `spdi_to_hgvs(...)` | HGVS ↔ VCF ↔ SPDI |
| | `c_to_p(c_hgvs, gene)` | Coding HGVS → protein HGVS |
| **Extraction** | `extract(ref, obs, acc?, prefix?)` | Diff two sequences → HGVS |
| **Liftover** | `liftover_g_to_assembly(hgvs, target?)` / `liftover_transcript(gene, c_hgvs, target?)` | Assembly liftover / transcript→genomic→liftover |
| **Transcription** | `c_to_g(c_hgvs, gene)` / `g_to_c(g_hgvs, gene)` / `get_cds_exon_map(gene)` | Coding ↔ genomic / exon structure |
| **Translation** | `translate_variant(c_hgvs, gene)` / `get_protein_effect(c_hgvs, gene)` | c.HGVS → p.HGVS / detailed effect dict |

## License

MIT License

> **Note on Naming:** This package is not affiliated with the read simulator "VarSim" (PMID: 25524895).