# VarSim API Reference

All public functions organized by category.

---

## Simulation

Generate all possible HGVS variant descriptions for a gene's MANE transcript.

| Function | Returns | Description |
|---|---|---|
| `cds(gene)` | `list[tuple]` | All CDS single-nucleotide variants. Each: `(c.HGVS, p.HGVS_1l, p.HGVS_3l)` |
| `missense(gene)` | `list[tuple]` | All codon-level variants with protein effect (missense/silent/start-loss) |
| `utr5(gene)` | `list[str]` | All 5-prime UTR SNVs (`c.-N`) |
| `utr3(gene)` | `list[str]` | All 3-prime UTR SNVs (`c.*N`) |
| `splice_site(gene)` | `list[str]` | Canonical GT-AG splice dinucleotide SNVs (`c.N±O`) |
| `aa_sub(gene)` | `list[tuple]` | All possible amino acid substitutions. Each: `(p.HGVS_1l, p.HGVS_3l)` |
| `codon_sub(gene)` | `list[str]` | All codon-to-codon substitutions as nucleotide HGVS strings |
| `frameshift(gene)` | `list[tuple]` | All 1bp deletions + 4×1bp insertions per CDS position. Each: `(c.HGVS, p.HGVS_1l, p.HGVS_3l)` |

---

## Parsing

| Function | Returns | Description |
|---|---|---|
| `parse(hgvs_str)` | `HGVSTag` | Parse any HGVS string into a structured NamedTuple |

**HGVSTag fields:** `acc`, `genomic_acc`, `prefix`, `start_pos`, `end_pos`, `start_offset`, `end_offset`, `ref`, `alt`, `variant_type`, `is_uncertain`, `fs_length`, `original`

---

## Validation

| Function | Returns | Description |
|---|---|---|
| `is_valid_syntax(hgvs_str)` | `bool` | Check if string can be parsed (syntax only) |
| `validate(hgvs_str)` | `list[dict]` | Full validation with issue list. Issues: `{severity, message}` |
| `validate_semantic(hgvs_str, ref_seq=None)` | `list[dict]` | Semantic validation (ref allele match, coordinate bounds) |
| `is_valid(hgvs_str, ref_seq=None)` | `bool` | Combined syntax + semantic validation |

---

## Normalization

| Function | Returns | Description |
|---|---|---|
| `normalize(hgvs_str, ref_seq=None)` | `str` | Full normalization: 3′ shift → ins→dup → allele minimization |
| `normalize_3prime_shift(hgvs_str, ref_seq)` | `str` | 3′ shift only |
| `ins_to_dup(hgvs_str, ref_seq)` | `str` | Convert insertion to duplication if applicable |

---

## Backtranslation

| Function | Returns | Description |
|---|---|---|
| `backtranslate(gene, p_hgvs)` | `list[str]` | Protein HGVS → list of possible c.HGVS strings (uses NCBI) |
| `backtranslate_protein(p_hgvs)` | `list[str]` | Protein HGVS → list of possible c.HGVS strings (genetic code only) |

---

## Conversion

| Function | Returns | Description |
|---|---|---|
| `hgvs_to_vcf(hgvs_str, chrom=None)` | `dict` | HGVS → VCF `{CHROM, POS, ID, REF, ALT}` |
| `vcf_to_hgvs(chrom, pos, ref, alt, acc, prefix)` | `str` | VCF → HGVS string |
| `hgvs_to_spdi(hgvs_str)` | `str` | HGVS → SPDI `ACC:POS:REF:ALT` (0-based) |
| `spdi_to_hgvs(spdi_str, prefix)` | `str` | SPDI → HGVS string |
| `c_to_p(hgvs_str, gene)` | `str` | Coding HGVS → protein HGVS |

---

## Extraction

| Function | Returns | Description |
|---|---|---|
| `extract(ref_seq, obs_seq, acc, prefix)` | `str` | Diff two sequences → minimal HGVS description |

---

## Liftover

| Function | Returns | Description |
|---|---|---|
| `liftover_g_to_assembly(hgvs_str, target_assembly)` | `str` | Lift genomic HGVS between assemblies (NCBI Remap API) |
| `liftover_transcript(gene, c_hgvs, target_assembly)` | `str` | Gene + c.HGVS → g.HGVS → liftover |

---

## Transcription

| Function | Returns | Description |
|---|---|---|
| `c_to_g(c_hgvs, gene)` | `str` | Coding HGVS → genomic HGVS |
| `g_to_c(g_hgvs, gene)` | `str` | Genomic HGVS → coding HGVS |
| `get_cds_exon_map(gene)` | `list[dict]` | Exon structure mapping for a gene |

---

## Translation

| Function | Returns | Description |
|---|---|---|
| `translate_variant(c_hgvs, gene)` | `str` | Coding HGVS → protein HGVS |
| `translate_variants(c_hgvs_list, gene)` | `list[str]` | Batch translation (single API fetch) |
| `get_protein_effect(c_hgvs, gene)` | `dict` | Detailed effect: `{p_hgvs_1letter, p_hgvs_3letter, effect_type, position, ref_aa, alt_aa}` |

**Effect types:** `missense`, `silent`, `nonsense`, `frameshift`, `stop_loss`, `start_loss`, `extension`, `non-coding`
