# VarSim — HGVS 变异模拟器与工具集

> **命名说明：** 本工具 VarSim 是一个用于生成 HGVS 命名的序列变异*模拟器*。它与读长模拟器"VarSim"（PMID: 25524895）无关，请勿混淆。

VarSim 是一个全面的 HGVS 变异命名工具集——可为 MANE 转录本模拟所有可能的 SNV 和移码变异，同时提供解析、验证、标准化、回译、格式转换、提取、Liftover、转录和翻译功能。全部依托 NCBI Entrez 实现。

## 安装

```powershell
pip install varsim
```

## 配置

查询 NCBI Entrez 数据库需设置两个环境变量：

- **`EMAIL`** —（必需）有效的邮箱地址，以便 NCBI 在查询出现问题时与您联系。
- **`API_KEY`** —（推荐）NCBI API 密钥，可大幅提高查询速率。请从 [NCBI 账户设置](https://www.ncbi.nlm.nih.gov/account/) 获取。

**Linux/macOS：**
```bash
export EMAIL="your.email@example.com"
export API_KEY="your_api_key_here"
```

**Windows (PowerShell)：**
```powershell
$env:EMAIL="your.email@example.com"
$env:API_KEY="your_api_key_here"
```

## 使用说明

```python
import varsim
```

### 1. 变异模拟

为 MANE 转录本生成所有可能的单核苷酸变异。在适用的情况下，结果包含核苷酸和蛋白质 HGVS 表示法。

| 函数 | 说明 |
|---|---|
| `cds(gene)` | 所有 CDS SNV → (c.HGVS, p.HGVS¹, p.HGVS³) |
| `utr5(gene)` / `utr3(gene)` | 所有 5′UTR / 3′UTR SNV → c.HGVS 字符串列表 |
| `splice_site(gene)` | 经典剪接位点 SNV（±1, ±2） |
| `aa_sub(gene)` | 所有氨基酸替换 → (p¹, p³) |
| `codon_sub(gene)` | 所有密码子级别替换 → c.HGVS 字符串列表 |
| `missense(gene)` | 含蛋白质效应的密码子变异（错义 / 同义） |
| `frameshift(gene)` | 所有 1bp 缺失和插入的移码变异 |

```python
>>> varsim.cds("INS")
[('NM_000207.3:c.1A>G', 'NP_000198.1:p.(M1?)', 'NP_000198.1:p.(Met1?)'), ...]

>>> varsim.splice_site("INS")
['NC_000011.10(NM_000207.3):c.187+1G>A', 'NC_000011.10(NM_000207.3):c.187+1G>T', ...]
```

### 2. HGVS 解析与验证

将 HGVS 字符串解析为结构化对象，或验证其语法和语义。

| 函数 | 说明 |
|---|---|
| `parse(hgvs)` | 解析 HGVS → `HGVSTag`（`.acc`、`.prefix`、`.variant_type`、`.ref`、`.alt`、`.start_pos` 等） |
| `is_valid(hgvs, ref_seq=None)` | 语法（及可选语义）通过则返回 `True` |
| `validate(hgvs)` | 详细验证 → `{"severity", "message"}` 字典列表 |

```python
>>> tag = varsim.parse("NM_000207.3:c.1A>G")
>>> tag.variant_type, tag.ref, tag.alt
('substitution', 'A', 'G')
>>> varsim.validate("NM_000207.3:c.1A>G")
[]
>>> varsim.is_valid("NM_000207.3:c.1A>G", ref_seq="ATGCGTACG...")
True
```

### 3. HGVS 标准化

按照 HGVS 规范将变异转换为标准形式。

| 函数 | 说明 |
|---|---|
| `normalize(hgvs, ref_seq=None)` | 完整流程：3′位移、ins→dup、等位基因最小化、范围标准化 |
| `normalize_3prime_shift(hgvs, ref_seq)` | 将变异尽量向 3′ 端移动 |
| `ins_to_dup(hgvs, ref_seq)` | 在适当情况下将插入转换为重复 |

```python
>>> varsim.normalize("c.4A>G", ref_seq="AAGC")
'c.2A>G'

>>> varsim.ins_to_dup("NM_000207.3:c.4_5insA", ref_seq="TAAA")
'NM_000207.3:c.3dup'
```

### 4. 回译

确定哪些核苷酸变化可能产生给定的蛋白质变异。

| 函数 | 说明 |
|---|---|
| `backtranslate(gene, p_hgvs)` | 蛋白质→核苷酸回译（使用真实的 MANE CDS） |
| `backtranslate_protein(p_hgvs)` | 纯密码子表回译（无需获取基因） |

```python
>>> varsim.backtranslate_protein("p.(V42G)")
['c.125T>G']
>>> varsim.backtranslate("G6PD", "p.(V42G)")  # 对照真实 CDS 验证
['NM_001360016.2:c.125T>G']
```

### 5. 格式转换

在 HGVS、VCF 和 SPDI 格式之间相互转换。

| 函数 | 说明 |
|---|---|
| `hgvs_to_vcf(hgvs, chrom=None)` | HGVS g./c. → VCF 字典 `{CHROM, POS, REF, ALT}` |
| `vcf_to_hgvs(chrom, pos, ref, alt, acc=None)` | VCF 记录 → HGVS 字符串 |
| `hgvs_to_spdi(hgvs)` | HGVS → SPDI 字符串 |
| `spdi_to_hgvs(spdi, prefix="g.")` | SPDI → HGVS 字符串 |
| `c_to_p(c_hgvs, gene)` | 编码 HGVS → 蛋白质 HGVS |

```python
>>> varsim.hgvs_to_vcf("NC_000023.11:g.123456A>G")
{'CHROM': 'NC_000023.11', 'POS': 123456, 'ID': '.', 'REF': 'A', 'ALT': 'G'}
>>> varsim.vcf_to_hgvs("X", 123456, "A", "G", acc="NC_000023.11")
'NC_000023.11:g.123456A>G'
```

### 6. 变异提取

对比两条序列并生成最精简的 HGVS 描述。

| 函数 | 说明 |
|---|---|
| `extract(ref_seq, obs_seq, acc="NM_000207.3", prefix="c.")` | 序列比对→HGVS 字符串 |

```python
>>> varsim.extract("ATGC", "ATTC", prefix="c.")
'NM_000207.3:c.3G>T'
>>> varsim.extract("ATGC", "ATC", prefix="c.")
'NM_000207.3:c.3del'
```

### 7. Liftover（基因组坐标转换）

通过 NCBI Remap API 在不同基因组组装版本之间映射变异。

| 函数 | 说明 |
|---|---|
| `liftover_g_to_assembly(hgvs, target_assembly="GRCh38")` | 在组装版本间转换 g.HGVS |
| `liftover_transcript(gene, c_hgvs, target_assembly="GRCh38")` | 转录本→基因组→Liftover 流水线 |

```python
>>> varsim.liftover_g_to_assembly("NC_000001.10:g.12345A>G", "GRCh38")
'NC_000001.11:g.12345A>G'
>>> varsim.liftover_transcript("G6PD", "c.1A>G", "GRCh38")
'NC_000023.11:g.153760607T>C'
```

### 8. 转录

利用外显子结构在编码坐标和基因组坐标之间转换。

| 函数 | 说明 |
|---|---|
| `c_to_g(c_hgvs, gene)` | 编码 (c.) → 基因组 (g.) 坐标 |
| `g_to_c(g_hgvs, gene)` | 基因组 (g.) → 编码 (c.) 坐标 |
| `get_cds_exon_map(gene)` | 外显子结构映射（cDNA + 基因组坐标） |

```python
>>> varsim.c_to_g("NM_001360016.2:c.1A>G", "G6PD")
'NC_000023.11:g.153760607A>G'
>>> varsim.get_cds_exon_map("G6PD")
[{'exon': 1, 'cds_start': 0, 'cds_end': 138, 'genomic_start': ..., 'strand': -1}, ...]
```

### 9. 翻译

将编码变异翻译为其蛋白质后果。

| 函数 | 说明 |
|---|---|
| `translate_variant(c_hgvs, gene)` | 编码 → 蛋白质 HGVS 字符串 |
| `translate_variants(c_hgvs_list, gene)` | 批量翻译多个 c.HGVS 字符串 |
| `get_protein_effect(c_hgvs, gene)` | 效应字典：`effect_type`、`position`、`ref_aa`、`alt_aa`、单字母和三字母 p.HGVS |

```python
>>> varsim.translate_variant("NM_000207.3:c.1A>G", "INS")
'NP_000198.1:p.(M1?)'
>>> eff = varsim.get_protein_effect("NM_000207.3:c.4A>G", "INS")
>>> eff["effect_type"]
'missense'
```

## API 参考

| 类别 | 函数 | 简介 |
|---|---|---|
| **模拟** | `cds(gene)` / `utr5(gene)` / `utr3(gene)` | CDS、5′UTR、3′UTR 的 SNV → c.HGVS + p.HGVS |
| | `splice_site(gene)` / `aa_sub(gene)` / `codon_sub(gene)` | 剪接位点 SNV / 氨基酸替换 / 密码子替换 |
| | `missense(gene)` / `frameshift(gene)` | 带蛋白质效应的密码子变异 / 移码 indel |
| **解析** | `parse(hgvs)` / `validate(hgvs)` / `is_valid(hgvs, ref_seq?)` | 解析→HGVSTag / 详细问题列表 / 布尔检查 |
| **标准化** | `normalize(hgvs, ref_seq?)` / `normalize_3prime_shift(...)` / `ins_to_dup(...)` | 完整标准化 / 3′位移 / ins→dup |
| **回译** | `backtranslate(gene, p_hgvs)` / `backtranslate_protein(p_hgvs)` | 通过 CDS / 密码子表实现蛋白质→核苷酸 |
| **转换** | `hgvs_to_vcf(...)` / `vcf_to_hgvs(...)` / `hgvs_to_spdi(...)` / `spdi_to_hgvs(...)` | HGVS ↔ VCF ↔ SPDI |
| | `c_to_p(c_hgvs, gene)` | 编码 HGVS → 蛋白质 HGVS |
| **提取** | `extract(ref, obs, acc?, prefix?)` | 双序列比对→HGVS |
| **Liftover** | `liftover_g_to_assembly(hgvs, target?)` / `liftover_transcript(gene, c_hgvs, target?)` | 组装版本转换 / 转录本→基因组→Liftover |
| **转录** | `c_to_g(c_hgvs, gene)` / `g_to_c(g_hgvs, gene)` / `get_cds_exon_map(gene)` | 编码↔基因组 / 外显子结构 |
| **翻译** | `translate_variant(c_hgvs, gene)` / `get_protein_effect(c_hgvs, gene)` | c.HGVS→p.HGVS / 详细效应字典 |

## 许可证

MIT License

> **命名说明：** 本软件包与读长模拟器"VarSim"（PMID: 25524895）无关。
