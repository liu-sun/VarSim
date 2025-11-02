# VarSim

> **Note on Naming:** This tool, VarSim, is a sequence variant *simulator* for generating HGVS nomenclature. It is not affiliated with, and should not be confused with, the "VarSim" read simulator (PMID: 25524895).

VarSim generates simulations for all possible single nucleotide variants (SNVs) for the Matched Annotation from NCBI and EMBL-EBI (MANE) transcript, along with the corresponding protein using HGVS notation.

**INSTALLATION**

```powershell
pip install varsim
```
CONFIGURATION (Mandatory)

Before you can use VarSim, you must set two environment variables to query the NCBI Entrez database.

    EMAIL: (Required by NCBI) A valid email address so NCBI can contact you if there's an issue with your queries.

    API_KEY: (Recommended) An NCBI API key, which allows for a much higher query rate. You can obtain one from your NCBI account settings.

Linux/macOS:
```Bash

export EMAIL="your.email@example.com"
export API_KEY="your_api_key_here"
```
Windows (PowerShell):

```powershell
set EMAIL="your.email@example.com"
set API_KEY="your_api_key_here"
```
USAGE

Variant Simulator
```Python

import varsim
# Assumes EMAIL and API_KEY are set
varsim.cds("INS")
```

```python
[('NM_000207.3:c.1A>G', 'NP_000198.1:p.(M1?)', 'NP_000198.1:p.(Met1?)'), ('NM_000207.3:c.1A>T', 'NP_000198.1:p.(M1?)', 'NP_000198.1:p.(Met1?)'), ('NM_000207.3:c.1A>C', 'NP_000198.1:p.(M1?)', 'NP_000198.1:p.(Met1?)'), ('NM_000207.3:c.2T>G', 'NP_000198.1:p.(M1?)', 'NP_000198.1:p.(Met1?)'), ('NM_000207.3:c.2T>A', 'NP_000198.1:p.(M1?)', 'NP_000198.1:p.(Met1?)'), ('NM_000207.3:c.2T>C', 'NP_000198.1:p.(M1?)', 'NP_000198.1:p.(Met1?)'), ('NM_000207.3:c.3G>A', 'NP_000198.1:p.(M1?)', 'NP_000198.1:p.(Met1?)'), ('NM_000207.3:c.3G>T', 'NP_000198.1:p.(M1?)', 'NP_000198.1:p.(Met1?)'), ('NM_000207.3:c.3G>C', 'NP_000198.1:p.(M1?)', 'NP_000198.1:p.(Met1?)'), ... ('NM_000207.3:c.328A>T', 'NP_000198.1:p.(N110Y)', 'NP_000198.1:p.(Asn110Tyr)'), ('NM_000207.3:c.328A>C', 'NP_000198.1:p.(N110H)', 'NP_000198.1:p.(Asn110His)'), ('NM_000207.3:c.329A>T', 'NP_000198.1:p.(N110I)', 'NP_000198.1:p.(Asn110Ile)'), ('NM_000207.3:c.329A>C', 'NP_000198.1:p.(N110T)', 'NP_000198.1:p.(Asn110Thr)'), ('NM_000207.3:c.330C>T', 'NP_000198.1:p.(N110=)', 'NP_000198.1:p.(Asn110=)'), ('NM_000207.3:c.330C>A', 'NP_000198.1:p.(N110=)', 'NP_000198.1:p.(Asn110=)'), ('NM_000207.3:c.330C>G', 'NP_000198.1:p.(N110=)', 'NP_000198.1:p.(Asn110=)')]
```
```Python

varsim.utr5("INS")
```
```python
['NM_000207.3:c.-59A>G', 'NM_000207.3:c.-59A>T', 'NM_000207.3:c.-59A>C', 'NM_000207.3:c.-58G>A', 'NM_000207.3:c.-58G>T', ... 'NM_000207.3:c.-2C>A', 'NM_000207.3:c.-2C>T', 'NM_000207.3:c.-1C>G', 'NM_000207.3:c.-1C>A', 'NM_000207.3:c.-1C>T']
```
```Python

varsim.utr3("INS")
```
```python
['NM_000207.3:c.*1A>G', 'NM_000207.3:c.*1A>T', 'NM_000207.3:c.*1A>C', 'NM_000207.3:c.*2C>G', 'NM_000207.3:c.*2C>A', ... 'NM_000207.3:c.*72G>T', 'NM_000207.3:c.*72G>C', 'NM_000207.3:c.*73C>G', 'NM_000207.3:c.*73C>A', 'NM_000207.3:c.*73C>T']
```
```Python

varsim.splice_site("INS")
```
```python
['NC_000011.10(NM_000207.3):c.187+1G>A', 'NC_000011.10(NM_000207.3):c.187+1G>T', 'NC_000011.10(NM_000207.3):c.187+1G>C', 'NC_000011.10(NM_000207.3):c.187+2T>A', 'NC_000011.10(NM_000207.3):c.187+2T>G', 'NC_000011.10(NM_000207.3):c.187+2T>C', 'NC_000011.10(NM_000207.3):c.188-2A>G', 'NC_000011.10(NM_000207.3):c.188-2A>T', 'NC_000011.10(NM_000207.3):c.188-2A>C', 'NC_000011.10(NM_000207.3):c.188-1G>A', 'NC_000011.10(NM_000207.3):c.188-1G>T', 'NC_000011.10(NM_000207.3):c.188-1G>C']
```