# VarSim
VarSim generates simulations for all possible single nucleotide variants (SNVs) for the Matched Annotation from NCBI and EMBL-EBI (MANE) transcript, along with the corresponding protein using HGVS notation.

**INSTALLATION**
```powershell
pip install varsim
```
Variation Simulation

```python
cds("INS")
```
NM_000207.3:c.1A>G	NP_000198.1:p.M1V	NP_000198.1:p.Met1Val

NM_000207.3:c.2T>G	NP_000198.1:p.M1R	NP_000198.1:p.Met1Arg

M_000207.3:c.2T>A	NP_000198.1:p.M1K	NP_000198.1:p.Met1Lys

NM_000207.3:c.3G>A	NP_000198.1:p.M1I	NP_000198.1:p.Met1Ile

NM_000207.3:c.1A>T	NP_000198.1:p.M1L	NP_000198.1:p.Met1Leu

...	...	...	...

NM_000207.3:c.328A>T	NP_000198.1:p.N110Y	NP_000198.1:p.Asn110Tyr

NM_000207.3:c.329A>T	NP_000198.1:p.N110I	NP_000198.1:p.Asn110Ile

NM_000207.3:c.330C>T	NP_000198.1:p.N110=	NP_000198.1:p.Asn110=

NM_000207.3:c.328A>C	NP_000198.1:p.N110H	NP_000198.1:p.Asn110His

NM_000207.3:c.329A>C	NP_000198.1:p.N110T	NP_000198.1:p.Asn110Thr
```python
utr5("INS")
```
NM_000207.3:c.-59A>G

NM_000207.3:c.-59A>T

NM_000207.3:c.-59A>C

NM_000207.3:c.-58G>A

NM_000207.3:c.-58G>T

...	...	...	...

NM_000207.3:c.-2C>A

NM_000207.3:c.-2C>T

NM_000207.3:c.-1C>G

NM_000207.3:c.-1C>A

NM_000207.3:c.-1C>T
```python
utr3("INS")
```
NM_000207.3:c.*1A>G

NM_000207.3:c.*1A>T

NM_000207.3:c.*1A>C

NM_000207.3:c.*2C>G

NM_000207.3:c.*2C>A

...	...	...	...

NM_000207.3:c.*72G>T

NM_000207.3:c.*72G>C

NM_000207.3:c.*73C>G

NM_000207.3:c.*73C>A

NM_000207.3:c.*73C>T
```python
splicing("INS")
```
NM_000207.3:c.187+2T>

NM_000207.3:c.188-2A>G

NM_000207.3:c.187+1G>A

NM_000207.3:c.187+2T>A

NM_000207.3:c.188-1G>A

NM_000207.3:c.187+1G>T

NM_000207.3:c.188-2A>T

NM_000207.3:c.188-1G>T

NM_000207.3:c.187+1G>C

NM_000207.3:c.187+2T>C

NM_000207.3:c.188-2A>C

NM_000207.3:c.188-1G>C