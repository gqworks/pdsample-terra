# Parkinson's PRS weights (PGS000902)

`pgs000902_nalls2019.tsv` is the PGS Catalog scoring file **PGS000902** (PRS90_PD)
harmonized to GRCh38.

| Field | Value |
| --- | --- |
| PGS Catalog | [PGS000902](https://www.pgscatalog.org/score/PGS000902/) |
| Publication | Nalls MA et al. *Lancet Neurol.* 2019;18:1091-1102. doi:10.1016/S1474-4422(19)30320-5 |
| PGS Catalog publication | PGP000235 |
| Weight type | beta (log-odds) |
| Variants | 90 genome-wide significant PD SNPs |

`pdvar` scores **`panel_prs`** as the weighted allele-count sum over **all 90
rows**. `CallSummaryLoci` genotypes those sites (plus APOE ε SNPs) from the WGS
BAM, including hom-ref. Missing genotypes are not imputed as hom-ref.

`score_in_panel_prs` is annotation only: `true` for SNPs inside a MonoPDAus
gene span or within 20 kb (SNCA cis tags rs356182 / rs5019538, plus COASY cis
rs12951632). It does not filter the default score.

Do not resurrect an empty weights stub. Every row carries the citation above.
The two APOE ε SNPs are **not** in PGS000902; they are genotyped separately.
