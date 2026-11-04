# GBA1-PD Browser lookup

`gba1_pd_browser.tsv` is a single-variant extract of the GBA1-PD Browser export
(Parlar et al., *Mov Disord* 2023; https://pdgenetics.shinyapps.io/GBA1Browser/).

| Column | Meaning |
|--------|---------|
| `cdna` | HGVS c. notation (join key for Gauchian `other_variants`) |
| `protein` | Full-length protein HGVS when present |
| `pd_class` | `severe`, `mild`, `risk`, or `unknown` |
| `gd_acmg` | ACMG-style label in Gaucher disease from the browser |
| `pos_hg38` | GRCh38 position when provided |
| `allele_frequency` | Browser gnomAD AF |
| `references` / `comments` | Browser provenance |

Multi-allele / haplotype rows from the browser export are dropped. `pdvar`
joins this table onto **Gauchian** GBA1 SNVs only and applies a configurable
score bump (`scoring.points.gba1_pd_*`). It is triage annotation, not an ACMG
classification.
