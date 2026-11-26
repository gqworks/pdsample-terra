# GBA1 resources

## GBA1-PD Browser lookup

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
joins this table onto **Gauchian** and **gba1_paralog** GBA1 SNVs (and VEP
GBA1 rows) matched by c. HGVS and applies a configurable score bump
(`scoring.points.gba1_pd_*`). It is triage annotation, not an ACMG
classification.

## Homology hotspot force-call

`gba1_homology_discriminators.tsv`, `gba1_homology_hotspots.tsv`, and
`gbap1_grch38_patches.tsv` drive `workflow/scripts/gba1_paralog_forcecall.py`
(WDL `Gba1ParalogHotspots`). Window reads are scored against the GRCh38 GBA1
template and a GBAP1 template with GRCh38’s GBA-like minor alleles patched to
the major GBAP1 bases (Orimo et al. 2024 / Woo et al. 2022). Gene-better
fragments are piled at L444P / A495P / D448H so mis-mapping to GBAP1 does not
silently drop those alleles when Gauchian returns a non-carrier.
