# Data dictionary

Outputs are one sample. `candidates.tsv` drops rows whose `inheritance_status` is `not_qualifying`.

## Files

| File | Contents |
| --- | --- |
| `candidates.tsv` | Ranked variants for the Terra row |
| `report.md`, `report.html` | The same triage narrative |
| `qc_summary.json` | Contamination, sex, coverage, reference check, and per-module status |
| `liability.tsv` | `solved`, `partially_explained`, or `unexplained`, plus a null explained-liability until penetrance is cited |
| `followup.tsv` | Recessive genes that still have a single heterozygous hit |
| `sample.json` | Paths the Assess task passed in |

## Candidate columns

| Column | Meaning |
| --- | --- |
| `rank` | Order inside the sample. Solved recessive statuses come before incidental |
| `gene`, `moi` | Panel symbol after aliasing, and AD / AR / XL / RISK / UNKNOWN |
| `variant_id` | `chrom-pos-ref-alt` with no `chr` prefix |
| `consequence`, `zygosity`, `genotype` | VEP consequence and the carried genotype. Male homozygous X small variants are hemizygous |
| `tier`, `evidence_score`, `evidence_reasons` | Point triage. Tier 1 is the highest score band |
| `inheritance_status` | `AR_comphet_candidate`, `AR_hom`, `AD_candidate`, `risk_biallelic`, `mtdna_variant`, `str_expanded`, `incidental`, and the other statuses in `src/pdvar/report.py` |
| `hit_source` | `small_variant`, `expansionhunter`, `gauchian`, `mtdna`, `sv`, or `introme_lookup` |
| `gnomad_af`, `gnomad_popmax_af` | From the VEP custom annotation when present. Missing means rare |
| `sv_panel_af` | Allele frequency on the GATK-SV reference panel, when that VCF ran |
| `clinvar_class` | `P/LP`, `B/LB`, `conflicting`, `VUS`, `other`, or `none` |
| `cadd`, `revel`, `spliceai` | Present when VEP plugins or Introme wrote them |
| `introme_score`, `introme_class` | Introme v2 score and high / moderate / low |
| `str_locus`, `str_expanded` | ExpansionHunter locus. `str_expanded` is null until a cited repeat cutoff exists |
| `heteroplasmy` | mtDNA allele fraction when the caller wrote `AF` |
| `incidental`, `disposition` | `TTN_truncation` / `incidental` for TTN null alleles. ClinVar is not overwritten |
| `lof_weight` | 1 while expression cutoffs are null |

`qc_summary.json` `modules` records `called` or `not_assessed` for small variants, ExpansionHunter, Gauchian, mtDNA, Introme, and VEP. GATK-SV is `called` when the SV VCF is present. A missing `sv_vcf` raises `FileNotFoundError`. `sv_status` on the QC object is `called`. The `sv` array still uses `no SV`, `SV carrier`, or `SV not assessed` for whether that sample is in the VCF header and carries a panel SV. `relatedness` is `not_run`.
