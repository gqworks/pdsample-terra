# Data dictionary

Outputs are one sample. `candidates.tsv` drops rows whose `inheritance_status` is `not_qualifying`.

## Files

Terra writes these Assess files into one `candidates` array column, in this order: `report.html`, `qc_summary.json`, `liability.tsv`, `followup.tsv`, `candidates.tsv`.

| File | Contents |
| --- | --- |
| `candidates.tsv` | Ranked variants for the Terra row |
| `report.html` | Per-sample triage report |
| `qc_summary.json` | Contamination, sex, coverage, reference check, and per-module status |
| `liability.tsv` | `solved`, `partially_explained`, or `unexplained`, plus a null explained-liability until penetrance is cited |
| `followup.tsv` | Recessive genes that still have a single heterozygous hit |

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
| `sv_panel_af` | Allele frequency on the GATK-SV reference panel, when that VCF ran |
| `clinvar_class` | `P/LP`, `B/LB`, `conflicting`, `VUS`, `other`, or `none` (from dbNSFP `clinvar_clnsig`) |
| `gnomad_af`, `gnomad_popmax_af` | From dbNSFP `gnomAD4.1_joint_AF` / `gnomAD4.1_joint_POPMAX_AF` (lifted to INFO for Introme). Missing means rare |
| `cadd`, `revel` | From academic dbNSFP 5.4a (`CADD_phred` / `REVEL_score`) or optional CADD/REVEL plugins |
| `sift`, `sift_pred`, `polyphen`, `polyphen_pred` | From dbNSFP (`SIFT_*` / `Polyphen2_HDIV_*`); optional VEP cache `--sift/--polyphen` if `vep_sift_polyphen` is enabled |
| `eigen`, `muttaster`, `muttaster_pred`, `fathmm`, `fathmm_pred` | From dbNSFP when set; 5.4a uses `fathmm-XF_coding_*` (maps to `fathmm` / `fathmm_pred`) |
| `eigen_high` / `sift_deleterious` / … | Evidence-reason flags when configured score/prediction thresholds hit |
| `avi_phred`, `avi_raw`, `avi_quantile`, `avi_url` | AlphaGenome AVI for SNVs from the precomputed table when `avi_table` was set. `avi_raw` and `avi_phred` are copied from the table. `avi_quantile` is derived as `1 - 10^(-avi_phred/10)`, so it is only as precise as the published Phred. Phred 20 is the top 1% of genome-wide SNVs. Empty when `avi_table` was unset, the allele is not an SNV, or the table had no row |
| `introme_score`, `introme_class` | Introme v2 score and high / moderate / low |
| `str_locus`, `str_expanded` | ExpansionHunter locus. `str_expanded` is null until a cited repeat cutoff exists |
| `heteroplasmy` | mtDNA allele fraction when the caller wrote `AF` |
| `incidental`, `disposition` | `TTN_truncation` / `incidental` for TTN null alleles. ClinVar is not overwritten |
| `lof_weight` | 1 while expression cutoffs are null |
| `age_onset`, `symptoms`, `ethnicity` | Patient fields copied from the sample row when they were provided |
| `onset_match`, `matched_features`, `phenotype_score` | Concordance against `config/gene_phenotypes.yaml`. Onset in range is +2, outside is -2, and each matched feature or enriched ancestry is +1 |
| `clinical_match` | `matched`, `unmatched`, `uninformative`, `no_profile`, or `not_assessed`. Risk rows are repeated under Phenotype concordance in the report |
| `levodopa`, `age_sampling`, `comorbidities` | From a phenotype report when one was supplied. Comorbidities are listed in the report. Levodopa use is a biomarker note |

`qc_summary.json` `modules` records `called` or `not_assessed` for small variants, ExpansionHunter, Gauchian, mtDNA, Introme, VEP, and AVI. AVI is `called` only when the Atlas table has at least one score. GATK-SV is `called` when the SV VCF came from `GATKSVPipelineSingleSample`, or `joint` when it was subset from optional `sv_joint_vcf`. A missing `sv_vcf` raises `FileNotFoundError`. `sv_status` on the QC object mirrors that value. The `sv` array still uses `no SV`, `SV carrier`, or `SV not assessed` for whether that sample is in the VCF header and carries a panel SV. `relatedness` is `not_run`.
