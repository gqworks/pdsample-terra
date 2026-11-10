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
| `rank` | Order inside the sample. Solved recessive statuses come first |
| `gene`, `moi` | Panel symbol after aliasing, and AD / AR / XL / RISK / UNKNOWN |
| `mechanism`, `mechanism_notes` | AD GoF/LoF/DN triage from the panel CSV (`gof`, `lof`, `dn`, `mixed`, `unknown`). Used to demote intron hits in GoF/DN genes to `AD_gof_intronic` |
| `variant_id` | `chrom-pos-ref-alt` with no `chr` prefix |
| `hgvsc`, `hgvsp` | VEP HGVS for the picked transcript, written next to `variant_id`. Coding change is `transcript:c.` and protein change is `transcript:p.` |
| `svtype`, `svlen` | Structural-variant type and length. `svlen` is INFO/SVLEN, sign kept, so a spec-compliant deletion is negative. When the tag is missing, the value is `end - pos` |
| `consequence`, `zygosity`, `genotype` | VEP consequence and the carried genotype. Male homozygous X small variants are hemizygous |
| `tier`, `evidence_score`, `evidence_reasons` | Point triage. Tier 1 is the highest score band |
| `inheritance_status` | `AR_comphet_candidate`, `AR_hom`, `AD_candidate`, `AD_gof_intronic` (weak intron in GoF/DN AD gene), `risk_biallelic`, `mtdna_pathogenic_candidate`, `mtdna_variant`, `str_expanded`, and the other statuses in `src/pdvar/report.py` |
| `hit_source` | `small_variant`, `expansionhunter`, `gauchian`, `mtdna`, or `sv` |
| `sv_panel_af` | True GATK-SV reference-panel AF when available (hard-filtered like gnomAD); not the joint-callset INFO/AF |
| `cohort_af` | In-callset SV allele frequency from INFO/AF or CN_NONREF_FREQ (annotation only; not a hard filter) |
| `clinvar_class` | `P/LP`, `B/LB`, `conflicting`, `VUS`, `other`, or `none` (from dbNSFP `clinvar_clnsig`) |
| `gnomad_af`, `gnomad_popmax_af` | Nuclear: from dbNSFP `gnomAD4.1_joint_AF` / `gnomAD4.1_joint_POPMAX_AF` (lifted to INFO for Introme). mtDNA: `gnomad_af` is MitImpact `AF_hom_gnomad` (gnomAD v3.1). Missing means rare |
| `gnomad_af_het` | mtDNA only: MitImpact `AF_het_gnomad` |
| `mitomap_status`, `mitomap_gb_freq`, `apogee2_class`, `apogee2_score`, `mitotip_score`, `mitotip_quartile`, `mlc_score`, `mitimpact_hgvs` | MitImpact annotations on kept `hit_source=mtdna` rows. `mlc_score` ≥ 0.5 adds `mlc_high`. Rows with heteroplasmy below 1% or a blank `MITOMAP_Disease_GenBank_Freq` are dropped before scoring |
| `cadd`, `revel` | From academic dbNSFP 5.4a (`CADD_phred` / `REVEL_score`) or optional CADD/REVEL plugins |
| `sift`, `sift_pred`, `polyphen`, `polyphen_pred` | From dbNSFP (`SIFT_*` / `Polyphen2_HDIV_*`); optional VEP cache `--sift/--polyphen` if `vep_sift_polyphen` is enabled |
| `eigen`, `muttaster`, `muttaster_pred`, `fathmm`, `fathmm_pred` | From dbNSFP when set; 5.4a uses `fathmm-XF_coding_*` (maps to `fathmm` / `fathmm_pred`) |
| `eigen_high` / `sift_deleterious` / … | Evidence-reason flags when configured score/prediction thresholds hit |
| `avi_phred`, `avi_raw`, `avi_quantile`, `avi_url` | AlphaGenome AVI for SNVs from the precomputed table when `avi_table` was set. `avi_raw` and `avi_phred` are copied from the table. `avi_quantile` is derived as `1 - 10^(-avi_phred/10)`, so it is only as precise as the published Phred. Phred 20 is the top 1% of genome-wide SNVs. Empty when `avi_table` was unset, the allele is not an SNV, or the table had no row |
| `introme_score`, `introme_class` | Introme v2 score and high / moderate / low |
| `splice_rule`, `splice_subgroup`, `splice_sav_pct`, `splice_heuristic_n` | Sullivan et al. 2025 donor (Figure 3) or acceptor (Figure 4) splice-disruption call for an SNV in those windows. `splice_rule` is `DD4` / `DA2` and the like, `splice_subgroup` is `standard`, `context`, or `auxiliary`, and the percent and n are the published spliceogenicity. Blank for indels, alleles outside the windows, and when the reference FASTA or GTF was not available. No evidence points |
| `str_locus`, `str_expanded` | ExpansionHunter locus. `str_expanded` is true when `max(REPCN)` meets the cited `str_loci` floor |
| `heteroplasmy` | mtDNA allele fraction from FORMAT `AF`. Alleles below 1% are not ingested |
| `age_onset`, `symptoms`, `ethnicity` | Patient fields parsed from the optional `phenotype_report` text when supplied |
| `onset_match`, `matched_features`, `phenotype_score` | Concordance against `config/gene_phenotypes.yaml`. Onset in range is +2, outside is -2, and each matched feature or enriched ancestry is +1 |
| `clinical_match` | `matched`, `unmatched`, `uninformative`, `no_profile`, or `not_assessed`. Risk rows are repeated under Phenotype concordance in the report |
| `levodopa`, `age_sampling`, `comorbidities` | From the phenotype report when one was supplied. Comorbidities are listed in the report. Levodopa use is a biomarker note |

`qc_summary.json` `modules` records `called` or `not_assessed` for small variants, ExpansionHunter, Gauchian, mtDNA, Introme, VEP, AVI, `splice_heuristics`, and `mtdna_annotation` (MitImpact join). `splice_heuristics` is `called` when the reference FASTA and GENCODE GTF were both available. AVI is `called` only when the Atlas table has at least one score. GATK-SV is `called` when the SV VCF came from `GATKSVPipelineSingleSample` (including fallback when `sample_id` is absent from optional `sv_joint_vcf`), or `joint` when it was subset from that joint VCF. A missing `sv_vcf` raises `FileNotFoundError`. `sv_status` on the QC object mirrors that value. The `sv` array still uses `no SV`, `SV carrier`, or `SV not assessed` for whether that sample is in the VCF header and carries a panel SV. `relatedness` is `not_run`.
