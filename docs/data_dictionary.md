# Data dictionary

Outputs are one sample. `candidates.tsv` drops rows whose `inheritance_status` is `not_qualifying` unless `followup_reason` is set. `report.html` omits `not_qualifying` rows even when they are in the TSV. How ranking and shortlisting work: [README § Variant prioritisation](../README.md#variant-prioritisation).

## Files

Terra writes these Assess files into one `candidates` array column, in this order: `{sample_id}_report.html`, `{sample_id}_candidates.tsv`, `{sample_id}_summary.tsv`. Sample-level genotypes and scores are the summary TSV only (not extra shortlist rows). Unique QC JSON fields (sex, contamination, coverage, reference, relatedness, pass/failures) are copied into the HTML report.

| File | Contents |
| --- | --- |
| `{sample_id}_report.html` | Per-sample triage report. Prints sample classification (`solved` / `partially_explained` / `unexplained`), explained-liability text, QC (sex, FREEMIX, coverage, reference, relatedness, pass/failures), modules, and the ranked table. Ranked candidates (in order): Gene name, Inheritance status, CPRA (`variant_id`), consequence, variant type, Zygosity, hgvsc, hgvsp, ClinVar, gnomad AF, gnomad popmax AF, introme score, AVI (`avi_phred`), heteroplasmy, mitomap_clinical_info, Spliceogenicity (`splice_sav_pct`), SV length, clinical_match, phenotype score, tier, evidence score, evidence reason. Omits `not_qualifying` rows |
| `{sample_id}_candidates.tsv` | Ranked shortlist for the Terra row. AR single-het / review-conflict follow-up alleles, unphased comp-hets, rescued second hits (including `not_qualifying`), and empty-shortlist module stubs are force-kept (after Mendelian solves) and carry `followup_reason` |
| `{sample_id}_summary.tsv` | One row per sample: APOE ε genotype, m.2158T>C, Haplogrep haplogroup, PGS000902 PRS (90 SNPs), mean MLC, MSS |

## Candidate columns

| Column | Meaning |
| --- | --- |
| `rank` | Order inside the sample after inheritance priority, ClinVar preference, tier, score, AVI, and phenotype tie-breaks. Mendelian / AD / risk / high-evidence mtDNA statuses lead; `AR_single_het` follow-up rows sort lower |
| `gene`, `moi` | Panel symbol after aliasing, and AD / AR / XL / RISK / UNKNOWN |
| `mechanism`, `mechanism_notes` | AD GoF/LoF/DN triage from the panel CSV (`gof`, `lof`, `dn`, `mixed`, `splice_gof`, `unknown`). Used to demote intron hits in GoF/DN genes to `AD_gof_intronic`; MAPT `splice_gof` skips that demotion |
| `variant_id` | `chrom-pos-ref-alt` with no `chr` prefix |
| `hgvsc`, `hgvsp` | VEP HGVS for the picked transcript, written next to `variant_id`. Coding change is `transcript:c.` and protein change is `transcript:p.` |
| `svtype`, `svlen` | Structural-variant type and length. `svlen` is INFO/SVLEN, sign kept, so a spec-compliant deletion is negative. When the tag is missing, the value is `end - pos` |
| `copy_number`, `expected_cn` | Sample copy number from FORMAT `CN`, else GATK-SV `RD_CN`, and expected copy number (`ECN`). Used by `pesr_caller_qc` and dosage scoring |
| `consequence`, `zygosity`, `genotype` | VEP consequence and the carried genotype. Male homozygous X small variants are hemizygous |
| `tier`, `evidence_score`, `evidence_reasons` | Point triage. Tier 1 is the highest score band |
| `inheritance_status` | `AR_comphet_candidate`, `AR_hom`, `AD_candidate`, `AD_gof_intronic` (weak intron in GoF/DN AD gene), `risk_biallelic`, `mtdna_pathogenic_candidate`, `mtdna_variant`, `str_expanded`, `module_not_assessed`, and the other statuses in `src/pdvar/report.py` |
| `str_n_expanded_alleles`, `str_motif`, `str_motif_pathogenic`, `str_moi` | Count of alleles at or above the cited floor, observed motif when parsed, whether that motif is the pathogenic one, and STR inheritance |
| `high_sr_background`, `freq_het`, `sv_artifact` | GATK-SV HIGH_SR_BACKGROUND flag, in-cohort FREQ_HET, and depth-only artifact tag (unscored; not a hard drop) |
| `hit_source` | `small_variant`, `expansionhunter`, `gauchian`, `mtdna`, or `sv` |
| `sv_panel_af` | True GATK-SV reference-panel AF when available (hard-filtered like gnomAD); not the joint-callset INFO/AF |
| `cohort_af` | In-callset SV allele frequency from INFO/AF or CN_NONREF_FREQ (annotation only; not a hard filter) |
| `sv_algorithms` | INFO/ALGORITHMS from GATK-SV |
| `sv_evidence` | FORMAT EV evidence classes (RD / PE / SR) |
| `rd_gq` | FORMAT RD_GQ when present |
| `ogq` | FORMAT OGQ when present |
| `clinvar_class` | `P/LP`, `B/LB`, `conflicting`, `VUS`, `other`, or `none` (from panel ClinVar P/LP `ClinVar_CLNSIG` or dbNSFP `clinvar_clnsig`) |
| `gnomad_af`, `gnomad_popmax_af` | Nuclear: prefer panel-sliced gnomAD v4.1 joint sites (`AF_joint` / `AF_grpmax_joint` → INFO); else dbNSFP `gnomAD4.1_joint_*` via CSQ lift. mtDNA: `gnomad_af` is MitImpact `AF_hom_gnomad` (gnomAD v3.1). Missing means rare |
| `gnomad_af_het` | mtDNA only: MitImpact `AF_het_gnomad` |
| `mitomap_status`, `mitomap_clinical_info`, `mitomap_gb_freq`, `apogee2_class`, `apogee2_score`, `mitotip_score`, `mitotip_quartile`, `mlc_score`, `mitimpact_hgvs` | MitImpact annotations on kept `hit_source=mtdna` rows. `mitomap_clinical_info` is `MITOMAP_Disease_Clinical_info`. `mlc_score` ≥ 0.5 adds `mlc_high`. Rows with heteroplasmy below 1% or a blank `MITOMAP_Disease_GenBank_Freq` are dropped before scoring |
| `cadd`, `revel` | From academic dbNSFP 5.4a (`CADD_phred` / `REVEL_score`) or optional CADD/REVEL plugins |
| `sift`, `sift_pred`, `polyphen`, `polyphen_pred` | From dbNSFP (`SIFT_*` / `Polyphen2_HDIV_*`); optional VEP cache `--sift/--polyphen` if `vep_sift_polyphen` is enabled |
| `eigen`, `muttaster`, `muttaster_pred`, `fathmm`, `fathmm_pred` | From dbNSFP when set; 5.4a uses `fathmm-XF_coding_*` (maps to `fathmm` / `fathmm_pred`) |
| `eigen_high` / `sift_deleterious` / … | Evidence-reason flags when configured score/prediction thresholds hit |
| `avi_phred`, `avi_raw`, `avi_quantile`, `avi_url` | AlphaGenome AVI for SNVs from the precomputed table when `avi_table` was set. `avi_raw` and `avi_phred` are copied from the table. `avi_quantile` is derived as `1 - 10^(-avi_phred/10)`, so it is only as precise as the published Phred. Phred 20 is the top 1% of genome-wide SNVs. Empty when `avi_table` was unset, the allele is not an SNV, or the table had no row |
| `introme_score`, `introme_class` | Introme v2 score and high / moderate / low. High is ≥ 0.9 and takes `introme_high` points. Gene-scoped second-hit fallback at 0.7 is moderate |
| `splice_rule`, `splice_subgroup`, `splice_sav_pct`, `splice_heuristic_n` | Sullivan et al. 2025 donor (Figure 3) or acceptor (Figure 4) splice-disruption call for an SNV in those windows. `splice_rule` is `DD4` / `DA2` and the like, `splice_subgroup` is `standard`, `context`, or `auxiliary`, and the percent and n are the published spliceogenicity. Blank for indels, alleles outside the windows, and when the reference FASTA or GTF was not available. No evidence points |
| `str_locus`, `str_expanded` | ExpansionHunter locus. `str_expanded` is true when `max(REPCN)` meets the cited `str_loci` floor |
| `heteroplasmy` | mtDNA allele fraction from FORMAT `AF`. Alleles below 1% are not ingested |
| `age_onset`, `symptoms`, `ethnicity` | Patient fields parsed from the optional `phenotype_report` text when supplied |
| `onset_match`, `matched_features`, `phenotype_score` | Concordance against `config/gene_phenotypes.yaml`. Onset in range is +2, outside is -2, and each matched feature or enriched ancestry is +1 |
| `clinical_match` | `matched`, `unmatched`, `uninformative`, `no_profile`, or `not_assessed`. Risk rows are repeated under Phenotype concordance in the report |
| `levodopa`, `age_sampling`, `comorbidities` | From the phenotype report when one was supplied. Comorbidities are listed in the report. Levodopa use is a biomarker note |
| `followup_reason` | Last column. Non-empty for alleles in an AR gene that still has only a single heterozygous hit (or a review_conflict), unphased comp-het/hemizygous pairs, rescued second-hit search rows (including Introme 0.7 fallback), early-onset recessive nulls, incomplete AR STR expansions, and empty-shortlist module stubs. Used to force-keep those rows |

## QC in `{sample_id}_report.html`

These fields used to live only in `qc_summary.json`. They are copied into the HTML report so that file is not a Terra output:

| Field | Meaning |
| --- | --- |
| `pass` / `failures` | QcGate hard-fail result. Failures list FREEMIX, coverage, or sex-check messages when those gates trip |
| `sex` | `value`, `source` (`provided` / `somalier` / `unknown`), `somalier`, and `provided` |
| `contamination` | VerifyBamID2 `FREEMIX`, the configured max, and whether it passed |
| `coverage` | mosdepth mean depth, the configured min, and whether it passed |
| `reference` | Whether BAM/CRAM contigs match the reference, plus contig counts |
| `relatedness` | Always `not_run` in this per-sample workflow |

Any other QC JSON key that is not already in Modules, Callsets, or Phenotype is printed as `other`.

## Sample summary (`{sample_id}_summary.tsv`)

These fields describe the sample. They are not extra `candidates.tsv` rows.

| Field | Meaning |
| --- | --- |
| `apoe_genotype` | Unphased ε2/ε3/ε4 from rs429358 (GRCh38 `chr19:44908684` T>C, ε4) and rs7412 (`chr19:44908822` C>T, ε2). `e2/e4_or_e1/e3` when the two SNPs are both heterozygous (phase unknown). APOE is **not** on the 133-gene Mendelian panel; `CallSummaryLoci` genotypes these two sites from the BAM. Citation: NCBI dbSNP / standard APOE ε2/ε3/ε4 definition |
| `apoe_rs429358_gt` / `apoe_rs7412_gt` | Sorted allele strings at those sites (`C/T`, `T/T`, …) |
| `apoe_status` | `called`, `partial` (one SNP missing), or `not_assessed` |
| `mt_2158_genotype` | m.2158T>C (SHLP2 p.Lys4Arg, rs41349444) from the mt VCF. `T` = reference (variant absent after a successful mito call), `C` = heteroplasmy ≥ 0.95, `T/C` = 0.01–0.95. Hudson et al. *Neurology* 2013; Kim et al. *Mol Psychiatry* 2024 |
| `mt_2158_heteroplasmy` | FORMAT `AF` at m.2158, or `0` when the allele is absent |
| `mt_2158_status` | `called` or `not_assessed` |
| `haplogroup` | Haplogrep 3 best hit on the mt VCF (`phylotree-rcrs@17.2`). Container: `quay.io/biocontainers/haplogrep3:3.2.2--hdfd78af_1`. Phylotree 17 rCRS is bundled in that image; `install-tree` is a no-op fallback |
| `haplogroup_quality` | Haplogrep quality score |
| `haplogroup_status` | `called` or `not_assessed` |
| `panel_prs` | Weighted allele-count sum of **all 90** PGS Catalog **PGS000902** SNPs (Nalls et al. *Lancet Neurol* 2019). `CallSummaryLoci` genotypes every site from the WGS BAM (hom-ref included). Missing sites are not imputed as hom-ref. `score_in_panel_prs` in the weights TSV is annotation only. Weights: `resources/prs/pgs000902_nalls2019.tsv` |
| `prs_sites_used` / `prs_sites_missing` / `prs_sites_total` | How complete that 90-SNP score is (`prs_sites_total` is 90 when the shipped weights file is used) |
| `prs_status` | `called` when at least one PGS SNP was genotyped; `not_assessed` otherwise |
| `mean_mlc` | Mean MitImpact `MI_MLC_score` across mt SNVs in the annotated mt VCF, **excluding this sample's Haplogrep FoundPolys** (haplogroup-defining / haplotagging variants present in that sample). Lake et al. *Nature* 2024 via MitImpact. There is no static in-repo haplotagging list |
| `mss` | MLC score sum (MSS; Lake et al. *Nature* 2024) of MitImpact `MI_MLC_score` over the **same** mt SNVs as `mean_mlc` (FoundPolys excluded). Lake's UKB MSS summed heteroplasmies; this column uses the mito-callset SNVs after haplotagging exclusion |
| `mlc_n_variants` / `mlc_n_excluded_haplotag` | SNVs considered and haplogroup-defining sites dropped |
| `mlc_haplotag_exclusion` | `haplogrep_found_polys`, `empty_found_polys`, or `not_applied` |
| `mlc_status` | `called` when at least one scored SNV remains, else `not_assessed` |

The HTML report QC/modules section records `called` or `not_assessed` for small variants, ExpansionHunter, Gauchian, mtDNA, Introme, VEP, AVI, `splice_heuristics`, `mtdna_annotation` (MitImpact join), `haplogrep`, `apoe`, `mt_2158`, `panel_prs`, and `mean_mlc` (covers MSS). `splice_heuristics` is `called` when the reference FASTA and GENCODE GTF were both available. AVI is `called` only when the Atlas table has at least one score. GATK-SV is `called` when the SV VCF came from `GATKSVPipelineSingleSample` (including fallback when `sample_id` is absent from optional `sv_joint_vcf`), `joint` when it was subset from that joint VCF, or `prior` when `this.sv_vcf` was reused. Introme is `prior` when a real non-empty `this.introme` was reused (empty/placeholder cells call Introme). A real non-empty `this.vep_vcf` skips the VEP annotation chain; Assess still runs. The same `ProbePriorFile` failsafe applies to `this.small_vcf`, `this.sv_vcf`, `this.expansionhunter_vcf`, `this.gauchian_json`, `this.mtdna_vcf`, and `this.avi_tsv`: empty/missing/placeholder/zero-byte cells run that producer. `mtdna_annotation` is `prior` when `this.mtdna_vcf` was reused. A missing `sv_vcf` raises `FileNotFoundError`. `sv_status` on the QC object mirrors that value. The `sv` array still uses `no SV`, `SV carrier`, or `SV not assessed` for whether that sample is in the VCF header and carries a panel SV. `relatedness` is `not_run`.
