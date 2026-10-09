# Methods

`PdSample` interprets one genome. Every call is made from that sample's CRAM. Scores are additive points from `config/config.yaml`. They are a triage rank in the style of ACMG evidence categories, and they are not an ACMG/AMP classification.

## Callers

Small variants are single-sample HaplotypeCaller GVCFs, genotyped with GenotypeGVCFs and filtered with the GATK single-sample hard-filter expressions documented in the README. VQSR is a cohort fit, and no CNN model is bundled, so hard filters are the filter used here. The default interval list is the panel gene-span BED.

GATK-SV single-sample (v0.29-beta) runs on every sample. The reference panel and hg38 resource files are required inputs, taken from the public 1KG bundle that tag defines or from a matched panel in workspace data. Overlap uses the gene and exon BEDs under `resources/panel/`. Reference-panel allele frequency on an SV is stored as `sv_panel_af` and filtered with the same ceiling as gnomAD AF. `pdvar assess` requires the SV VCF; it does not treat a missing file as a skipped call. A sample that is absent from an existing SV VCF header is still `SV not assessed` on that row of the coverage table.

ExpansionHunter uses the vendored 19-locus catalog. A repeat count is a disease call only when `str_loci.<locus>.pathogenic_min_repeats` and `citation` are both set. Until then the row can be kept as an STR carrier and is not labelled expanded.

Gauchian copy number 0, or an explicit deletion, becomes a homozygous `transcript_ablation` in `GBA1`. A recombinant allele stays a structural variant unless the Gauchian file itself marks it pathogenic.

MitochondriaPipeline 4.7.0.0 calls chrM. Carried mtDNA variants are labelled `mtdna_variant`. That status is visible on the candidate list and does not mark the sample solved.

## Annotation and ranking

VEP 116 `--pick` supplies consequence, symbol, and HGVS. Optional custom gnomAD and ClinVar INFO fields override CSQ frequency and clinical significance. CADD, REVEL, and SpliceAI are read from CSQ when a plugin added them. Missing gnomAD AF is treated as rare.

Introme v2 (`d6da48bc`) scores splice consequences on the panel BED. When VEP wrote a gnomAD custom annotation, `INFO/gnomAD_AF` is available to Introme's rare filter. A second Introme allele is added inside the same sample for an AR gene that already has a qualifying hit, and only when `gnomad_af` is missing or at or below `lookup_max_gnomad_af` (0.01).

Inheritance uses the panel CSV (`GBA` aliased to `GBA1`). Two qualifying heterozygous alleles in an AR gene are `AR_comphet_candidate` unless both genotypes are phased (`|`) and share a `PS` tag, in which case the pair is `AR_comphet_phased` or collapsed to a single het. Statistical phasing is disabled. A biallelic finding in a risk gene such as GBA1 is `risk_biallelic` and counts as partially explained. TTN null alleles are tagged `incidental` (`TTN_truncation`). ClinVar is left as reported, and the row is excluded from solved and partial liability.

Isoform and expression weights stay neutral while `expression.min_brain_pext`, `min_gtex_tpm`, and `unexpressed_weight` are null. Biomarker multipliers stay 1 while their citations are null. `explained_liability` is numeric only when `liability.penetrance` has a value and a citation. The shipped penetrance map is empty.
