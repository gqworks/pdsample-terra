# AD gene mechanism sources (GoF vs LoF)

Labels on `monopdaus-gene-inheritance.csv` (`mechanism`, `mechanism_notes`) are a
literature triage synthesis for whether **splice / intronic LoF-like alleles** are
expected to explain autosomal-dominant disease in that gene.

| Label | Meaning for Introme / intron triage |
| --- | --- |
| `gof` | Toxic or enzymatic gain of function; splice LoF is not the Mendelian mechanism |
| `dn` | Dominant-negative missense (or equivalent); splice LoF usually does not recapitulate disease |
| `lof` | Haploinsufficiency / LoF; splice-disrupting alleles can be relevant |
| `mixed` | Both GoF and LoF evidence; keep splice relevance (conservative) |
| `unknown` | Insufficient evidence |

Pipeline behaviour: when `inheritance.demote_gof_ad_intronic` is true, weak intron /
Introme-only hits in AD genes labelled `gof` or `dn` are demoted to
`AD_gof_intronic` (below Mendelian shortlist priority). ClinVar P/LP and coding
alleles are not demoted. `lof` and `mixed` keep the existing AD intron demotion
(priority 1) only.

## Core Mendelian PD

| Gene | Label | Key sources |
| --- | --- | --- |
| SNCA | gof | GeneReviews NBK1223; Singleton 2003 Science (duplication); Chartier-Harlin 2004 Lancet |
| LRRK2 | gof | GeneReviews NBK1208; Luzón-Toro 2007 Hum Mol Genet (G2019S kinase GoF) |
| VPS35 | dn | Trinh 2018 MDSGene Mov Disord; Sassone 2021 / Williams 2022 retromer DN reviews |
| CHCHD2 | mixed | Funayama 2015 Lancet Neurol; Frontiers Neurosci 2025 CHCHD2 review (LoF + T61I toxic GoF) |
| RAB32 | gof | Gustavsson 2024 Lancet Neurol (Ser71Arg activates LRRK2) |

## Atypical parkinsonism / dystonia differentials (selected)

| Gene | Label | Key sources |
| --- | --- | --- |
| DCTN1 | dn | Wittke 2021 MDSGene atypical PD; Perry syndrome CAP-Gly missense |
| ATP1A3 | gof | GeneReviews ATP1A3; Kaneko / ionic-leak GoF literature |
| ANO3 | gof | GeneReviews dystonia; Brain 2023 awad412 (channel/scramblase) |
| GCH1 | lof | GeneReviews DRD; OMIM 600225 (whole-gene deletions → HI) |
| KMT2B | lof | GeneReviews dystonia |
| TUBB4A | dn | GeneReviews DYT4 / H-ABC |
| EIF2AK2 | gof | Mao 2020 Am J Hum Genet |
| PRKCG | gof | GeneReviews SCA14 |
| PPP2R5D | gof | OMIM PPP2R5D / Jordan syndrome missense |

## Dementia / ALS / PFBC / mitochondrial (selected)

| Gene | Label | Key sources |
| --- | --- | --- |
| APP, PSEN1, PSEN2 | gof | GeneReviews Alzheimer Disease Overview |
| MAPT | gof | GeneReviews MAPT |
| GRN | lof | GeneReviews GRN-FTD |
| PRNP | gof | GeneReviews genetic prion |
| TBK1 | lof | Freischmidt 2015 Nat Neurosci |
| SQSTM1 | lof | GeneReviews ALS |
| FUS | mixed | GeneReviews ALS (toxic missense + some LoF) |
| CSF1R | lof | GeneReviews ALSP; Konno 2014 Neurology (HI / splice) |
| SLC20A2, PDGFB, PDGFRB, XPR1 | lof | Balck 2021 MDSGene PFBC; GeneReviews PFBC |
| OPA1 | lof | GeneReviews DOA; PMC10645107 (primarily HI) |
| AFG3L2 | lof | GeneReviews SCA28; haploinsufficiency reports |
| TWNK, SLC25A4 | dn | GeneReviews / OMIM adPEO dominant-negative |
| VCP | dn | GeneReviews VCP |
| FTL | gof | GeneReviews neuroferritinopathy |
| DNAJC5 | gof | GeneReviews NCL; PMC9809949 |
| STXBP1, SCN1A, FOXG1 | lof | GeneReviews (haploinsufficiency encephalopathies) |

This table is a starter synthesis for triage, not a systematic review. Expert review
is required before use beyond research ranking.
