# Sources for `gene_phenotypes.yaml`

Starter literature synthesis behind each profile in `config/gene_phenotypes.yaml`. It is not a
formal systematic review. Onset ranges and keywords were chosen for triage concordance and must
be reviewed by a clinician before they inform reporting.

Profile confidence:

- **high**: established gene-disease relationship with a GeneReviews / MDSGene / OMIM clinical
  description that gives onset and core features.
- **medium**: established or probable relationship, but few families or a broad phenotype, so
  onset and keywords are approximate.
- **low**: unreplicated, contested, or phenotype-poor. The profile is deliberately minimal
  (usually no `onset_range` and no keywords) so it cannot drive a match.

Citations were checked against live sources for the newer or thinly reported genes (RAB32, PSMF1,
DAGLB, PTPA, CHCHD2, ATP6AP2, the MDSGene atypical and PFBC reviews). Other GeneReviews and
OMIM pointers are standard references and should be spot-checked during clinical review.

Main references used across genes:

- GeneReviews Parkinson Disease Overview: https://www.ncbi.nlm.nih.gov/books/NBK1223/
- MDSGene reviews: https://www.mdsgene.org (Trinh 2018 doi:10.1002/mds.27527; Kasten 2018
  doi:10.1002/mds.27352; Wittke 2021 doi:10.1002/mds.28517; Balck 2021 doi:10.1002/mds.28753)
- OMIM clinical synopses: https://omim.org
- Lange et al. 2022, MDS Task Force nomenclature update (Mov Disord)


## Mendelian PD with a typical (levodopa-responsive) picture

| Gene | MOI | Onset | Confidence | Sources |
| --- | --- | --- | --- | --- |
| SNCA | AD | 20-65 | high | Trinh et al. 2018 MDSGene SNCA/LRRK2/VPS35 review, Mov Disord, doi:10.1002/mds.27527; GeneReviews Parkinson Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1223/) |
| LRRK2 | AD | 35-80 | high | GeneReviews LRRK2 Parkinson Disease (https://www.ncbi.nlm.nih.gov/books/NBK1208/); Updated MDSGene LRRK2 review 2025, npj Parkinsons Dis (s41531-025-00881-9); median onset ~56 y; Trinh et al. 2018 MDSGene SNCA/LRRK2/VPS35 review, Mov Disord, doi:10.1002/mds.27527 |
| VPS35 | AD | 30-70 | high | Trinh et al. 2018 MDSGene SNCA/LRRK2/VPS35 review, Mov Disord, doi:10.1002/mds.27527; GeneReviews Parkinson Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1223/) |
| PRKN | AR | 10-45 | high | Kasten et al. 2018 MDSGene Parkin/PINK1/DJ-1 review, Mov Disord, doi:10.1002/mds.27352; GeneReviews Parkin Type of Early-Onset PD (https://www.ncbi.nlm.nih.gov/books/NBK1478/) |
| PINK1 | AR | 20-50 | high | Kasten et al. 2018 MDSGene Parkin/PINK1/DJ-1 review, Mov Disord, doi:10.1002/mds.27352; GeneReviews Parkinson Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1223/) |
| PARK7 | AR | 15-45 | high | Kasten et al. 2018 MDSGene Parkin/PINK1/DJ-1 review, Mov Disord, doi:10.1002/mds.27352; GeneReviews Parkinson Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1223/) |
| CHCHD2 | AD | 35-70 | medium | OMIM #616710 PARK22 (https://omim.org/entry/616710); mean onset ~56 y in Japanese families; GeneReviews Parkinson Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1223/); Replication uneven outside East Asia |
| RAB32 | AD | 30-80 | high | GeneReviews RAB32-Related PD (https://www.ncbi.nlm.nih.gov/books/NBK620406/); median onset 54 y (31-82); Gustavsson et al. 2024 Lancet Neurol, doi:10.1016/S1474-4422(24)00121-2 (p.Ser71Arg only) |
| DAGLB | AR | 15-40 | medium | Liu et al. 2022 Nat Commun, doi:10.1038/s41467-022-31168-9 (6 Chinese patients, onset <=40 y) |

## Early-onset / atypical parkinsonism (PARK genes with distinctive signs)

| Gene | MOI | Onset | Confidence | Sources |
| --- | --- | --- | --- | --- |
| ATP13A2 | AR | 10-30 | high | Wittke et al. 2021 MDSGene atypical parkinsonism review, Mov Disord, doi:10.1002/mds.28517 (PMC9070562); GeneReviews Parkinson Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1223/) |
| DNAJC6 | AR | 5-45 | high | Wittke et al. 2021 MDSGene atypical parkinsonism review, Mov Disord, doi:10.1002/mds.28517 (PMC9070562); median onset 11 y (IQR 10-29) |
| SYNJ1 | AR | 10-40 | high | Wittke et al. 2021 MDSGene atypical parkinsonism review, Mov Disord, doi:10.1002/mds.28517 (PMC9070562); median onset 22 y |
| FBXO7 | AR | 10-30 | high | Wittke et al. 2021 MDSGene atypical parkinsonism review, Mov Disord, doi:10.1002/mds.28517 (PMC9070562); OMIM #260300 PARK15 |
| VPS13C | AR | 20-45 | high | Wittke et al. 2021 MDSGene atypical parkinsonism review, Mov Disord, doi:10.1002/mds.28517 (PMC9070562); Lesage et al. 2016 Am J Hum Genet |
| PLA2G6 | AR | 0-40 | high | GeneReviews PLA2G6-Associated Neurodegeneration (https://www.ncbi.nlm.nih.gov/books/NBK1675/); GeneReviews NBIA Overview (https://www.ncbi.nlm.nih.gov/books/NBK121988/) |
| DCTN1 | AD | 35-70 | high | Wittke et al. 2021 MDSGene atypical parkinsonism review, Mov Disord, doi:10.1002/mds.28517 (PMC9070562); Perry syndrome median onset 49 y; Lange et al. 2022 MDS Task Force nomenclature update (Mov Disord) |
| PSMF1 | AR | 0-50 | medium | Magrinelli et al. 2024, PMC11326324 (22 patients; missense -> EOPD, splice/LoF -> childhood/perinatal) |
| PTPA | AR | 10-30 | medium | Fevga et al. 2023 Brain, doi:10.1093/brain/awac326 (PMC10115167; 2 families, onset 11-16 y) |
| PTRHD1 | AR | 15-40 | medium | Khodadadi et al. 2017 Mov Disord; OMIM 620747 |
| RAB39B | XL | 10-50 | high | GeneReviews Parkinson Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1223/); OMIM: https://omim.org/search?search=RAB39B |
| ATP6AP2 | XL | 14-58 | medium | OMIM #300911 XPDS (https://omim.org/entry/300911); onset 14-58 y; Korvatska et al. 2013 Hum Mol Genet, doi:10.1093/hmg/ddt180 |
| WDR45 | XL | 0-40 | high | GeneReviews BPAN (https://www.ncbi.nlm.nih.gov/books/NBK424403/); GeneReviews NBIA Overview (https://www.ncbi.nlm.nih.gov/books/NBK121988/) |

## Risk genes (moi: risk). Onset is not scored unless the literature gives an anchor.

| Gene | MOI | Onset | Confidence | Sources |
| --- | --- | --- | --- | --- |
| GBA1 | risk | 40-75 | high | GeneReviews Parkinson Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1223/); Sidransky et al. 2009 NEJM, doi:10.1056/NEJMoa0901281; GeneReviews Gaucher Disease (https://www.ncbi.nlm.nih.gov/books/NBK1269/) |
| PSAP | risk | not scored | medium | OMIM #619491 PARK24 (https://omim.org/entry/619491); GeneReviews Parkinson Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1223/); Phenotype data sparse; no onset anchor |
| SMPD1 | risk | not scored | medium | Gan-Or et al. 2013 Neurology (p.L302P Ashkenazi founder); OMIM: https://omim.org/search?search=SMPD1 |
| TMEM175 | risk | not scored | medium | Jinn et al. 2019 (PMC6859430); GWAS risk p.M393T; Not Mendelian; no distinctive phenotype |
| HEXA | AR | 0-40 | high | GeneReviews HEXA Disorders (https://www.ncbi.nlm.nih.gov/books/NBK1218/); Late-onset Tay-Sachs presents teens-30s; AJ founder alleles studied as PD risk |

## Low-evidence / unreplicated PARK candidates. Kept minimal on purpose: no onset scoring, few or no keywords.

| Gene | MOI | Onset | Confidence | Sources |
| --- | --- | --- | --- | --- |
| DNAJC13 | unknown | not scored | low | GeneReviews Parkinson Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1223/); Farhan et al. 2020 Neurobiol Aging (not supported) |
| EIF4G1 | unknown | not scored | low | GeneReviews Parkinson Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1223/); Large case-control cohorts do not support PARK18 |
| GIGYF2 | unknown | not scored | low | GeneReviews Parkinson Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1223/); PARK11 not supported as Mendelian PD |
| TMEM230 | unknown | not scored | low | GeneReviews Parkinson Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1223/); Deng et al. 2016 contested; unreplicated |
| RIC3 | unknown | not scored | low | Single-family report; unreplicated (Puschmann 2017 review) |
| UQCRC1 | unknown | not scored | low | Lin et al. 2020 Brain (Taiwanese AD parkinsonism + polyneuropathy); not replicated in Europeans |
| WASL | unknown | not scored | low | medRxiv 2020 single EOPD pedigree; unreplicated |
| PODXL | unknown | not scored | low | Single AR juvenile PD pedigree; awaiting confirmation |
| SOX6 | unknown | not scored | low | No established Mendelian PD phenotype (mechanistic literature only) |
| NR4A2 | unknown | not scored | low | OMIM: https://omim.org/search?search=NR4A2; Scattered DRD/parkinsonism-with-ID reports; not in GeneReviews PD tables |
| ANG | unknown | not scored | low | OMIM: https://omim.org/search?search=ANG; ALS/PD association contested |
| SCP2 | AR | not scored | low | OMIM: https://omim.org/search?search=SCP2; Single adult case (Ferdinandusse et al. 2006) |
| PDE8B | AD | not scored | low | OMIM: https://omim.org/search?search=PDE8B; Appenzeller et al. 2010 single ADSD family |
| ALPL | AD | not scored | low | GeneReviews Hypophosphatasia (https://www.ncbi.nlm.nih.gov/books/NBK1150/); No parkinsonism-specific phenotype data; adult HPP skeletal/dental keywords only |

## Dystonia, dystonia-parkinsonism and monoamine-neurotransmitter disorders

| Gene | MOI | Onset | Confidence | Sources |
| --- | --- | --- | --- | --- |
| ANO3 | AD | 5-45 | high | GeneReviews Hereditary Dystonia Overview (https://www.ncbi.nlm.nih.gov/books/NBK1155/); OMIM: https://omim.org/search?search=ANO3 |
| AOPEP | AR | 5-40 | high | Zech et al. 2022 Mov Disord (AOPEP biallelic dystonia); OMIM: https://omim.org/search?search=AOPEP |
| ATP1A3 | AD | 0-45 | high | GeneReviews ATP1A3-Related Neurologic Disorders (https://www.ncbi.nlm.nih.gov/books/NBK1115/) |
| KMT2B | AD | 0-25 | high | GeneReviews Hereditary Dystonia Overview (https://www.ncbi.nlm.nih.gov/books/NBK1155/); OMIM: https://omim.org/search?search=KMT2B |
| PRKRA | AR | 2-30 | high | GeneReviews Hereditary Dystonia Overview (https://www.ncbi.nlm.nih.gov/books/NBK1155/); OMIM: https://omim.org/search?search=PRKRA |
| TUBB4A | AD | 0-40 | high | GeneReviews Hereditary Dystonia Overview (https://www.ncbi.nlm.nih.gov/books/NBK1155/); OMIM: https://omim.org/search?search=TUBB4A |
| GCH1 | AD | 0-50 | high | GeneReviews GTP Cyclohydrolase 1-Deficient DRD (https://www.ncbi.nlm.nih.gov/books/NBK1508/); Adult-onset parkinsonism and PD risk reported (Mencacci et al. 2014 Brain) |
| TH | AR | 0-12 | high | GeneReviews Tyrosine Hydroxylase Deficiency (https://www.ncbi.nlm.nih.gov/books/NBK1437/) |
| SPR | AR | 0-10 | high | GeneReviews Sepiapterin Reductase Deficiency (https://www.ncbi.nlm.nih.gov/books/NBK304122/) |
| PTS | AR | 0-5 | high | OMIM: https://omim.org/search?search=PTS; BH4 deficiency reviews (Opladen et al. 2020 Orphanet J Rare Dis) |
| QDPR | AR | 0-5 | high | OMIM: https://omim.org/search?search=QDPR; Opladen et al. 2020 Orphanet J Rare Dis (BH4 guideline) |
| DDC | AR | 0-2 | high | OMIM: https://omim.org/search?search=DDC; Wassenberg et al. 2017 Orphanet J Rare Dis (AADC consensus) |
| DNAJC12 | AR | 0-40 | high | Anikster et al. 2017 Am J Hum Genet; OMIM: https://omim.org/search?search=DNAJC12 |
| SLC6A3 | AR | 0-30 | high | OMIM: https://omim.org/search?search=SLC6A3; Kurian et al. 2011 Lancet Neurol (DTDS incl. atypical later onset) |
| SLC18A2 | AR | 0-5 | high | OMIM: https://omim.org/search?search=SLC18A2; Rilstone et al. 2013 NEJM |
| SLC30A10 | AR | 0-40 | high | OMIM: https://omim.org/search?search=SLC30A10; Tuschl et al. 2012 Am J Hum Genet; adult parkinsonism reported |
| SLC39A14 | AR | 0-15 | high | OMIM: https://omim.org/search?search=SLC39A14; Tuschl et al. 2016 Nat Commun |
| EIF2AK2 | AD | 0-30 | medium | OMIM: https://omim.org/search?search=EIF2AK2; Mao et al. 2020 Am J Hum Genet; Kuipers et al. 2021 (dystonia) |

## NBIA and metal-accumulation disorders

| Gene | MOI | Onset | Confidence | Sources |
| --- | --- | --- | --- | --- |
| PANK2 | AR | 0-40 | high | GeneReviews PKAN (https://www.ncbi.nlm.nih.gov/books/NBK1490/); GeneReviews NBIA Overview (https://www.ncbi.nlm.nih.gov/books/NBK121988/) |
| C19ORF12 | AR | 3-40 | high | GeneReviews MPAN (https://www.ncbi.nlm.nih.gov/books/NBK185329/); GeneReviews NBIA Overview (https://www.ncbi.nlm.nih.gov/books/NBK121988/) |
| COASY | AR | 0-20 | high | GeneReviews NBIA Overview (https://www.ncbi.nlm.nih.gov/books/NBK121988/); OMIM: https://omim.org/search?search=COASY |
| FA2H | AR | 0-15 | high | GeneReviews NBIA Overview (https://www.ncbi.nlm.nih.gov/books/NBK121988/); OMIM: https://omim.org/search?search=FA2H |
| FTL | AD | 15-60 | high | GeneReviews Neuroferritinopathy (https://www.ncbi.nlm.nih.gov/books/NBK1141/) |
| CP | AR | 25-60 | high | GeneReviews Aceruloplasminemia (https://www.ncbi.nlm.nih.gov/books/NBK1493/) |
| ATP7B | AR | 5-45 | high | GeneReviews Wilson Disease (https://www.ncbi.nlm.nih.gov/books/NBK1512/) |

## Primary familial brain calcification (PFBC)

| Gene | MOI | Onset | Confidence | Sources |
| --- | --- | --- | --- | --- |
| SLC20A2 | AD | 20-70 | high | Balck et al. 2021 MDSGene PFBC review, Mov Disord, doi:10.1002/mds.28753; median onset 47 y; GeneReviews Primary Familial Brain Calcification (https://www.ncbi.nlm.nih.gov/books/NBK1421/) |
| PDGFB | AD | 5-55 | high | Balck et al. 2021 MDSGene PFBC review, Mov Disord, doi:10.1002/mds.28753; median onset 30 y; GeneReviews Primary Familial Brain Calcification (https://www.ncbi.nlm.nih.gov/books/NBK1421/) |
| PDGFRB | AD | 20-70 | high | Balck et al. 2021 MDSGene PFBC review, Mov Disord, doi:10.1002/mds.28753; lowest penetrance; GeneReviews Primary Familial Brain Calcification (https://www.ncbi.nlm.nih.gov/books/NBK1421/) |
| XPR1 | AD | 20-70 | high | Balck et al. 2021 MDSGene PFBC review, Mov Disord, doi:10.1002/mds.28753; median onset 44 y; GeneReviews Primary Familial Brain Calcification (https://www.ncbi.nlm.nih.gov/books/NBK1421/) |
| MYORG | AR | 20-70 | high | Balck et al. 2021 MDSGene PFBC review, Mov Disord, doi:10.1002/mds.28753; median onset 46 y; GeneReviews Primary Familial Brain Calcification (https://www.ncbi.nlm.nih.gov/books/NBK1421/) |
| JAM2 | AR | 5-45 | high | Balck et al. 2021 MDSGene PFBC review, Mov Disord, doi:10.1002/mds.28753; few onset data (20-38 y); PMC10342144 |

## Differential: Alzheimer, FTD, prion, ALS spectrum

| Gene | MOI | Onset | Confidence | Sources |
| --- | --- | --- | --- | --- |
| APP | AD | 30-65 | high | GeneReviews Alzheimer Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1161/); OMIM: https://omim.org/search?search=APP |
| PSEN1 | AD | 25-65 | high | GeneReviews Alzheimer Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1161/); OMIM: https://omim.org/search?search=PSEN1 |
| PSEN2 | AD | 40-80 | high | GeneReviews Alzheimer Disease Overview (https://www.ncbi.nlm.nih.gov/books/NBK1161/); OMIM: https://omim.org/search?search=PSEN2 |
| MAPT | AD | 25-65 | high | GeneReviews MAPT-Related Disorders (https://www.ncbi.nlm.nih.gov/books/NBK1505/) |
| GRN | AD | 35-80 | high | GeneReviews GRN Frontotemporal Dementia (https://www.ncbi.nlm.nih.gov/books/NBK1371/) |
| PRNP | AD | 20-80 | high | GeneReviews Genetic Prion Disease (https://www.ncbi.nlm.nih.gov/books/NBK1229/) |
| CHCHD10 | AD | 20-70 | medium | OMIM: https://omim.org/search?search=CHCHD10; Bannwarth et al. 2014 Brain |
| FUS | AD | 15-60 | high | GeneReviews Amyotrophic Lateral Sclerosis Overview (https://www.ncbi.nlm.nih.gov/books/NBK1450/); OMIM: https://omim.org/search?search=FUS |
| TBK1 | AD | 40-80 | high | GeneReviews Amyotrophic Lateral Sclerosis Overview (https://www.ncbi.nlm.nih.gov/books/NBK1450/); OMIM: https://omim.org/search?search=TBK1 |
| SQSTM1 | AD | 40-80 | high | GeneReviews Amyotrophic Lateral Sclerosis Overview (https://www.ncbi.nlm.nih.gov/books/NBK1450/); OMIM: https://omim.org/search?search=SQSTM1 |
| VCP | AD | 30-70 | high | GeneReviews VCP Multisystem Proteinopathy (https://www.ncbi.nlm.nih.gov/books/NBK1476/) |

## Differential: ataxia, HSP and neuropathy

| Gene | MOI | Onset | Confidence | Sources |
| --- | --- | --- | --- | --- |
| AFG3L2 | AD | 10-60 | high | GeneReviews SCA28 (https://www.ncbi.nlm.nih.gov/books/NBK54582/) |
| PRKCG | AD | 10-60 | high | GeneReviews SCA14 (https://www.ncbi.nlm.nih.gov/books/NBK1399/) |
| STUB1 | AR | 10-60 | medium | OMIM: https://omim.org/search?search=STUB1; SCA48 (AD, later onset with cognitive/psychiatric signs) and SCAR16 (AR) |
| SPG11 | AR | 5-35 | high | GeneReviews SPG11 (https://www.ncbi.nlm.nih.gov/books/NBK1210/); GeneReviews Hereditary Spastic Paraplegia Overview (https://www.ncbi.nlm.nih.gov/books/NBK1509/) |
| SPG7 | AR | 20-65 | high | GeneReviews SPG7 (https://www.ncbi.nlm.nih.gov/books/NBK1107/); GeneReviews Hereditary Spastic Paraplegia Overview (https://www.ncbi.nlm.nih.gov/books/NBK1509/) |
| ZFYVE26 | AR | 5-30 | high | GeneReviews Hereditary Spastic Paraplegia Overview (https://www.ncbi.nlm.nih.gov/books/NBK1509/); OMIM: https://omim.org/search?search=ZFYVE26 |
| KIF5A | AD | 5-60 | high | GeneReviews Hereditary Spastic Paraplegia Overview (https://www.ncbi.nlm.nih.gov/books/NBK1509/); GeneReviews Amyotrophic Lateral Sclerosis Overview (https://www.ncbi.nlm.nih.gov/books/NBK1450/) |
| PNPLA6 | AR | 0-40 | high | OMIM: https://omim.org/search?search=PNPLA6; Synofzik et al. 2014 Brain (PNPLA6 spectrum) |
| FIG4 | AR | 0-40 | medium | OMIM: https://omim.org/search?search=FIG4; GeneReviews CMT overview (CMT4J) |
| SLC25A46 | AR | 0-30 | medium | OMIM: https://omim.org/search?search=SLC25A46; Abrams et al. 2015 Nat Genet |
| FLVCR1 | AR | 0-30 | medium | OMIM: https://omim.org/search?search=FLVCR1 |

## Differential: mitochondrial disease

| Gene | MOI | Onset | Confidence | Sources |
| --- | --- | --- | --- | --- |
| POLG | AR | 10-60 | medium | GeneReviews POLG-Related Disorders (https://www.ncbi.nlm.nih.gov/books/NBK26471/) |
| TWNK | AD | 20-70 | medium | OMIM: https://omim.org/search?search=TWNK; adPEO with parkinsonism reported (Baloh et al. 2007) |
| SLC25A4 | AD | 20-60 | medium | OMIM: https://omim.org/search?search=SLC25A4 |
| OPA1 | AD | 0-50 | high | GeneReviews Optic Atrophy Type 1 (https://www.ncbi.nlm.nih.gov/books/NBK1248/); Carelli et al. 2015 Brain (OPA1 parkinsonism/dementia) |
| OPA3 | AD | 0-40 | medium | OMIM: https://omim.org/search?search=OPA3 |
| SERAC1 | AR | 0-15 | high | OMIM: https://omim.org/search?search=SERAC1; GeneReviews MEGDEL syndrome |
| WARS2 | AR | 0-40 | medium | OMIM: https://omim.org/search?search=WARS2; Burke et al. 2018 Mov Disord (infantile parkinsonism) |

## Differential: lysosomal, metabolic and leukodystrophy

| Gene | MOI | Onset | Confidence | Sources |
| --- | --- | --- | --- | --- |
| ARSA | AR | 0-60 | high | GeneReviews Arylsulfatase A Deficiency (https://www.ncbi.nlm.nih.gov/books/NBK1130/) |
| CLN3 | AR | 3-15 | high | GeneReviews Neuronal Ceroid-Lipofuscinoses (https://www.ncbi.nlm.nih.gov/books/NBK1428/) |
| DNAJC5 | AD | 20-50 | high | GeneReviews Neuronal Ceroid-Lipofuscinoses (https://www.ncbi.nlm.nih.gov/books/NBK1428/); OMIM: https://omim.org/search?search=DNAJC5 |
| TPP1 | AR | 1-10 | high | GeneReviews Neuronal Ceroid-Lipofuscinoses (https://www.ncbi.nlm.nih.gov/books/NBK1428/) |
| GLB1 | AR | 0-30 | high | GeneReviews GLB1-Related Disorders (https://www.ncbi.nlm.nih.gov/books/NBK164500/) |
| NPC1 | AR | 0-60 | high | GeneReviews Niemann-Pick Disease Type C (https://www.ncbi.nlm.nih.gov/books/NBK1296/) |
| NPC2 | AR | 0-60 | high | GeneReviews Niemann-Pick Disease Type C (https://www.ncbi.nlm.nih.gov/books/NBK1296/) |
| LYST | AR | 15-50 | high | GeneReviews Chediak-Higashi Syndrome (https://www.ncbi.nlm.nih.gov/books/NBK5188/) |
| CYP27A1 | AR | 10-50 | high | OMIM: https://omim.org/search?search=CYP27A1; GeneReviews Cerebrotendinous Xanthomatosis |
| CSF1R | AD | 20-70 | high | GeneReviews CSF1R-Related Leukoencephalopathy (https://www.ncbi.nlm.nih.gov/books/NBK100239/) |
| POLR3A | AR | 0-40 | high | GeneReviews POLR3-Related Leukodystrophy (https://www.ncbi.nlm.nih.gov/books/NBK99167/) |
| SLC19A3 | AR | 0-30 | high | GeneReviews Biotin-Thiamine-Responsive Basal Ganglia Disease (https://www.ncbi.nlm.nih.gov/books/NBK169615/) |
| VPS13A | AR | 15-50 | high | GeneReviews Chorea-Acanthocytosis (https://www.ncbi.nlm.nih.gov/books/NBK1387/) |
| VAC14 | AR | 0-15 | medium | OMIM: https://omim.org/search?search=VAC14; Lenk et al. 2016 Am J Hum Genet |
| PGK1 | XL | 0-30 | medium | OMIM: https://omim.org/search?search=PGK1; Juvenile parkinsonism case reports (Sotiriou et al. 2010) |
| ADAR | AD | 0-15 | medium | GeneReviews Aicardi-Goutieres Syndrome (https://www.ncbi.nlm.nih.gov/books/NBK1475/); OMIM: https://omim.org/search?search=ADAR |

## Differential: neurodevelopmental and epilepsy genes. onset_range spans the neurological presentation and documented later parkinsonism.

| Gene | MOI | Onset | Confidence | Sources |
| --- | --- | --- | --- | --- |
| DHDDS | AD | 0-35 | medium | OMIM: https://omim.org/search?search=DHDDS; Galosi et al. 2022 (DHDDS/NUS1 parkinsonism) |
| FOXG1 | AD | 0-5 | high | OMIM: https://omim.org/search?search=FOXG1; GeneReviews FOXG1 Syndrome |
| FRRS1L | AR | 0-5 | high | OMIM: https://omim.org/search?search=FRRS1L |
| HIVEP2 | AD | 0-10 | medium | OMIM: https://omim.org/search?search=HIVEP2 |
| MECP2 | XL | 0-40 | high | GeneReviews MECP2 Disorders (https://www.ncbi.nlm.nih.gov/books/NBK1497/) |
| NAA15 | AD | 0-10 | medium | OMIM: https://omim.org/search?search=NAA15 |
| NUS1 | AD | 0-45 | medium | OMIM: https://omim.org/search?search=NUS1; Guo et al. 2018 (NUS1 early-onset PD); Galosi et al. 2022 |
| PPP2R5D | AD | 0-45 | high | OMIM: https://omim.org/search?search=PPP2R5D; Kim et al. 2020 (E200K early-onset parkinsonism) |
| SCN1A | AD | 0-45 | high | GeneReviews SCN1A Seizure Disorders (https://www.ncbi.nlm.nih.gov/books/NBK1318/); Fasano et al. 2014 (parkinsonism in adult Dravet) |
| SLC9A6 | XL | 0-40 | high | OMIM: https://omim.org/search?search=SLC9A6; GeneReviews Christianson Syndrome |
| STXBP1 | AD | 0-45 | high | OMIM: https://omim.org/search?search=STXBP1; Adult STXBP1 parkinsonism reports (Stamberger 2016 review) |
| TBC1D24 | AR | 0-20 | medium | OMIM: https://omim.org/search?search=TBC1D24 |
| UBTF | AD | 0-15 | medium | OMIM: https://omim.org/search?search=UBTF; Edvardson et al. 2017 Am J Hum Genet (E210K) |
| EPM2A | AR | 8-20 | high | GeneReviews Progressive Myoclonus Epilepsy, Lafora Type (https://www.ncbi.nlm.nih.gov/books/NBK1389/) |
| NHLRC1 | AR | 8-25 | high | GeneReviews Progressive Myoclonus Epilepsy, Lafora Type (https://www.ncbi.nlm.nih.gov/books/NBK1389/) |
| DCAF17 | AR | 10-30 | high | OMIM: https://omim.org/search?search=DCAF17; GeneReviews Woodhouse-Sakati Syndrome |
