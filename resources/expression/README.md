# Brain expression table

`pdvar isoform` reads a TSV with these columns:

| Column | Meaning |
|--------|---------|
| chrom, start, end | Interval. Empty for a gene-level GTEx row. |
| gene | HGNC symbol |
| gene_id | Ensembl gene id when known |
| brain_pext | gnomAD proportion expressed across transcripts, if you exported one |
| gtex_brain_median_tpm | Median of substantia nigra, putamen, caudate, and frontal cortex |
| expression_resolution | `gene` or `interval` |
| source | Where the row came from |

The shipped config leaves `expression.resource`, `min_brain_pext`, `min_gtex_tpm`, and `unexpressed_weight` null. Weights stay 1 until a cutoff and `unexpressed_weight` are both set.

GTEx v8 gene-median TPM (GRCh38, gene-level, coarse relative to an exon):

`https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`

gnomAD pext is not a flat HTTPS file. v2 base-level data is a requester-pays Hail table (`gs://gnomad-public-requester-pays/papers/2019-tx-annotation/gnomad_browser/all.baselevel.021620.ht`, GRCh37 / GTEx v7). v4 pext uses GTEx v10 Hail tables. Export `chrom, start, end, gene, brain_pext` yourself. See https://gnomad.broadinstitute.org/help/pext.

```bash
python scripts/fetch_expression_resources.py print-urls
python scripts/fetch_expression_resources.py gtex --input resources/expression/downloads/gtex.gct.gz --output resources/expression/brain_expression.tsv
```

Put downloads under `resources/expression/downloads/`. That directory is gitignored.
