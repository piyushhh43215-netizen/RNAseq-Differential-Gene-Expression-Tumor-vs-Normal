# Data and reproducibility notes

## Source data

The analysis uses the GEO series **GSE81089** and expects a gene-level raw read-count matrix. The raw matrix is intentionally excluded from this GitHub package; retrieve it from GEO/the original source and keep the source attribution with your project.

## Input format expected

The R script expects:
- one gene identifier column named `Ensembl_gene_id`;
- one column per sample containing integer raw read counts;
- sample IDs compatible with the Tumor/Normal and paired-patient naming assumptions in the script.

## Before publishing final scientific conclusions

- Run the R script from start to finish in a clean R session.
- Verify the sample labels and patient pairing against the GEO metadata.
- Confirm that the design formula is appropriate for the selected samples.
- Check package versions and generated `sessionInfo.txt`.
- Interpret adjusted p-values and log2 fold changes alongside biological context.

The CSV files in `results/` are prior output files retained from the earlier workflow. The reconstructed script writes new files to a separate output folder so the older outputs are not overwritten.
