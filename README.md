# RNA-seq Differential Gene Expression: Tumor vs Normal

A bioinformatics project analyzing human RNA-seq read counts from **GSE81089** to identify differentially expressed genes (DEGs) between tumor and normal samples, followed by functional enrichment analysis.

## Project overview

This project demonstrates an RNA-seq analysis workflow using R and Bioconductor:

1. Import a raw gene-level read-count matrix.
2. Prepare sample metadata and identify paired tumor/normal samples.
3. Filter genes with very low read counts.
4. Run differential expression analysis with **DESeq2**.
5. Annotate Ensembl gene IDs using **AnnotationDbi** and **org.Hs.eg.db**.
6. Perform Gene Ontology Biological Process and KEGG pathway enrichment using **clusterProfiler**.
7. Export tables and visualization files for interpretation.

## Repository structure

```text
RNAseq_Tumor_vs_Normal/
├── README.md
├── .gitignore
├── LICENSE
├── scripts/
│   └── RNAseq_Tumor_vs_Normal.R
├── results/
│   ├── Annotated_DEGs_Tumor_vs_Normal.csv
│   ├── DESeq2_significant_genes.csv
│   ├── DESeq2_significant_genes_annotated.csv
│   ├── GO_Biological_Process_Results.csv
│   └── KEGG_Pathway_Results.csv
└── docs/
    └── DATA_AND_REPRODUCIBILITY.md
```

## Dataset

- **Accession:** [GSE81089](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81089)
- **Input type:** gene-level raw read counts (featureCounts-style table)
- **Comparison:** Tumor vs Normal
- **Input file expected by the script:** `GSE81089_readcounts_featurecounts.tsv.gz`

The raw dataset is not included in this repository. Download it from the source and update `input_file` in the R script to point to its local path. Use raw counts for DESeq2; do not use FPKM values as DESeq2 input.

## Requirements

- R
- Bioconductor packages: `DESeq2`, `AnnotationDbi`, `org.Hs.eg.db`, `clusterProfiler`, `enrichplot`
- CRAN packages: `ggplot2`, `pheatmap`

The script checks for and installs missing packages.

## How to run

1. Download the raw count file from the GEO project/source and save it locally.
2. Open `scripts/RNAseq_Tumor_vs_Normal.R` in RStudio.
3. Edit the `input_file` path near the top of the script so it points to your downloaded file.
4. Review the sample metadata and paired-sample assumptions printed by the script.
5. Run the script from RStudio.
6. Results are written to the `results_rebuilt` folder configured in the script. Update `output_dir` if you prefer a different location.

## Main methods

- Differential expression: DESeq2
- Multiple-testing correction: Benjamini–Hochberg adjusted p-values
- Significance criteria in the script: `padj < 0.05` and `|log2FoldChange| > 1`
- Functional enrichment: GO Biological Process and KEGG

## Results included

The `results/` folder contains CSV exports from the earlier analysis workflow, including DEG tables and GO/KEGG enrichment tables. The R script writes a separate set of rebuilt outputs to avoid overwriting these files.

**Reproducibility note:** the script has been reconstructed from the analysis history and is provided as the project workflow. Run it end-to-end and review its metadata/design checks before treating the rebuilt outputs as final. Existing CSV exports are retained as prior outputs, not represented as proof that this exact script has already completed successfully in the current environment.

## Skills demonstrated

`R` · `RNA-seq` · `DESeq2` · `Differential Gene Expression` · `Bioconductor` · `Gene Annotation` · `GO Enrichment` · `KEGG Pathway Analysis` · `Data Visualization`

## Author

**Piyush Kumar**  
Bioinformatics & Computational Biology learner
