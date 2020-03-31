# DA_grafting_analysis
Analysis of single cell RNA-seq data from grafting of human cells in rat brains. 



All code required to reproduce the main analysis figures of the paper “Single Cell Gene Expression Analysis Reveals Human Stem Cell-Derived Graft Composition in a Cell Therapy Model of Parkinson’s Disease" by Tiklova et al. 

The analysis has been performed in two steps:

* All Smart-seq2 data analysis was performed with Seurat v2.3.4. All scripts for that analysis is in folder [bin][bin/].
* All 10x and integrated analysis was performed with Seurat v3.0.0. All scripts for that analysis is in folder [bin_10x][bin_10x/].

SS2 data is deposited at GEO GSE118412, expression matrices are provided as gz files in folder [data][data/]

10x data is deposited at GEO GSE132758.

## Code for SS2 analysis

The scripts used to produce all manuscript figures are:

* Seurat analysis for the two sets of cells ( before grafting and after grafting ) and cell cycle scoring [Seurat_analysis](bin/Seurat_analysis.md).
* DE detection with SAMseq [run_samseq](bin/run_samseq.md).
* Comparison to the Mouse Brain Atlas clusters [compare_to_linnarsson_atlas](bin/compare_to_linnarsson_atlas.md).
* Plotting gene expression in violin plots, heatmaps and onto tsne [plot_genes](plot_genes.md)

## Code for 10x analysis

The scripts used to produce all figures are:

* Elimination of rat contaminating cells [Check_rno_vs_hsa](bin_10x/Check_rno_vs_hsa.md)
* QC filtering [Summarize_QC_stats](bin_10x/Summarize_QC_stats.md)
* Seurat analysis with one rat at a time [Seurat_per_rat](bin_10x/Seurat_per_rat.md)
