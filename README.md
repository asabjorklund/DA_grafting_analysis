# DA_grafting_analysis
Analysis of single cell RNA-seq data from grafting of human cells in rat brains. 

All code required to reproduce the main analysis figures of the paper “Unbiased Gene Expression Analysis Reveals Human Stem Cell-Derived Graft Composition in a Cell Therapy Model of Parkinson’s Disease" by Tiklova et al. 

The scripts used to produce all figures are:

* Seurat analysis for the two sets of cells ( before grafting and after grafting ) and cell cycle scoring [Seurat_analysis](bin/Seurat_analysis.md).

* DE detection with SAMseq [run_samseq](bin/run_samseq.md).

* Comparison to the Mouse Brain Atlas clusters [compare_to_linnarsson_atlas](bin/compare_to_linnarsson_atlas.md).

* Plotting gene expression in violin plots, heatmaps and onto tsne [plot_genes](plot_genes.md)