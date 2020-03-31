For each rat, run all of seurat pipeline for

-   dimensionality reduction,
-   clustering,

``` r
suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
suppressMessages(require(gridExtra))
suppressMessages(require(ggplot2))
suppressMessages(require(pheatmap))
```

Function to run all Seurat steps.

``` r
sel.pc <- 1:30
npcs <- 30

run_seurat <- function(sample, dataA){

  # create seurat object for that sample
  data <- subset(x = dataA, cells = colnames(dataA)[dataA@meta.data$Rat == sample])

  data <- NormalizeData(object = data, verbose = F)
  # find variable genes, default is top 2K variable genes.
  data <- FindVariableFeatures(object = data, verbose = F)
  data <- ScaleData(object = data, vars.to.regress = c("nGene","percent.mito"), verbose = F)
  data <- RunPCA(object = data, features = VariableFeatures(object = data), verbose=F)
  data <- JackStraw(data, verbose=F)
  data <- ScoreJackStraw(data, dims = 1:20)
  sel.pc <- which(JS(object = data[['pca']], slot = 'overall')[,2] < 1e-3)
  
  # cluster
  set.seed(1)
  data <- FindNeighbors(object = data, dims = sel.pc, verbose = F)
  data <- FindClusters(object = data, resolution = 0.6, verbose = F)

  # run tSNE
  set.seed(1)
  data <- RunTSNE(object = data, dims = sel.pc)

  # runUMAP
  set.seed(1)
  data <- RunUMAP(object = data, dims = sel.pc)
  
  return(data)
}
```

### Plot overview of the data

First plot PCA for first PCs, and the top loading genes for the PCs.

Then tSNE and UMAP with clustering + detected genes

``` r
plot_seurat <- function(data){
  
  # first plot some pca plots
  p1 <- DimPlot(data, reduction = "pca")
  p2 <- DimPlot(data, reduction = "pca", dims = 3:4)
  p3 <- DimPlot(data, reduction = "pca", dims = 5:6)
  p4 <- DimPlot(data, reduction = "pca", dims = 7:8)
  grid.arrange(p1,p2,p3,p4, nrow=2)
  # plot gene loadings to first 12 pcs.
  print(DimHeatmap(data, dims=1:6, nfeatures = 20))
  print(DimHeatmap(data, dims=7:12, nfeatures = 20))

  
  # then tSNE + UMAP
  small.leg <- theme(legend.text = element_text(size=6))
  p1 <- DimPlot(data, label = T, reduction = "tsne") + NoAxes() + small.leg
  p2 <- FeaturePlot(data, features = "nFeature_RNA", reduction = "tsne") + NoAxes() + small.leg
  p3 <- FeaturePlot(data, features = "percent.mito", reduction = "tsne") + NoAxes() + small.leg
  grid.arrange(p1,p2,p3, nrow=2)  
  
  # plot split by sample
  print(DimPlot(data, split.by = "Sample", reduction = "tsne"))
  
  
  p1 <- DimPlot(data, label = T, reduction = "umap") + NoAxes() + small.leg
  p2 <- FeaturePlot(data, features = "nFeature_RNA", reduction = "umap") + NoAxes() + small.leg
  p3 <- FeaturePlot(data, features = "percent.mito", reduction = "umap") + NoAxes() + small.leg
  grid.arrange(p1,p2,p3, nrow=2)  

  # plot split by sample
  print(DimPlot(data, split.by = "Sample", reduction = "umap"))  
  
  # plot some marker genes.
  plot.genes <- c("TH","SLC6A3","SNAP25","GAP43", "NEFL","DLK1","OLIG1","GFAP","AQP4","COL1A1","DCN","FBLN1")
  p1 <- FeaturePlot(data, features = plot.genes, reduction = "tsne", pt.size = 0.2, combine = F, cols = c("yellow","red","black"))
  p1 <- lapply(p1, function(x) x + NoAxes() + NoLegend())
  grid.arrange(grobs=p1, ncol=4)

}
```

``` r
force <- FALSE

# use the clustered seurat object with all cells, 
# will give coloring according to cluster in that analysis
load("../data/processed2/filtered_seurat_object.Rdata")
samples <- unique(dataA@meta.data$Rat)

for (sample in samples) {
  print(sample)
  plotfile = sprintf("../data/processed2/per_sample/Seurat_rat_%s.pdf", sample)
  save.file = sprintf("../data/processed2/per_sample/seurat_object_%s.Rdata",sample)
  
  if (file.exists(save.file) & !force){
    load(save.file)
  }else{
    ss <- run_seurat(sample, dataA)
    save(ss, file=save.file)
  }
  
  pdf(plotfile)
  plot_seurat(ss)
  dev.off()
}
```

    ## [1] "rat45"

    ## NULL

    ## NULL

    ## [1] "rat11"

    ## NULL

    ## NULL

    ## [1] "rat39"

    ## NULL

    ## NULL

``` r
sessionInfo()
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-apple-darwin13.4.0 (64-bit)
    ## Running under: macOS  10.15.2
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /Users/asbj/miniconda3/envs/seurat3/lib/R/lib/libRblas.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] pheatmap_1.0.12 ggplot2_3.2.1   gridExtra_2.3   Matrix_1.2-17  
    ## [5] Seurat_3.0.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tsne_0.1-3          nlme_3.1-141        bitops_1.0-6       
    ##  [4] RColorBrewer_1.1-2  httr_1.4.1          rprojroot_1.3-2    
    ##  [7] sctransform_0.2.0   tools_3.5.1         backports_1.1.5    
    ## [10] R6_2.4.0            irlba_2.3.3         KernSmooth_2.23-15 
    ## [13] lazyeval_0.2.2      colorspace_1.4-1    withr_2.1.2        
    ## [16] npsurv_0.4-0        tidyselect_0.2.5    compiler_3.5.1     
    ## [19] plotly_4.9.0        labeling_0.3        caTools_1.17.1.2   
    ## [22] scales_1.0.0        lmtest_0.9-37       ggridges_0.5.1     
    ## [25] pbapply_1.4-2       stringr_1.4.0       digest_0.6.22      
    ## [28] rmarkdown_1.10      R.utils_2.9.0       pkgconfig_2.0.3    
    ## [31] htmltools_0.4.0     bibtex_0.4.2        htmlwidgets_1.5.1  
    ## [34] rlang_0.4.1         zoo_1.8-6           jsonlite_1.6       
    ## [37] ica_1.0-2           gtools_3.8.1        dplyr_0.8.3        
    ## [40] R.oo_1.22.0         magrittr_1.5        Rcpp_1.0.2         
    ## [43] munsell_0.5.0       ape_5.3             reticulate_1.13    
    ## [46] lifecycle_0.1.0     R.methodsS3_1.7.1   stringi_1.4.3      
    ## [49] yaml_2.2.0          gbRd_0.4-11         MASS_7.3-51.4      
    ## [52] gplots_3.0.1.1      Rtsne_0.15          plyr_1.8.4         
    ## [55] grid_3.5.1          parallel_3.5.1      gdata_2.18.0       
    ## [58] listenv_0.7.0       ggrepel_0.8.1       crayon_1.3.4       
    ## [61] lattice_0.20-38     cowplot_1.0.0       splines_3.5.1      
    ## [64] SDMTools_1.1-221.1  zeallot_0.1.0       knitr_1.20         
    ## [67] pillar_1.4.2        igraph_1.2.4.1      future.apply_1.3.0 
    ## [70] reshape2_1.4.3      codetools_0.2-16    glue_1.3.1         
    ## [73] evaluate_0.14       lsei_1.2-0          metap_1.1          
    ## [76] data.table_1.12.6   vctrs_0.2.0         png_0.1-7          
    ## [79] Rdpack_0.11-0       gtable_0.3.0        RANN_2.6.1         
    ## [82] purrr_0.3.3         tidyr_1.0.0         future_1.14.0      
    ## [85] assertthat_0.2.1    rsvd_1.0.2          survival_2.44-1.1  
    ## [88] viridisLite_0.3.0   tibble_2.1.3        cluster_2.1.0      
    ## [91] globals_0.12.4      fitdistrplus_1.0-14 ROCR_1.0-7
