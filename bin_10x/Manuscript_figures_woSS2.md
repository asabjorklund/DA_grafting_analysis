Plot figures for manuscript. Use the integrated object from
Data\_integration.Rmd (file
seurat\_integrated\_rat\_scale\_res0.1.Rdata)

Plot UMAP with:

-   Clusters (4 main celltypes)
-   FACS/NotFACS
-   SS2/10x

Violin plots for the 4 clusters - same markers as in SS2 data only.

``` r
suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
suppressMessages(require(gridExtra))
suppressMessages(require(ggplot2))
suppressMessages(require(pheatmap))
```

``` r
load(file = "../data/processed2/integrated/seurat_integrated_rat_scale_res0.1.Rdata")

# reassing cluster names to celltypes
celltypes <- c("Astrocyte","FBLC","Neuron","Cycling")
names(celltypes) <- as.character(0:3)

ct <- celltypes[data.integrated@active.ident]
data.integrated <- AddMetaData(data.integrated, ct, col.name = "CelltypeCluster")
```

### plot UMAPs

``` r
# make one version with labels in plot
scaleC <- scale_colour_manual(values = c("blue1", "gray","green3", "red1"))

# Per cluster with legend
p1 <- DimPlot(data.integrated, group.by = "CelltypeCluster") + NoAxes() +scaleC

# one with FACS/notFACS
p2 <- DimPlot(data.integrated, group.by = "Type", order = c("NotFACS","FACS")) + NoAxes() + scale_colour_manual(values = c("cyan", "magenta"))

# one with Rat
p3 <- DimPlot(data.integrated, group.by = "Rat") + NoAxes()

pdf("../data/figures2/umap_integrated_woSS2.pdf")
grid.arrange(p1,p2,p3, ncol=2)

# all same plots without legends.
p1 <- DimPlot(data.integrated, group.by = "CelltypeCluster") + NoAxes() + 
  scaleC + NoLegend() + ggtitle("Clusters")
p2 <- DimPlot(data.integrated, group.by = "Type", order = c("NotFACS","FACS")) + NoAxes() + NoLegend() + 
  scale_colour_manual(values = c("cyan", "magenta")) + ggtitle("FACS type")
p3 <- DimPlot(data.integrated, group.by = "Rat") + NoAxes() + NoLegend() +  ggtitle("Rat")
grid.arrange(p1,p2,p3, ncol=2)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### gene violin plots as in SS2 figure.

``` r
genes <- list()
genes$Astrocyte <-  c("AQP4", "GFAP", "SLC1A3",  "GJA1", "EDNRB", "SLC4A4")
genes$Oligodentrocyte <- c("OLIG1", "OLIG2","NKX2-2",  "SOX10",  "PLP1",  "PMP2")
genes$PanNeuronal <- c("GAP43", "RBFOX3","NSG2","SNAP25", "NEFL", "SYN1")
genes$DopamineNeuron <-c("TH","NR4A2","SLC18A2", "DDC",  "RET", "GFRA1")
genes$FBLC <- c("PDGFRA",   "COL1A1",  "COL1A2",  "LUM",     "DCN",     "FBLN1")

temp <- data.frame(t(data.integrated@assays$RNA@data[unlist(genes),]), check.names = F)
temp$cluster <- data.integrated@meta.data$CelltypeCluster
short <- temp$cluster
short[short=="Cycling"] <- "CC"
short[short=="Neuron"] <- "N"
short[short=="Astrocyte"] <- "AC"
temp$cluster <- short

pl <- list()
noleg <- theme(legend.position="none")
scaleF <- scale_fill_manual(values = c("blue1", "gray","green3", "red1"))
for (gene in unlist(genes)){
  pl[[gene]] <- ggplot() + geom_violin(data=temp,aes_(x=as.name("cluster"),y=as.name(gene),fill=as.name("cluster"), colour = as.name("cluster") ), scale="width") + 
    theme_void() + noleg + scaleF + scaleC + theme(title =element_text(size=6, face='bold')) + ggtitle(gene)
}

# OBS! to get it filled by row instead of column, rearrange the plots
pl2 <- pl[as.vector(t(matrix(names(pl),nrow=6)))]

pdf("../data/figures2/violins_integrated_woSS2.pdf", width = 14, height = 7)
grid.arrange(grobs=pl2, nrow=6)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### Umap plots of gene list

Gene list from Katarina with genes (gene\_list\_NatureMedicine.csv)

``` r
genes2 <- read.table("../data/figures2/gene_list_NatureMedicine.csv", header = T)

genes2$symbol <- unlist(lapply(strsplit(as.character(genes2$Gene),":"), function(x) x[1]))

# is all the same genes as above except 5
# [1] "PCDH15"  "SLC17A6" "SLC17A8" "GAD2"    "SLC32A1"
m <- match(genes2$symbol, unlist(genes))
genes2[which(is.na(m)),"symbol"]
```

    ## [1] "PCDH15"  "SLC17A6" "SLC17A8" "GAD2"    "SLC32A1"

``` r
data.integrated@active.assay <- "RNA"
plots <- FeaturePlot(data.integrated, features = genes2$symbol, reduction = "umap", combine=F, cols = c("yellow","red","black"))
plots <- lapply(plots, function(x) x + NoAxes() + theme(legend.text = element_text(size=6)))

pdf("../data/figures2/umap_gene_plots_woSS2.pdf", width = 10, height = 7)
for (i in 1:4){
  sel <- (i*9-8):min(i*9,length(plots))
  grid.arrange(grobs=plots[sel])
}
dev.off()
```

    ## quartz_off_screen 
    ##                 2

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
