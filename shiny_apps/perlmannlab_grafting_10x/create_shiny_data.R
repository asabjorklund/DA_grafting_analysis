library(Seurat)
library(Matrix)
library(data.table)

infile <- "/Users/asbj/projects/Perlmann/analysis/grafting/10x/data/processed2/integrated/seurat_integrated_rat_scale_res0.1.Rdata"
load(infile)


# add in celltype info with cell names
# add in FACS/notFACS


# reassing cluster names to celltypes
celltypes <- c("Astrocyte","VLMC","Neuron","Cycling")
names(celltypes) <- as.character(0:3)
ct <- celltypes[data.integrated@active.ident]
data.integrated <- AddMetaData(data.integrated, ct, col.name = "CelltypeCluster")


# need to save metadata + expression.

meta <- data.integrated@meta.data
umap <- data.integrated@reductions$umap@cell.embeddings
tsne <- data.integrated@reductions$tsne@cell.embeddings
expr <- data.integrated@assays$RNA@data

meta$tsne_x <- tsne[,1]
meta$tsne_y <- tsne[,2]
meta$umap_x <- umap[,1]
meta$umap_y <- umap[,2]


# select only expressed genes for drop-down menu (expressed in >1 cell)
rs <- rowSums(data.integrated@assays$RNA@counts)
Gene_names <- names(rs)[rs>1]
genes <- as.data.table(Gene_names)

expr <- expr[Gene_names,]


savefile <- "integrated_rat_scale_res0.1.Rdata"
save(meta,umap,tsne,expr,genes, file = savefile)


# make one version with labels in plot
#scaleC <- scale_colour_manual(values = c("blue1", "gray","green3", "red1"))
