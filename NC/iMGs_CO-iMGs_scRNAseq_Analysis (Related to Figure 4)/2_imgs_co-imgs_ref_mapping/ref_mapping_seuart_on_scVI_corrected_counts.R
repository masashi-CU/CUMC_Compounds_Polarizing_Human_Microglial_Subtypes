##########################################################################|
# Mapping CO-iMGs treated with DMSO or Camptothecin (0.1 uM and 0.5 uM)
# onto iMGs reference
# Natacha Comandante-Lou (nc3018 at columbia dot edu)
##########################################################################|
setwd("analysis_path")

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(zellkonverter)
library(scCustomize)
library(ggpubr)
library(logger)

source("seurat_reference_mapping.R")
# Load Data ---------------------------
#ref
ad <- readH5AD('scvi_output/scvi_corrected_counts.h5ad')
img<- as.Seurat(ad, counts = "X", data = "batch_corrected_counts") # convert to a seurat obj
saveRDS(img, file = "scvi_output/scvi_corrected_counts.rds")

#query
mic = LoadH5Seurat("data/mic_sct_clus.h5Seurat")
DimPlot_scCustom(mic)



# Project CO-iMGs on iMGs data(downsampled) ---------------------------
set.seed(124)
img.sub = img[, sample(colnames(img), size = 20000, replace=F)]

mic.sub = subset(x = mic, subset = treatment %in% c("DMSO","Camptothecin"))

nfeatures =1000
npcs = 15

output.dir = sprintf("ref_mapping_output/bcc_ref_mapping_all_n=20000/Camptothecin nhvgs = %d, npcs = %d",nfeatures, npcs)

output =seurat_reference_mapping_pipeline(ref = img.sub ,query = mic.sub,
                                          ref.assay = "originalexp",
                                          query.assay = "RNA",
                                          nfeatures = nfeatures,
                                          npcs = npcs,
                                          pred.var = "originalexp_snn_res.0.7",
                                          resolution = 0.7,
                                          reference.reduction = "pca",
                                          normalization.method = "LogNormalize",
                                          reduction.model = "umap",
                                          plot = TRUE,
                                          plot.reduction = "umap",
                                          plot.split.by = "treatment",
                                          output.dir = output.dir,
                                          save = TRUE)
# refquery = output$refquery
# saveRDS(refquery,file = file.path(output.dir, "refquery.rds"))
p = DimPlot_scCustom(subset(output$refquery, subset=treatment%in%c("DMSO","Camptothecin")),reduction = "umap", split.by = "id",group.by = "treatment")
ggsave(plot = p, filename = file.path(output.dir,"dim_reduc","umap_refquery_dmso-campto_only.png"),width = 16, height = 7)

# save as h5ad for meld
refquery = output$refquery
SaveH5Seurat(refquery, filename = file.path(output.dir,"refquery.h5Seurat"),overwrite = TRUE)
Convert(file.path(output.dir,"refquery.h5Seurat"), dest = "h5ad",overwrite = TRUE)



