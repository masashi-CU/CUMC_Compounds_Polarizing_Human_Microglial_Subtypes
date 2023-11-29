###############################################################################|
# Run MELD analysis on the Mouse microglia data (after batch correction by sex)
# Natacha Comandante-Lou (nc3018 at columbia dot edu)
###############################################################################|
setwd("analysis_path")
library(tidyverse)
library(reticulate)
library(Seurat)
library(SeuratDisk)
library(scCustomize)
library(dplyr)
library(ggplot2)
library(zellkonverter)
library(scCustomize)
library(ggpubr)
library(rstatix)
library(ggbeeswarm)
library(RColorBrewer)
library(schex)
library(logger)
use_condaenv("/home/nc3018/conda/my-envs/sctools-3.8")


pd = import("pandas")
np = import("numpy")
sc = import("scanpy")
os = import("os")
scprep = import("scprep")

###############################################################################|
# MELD PARAMETERS ------
###############################################################################|
# Dataset options: Male (M), Female (F)
opt = "M" 

# Meld Params:
npcs = 15
batch_key = "replicate" 
ctrl_key = "DMSO" #control label

meld_output_dir = file.path("meld_analysis_output_downsampled_mouse",opt)
###############################################################################|
# Load Data and Assemble Seurat Input ------
###############################################################################|
# Assemble Seurat Object using batch corrected counts

scvi_data_path = "scvi_output_downsampled_M_mouse_LogNormalize/"
ad <- readH5AD(os$path$join(scvi_data_path, 'scvi_corrected_counts.h5ad'))
adata_Seurat <- as.Seurat(ad, counts = "X", data = "batch_corrected_counts")

p = DimPlot_scCustom(adata_Seurat,group.by = "identity",reduction = "X_umap",split.by = "treatment")
ggsave(p, filename = file.path(meld_output_dir,"dim_reduc","umap_before_bcc_split.by=treatment_group.by=identity.png"), width = 12,height = 3) 

ndim = 20
sobj = RunUMAP(adata_Seurat, dims = 1:ndim, reduction = "X_scVI")
p = DimPlot_scCustom(sobj,group.by = "identity",reduction = "umap",split.by = "treatment")
ggsave(p, filename = file.path(meld_output_dir,"dim_reduc",sprintf("umap_after_bcc_split.by=treatment_group.by=identity_dim=%d.png",ndim)), width = 12,height = 3) 



###############################################################################|
# PCA. -------
###############################################################################|

# read h5ad
adata = sc$read_h5ad(file.path(scvi_data_path,'scvi_corrected_counts.h5ad'))

# modify treatment label
metadata = adata$obs
colnames(metadata)[colnames(metadata)=="gender"] = "sex"
metadata =dplyr::mutate(metadata, replicate = paste(orig.ident, sex, sep = "."))
adata$obs = metadata
metadata$orig.ident = as.factor(metadata$orig.ident)

layer = 'batch_corrected_counts'
log_info("Using layer: ", layer)
data_bcc = adata$layers[layer]
data_sqrt = np$sqrt(data_bcc)

data_pca = scprep$reduce$pca(data_sqrt,n_components=as.integer(100))
adata$obsm['bcc_pca'] = data_pca #store pca calculated from batch corrected count back to the anndata obj




feature.list_mouse = readRDS("module_lists_mouse.rds")
#######################################################################|
# Run MELD Pipeline ---------
#######################################################################|

source('meld_analysis_pipeline.R')

condition_key = "identity"
# knn = 5
# beta = 44

mypal_drug = c("#000000","#FFBF80","#A2D1DE","#9180ED","#EBE76C","#B71375","#FF8080")

adata_Seurat = meld_pipeline(adata,
                             sobj,
                             
                             # MELD params
                             meld = TRUE,
                             npcs = npcs,
                             condition_key = condition_key,
                             batch_key = batch_key,
                             ctrl_key = ctrl_key, #control label
                             param_sweep = TRUE,
                             # knn = knn,
                             # beta = beta,
                             metadata = metadata,
                             pca_key = "bcc_pca",
                             
                             # MELD post-process params
                             meld_postproc = TRUE,
                             meld_output_dir = meld_output_dir,
                             
                             # Params for seuart_pipeline
                             run_seurat_pipeline = FALSE,
                             normalization = "LogNormalize",
                             n.var.features = 3000,
                             umap = FALSE, #recalculate umap of adata_Seurat
                             umap.assay = "originalexp",
                             umap.ndim = NULL,
                             umap.reduction = "pca",
                             umap.nn = 30L,
                             umap.md = 0.3,
                             umap.metric = "cosine",
                             umap.plot = TRUE,
                             save.obj = FALSE,
                             
                             # Params for modules calculations
                             module.score = TRUE,
                             module.assay = "originalexp",
                             feature.list = feature.list_mouse,
                             module.method = "seurat",
                             module.plot.reduction = "umap",

                             # Params for umap plot
                             plot.group.by = "orig.ident",
                             plot.assay = "originalexp",
                             plot.colors.use = mypal_drug,
                             schex = TRUE,
                             save.result = TRUE
)




