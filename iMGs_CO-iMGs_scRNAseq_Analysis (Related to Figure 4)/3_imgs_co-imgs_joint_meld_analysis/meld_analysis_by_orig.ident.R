###############################################################################|
# Run Meld analysis on the integrated CO-iMGs and iMGs data
# Source: https://github.com/KrishnaswamyLab/MELD/blob/main/notebooks/MELD_Quickstart.ipynb
# Natacha Comandante-Lou (nc3018 at columbia dot edu)
###############################################################################|
setwd("analysis_path")
library(reticulate)
library(logger)
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(scCustomize)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggbeeswarm)
use_condaenv("/home/nc3018/conda/my-envs/sctools-3.8")

###############################################################################|
# Load python functions
###############################################################################|
log_info("Loading Python Functions...")
os = import("os")
ad = import("anndata")
sc = import("scanpy")
pd = import("pandas")
np = import("numpy")
scprep = import("scprep")

source_python("run_meld.py")

###############################################################################|
# Set up directories ------
###############################################################################|

# Meld output directories
opt = 'condition_key=orig.ident'
condition_key = 'orig.ident'
output_dir = file.path("meld_analysis_output",opt,"meld_output")
figure_dir  = file.path("meld_analysis_output",opt,"meld_figures")


# Meld output subfolders
npcs = as.integer(15) #number of pcs used as input to meld
dat_dir = file.path(output_dir,sprintf("npc=%d",npcs))
fig_dir = file.path(figure_dir,sprintf("npc=%d",npcs))

if(!dir.exists(dat_dir)){dir.create(dat_dir, recursive = "TRUE")}
if(!dir.exists(fig_dir)){dir.create(fig_dir, recursive = "TRUE")}
###############################################################################|
# Log -------
###############################################################################|

log_layout(layout_simple)
log_file = file.path(dat_dir , "meld.log")
if(file.exists(log_file)){
  log_warn("Overwriting previous log: ", log_file)
  file.remove(log_file)
}

log_appender(appender_tee(file = log_file))

log_info("#--------------------------------------------------------------------#")
log_info("#                        MELD ANALYSIS                               #")
log_info("#--------------------------------------------------------------------#")
###############################################################################|
# Load Data ------------------
###############################################################################|
start = Sys.time()
## Mapped data source -----
mapped_data_path = sprintf("ref_mapping_output/bcc_ref_mapping_all_n=20000/Camptothecin nhvgs = 1000, npcs = %d", npcs)
log_info("Reading Data from: ", file.path(mapped_data_path,"refquery.h5ad"))
refquery = readRDS(file.path(mapped_data_path,"refquery.rds"))

## Read adata ------------------
# convert to adata
SaveH5Seurat(refquery, filename = file.path(mapped_data_path,"refquery.h5Seurat"),overwrite = TRUE)
Convert(file.path(mapped_data_path,"refquery.h5Seurat"), dest = "h5ad",overwrite = TRUE)
adata = sc$read_h5ad(file.path(mapped_data_path,"refquery.h5ad"))

# modify treatment label
#treatment = adata$obs['treatment']
metadata = adata$obs
metadata$batch[metadata$batch=="NA"] = "3"
metadata = dplyr::mutate(metadata, treatment = case_when(orig.ident %in% c("PV028-Camptothecin","PV011-Camptothecin") ~ "Camptothecin (0.1 uM)",
                                                  orig.ident %in% c("PV029-Camptothecin") ~ "Camptothecin (0.5 uM)",
                                                  orig.ident %in% c("PV020-Luotonin A") ~ "LuotoninA",
                                                  TRUE ~ treatment))

## Run Meld Pipeline ---------
log_info("Running Meld Pipeline...")


log_info("Number of PCs selected as MELD input: ", npcs)


data_pca = adata$obsm[["X_pca"]][,1:npcs]

log_info("Parameter Sweep... ")

top_result = meld_param_sweep(data = data_pca ,knn_range = as.integer(np$arange(1,28)),beta_range = as.integer(np$arange(1,200)),output_dir = fig_dir)

knn = as.integer(top_result['knn'])
beta = as.integer(top_result['beta'])

adata_meld = run_meld(data = data_pca, knn = knn, beta = beta , metadata=metadata, output_dir = dat_dir, fig_dir = fig_dir ,
                      condition_key = condition_key,
                      batch_key = "batch", ctrl_key = "DMSO")
#-------------------------#
end = Sys.time()
log_success("Elapsed Time (Meld Analysis): ", round(end-start, 2) ," ", units(end-start))

