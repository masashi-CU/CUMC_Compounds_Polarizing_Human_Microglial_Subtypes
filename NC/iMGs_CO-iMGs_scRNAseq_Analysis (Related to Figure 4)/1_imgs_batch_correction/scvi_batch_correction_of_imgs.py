###############################################################################|
# Pipeline for correcting batch effects across iMG samples using a scVI model
#             Natacha Comandante-Lou (nc3018 at columbia dot edu)
# Input: h5ad file of the raw count data
# Output: AnnData object with the "batch_corrected_counts" layer
# Source: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/api_overview.html
##############################################################################|
import os
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import scprep
import scvi

def scvi_batch_correct_pipeline(adata, n_top_genes = 3000, batch_key = "batch", seed = 888,n_layers = 2, n_latent =30, model_init = None, outdir = "output", save = True):
  ###############################################################################|
  # PREPROCESSING
  ##############################################################################|
  # Normalize the data ---
  # save the raw counts to a layer
  adata.layers["counts"] = adata.raw.X.copy()
  sc.pp.normalize_total(adata, target_sum=1e4)
  sc.pp.log1p(adata)
  adata.raw = adata  # freeze the state in `.raw`
  
  # Highly variable genes ---
  sc.pp.highly_variable_genes(adata, flavor = 'seurat', n_top_genes = n_top_genes, batch_key=batch_key)
  
  ##############################################################################|
  # BATCH CORRECTION
  ##############################################################################|
  # Set-up
  
  scvi.settings.seed = seed
  
  use_gpu = True
  
  if not use_gpu:
    n_cores = 1
    scvi.settings.num_threads = 1
    scvi.settings.dl_num_workers = n_cores
  
  if model_init is None:
    # Load processed data
    raw_counts_layer = 'counts' # raw UMI counts layer
    batch_col = batch_key # metadata colum specifying batches
    scvi.model.SCVI.setup_anndata(adata, layer = raw_counts_layer, batch_key = batch_col)
    model = scvi.model.SCVI(adata, n_layers = n_layers, n_latent = n_latent, gene_likelihood= "nb")   
    model.train()
    
  else:
    
    model = SCVI.load(model_init)
    
  
  corrected_data = model.get_normalized_expression(transform_batch = sorted(adata.obs[batch_col].unique()),
                                                   library_size = 1e4)
  adata.layers['batch_corrected_counts'] = corrected_data
  
  corrected_data.iloc[:,:] = np.log1p(corrected_data.values)
  adata.layers['log_batch_corrected_counts'] = corrected_data
  
  latent = model.get_latent_representation()
  adata.obsm["X_scVI"] = latent
  
  if save:
    #adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    adata.write_h5ad(os.path.join(outdir, 'scvi_corrected_counts.h5ad'), compression='gzip')
  
  model.save(outdir+"scvi_model/")
  return adata





##################
# RUN
##################
os.chdir('analysis_path')
input_data_path = 'data/h5ad/'
output_data_path  = 'scvi_output/'

# Load data-----
adata = sc.read_h5ad(os.path.join(input_data_path, 'img_filtered_rm_dbl.h5ad'))
adata.obs['orig.ident'] = adata.obs['orig.ident'].astype("str") +"-"+ adata.obs['treatment'].astype("str")
adata.obs['orig.ident'] = adata.obs['orig.ident'].astype("category")
adata.obs["batch"] = adata.obs["batch"].astype("category")

adata = scvi_batch_correct_pipeline(adata, n_top_genes = 3000, batch_key = "batch", seed = 888,n_layers = 2, n_latent =30, model_init = None, outdir = output_data_path, save = True)
