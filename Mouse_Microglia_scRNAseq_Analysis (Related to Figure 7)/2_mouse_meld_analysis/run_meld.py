###############################################################################|
# Wrapper function for MELD
# Source: https://github.com/KrishnaswamyLab/MELD/blob/main/notebooks/MELD_Quickstart.ipynb
# Natacha Comandante-Lou (nc3018 at columbia dot edu)
###############################################################################|
import os as os
import numpy as np
import pandas as pd
import graphtools
import graphtools as gt
import matplotlib.pyplot as plt
import scprep
import scanpy as sc
import phate
import sklearn
import meld
import tasklogger
import collections
import warnings
from collections import defaultdict
import anndata as ad
import sklearn.mixture
from scipy.spatial.distance import pdist, cdist, squareform
warnings.simplefilter("ignore")
import pickle
np.__version__

##############################################################################|
# Run MELD ---------------
##############################################################################|

def run_meld(data, knn, beta, metadata, graph=None, output_dir = "data",fig_dir = "figures", condition_key = "orig.ident", batch_key = "batch", ctrl_key = "DMSO" ):
  # Build graph - learned from all cells from all samples. 
  if graph is None:
    G = gt.Graph(data, knn=int(knn), use_pygsp=True)
  else:
    G = pickle.load(open(graph,"rb"))
  
  # Meld estimate the density of each sample.
  meld_op = meld.MELD(beta=beta)
  sample_densities = meld_op.fit_transform(G, sample_labels=metadata[condition_key])
  
  pickle.dump(G, open(output_dir+'/graph_sqrt.pkl', 'wb'))
  
  ## Normalization ------
  # normalize to control within each batch
  def ctrl_normalize_densities_by_replicate(sample_densities, ctrl_label, replicate):
    sample_likelihoods = sample_densities.copy()
    replicates = np.unique(replicate)
    for rep in replicates:
      conditions = np.unique(metadata[replicate==rep][condition_key])
      for con in conditions:
        ctrl = conditions[pd.Series(conditions).str.contains(ctrl_label)].all()
        curr_cols = [con,ctrl]
        print("ctrl normalize --------")
        print(curr_cols)
        sample_likelihoods[con] = sklearn.preprocessing.normalize(sample_densities[curr_cols], norm='l1')[:,0]
        #print( sample_likelihoods[curr_cols])
    return sample_likelihoods
  
  ctrl_normalized_sample_likelihoods = ctrl_normalize_densities_by_replicate(sample_densities,ctrl_key, metadata[batch_key])
  
  # normalize within each batch
  def normalize_densities_by_replicate(sample_densities, replicate):
    sample_likelihoods = sample_densities.copy()
    replicates = np.unique(replicate)
    for rep in replicates:
      curr_cols = np.unique(metadata[replicate==rep][condition_key])
      print("normalize -------------")
      print(curr_cols)
      sample_likelihoods[curr_cols] = sklearn.preprocessing.normalize(sample_densities[curr_cols], norm='l1')
      #print( sample_likelihoods[curr_cols])
    return sample_likelihoods
  
  normalized_sample_likelihoods = normalize_densities_by_replicate(sample_densities,metadata[batch_key])
  
  adata_meld = ad.AnnData(data)
  adata_meld.uns['meld_densities'] = sample_densities
  adata_meld.uns['meld_ctrl_normalized_densities'] = ctrl_normalized_sample_likelihoods
  adata_meld.uns['meld_normalized_densities'] = normalized_sample_likelihoods
  
  def gate_distr(density_name, condition_oi, plot):
# Fit meld score distribution as a gaussian mixture model
    mixture_model = sklearn.mixture.GaussianMixture(n_components=3)
    classes = mixture_model.fit_predict(adata_meld.uns[density_name][condition_oi].values.reshape(-1,1))
    # Reorder the class labels so that 0 has the lowest average likelihood and 2 has the highest.
    classes = scprep.utils.sort_clusters_by_values(classes, adata_meld.uns[density_name][condition_oi])
    if plot:
      density_fig_dir = fig_dir+'/meld_density_gate_ncomp=3'
      if not os.path.exists(density_fig_dir):
        os.makedirs(density_fig_dir)
      # Plot
      fig, ax = plt.subplots(1, figsize=(5,4))
      for c in np.unique(classes):
          ax.hist(adata_meld.uns[density_name][condition_oi][classes == c], range=(0,1), bins=100,
                 color=meld.utils.get_meld_cmap()([0,255//2, 255])[c])
      
      fig_title = density_name + "\n"+condition_oi
      fig_name = density_fig_dir+'/'+density_name +'_'+condition_oi +'.png'
      plt.title(fig_title)
      fig.tight_layout()
      fig.savefig(fig_name.encode('utf-8'), dpi=300)
    return classes
  
  condition_list = adata_meld.uns['meld_densities'].columns
  
  density_name = 'meld_ctrl_normalized_densities'
  res = pd.DataFrame([gate_distr(density_name, c, True) for c in condition_list],condition_list).transpose()
  adata_meld.uns[density_name+"_class"] = res
  
  density_name = 'meld_normalized_densities'
  res = pd.DataFrame([gate_distr(density_name, c, True) for c in condition_list],condition_list).transpose()
  adata_meld.uns[density_name+"_class"] = res
  adata_meld.write_h5ad(os.path.join(output_dir, 'meld_output.h5ad'))
  return adata_meld

##############################################################################|
# MELD PARAMETER SWEEP ----------
##############################################################################|
def meld_param_sweep(data,knn_range ,beta_range ,output_dir = "figures"):
  #os.mkdir(output_dir)
  benchmarker = meld.Benchmarker()
  
  # 3D PHATE components are used to create the ground truth PDF
  benchmarker.fit_phate(data);
  
  
  from joblib import Parallel, delayed
  
  def simulate_pdf_calculate_likelihood(benchmarker, seed, beta):
    benchmarker.set_seed(seed)
    benchmarker.generate_ground_truth_pdf()
    
    benchmarker.generate_sample_labels()
    benchmarker.calculate_MELD_likelihood(beta=beta)
    MELD_mse = benchmarker.calculate_mse(benchmarker.expt_likelihood)
    return MELD_mse, seed, beta, benchmarker.graph.knn
  
  results = []
  adata = ad.AnnData(data)
  
  with Parallel(n_jobs=36) as p:
      for knn in knn_range:
        # doing this outside the parallel loop because building the graph takes the longest
        benchmarker.fit_graph(adata.X, knn=knn)
        print(knn)
        curr_results = p(delayed(simulate_pdf_calculate_likelihood)(benchmarker, seed, beta) \
                                       for seed in range(25) for beta in beta_range)
        curr_results = pd.DataFrame(curr_results, columns = ['MSE', 'seed', 'beta', 'knn'])
        results.append(curr_results)
  
  results = pd.concat(results, axis=0)
  results.to_pickle(os.path.join(output_dir, "meld_param_sweep.pkl"))
  
  # 
  # # Load results 
  # results= pd.read_pickle('data_test/meld_param_sweep_knn=40.pkl')
  # 
  # 
  # We want to take the average of each set of random seeds for each combination of beta and knn values
  results_wide = results.groupby(['beta', 'knn']).mean().sort_values(by='MSE').reset_index()
  
  scprep.plot.scatter(results_wide['beta'], results_wide['knn'],
                           s=1, c=results_wide['MSE'], vmax=0.01, vmin=0, cmap='inferno_r')
  plt.show()
  
  # Highlight the top performing combination with a large red dot
  top_result = results_wide.sort_values('MSE').iloc[0]
  plt.scatter(top_result['beta'], top_result['knn'], c='r', s=40, linewidth=1, edgecolor='k')
  plt.show()
  top_result
  
  fig = plt.gcf()
  fig.set_size_inches(5,4)
  fig.savefig(os.path.join(output_dir, "meld_param_sweep.png"), dpi=300)
  
  return top_result


