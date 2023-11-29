#' MELD pipline:
#'  - R wrapper for run_meld.py
#'  - MELD postanalysis
#'  Natacha Comandante-Lou (nc3018 at columbia dot edu)
require(reticulate)
require(tidyverse)
require(logger)
require(dplyr)
require(Seurat)
require(SeuratDisk)
require(scCustomize)
require(RColorBrewer)
require(ggplot2)
require(ggpubr)
require(rstatix)
require(ggbeeswarm)
require(zellkonverter)
require(schex)
###############################################################################|
# Load python functions
###############################################################################|
log_info("Loading Python Functions...")
use_condaenv("/home/nc3018/conda/my-envs/sctools-3.8")
os = import("os")
ad = import("anndata")
sc = import("scanpy")
pd = import("pandas")
np = import("numpy")
scprep = import("scprep")

source_python("run_meld.py")


###############################################################################|
# MELD PIPELINE-----------
meld_pipeline = function(adata,
                         adata_Seurat,
                         
                         # MELD params
                         meld = TRUE,
                         npcs = 15,
                         condition_key,
                         batch_key,
                         ctrl_key, #control label
                         param_sweep = FALSE,
                         knn,
                         beta,
                         metadata,
                         pca_key = "X_pca", #key in adata.obsm storing pca
                         # MELD post-process params
                         meld_postproc = TRUE,
                         meld_output_dir,
                         # Params for seuart_pipeline
                         run_seurat_pipeline = TRUE,
                         normalization = "LogNormalize",
                         n.var.features = 3000,
                         umap = TRUE, #recalculate umap of adata_Seurat
                         umap.assay = "RNA",
                         umap.ndim = NULL,
                         umap.reduction = "pca",
                         umap.nn = 30L,
                         umap.md = 0.3,
                         umap.metric = "cosine",
                         umap.plot = TRUE,
                         save.obj = FALSE,
                         
                         # Params for modules calculations
                         module.score = FALSE,
                         feature.list,
                         module.method = "seurat",
                         module.plot.reduction = "umap",
                         module.assay = "RNA",
                         
                         # Params for umap plot
                         plot.group.by = "orig.ident",
                         plot.assay = "RNA",
                         plot.colors.use,
                         schex = TRUE,
                         save.result = FALSE
){
  #-----------------------------------------------------------------------------|
  # Set up directories ------
  #-----------------------------------------------------------------------------|
  
  # Meld output directories
  opt = sprintf("batch_key = %s, condition_key = %s", batch_key, condition_key)
  
  output_dir = file.path(meld_output_dir,opt,"meld_output")
  figure_dir  = file.path(meld_output_dir,opt,"meld_figures")
  
  
  # Meld output subfolders
  npcs = as.integer(npcs) #number of pcs used as input to meld
  dat_dir = file.path(output_dir,sprintf("npc=%d",npcs))
  fig_dir = file.path(figure_dir,sprintf("npc=%d",npcs))
  
  if(!dir.exists(dat_dir)){dir.create(dat_dir, recursive = "TRUE")}
  if(!dir.exists(fig_dir)){dir.create(fig_dir, recursive = "TRUE")}
  
  #-----------------------------------------------------------------------------|
  # Log -------
  #-----------------------------------------------------------------------------|
  
  log_layout(layout_simple)
  log_file = file.path(dat_dir , "meld.log")
  if(file.exists(log_file)){
    log_warn("Overwriting previous log: ", log_file)
    file.remove(log_file)
  }
  
  log_appender(appender_tee(file = log_file))
  
  
  #-----------------------------------------------------------------------------|
  # Run MELD -------
  #-----------------------------------------------------------------------------|
  if (meld){
    start = Sys.time()
    log_info("#--------------------------------------------------------------------#")
    log_info("#                        MELD ANALYSIS                               #")
    log_info("#--------------------------------------------------------------------#")
    
    ## Run Meld Pipeline ---------
    log_info("Running Meld Pipeline...")
    
    
    log_info("Number of PCs selected as MELD input: ", npcs)
    
    
    data_pca = adata$obsm[[pca_key]][,1:npcs]
    
    if (param_sweep){
      log_info("Parameter Sweep... ")
      top_result = meld_param_sweep(data = data_pca ,knn_range = as.integer(np$arange(1,28)),beta_range = as.integer(np$arange(1,200)),output_dir = fig_dir)
      
      knn = as.integer(top_result['knn'])
      beta = as.integer(top_result['beta'])
      
      log_info("Top Results:")
      log_info("knn = ", knn)
      log_info("beta = ", beta)
      
    }else{
      log_info("No parameter sweep, using the following parameters:")
      
      knn = as.integer(knn)
      beta = as.integer(beta)
      
      log_info("knn = ", knn)
      log_info("beta = ", beta)
    }
    
    
    adata_meld = run_meld(data = data_pca, knn = knn, beta = beta , metadata=metadata, output_dir = dat_dir, fig_dir = fig_dir ,
                          condition_key = condition_key,
                          batch_key = batch_key, ctrl_key = ctrl_key)
    #-------------------------#
    end = Sys.time()
    log_success("Elapsed Time (Meld Analysis): ", round(end-start, 2) ," ", units(end-start))
    
  }
  
  
  #-----------------------------------------------------------------------------|
  # Meld Post-processing ------------------
  #-----------------------------------------------------------------------------|
  if (meld_postproc){
    meld_data_path = dat_dir
    meld_fig_path = fig_dir
    
    adata_Seurat = meld_postprocessing(adata_Seurat,
                                       meld_data_path, #path to meld_output.h5ad
                                       meld_fig_path,
                                       
                                       # Params for seuart_pipeline
                                       run_seurat_pipeline = run_seurat_pipeline,
                                       normalization = normalization ,
                                       n.var.features = n.var.features,
                                       umap = umap, #recalculate umap of adata_Seurat
                                       umap.assay = umap.assay,
                                       umap.ndim =umap.ndim,
                                       umap.reduction = umap.reduction,
                                       umap.nn = umap.nn,
                                       umap.md = umap.md,
                                       umap.metric = umap.metric,
                                       umap.plot = umap.plot,
                                       save.obj = save.obj,
                                       
                                       # Params for modules calculations
                                       module.score = module.score,
                                       feature.list = feature.list,
                                       module.method = module.method,
                                       module.plot.reduction = module.plot.reduction,
                                       module.assay = module.assay,
                                       # Params for umap plot
                                       plot.group.by = plot.group.by,
                                       plot.assay = plot.assay,
                                       plot.colors.use = plot.colors.use,
                                       
                                       schex = schex,
                                       # Params for Meld Visualizations
                                       meld.assay = c("meld_ctrl_norm_density","meld_norm_density"))
    
    
  }
  
  
  if(save.result){
    saveRDS(adata_Seurat, file = file.path(dat_dir,"meld_res_sobj.rds"))
  }
  return( adata_Seurat)
}





###############################################################################|
# MELD POST-PROCESSING ------------------
###############################################################################|
meld_postprocessing = function(adata_Seurat,
                               meld_data_path, #path to meld_output.h5ad
                               meld_fig_path,
                               
                               # Params for seuart_pipeline
                               run_seurat_pipeline = TRUE,
                               normalization = "LogNormalize",
                               n.var.features = 3000,
                               umap = TRUE, #recalculate umap of adata_Seurat
                               umap.assay = "RNA",
                               umap.ndim = NULL,
                               umap.reduction = "pca",
                               umap.nn = 30L,
                               umap.md = 0.3,
                               umap.metric = "cosine",
                               umap.plot = TRUE,
                               save.obj = FALSE,
                               
                               # Params for modules calculations
                               module.score = TRUE,
                               module.assay = "RNA",
                               feature.list = feature.list,
                               module.method = "seurat",
                               module.plot.reduction = "umap",
                               
                               # Params for umap plot
                               plot.group.by = "orig.ident",
                               plot.assay = "RNA",
                               plot.colors.use,
                               schex = FALSE,
                               # Params for Meld Visualizations
                               meld.assay = c("meld_ctrl_norm_density","meld_norm_density")
                               
){
  
  log_info("#------------------------------------------------------------------#")
  log_info("#                     MELD POST-PROCESSING                         #")
  log_info("#------------------------------------------------------------------#")
  
  
  log_info("Meld Post-processing from:\n ", meld_data_path)
  log_info("Meld Figures to:\n ",meld_fig_path)
  #-----------------------------------------------------------------------------|
  ## Load and Add MELD Output ------------------
  #-----------------------------------------------------------------------------|
  ad <- readH5AD(os$path$join(meld_data_path, 'meld_output.h5ad'))
  
  log_info("Adding meld data to input seurat objects...")
  # add object
  meld_mat = as.matrix(ad@metadata[["meld_densities"]])
  rownames(meld_mat) = rownames(adata_Seurat@meta.data)
  adata_Seurat[["meld_density"]] <- CreateAssayObject(data  = t(meld_mat))
  
  # add object
  meld_mat = as.matrix(ad@metadata[["meld_normalized_densities"]])
  rownames(meld_mat) = rownames(adata_Seurat@meta.data)
  adata_Seurat[["meld_norm_density"]] <- CreateAssayObject(data  = t(meld_mat)) #normalized density within the same batch
  
  # add object
  meld_mat = as.matrix(ad@metadata[["meld_ctrl_normalized_densities"]])
  rownames(meld_mat) = rownames(adata_Seurat@meta.data)
  adata_Seurat[["meld_ctrl_norm_density"]] <- CreateAssayObject(data  = t(meld_mat)) #normalized density with ctrl within the same batch
  
  # add object
  # meld_densities = adata_Seurat@assays[["meld_ctrl_norm_density"]]@data%>%t()%>%as.data.frame()
  # treatment_list = unique(adata_Seurat@meta.data$treatment)
  # meld_densities_mean = sapply(treatment_list, function(x) rowMeans(meld_densities[, grep(x, names(meld_densities))]))
  # 
  # adata_Seurat[["mean_meld_ctrl_norm_density"]] = CreateAssayObject(data = t(meld_densities_mean))
  
  
  
  # add metadata
  # three class
  meld_class_mat = as.data.frame(ad@metadata[["meld_ctrl_normalized_densities_class"]])
  rownames(meld_class_mat) = rownames(adata_Seurat@meta.data)
  adata_Seurat = AddMetaData(adata_Seurat,meld_class_mat,col.name = paste(colnames(meld_class_mat),"-three_ctrl_norm_class",sep = "") )
  # two class (combined low-med)
  col_oi = colnames(adata_Seurat@meta.data%>%select(matches("three_ctrl_norm_class")))
  for(class_oi in col_oi){
    new_col = stringr::str_replace(class_oi, "three","two")
    adata_Seurat@meta.data[[new_col]] = adata_Seurat@meta.data[[class_oi]]
    adata_Seurat@meta.data[[new_col]][adata_Seurat@meta.data[[new_col]]==1] =0
  }
  
  
  # add metadata
  #three class
  meld_class_mat = as.data.frame(ad@metadata[["meld_normalized_densities_class"]])
  rownames(meld_class_mat) = rownames(adata_Seurat@meta.data)
  adata_Seurat = AddMetaData(adata_Seurat,meld_class_mat,col.name = paste(colnames(meld_class_mat),"-three_norm_class",sep = "") )
  # two class (combined low-med)
  col_oi = colnames(adata_Seurat@meta.data%>%select(matches("three_norm_class")))
  for(class_oi in col_oi){
    new_col = stringr::str_replace(class_oi, "three","two")
    adata_Seurat@meta.data[[new_col]] = adata_Seurat@meta.data[[class_oi]]
    adata_Seurat@meta.data[[new_col]][adata_Seurat@meta.data[[new_col]]==1] =0
  }
  
  #-----------------------------------------------------------------------------|
  ## UMAP Viz ------------------
  #-----------------------------------------------------------------------------|
  log_info("Calculating UMAP for visualizations..")
  tryCatch({
    if (run_seurat_pipeline){
  
      log_info("Running Seurat Pipeline on adata_Seurat...")
      adata_Seurat  = seurat_pipeline(data = adata_Seurat,
                                      normalization = normalization,
                                      normalization.assay = "RNA",
                                      vars.to.regress = NULL,
                                      n.var.features = n.var.features,
                                      cluster = FALSE,
                                      
                                      pca = TRUE,
                                      pca.plot = TRUE,
                                      max_npcs = 50, 
                                      pct.sd.thres = 0.1,
                                      
                                      umap = umap,
                                      umap.assay = umap.assay,
                                      umap.ndim = umap.ndim,
                                      umap.reduction = umap.reduction,
                                      umap.nn = umap.nn,
                                      umap.md = umap.md,
                                      umap.metric = umap.metric,
                                      umap.plot = umap.plot,
                                      
                                      plot.group.by = plot.group.by,
                                      plot.colors.use = plot.colors.use,
                                      plot.assay = plot.assay,
                                      plot.type = "DimPlot",
                                      
                                      add.module = FALSE,
                                      
                                      out.dat.dir = file.path(meld_data_path,"seurat_pipeline"),
                                      out.fig.dir = file.path(meld_fig_path,"seurat_pipeline"),
                                      
                                      save.obj = save.obj,
                                      verbose = TRUE)
    }
    
    
    p = DimPlot_scCustom(adata_Seurat,reduction = "umap", group.by = plot.group.by,colors_use = plot.colors.use)
    ggsave(p, filename = file.path(meld_fig_path,"seurat_pipeline","dim_reduc",sprintf("umap_%s.png",plot.group.by)))
  
  },error = function(e){
    log_error("Error in seurat_pipeline: ",message(e))
    return(adata_Seurat)
    })
  #-----------------------------------------------------------------------------|
  # Calculate Module Score -----
  #-----------------------------------------------------------------------------|
  tryCatch({
    if (module.score){
      source("/mnt/mfs/ctcn/team/natacha/RScript/seurat_get_summary_scores.R")
      log_info("Calculating module scores...")
      adata_Seurat = seurat_get_summary_scores(data = adata_Seurat, feature.list = feature.list,
                                               module.key = "MODULE",
                                               method = module.method,
                                               assay = module.assay,
                                               out.fig.dir = file.path(meld_fig_path,"module_scores"),
                                               plot = TRUE,
                                               plot.col.limits = NULL,
                                               plot.reduction = module.plot.reduction)
    }
  },error = function(e){log_error("Error in seurat_get_summary_scores: ",message(e))
    return(adata_Seurat)
    })
  #-----------------------------------------------------------------------------|
  # Visualize Meld Density --------
  #-----------------------------------------------------------------------------|
  tryCatch({
    log_info("Visualize Meld Density:")
    
    viz_meld_umap = function(feature, sobj, assay, save = TRUE, save.dir, device = "pdf"){
      if (!dir.exists(save.dir)){dir.create(save.dir, recursive = TRUE)}
      tmp  = DefaultAssay(sobj)
      DefaultAssay(sobj) = assay
      
      #p =FeaturePlot_scCustom(sobj, features = feature,reduction = "umap")
      p = FeaturePlot(sobj, features = feature,reduction = "umap")+
        scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu")), limits = c(0,1))
      DefaultAssay(sobj) = tmp
      
      if(save){
        ggsave(p , file = file.path(save.dir,sprintf("%s_%s.%s",assay, stringr::str_replace(feature," ",""),device)), width = 4, height = 3.5)
      }
      return(p)
    }
    
    
    
    
    viz_meld_umap_schex = function(feature, sobj, assay, save = TRUE, save.dir, device = "pdf"){
      if (!dir.exists(save.dir)){dir.create(save.dir, recursive = TRUE)}
      tmp  = DefaultAssay(sobj)
      DefaultAssay(sobj) = assay

      sobj_small <- make_hexbin(sobj, nbins = 30, dimension_reduction = "UMAP")
      
      p <- schex::plot_hexbin_feature(sobj_small, type="data", feature=feature, 
                                      action="mean", xlab="UMAP1", ylab="UMAP2", 
                                      title=paste0("Mean of ", feature))+
        scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu")), limits = c(0,1))
      
      DefaultAssay(sobj) = tmp
      
      if(save){
        ggsave(p, file = file.path(save.dir,sprintf("schex_%s_%s.%s",assay, stringr::str_replace(feature," ",""),device)), width = 4, height = 3.5)
      }
      return(p)
    }
    
    
    ## Meld Treatment-associated Density Normalized to Ctrl within Batch -----------
    
    for (assay_oi in meld.assay){
      log_info(assay_oi)
      # single-cell
      p_list1 = lapply(rownames(adata_Seurat@assays[[assay_oi]]),
                      viz_meld_umap, adata_Seurat, assay = assay_oi,
                      save = TRUE, save.dir = file.path(meld_fig_path,"meld_umap"),device = "pdf")
      fig1 = ggarrange(plotlist = p_list1, nrow = 4, ncol = 4)
      ggsave(fig1 , file = file.path(meld_fig_path,"meld_umap",sprintf("%s.png", assay_oi)), width = 15, height = 10)
      
      if(schex){
        # schex
        p_list2 = lapply(rownames(adata_Seurat@assays[[assay_oi]]),
                         viz_meld_umap_schex, adata_Seurat, assay = assay_oi,
                         save = TRUE, save.dir = file.path(meld_fig_path,"meld_umap"),device = "pdf")
        fig2 = ggarrange(plotlist = p_list2, nrow = 4, ncol = 4)
        ggsave(fig2 , file = file.path(meld_fig_path,"meld_umap",sprintf("schex_%s.png", assay_oi)), width = 15, height = 10)
      }

    }
    
    },error = function(e){
      log_error("Error in viz_meld_umap",message(e))
      return(adata_Seurat)})
  
  
  
  
  
  
  
  return(adata_Seurat)
  
  
}
