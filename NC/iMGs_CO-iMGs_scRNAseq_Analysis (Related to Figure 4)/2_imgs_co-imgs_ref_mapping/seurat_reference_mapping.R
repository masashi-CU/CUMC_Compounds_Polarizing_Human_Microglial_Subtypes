#' Seurat Reference Mapping Pipeline
#' Source: https://satijalab.org/seurat/articles/integration_mapping
#' Natacha Comandante-Lou (nc3018 at columbia dot edu)
#' @import Seurat
#' @import logger
#' @import ggplot2
#' @import scCustomize
#' @export

#############################################################################|
# Preprocess Data ----------------------
#############################################################################|
my_cols = c(`1`="#7FCCDD",`2`="#4559A4", `3`="#8CCF96", `4`="#FFAC53",`5`="#F4EA3B",
            `6`="#7E66F6",`7`="#8056A0",`8`="#FF6D00", `9`="#E0312B", `10`="#2F6539",`11`='#FF00F6', `12`="#000000",
            `13` = "#3283FEFF",`14` = "#B00068FF",`15` = "#1CFFCEFF",`16` = "#90AD1CFF",`17` = "#2ED9FFFF",`18`= "#85660DFF")
prep_transfer = function(sobj, nhvg = 3000, assay = "SCT", resolution = 0.7){
  
  log_info("#-----------------------------------------------------------------#")
  log_info("#                 PREPARE DATA FOR TRANSFER                       #")
  log_info("#-----------------------------------------------------------------#")
  log_info("* Assay: ", assay)
  log_info("* Number of top variable genes: ", nhvg)
  log_info("* Clustering resolution: ", resolution)
  
  DefaultAssay(sobj) = assay
  if (assay == "RNA"){
    sobj = NormalizeData(sobj, verbose = FALSE)%>%
      FindVariableFeatures(selection.method = "vst", nfeatures = nhvg,verbose = FALSE)%>%
      ScaleData()%>%RunPCA(assay = assay)
    
    sobj = RunUMAP(sobj, assay = assay, reduction = "pca", dims = 1:sobj@misc[["optimalPC"]], verbose = FALSE,return.model = TRUE)%>%
      FindNeighbors(reduction = "pca", dims = 1:sobj@misc[["optimalPC"]],verbose = FALSE) %>%
      FindClusters(resolution = resolution, verbose = FALSE)
  }
  if (assay == "SCT"){
    sobj = SCTransform(sobj, vst.flavor = "v2", verbose = FALSE) %>%
      RunPCA(assay = assay)
    sobj = RunUMAP(sobj,assay = assay, reduction = "pca", dims = 1:sobj@misc[["optimalPC"]], verbose = FALSE,return.model = TRUE)%>%
      FindNeighbors(reduction = "pca", dims = 1:sobj@misc[["optimalPC"]],verbose = FALSE) %>%
      FindClusters(resolution = resolution, verbose = FALSE)
    
    
    # features <- SelectIntegrationFeatures(object.list = sobj, nfeatures = nhvg)
    # sobj <- PrepSCTIntegration(object.list = sobj, anchor.features = features)
    # 
  }
  
  if (assay == "originalexp"){
    # batched corrected count in scVI was already log normalized 
    # (see https://ccc-protocols.readthedocs.io/en/latest/notebooks/ccc_R/S1_Batch_Correction.html)
    sobj = FindVariableFeatures(sobj,selection.method = "vst", nfeatures =nfeatures,verbose = FALSE)%>%
      ScaleData()%>%RunPCA(assay = assay)
    
    nPCs = dim(sobj@reductions[["X_scVI"]])[2]
    sobj = RunUMAP( sobj, dims = 1:nPCs, reduction = "X_scVI", n.components = 2,assay = "originalexp",slot = "data",return.model = TRUE)%>%
      FindNeighbors(dims = 1:nPCs, reduction = "X_scVI",verbose = FALSE)%>%
      FindClusters(resolution = resolution, verbose = FALSE)
    
  }  
  return(sobj)
}

#############################################################################|
# Classification of query-----------
#############################################################################|

map_query = function(ref,query,
                     npcs = NULL,
                     normalization.method = "SCT",
                     pred.var,
                     reference.reduction = "pca",
                     reduction.model = "umap",
                     plot = TRUE,
                     plot.reduction = "umap",
                     plot.split.by,
                     output.dir = "figures"
                     ){

  log_info("#-----------------------------------------------------------------#")
  log_info("#                        MAP QUERY                                #")
  log_info("#-----------------------------------------------------------------#")
  log_info("* Normalization Method: ", normalization.method)
  log_info("* Reference Reduction: ", reference.reduction)
  log_info("* Predicting variable: ", pred.var)
  
  if (is.null(npcs)){nPCs = ref@misc[["optimalPC"]]}else{nPCs = npcs}
  
  reference.reduction = "pca"
  
  ref.levels = sort(unique(ref[[pred.var]][,1]))
  ref[[pred.var]][,1] = factor(ref[[pred.var]][,1], levels = ref.levels)
  Idents(ref) = pred.var
  
  anchors <- FindTransferAnchors(
    reference = ref,
    query = query,
    normalization.method =  normalization.method,
    reference.reduction = reference.reduction,
    dims = 1:nPCs
  )
  
  query<- MapQuery(
    anchorset = anchors,
    query =query,
    reference = ref,
    refdata = list(cell.label = pred.var),
    reference.reduction = reference.reduction, 
    reduction.model = reduction.model
  )
  
  query$predicted.cell.label =factor(query$predicted.cell.label, levels = ref.levels)
  
  
  
  if (plot){
    p1 = DimPlot_scCustom(ref, reduction = plot.reduction, split.by = plot.split.by, label = FALSE, label.size = 3,colors_use = my_cols,
                          repel = TRUE) + NoLegend()+ggtitle("Reference annotations")+xlim(c(-10,11))+ylim(c(-10,6))
    
    
    p2 = DimPlot_scCustom(query, reduction = plot.reduction, split.by =  plot.split.by,group.by = "predicted.cell.label", label = FALSE,label.size = 3,
                          colors_use = my_cols,repel = TRUE) + NoLegend()+ ggtitle("Query Predicted labels")+xlim(c(-10,11))+ylim(c(-10,6))
    
    p = ggpubr::ggarrange(p1,p2, nrow = 2)
    ggsave(p, filename = file.path(output.dir,sprintf("dim_reduc/%s_ref_vs_query_%s.png", plot.reduction, normalization.method)), width = 16, height = 15)
  }
  
  return(query)
  
}
 
#############################################################################|
# Run Pipeline -----------
#############################################################################|
seurat_reference_mapping_pipeline = function(ref,query,
                                             ref.assay = "SCT",
                                             query.assay = "SCT",
                                             nfeatures = 3000,
                                             pred.var,
                                             resolution = 0.7,
                                             npcs = NULL,
                                             normalization.method = "SCT",
                                             reference.reduction = "pca",
                                             reduction.model = "umap",
                                             plot = TRUE,
                                             plot.reduction = "umap",
                                             plot.split.by,
                                             output.dir = "figures",
                                             save=TRUE){
  
  if(!dir.exists(output.dir)){dir.create(output.dir, recursive = TRUE)}
   # Log -------
  log_layout(layout_simple)
  log_file = file.path(output.dir , "seurat_ref_mapping.log")
  if(file.exists(log_file)){
    log_warn("Overwriting previous log: ", log_file)
    file.remove(log_file)
  }
  
  log_appender(appender_tee(file = log_file))
  # Run ------
  
 
  start = Sys.time()
  
  ref = prep_transfer(ref, nhvg = nfeatures,assay = ref.assay)
  query = prep_transfer(query, nhvg = nfeatures,assay = query.assay)
  
  query = map_query(ref,query,
                    npcs = npcs,
                    normalization.method =  normalization.method,
                    pred.var = pred.var,
                    reference.reduction = reference.reduction,
                    reduction.model = reduction.model,
                    plot = plot,
                    plot.reduction = plot.reduction,
                    plot.split.by = plot.split.by,
                    output.dir = output.dir)
  
  
  # merge reference and query
  log_info("Merging reference and query...")
  ref$id <- 'reference'
  nPCs = ref@misc[["optimalPC"]]
  query$id <- 'query'
  refquery <- merge(ref, query)
  refquery[["pca"]] <- merge(ref[["pca"]], query[["pca"]])
  refquery <- RunUMAP(refquery, reduction = 'pca', dims = 1:nPCs)
  if (plot){
    p = DimPlot_scCustom(refquery, group.by = 'treatment', split.by = "id", shuffle = TRUE,pt.size = 0.5)
    ggsave(p, filename = file.path(output.dir,sprintf("dim_reduc/%s_refquery_merged_%s.png", plot.reduction,normalization.method)), width = 16, height = 7)
  }
  
  end = Sys.time()
  log_success("Elapsed Time: ", round(end-start, 2) ," ", units(end-start))
  
  output = list()
  output$ref = ref
  output$query = query
  output$refquery = refquery
  if(save){
    log_info("Saving...")
    
    saveRDS(refquery,file = file.path(output.dir, "refquery.rds"))
  }

  
  
  log_appender(appender_console)
  return(output)
}


