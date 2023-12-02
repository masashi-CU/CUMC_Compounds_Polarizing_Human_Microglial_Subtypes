###############################
#         jft2113             #
#         5.07.23             #
#         R4.1                #
############################### 

#load in packages for analysis#
library(Matrix)
library(Seurat)
library(demuxmix)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(viridis)
library(topGO)
library(org.Hs.eg.db)
library(grid)
library(ggtree)
library(tidyverse)
library(DropletUtils)
library(pheatmap)
library(reshape2)
library(forcats)
library(writexl)
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(patchwork)
library(data.table)
library(magrittr)
library(qvalue)
library(rstatix)
library(ggpubr)
library(openxlsx)

####set up workspace and load in utilities/wrapper functions####
setwd("Microglia_project/Analysis_2020/Full_uglia_analysis_2020/")
#load demuxmix and helper functions
source("full_functions.R")
#now, create a function for chooseR to run on the dataset for clustering optimization#
source("ChooseR/helper_functions.R")
source("ChooseR/pipeline.R")
chooseR <- function(obj, npcs = 20, resolutions = c(0.3,0.6,0.9,1.2,1.5,2,4,6), assay = "SCT", reduction = "pca", results_path = "results/"){
  # Run pipeline
  for (res in resolutions) {
    message(paste0("Clustering ", res, "..."))
    message("\tFinding ground truth...")
    
    # "Truths" will be stored at glue::glue("{reduction}.{assay}_res.{res}")
    obj <- find_clusters(
      obj,
      reduction = reduction,
      assay = assay,
      resolution = res,
      npcs = npcs
    )
    clusters <- obj[[glue::glue("{reduction}.{assay}_res.{res}")]]
    
    # Now perform iterative, sub-sampled clusters
    results <- multiple_cluster(
      obj,
      n = 50,
      size = 0.8,
      npcs = npcs,
      res = res,
      reduction = reduction,
      assay = assay
    )
    
    # Now calculate the co-clustering frequencies
    message(paste0("Tallying ", res, "..."))
    # This is the more time efficient vectorisation
    # However, it exhausts vector memory for (nearly) all datasets
    # matches <- purrr::map(columns, find_matches, df = results)
    # matches <- purrr::reduce(matches, `+`)
    columns <- colnames(dplyr::select(results, -cell))
    mtchs <- matrix(0, nrow = dim(results)[1], ncol = dim(results)[1])
    i <- 1 # Counter
    for (col in columns) {
      message(paste0("\tRound ", i, "..."))
      mtchs <- Reduce("+", list(
        mtchs,
        find_matches(col, df = results)
      ))
      i <- i + 1
    }
    
    message(paste0("Scoring ", res, "..."))
    mtchs <- dplyr::mutate_all(
      dplyr::as_tibble(mtchs),
      function(x) dplyr::if_else(Re(x) > 0, percent_match(x), 0)
    )
    
    # Now calculate silhouette scores
    message(paste0("Silhouette ", res, "..."))
    sil <- cluster::silhouette(
      x = as.numeric(as.character(unlist(clusters))),
      dmatrix = (1 - as.matrix(mtchs))
    )
    saveRDS(sil, paste0(results_path, "silhouette_", res, ".rds"))
    
    # Finally, calculate grouped metrics
    message(paste0("Grouping ", res, "..."))
    grp <- group_scores(mtchs, unlist(clusters))
    saveRDS(grp, paste0(results_path, "frequency_grouped_", res, ".rds"))
    sil <- group_sil(sil, res)
    saveRDS(sil, paste0(results_path, "silhouette_grouped_", res, ".rds"))
  }
  return(obj)
}
chooseR_viz <- function(obj, npcs = 20, resolutions = c(0.3,0.6,0.9,1.2,1.5,2,4,6), assay = "SCT", reduction = "pca", results_path = "results/"){
  scores <- purrr::map(
    paste0(results_path, "silhouette_grouped_", resolutions, ".rds"),
    readRDS
  )
  scores <- dplyr::bind_rows(scores) %>%
    dplyr::group_by(res) %>%
    dplyr::mutate("n_clusters" = dplyr::n()) %>%
    dplyr::ungroup()
  meds <- scores %>%
    dplyr::group_by(res) %>%
    dplyr::summarise(
      "boot" = list(boot_median(avg_sil)),
      "n_clusters" = mean(n_clusters)
    ) %>%
    tidyr::unnest_wider(boot)
  
  writexl::write_xlsx(meds, paste0(results_path, "median_ci.xlsx"))
  
  # Find thresholds
  threshold <- max(meds$low_med)
  choice <- as.character(
    meds %>%
      dplyr::filter(med >= threshold) %>%
      dplyr::arrange(n_clusters) %>%
      tail(n = 1) %>%
      dplyr::pull(res)
  )
  
  # And plot!
  ggplot(meds, aes(factor(res), med)) +
    geom_crossbar(
      aes(ymin = low_med, ymax = high_med),
      fill = "grey",
      size = 0.25
    ) +
    geom_hline(aes(yintercept = threshold), colour = "blue") +
    geom_vline(aes(xintercept = choice), colour = "red") +
    geom_jitter(
      data = scores,
      aes(factor(res), avg_sil),
      size = 0.35,
      width = 0.15
    ) +
    scale_x_discrete("Resolution") +
    scale_y_continuous(
      "Silhouette Score",
      expand = c(0, 0),
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.25),
      oob = scales::squish
    ) +
    cowplot::theme_minimal_hgrid() +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
    )
  
  ggsave(
    filename = paste0(results_path, "silhouette_distribution_plot.png"),
    dpi = 300,
    height = 3.5,
    width = 3.5,
    units = "in"
  )
  
  # Finally, a dot plot of silhouette scores to help identify less robust clusters
  # The initial pipe is to order the clusters by silhouette score
  scores %>%
    dplyr::filter(res == choice) %>%
    dplyr::arrange(dplyr::desc(avg_sil)) %>%
    dplyr::mutate_at("cluster", ordered, levels = .$cluster) %>%
    ggplot(aes(factor(cluster), avg_sil)) +
    geom_point() +
    scale_x_discrete("Cluster") +
    scale_y_continuous(
      "Silhouette Score",
      expand = c(0, 0),
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.25),
      oob = scales::squish
    ) +
    cowplot::theme_minimal_grid() +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
    )
  
  ggsave(
    filename = paste0(results_path, "silhouette_point_plot_", choice, ".png"),
    dpi = 300,
    height = 3.5,
    width = 3.5,
    units = "in"
  )
  return(choice)
}
chooseR_viz_extended <- function(obj, assay = "SCT", reduction = "pca", results_path = "results/", choice = 0.5){
  # First is a cluster average co-clustering heatmap
  # Read the data
  grp <- readRDS(paste0(results_path, "frequency_grouped_", choice, ".rds"))
  
  # As the data is symmetrical, we do not need the upper triangle
  grp <- grp %>%
    pivot_wider(names_from = "cell_2", values_from = "avg_percent") %>%
    dplyr::select(str_sort(colnames(.), numeric = T)) %>%
    column_to_rownames("cell_1")
  grp[lower.tri(grp)] <- NA
  grp <- grp %>%
    as_tibble(rownames = "cell_1") %>%
    pivot_longer(-cell_1, names_to = "cell_2", values_to = "avg_percent") %>%
    mutate_at("cell_2", ordered, levels = unique(.$cell_1)) %>%
    mutate_at("cell_1", ordered, levels = unique(.$cell_1))
  
  # And plot!
  plot <- ggplot(grp, aes(factor(cell_1), cell_2, fill = avg_percent)) +
    geom_tile() +
    scale_x_discrete("Cluster", expand = c(0, 0)) +
    scale_y_discrete(
      "Cluster",
      limits = rev(levels(grp$cell_2)),
      expand = c(0, 0)
    ) +
    scale_fill_distiller(
      " ",
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      palette = "RdYlBu",
      na.value = "white"
    ) +
    coord_fixed() +
    theme(
      axis.ticks = element_line(colour = "black"),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.position = c(0.9, 0.9)
    ) +
    guides(fill = guide_colorbar(barheight = 3, barwidth = 1))
  
  ggsave(
    plot = plot,
    filename = paste0(results_path, "coclustering_heatmap_", choice, ".png"),
    dpi = 300,
    height = 3.5,
    width = 3.5,
    units = "in"
  )
  
  # Let's add the silhouette scores to the Seurat object!
  sil_scores <- readRDS(paste0(results_path, "silhouette_", choice, ".rds"))
  sil_scores <- as.data.frame(sil_scores[, 3], row.names = Seurat::Cells(obj))
  colnames(sil_scores) <- c("sil_score")
  obj <- AddMetaData(obj, metadata = sil_scores)
  
  # Let's visualise the selected cluster
  # If your data has known  clusters, you could also visualise those!
  # Remember, truths are in "glue({reduction}.{assay}_res.{choice})"
  # Seurat Changes color scheme if you order your data, so we provide
  # the following helper function to restore defaults
  gg_color <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    colours <- hcl(h = hues, c = 100, l = 65)[1:n]
    return(colours)
  }
  
  plot <- DimPlot(
    obj,
    reduction = "umap",
    group.by = glue::glue("{reduction}.{assay}_res.{choice}"),
    pt.size = 0.5,
    # cols = gg_color(6) # Only necessary if you have ordered your clusters
  )
  
  ggsave(
    plot = plot,
    filename = paste0(results_path, choice, "_cluster_umap.png"),
    dpi = 300,
    height = 5,
    width = 5,
    units = "in"
  )
  
  # We also find it useful to visualise the silhouette scores on the UMAP!
  plot <- FeaturePlot(
    obj,
    "sil_score",
    reduction = "umap",
    pt.size = 0.5,
    min.cutoff = -1,
    max.cutoff = 1
  ) +
    scale_colour_distiller(
      palette = "RdYlBu",
      labels = c(-1, 0, 1),
      breaks = c(-1, 0, 1),
      limits = c(-1, 1)
    )
  
  ggsave(
    plot = plot,
    filename = paste0(results_path, choice, "_silhouette_umap.png"),
    dpi = 300,
    height = 5,
    width = 5,
    units = "in"
  )
  ########USABLE ONLY IF THERE ARE EXPERT ANNOTATIONS AVAILABLE FOR THE DATASET##############
  # # This dataset does have annotated cell types, so we can create a heatmap
  # # to see how well they line up. Remember, "truths" will always be in
  # # glue::glue("{reduction}.{assay}_res.{choice}"), but annotated cell types
  # # might be in different locations. Be sure to change the names as needed.
  # ids <- as_tibble(
  #   obj[[c(glue::glue("{reduction}.{assay}_res.{choice}"), "CellType")]]
  # ) %>%
  #   mutate_at(
  #     c(glue::glue("{reduction}.{assay}_res.{choice}"), "CellType"),
  #     as.factor
  #   ) %>%
  #   group_by(pca.SCT_res.1.6, CellType) %>%
  #   summarise("count" = n()) %>%
  #   ungroup() %>%
  #   complete(pca.SCT_res.1.6, CellType, fill = list("count" = 0)) %>%
  #   # Count suggested clusters
  #   group_by(pca.SCT_res.1.6) %>%
  #   mutate("n_suggested" = sum(count)) %>%
  #   ungroup() %>%
  #   # Count known clusters
  #   group_by(CellType) %>%
  #   mutate("n_known" = sum(count)) %>%
  #   ungroup() %>%
  #   # Calculate statistics
  #   mutate("dice" = (2 * count) / (n_suggested + n_known)) %>%
  #   mutate(
  #     "pca.SCT_res.1.6" = fct_reorder(pca.SCT_res.1.6, dice, max, .desc = TRUE),
  #     "CellType" = fct_reorder(CellType, dice, max, .desc = TRUE)
  #   )
  # 
  # plot <- ggplot(ids, aes(pca.SCT_res.1.6, CellType, fill = dice)) +
  #   geom_tile() +
  #   scale_x_discrete("Suggested Clusters", expand = c(0, 0)) +
  #   scale_y_discrete(
  #     "Known Clusters",
  #     limits = rev(levels(ids$CellType)),
  #     expand = c(0, 0)
  #   ) +
  #   coord_fixed() +
  #   scale_fill_distiller(
  #     NULL,
  #     limits = c(0, 1),
  #     breaks = c(0, 0.5, 1),
  #     palette = "RdYlBu"
  #   ) +
  #   theme(
  #     axis.ticks = element_line(colour = "black"),
  #     axis.text = element_text(size = 6),
  #     axis.title = element_text(size = 8),
  #     legend.text = element_text(size = 7),
  #   ) +
  #   guides(fill = guide_colorbar(barheight = 3, barwidth = 0.5))
  # 
  # ggsave(
  #   plot = plot,
  #   filename = paste0(results_path, "cluster_relation_heatmap.png"),
  #   dpi = 300,
  #   height = 3.5,
  #   width = 3.5,
  #   units = "in"
  # )
}
GO_topGO <- function(total_terms = 100, top_x = 30, ont = "BP", bg_genes, diffex_genes, condition_name, color = "orange", outdir, gene_ID = "symbol"){
  DEgenes = factor(as.integer(bg_genes %in% diffex_genes))
  names(DEgenes) = bg_genes
  GOdata <- new("topGOdata",
                ontology = ont, # use biological process ontology
                allGenes = DEgenes,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.org, mapping = "org.Hs.eg.db", ID = gene_ID)
  test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(GOdata, test.stat)
  allRes = GenTable(GOdata, Fisher = resultFisher, topNodes = total_terms, numChar = 100, orderBy = "Fisher")
  allRes$Fisher = as.numeric(allRes$Fisher)
  allRes$Fisher_corrected = p.adjust(allRes$Fisher, "BH")
  allRes_trim = allRes[1:top_x,]
  
  ggplot(allRes_trim, aes(x = seq(1:top_x), y = abs(log(Fisher_corrected)))) +
    geom_segment(aes(x=seq(1:top_x), xend=seq(1:top_x), y=0, yend=-log(Fisher_corrected)), color="grey15", size = 2.5) +
    geom_point(color=color, size=15, shape = 21, fill = alpha(color, 0.97), alpha = 0.9, stroke = 2) + 
    geom_hline(yintercept = -log(0.05,10), linetype="dashed", color = "black", size=1)+
    theme_minimal() +
    coord_flip() +
    theme(axis.text.y = element_text(face="bold", color="black", size = 20), axis.text.x = element_text(color="black", size = 15),
          axis.title.y = element_blank(), axis.title.x = element_text(color="black", size = 25),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(colour = "black", size = 1, fill = NA)
    ) +
    scale_x_continuous(breaks = c(1:top_x), labels = c(allRes_trim$Term), expand = expansion(add = c(0.85, 0.85)))+ 
    scale_y_continuous(expand = expansion(add = c(0.03,1))) + 
    xlab("") + ylab("-log(FDR-corrected p-value)") + ggtitle(paste(condition_name, ": top ", top_x, " ", ont, sep = ""))
  ggsave(paste(outdir, condition_name, "_top_", top_x, "_", ont, ".png", sep = ""), units = "mm", width = 450, height = 300, dpi = 300)
  return(allRes)
}

#markers from orignal dataset for comparison#
final_156 = c("MEF2A", "QKI", "ITPR2", "DENND3", "SRGAP2", "MAF", "ATM") ###IL6ST for cluster 5, PNISR, DDX17 in place of ATM would be fine too ###
final_2349 = c("CX3CR1", "FCGR1A", "TYROBP", "ITM2B", "AIF1", "GPX1","FCER1G", "FTL", "TPT1", "APOC1")
final_81011 = c("CXCR4", "SRGN", "CD74", "HLA-C", "CYBA", "SPP1", "LGALS1", "CD9") ###FOLR2 also works for C10, but is it a BAM marker...? CPVL, FPR1 are also good. GLDN, LIPA, APOC1 for C11###
final_12 = c("H2AFZ", "MKI67", "PCNA")
final = rev(unique(c(final_156, final_2349, final_81011, final_12)))

#requires a seurat object with a computed UMAP, a module_list with w/ 2 columns ("Gene", "Pathway"), final output destination#
#first one plots module score on the UMAP for the data#
plot_module_enrichment <- function(seurat, module_list, output_dir = "../Compound_work/HMC3 scRNAseq/Viz/Metabolism/") {
  colnames(module_list) <- c("Gene", "Pathway")
  #subset to those genes in the dataset#
  module_list = module_list[module_list$Gene %in% rownames(seurat),]
  #plot the various modules overlaid on the data#
  all_modules = unique(module_list$Pathway)
  for(module in all_modules){
    all_genes = list(module_list$Gene[module_list$Pathway == module])
    seurat = AddModuleScore(object = seurat, features = all_genes, name = gsub(pattern = " ", x = module, replacement =  "_"))
  }
  pdf(paste0(output_dir, "pathway_modules.pdf"), width = 11, height = 8.5)
  for(module in all_modules){
    col_nam <-
      gsub(pattern = "-", replacement = ".", x = module) %>%
      gsub(pattern = " ", replacement =  "_", .) %>%
      gsub(pattern = ",", replacement =  ".", .)
    print(FeaturePlot(seurat, features = paste0(col_nam, 1), pt.size = 1) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), 
                                                                                                    guide = guide_colorbar(label = TRUE,
                                                                                                                           draw.ulim = TRUE, 
                                                                                                                           draw.llim = TRUE,
                                                                                                                           frame.colour = "black", ticks = TRUE, 
                                                                                                                           label.position = "bottom",
                                                                                                                           barwidth = 18,
                                                                                                                           barheight = 1.3, 
                                                                                                                           direction = 'horizontal')) +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position = "bottom", text=element_text(size = 20)) + NoAxes() + ggtitle(module))
  }
  dev.off()
}
#second one plots module scores per cluster#
plot_cluster_module_enrichment <- function(seurat, module_list, allgenes, markers, markers_down, FDR_thresh = 0.01, output_dir){
  #set up the modules#
  colnames(module_list) <- c("Gene", "Pathway")
  all_modules = unique(module_list$Pathway)
  all_assoc = list()
  for(module in all_modules){
    all_genes = module_list$Gene[module_list$Pathway == module]
    all_assoc[[gsub(pattern = " ", x = module, replacement =  "_")]] = all_genes
  }
  #we assemble our data frame for testing#
  markers <- subset(markers, num_down_types>=3)
  markers$cluster <- paste("MG_",markers$up_type,sep="")
  df<-data.frame(matrix(NA,nrow=length(unique(markers$cluster)),ncol=6))
  rownames(df)<-unique(markers$cluster)
  colnames(df)<-c("n.genes","N.genes","n.cellTypeMarker","N.all","p.genes","common_genes")
  
  markers_down <- subset(markers_down, num_up_types >= 3)
  markers_down$cluster <- paste("MG_",markers_down$down_type,sep="")
  df_down<-data.frame(matrix(NA,nrow=length(unique(markers_down$cluster)),ncol=6))
  rownames(df_down)<-unique(markers_down$cluster)
  colnames(df_down)<-c("n.genes","N.genes","n.cellTypeMarker","N.all","p.genes","common_genes")
  #subset everything for testing#
  geneInfo1<-subset(markers,select=c(gene,cluster))
  geneInfo1<-geneInfo1[which(!is.na(geneInfo1$gene)),]
  geneInfo1<-geneInfo1[which(geneInfo1$gene!=""),]
  geneInfo1<-geneInfo1[!duplicated(geneInfo1),]
  
  geneInfo2<-subset(markers_down,select=c(gene,cluster))
  geneInfo2<-geneInfo2[which(!is.na(geneInfo2$gene)),]
  geneInfo2<-geneInfo2[which(geneInfo2$gene!=""),]
  geneInfo2<-geneInfo2[!duplicated(geneInfo2),]
  
  ##pull the list of interest for upregulated genes##
  overall = data.frame(matrix(ncol = length(unique(markers$up_type)), nrow = length(all_assoc)))
  rownames(overall) = sapply(strsplit(names(all_assoc), "[.]"), "[", 1)
  colnames(overall) = c(1:length(unique(markers$up_type)))
  overall_number = data.frame(matrix(ncol = length(unique(markers$up_type)), nrow = length(all_assoc)))
  rownames(overall_number) = rownames(overall); colnames(overall_number) = colnames(overall)
  
  ##pull the list of interest for downregulated genes##
  overall_down = data.frame(matrix(ncol = length(unique(markers$up_type)), nrow = length(all_assoc)))
  rownames(overall_down) = sapply(strsplit(names(all_assoc), "[.]"), "[", 1)
  colnames(overall_down) = c(1:length(unique(markers$up_type)))
  
  for(j in 1:length(all_assoc)){
    disease = names(all_assoc)[j]
    d_genes = data.table(all_assoc[[j]])
    colnames(d_genes) = "gene"
    d_genes = unique(d_genes)
    
    for (i in 1:length(unique(geneInfo1$cluster))){
      print(i)
      print(unique(geneInfo1$cluster)[i])
      
      marker.list <- geneInfo1[which(geneInfo1$cluster==unique(geneInfo1$cluster)[i]),]
      
      a <- dim(merge(d_genes,marker.list,by="gene"))[1]
      b <- dim(merge(d_genes,allgenes,by="gene"))[1]
      c <- dim(marker.list)[1]
      d <- dim(allgenes)[1]-dim(marker.list)[1]
      
      gene.list <- merge(d_genes,marker.list,by="gene")
      gene.list <-paste0( gene.list$gene, collapse=",")
      
      df[i,1]<-a
      df[i,2]<-b
      df[i,3]<-c
      df[i,4]<-c+d
      df[i,5]<-fisher.test(matrix(c(a, b-a, c-a, d-b+a), 2, 2), alternative='greater')$p.value
      
      df[i,6]<-gene.list
    }
    ###fdr correction of p values###
    df$fdr.p.genes<-qvalue(df$p.genes)$qvalue
    df$Sample<-"Microglia"
    df$cluster<-rownames(df)
    write.csv(df,paste(output_dir, "up_", disease, ".csv",sep =""),row.names = F)
    #load these in, remove, 13, and add to aggregate frames#
    enrichment <- read.csv(paste(output_dir, "up_", disease, ".csv",sep =""))
    overall[disease,] = t(enrichment$fdr.p.genes)
    overall_number[disease,] = t(enrichment$n.genes)
  }
  
  ##now for downreg list##
  for(j in 1:length(all_assoc)){
    disease = names(all_assoc)[j]
    d_genes = data.table(all_assoc[[j]])
    colnames(d_genes) = "gene"
    d_genes = unique(d_genes)
    
    for (i in 1:length(unique(geneInfo2$cluster))){
      print(i)
      print(unique(geneInfo2$cluster)[i])
      
      marker.list <- geneInfo2[which(geneInfo2$cluster==unique(geneInfo2$cluster)[i]),]
      
      a <- dim(merge(d_genes,marker.list,by="gene"))[1]
      b <- dim(merge(d_genes,allgenes,by="gene"))[1]
      c <- dim(marker.list)[1]
      d <- dim(allgenes)[1]-dim(marker.list)[1]
      
      gene.list <- merge(d_genes,marker.list,by="gene")
      gene.list <-paste0( gene.list$gene, collapse=",")
      
      df_down[i,1]<-a
      df_down[i,2]<-b
      df_down[i,3]<-c
      df_down[i,4]<-c+d
      df_down[i,5]<-fisher.test(matrix(c(a, b-a, c-a, d-b+a), 2, 2), alternative='greater')$p.value
      
      df_down[i,6]<-gene.list
    }
    ###fdr correction of p values###
    df_down$fdr.p.genes<-qvalue(df_down$p.genes)$qvalue
    df_down$Sample<-"Microglia"
    df_down$cluster<-rownames(df_down)
    write.csv(df_down,paste(output_dir, "down_", disease, ".csv",sep =""),row.names = F)
    # load these in, remove, 13, and add to aggregate frames#
    enrichment <- read.csv(paste(output_dir, "down_", disease, ".csv",sep =""))
    overall_down[disease,] = t(enrichment$fdr.p.genes)
  }
  
  ###log-transform the matrices for convenience###
  overall = -log(overall, base = 10)
  overall[overall < -log(FDR_thresh,10)] = 0
  overall = trunc(overall*100)/100
  
  range <- max(overall)
  overall[overall == 0] = NA
  #subset out those rows that have something in them#
  overall = overall[apply(overall,MARGIN = 1, FUN = function(x){sum(is.na(x)) < 12}),]
  
  overall_2 = overall 
  overall_2[is.na(overall_2)] = ""
  
  ##heatmap to visualize + association##
  png(paste0(output_dir, "FDR0.01_metabo_up_pheatmap.png"), units = "mm", width = 400, height = 600, res = 300)
  pheatmap(overall,cluster_cols = F, scale = "none", color=colorRampPalette(c("white","red"))(90), cluster_rows = F, angle_col = 0, border_color = "black", cellwidth = 50, cellheight = 50, 
           display_numbers = overall_2, number_format = "%s", number_color = "black", breaks = seq(0, range*1.1, length.out = 100), fontsize = 20, fontsize_number = 15, na_col = "white")
  dev.off()
  
  overall_down = -log(overall_down, base = 10)
  overall_down[overall_down < -log(FDR_thresh,10)] = 0
  overall_down = trunc(overall_down*100)/100
  
  range <- max(overall_down)
  overall_down[overall_down == 0] = NA
  overall_down = overall_down[apply(overall_down,MARGIN = 1, FUN = function(x){sum(is.na(x)) < 12}),]
  
  overall_2 = overall_down 
  overall_2[is.na(overall_2)] = ""
  ##heatmap to visualize - association##
  png(paste0(output_dir, "FDR0.01_metabo_down_pheatmap.png"), units = "mm", width = 400, height = 600, res = 300)
  pheatmap(overall_down,cluster_cols = F, scale = "none", color=colorRampPalette(c("white","blue"))(90), cluster_rows = F, angle_col = 0, border_color = "black", cellwidth = 50, cellheight = 50, 
           display_numbers = overall_2, number_format = "%s", number_color = "white", breaks = seq(0, range*1.1, length.out = 100), fontsize = 20, fontsize_number = 15, na_col = "white")
  dev.off()
}

####retrieve location of scRNA-seq and read in the data####
data_directory <- "../Data/Drug_Tx/"
drug_files = list.files(data_directory)[1]

#for now, only one dataset to pull out#
seuratmat = list()
hashtag_frame = data.frame(Matrix(nrow = 0, ncol = 27))
#filtering parameters: min UMIs of 500, essentially no upper border, filter any cells w/o HTO counts since the sample was hashed, 
#filter out ribosomal and mitochondrial transcripts 
umimin = 500; umimax = 1e10; HTO_filter = 0; RP_filt = T; MT_filt = T
#read in the HTO mapping#
HTO_mapping = data.frame(read.csv("../Compound_work/HMC3 scRNAseq/HTO_mapping.csv"))
colnames(HTO_mapping) = c("Name", "Condition", paste0("Hashtag0",4:9), "Hashtag10")
#loop to read in different files and perform basic pre-processing#
for(file in drug_files){
  #read in initial file#
  tempmat=Read10X_h5(file=paste(data_directory, file, "/filtered_feature_bc_matrix.h5", sep =""), use.names = T,unique.features = T)
  hashmat=tempmat[[2]]
  rownames(hashmat) = sapply(strsplit(rownames(hashmat), "_"), "[", 1)
  cellmat=tempmat[[1]]
  joint.bcs <- intersect(colnames(cellmat), colnames(hashmat))
  cellmat = cellmat[,joint.bcs]; hashmat = hashmat[,joint.bcs]
  ###mito filtering: choose the more $stringent$ of median or IQR filtering###
  if(MT_filt == T){
    mtpct=colSums(cellmat[grep("^MT-",rownames(cellmat)),])/colSums(cellmat)*100
    med_cutoff = median(mtpct) + 2*mad(mtpct)
    IQR_cutoff = as.numeric(quantile(mtpct)[4] + IQR(mtpct)*1.5)
    mitofilter = min(med_cutoff, IQR_cutoff)
    mitofilter = max(mitofilter, 10)
    print(paste(file, " mito: med-", med_cutoff, ", IQR-", IQR_cutoff, sep = ""))
    cellmat=cellmat[,which(mtpct<mitofilter)]  ###remove mito-heavy cells
    hashmat=hashmat[,which(mtpct<mitofilter)]
  }
  ###filter genes unlikely to be biologically relevant out, then filter on UMI count
  if(RP_filt == T){
    cellmat=cellmat[grep("^MT-|^RP[0-9]|^BC[0-9]|^RPL|^RPS|^MTRNR|-PS",rownames(cellmat),invert=T),] ###variant to remove all RP, MT, and pseudogenes
  }else if (MT_filt == T){
    cellmat=cellmat[grep("^MT-|^BC[0-9]|^MTRNR|-PS",rownames(cellmat),invert=T),] ###variant to remove MT + PS
  } else{
    cellmat=cellmat[grep("^BC[0-9]|-PS",rownames(cellmat),invert=T),] ###variant to remove MT + PS
  }
  #what is the distribution of UMIs here#
  quantile(colSums(cellmat))
  hist(colSums(cellmat))
  
  keepcells=which(colSums(cellmat) > umimin & colSums(cellmat) < umimax & colSums(hashmat) >= HTO_filter)     ###Do we want this now or after hashtag deconvolution?###
  hashtag_frame = rbind(hashtag_frame, t(rowSums(hashmat[,keepcells])))
  rownames(hashtag_frame)[dim(hashtag_frame)[1]] = file
  
  #summed = rowsum(tempmat, row.names(tempmat))
  dat_seurat = CreateSeuratObject(counts=cellmat[,keepcells])
  hashmat = hashmat[,keepcells]
  dat_seurat[["HTO"]]=CreateAssayObject(counts=hashmat)
  tempcol=dat_seurat$orig.ident
  barcodes = names(tempcol)
  names(tempcol)=paste(file,"-",1:length(tempcol),sep="")
  dat_seurat = RenameCells(dat_seurat, new.names = names(tempcol))
  tempcol = rep(file, length(tempcol))
  dat_seurat$orig.ident = tempcol
  names(barcodes) = names(tempcol)
  dat_seurat = AddMetaData(dat_seurat, metadata = barcodes, col.name = "barcodes")
  #run demuxmix here#
  rna = colSums(dat_seurat@assays$RNA@counts)
  p.acpt = 0.9
  dmm <- demuxmix(hto = dat_seurat@assays$HTO@counts, 
                  rna=rna,
                  p.acpt = p.acpt^nrow(hashmat))
  #QC metrics#
  plotDmmHistogram(dmm)
  dmmOverlap(dmm)
  plotDmmPosteriorP(dmm)
  # plotDmmScatter(dmm)
  classLabels <- dmmClassify(dmm)
  colnames(classLabels)[1] = "HTO_assign"
  dat_seurat = AddMetaData(dat_seurat, metadata = classLabels)
  #reassign associated compound phenotype#
  type_col = dat_seurat$Type
  for(feature in rownames(hashmat)){
    type_col[dat_seurat$HTO_assign == feature] = HTO_mapping[HTO_mapping$Name == file,colnames(HTO_mapping) == feature]  
  }
  dat_seurat = AddMetaData(dat_seurat, metadata = type_col, col.name = "identity")
  seuratmat[[file]] = dat_seurat
}

####Now, analyze the data and identify doublets for removal####
#SCTransform workflow - could probably also use ChooseR, but this is not the final clustering#
dat_seurat = dat_seurat %>% SCTransform() %>%  RunPCA()
ElbowPlot(dat_seurat)
dat_seurat = dat_seurat %>% FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = c(0.3,0.5,0.7,1)) %>%
  RunUMAP(dims = 1:20)
DimPlot(dat_seurat, group.by = "identity", pt.size = 2) + ggtitle("") + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 40)) + NoAxes()
ggsave("../Compound_work/HMC3 scRNAseq/Viz/SCT_UMAP.png", units = "mm", width = 450, height = 350, dpi = 300)

dat_seurat = SetIdent(dat_seurat, value = "identity")
SCT_all_markers = FindAllMarkers(dat_seurat, test.use = "MAST",
                                 min.pct = 0.05, only.pos = F, assay = "RNA", slot = "data", logfc.threshold = 0.25)
SCT_all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> std_markers
DotPlot(dat_seurat, assay = "RNA", features = unique(std_markers$gene), dot.min = 0.1, cluster.idents = F, scale = T) + coord_flip() +
  scale_colour_viridis(option="magma") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.title = element_blank()) + scale_y_discrete(position = "right")
ggsave("../Compound_work/HMC3 scRNAseq/Viz/SCT_dotplot.png", units = "mm", width = 300, height = 250, dpi = 300)

#now, run DoubletFidner
dat_seurat = dat_seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(tmp)
nPC <- 15

seed = dim(dat_seurat)[2]
set.seed(seed)
sweep_list <- DoubletFinder::paramSweep_v3(dat_seurat, PCs = 1:nPC, sct = F)
sweep_stats <- DoubletFinder::summarizeSweep(sweep_list, GT = FALSE)

bcmvn = DoubletFinder::find.pK(sweep_stats) %>%
  mutate(pK = as.numeric(as.character(pK)))

max_pK = bcmvn %>% filter(MeanBC == max(bcmvn$MeanBC)) %>% pull(pK)
max_pK = as.numeric(as.character(max_pK))

expected_rate = 0.2
dat_seurat <- doubletFinder_v3(seu = dat_seurat,
                               PCs = 1:nPC,
                               pN = 0.25, # default
                               pK = max_pK,
                               nExp = expected_rate * dim(dat_seurat)[2]) # this is the value from the 10X multiplet rate table depending on v2 or v3 **multiplied by the number of pre QC cells**

DF.name = colnames(dat_seurat@meta.data)[grepl("DF.classification", colnames(dat_seurat@meta.data))]

#visualize comparative results of DF and demux
DimPlot(res, group.by = "Type", pt.size = 2) + ggtitle("") + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 40)) + NoAxes()
ggsave("../Compound_work/HMC3 scRNAseq/Viz/std_UMAP_type.png", units = "mm", width = 450, height = 350, dpi = 300)
DimPlot(res, group.by = DF.name, pt.size = 2) + ggtitle("") + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 40)) + NoAxes()
ggsave("../Compound_work/HMC3 scRNAseq/Viz/std_UMAP_DF_results.png", units = "mm", width = 450, height = 350, dpi = 300)
VlnPlot(res, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1) + ggtitle("DoubletFinder")
ggsave("../Compound_work/HMC3 scRNAseq/Viz/DF_violin.png", units = "mm", width = 300, height = 200, dpi = 300)
VlnPlot(res, features = "nFeature_RNA", group.by = "Type", pt.size = 0.1) + ggtitle("demuxmix")
ggsave("../Compound_work/HMC3 scRNAseq/Viz/demux_violin.png", units = "mm", width = 300, height = 200, dpi = 300)
VlnPlot(res, features = "nFeature_RNA", group.by = "identity", pt.size = 0.1) + ggtitle("demuxmix")
ggsave("../Compound_work/HMC3 scRNAseq/Viz/demux_ident_violin.png", units = "mm", width = 400, height = 200, dpi = 300)
#time to remove the doublets#
dat_seurat = dat_seurat[,res$`DF.classifications_0.25_5e-04_8537.75` == "Singlet" & res$Type == "singleton"]
dat_seurat = DietSeurat(dat_seurat, counts = T, scale.data = F, data = T, assays = c("RNA", "HTO", "SCT"),
                        dimreducs = "umap")
save(dat_seurat, file = "../Compound_work/HMC3 scRNAseq/HMC3_scSTIM_doubletrm.rda")

DimPlot(dat_seurat, group.by = "identity", pt.size = 2) + ggtitle("") + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 40)) + NoAxes()
ggsave("../Compound_work/HMC3 scRNAseq/Viz/SCT_UMAP_postdoublet.png", units = "mm", width = 450, height = 350, dpi = 300)

#####load the doublet-removed data back in and conduct a series of different analyses.
#first, load in the HMC3 dataset#
load("../Compound_work/HMC3 scRNAseq/HMC3_scSTIM_doubletrm.rda")

#now, let's run chooseR on this dataset#
dat_seurat_chooseR <- chooseR(obj = dat_seurat, npcs = 20, resolutions = c(0.3,0.6,0.9,1.2,1.5,2,4,6), assay = "SCT", reduction = "pca", results_path = "../Compound_work/HMC3 scRNAseq/ChooseR/")
rm(dat_seurat); gc()

# Save original data, with calculated clusterings#
# saveRDS(dat_seurat_chooseR, paste0("../Compound_work/HMC3 scRNAseq/ChooseR/", "chooseR_data.rds"))
dat_seurat_chooseR <- readRDS(paste0("../Compound_work/HMC3 scRNAseq/ChooseR/", "chooseR_data.rds"))
#create chooseR visualizations#
choice <- chooseR_viz(obj = dat_seurat_chooseR, npcs = 20, resolutions = c(0.3,0.6,0.9,1.2,1.5,2,4,6), assay = "SCT", reduction = "pca", results_path = "../Compound_work/HMC3 scRNAseq/ChooseR/")
#in this case, optimal choice was 0.6#
dat_seurat_chooseR <- SetIdent(object = dat_seurat_chooseR, value = paste0("pca.SCT_res.", choice))
new.ids = 1:length(unique(dat_seurat_chooseR$pca.SCT_res.0.6))
names(new.ids) <- 0:(length(unique(dat_seurat_chooseR$pca.SCT_res.0.6))-1)
dat_seurat_chooseR <- RenameIdents(dat_seurat_chooseR, new.ids)
DimPlot(dat_seurat_chooseR, pt.size = 2, label = T) + ggtitle("ChooseR optimal clustering v1") + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 40)) + NoAxes()
ggsave("../Compound_work/HMC3 scRNAseq/Viz/chooseR_optimalv1.png", units = "mm", width = 450, height = 350, dpi = 300)

#perform pairwise differential expression analysis#
out_stub = "../Compound_work/HMC3 scRNAseq/pairwise_diffex/"
all_pairwise_diffex(dat_seurat_chooseR, outdir = paste(out_stub, "RNA_MAST/RNA_MAST_", sep =""),
                    assay = "RNA", slot = "data", test = "MAST")
markers_up = pairwise_genes(DEG_file = paste(out_stub, "RNA_MAST/RNA_MAST_alldegenes.csv", sep =""), allpops = as.numeric(levels(Idents(dat_seurat_chooseR))), upreg = T)
write.csv(markers_up, paste(out_stub, "RNA_MAST/microglia_pairwise_up.csv", sep = ""))
markers_down = pairwise_genes(DEG_file = paste(out_stub, "RNA_MAST/RNA_MAST_alldegenes.csv", sep =""), allpops = as.numeric(levels(Idents(dat_seurat_chooseR))), upreg = F)
write.csv(markers_down, paste(out_stub, "RNA_MAST/microglia_pairwise_down.csv", sep = ""))

#examine association of original marker gene sets with these clusters#
DotPlot(dat_seurat_chooseR, assay = "RNA", features = final, dot.min = 0.1, cluster.idents = F, scale = T) + coord_flip() +
  scale_colour_viridis(option="magma", guide = guide_colorbar(frame.colour = "black", ticks = TRUE)) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.title = element_blank()) + scale_y_discrete(position = "right") +
  geom_vline(xintercept = 3.5, linetype = "dotted") + geom_vline(xintercept = 6.5, linetype = "dotted") + 
  geom_vline(xintercept = 9.5, linetype = "dotted") + geom_vline(xintercept = 11.5, linetype = "dotted") + 
  geom_vline(xintercept = 20.5, linetype = "dotted")
ggsave("../Compound_work/HMC3 scRNAseq/Viz/newclust_stdmarkers_dotplot.png", units = "mm", width = 150, height = 250, dpi = 300)

#annotation of the diffex genes for these clusters#
total_num = 100
FDR_thresh = 0.05
genes = rownames(dat_seurat_chooseR)
entrez_names = mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")
out_stub = "../Compound_work/HMC3 scRNAseq/Viz/Anno/"

#first, reactome#
genes_grouped_up = list()
genes_grouped_down = list()
markers_up <- markers_up[markers_up$num_down_types > 3,]
markers_down <- markers_down[markers_down$num_up_types > 3,]
for(i in unique(markers_down$down_type)){
  #up#
  module_df = data.frame(Entrez = entrez_names, group = 20)
  module_df$group[genes %in% markers_up$gene[markers_up$up_type == unique(markers_up$up_type)[i]][1:total_num]] = i
  module_df = module_df[!(is.na(module_df$Entrez)),]
  genes_grouped_up[[i]] = module_df$Entrez[module_df$group == i]
  #down#
  module_df = data.frame(Entrez = entrez_names, group = 20)
  module_df$group[genes %in% markers_down$gene[markers_down$down_type == unique(markers_down$down_type)[i]][1:total_num]] = i
  module_df = module_df[!(is.na(module_df$Entrez)),]
  genes_grouped_down[[i]] = module_df$Entrez[module_df$group == i]
}
names(genes_grouped_up) = as.numeric(levels(Idents(dat_seurat_chooseR))); names(genes_grouped_down) = as.numeric(levels(Idents(dat_seurat_chooseR)))

upsamp <- compareCluster(geneCluster = genes_grouped_up, fun ="enrichPathway", )
dotplot(upsamp, color = "qvalue", label_format = 30) + geom_point(aes(size=GeneRatio), shape = 21, colour="black", stroke=0.5) + 
  theme(legend.key = element_rect(fill = "white", colour = "black")) + scale_fill_gradient(low = "blue", high = "red", guide = guide_colorbar(frame.colour = "black"))
ggsave(paste0(out_stub, "up_FDR0.05_Reactome_dotplot_all.png"), units = "mm", width = 400, height = 350, dpi = 300)

downsamp <- compareCluster(geneCluster = genes_grouped_down, fun ="enrichPathway")
dotplot(downsamp, color = "qvalue", label_format = 30) + geom_point(aes(size=GeneRatio), shape = 21, colour="black", stroke=0.5) + 
  theme(legend.key = element_rect(fill = "white", colour = "black")) + scale_fill_gradient(low = "blue", high = "red", guide = guide_colorbar(frame.colour = "black"))
ggsave(paste0(out_stub, "down_FDR0.05_Reactome_dotplot_all.png"), units = "mm", width = 400, height = 350, dpi = 300)

#then GO#
bg_genes <- rownames(dat_seurat_chooseR)
BP_up = list()
BP_down = list()
for(ident in unique(markers_down$down_type)){
  upreg_genes = markers_up$gene[markers_up$up_type == ident][1:total_num]
  if(!(identical(upreg_genes, character(0)))){
    BP_up[[ident]] <- GO_topGO(top_x = 30, ont = "BP", bg_genes = bg_genes, diffex_genes = upreg_genes, 
                               condition_name = paste0(ident, "_up"), 
                               outdir = "../Compound_work/HMC3 scRNAseq/Viz/Anno/", gene_ID = "symbol")
  }
  downreg_genes = markers_down$gene[markers_down$down_type == ident][1:total_num]
  if(!(identical(downreg_genes, character(0)))){
    BP_down[[ident]] <- GO_topGO(top_x = 30, ont = "BP", bg_genes = bg_genes, diffex_genes = downreg_genes, 
                                 condition_name = paste0(ident, "_down"), 
                                 outdir = "../Compound_work/HMC3 scRNAseq/Viz/Anno/", gene_ID = "symbol")
  }
}

#here, I looked at marker genes from these clusters in this data vs. original data to look at similarity between original and new clusters#
c1_markers <- c("LDHA", "TMSB10", "GAPDH", "LGALS1", "TSPO")
c2_markers <- c("RRM2", "ZWINT", "CLSPN", "EIF6", "NIFK")
c3_markers <- c("ATF4", "GAS5", "IGBP1", "CD44", "LGALS3")
c4_markers <- c("CDKN3", "CENPA", "ACTG1", "HMGB2", "DLGAP5")
c5_markers <- c("CENPX", "B2M", "IFI30", "HLA-B", "TIMP1")
c6_markers <- c("NPPB", "CCN2", "TAGLN", "CALD1", "TMSB4X")
c7_markers <- c("ADAMTS1", "BOD1L1", "DDX17", "FUS", "NOTCH2")
c8_markers <- c("AP1S1", "HRAS", "LAMTOR2", "CD81", "S100A11")
c9_markers <- c("SNAPC1", "GADD45B", "MAFG", "CCNL1", "TOP1")
c10_markers <- c("ALCAM", "COL4A1","HSP90B1", "ITGB1", "TGFB2")
c11_markers <- c("CDKN1A", "PLTP", "IGFBP7", "CD82", "CAV1")
c12_markers <- c("IFIT3", "FLOT1", "VIM", "ISG20", "SRRM1")
c13_markers <- c("CDK1", "UBE2C", "EIF4A1", "TPX2", "FGF5")

#genes with the new clusters#
final = rev(unique(c(c1_markers, c2_markers, c3_markers, c4_markers, c5_markers, c6_markers, c7_markers, c8_markers, c9_markers, c10_markers, c11_markers, c12_markers, c13_markers)))
DotPlot(dat_seurat_chooseR, assay = "RNA", features = final, dot.min = 0.1, cluster.idents = F, scale = T) + coord_flip() +
  scale_colour_viridis(option="magma", guide = guide_colorbar(frame.colour = "black", ticks = TRUE)) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.title = element_blank()) + scale_y_discrete(position = "right") +
  geom_vline(xintercept = 5.5, linetype = "dotted") + geom_vline(xintercept = 10.5, linetype = "dotted") + geom_vline(xintercept = 15.5, linetype = "dotted") + geom_vline(xintercept = 20.5, linetype = "dotted") + 
  geom_vline(xintercept = 25.5, linetype = "dotted") + geom_vline(xintercept = 30.5, linetype = "dotted") + geom_vline(xintercept = 35.5, linetype = "dotted") + geom_vline(xintercept = 40.5, linetype = "dotted") +
  geom_vline(xintercept = 45.5, linetype = "dotted") + geom_vline(xintercept = 50.5, linetype = "dotted") + geom_vline(xintercept = 55.5, linetype = "dotted") +  geom_vline(xintercept = 60.5, linetype = "dotted")
ggsave("../Compound_work/HMC3 scRNAseq/Viz/newclust_newmarkers_dotplot.png", units = "mm", width = 200, height = 500, dpi = 300)

#genes with the original microglial dataset#
DotPlot(uglia_only, assay = "RNA", features = final, dot.min = 0.1, cluster.idents = F, scale = T) + coord_flip() +
  scale_colour_viridis(option="magma", guide = guide_colorbar(frame.colour = "black", ticks = TRUE)) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.title = element_blank()) + scale_y_discrete(position = "right") +
  geom_vline(xintercept = 5.5, linetype = "dotted") + geom_vline(xintercept = 10.5, linetype = "dotted") + geom_vline(xintercept = 15.5, linetype = "dotted") + geom_vline(xintercept = 20.5, linetype = "dotted") + 
  geom_vline(xintercept = 25.5, linetype = "dotted") + geom_vline(xintercept = 30.5, linetype = "dotted") + geom_vline(xintercept = 35.5, linetype = "dotted") + geom_vline(xintercept = 40.5, linetype = "dotted") +
  geom_vline(xintercept = 45.5, linetype = "dotted") + geom_vline(xintercept = 50.5, linetype = "dotted") + geom_vline(xintercept = 55.5, linetype = "dotted") +  geom_vline(xintercept = 60.5, linetype = "dotted")
ggsave("../Compound_work/HMC3 scRNAseq/Viz/oldclust_newmarkers_dotplot.png", units = "mm", width = 200, height = 500, dpi = 300)

##now perform mitotyping analysis in the compound-stimulated data at the overall and cluster levels##
#load in the mitotype list#
mito_type_list = read.csv("../Metabolism/Gene-to-pathway-list_hm.csv")
plot_module_enrichment(seurat = dat_seurat_chooseR, module_list = mito_type_list, output_dir = "../Compound_work/HMC3 scRNAseq/Viz/Metabolism/Mito_")

###load in the separated cluster markers###
out_stub = "../Compound_work/HMC3 scRNAseq/pairwise_diffex/"
markers = read.csv(paste(out_stub, "RNA_MAST/microglia_pairwise_up.csv", sep = ""))
markers_down = read.csv(paste(out_stub, "RNA_MAST/microglia_pairwise_down.csv", sep = ""))
FDR_thresh = 0.01

###read in all genes present in the microglia-only dataset for background###
genes <- rownames(dat_seurat_chooseR)
allgenes = data.table(genes)
colnames(allgenes) = "gene"

plot_cluster_module_enrichment(seurat = dat_seurat_chooseR, module_list = mito_type_list, allgenes = allgenes, markers = markers, 
                               markers_down = markers_down, FDR_thresh = 0.01, output_dir = "../Compound_work/HMC3 scRNAseq/Viz/Metabolism/")

#plot stacked barplots for each of the clusters showing representation in different compound treatment situations
cols = c( "#EE7600", "gray30",  "#91CADB", "#7A67EE", "gray90")

region_by_type = prop.table(table(Idents(dat_seurat_chooseR), dat_seurat_chooseR$identity), margin = 1)
region_by_type <- region_by_type[,c("Untreated", "DMSO", "Narciclasine", "Torin2", "Camptothecin")]
stack_frame = data.frame(Matrix(ncol = 5, nrow = 0))
for(i in 1:dim(region_by_type)[[2]]){
  tempcol = region_by_type[,i]
  temp_df = cbind(tempcol, tempcol, tempcol)
  temp_df[,1] = as.numeric(tempcol)
  temp_df[,2] = names(tempcol)
  temp_df[,3] = colnames(region_by_type)[i]
  stack_frame = rbind(stack_frame, temp_df)
}
colnames(stack_frame) = c("proportion", "clust", "treatment")
stack_frame$clust = factor(stack_frame$clust, levels = c(1:13))
stack_frame$proportion = as.numeric(stack_frame$proportion)
stack_frame$treatment = factor(stack_frame$treatment,levels = unique(stack_frame$treatment))
ggplot(stack_frame, aes(fill = treatment, y=proportion, x=clust)) + geom_bar(position = "stack", stat = "identity", colour = "black") + theme_classic() +
  scale_y_continuous(expand = c(0,0))+ ylab("Proportion of each cluster") + scale_fill_manual(values=cols) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.ticks.x = element_blank(), axis.text.y = element_text(size = 30, color = "black"), plot.margin = margin(15, 10, 10, 10),
        legend.key.size = unit(0.5, 'in'), legend.title = element_text(size=15),legend.text = element_text(size=15))
ggsave("../Compound_work/HMC3 scRNAseq/Viz/prop_tx.png", units = "mm", width = 300, height = 200, dpi = 300)

# DimPlot(dat_seurat_chooseR, pt.size = 2) + ggtitle("ChooseR optimal clustering v1") + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 40)) + NoAxes()
DimPlot(dat_seurat_chooseR, pt.size = 2, group.by = "identity", cols = cols) + ggtitle("Compound identity") + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 30)) + NoAxes()
ggsave("../Compound_work/HMC3 scRNAseq/Viz/compound_identity.png", units = "mm", width = 450, height = 350, dpi = 300)

#specifically, plot the module scores of several clusters on this data, starting with cluster 10#
old_markers <- read.csv("2020_Analysis_Outs/pairwise_diffex/downsampled/split/RNA_MAST/microglia_pairwise_up.csv")
c10_markers <- old_markers$gene[old_markers$up_type == 10][1:50]
c10_markers <- list(c10_markers)
dat_seurat_chooseR = AddModuleScore(object = dat_seurat_chooseR, features = c10_markers, name = "10.module")
print(FeaturePlot(dat_seurat_chooseR, features = "HLA.high.score1", pt.size = 1) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), 
                                                                                                           guide = guide_colorbar(label = TRUE,
                                                                                                                                  draw.ulim = TRUE, 
                                                                                                                                  draw.llim = TRUE,
                                                                                                                                  frame.colour = "black", ticks = TRUE, 
                                                                                                                                  label.position = "bottom",
                                                                                                                                  barwidth = 18,
                                                                                                                                  barheight = 1.3, 
                                                                                                                                  direction = 'horizontal')) +
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position = "bottom", text=element_text(size = 20)) + NoAxes() + ggtitle("Cluster 10 score"))
ggsave("../Compound_work/HMC3 scRNAseq/Viz/c10_module_score.png", units = "mm", width = 300, height = 300, dpi = 300)

c16_markers <- intersect(old_markers$gene[old_markers$up_type == 1], old_markers$gene[old_markers$up_type == 6])[1:100]
c16_markers <- list(c16_markers)
dat_seurat_chooseR = AddModuleScore(object = dat_seurat_chooseR, features = c16_markers, name = "1.6.module")
print(FeaturePlot(dat_seurat_chooseR, features = "X1.6.module1", pt.size = 1) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), 
                                                                                                        guide = guide_colorbar(label = TRUE,
                                                                                                                               draw.ulim = TRUE, 
                                                                                                                               draw.llim = TRUE,
                                                                                                                               frame.colour = "black", ticks = TRUE, 
                                                                                                                               label.position = "bottom",
                                                                                                                               barwidth = 18,
                                                                                                                               barheight = 1.3, 
                                                                                                                               direction = 'horizontal')) +
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position = "bottom", text=element_text(size = 20)) + NoAxes() + ggtitle("Cluster 1/6 score"))
ggsave("../Compound_work/HMC3 scRNAseq/Viz/c1+6_module_score.png", units = "mm", width = 300, height = 300, dpi = 300)
for(clust in 1:12){
  clust_markers <- old_markers$gene[old_markers$up_type == clust][1:50]
  clust_markers <- list(clust_markers)
  dat_seurat_chooseR = AddModuleScore(object = dat_seurat_chooseR, features = clust_markers, name = paste0(clust, "_module"))
  
  print(FeaturePlot(dat_seurat_chooseR, features = paste0("X",clust, "_module1"), pt.size = 1) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), 
                                                                                                                         guide = guide_colorbar(label = TRUE,
                                                                                                                                                draw.ulim = TRUE, 
                                                                                                                                                draw.llim = TRUE,
                                                                                                                                                frame.colour = "black", ticks = TRUE, 
                                                                                                                                                label.position = "bottom",
                                                                                                                                                barwidth = 18,
                                                                                                                                                barheight = 1.3, 
                                                                                                                                                direction = 'horizontal')) +
          theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position = "bottom", text=element_text(size = 20)) + NoAxes() + ggtitle(paste0("Cluster ", clust, " score")))
  ggsave(paste0("../Compound_work/HMC3 scRNAseq/Viz/Cluster_module_scores/", clust, "_module_score.png"), units = "mm", width = 300, height = 300, dpi = 300)
}

#assemble and plot signature scores for different perturbations
old_markers <- read.csv("2020_Analysis_Outs/pairwise_diffex/downsampled/split/RNA_MAST/microglia_pairwise_up.csv")

c16_markers <- intersect(old_markers$gene[old_markers$up_type == 1], old_markers$gene[old_markers$up_type == 6])[1:100]
c16_markers <- list(c16_markers)
dat_seurat_chooseR = AddModuleScore(object = dat_seurat_chooseR, features = c16_markers, name = "1.6.module")
c49_markers <- intersect(old_markers$gene[old_markers$up_type == 4], old_markers$gene[old_markers$up_type == 9])[1:100]
c49_markers <- list(c49_markers)
dat_seurat_chooseR = AddModuleScore(object = dat_seurat_chooseR, features = c49_markers, name = "4.9.module")
c8_markers <- old_markers$gene[old_markers$up_type == 8][1:100]
c8_markers <- list(c8_markers)
dat_seurat_chooseR = AddModuleScore(object = dat_seurat_chooseR, features = c8_markers, name = "8.module")
c10_markers <- old_markers$gene[old_markers$up_type == 10][1:100]
c10_markers <- list(c10_markers)
dat_seurat_chooseR = AddModuleScore(object = dat_seurat_chooseR, features = c10_markers, name = "10.module")
c11_markers <- old_markers$gene[old_markers$up_type == 11][1:100]
c11_markers <- list(c11_markers)
dat_seurat_chooseR = AddModuleScore(object = dat_seurat_chooseR, features = c11_markers, name = "11.module")
c12_markers <- old_markers$gene[old_markers$up_type == 12][1:100]
c12_markers <- list(c12_markers)
dat_seurat_chooseR = AddModuleScore(object = dat_seurat_chooseR, features = c12_markers, name = "12.module")

#visualize the results per gene target#
df <- data.frame(Idents(dat_seurat_chooseR), dat_seurat_chooseR$identity, dat_seurat_chooseR@meta.data[,grep("module", colnames(dat_seurat_chooseR@meta.data))])
colnames(df) <- c("cluster", "treatment", paste0("module", c("1/6", "4/9", "8", "10", "11", "12")))
df$cluster <- factor(df$cluster, levels = c(1:13))
cols = c( "#EE7600", "gray30",  "#91CADB", "#7A67EE", "gray90")

plot_list <- lapply(colnames(df)[grep("module", colnames(df))], FUN = function(x) {
  ggplot(df, aes(fill=cluster, y=df[,x], x=cluster)) + 
    geom_violin(position="dodge", alpha=0.8) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.title.x = element_blank(), legend.position = "none") + ylab(x)
})
wrap_plots(plot_list, ncol = 1, nrow = 6)
ggsave("../Compound_work/HMC3 scRNAseq/Viz/cluster_module_scores_by_cluster.png", units = "mm", width = 300, height = 225, dpi = 300)

plot_list <- lapply(colnames(df)[grep("module", colnames(df))], FUN = function(x) {
  ggplot(df, aes(fill=treatment, y=df[,x], x=treatment)) + 
    geom_violin(position="dodge", alpha=0.8) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.title.x = element_blank(), legend.position = "none") + ylab(x) + scale_fill_manual(values = cols)
})
wrap_plots(plot_list, ncol = 1, nrow = 6)
ggsave("../Compound_work/HMC3 scRNAseq/Viz/cluster_module_scores_by_treatment.png", units = "mm", width = 150, height = 300, dpi = 300)

#save the p values for the above analyses#
treatment_p <- list()
for(i in paste0("module", c("1/6", "4/9", "10", "8", "11", "12"))){
  stat.test <- df %>% pairwise_wilcox_test(as.formula(paste0(i, " ~ treatment")), 
                                           p.adjust.method = "holm", exact = FALSE, 
                                           detailed = F, conf.level = 0.999)
  treatment_p[[i]] <- stat.test
}
tx_wb <- createWorkbook()
lapply(seq_along(treatment_p), function(x){
  addWorksheet(wb = tx_wb, sheetName = gsub(names(treatment_p)[x], pattern = "/", replacement = "."))
  writeData(wb = tx_wb, sheet = x, treatment_p[[x]])
})
saveWorkbook(tx_wb, "../Compound_work/HMC3 scRNAseq/Viz/cluster_module_scores_by_treatment.xlsx", overwrite = TRUE)

module_p <- list()
for(i in paste0("module", c("1/6", "4/9", "10", "8", "11", "12"))){
  stat.test <- df %>% pairwise_wilcox_test(as.formula(paste0(i, " ~ cluster")), 
                                           p.adjust.method = "holm", exact = FALSE,
                                           detailed = F, conf.level = 0.999)
  module_p[[i]] <- stat.test
}
clust_wb <- createWorkbook()
lapply(seq_along(module_p), function(x){
  addWorksheet(wb = clust_wb, sheetName = gsub(names(module_p)[x], pattern = "/", replacement = "."))
  writeData(wb = clust_wb, sheet = x, module_p[[x]])
})
saveWorkbook(clust_wb, "../Compound_work/HMC3 scRNAseq/Viz/cluster_module_scores_by_cluster.xlsx", overwrite = TRUE)


##GSEA for more quantitative assessment of annotation scores per treatment condition, etc.##
#load in the separated cluster markers#
diffex_stub = "2020_Analysis_Outs/pairwise_diffex/downsampled/split/"
old_markers = read.csv(paste(diffex_stub, "RNA_MAST/microglia_pairwise_up.csv", sep = ""))
num_genes = 100

#set up lists for testing#
c16_markers <- intersect(old_markers$gene[old_markers$up_type == 1], old_markers$gene[old_markers$up_type == 6])[1:num_genes]
c49_markers <- intersect(old_markers$gene[old_markers$up_type == 4], old_markers$gene[old_markers$up_type == 9])[1:num_genes]
c10_markers <- old_markers$gene[old_markers$up_type == 10][1:num_genes]
c8_markers <- old_markers$gene[old_markers$up_type == 8][1:num_genes]
c11_markers <- old_markers$gene[old_markers$up_type == 11][1:num_genes]
c12_markers <- old_markers$gene[old_markers$up_type == 12][1:num_genes]
module_table <- data.frame(term = as.character(sapply(c("1/6", "4/9", "8", "10", "11", "12"), FUN = function(x) {rep(x, num_genes)})),
                           gene = c(c16_markers, c49_markers, c8_markers,c10_markers,c11_markers,c12_markers))
#collapse and aggregate control conditions#
dat_seurat_chooseR$simplified_tx <- fct_collapse(
  factor(dat_seurat_chooseR$identity, levels = c("DMSO", "Untreated", "Camptothecin", "Torin2", "Narciclasine")), 
  Ctrl = c("DMSO", "Untreated")
)
dat_seurat_chooseR <- NormalizeData(dat_seurat_chooseR)
uglia_colors = c("#91CADB", "turquoise3", "darkorange2", "#3D643D", "magenta2", "black")
names(uglia_colors) <- c("1/6", "4/9", "8", "10", "11", "12")

#perform GSEA analysis and visualization for each of the singular treatment conditions#
for(drug in c("Camptothecin", "Torin2", "Narciclasine")){
  drug_lfc <- FoldChange(object = dat_seurat_chooseR, 
                         group.by = "simplified_tx",
                         ident.1 = drug,
                         ident.2 = "Ctrl",
                         assay = "RNA", 
                         slot = "data")
  drug_lfc <- drug_lfc[order(drug_lfc$avg_log2FC,decreasing = T),]
  gene_list <- drug_lfc$avg_log2FC
  names(gene_list) <- rownames(drug_lfc)
  gsea_result <- GSEA(geneList = gene_list,
                      exponent = 1,
                      minGSSize = 10,
                      maxGSSize = 500,
                      eps = 0,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "BH",
                      TERM2GENE = module_table,
                      verbose = T,
                      by = "fgsea"
  )
  for(j in unique(gsea_result@result$ID)){
    print(gseaplot(gsea_result, geneSetID = j, title = paste0("Cluster ", j, " signature"), color.line = uglia_colors[j]))
    ggsave(filename = paste0("../Compound_work/HMC3 scRNAseq/Viz/GSEA/", drug, "_", gsub(x = j, pattern = "/", replacement = "."), "_GSEA.png"), 
           units = "mm", width = 200, height = 350, dpi = 300)
  }
  if("1/6" %in% unique(gsea_result@result$ID)){
    tmp <- gseaplot2(gsea_result, geneSetID = c("1/6", "10"), subplots = c(1),
                     title = paste0(drug,": all signatures"), color = uglia_colors[c("1/6", "10")], pvalue_table = F, base_size = 20) + 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    tmp2 <- gseaplot2(gsea_result, geneSetID = c("1/6", "10"), subplots = c(2), color = uglia_colors[ c("1/6", "10")], base_size = 20) + 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
    tmp/tmp2 + plot_layout(heights = c(4,1))
    ggsave(filename = paste0("../Compound_work/HMC3 scRNAseq/Viz/GSEA/", drug, "_10vs1.6_GSEA.png"), 
           units = "mm", width = 200, height = 300, dpi = 300)
    #manually add a  p-value table back in#
    pd <- gsea_result[c("1/6", "10"), c("Description", "pvalue", "p.adjust")]
    rownames(pd) <- pd$Description
    pd <- pd[,-1]
    pd <- round(pd, 4)
    tp <- ggtable(pd, tmp)
    print(tp)
    ggsave(filename = paste0("../Compound_work/HMC3 scRNAseq/Viz/GSEA/", drug, "_10vs1.6_table.png"), 
           units = "mm", width = 50, height = 50, dpi = 300)
  }
  
  if(length(unique(gsea_result@result$ID)) > 0){
    tmp <- gseaplot2(gsea_result, geneSetID = unique(gsea_result@result$ID), subplots = c(1),
                     title = paste0(drug,": all signatures"), color = uglia_colors[unique(gsea_result@result$ID)], pvalue_table = F, base_size = 20) + 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    tmp2 <- gseaplot2(gsea_result, geneSetID = unique(gsea_result@result$ID), subplots = c(2), color = uglia_colors[unique(gsea_result@result$ID)], base_size = 20) + 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
    tmp/tmp2 + plot_layout(heights = c(4,1))
    ggsave(filename = paste0("../Compound_work/HMC3 scRNAseq/Viz/GSEA/", drug, "_all_GSEA.png"), 
           units = "mm", width = 200, height = 300, dpi = 300)
    #manually add a  p-value table back in#
    pd <- gsea_result[unique(gsea_result@result$ID), c("Description", "pvalue", "p.adjust")]
    rownames(pd) <- pd$Description
    pd <- pd[,-1]
    pd <- round(pd, 4)
    tp <- ggtable(pd, tmp)
    print(tp)
    ggsave(filename = paste0("../Compound_work/HMC3 scRNAseq/Viz/GSEA/", drug, "_all_table.png"), 
           units = "mm", width = 50, height = 50, dpi = 300)
  }
}

#GSEA for mitotypes - same as above, but for mitotype scores instead of for original cluster module scores#
mito_type_list = read.csv("../Metabolism/Gene-to-pathway-list_hm.csv")
colnames(mito_type_list) <- c("gene", "term")
mito_type_list <- mito_type_list[,c("term", "gene")]
for(drug in c("Camptothecin", "Torin2", "Narciclasine")){
  drug_lfc <- FoldChange(object = dat_seurat_chooseR,
                         group.by = "simplified_tx",
                         ident.1 = drug,
                         ident.2 = "Ctrl",
                         assay = "RNA", 
                         slot = "data")
  drug_lfc <- drug_lfc[order(drug_lfc$avg_log2FC,decreasing = T),]
  gene_list <- drug_lfc$avg_log2FC
  names(gene_list) <- rownames(drug_lfc)
  gsea_result <- GSEA(geneList = gene_list,
                      exponent = 1,
                      minGSSize = 5,
                      maxGSSize = 500,
                      eps = 0,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "BH",
                      TERM2GENE = mito_type_list,
                      verbose = T,
                      by = "fgsea"
  )
  for(j in unique(gsea_result@result$ID)){
    print(gseaplot(gsea_result, geneSetID = j, title = paste0(j, " signature")))
    ggsave(filename = paste0("../Compound_work/HMC3 scRNAseq/Viz/GSEA/Metabolism/", drug, "_", gsub(x = j, pattern = "/", replacement = "."), "_plot.png"), 
           units = "mm", width = 200, height = 350, dpi = 300)
  }
  if(length(unique(gsea_result@result$ID)) > 0){
    tmp <- gseaplot2(gsea_result, geneSetID = unique(gsea_result@result$ID), subplots = c(1),
                     title = paste0(drug,": all signatures"), pvalue_table = F, base_size = 20) + 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
    tmp2 <- gseaplot2(gsea_result, geneSetID = unique(gsea_result@result$ID), subplots = c(2), base_size = 20) + 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
    tmp/tmp2 + plot_layout(heights = c(4,1))
    ggsave(filename = paste0("../Compound_work/HMC3 scRNAseq/Viz/GSEA/Metabolism/", drug, "_all.png"), 
           units = "mm", width = 200, height = 300, dpi = 300)
    #manually add a  p-value table back in#
    pd <- gsea_result[unique(gsea_result@result$ID), c("Description", "pvalue", "p.adjust")]
    rownames(pd) <- pd$Description
    pd <- pd[,-1]
    pd <- round(pd, 4)
    tp <- ggtable(pd, tmp)
    print(tp)
    ggsave(filename = paste0("../Compound_work/HMC3 scRNAseq/Viz/GSEA/Metabolism/", drug, "all_table.png"), 
           units = "mm", width = 150, height = 150, dpi = 300)
  }
}