#revised GBM slice analysis#
library(Matrix)
require(Seurat)
require(dplyr)
require(SeuratWrappers)
library(ggplot2)
library(harmony)
library(viridis)
library(grid)
library(ggtree)
library(tidyverse)
library(org.Hs.eg.db)
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(RColorBrewer)
library(patchwork)
library(EnhancedVolcano)
library(pheatmap)
library(data.table)

setwd("Microglia_project/Analysis_2020/Full_uglia_analysis_2020/")
#color_palette
c26 <- c(
  "gray70",
  "gray20",
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown", "black"
  
)
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

#read in files
allmat=list()
allgen=c()
allcells=c()
filenames = list.files(path = "../Data/Wenting_slice_myeloid/Data/")
# umifilter=1000  ###UMI threshold to keep cells/nuclei. Put at 1000 for uglia.
# mito_pct_filter=10 ##mitochondrial read percentage filter
umimin = 0
umimax = 12500
#generate individual matrices per dataset
pdf(file =paste0("../Data/Wenting_slice_myeloid/cell_metrics.pdf"), onefile = T, paper = "a4r", width = 8.5, height = 11)
for (i in 1:length(filenames)) {
  tempname=paste("../Data/Wenting_slice_myeloid/Data/", filenames[i], sep ="")
  temp=read.table(file=tempname, sep="\t",skip=2,as.is=T)
  temp = temp[,-1]
  ###drop the ENS IDs, zero rows, and sum duplicate genes
  temp <- temp[order(temp[,1]),]
  temp = temp[rowSums(temp[, -1])>0, ]
  summed = rowsum(temp[,-1], temp[,1])
  genes=rownames(summed)
  ###filter RP genes and apply a lax UMI filter###
  filter_genes = grep("^RP[0-9]|^BC[0-9]|^RPL|^RPS|-PS",genes) ###variant to remove all RPs
  summed = summed[-filter_genes,]
  names = strsplit(filenames[[i]], "[.]")[[1]][1]
  names = gsub("-","_",names)
  colnames(summed)=paste(names,"-",1:ncol(summed),sep="")
  umicount=colSums(summed)
  hist(colSums(summed), xlim = c(0,20000), breaks = 100, main = names)
  
  keepcols=which(umicount>umimin)
  summed=summed[,keepcols]
  umicount=colSums(summed)
  keepcols=which(umicount<umimax)
  summed=summed[,keepcols]
  tempmat=as.matrix(summed)
  allmat[[i]]=Matrix(tempmat,sparse=T)
  allgen=c(allgen,rownames(tempmat))
  allcells=c(allcells,colnames(tempmat))
  print(names)
}
dev.off()

# #merge all matrices
allgen=unique(allgen)
alldat=Matrix(0,nrow=length(allgen),ncol=length(allcells),sparse=T)
rownames(alldat)=allgen
colnames(alldat)=allcells
for (ii in 1:length(allmat)) {
  alldat[rownames(allmat[[ii]]),colnames(allmat[[ii]])]=allmat[[ii]]
}
rm(genes, summed, temp, tempmat, allcells, allgen, allmat, outval); gc()
saveRDS(alldat, "GBMmat_0+12500_filterRP.rds")
#reload and process#
alldat = readRDS("GBMmat_0+12500_filterRP.rds")
GBM.seurat = SeuratObject::CreateSeuratObject(counts = alldat, min.cells = 3, names.delim = "-")
GBM.seurat$orig.ident <- sapply(strsplit(names(GBM.seurat$orig.ident), "-"), "[", 1)

###add metadata###
##fix up the naming of the drug metabata and then add to seurat object##
drug_pert_metadata = data.frame(read.csv("../Data/Wenting_slice_myeloid/Name_Drug_matching.csv", header = T))[,1:3]
colnames(drug_pert_metadata) = c("Name", "Drug", "TB_num")
drug_pert_metadata$Name = gsub(drug_pert_metadata$Name, pattern = "-", replacement = "_")
#not all drugs are present in the data I have currently#
# drug_pert_metadata = drug_pert_metadata[drug_pert_metadata$Name %in% unique(GBM.seurat$orig.ident),]
#create a metadata column that can be used to add to the seurat object#
final_pert_metadata = unlist(sapply(unique(GBM.seurat$orig.ident), 
                                    FUN = function(x){rep(drug_pert_metadata$Drug[drug_pert_metadata$Name == x], table(GBM.seurat$orig.ident)[x])}))
names(final_pert_metadata) = colnames(GBM.seurat)
GBM.seurat = AddMetaData(object = GBM.seurat, metadata = final_pert_metadata, col.name = "Drug_Tx")
#metadata for sequencing run #
seq_id = sapply(strsplit(colnames(GBM.seurat), split = "_"), "[", 1)
names(seq_id) = colnames(GBM.seurat)
GBM.seurat = AddMetaData(object = GBM.seurat, metadata = seq_id, col.name = "seq_id")
#metadata for tumor bank ID#
TB_metadata = unlist(sapply(unique(GBM.seurat$orig.ident), 
                            FUN = function(x){rep(drug_pert_metadata$TB_num[drug_pert_metadata$Name == x], table(GBM.seurat$orig.ident)[x])}))
names(TB_metadata) = colnames(GBM.seurat)
GBM.seurat = AddMetaData(object = GBM.seurat, metadata = TB_metadata, col.name = "TB_ID")

#Seurat integration#
split.list <- SplitObject(GBM.seurat, split.by = "TB_ID")
split.list <- lapply(X = split.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = split.list, nfeatures = 3000)
split.list <- PrepSCTIntegration(object.list = split.list, anchor.features = features)
GBM.anchors <- FindIntegrationAnchors(object.list = split.list, normalization.method = "SCT",
                                      anchor.features = features)
GBM.sct.combined <- IntegrateData(anchorset = GBM.anchors, normalization.method = "SCT")
GBM.sct.combined <- RunPCA(GBM.sct.combined, verbose = FALSE)
GBM.sct.combined <- RunUMAP(GBM.sct.combined, reduction = "pca", dims = 1:30)
GBM.sct.combined <- FindNeighbors(GBM.sct.combined, reduction = "pca", dims = 1:30)
GBM.sct.combined <- FindClusters(GBM.sct.combined, resolution = 0.5)
DimPlot(GBM.sct.combined)
ggsave("../Data/Wenting_slice_myeloid/Viz/final_UMAP.png", units = "in", width = 4.5, height = 4.5, dpi = 300)
DimPlot(GBM.sct.combined, group.by = "TB_ID")
ggsave("../Data/Wenting_slice_myeloid/Viz/final_TB_ID.png", units = "in", width = 4.5, height = 4.5, dpi = 300)
DimPlot(GBM.sct.combined, group.by = "Drug_Tx")
ggsave("../Data/Wenting_slice_myeloid/Viz/final_Drug_ID.png", units = "in", width = 4.5, height = 4.5, dpi = 300)
DimPlot(GBM.sct.combined, group.by = "Drug_Tx", split.by = "TB_ID")
ggsave("../Data/Wenting_slice_myeloid/Viz/final_split_TB_Drug.png", units = "in", width = 12, height = 4.5, dpi = 300)

#---------------------------------------------------------------------------------#
#load back in and restart analysis#
load("../Data/Wenting_slice_myeloid/final_GBM_slicemyeloid.rda")
GBM.sct.combined <- SetIdent(GBM.sct.combined, value =  "TB_ID")
DefaultAssay(GBM.sct.combined) <- "RNA"

c19 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "steelblue4",
  "darkturquoise", "yellow4",
  "brown"
)

DimPlot(GBM.sct.combined, split.by = "Drug_Tx", group.by = "TB_ID", pt.size = 0.5, ncol = 7, cols = c19) + ggtitle("Drug treatments") + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 15)) + NoAxes()
ggsave("../Data/Wenting_slice_myeloid/Viz/final_drug_tx_TBID.png", units = "in", width = 14, height = 10, dpi = 300)

#let's compute differential expression for each TBID in comparison to control here
diffex_list <- list()
# GBM.sct.combined <- subset(GBM.sct.combined, idents = drug_TB_IDs[drug_TB_IDs %in% unique(Idents(GBM.sct.combined))], invert = T)
for(TB_id in unique(GBM.sct.combined$TB_ID)){
  sub_seurat <- subset(GBM.sct.combined, idents = TB_id)
  unique_drugs <- unique(sub_seurat$Drug_Tx)
  sub_seurat <- NormalizeData(sub_seurat)
  diffex_list[[TB_id]] <- lapply(unique_drugs[unique_drugs != "Ctrl"], FUN = function(x){
    FindMarkers(object = sub_seurat, test.use = "MAST", ident.1 = "Ctrl", ident.2 = x, min.pct = 0.05, 
                only.pos = F, assay = "RNA", slot = "data", logfc.threshold = 0.25, group.by = "Drug_Tx")
  })
  names(diffex_list[[TB_id]]) <- unique_drugs[unique_drugs != "Ctrl"]
  dir.create(paste0("../Data/Wenting_slice_myeloid/comparative_diffex/", trimws(TB_id)))
  for(name in names(diffex_list[[TB_id]])){
    write.csv(x = diffex_list[[TB_id]][[name]], file = paste0("../Data/Wenting_slice_myeloid/comparative_diffex/", trimws(TB_id), "/Ctrl_vs_", name, "_markers.csv"))
  }
}

#---------------------------------------------to run everytime you may restart the analysis with different drug here---------------------------------------------------------------#
#load back in and restart analysis#
load("../Data/Wenting_slice_myeloid/final_GBM_slicemyeloid.rda")
GBM.sct.combined <- SetIdent(GBM.sct.combined, value =  "TB_ID")
DefaultAssay(GBM.sct.combined) <- "RNA"

#load in the diffex lists#
diffex_list <- readRDS("../Data/Wenting_slice_myeloid/comparative_diffex/diffex.rds")
#let's do a basic analysis with topotecan first, for example#
drug1 <- "Topotecan" #change for a different drug target#


#
npcs <- 20
names(c26) <- unique(GBM.sct.combined$Drug_Tx)
mito_type_list = read.csv("../Metabolism/Gene-to-pathway-list_hm.csv")
#
uglia_colors = c("#91CADB", "turquoise3", "darkorange2", "#3D643D", "magenta2", "black")
names(uglia_colors) <- c("1/6", "4/9", "8", "10", "11", "12")
rwb <- colorRampPalette(colors = c("blue", "white", "red"))

#
drug_pert_metadata = data.frame(read.csv("../Data/Wenting_slice_myeloid/Name_Drug_matching.csv", header = T))[,1:3]
# drug_pert_metadata$Drug_Tx <- trimws(drug_pert_metadata$Drug_Tx)
drug_TB_IDs <- drug_pert_metadata$X.TB[drug_pert_metadata$Drug_Tx == drug1]
diffex_subset <- diffex_list[names(diffex_list) %in% drug_TB_IDs]
diffex_matrices <- lapply(diffex_subset, FUN = function(x){return(x[[drug1]])})

#set up prior cluster module scores#
diffex_stub = "2020_Analysis_Outs/pairwise_diffex/downsampled/split/"
old_markers = read.csv(paste(diffex_stub, "RNA_MAST/microglia_pairwise_up.csv", sep = ""))
num_genes = 100

#set up lists for testing with GSEA#
c16_markers <- intersect(old_markers$gene[old_markers$up_type == 1], old_markers$gene[old_markers$up_type == 6])[1:num_genes]
c49_markers <- intersect(old_markers$gene[old_markers$up_type == 4], old_markers$gene[old_markers$up_type == 9])[1:num_genes]
c8_markers <- old_markers$gene[old_markers$up_type == 8][1:num_genes]
c10_markers <- old_markers$gene[old_markers$up_type == 10][1:num_genes]
c11_markers <- old_markers$gene[old_markers$up_type == 11][1:num_genes]
c12_markers <- old_markers$gene[old_markers$up_type == 12][1:num_genes]
module_table <- data.frame(term = as.character(sapply(c("1/6", "4/9", "8", "10", "11", "12"), FUN = function(x) {rep(x, num_genes)})),
                           gene = c(c16_markers, c49_markers, c8_markers,c10_markers,c11_markers,c12_markers))
#for input into seurat module scoring#
marker_list <- list() 
marker_list[["1/6"]] <- c16_markers
marker_list[["4/9"]] <- c49_markers
marker_list[["8"]] <- c8_markers
marker_list[["10"]] <- c10_markers
marker_list[["11"]] <- c11_markers
marker_list[["12"]] <- c12_markers

#--------------------------------module scores-------------------------------------------------#
sub_seurat_list <- list()
for(TB_id in drug_TB_IDs){
  sub_seurat <- subset(GBM.sct.combined, idents = TB_id)
  sub_seurat <- NormalizeData(sub_seurat)
  sub_seurat <- sub_seurat %>% NormalizeData %>% SCTransform() %>%  RunPCA() %>% FindNeighbors(dims = 1:npcs) %>% RunUMAP(dims = 1:npcs)
  DimPlot(sub_seurat, group.by = "Drug_Tx", pt.size = 2, cols = c26[unique(sub_seurat$Drug_Tx)]) + ggtitle(paste0(trimws(TB_id), " drug treatments")) + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 15)) + NoAxes()
  ggsave(paste0("../Data/Wenting_slice_myeloid/comparative_viz/", trimws(TB_id), "/drug_tx_UMAP.png"), units = "in", width = 14, height = 10, dpi = 300)
  for(i in 1:length(marker_list)){
    module_name <- paste0(x= gsub(names(marker_list)[i], pattern = "/", replacement = "." ), ".module")
    sub_seurat = AddModuleScore(object = sub_seurat, features = marker_list[i], name = module_name)
    print(FeaturePlot(sub_seurat, features = paste0("X",module_name,"1")), pt.size = 4) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), 
                                                                                                                  guide = guide_colorbar(label = TRUE,
                                                                                                                                         draw.ulim = TRUE, 
                                                                                                                                         draw.llim = TRUE,
                                                                                                                                         frame.colour = "black", ticks = TRUE, 
                                                                                                                                         label.position = "bottom",
                                                                                                                                         barwidth = 18,
                                                                                                                                         barheight = 1.3, 
                                                                                                                                         direction = 'horizontal')) +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position = "bottom", text=element_text(size = 20)) + NoAxes() + 
      ggtitle(paste0(trimws(TB_id), " cluster ",names(marker_list)[i]," score"))
    ggsave(paste0("../Data/Wenting_slice_myeloid/comparative_viz/",trimws(TB_id), "/",module_name, "_score.png"), units = "in", width = 14, height = 11, dpi = 300)
  }
  df <- data.frame(sub_seurat$Drug_Tx, sub_seurat@meta.data[,grep("module", colnames(sub_seurat@meta.data))])
  colnames(df) <- c("drug", paste0("module", c("1/6", "4/9", "8", "10", "11", "12")))
  df$drug <- factor(df$drug, levels = c("Ctrl", unique(df$drug)[-c(unique(df$drug)=="Ctrl")][order(unique(df$drug)[-c(unique(df$drug)=="Ctrl")])]))
  plot_list <- lapply(colnames(df)[grep("module", colnames(df))], FUN = function(x) {
    ggplot(df, aes(y=df[,x], x=drug, fill = drug)) + 
      geom_violin(position="dodge", alpha=0.75) +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.title.x = element_blank()) + 
      ylab(x) + scale_fill_manual(values = c26[unique(sub_seurat$Drug_Tx)])
    
  })
  wrap_plots(plot_list, ncol = 1, nrow = 6)
  if(length(unique(sub_seurat$Drug_Tx)) > 3){
    ggsave(paste0("../Data/Wenting_slice_myeloid/comparative_viz/",trimws(TB_id), "/module_score_by_drug.png"), units = "in", width = 12, height = 12, dpi = 300)
  }
  else{
    ggsave(paste0("../Data/Wenting_slice_myeloid/comparative_viz/",trimws(TB_id), "/module_score_by_drug.png"), units = "in", width = 6, height = 12, dpi = 300)
  }
  dir.create(paste0("../Data/Wenting_slice_myeloid/comparative_viz/", trimws(TB_id), "/Metabolism"))
  plot_module_enrichment(seurat = sub_seurat, module_list = mito_type_list, output_dir = paste0("../Data/Wenting_slice_myeloid/comparative_viz/",trimws(TB_id), "/Metabolism/Mito_"))
  sub_seurat_list[[TB_id]] <- sub_seurat
}

#---------------------------------GSEA------------------------------------------------#
lfc_lists <- list()
gsea_lists <- list()
for(TB_id in drug_TB_IDs){
  dir.create(paste0("../Data/Wenting_slice_myeloid/comparative_viz/", trimws(TB_id), "/", trimws(drug1), "/GSEA"))
  # sub_seurat_list[[TB_id]] <- SetIdent(sub_seurat_list[[TB_id]], value =  "Drug_Tx")
  drug_lfc <- FoldChange(object = sub_seurat_list[[TB_id]], 
                         group.by = "Drug_Tx",
                         ident.1 = drug1,
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
                      pvalueCutoff = 0.5,
                      pAdjustMethod = "BH",
                      TERM2GENE = module_table,
                      verbose = T,
                      by = "fgsea"
  )
  lfc_lists[[TB_id]] <- drug_lfc
  gsea_lists[[TB_id]] <- gsea_result
}

colnames(mito_type_list) <- c("gene", "term")
mito_type_list <- mito_type_list[,c("term", "gene")]
metabo_enrich <- list()
for(TB_id in drug_TB_IDs){
  dir.create(paste0("../Data/Wenting_slice_myeloid/comparative_viz/", trimws(TB_id), "/", trimws(drug1), "/GSEA_metabolism/"))
  # sub_seurat_list[[TB_id]] <- SetIdent(sub_seurat_list[[TB_id]], value =  "Drug_Tx")
  drug_lfc <- FoldChange(object = sub_seurat_list[[TB_id]], 
                         group.by = "Drug_Tx",
                         ident.1 = drug1,
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
    ggsave(filename = paste0("../Data/Wenting_slice_myeloid/comparative_viz/", trimws(TB_id), "/", trimws(drug1), "/GSEA_metabolism/", gsub(x = j, pattern = "/", replacement = "."), "_plot.png"), 
           units = "mm", width = 200, height = 350, dpi = 300)
  }
  if(length(unique(gsea_result@result$ID)) > 0){
    tmp <- gseaplot2(gsea_result, geneSetID = unique(gsea_result@result$ID), subplots = c(1),
                     title = paste0(drug1,": all signatures"), pvalue_table = F, base_size = 20) + 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
    tmp2 <- gseaplot2(gsea_result, geneSetID = unique(gsea_result@result$ID), subplots = c(2), base_size = 20) + 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
    tmp/tmp2 + plot_layout(heights = c(4,1))
    ggsave(filename = paste0("../Data/Wenting_slice_myeloid/comparative_viz/", trimws(TB_id), "/", trimws(drug1), "/GSEA_metabolism/all_GSEA.png"), 
           units = "mm", width = 200, height = 300, dpi = 300)
    #manually add a  p-value table back in#
    pd <- gsea_result[unique(gsea_result@result$ID), c("Description", "pvalue", "p.adjust")]
    rownames(pd) <- pd$Description
    pd <- pd[,-1]
    pd <- round(pd, 4)
    tp <- ggtable(pd, tmp)
    print(tp)
    ggsave(filename = paste0("../Data/Wenting_slice_myeloid/comparative_viz/", trimws(TB_id), "/", trimws(drug1), "/GSEA_metabolism/all_table.png"), 
           units = "mm", width = 150, height = 150, dpi = 300)
  }
  metabo_enrich[[TB_id]] <- gsea_result
}
saveRDS(metabo_enrich, file = paste0("../Data/Wenting_slice_myeloid/comparative_anno/", trimws(drug1), "/GSEA_metabo_enrich.rds"))

#----------------------------------aggregate scored heatmap-----------------------------------------------#
#recreate sub-seurat and lfc lists#
sub_seurat_list <- list()
lfc_lists <- list()
for(TB_id in drug_TB_IDs){
  sub_seurat <- subset(GBM.sct.combined, idents = TB_id)
  sub_seurat <- NormalizeData(sub_seurat)
  sub_seurat <- sub_seurat %>% NormalizeData %>% SCTransform() %>%  RunPCA() %>% FindNeighbors(dims = 1:npcs) %>% RunUMAP(dims = 1:npcs)
  sub_seurat_list[[TB_id]] <- sub_seurat
  drug_lfc <- FoldChange(object = sub_seurat_list[[TB_id]], 
                         group.by = "Drug_Tx",
                         ident.1 = drug1,
                         ident.2 = "Ctrl",
                         assay = "RNA", 
                         slot = "data")
  drug_lfc <- drug_lfc[order(drug_lfc$avg_log2FC,decreasing = T),]
  lfc_lists[[TB_id]] <- drug_lfc
}
#create heatmaps#
dev.off()
agg_data <- data.frame(unified_gene_list = unique(unlist(sapply(diffex_matrices, FUN = function(x){x$X}))))
for(i in 1:length(diffex_matrices)){
  lfc_col <- sapply(agg_data$unified_gene_list, FUN = function(x) {
    val = lfc_lists[[i]]$avg_log2FC[rownames(lfc_lists[[i]]) == x]
    return(val)
  })
  p_val_col <- sapply(agg_data$unified_gene_list, FUN = function(x) {
    val = diffex_matrices[[i]]$p_val_adj[diffex_matrices[[i]]$X == x]
    if(identical(val, numeric(0))){
      return(1)
    } else{
      return(val)
    }
  })
  tmp_frame <- data.frame(as.numeric(lfc_col), as.numeric(p_val_col))
  colnames(tmp_frame) <- c(paste0(names(diffex_matrices)[i], "_lfc"), paste0(names(diffex_matrices)[i], "_p_val"))
  agg_data <- cbind(agg_data, tmp_frame)
}
agg_data <- agg_data[,c(1, grep("lfc", colnames(agg_data)), grep("p_val", colnames(agg_data)))]
agg_data$lfc_avg <- rowMeans(replace(agg_data[,grep("lfc", colnames(agg_data))], agg_data[,grep("lfc", colnames(agg_data))] == 0, NA), na.rm = TRUE)
agg_data$num_sig <- apply(agg_data[,grep("p_val", colnames(agg_data))], MARGIN = 1, FUN = function(x){sum(x < 0.05)})
agg_data <- agg_data[order(abs(agg_data$lfc_avg)),]
row.names(agg_data) <- agg_data$unified_gene_list

#add aggregated testing#
sub_seurat <- subset(GBM.sct.combined, idents = drug_TB_IDs)
sub_seurat <- NormalizeData(sub_seurat)
sub_seurat <- sub_seurat %>% NormalizeData %>% SCTransform() %>%  RunPCA() %>% FindNeighbors(dims = 1:npcs) %>% RunUMAP(dims = 1:npcs)
drug_lfc <- FoldChange(object = sub_seurat, 
                       group.by = "Drug_Tx",
                       ident.1 = drug1,
                       ident.2 = "Ctrl",
                       assay = "RNA", 
                       slot = "data")
drug_lfc <- drug_lfc[order(drug_lfc$avg_log2FC,decreasing = T),]
diffex_markers <- FindMarkers(object = sub_seurat, test.use = "MAST", ident.1 = drug1, ident.2 = "Ctrl", min.pct = 0.05, 
                              only.pos = F, assay = "RNA", slot = "data", logfc.threshold = 0.25, group.by = "Drug_Tx")
#add these into the agg_data frame#
lfc_col <- sapply(agg_data$unified_gene_list, FUN = function(x) {
  val = drug_lfc$avg_log2FC[rownames(drug_lfc) == x]
  return(val)
})
p_val_col <- sapply(agg_data$unified_gene_list, FUN = function(x) {
  val = diffex_markers$p_val_adj[rownames(diffex_markers) == x]
  if(identical(val, numeric(0))){
    return(1)
  } else{
    return(val)
  }
})
agg_data$merged_lfc <- lfc_col
agg_data$merged_p_val <- p_val_col

#alternate plotting of a cluster (i.e. cluster 10)#
for(i in 1:length(marker_list)){
  sub_agg_dat <- agg_data[row.names(agg_data) %isetwdn% marker_list[[i]],]
  # row_anno <- data.frame(sub_agg_dat$num_sig, row.names = row.names(sub_agg_dat))
  # colnames(row_anno) <- "p<0.05"
  anno_mat <- matrix("", nrow = dim(sub_agg_dat)[1], ncol = length(unique(drug_TB_IDs)) + 1)
  anno_mat[sub_agg_dat[,grep("_p_", colnames(agg_data))] < 0.05] <- "*"
  
  png(paste0("../Data/Wenting_slice_myeloid/comparative_anno/", trimws(drug1), "/", gsub(x = names(marker_list)[i], pattern = "/", replacement = "."), "_lfc_heatmap_custom.png"), 
      units = "mm", width = 475, height = 200, res = 300)
  myBreaks <- c(seq(min(sub_agg_dat[,grep("_lfc", colnames(agg_data))]), 0, length.out=ceiling(100/2) + 1), 
                seq(max(sub_agg_dat[,grep("_lfc", colnames(agg_data))])/100, max(sub_agg_dat[,grep("_lfc", colnames(agg_data))]), length.out=floor(100/2)))
  pheatmap(t(sub_agg_dat[,grep("_lfc", colnames(agg_data))]), color = rwb(100), border_color = "black", cluster_rows = F,
           labels_row = sapply(strsplit(colnames(sub_agg_dat)[grep("_lfc", colnames(sub_agg_dat))], "_"), FUN = "[", 1), breaks = myBreaks, scale = "none", 
           display_numbers = t(anno_mat), number_color = "black", fontsize_number = 14, cellwidth = 12, cellheight = 12)
  dev.off()
}

#-------------------------------aggregated GSEA---------------------------------------------------#
#GSEA on the aggregated dataset#
dir.create(paste0("../Data/Wenting_slice_myeloid/comparative_anno/", trimws(drug1), "/GSEA"))

drug_lfc <- drug_lfc[order(drug_lfc$avg_log2FC,decreasing = T),]
gene_list <- drug_lfc$avg_log2FC
names(gene_list) <- rownames(drug_lfc)
gsea_result <- GSEA(geneList = gene_list,
                    exponent = 1,
                    minGSSize = 10,
                    maxGSSize = 500,
                    eps = 0,
                    pvalueCutoff = 0.5,
                    pAdjustMethod = "BH",
                    TERM2GENE = module_table,
                    verbose = T,
                    by = "fgsea"
)

gsea_lists[["merge"]] <- gsea_result


for(j in unique(gsea_result@result$ID)){
  print(gseaplot(gsea_result, geneSetID = j, title = paste0("Cluster ", j, " signature"), color.line = uglia_colors[j]))
  ggsave(filename = paste0("../Data/Wenting_slice_myeloid/comparative_anno/", trimws(drug1), "/GSEA/", gsub(x = j, pattern = "/", replacement = "."), "_plot.png"), 
         units = "mm", width = 200, height = 350, dpi = 300)
}

if(length(unique(gsea_result@result$ID)) > 0){
  tmp <- gseaplot2(gsea_result, geneSetID = unique(gsea_result@result$ID), subplots = c(1),
                   title = paste0(drug1,": all signatures"), color = uglia_colors[unique(gsea_result@result$ID)], pvalue_table = F, base_size = 20) + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  tmp2 <- gseaplot2(gsea_result, geneSetID = unique(gsea_result@result$ID), subplots = c(2), color = uglia_colors[unique(gsea_result@result$ID)], base_size = 20) + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  tmp/tmp2 + plot_layout(heights = c(4,1))
  ggsave(filename = paste0("../Data/Wenting_slice_myeloid/comparative_anno/", trimws(drug1), "/GSEA/all_GSEA.png"), 
         units = "mm", width = 200, height = 300, dpi = 300)
  #manually add a  p-value table back in#
  pd <- gsea_result[unique(gsea_result@result$ID), c("Description", "pvalue", "p.adjust")]
  rownames(pd) <- pd$Description
  pd <- pd[,-1]
  pd <- round(pd, 4)
  tp <- ggtable(pd, tmp)
  print(tp)
  ggsave(filename = paste0("../Data/Wenting_slice_myeloid/comparative_anno/", trimws(drug1), "/GSEA/all_table.png"), 
         units = "mm", width = 50, height = 50, dpi = 300)
}

#make heatmap from aggregated gsea lists#
enrich_results <- data.table()
p_results <- data.table()
for(sample in c(drug_TB_IDs, "merge")){
  gsea_result <- gsea_lists[[sample]]@result
  gsea_result <- gsea_result[order(gsea_result$ID),]
  enrich_results <- cbind(enrich_results, gsea_result$enrichmentScore)
  p_results <- cbind(p_results, gsea_result$p.adjust)
  # enrich_tab <- cbind(gsea_result$ID, sample, gsea_result$enrichmentScore, gsea_result$p.adjust)
  # enrich_results <- rbind(enrich_results, enrich_tab)
}
colnames(enrich_results) <- c(drug_TB_IDs, "merge"); colnames(p_results) <- c(drug_TB_IDs, "merge")
rownames(enrich_results) <- gsea_result$ID; rownames(p_results) <- gsea_result$ID

#pivot the table for heatmap plotting
# enrich_results <- t(enrich_results)
# colnames(enrich_results) <- gsea_result$ID
anno_mat <- matrix("", nrow = dim(enrich_results)[1], ncol = length(unique(drug_TB_IDs)) + 1)
anno_mat[p_results<0.05] <- "*"; anno_mat[p_results < 0.01] <- "**"; anno_mat[p_results < 0.001] <- "***"
rownames(anno_mat) <- gsea_result$ID;colnames(anno_mat) <- c(drug_TB_IDs, "merge")

png(paste0("../Data/Wenting_slice_myeloid/comparative_anno/", trimws(drug1), "/GSEA/", "GSEA_heatmap.png"), units = "mm", width = 450, height = 500, res = 300)
pheatmap(enrich_results, scale = "none", color = rwb(100), border_color = "black", cluster_cols = F, cluster_rows = T, cellwidth = 50, cellheight = 50, 
         display_numbers = anno_mat, number_color = "black", fontsize_number = 27)
dev.off()