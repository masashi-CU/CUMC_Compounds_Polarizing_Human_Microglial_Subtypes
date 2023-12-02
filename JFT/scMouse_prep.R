###############################
#         jft2113             #
#         07/27/22            #
#         R4.1                #
###############################
#  analyze VCH mouse data     #

####-------------continued analysis on the cluster-------------------------------------####
library(Matrix)
library(Seurat)
library(demuxmix)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(viridis)
library(grid)
library(ggtree)
library(tidyverse)
library(forcats)
library(MAST)
library(DropletUtils)
library(harmony)
library(patchwork)
library(schex)
library(org.Mm.eg.db)
library(biomaRt)
library(RColorBrewer)
library(rstatix)
library(clusterProfiler)
library(enrichplot)

source("~/ChooseR/helper_functions.R")
source("~/ChooseR/pipeline.R")
source("~/full_uglia_v2.0/full_functions.R")
harmony_pipeline <- function(seurat, var_genes = 3000, PCs = 20, rez = c(0.35, 0.5, 0.75), harm_regress = c("orig.ident"), 
                             SCT_regress = c("nCount_RNA"), apply_SCT = F, lambda = 1, theta = 2){
  if(apply_SCT == T){
    train_harm <- SCTransform(seurat, verbose = FALSE, variable.features.n = var_genes, vars.to.regress = SCT_regress)
  }
  else{
    train_harm <- FindVariableFeatures(seurat, verbose = FALSE, nfeatures = var_genes)
    train_harm <- ScaleData(train_harm, verbose = FALSE)
  }
  train_harm <- RunPCA(train_harm, npcs = PCs, verbose = FALSE)
  train_harm <- RunHarmony(train_harm, group.by.vars = harm_regress, lambda = lambda, theta = theta)
  train_harm <- RunUMAP(train_harm, reduction = "harmony", dims = 1:PCs)
  train_harm <- FindNeighbors(train_harm, reduction = "harmony", dims = 1:PCs)
  train_harm <- FindClusters(train_harm, resolution = rez)
  return(train_harm)
}
schex_plots <- function(dataset, gene_list, set_lim = T){
  if(set_lim == T){
    gg_list = lapply(gene_list, FUN = function(x){
      y = plot_hexbin_gene(dataset, type = "data", gene = x, action = "mean", title = x) + xlim(-8.5, 7.5) + ylim(-5.25, 5) + 
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 15), axis.title = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 20))+ scale_fill_viridis_c(guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
    })
  }
  else{
    gg_list = lapply(gene_list, FUN = function(x){
      y = plot_hexbin_gene(dataset, type = "data", gene = x, action = "mean", title = x) + 
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 15), axis.title = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 20))+ scale_fill_viridis_c(guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
    })
  }
  return(gg_list)
}
convertHumanGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
plot_module_enrichment <- function(seurat, module_list, output_dir = "../Compound_work/HMC3 scRNAseq/Viz/Metabolism/") {
  colnames(module_list) <- c("Pathway", "Gene")
  #subset to those genes in the dataset#
  module_list = module_list[module_list$Gene %in% rownames(seurat),]
  #plot the various modules overlaid on the data#
  all_modules = unique(module_list$Pathway)
  for(module in all_modules){
    all_genes = list(module_list$Gene[module_list$Pathway == module])
    seurat = AddModuleScore(object = seurat, features = all_genes, name = gsub(pattern = " ", x = module, replacement =  "_"))
  }
  schex <- make_hexbin(seurat, nbins = 100, dimension_reduction = "UMAP")
  
  pdf(paste0(output_dir, "pathway_modules.pdf"), width = 11, height = 8.5)
  for(module in all_modules){
    col_nam <-
      gsub(pattern = "-", replacement = ".", x = module) %>%
      gsub(pattern = " ", replacement =  "_", .) %>%
      gsub(pattern = ",", replacement =  ".", .)
    
    print(plot_hexbin_meta(schex, col = paste0(col_nam, 1), action = "mean") +  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), 
                                                                                                     guide = guide_colorbar(label = TRUE,
                                                                                                                            draw.ulim = TRUE, 
                                                                                                                            draw.llim = TRUE,
                                                                                                                            frame.colour = "black", ticks = TRUE, 
                                                                                                                            label.position = "bottom",
                                                                                                                            barwidth = 18,
                                                                                                                            barheight = 1.3, 
                                                                                                                            direction = 'horizontal')) +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position = "bottom", text=element_text(size = 20),
                  plot.title = element_text(hjust = 0.5, size = 30)) + NoAxes() + ggtitle(module))
  }
  dev.off()
}

#starting with the compound single-cell data#
data_directory <- "Data/"
drug_files <- list.files(data_directory)[grep("PV", list.files(data_directory))]

#cols for compound treatment and other choices#
compound_cols <- c( "#EE7600", "gray30",  "#91CADB")
gender_cols <- c("#CE4969", "#5887C2")

#------------------------------------prep the datasets ---------------------------------------------# 
seuratmat <- list()
hashtag_frame <- data.frame(Matrix(nrow = 0, ncol = 27))
#filtering parameters#
umimin = 500; umimax = 1e10; HTO_filter = 0; RP_filt = T; MT_filt = T
#read in the HTO mappint#
HTO_mapping = data.frame(read.csv("HTO_mapping.csv"))
colnames(HTO_mapping) = c("Name", "Condition", paste0("Hashtag0",1:4))
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
    mtpct=colSums(cellmat[grep("^mt-",rownames(cellmat)),])/colSums(cellmat)*100
    med_cutoff = median(mtpct) + 2*mad(mtpct)
    IQR_cutoff = as.numeric(quantile(mtpct)[4] + IQR(mtpct)*1.5)
    mitofilter = min(med_cutoff, IQR_cutoff)
    mitofilter = max(mitofilter, 5)
    print(paste(file, " mito: med-", med_cutoff, ", IQR-", IQR_cutoff, sep = ""))
    cellmat=cellmat[,which(mtpct<mitofilter)]  ###remove mito-heavy cells
    hashmat=hashmat[,which(mtpct<mitofilter)]
  }
  ###filter genes unlikely to be biologically relevant out, then filter on UMI count
  if(RP_filt == T){
    cellmat=cellmat[grep("mt-|^Rp[0-9]|^Bc[0-9]|^Rpl|^Rps|^mtrnr|-ps|^Gm[0-9]+",rownames(cellmat),invert=T),] ###variant to remove all RP, MT, and pseudogenes
  }else if (MT_filt == T){
    cellmat=cellmat[grep("^mt-|^bc[0-9]|^mtrnr|-ps|^Gm[0-9]+",rownames(cellmat),invert=T),] ###variant to remove MT + PS
  } else{
    cellmat=cellmat[grep("^BC[0-9]|-PS",rownames(cellmat),invert=T),] ###variant to remove MT + PS
  }
  
  #what is the distribution of UMIs here#
  quantile(colSums(cellmat))

  
  keepcells=which(colSums(cellmat) > umimin & colSums(cellmat) < umimax & colSums(hashmat) >= HTO_filter)     ###Do we want this now or after hashtag deconvolution?###
  hashtag_frame = rbind(hashtag_frame, t(rowSums(hashmat[,keepcells])))
  rownames(hashtag_frame)[dim(hashtag_frame)[1]] = file
  
  #summed = rowsum(tempmat, row.names(tempmat))
  dat_seurat = CreateSeuratObject(counts=cellmat[,keepcells])
  #remove the other hashes not present in the sample#
  hashmat <- hashmat[rownames(hashmat) %in% colnames(HTO_mapping)[!is.na(HTO_mapping[1,])][-1],keepcells]
  dat_seurat[["HTO"]]=CreateAssayObject(counts=hashmat)
  tempcol=dat_seurat$orig.ident
  barcodes = names(tempcol)
  names(tempcol)=paste(file,"-",1:length(tempcol),sep="")
  dat_seurat = RenameCells(dat_seurat, new.names = names(tempcol))
  tempcol = rep(file, length(tempcol))
  dat_seurat$orig.ident = tempcol
  names(barcodes) = names(tempcol)
  dat_seurat = AddMetaData(dat_seurat, metadata = barcodes, col.name = "barcodes")
  seuratmat[[file]] = dat_seurat
}
rm(tempmat, cellmat, dat_seurat); gc()
load("mouse_demux_df.rda")
for(file in drug_files){
  seuratmat[[file]]@meta.data <- meta_info[[file]]
}
saveRDS(seuratmat, file = "combined_doublet_files.rds")

#--------------------------------------------------------load back in----------------------------------------------------------------------------------------#
seuratmat <- readRDS("combined_doublet_files.rds")
#add in the DF classifications#
for(file in names(seuratmat)){
  seuratmat[[file]]$DF_class <- seuratmat[[file]]@meta.data[,19]
}
#merge
combined_mouse_dat <- merge(seuratmat[[1]], seuratmat[-1])
rm(seuratmat); gc()
#look at DF vs. demuxmix
table(combined_mouse_dat$Type, combined_mouse_dat$DF_class)
combined_mouse_dat <- combined_mouse_dat[,(combined_mouse_dat$Type == "singleton" & combined_mouse_dat$DF_class == "Singlet")]
#add some metadata columns#
#compound
compound_col <- combined_mouse_dat$identity
compound_col[grep("DMSO", x = compound_col)] <- "DMSO"
compound_col[grep("Campto", x = compound_col)] <- "Camptothecin"
compound_col[grep("Narcic", x = compound_col)] <- "Narciclasine"
#gender
gender_col <- combined_mouse_dat$identity
gender_col[grep("F", x = gender_col)] <- "F"
gender_col[grep("M", x = gender_col)] <- "M"

combined_mouse_dat <- combined_mouse_dat %>% 
  AddMetaData(metadata = compound_col, col.name = "treatment") %>%
  AddMetaData(metadata = gender_col, col.name = "gender")
#run basic dim reduc pipeline#
set.seed(dim(combined_mouse_dat)[2])
combined_mouse_dat <- combined_mouse_dat %>% SCTransform() %>%  RunPCA()
ElbowPlot(combined_mouse_dat, ndims = 50)
combined_mouse_dat <- combined_mouse_dat %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = c(0.3, 0.5, 0.7, 1)) %>% RunUMAP(dims = 1:30) %>% 
  NormalizeData(assay = "RNA")

saveRDS(combined_mouse_dat, file = "combined_mouse_std.rds")
# combined_mouse_dat <- readRDS("combined_mouse_std.rds")

#plot basic visualization
DimPlot(combined_mouse_dat, pt.size = 1, label = F, group.by = "orig.ident", raster = F) + ggtitle("orig.ident") + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 40)) + NoAxes()
ggsave(paste0("Viz/combined_orig.png"), units = "mm", width = 450, height = 350, dpi = 300)
DimPlot(combined_mouse_dat, pt.size = 1, label = F, group.by = "identity", raster = F) + ggtitle("Identity") + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 40)) + NoAxes()
ggsave(paste0("Viz/combined_ident.png"), units = "mm", width = 650, height = 350, dpi = 300)
DimPlot(combined_mouse_dat, pt.size = 2, group.by = "treatment", cols = compound_cols, raster = F) + ggtitle("Treatment") + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 30)) + NoAxes()
ggsave("Viz/combined_treatment.png", units = "mm", width = 450, height = 350, dpi = 300)
DimPlot(combined_mouse_dat, pt.size = 2, group.by = "gender", cols = gender_cols, raster = F) + ggtitle("Gender") + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 30)) + NoAxes()
ggsave("Viz/combined_gender.png", units = "mm", width = 450, height = 350, dpi = 300)

#for example, what does this clustering structure look like#
combined_mouse_dat <- SetIdent(combined_mouse_dat, value = "SCT_snn_res.0.5")
new.ids = 1:length(unique(combined_mouse_dat$SCT_snn_res.0.5))
names(new.ids) <- 0:(length(unique(combined_mouse_dat$SCT_snn_res.0.5))-1)
combined_mouse_dat <- RenameIdents(combined_mouse_dat, new.ids)
DimPlot(combined_mouse_dat, pt.size = 1, label = F, raster = F) + ggtitle("Identity") + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), text=element_text(size = 40)) + NoAxes()
ggsave(paste0("Viz/combined_basic_clustering.png"), units = "mm", width = 550, height = 350, dpi = 300)

#characteristic marker genes#
classicgenes = c("Snap25", "Eno2", "Gad1", "Gad2", "Slc32a1", "Dlx5", "Lhx6", "Sst", "Pvalb", "Vip", "Htr3a",
                 "Ndnf", "Lamp5", "Slc17a6", "Slc17a7", "Enc1", "Cux2", "Cbln4", "Scnn1a", "Rorb", "Fezf2",
                 "Foxp2", "Ntsr1", "Ctgf", "Aqp4", "Fgfr3", "Gfap", "Olig1", "Olig2", "Opalin", "Pdgfra", 
                 "Aif1", "C1qa", "Siglech", "Ctss", "Kcnj8", "Acta2", "Fcn1", "Vcan", "Cd3e", "Gzmb", "Ccr7",
                 "Cd79a", "Cd8a")
DotPlot(combined_mouse_dat, assay = "RNA", features = classicgenes, dot.min = 0, cluster.idents = F, scale = T) + coord_flip() +
  scale_colour_viridis(option="magma", guide = guide_colorbar(frame.colour = "black", ticks = TRUE)) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.title = element_blank()) + scale_y_discrete(position = "right")
ggsave("Viz/celltype_markers_v1.png", units = "in", width = 7, height = 10, dpi = 300)

#schex for slightly clearer visualization#
full <- make_hexbin(combined_mouse_dat, nbins = 100, dimension_reduction = "UMAP")
png("Viz/schex_density.png", units = "mm", width = 350, height = 350, res = 300)
plot_hexbin_density(full)
dev.off()

#plot metadata#
full$assigned_idents = factor(Idents(combined_mouse_dat))
png("Viz/schex_clustering.png", units = "mm",  width = 350, height = 350, res = 300)
plot_hexbin_meta(full, col = "assigned_idents", action = "majority")
dev.off()
png("Viz/schex_orig.png", units = "mm",  width = 350, height = 350, res = 300)
plot_hexbin_meta(full, col = "orig.ident", action = "majority")
dev.off()
png("Viz/schex_ident.png", units = "mm",  width = 350, height = 350, res = 300)
plot_hexbin_meta(full, col = "treatment", action = "majority")
dev.off()
png("Viz/schex_gender.png", units = "mm",  width = 350, height = 350, res = 300)
plot_hexbin_meta(full, col = "gender", action = "majority")
dev.off()

#plot some marker genes on the schex#
final_features <- c("Snap25", "Gad1", "Pvalb", "Cux2", "Aqp4", "Gfap", "Olig1", "Pdgfra", "Aif1", "C1qa", "Siglech", "Ctss", "Acta2", "Vcan", "Cd3e", "Cd8a")
gg_markers = schex_plots(dataset = full, gene_list = final_features, set_lim = F)
wrap_plots(gg_markers, ncol = 4, nrow = 4)
ggsave("Viz/feature_schex.png", units = "mm",width = 550, height = 450, dpi = 300)

#-------------------------------------compute mitotypes on Natacha's objects------------------------------------------------------#
mito_type_list = read.csv("MouseMitoCarta_Gene_to_Pathway.csv")
colnames(mito_type_list) <- c("term", "gene")

f_mice <- readRDS("/mnt/mfs/hgrcgrid/homes/nc3018/For_John/meld_res_sobj_F.rds")
plot_module_enrichment(seurat = f_mice, module_list = mito_type_list, output_dir = "Viz/Metabolism/revised/F_Mito_")

m_mice <- readRDS("/mnt/mfs/hgrcgrid/homes/nc3018/For_John/meld_res_sobj_M.rds")
plot_module_enrichment(seurat = m_mice, module_list = mito_type_list, output_dir = "Viz/Metabolism/revised/M_Mito_")

mito_type_list = read.csv("Gene-to-pathway-list_hm.csv")
mito_type_list <- mito_type_list[,c(2,1)]
colnames(mito_type_list) <- c("term", "gene")
img <- readRDS("/mnt/mfs/hgrcgrid/homes/nc3018/For_John/meld_res_sobj_img.rds")
plot_module_enrichment(seurat = img, module_list = mito_type_list, output_dir = "Viz/Metabolism/revised/iMG_Mito_")
