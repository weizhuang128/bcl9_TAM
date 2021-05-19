

#library(SingleR)
library(Seurat)
library(dplyr)
library(reshape2)
library(tidyr)
library(randomcoloR)
library(RColorBrewer)
library(monocle)



############Start here

####################################

setwd("E:/BaiduNetdiskDownload/Manuscript of scRNAseq BCL9 KD in mouse TIME/revised/TAM")

hsBCL9_pbmc_part<-readRDS("BCL9_scRNA_TAM.rds")


#hsBCL9_pbmc_part<- hsBCL9_pbmc_part[,which((hsBCL9_pbmc_part$seurat_clusters %in% c(0:9)))]
#hsBCL9_pbmc_part<- hsBCL9_pbmc_part[,which(!(hsBCL9_pbmc_part$seurat_clusters %in% c(6:8)))]

hsBCL9_pbmc_part <- ScaleData(hsBCL9_pbmc_part)
hsBCL9_pbmc_part <- ScaleData(hsBCL9_pbmc_part, vars.to.regress = c("nFeature_RNA","nCount_RNA", "percent.mt"), features = rownames(hsBCL9_pbmc_part))
hsBCL9_pbmc_part <- FindVariableFeatures(hsBCL9_pbmc_part, selection.method = "vst", nfeatures = 2000)
#hsBCL9_pbmc_part <- NormalizeData(hsBCL9_pbmc_part, normalization.method = "LogNormalize", scale.factor = 10000)

#####


# Identify the 10 most highly variable genes


# plot variable features with and without labels
#all.genes <- rownames(hsBCL9_pbmc_part)
#hsBCL9_pbmc_part <- ScaleData(hsBCL9_pbmc_part, features = rownames(hsBCL9_pbmc_part),block.size = 100,min.cells.to.block = 1000)

hsBCL9_pbmc_part <- RunPCA(hsBCL9_pbmc_part)
#VizDimLoadings(hsBCL9_pbmc_part, dims = 1:20, reduction = "pca")

ElbowPlot(object = hsBCL9_pbmc_part)

print(DimPlot(hsBCL9_pbmc_part, reduction = "pca"))

setwd("E:/BaiduNetdiskDownload/Manuscript of scRNAseq BCL9 KD in mouse TIME/revised/TAM")


pdf(paste0("TAM_hsBCL9_pbmc_part_DimHeatmap_",Sys.Date(),".pdf"))

DimHeatmap(hsBCL9_pbmc_part, dims = 1, cells = 500, balanced = TRUE)

dev.off()





hsBCL9_pbmc_part <- RunUMAP(hsBCL9_pbmc_part, dims = 1:20)#20

hsBCL9_pbmc_part <- RunTSNE(hsBCL9_pbmc_part, dims = 1:20)#20


hsBCL9_pbmc_part <- FindNeighbors(hsBCL9_pbmc_part, dims = 1:5)
hsBCL9_pbmc_part <- FindClusters(hsBCL9_pbmc_part, resolution = 0.2)#0.5





palette_seurat_clusters<-distinctColorPalette(length(unique(hsBCL9_pbmc_part$seurat_clusters)))
#palette_seurat_clusters <- brewer.pal( length(unique(hsBCL9_pbmc_part$seurat_clusters)),name ="Set3")

pdf(paste0("TAM_hsBCL9_pbmc_part_umap_cluster_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part, reduction = "umap",label = TRUE ,cols = palette_seurat_clusters ))
dev.off()

pdf(paste0("TAM_hsBCL9_pbmc_part_tsne_cluster_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part, reduction = "tsne",label = TRUE ,cols = palette_seurat_clusters))
dev.off()

palette_Sample <- brewer.pal( length(unique(hsBCL9_pbmc_part$Sample)),name ="Set3")
pdf(paste0("TAM_hsBCL9_pbmc_part_umap_Sample_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part, reduction = "umap", group.by = "Sample",label = TRUE,cols = palette_Sample ))
dev.off()

pdf(paste0("TAM_hsBCL9_pbmc_part_tsne_Sample_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part, reduction = "tsne", group.by = "Sample",label = TRUE ,cols = palette_Sample))
dev.off()

palette_Group <- brewer.pal( length(unique(hsBCL9_pbmc_part$Group)),name ="Dark2")
pdf(paste0("TAM_hsBCL9_pbmc_part_umap_Group_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part, reduction = "umap", group.by = "Group",cols = palette_Group ))
dev.off()
pdf(paste0("TAM_hsBCL9_pbmc_part_tsne_Group_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part, reduction = "tsne", group.by = "Group" ,cols = palette_Group))
dev.off()


palette_group2and3 <- brewer.pal( length(unique(hsBCL9_pbmc_part$Group)),name ="Accent")
pdf(paste0("TAM_hsBCL9_pbmc_part_umap_group2and3_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part, reduction = "umap", group.by = "group2and3",cols = palette_group2and3))
dev.off()

pdf(paste0("TAM_hsBCL9_pbmc_part_tsne_group2and3_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part, reduction = "tsne", group.by = "group2and3",cols =palette_group2and3 ))
dev.off()





endo_list <- c("Cxcl2","Cxcl10","Il1b","Ccl8","Cbr2","Folr2","C1qa","Cd74","Apoe")

endo_list<-intersect(rownames(hsBCL9_pbmc_part),endo_list)

for (plotgene in endo_list)
  
{
  
  #plotgene<-"H2???Eb1"
  pdf(paste0("TAM_hsBCL9_pbmc_partM1M2_partgenes_umap_",plotgene,Sys.Date(),".pdf"))
  
  print(FeaturePlot(hsBCL9_pbmc_part, features = plotgene,cols = c("grey","red") ))
  dev.off()
  pdf(paste0("TAM_hsBCL9_pbmc_partM1M2_partgenes_tsne_",plotgene,Sys.Date(),".pdf"))
  
  print(FeaturePlot(hsBCL9_pbmc_part, reduction = "tsne",features = plotgene,cols = c("grey","red") ))
  dev.off()
  
}









pbmc.markers <- FindAllMarkers(hsBCL9_pbmc_part, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, num.cores = 64)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


pdf(paste0("TAM_hsBCL9_pbmc_part_Seurat_DoHeatmap_",Sys.Date(),".pdf"),10,40)
print(DoHeatmap(hsBCL9_pbmc_part, features = top10$gene) + NoLegend())
dev.off()


features<- pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


features_gene <- c("Cxcl2","Cxcl10","Il1b","Ccl8","Cbr2","Folr2","C1qa","Cd74","Apoe")
pdf(paste0("TAM_clustering_M1M2_Seurat_DotPlot_3",Sys.Date(),".pdf"))
print(DotPlot(hsBCL9_pbmc_part, features = features_gene) + RotatedAxis())

dev.off()




seurat_clusters_count<-data.frame(seurat_clusters=hsBCL9_pbmc_part$seurat_clusters,Group=hsBCL9_pbmc_part$Group)
seurat_clusters_count<-t(table(seurat_clusters_count))
barpt<-ggplot(melt(seurat_clusters_count), aes(fill=Group, y=value, x=as.character(seurat_clusters))) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = palette_Group)+
  geom_text(aes(label=value),stat='identity',position=position_fill(vjust=0.5))+
  coord_flip()

pdf(paste0("TAM_seurat_clusters_barptM1M2_",Sys.Date(),".pdf"))
barpt
dev.off()





patient_name <-"TAM"



scRNA <- hsBCL9_pbmc_part

library(SingleR)
dir.create("CellType")
refdata <- MonacoImmuneData()
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters

cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main,
                    method = "cluster", clusters = clusters,
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"CellType/celltype_Monaco.csv",row.names = F)
scRNA@meta.data$celltype_Monaco = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_Monaco'] <- celltype$celltype[i]}

pdf(file=paste0(patient_name,Sys.Date(),"BCell_cell_Type_TSNE.pdf"))

DimPlot(scRNA, group.by="celltype_Monaco", repel=T, label=T, label.size=5, reduction='tsne')
dev.off()

pdf(file=paste0(patient_name,Sys.Date(),"BCell_cell_Type_umap.pdf"))

DimPlot(scRNA, group.by="celltype_Monaco", repel=T, label=T, label.size=5, reduction='umap')
dev.off()




refdata <- DatabaseImmuneCellExpressionData()
# load('~/database/SingleR_ref/ref_DICE_1561s.RData')
# refdata <- ref_DICE
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main,
                    method = "cluster", clusters = clusters,
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"CellType/celltype_DICE.csv",row.names = F)
scRNA@meta.data$celltype_DICE = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_DICE'] <- celltype$celltype[i]}

pdf(file=paste0(patient_name,Sys.Date(),"BCell_cell_Type_TSNE_DICE.pdf"))
DimPlot(scRNA, group.by="celltype_DICE", repel=T, label=T, label.size=5, reduction='tsne')
dev.off()
pdf(file=paste0(patient_name,Sys.Date(),"BCell_cell_Type_umap_DICE.pdf"))
DimPlot(scRNA, group.by="celltype_DICE", repel=T, label=T, label.size=5, reduction='umap')
dev.off()


hsBCL9_pbmc_part@meta.data[["Perturbed"]] <-hsBCL9_pbmc_part$Group %in% c("veh","NT")




#make GSEA table
##################


hsBCL9_pbmc_partMX<- hsBCL9_pbmc_part[,which(hsBCL9_pbmc_part$seurat_clusters ==4)]

matrix_temp_scale <- as.matrix(hsBCL9_pbmc_partMX@assays[["RNA"]]@data)  #scale.data

Group <- hsBCL9_pbmc_partMX$Perturbed


#convert mouse to human
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=F)
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(genesV2)
}



genegenegene<-convertMouseGeneList(rownames(matrix_temp_scale))

human_gene_list<-genegenegene[match(rownames(matrix_temp_scale),genegenegene$MGI.symbol),2]
gene_name_list <-human_gene_list



sink(paste0(patient_name,"_High_BCL9_TAMH1H2_Score_GSEA_culster4.cls"))
cat(as.numeric(length(colnames(matrix_temp_scale))),as.numeric(length(unique(as.character(Group)))),"1",sep = "\t","\n")

cat("#",unique(as.character(Group)),sep = "\t","\n")

cat(as.character(Group),sep = "\t","\n")

sink()


rm(pbmc)

ready_subset<-cbind(NAME=c(gene_name_list),Description="na",matrix_temp_scale)
ready_subset[1:10,1:10]


sink(paste0(patient_name,"_High_BCL9_endo_Score_GSEA_culster4.gct"))

cat("#1.2",sep = "\t","\n")
cat(dim(matrix_temp_scale),sep = "\t","\n")
write.table(ready_subset,file=paste0(patient_name,"_High_BCL9_endo_Score_GSEA_culster4.gct"),quote = F,sep = "\t",col.names=T,row.names = F,append = T)
sink()




#plot GSEA
#########
REST<-read.delim("E:/BaiduNetdiskDownload/Manuscript of scRNAseq BCL9 KD in mouse TIME/CT26 tumor shBCL9 scRNAseq results_out3/my_analysis.Gsea.1597976157032/gsea_report_for_FALSE_1597976157032.tsv",header = T)
test_group<-read.delim("E:/BaiduNetdiskDownload/Manuscript of scRNAseq BCL9 KD in mouse TIME/CT26 tumor shBCL9 scRNAseq results_out3/my_analysis.Gsea.1597976157032/gsea_report_for_TRUE_1597976157032.tsv",header = T)




REST_20<-subset(REST, GS.DETAILS!="")
REST_20$Group<-c("FALSE")
test_group_20<-subset(test_group, GS.DETAILS!="")
test_group_20$Group<-c("TRUE")
All<-rbind(REST_20, test_group_20)

bar_all_KEGG<-ggplot(All, aes(reorder(NAME, NES), NES, fill=Group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c(palette_group2and3))+#
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA (Endothelial cells)")

pdf(file="E:/BaiduNetdiskDownload/Manuscript of scRNAseq BCL9 KD in mouse TIME/CT26 tumor shBCL9 scRNAseq results_out3/Pathways_NES_from_GSEA.pdf",width=10, height=10)
bar_all_KEGG
dev.off()




###########################################





hsBCL9_pbmc_part0_2_4_temp<- hsBCL9_pbmc_part[,which(hsBCL9_pbmc_part$seurat_clusters %in% c(0,2,4))]

matrix_temp <- hsBCL9_pbmc_part0_2_4_temp@assays[["RNA"]]@counts
matrix_temp_scale <- hsBCL9_pbmc_part0_2_4_temp@assays[["RNA"]]@data

gc()

HSMM <- newCellDataSet(as.matrix(matrix_temp))

HSMM@phenoData@data[["Sample"]] <- hsBCL9_pbmc_part0_2_4_temp@meta.data[["Sample"]]
HSMM@phenoData@data[["Group"]] <- hsBCL9_pbmc_part0_2_4_temp@meta.data[["Group"]]
HSMM@phenoData@data[["Perturbed"]] <- hsBCL9_pbmc_part0_2_4_temp@meta.data[["Perturbed"]]



#rpc_matrix <- relative2abs(HSMM, method = "num_genes")
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
length(expressed_genes)
print(head(pData(HSMM)))

HSMM <- detectGenes(HSMM, min_expr = 0.1)

L <- log(exprs(HSMM[expressed_genes,]))

melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
gc()
qplot(value, geom="density", data=melted_dens_df) +  stat_function(fun = dnorm, size=0.5, color='red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")




disp_table <- dispersionTable(HSMM)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.01 &
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id


HSMM <- setOrderingFilter(HSMM, ordering_genes)
print(plot_ordering_genes(HSMM))
print(plot_pc_variance_explained(HSMM, return_all = F)) # norm_method = 'log',

gc()
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 20,
                        reduction_method = 'tSNE', verbose = T)


HSMM <- clusterCells(HSMM, num_clusters=2)


print(plot_cell_clusters(HSMM))

gc()
HSMM_myo <- estimateDispersions(HSMM)

disp_table <- dispersionTable(HSMM_myo)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.01 &
                           dispersion_empirical >=1 * dispersion_fit)$gene_id


HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
print(plot_ordering_genes(HSMM_myo))

gc()
HSMM_myo <- reduceDimension(HSMM_myo, max_components=2)
HSMM_myo <- orderCells(HSMM_myo)


patient_name<-"BCL9" 
LSEC_population <- "TAM"
gc()


pdf(file=paste0("hsBCL9_pbmc_partM1M20_2_4",patient_name,LSEC_population,Sys.Date(),"_singlemonocleplot.pdf"))
print(plot_cell_trajectory(HSMM_myo))
dev.off()
pdf(file=paste0("hsBCL9_pbmc_partM1M20_2_4",patient_name,LSEC_population,Sys.Date(),"_singlemonocleplotPseudotime.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by="Pseudotime"))
dev.off()

pdf(file=paste0("hsBCL9_pbmc_partM1M20_2_4",patient_name,LSEC_population,Sys.Date(),"_singlemonocleplotSample.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "Sample")+ scale_color_manual(breaks = waiver(),values=palette_Sample))
dev.off()
pdf(file=paste0("hsBCL9_pbmc_partM1M20_2_4",patient_name,LSEC_population,Sys.Date(),"_singlemonocleplotGroup.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "Group")+ scale_color_manual(breaks = waiver(),values=palette_Group))
dev.off()
pdf(file=paste0("hsBCL9_pbmc_partM1M20_2_4",patient_name,LSEC_population,Sys.Date(),"_singlemonocleplotseurat_Step3_clusters.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "Step3_clusters")+ scale_color_manual(breaks = waiver(),values=palette_seurat_clusters))
dev.off()
pdf(file=paste0("hsBCL9_pbmc_partM1M20_2_4",patient_name,LSEC_population,Sys.Date(),"_singlemonocleplotseurat_Step1_clusters.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "Step1_clusters")+ scale_color_manual(breaks = waiver(),values=palette_group2and3))
dev.off()



#plot_gene <- "C1qb"
pdf(file=paste0("hsBCL9_pbmc_partM1M20_2_4",patient_name,LSEC_population,Sys.Date(),plot_gene,"_",patient_name,LSEC_population,Sys.Date(),"_plot_complex_cell_trajectorysinglemonocle_plotgene.pdf"))

print(plot_cell_trajectory(HSMM_myo, color_by =(matrix_temp_scale[plot_gene,]))+scale_colour_gradientn(colours = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b"),limits=c(min(matrix_temp_scale[plot_gene,]),max(matrix_temp_scale[plot_gene,])))+labs(title=paste("Cluster",LSEC_population,plot_gene,patient_name)))
dev.off()







expressed_genes <- row.names(subset(fData(HSMM_myo), num_cells_expressed >= 100))

diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)",cores=8)
sig_gene_names <- row.names(subset(diff_test_res, (qval < 10**-100 & use_for_ordering=="TRUE")))




pdf(file=paste0("hsBCL9_pbmc_partM1M20_2_4",LSEC_population,patient_name,"_finishImages.pdf"),10,20)

print(plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                              num_clusters = 6,
                              cores = 8,
                              show_rownames = T))

dev.off()



###############################################



hsBCL9_pbmc_part2_9_12_temp<- hsBCL9_pbmc_part[,which(hsBCL9_pbmc_part$seurat_clusters %in% c(2,9,12))]


pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_temp_umap_cluster_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part2_9_12_temp, reduction = "umap",label = TRUE ,cols = palette_seurat_clusters ))
dev.off()

pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_temp_tsne_cluster_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part2_9_12_temp, reduction = "tsne",label = TRUE ,cols = palette_seurat_clusters))
dev.off()



hsBCL9_pbmc_part2_9_12<-hsBCL9_pbmc_part2_9_12_temp


hsBCL9_pbmc_part2_9_12@meta.data[["Step1_clusters"]] <-hsBCL9_pbmc_part2_9_12_temp$seurat_clusters


hsBCL9_pbmc_part2_9_12 <- ScaleData(hsBCL9_pbmc_part2_9_12)

hsBCL9_pbmc_part2_9_12 <- FindVariableFeatures(hsBCL9_pbmc_part2_9_12, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(hsBCL9_pbmc_part2_9_12)
hsBCL9_pbmc_part2_9_12 <- ScaleData(hsBCL9_pbmc_part2_9_12, features = rownames(hsBCL9_pbmc_part2_9_12),block.size = 100,min.cells.to.block = 1000)

hsBCL9_pbmc_part2_9_12 <- RunPCA(hsBCL9_pbmc_part2_9_12)
#VizDimLoadings(hsBCL9_pbmc_part2_9_12, dims = 1:20, reduction = "pca")

ElbowPlot(object = hsBCL9_pbmc_part2_9_12)

print(DimPlot(hsBCL9_pbmc_part2_9_12, reduction = "pca"))

setwd("E:/BaiduNetdiskDownload/Manuscript of scRNAseq BCL9 KD in mouse TIME/revised/TAM")


pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_DimHeatmap_",Sys.Date(),".pdf"))

DimHeatmap(hsBCL9_pbmc_part2_9_12, dims = 1, cells = 500, balanced = TRUE)

dev.off()

hsBCL9_pbmc_part2_9_12 <- FindNeighbors(hsBCL9_pbmc_part2_9_12, dims = 1:10)



hsBCL9_pbmc_part2_9_12 <- RunUMAP(hsBCL9_pbmc_part2_9_12, dims = 1:20)

hsBCL9_pbmc_part2_9_12 <- RunTSNE(hsBCL9_pbmc_part2_9_12, dims = 1:20)

hsBCL9_pbmc_part2_9_12 <- FindClusters(hsBCL9_pbmc_part2_9_12, resolution = 0.5)


palette_seurat_clusters<-distinctColorPalette(length(unique(hsBCL9_pbmc_part2_9_12$seurat_clusters)))
#palette_seurat_clusters <- brewer.pal( length(unique(hsBCL9_pbmc_part2_9_12$seurat_clusters)),name ="Set3")

pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_umap_cluster_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part2_9_12, reduction = "umap",label = TRUE ,cols = palette_seurat_clusters ))
dev.off()

pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_tsne_cluster_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part2_9_12, reduction = "tsne",label = TRUE ,cols = palette_seurat_clusters))
dev.off()

palette_Sample <- brewer.pal( length(unique(hsBCL9_pbmc_part2_9_12$Sample)),name ="Set3")
pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_umap_Sample_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part2_9_12, reduction = "umap", group.by = "Sample",label = TRUE,cols = palette_Sample ))
dev.off()

pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_tsne_Sample_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part2_9_12, reduction = "tsne", group.by = "Sample",label = TRUE ,cols = palette_Sample))
dev.off()

palette_Group <- brewer.pal( length(unique(hsBCL9_pbmc_part2_9_12$Group)),name ="Dark2")
pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_umap_Group_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part2_9_12, reduction = "umap", group.by = "Group",cols = palette_Group ))
dev.off()
pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_tsne_Group_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part2_9_12, reduction = "tsne", group.by = "Group" ,cols = palette_Group))
dev.off()


palette_group2and3 <- brewer.pal( length(unique(hsBCL9_pbmc_part2_9_12$Group)),name ="Accent")
pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_umap_Step1_clusters_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part2_9_12, reduction = "umap", group.by = "Step1_clusters",cols = palette_group2and3))
dev.off()

pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_tsne_Step1_clusters_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part2_9_12, reduction = "tsne", group.by = "Step1_clusters",cols =palette_group2and3 ))
dev.off()




hsBCL9_pbmc_part2_9_12_temp@meta.data[["Step3_clusters"]] <-hsBCL9_pbmc_part2_9_12$seurat_clusters



palette_seurat_clusters<-distinctColorPalette(length(unique(hsBCL9_pbmc_part2_9_12$seurat_clusters)))
#palette_seurat_clusters <- brewer.pal( length(unique(hsBCL9_pbmc_part2_9_12$seurat_clusters)),name ="Set3")

pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_temp_umap_Step1_clusters_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part2_9_12_temp, reduction = "umap", group.by = "Step3_clusters",label = TRUE ,cols = palette_seurat_clusters ))
dev.off()

pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_temp_tsne_Step1_clusters_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part2_9_12_temp, reduction = "tsne", group.by = "Step3_clusters",label = TRUE ,cols = palette_seurat_clusters))
dev.off()



seurat_clusters_count<-data.frame(seurat_clusters=hsBCL9_pbmc_part2_9_12$seurat_clusters,Group=hsBCL9_pbmc_part2_9_12$Group)
seurat_clusters_count<-t(table(seurat_clusters_count))
barpt<-ggplot(melt(seurat_clusters_count), aes(fill=Group, y=value, x=as.character(seurat_clusters))) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = palette_Group)+
  geom_text(aes(label=value),stat='identity',position=position_fill(vjust=0.5))+
  coord_flip()

pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_seurat_clusters_barpt_",Sys.Date(),".pdf"))
barpt
dev.off()



pbmc.markers <- FindAllMarkers(hsBCL9_pbmc_part2_9_12, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, num.cores = 64)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)





pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_Seurat_DoHeatmap_",Sys.Date(),".pdf"),10,40)
print(DoHeatmap(hsBCL9_pbmc_part2_9_12, features = top10$gene) + NoLegend())
dev.off()


endo_list<-c("Hbegf","Cav2","Phgdh","Esm1","Hspa5","Lars2","Gm26917","Cebpd","Serping1","Id3","Ogn","Serpinf1","Mgp","Dcn","Cp","Tmem176b","Aspn","Gzmc","Gzma","Ccl5","Cd52","AW112010","Cd3e","Ltb","Gzmb","Il2rb","Nkg7","Ubc","Rad21","H1f0","Col1a1","Hist1h1e","Col3a1","Hist1h1b","Gm42418","AY036118","Hist1h2af","Cxcl10","Ly6a","Il1b","Ifit1","Dusp1","Isg15","Irf7","Rsad2","Ly6c2","Plac8","Ccl8","Wfdc17","Pf4","Ccl6","Ccl12","Lyz2","Apoe","C1qb","C1qc","C1qa","Nppb","Krt20","Qrfp","2200002D01Rik","Lxn","Igfbp6","Prkg2","Timp1","Hmga1","Lgals7","Tuba1c","Cdc20","Cenpa","Cenpf","Cks2","Tubb4b","Top2a","Birc5","Ube2c","Hmgb2")


endo_list<-intersect(rownames(hsBCL9_pbmc_part2_9_12),endo_list)

for (plotgene in endo_list)
  
{
  
  #plotgene<-"H2???Eb1"
  pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_partgenes_umap_",plotgene,Sys.Date(),".pdf"))
  
  print(FeaturePlot(hsBCL9_pbmc_part2_9_12, features = plotgene,cols = c("grey","red") ))
  dev.off()
  pdf(paste0("TAM_hsBCL9_pbmc_part2_9_12_partgenes_tsne_",plotgene,Sys.Date(),".pdf"))
  
  print(FeaturePlot(hsBCL9_pbmc_part2_9_12, reduction = "tsne",features = plotgene,cols = c("grey","red") ))
  dev.off()
  
}




#hsBCL9_pbmc_part4_5_7_temp<- hsBCL9_pbmc_part2_9_12[,which(hsBCL9_pbmc_part2_9_12$seurat_clusters %in% c(4,5,7))]

hsBCL9_pbmc_part4_5_7_temp<- hsBCL9_pbmc_part2_9_12
matrix_temp <- hsBCL9_pbmc_part4_5_7_temp@assays[["RNA"]]@counts
matrix_temp_scale <- hsBCL9_pbmc_part4_5_7_temp@assays[["RNA"]]@data

gc()

HSMM <- newCellDataSet(as.matrix(matrix_temp))

HSMM@phenoData@data[["Sample"]] <- hsBCL9_pbmc_part4_5_7_temp@meta.data[["Sample"]]
HSMM@phenoData@data[["Group"]] <- hsBCL9_pbmc_part4_5_7_temp@meta.data[["Group"]]
HSMM@phenoData@data[["Step3_clusters"]] <- hsBCL9_pbmc_part4_5_7_temp@meta.data[["seurat_clusters"]]
HSMM@phenoData@data[["Step1_clusters"]] <- hsBCL9_pbmc_part4_5_7_temp@meta.data[["Step1_clusters"]]


#rpc_matrix <- relative2abs(HSMM, method = "num_genes")
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
length(expressed_genes)
print(head(pData(HSMM)))

HSMM <- detectGenes(HSMM, min_expr = 0.1)

L <- log(exprs(HSMM[expressed_genes,]))

melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
gc()
qplot(value, geom="density", data=melted_dens_df) +  stat_function(fun = dnorm, size=0.5, color='red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")




disp_table <- dispersionTable(HSMM)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.01 &
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id


HSMM <- setOrderingFilter(HSMM, ordering_genes)
print(plot_ordering_genes(HSMM))
print(plot_pc_variance_explained(HSMM, return_all = F)) # norm_method = 'log',

gc()
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 20,
                        reduction_method = 'tSNE', verbose = T)


HSMM <- clusterCells(HSMM, num_clusters=2)


print(plot_cell_clusters(HSMM))

gc()
HSMM_myo <- estimateDispersions(HSMM)

disp_table <- dispersionTable(HSMM_myo)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.01 &
                           dispersion_empirical >=1 * dispersion_fit)$gene_id


HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
print(plot_ordering_genes(HSMM_myo))

gc()
HSMM_myo <- reduceDimension(HSMM_myo, max_components=2)
HSMM_myo <- orderCells(HSMM_myo)


patient_name<-"BCL9" 
LSEC_population <- "TAM"
gc()


pdf(file=paste0("hsBCL9_pbmc_part2_9_12",patient_name,LSEC_population,Sys.Date(),"_singlemonocleplot.pdf"))
print(plot_cell_trajectory(HSMM_myo))
dev.off()
pdf(file=paste0("hsBCL9_pbmc_part2_9_12",patient_name,LSEC_population,Sys.Date(),"_singlemonocleplotPseudotime.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by="Pseudotime"))
dev.off()

pdf(file=paste0("hsBCL9_pbmc_part2_9_12",patient_name,LSEC_population,Sys.Date(),"_singlemonocleplotSample.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "Sample")+ scale_color_manual(breaks = waiver(),values=palette_Sample))
dev.off()
pdf(file=paste0("hsBCL9_pbmc_part2_9_12",patient_name,LSEC_population,Sys.Date(),"_singlemonocleplotGroup.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "Group")+ scale_color_manual(breaks = waiver(),values=palette_Group))
dev.off()
pdf(file=paste0("hsBCL9_pbmc_part2_9_12",patient_name,LSEC_population,Sys.Date(),"_singlemonocleplotseurat_Step3_clusters.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "Step3_clusters")+ scale_color_manual(breaks = waiver(),values=palette_seurat_clusters))
dev.off()
pdf(file=paste0("hsBCL9_pbmc_part2_9_12",patient_name,LSEC_population,Sys.Date(),"_singlemonocleplotseurat_Step1_clusters.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "Step1_clusters")+ scale_color_manual(breaks = waiver(),values=palette_group2and3))
dev.off()



#plot_gene <- "C1qb"
pdf(file=paste0("hsBCL9_pbmc_part2_9_12",patient_name,LSEC_population,Sys.Date(),plot_gene,"_",patient_name,LSEC_population,Sys.Date(),"_plot_complex_cell_trajectorysinglemonocle_plotgene.pdf"))

print(plot_cell_trajectory(HSMM_myo, color_by =(matrix_temp_scale[plot_gene,]))+scale_colour_gradientn(colours = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b"),limits=c(min(matrix_temp_scale[plot_gene,]),max(matrix_temp_scale[plot_gene,])))+labs(title=paste("Cluster",LSEC_population,plot_gene,patient_name)))
dev.off()



expressed_genes <- row.names(subset(fData(HSMM_myo), num_cells_expressed >= 100))

diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)",cores=8)
sig_gene_names <- row.names(subset(diff_test_res, (qval < 10**-100 & use_for_ordering=="TRUE")))




pdf(file=paste0("hsBCL9_pbmc_part2_9_12",LSEC_population,patient_name,"_finishImages.pdf"),10,20)

print(plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                              num_clusters = 6,
                              cores = 8,
                              show_rownames = T))

dev.off()


total_DF<-data.frame(ID=colnames(pbmc))

part_DF<-data.frame(ID=colnames(hsBCL9_pbmc_part2_9_12),Step3_clusters=hsBCL9_pbmc_part2_9_12$seurat_clusters,Step1_clusters=hsBCL9_pbmc_part2_9_12$Step1_clusters)


merged_part_DF<-merge(total_DF,part_DF,by="ID",all.x=T)

merged_part_DF<-merged_part_DF[match(colnames(pbmc),merged_part_DF$ID),]

pbmc@meta.data[["Step3_clusters"]]<-merged_part_DF$Step3_clusters
pbmc@meta.data[["Step1_clusters"]]<-merged_part_DF$Step1_clusters



pdf(paste0("TAM_allcells_Step3_clusters_",Sys.Date(),".pdf"))
print(DimPlot(pbmc, reduction = "umap", group.by = "Step3_clusters",cols = palette_Sample))
dev.off()
pdf(paste0("TAM_allcells_Step1_clusters_",Sys.Date(),".pdf"))
print(DimPlot(pbmc, reduction = "umap", group.by = "Step1_clusters",cols = palette_Sample))
dev.off()
