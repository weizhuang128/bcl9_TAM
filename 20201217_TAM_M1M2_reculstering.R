


#library(SingleR)
library(Seurat)
library(dplyr)
library(reshape2)
library(tidyr)
library(randomcoloR)
library(RColorBrewer)
library(monocle)

################
setwd("E:/BaiduNetdiskDownload/Manuscript of scRNAseq BCL9 KD in mouse TIME/revised")

pbmc<-readRDS("bcl9_all_cells_data.rds") #reference https://doi.org/10.3389/fonc.2020.603702

hsBCL9_pbmc_part<- pbmc[,which(pbmc$seurat_clusters %in% c(3))]

hsBCL9_pbmc_part <- ScaleData(hsBCL9_pbmc_part)
hsBCL9_pbmc_part <- FindVariableFeatures(hsBCL9_pbmc_part, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(hsBCL9_pbmc_part)
hsBCL9_pbmc_part <- ScaleData(hsBCL9_pbmc_part, features = rownames(hsBCL9_pbmc_part),block.size = 100,min.cells.to.block = 1000)
hsBCL9_pbmc_part <- RunPCA(hsBCL9_pbmc_part)

hsBCL9_pbmc_part <- RunUMAP(hsBCL9_pbmc_part, dims = 1:20)#20
hsBCL9_pbmc_part <- RunTSNE(hsBCL9_pbmc_part, dims = 1:20)#20
hsBCL9_pbmc_part <- FindNeighbors(hsBCL9_pbmc_part, dims = 1:10)
hsBCL9_pbmc_part <- FindClusters(hsBCL9_pbmc_part, resolution = 0.2)#0.5
hsBCL9_pbmc_part<- hsBCL9_pbmc_part[,which(!(hsBCL9_pbmc_part$seurat_clusters %in% c(6:8)))]


hsBCL9_pbmc_part <- ScaleData(hsBCL9_pbmc_part)
hsBCL9_pbmc_part <- ScaleData(hsBCL9_pbmc_part, vars.to.regress = c("nFeature_RNA","nCount_RNA", "percent.mt"), features = rownames(hsBCL9_pbmc_part))
hsBCL9_pbmc_part <- FindVariableFeatures(hsBCL9_pbmc_part, selection.method = "vst", nfeatures = 2000)


hsBCL9_pbmc_part <- RunPCA(hsBCL9_pbmc_part)



hsBCL9_pbmc_part <- RunUMAP(hsBCL9_pbmc_part, dims = 1:20)#20

hsBCL9_pbmc_part <- RunTSNE(hsBCL9_pbmc_part, dims = 1:20)#20


hsBCL9_pbmc_part <- FindNeighbors(hsBCL9_pbmc_part, dims = 1:5)
hsBCL9_pbmc_part <- FindClusters(hsBCL9_pbmc_part, resolution = 0.2)#0.5

hsBCL9_pbmc_part@meta.data[["Perturbed"]] <-hsBCL9_pbmc_part$Group %in% c("veh","NT")

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


palette_group2and3 <- brewer.pal( length(unique(hsBCL9_pbmc_part$Perturbed)),name ="Accent")
pdf(paste0("TAM_hsBCL9_pbmc_part_umap_Perturbed_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part, reduction = "umap", group.by = "Perturbed",cols = palette_group2and3))
dev.off()

pdf(paste0("TAM_hsBCL9_pbmc_part_tsne_Perturbed_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part, reduction = "tsne", group.by = "Perturbed",cols =palette_group2and3 ))
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





seurat_clusters_count<-data.frame(seurat_clusters=hsBCL9_pbmc_part$seurat_clusters,Perturbed=hsBCL9_pbmc_part$Perturbed)
seurat_clusters_count<-t(table(seurat_clusters_count))
barpt<-ggplot(melt(seurat_clusters_count), aes(fill=Perturbed, y=value, x=as.character(seurat_clusters))) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = palette_group2and3)+
  geom_text(aes(label=value),stat='identity',position=position_fill(vjust=0.5))+
  coord_flip()

pdf(paste0("TAM_seurat_clusters_barptM1M2Perturbed_",Sys.Date(),".pdf"))
barpt
dev.off()


rm(pbmc)



gc()

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
                         mean_expression >= 0.00001 &
                           dispersion_empirical >=1 * dispersion_fit)$gene_id


HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
print(plot_ordering_genes(HSMM_myo))

gc()
HSMM_myo <- reduceDimension(HSMM_myo, max_components=2)
HSMM_myo <- orderCells(HSMM_myo)
#saveRDS(hsBCL9_pbmc_part, file = "hsBCL9_pbmc_part_TAM_M1M2.rds")
#save(HSMM_myo,file = "hsBCL9_HSMM_myo_part_TAM_M1M2")

print(plot_cell_trajectory(HSMM_myo))




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
pdf(file=paste0("hsBCL9_pbmc_partM1M20_2_4",patient_name,LSEC_population,Sys.Date(),"_singlemonocleplotseurat_Perturbed.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "Perturbed")+ scale_color_manual(breaks = waiver(),values=palette_group2and3))
dev.off()


endo_list <- c("Cxcl2","Cxcl10","Il1b","Ccl8","Cbr2","Folr2","C1qa","Cd74","Apoe")


for (plot_gene in endo_list)
  
{

#plot_gene <- "Apoe"
pdf(file=paste0("hsBCL9_pbmc_partM1M20_2_4",patient_name,LSEC_population,Sys.Date(),plot_gene,"_",patient_name,LSEC_population,Sys.Date(),"_plot_complex_cell_trajectorysinglemonocle_plotgene.pdf"))

print(plot_cell_trajectory(HSMM_myo, color_by =(matrix_temp_scale[plot_gene,]))+scale_colour_gradientn(colours = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b"),limits=c(min(matrix_temp_scale[plot_gene,]),max(matrix_temp_scale[plot_gene,])))+labs(title=paste("Cluster",LSEC_population,plot_gene,patient_name)))
dev.off()


}




expressed_genes <- row.names(subset(fData(HSMM_myo), num_cells_expressed >= 100))

diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)",cores=8)
sig_gene_names <- row.names(subset(diff_test_res, (qval < 10**-300 & use_for_ordering=="TRUE")))




pdf(file=paste0("hsBCL9_pbmc_partM1M20_2_4",LSEC_population,patient_name,"_finishImages.pdf"),10,20)

print(plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                              num_clusters = 6,
                              cores = 8,
                              show_rownames = T))

dev.off()




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


setwd("E:/BaiduNetdiskDownload/Manuscript of scRNAseq BCL9 KD in mouse TIME/revised/TAM")

#plot GSEA 0
#########
REST<-read.delim("culster0_C5_my_analysis.Gsea.1607311661915/gsea_report_for_FALSE_1607311661915.tsv",header = T)
test_group<-read.delim("culster0_C5_my_analysis.Gsea.1607311661915/gsea_report_for_TRUE_1607311661915.tsv",header = T)




REST_20<-subset(REST, GS.DETAILS!="")
REST_20$BCL9_Perturbation<-c("FALSE")
test_group_20<-subset(test_group, GS.DETAILS!="")
test_group_20$BCL9_Perturbation<-c("TRUE")
All<-rbind(REST_20, test_group_20)

bar_all_KEGG<-ggplot(All, aes(reorder(NAME, NES), NES, fill=BCL9_Perturbation)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c(palette_group2and3))+#
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA (M0 in MSigDB C5)")

pdf(file="M0_in_BCL9_Perturbation C5_Pathways_NES_from_GSEA.pdf",width=10, height=10)
bar_all_KEGG
dev.off()




###########################################

REST<-read.delim("culster0_C7_my_analysis.Gsea.1607310337087/gsea_report_for_FALSE_1607310337087.tsv",header = T)
test_group<-read.delim("culster0_C7_my_analysis.Gsea.1607310337087/gsea_report_for_TRUE_1607310337087.tsv",header = T)




REST_20<-subset(REST, GS.DETAILS!="")
REST_20$BCL9_Perturbation<-c("FALSE")
test_group_20<-subset(test_group, GS.DETAILS!="")
test_group_20$BCL9_Perturbation<-c("TRUE")
All<-rbind(REST_20, test_group_20)

bar_all_KEGG<-ggplot(All, aes(reorder(NAME, NES), NES, fill=BCL9_Perturbation)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c(palette_group2and3))+#
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA (M0 in MSigDB C7)")

pdf(file="M0_in_BCL9_Perturbation C7_Pathways_NES_from_GSEA.pdf",width=12, height=10)
bar_all_KEGG
dev.off()



#plot GSEA 2
#########

REST<-read.delim("culster2_C5_my_analysis.Gsea.1607310375587/gsea_report_for_FALSE_1607310375587.tsv",header = T)
test_group<-read.delim("culster2_C5_my_analysis.Gsea.1607310375587/gsea_report_for_TRUE_1607310375587.tsv",header = T)

test_group[which(test_group$NES=="---"),"NES"]<-min(na.omit(as.numeric(test_group$NES)))


REST_20<-subset(REST, GS.DETAILS!="")
REST_20$BCL9_Perturbation<-c("FALSE")
test_group_20<-subset(test_group, GS.DETAILS!="")
test_group_20$BCL9_Perturbation<-c("TRUE")
All<-rbind(REST_20, test_group_20)

All$NES<-as.numeric(All$NES)

bar_all_KEGG<-ggplot(All, aes(reorder(NAME, NES), NES, fill=BCL9_Perturbation)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c(palette_group2and3))+#
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA (M1 in MSigDB C5)")

pdf(file="M1_in_BCL9_Perturbation C5_Pathways_NES_from_GSEA.pdf",width=12, height=10)
bar_all_KEGG
dev.off()




###########################################


REST<-read.delim("culster2_C7_my_analysis.Gsea.1607312401719/gsea_report_for_FALSE_1607312401719.tsv",header = T)
test_group<-read.delim("culster2_C7_my_analysis.Gsea.1607312401719/gsea_report_for_TRUE_1607312401719.tsv",header = T)

test_group[which(test_group$NES=="---"),"NES"]<-min(na.omit(as.numeric(test_group$NES)))


REST_20<-subset(REST, GS.DETAILS!="")
REST_20$BCL9_Perturbation<-c("FALSE")
test_group_20<-subset(test_group, GS.DETAILS!="")
test_group_20$BCL9_Perturbation<-c("TRUE")
All<-rbind(REST_20, test_group_20)
All$NES<-as.numeric(All$NES)
bar_all_KEGG<-ggplot(All, aes(reorder(NAME, NES), NES, fill=BCL9_Perturbation)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c(palette_group2and3))+#
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA (M1 in MSigDB C7)")



pdf(file="M1_in_BCL9_Perturbation C7_Pathways_NES_from_GSEA.pdf",width=12, height=10)
bar_all_KEGG
dev.off()



#####################

#plot GSEA 4
#########

REST<-read.delim("culster4_C5_my_analysis.Gsea.1607316685244/gsea_report_for_FALSE_1607316685244.tsv",header = T)
test_group<-read.delim("culster4_C5_my_analysis.Gsea.1607316685244/gsea_report_for_TRUE_1607316685244.tsv",header = T)




REST_20<-subset(REST, GS.DETAILS!="")
REST_20$BCL9_Perturbation<-c("FALSE")
test_group_20<-subset(test_group, GS.DETAILS!="")
test_group_20$BCL9_Perturbation<-c("TRUE")
All<-rbind(REST_20, test_group_20)

bar_all_KEGG<-ggplot(All, aes(reorder(NAME, NES), NES, fill=BCL9_Perturbation)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c(palette_group2and3))+#
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA (M2 in MSigDB C5)")

pdf(file="M2_in_BCL9_Perturbation C5_Pathways_NES_from_GSEA.pdf",width=12, height=10)
bar_all_KEGG
dev.off()





###########################################


REST<-read.delim("culster4_C7_my_analysis.Gsea.1607317551424/gsea_report_for_FALSE_1607317551424.tsv",header = T)
test_group<-read.delim("culster4_C7_my_analysis.Gsea.1607317551424/gsea_report_for_TRUE_1607317551424.tsv",header = T)




REST_20<-subset(REST, GS.DETAILS!="")
REST_20$BCL9_Perturbation<-c("FALSE")
test_group_20<-subset(test_group, GS.DETAILS!="")
test_group_20$BCL9_Perturbation<-c("TRUE")
All<-rbind(REST_20, test_group_20)

bar_all_KEGG<-ggplot(All, aes(reorder(NAME, NES), NES, fill=BCL9_Perturbation)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c(palette_group2and3))+#
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA (M2 in MSigDB C7)")

pdf(file="M2_in_BCL9_Perturbation C7_Pathways_NES_from_GSEA.pdf",width=12, height=10)
bar_all_KEGG
dev.off()




##########################
###########################################


REST<-read.delim("culster2_C3_my_analysis.Gsea.1607394614083/gsea_report_for_FALSE_1607394614083.tsv",header = T)
test_group<-read.delim("culster2_C3_my_analysis.Gsea.1607394614083/gsea_report_for_TRUE_1607394614083.tsv",header = T)
test_group[which(test_group$NES=="---"),"NES"]<-min(na.omit(as.numeric(test_group$NES)))



REST_20<-subset(REST, GS.DETAILS!="")
REST_20$BCL9_Perturbation<-c("FALSE")
test_group_20<-subset(test_group, GS.DETAILS!="")
test_group_20$BCL9_Perturbation<-c("TRUE")
All<-rbind(REST_20, test_group_20)
All$NES<-as.numeric(All$NES)
bar_all_KEGG<-ggplot(All, aes(reorder(NAME, NES), NES, fill=BCL9_Perturbation)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c(palette_group2and3))+#
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Transcription factor (M1 in MSigDB C3 TF)")

pdf(file="M1_in_BCL9_Perturbation C3 TF_Pathways_NES_from_GSEA.pdf",width=12, height=10)
bar_all_KEGG
dev.off()


###########################################


REST<-read.delim("culster0_C3_my_analysis.Gsea.1607394273302/gsea_report_for_FALSE_1607394273302.tsv",header = T)
test_group<-read.delim("culster0_C3_my_analysis.Gsea.1607394273302/gsea_report_for_TRUE_1607394273302.tsv",header = T)




REST_20<-subset(REST, GS.DETAILS!="")
REST_20$BCL9_Perturbation<-c("FALSE")
test_group_20<-subset(test_group, GS.DETAILS!="")
test_group_20$BCL9_Perturbation<-c("TRUE")
All<-rbind(REST_20, test_group_20)

bar_all_KEGG<-ggplot(All, aes(reorder(NAME, NES), NES, fill=BCL9_Perturbation)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c(palette_group2and3))+#
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Transcription factor (M0 in MSigDB C3 TF)")

pdf(file="M0_in_BCL9_Perturbation C3 TF_Pathways_NES_from_GSEA.pdf",width=12, height=10)
bar_all_KEGG
dev.off()


###########################################


REST<-read.delim("culster4_C3_my_analysis.Gsea.1607394625031/gsea_report_for_FALSE_1607394625031.tsv",header = T)
test_group<-read.delim("culster4_C3_my_analysis.Gsea.1607394625031/gsea_report_for_TRUE_1607394625031.tsv",header = T)




REST_20<-subset(REST, GS.DETAILS!="")
REST_20$BCL9_Perturbation<-c("FALSE")
test_group_20<-subset(test_group, GS.DETAILS!="")
test_group_20$BCL9_Perturbation<-c("TRUE")
All<-rbind(REST_20, test_group_20)

bar_all_KEGG<-ggplot(All, aes(reorder(NAME, NES), NES, fill=BCL9_Perturbation)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c(palette_group2and3))+#
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Transcription factor (M2 in MSigDB C3 TF)")

pdf(file="M2_in_BCL9_Perturbation C3 TF_Pathways_NES_from_GSEA.pdf",width=12, height=10)
bar_all_KEGG
dev.off()



