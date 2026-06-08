# qlogin -l mem_req=10G,s_vmem=10G -pe def_slot 8

### TEST R4.0
# module load R/4.0.2
# export R_LIBS=~/R_PRJ/R4.0_def
# export HDF5_DIR=~/tools/local/hdf5;
# export LD_LIBRARY_PATH=${HDF5_DIR}/lib:$LD_LIBRARY_PATH;
# export PATH=~/tools/local/FIt-SNE/bin:/home/kuze/tools/local/FIt-SNE/:$PATH

library(dplyr)
library(Seurat)
library(ggplot2)
library(reticulate)
library(reshape2)
library(cowplot)
# library(ggpubr)

### Function
MakerPlot <- function( OBJ , GeneSet , prj = "prg" , pmethod ="umap",cluster = "seurat_clusters" , width=16, height=16 ){
	file = paste0(outputdir,"/MarkerSet_",prj,".pdf")
	pdf(file=file, paper="a4r",height = height ,width= width)
	p1 <- DimPlot(object = OBJ, group.by= cluster,reduction = pmethod,
		pt.size = 0.5 ,label=T )+ ggtitle(label = paste0(cluster,"_",pmethod))
	for(i in 1:length(GeneSet)){
		Fp <- FeaturePlot(object= OBJ , features=GeneSet[i] ,reduction=pmethod, pt.size = 0.5)
		Vp <- VlnPlot(object = OBJ , features=GeneSet[i])
		p <- plot_grid( Fp, Vp , p1 , ncol =2 , labels="auto", align="h" )
		print(p)
	}
	dev.off()
	print(paste0("finish ", prj ," plot!!" ) )
}


##################################################
# Read Data set
##################################################
sample_inf <- read.delim( "sample_list.txt"  ,header = T, stringsAsFactors=F ,sep="\t")
donar <-  commandArgs(trailingOnly=TRUE)[1]

EachSampleDir <- "EachSample/"
outputdir <- paste0( "Integrate_donar/", donar ,"/" )
dir.create( outputdir , recursive = T)

matching_samples <- sample_inf %>%
  filter(Donar == donar ) %>%
  pull(Sample)


##################################################
# Merge & Cluster
##################################################

obj.list <- NULL
for ( Sample in matching_samples ){
	rds <- readRDS( paste0( EachSampleDir , Sample , "/SeuratOBJ_GEX_",Sample,"_SingleR.RData") )
	rds[["Sample"]] <- Sample
	rds[["Donar"]] <- donar
	obj.list <- c( obj.list , rds )
}

saveRDS( obj.list  ,paste0( outputdir , "OBJlist.rds") )

names(obj.list) <- matching_samples
SeuratOBJ <- merge(x = obj.list[[1]] , y = obj.list[2:length(obj.list) ] )
saveRDS( SeuratOBJ , paste0( outputdir , "rawSeuratOBJ.rds") )


#SeuratOBJ[["RNA"]] <- split(SeuratOBJ[["RNA"]], f = SeuratOBJ$Sample , layers="counts" )
SeuratOBJ <- NormalizeData(SeuratOBJ)
SeuratOBJ <- FindVariableFeatures(SeuratOBJ)
SeuratOBJ <- ScaleData(SeuratOBJ)
SeuratOBJ <- RunPCA(SeuratOBJ)

# After preprocessing, we integrate layers with added parameters specific to Harmony:
SeuratOBJ <- IntegrateLayers(object = SeuratOBJ, method = HarmonyIntegration, orig.reduction = "pca",
  new.reduction = 'harmony', verbose = FALSE)

SeuratOBJ <- RunUMAP(SeuratOBJ, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")
SeuratOBJ <- FindNeighbors(SeuratOBJ, reduction = "harmony", dims = 1:20)
SeuratOBJ <- FindClusters(SeuratOBJ, resolution = 0.5, cluster.name = "harmony_cluster")

saveRDS( SeuratOBJ , paste0( outputdir , "Harmony_integrated_SeuratOBJ.rds")) 

#  reduction = "umap.cca"
##################################################
# Plot
#################################################

filePath= paste0(outputdir,"/Integrated_Cluster_",donar,".pdf" )
p<-DimPlot(object = SeuratOBJ , group.by="harmony_cluster",reduction = "umap.harmony", pt.size = 0.2 , label=T ) + ggtitle(label = "Cluster UMAP ")
ggsave(file = filePath, plot = p, dpi=100, width=5, height=5)

filePath=paste0(outputdir,"/Integrated_Sample_",donar,".pdf")
p<-DimPlot(object = SeuratOBJ , group.by="Sample",reduction = "umap.harmony", pt.size = 0.2 ,label=T ) + ggtitle(label = "Sample UMAP")
ggsave(file = filePath, plot = p, dpi=100, width=5, height=5)

filePath=paste0(outputdir,"/Integrated_CART_",donar,".pdf")
p<-DimPlot(object = SeuratOBJ , group.by="CART",reduction = "umap.harmony", pt.size = 0.2 ,label=T ) + ggtitle(label = "Sample UMAP")
ggsave(file = filePath, plot = p, dpi=100, width=5, height=5)

filePath=paste0(outputdir,"/Integrated_CellType_",donar,".pdf")
p<-DimPlot(object = SeuratOBJ , group.by="SingleR.monaco",reduction = "umap.harmony", pt.size = 0.2 ,label=T ) + ggtitle(label = "Sample UMAP")
ggsave(file = filePath, plot = p, dpi=100, width=16, height=5)




