# qlogin -l mem_req=10G,s_vmem=10G -pe def_slot 8

### TEST R4.0
# module load R/4.0.2
# export R_LIBS=~/R_PRJ/R4.0_def
# export LD_LIBRARY_PATH=${HDF5_DIR}/lib:$LD_LIBRARY_PATH;

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

id <- commandArgs(trailingOnly=TRUE)[1]
Sample <- sample_inf[ id,"Sample" ]
Donar  <- sample_inf[ id,"Donar" ]

outputdir <- paste0( "EachSample/",Sample ,"/" )
dir.create( outputdir , recursive = T)

SeuratOBJ.data <- Read10X(data.dir = paste0("CR/",Sample,"/per_sample_outs/" ,Sample, "/count/sample_filtered_feature_bc_matrix/" ) )
SeuratOBJ <- CreateSeuratObject(counts = SeuratOBJ.data, project = Sample)
SeuratOBJ[["percent.mt"]] <- PercentageFeatureSet(SeuratOBJ, pattern = "^MT-")
SeuratOBJ[["percent.rib"]] <- PercentageFeatureSet(SeuratOBJ,pattern="^RP[LS]") 
# TCR
VDJdir  <- paste0("CR/",Sample,"/per_sample_outs/" ,Sample, "/vdj_t/" )
TCRfile <- paste0( VDJdir , "/filtered_contig_annotations.csv" )
TCRdata <- read.delim(TCRfile  ,header = T, stringsAsFactors=F ,sep=",")
TCRdata <- TCRdata[ TCRdata$is_cell == "true" ,c("barcode","raw_clonotype_id")]

TCloneFile <- paste0(VDJdir,"/clonotypes.csv")
TClonaData <- read.delim(TCloneFile  ,header = T,stringsAsFactors=F ,sep=",")

uqTCRdata <- unique(TCRdata)
uqTCRdata <- uqTCRdata[uqTCRdata[,2] != "",]
rownames(uqTCRdata) <- uqTCRdata$barcode
colnames(uqTCRdata) <- c("barcode","clonotype_id")
TCRinf <- left_join(uqTCRdata, TClonaData, by = "clonotype_id")
TCRinf[is.na(TCRinf)] <- "None"
colnames(TCRinf) <- c("barcode" , "TCR_clonotype_id" , "TCR_frequency", "TCR_proportion","TCR_cdr3s_aa","TCR_cdr3s_nt","TCR_inkt_evidence","TCR_mait_evidence")

# CART
CART_file <- paste0("CART_data/CART_",Sample,"_cb_list.txt")
CARTdata <- read.delim(CART_file  ,header = T, stringsAsFactors=F ,sep=",")

#RNA seq と対応付け
RNAseq_bc <-  as.data.frame( colnames(SeuratOBJ), stringsAsFactors=F)
colnames(RNAseq_bc) <- "barcode"

RNAvdj <- left_join(RNAseq_bc, TCRinf, by = "barcode")
RNAvdj <- left_join(RNAvdj , CARTdata, by = "barcode")
#RNAvdj <- left_join(RNAseq_bc, VDJinf ,BCRinf, by = "barcode")

RNAvdj[is.na(RNAvdj)] <- "None"
rownames(RNAvdj) <- RNAvdj[,"barcode"]
SeuratOBJ <- AddMetaData(SeuratOBJ, metadata = RNAvdj)
 SeuratOBJ[["VDJ"]] <- if_else (SeuratOBJ$TCR_clonotype_id  == "None", "None" ,"VDJdetceted")

SeuratOBJ <- RenameCells( SeuratOBJ ,  add.cell.id = Sample)
saveRDS( SeuratOBJ , paste0(outputdir ,"/raw_SeuratOBJ.rds"))


###############################################
#
###############################################
# SeuratOBJ <-readRDS( paste0(outputdir ,"/raw_SeuratOBJ.rds"))

png( paste0(outputdir,"/VlnPlot_",Sample,"_GEX.png"), width = 1000, height = 500)
VlnPlot(SeuratOBJ, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rib"), ncol = 4, group.by ="VDJ")
dev.off()

pdf( paste0(outputdir,"/FeatureScatter_",Sample,"_GEX.pdf") , width = 16, height = 16)
plot1 <- FeatureScatter(SeuratOBJ, feature1 = "percent.rib", feature2 = "percent.mt")
plot2 <- FeatureScatter(SeuratOBJ, feature1 = "percent.rib", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(SeuratOBJ, feature1 = "percent.rib", feature2 = "nCount_RNA")
plot4 <- FeatureScatter(SeuratOBJ, feature1 = "percent.rib", feature2 = "CD3E", slot = "counts")
(plot1 + plot2) / ( plot3+ plot4)
dev.off()


#SeuratOBJ <- subset(SeuratOBJ, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000  & nCount_RNA > 2500 &  nCount_RNA < 10000 & percent.mt < 10 &  percent.rib < 5 )
SeuratOBJ <- subset(SeuratOBJ, subset = nFeature_RNA > 500 & percent.mt < 30 &  percent.rib > 1 )
png( paste0(outputdir,"/VlnPlot_",Sample,"_GEX_AfterQC.png"), width = 1000, height = 500)
VlnPlot(SeuratOBJ, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rib"), ncol = 4, group.by ="VDJ")
dev.off()


SeuratOBJ <- NormalizeData(SeuratOBJ)
SeuratOBJ <- FindVariableFeatures(SeuratOBJ, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(SeuratOBJ)
SeuratOBJ <- ScaleData(SeuratOBJ , features = all.genes)
SeuratOBJ <- RunPCA(SeuratOBJ, npcs = 10)



filePath= paste0(outputdir,"/VizDimLoadings.pdf")
p<-VizDimLoadings(SeuratOBJ, dims = 1:10, reduction = "pca")
ggsave(file = filePath, plot = p, dpi=50, width=16, height=16)

filePath= paste0(outputdir,"/PCAPlot.pdf")
p<-DimPlot(SeuratOBJ, reduction = "pca", group.by="orig.ident")
ggsave(file = filePath, plot = p, dpi=100, width=5, height=5)

filePath <- paste0(outputdir,"/PCheatmap.pdf")
p<-DimHeatmap(SeuratOBJ, dims = 1:10, cells = 500, balanced = TRUE,fast = FALSE)
ggsave(file = filePath, plot = p, dpi=100, width=24, height=24)

filePath= paste0(outputdir,"/VlnPlot_ElbowPlot.pdf")
p<-ElbowPlot(SeuratOBJ , ndims=20 )
ggsave(file = filePath, plot = p, dpi=100, width=8, height=8)

pcadim=c(1:10)
resol = 0.5

SeuratOBJ <- FindNeighbors(SeuratOBJ, reduction = "pca", dims = pcadim)
SeuratOBJ <- FindClusters(SeuratOBJ, resolution = resol)
SeuratOBJ <- RunUMAP(SeuratOBJ, reduction = "pca", dims = pcadim)


freq_table <- prop.table(x = table(Idents(SeuratOBJ), SeuratOBJ@meta.data[, "orig.ident"]), margin = 2)
freq_mmod<-melt(t(freq_table),id=c("sample","cluster","val"))
colnames(freq_mmod) <- c("sample","cluster","val")
freq_mmod$cluster <- as.factor( freq_mmod$cluster )

sink( paste0(outputdir,"/frek_table.txt") )
table(Idents(SeuratOBJ), SeuratOBJ@meta.data[, "orig.ident"])
sink()

sink(  paste0(outputdir,"/freq_tableprop.txt"))
freq_table
sink()

sink( paste0(outputdir,"/frek_table_VDJ.txt") )
table(Idents(SeuratOBJ), SeuratOBJ@meta.data[, "VDJ"])
sink()
sink( paste0(outputdir,"/frek_table_CART.txt") )
table(Idents(SeuratOBJ), SeuratOBJ@meta.data[, "CART"])
sink()


filePath = paste0(outputdir,"/dimreduction_cluster_umap.pdf")
p <- DimPlot(object = SeuratOBJ, group.by="seurat_clusters",reduction = "umap", pt.size = 0.2 ,
                 label=T )+ ggtitle(label = "UMAP cluster")
ggsave(file = filePath, plot = p, dpi=120, width=8, height=8)

filePath = paste0(outputdir,"/dimreduction_CART_umap.pdf")
p <- DimPlot(object = SeuratOBJ, group.by="CART",reduction = "umap", pt.size = 0.2 ,
                 label=T )+ ggtitle(label = "UMAP cluster")
ggsave(file = filePath, plot = p, dpi=120, width=8, height=8)

filePath = paste0(outputdir,"/dimreduction_VDJ_umap.pdf")
p <- DimPlot(object = SeuratOBJ, group.by="VDJ",reduction = "umap", pt.size = 0.2 ,
                 label=T )+ ggtitle(label = "UMAP cluster")
ggsave(file = filePath, plot = p, dpi=120, width=8, height=8)



SeuratOBJ.markers <- FindAllMarkers(SeuratOBJ, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(SeuratOBJ.markers, file= paste0(outputdir,"/all_Seurat_clsuterMarker.txt"), quote=F, col.names=NA , sep="\t")

SeuratOBJ.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- SeuratOBJ.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

filePath= paste0(outputdir,"/DEG_heatmap.pdf")
p<-DoHeatmap(SeuratOBJ, features = top10$gene,raster =FALSE, group.by="seurat_clusters" )
ggsave(file = filePath, plot = p, dpi=100, width=16, height=16)

filePath= paste0(outputdir,"/ClusterQC_Vinplot","",".pdf")
p=VlnPlot( SeuratOBJ , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file = filePath, plot = p, dpi=100, width=16, height=5)


cellMaker <- read.table("./PBMC_marker.txt" , sep="\t" , header= T )
cellMaker[ !cellMaker$Gene %in% rownames(SeuratOBJ)  , ]
for ( i in unique(cellMaker$CellType)){
	        MakerPlot( SeuratOBJ , cellMaker[ cellMaker$CellType == i, 1] , i )
}
saveRDS(SeuratOBJ, file =  paste0( outputdir, "/SeuratOBJ_GEX_",Sample,"_ana_OBJ.RData") )

########################################
# Singler Annotation 
########################################

library(celldex)
ref <- MonacoImmuneData()

library(SingleR)
DefaultAssay(SeuratOBJ) <- "RNA"
SeuratOBJ.sce <- as.SingleCellExperiment( SeuratOBJ)
pred.SeuratOBJ <- SingleR(test = SeuratOBJ.sce , ref = ref , assay.type.test=1 , labels =ref$label.fine)
SeuratOBJ[["SingleR.monaco"]] <- pred.SeuratOBJ$labels

max.scoreMat <- apply( pred.SeuratOBJ$scores , 1, function(x){ return( c( celltype=  names(which.max( x)), score = max(x)  ))} ) %>% t %>% as.data.frame
rownames(max.scoreMat) <- rownames(pred.SeuratOBJ)
max.scoreMat <- max.scoreMat %>% mutate( cell = rownames(.) , Sample = Sample , score = as.numeric(score))

SeuratOBJ[["SingleR.score"]] <- max.scoreMat$score
#score_low_val <- 0.2
#SeuratOBJ[["SingleR.filter"]] <-  SeuratOBJ@meta.data %>% 
#									select( "SingleR.monaco" ,"SingleR.score" ) %>% 
#									mutate( SingleR.filter = if_else( SingleR.score > score_low_val , SingleR.monaco , "Undef"  )) %>%
#									pull(SingleR.filter)

cellMaker2 <- read.table("./Yoshitake_Marker_v2.txt" , sep="\t" , header= T )
cellMaker2[ !cellMaker2$Gene %in% rownames(SeuratOBJ)  , ]

sink( paste0(outputdir,"/frek_table_Annotation_celltype_data.txt") )
table( SeuratOBJ@meta.data[, "SingleR.monaco"] , SeuratOBJ@meta.data[, "seurat_clusters"])
sink()

sink( paste0(outputdir,"/frek_table_Annotation_celltype_data_CART.txt") )
table( SeuratOBJ@meta.data[, "SingleR.monaco"] , SeuratOBJ@meta.data[, "CART"])
sink()

sink( paste0(outputdir,"/frek_table_Annotation_celltype_data_VDJ.txt") )
table( SeuratOBJ@meta.data[, "SingleR.monaco"] , SeuratOBJ@meta.data[, "VDJ"])
sink()


filePath= paste0(outputdir,"/SingleR_UMAP.pdf" )
p<-DimPlot(object = SeuratOBJ , group.by="SingleR.monaco",reduction = "umap", pt.size = 0.2 ,
         label=T )+ ggtitle(label = "UMAP sample")
ggsave(file = filePath, plot = p, dpi=100, width=16, height=6)
filePath=paste0(outputdir,"/SeuratCluster_UMAP.pdf")
p<-DimPlot(object = SeuratOBJ , group.by="seurat_clusters",reduction = "umap", pt.size = 0.2 ,label=T )+ ggtitle(label = "UMAP sample")
ggsave(file = filePath, plot = p, dpi=100, width=16, height=16)
filePath=paste0(outputdir,"/SeuratCluster_DotPlot.pdf")
p<-DotPlot(SeuratOBJ , features = unique(cellMaker2[,1]), dot.scale = 8,group.by="seurat_clusters") +RotatedAxis()
ggsave(file = filePath, plot = p, dpi=100, width=16, height=6)
filePath=paste0(outputdir, "/SingleR_DotPlot.pdf")
p<-DotPlot(SeuratOBJ  , features = unique(cellMaker2[,1]), dot.scale = 8,group.by="SingleR.monaco") +RotatedAxis()
ggsave(file = filePath, plot = p, dpi=100, width=16, height=6)

TMaker <- c( cellMaker2[cellMaker2[,2] == "Tcell",1] , ""  )
TMaker <- c( TMaker ,"GNLY","NKG7","FAS","CD28")
p1 <-DotPlot(SeuratOBJ, features = unique(TMaker), dot.scale = 8,group.by="seurat_clusters") +RotatedAxis()
filePath= paste0(outputdir,"/Dotplot_Cluster.pdf" )
ggsave(file = filePath, plot = p1, dpi=100, width=16, height=6)

p2 <-DotPlot(SeuratOBJ, features = unique(TMaker), dot.scale = 8,group.by="SingleR.monaco") +RotatedAxis()
filePath= paste0(outputdir,"/Dotplot_SingleR.pdf" )
ggsave(file = filePath, plot = p2, dpi=100, width=16, height=6)


saveRDS(SeuratOBJ, file =  paste0( outputdir, "/SeuratOBJ_GEX_",Sample,"_SingleR.RData") )
