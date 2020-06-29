library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)

# Merging the 4 datasets of single nuclei RNAseq

sol.data <- Read10X(data.dir = "./Raw_features/sol/")
sol <- CreateSeuratObject(counts = sol.data, project = "Sol")

quad.data <- Read10X(data.dir = "./Raw_features/quad/")
quad <- CreateSeuratObject(counts = quad.data, project = "quad")

Tib.data <- Read10X(data.dir = "./Raw_features/Tib/")
Tib <- CreateSeuratObject(counts = Tib.data, project = "Tib")

Fastmix.data <- Read10X(data.dir = "./Raw_features/Fast_mix/")
Fastmix <- CreateSeuratObject(counts = Fastmix.data, project = "Fastmix")


merge.data <- merge(sol, y=c(quad, Tib, Fastmix), add.cell.ids = c("sol", "quad", "Tib", "Fastmix"), project = "Full Merge Project")
table(merge.data$orig.ident)


merge.data[["percent.mt"]] <- PercentageFeatureSet(merge.data, pattern = "^mt-")
mergeData <- subset(merge.data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
table(mergeData$orig.ident)


mergeData.list <- SplitObject(object = mergeData, split.by = "orig.ident")

for (i in 1:length(x = mergeData.list)) {
    mergeData.list[[i]] <- NormalizeData(object = mergeData.list[[i]], verbose = FALSE)
    mergeData.list[[i]] <- FindVariableFeatures(object = mergeData.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}


mergeData.anchors <- FindIntegrationAnchors(object.list = mergeData.list, dims = 1:30)
mergeData.integrated <- IntegrateData(anchorset = mergeData.anchors, dims = 1:30)

DefaultAssay(mergeData.integrated) <- "integrated"


# Run the standard workflow for visualization and clustering

mergeData.integrated <- ScaleData(object = mergeData.integrated, verbose = FALSE)
mergeData.integrated <- RunPCA(object = mergeData.integrated, npcs = 30, verbose = FALSE)
mergeData.integrated <- RunUMAP(object = mergeData.integrated, reduction = "pca", dims = 1:30)

mergeData.integrated <- readRDS("./mergeDataIntegrated.rds")

name <- "BI19-010_scRNA_MergedData_Sol-Quad-Tib-Fastmix"

mergeData.integrated <- FindNeighbors(mergeData.integrated, dims = 1:20)
mergeData.integrated <- FindClusters(mergeData.integrated, resolution = 0.5)



# UMAP plot of the different cluster 

	DimPlot(object = mergeData.integrated, reduction = "umap", group.by = "orig.ident", label = FALSE)
	DimPlot(mergeData.integrated, reduction = "umap", label = TRUE)

# Names of the clusters (fig1b and sup 2b and 2c)

	new.cluster.ids <- c("Myh1+2 Cells -1", 
		"Myh4 Cells -1",
		"Myh1+2 Cells -2",
		"Myh4 Cells -2",
		"FAPS-1",
		"Myh4 Cells -3",
		"FAPS-2", 
		"Newly fused myonuclei",
		"FAPS-3",
		"Myh7 Cells",
		"MuSC",
		"Tenocyte",
		"Smooth muscular",
		"Endothelial",
		"B/T Cells",
		"Adipocyte", 
		"MTJ Myonuclei",
		"Pericyte",
		"NMJ Myonuclei")
	names(new.cluster.ids) <- levels(mergeData.integrated)
	mergeData.integrated <- RenameIdents(mergeData.integrated, new.cluster.ids)

# UMAP plot of the different cluster with name legends (fig1b and sup 2b and 2c)

png(filename=paste(name,"UMAP.png",sep="_"),width=10 ,height=10, units="in",res = 600 ) 
	DimPlot(mergeData.integrated, reduction = "umap", label = FALSE)
dev.off()

png(filename=paste(name,"UMAP_withClusterNames.png",sep="_"),width=10 ,height=10, units="in",res = 600 ) 
	DimPlot(mergeData.integrated, reduction = "umap", label = TRUE)
dev.off()

png(filename=paste(name,"UMAP_byOrig.png",sep="_"),width=10 ,height=10, units="in",res = 600 ) 
	DimPlot(mergeData.integrated, reduction = "umap", label = FALSE, group.by = "orig.ident")
dev.off()
   


###########################
# Differential expression###########
###########################

DefaultAssay(mergeData.integrated) <- "RNA"

# New colonne cluster_sample

	mergeData.integrated$cluster.sample <- paste(mergeData.integrated$orig.ident, mergeData.integrated$seurat_clusters, sep = "_")


# Specific markers of each clusters

	Idents(mergeData.integrated) <- "seurat_clusters"

	ncluster <- length(levels(mergeData.integrated$seurat_clusters)) - 1
	for (i in 0:ncluster) {
		conservedMarkers <- FindConservedMarkers(mergeData.integrated, ident.1 = i, grouping.var = "orig.ident", verbose = FALSE)
		currentCluster.Name <- new.cluster.ids[i+1]
		nameFile <- paste(name,"FindConservedMarkers_Cluster",currentCluster.Name,".txt",sep="_")
		write.table(conservedMarkers,file=nameFile,sep="\t")
	}
png(filename=paste(name,"FeaturePlot.png",sep=""),width=10 ,height=7, units="in",res = 600 ) 




# Feature plot of differents markers (sup 2a)

FeaturePlot(mergeData.integrated, features = c("Myh1","Myh2","Myh4","Myh7","Ttn","Pdgfra","Pax7","Pparg","Pecam1","Ptprc","Myh11","Col11a1","Col24a1","Colq"))
dev.off()




###########################
# Differential expression of MTJ myonuclei  (fig 2b)###########
###########################

	data.currentClusters <- mergeData.integrated

	new.cluster.ids <- c("Myonuclei", 
		"Myonuclei",
		"Myonuclei",
		"Myonuclei",
		"FAPS-1",
		"Myonuclei",
		"FAPS-2", 
		"Newly fused myonuclei",
		"FAPS-3",
		"Myonuclei",
		"MuSC",
		"Tenocyte",
		"Smooth muscular",
		"Endothelial",
		"B/T Cells",
		"Adipocyte", 
		"MTJ Myonuclei",
		"Pericyte",
		"NMJ Myonuclei")
	names(new.cluster.ids) <- levels(data.currentClusters)
	data.currentClusters <- RenameIdents(data.currentClusters, new.cluster.ids)
	
	MTJ.markers <- FindMarkers(data.currentClusters, ident.1 = "MTJ Myonuclei", ident.2 = "Myonuclei", min.pct = 0.25, only.pos = TRUE)
	MTJ.markers2 <- FindMarkers(data.currentClusters, ident.2 = "MTJ Myonuclei", ident.1 = "Myonuclei", min.pct = 0.25, only.pos = TRUE)

	
data.currentClusters <- subset(data.currentClusters, idents = c("MTJ Myonuclei","Myonuclei"))
data.currentClusters <- ScaleData(object = data.currentClusters, verbose = FALSE)
	
top30 <- head(MTJ.markers, n=30)
top30.2 <- head(MTJ.markers2, n=30)
png(filename=paste(name,"_Heatmap_MTJ-vs-Myonuclei_30genes.png",sep="_"),width=7 ,height=14, units="in",res = 600 ) 
DoHeatmap(data.currentClusters, features = c(rownames(top30),rownames(top30.2)))


DoHeatmap(subset(data.currentClusters, downsample = 100), features = c(rownames(top30),rownames(top30.2)))




###########################
# Differential expression of NMJ myonuclei  (fig 2b)###########
###########################

	data.currentClusters <- mergeData.integrated

	new.cluster.ids <- c("Myonuclei", 
		"Myonuclei",
		"Myonuclei",
		"Myonuclei",
		"FAPS-1",
		"Myonuclei",
		"FAPS-2", 
		"Newly fused myonuclei",
		"FAPS-3",
		"Myonuclei",
		"MuSC",
		"Tenocyte",
		"Smooth muscular",
		"Endothelial",
		"B/T Cells",
		"Adipocyte", 
		"MTJ Myonuclei",
		"Pericyte",
		"NMJ Myonuclei")
	names(new.cluster.ids) <- levels(data.currentClusters)
	data.currentClusters <- RenameIdents(data.currentClusters, new.cluster.ids)
	
	NMJ.markers <- FindMarkers(data.currentClusters, ident.1 = "NMJ Myonuclei", ident.2 = "Myonuclei", min.pct = 0.25, only.pos = TRUE)
	NMJ.markers2 <- FindMarkers(data.currentClusters, ident.2 = "NMJ Myonuclei", ident.1 = "Myonuclei", min.pct = 0.25, only.pos = TRUE)

	
data.currentClusters <- subset(data.currentClusters, idents = c("NMJ Myonuclei","Myonuclei"))
data.currentClusters <- ScaleData(object = data.currentClusters, verbose = FALSE)
	
top30 <- head(NMJ.markers, n=30)
top30.2 <- head(NMJ.markers2, n=30)
png(filename=paste(name,"_Heatmap_NMJ-vs-Myonuclei_30genes.png",sep="_"),width=7 ,height=14, units="in",res = 600 ) 
DoHeatmap(data.currentClusters, features = c(rownames(top30),rownames(top30.2)))
dev.off()

DoHeatmap(subset(data.currentClusters, downsample = 100), features = c(rownames(top30),rownames(top30.2)))


###########################
# Differential expression of Unknow (newly fused) myonuclei  (fig 2b)###########
###########################


	data.currentClusters <- mergeData.integrated

	new.cluster.ids <- c("Myonuclei", 
		"Myonuclei",
		"Myonuclei",
		"Myonuclei",
		"FAPS-1",
		"Myonuclei",
		"FAPS-2", 
		"Newly fused myonuclei",
		"FAPS-3",
		"Myonuclei",
		"MuSC",
		"Tenocyte",
		"Smooth muscular",
		"Endothelial",
		"B/T Cells",
		"Adipocyte", 
		"MTJ Myonuclei",
		"Pericyte",
		"NMJ Myonuclei")
	names(new.cluster.ids) <- levels(data.currentClusters)
	data.currentClusters <- RenameIdents(data.currentClusters, new.cluster.ids)
	
	Newlyfused.markers <- FindMarkers(data.currentClusters, ident.1 = "Newly fused myonuclei", ident.2 = "Myonuclei", min.pct = 0.25, only.pos = TRUE)
	Newlyfused.markers2 <- FindMarkers(data.currentClusters, ident.2 = "Newly fused myonuclei", ident.1 = "Myonuclei", min.pct = 0.25, only.pos = TRUE)

	
data.currentClusters <- subset(data.currentClusters, idents = c("Newly fused myonuclei","Myonuclei"))
data.currentClusters <- ScaleData(object = data.currentClusters, verbose = FALSE)
	
top30 <- head(Newlyfused.markers, n=30)
top30.2 <- head(Newlyfused.markers2, n=30)
png(filename=paste(name,"_Heatmap_NewlyFused-vs-Myonuclei_30genes.png",sep="_"),width=7 ,height=14, units="in",res = 600 ) 
DoHeatmap(data.currentClusters, features = c(rownames(top30),rownames(top30.2)))
dev.off()


DoHeatmap(subset(data.currentClusters, downsample = 100), features = c(rownames(top30),rownames(top30.2)))


##########################################"


# Number of cells per cluster and conditions (Sup 2d)

	tableData.nCells <- table(mergeData.integrated@active.ident, mergeData.integrated@meta.data$orig.ident)
	tableNOrig <- paste(name,"nCells_byCluster_bySample.txt",sep="_")
	write.table(tableData.nCells,file=tableNOrig,sep="\t")




# Myh expression analysis
	
	# Global expression of Myh
	myh.features <- c("Myh1","Myh2","Myh4","Myh7")
	FeaturePlot(mergeData.integrated, features = myh.features)
	RidgePlot(mergeData.integrated, features = myh.features, ncol = 2)
	
	# Table of Myh expression
	myh.cells <- subset(mergeData.integrated, idents = c("Myh2 Cells","Myh7 Cells","Myh1 Cells","Myh4 Cells"))
	table(myh.cells@active.ident, myh.cells@meta.data$orig.ident)
	
		
	data.currentClusters <- mergeData.integrated

	new.cluster.ids <- c("Myonuclei", 
		"Myonuclei",
		"Myonuclei",
		"Myonuclei",
		"FAPS-1",
		"Myonuclei",
		"FAPS-2", 
		"Newly fused myonuclei",
		"FAPS-3",
		"Slow Myonuclei",
		"MuSC",
		"Tenocyte",
		"Smooth muscular",
		"Endothelial",
		"B/T Cells",
		"Adipocyte", 
		"MTJ Myonuclei",
		"Pericyte",
		"NMJ Myonuclei")
	names(new.cluster.id
	names(new.cluster.ids) <- levels(data.currentClusters)
	data.currentClusters <- RenameIdents(data.currentClusters, new.cluster.ids)
	
	
# Coexpression of two Myh in myonuclei?	

		myh.cells <- subset(data.currentClusters, idents = c("Myonuclei","Slow Myonuclei"))
png(filename=paste(name,"FeatureScatter_Myh7vsMyh1.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh7", feature2 = "Myh1", slot = "counts")
dev.off()	
	
		myh.cells <- subset(data.currentClusters, idents = c("Myonuclei","Slow Myonuclei"))
png(filename=paste(name,"FeatureScatter_Myh7vsMyh4.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh7", feature2 = "Myh4", slot = "counts")
dev.off()	
	
			myh.cells <- subset(data.currentClusters, idents = c("Myonuclei","Slow Myonuclei"))
png(filename=paste(name,"FeatureScatter_Myh2vsMyh4.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh2", feature2 = "Myh4", slot = "counts")
dev.off()	


	myh.cells <- subset(data.currentClusters, idents = c("Myonuclei","Slow Myonuclei"))
png(filename=paste(name,"FeatureScatter_Myh7vsMyh2_new.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh7", feature2 = "Myh2", slot = "counts")
dev.off()	
	
	myh.cells <- subset(data.currentClusters, idents = c("Myonuclei","Slow Myonuclei"))
png(filename=paste(name,"FeatureScatter_Myh1vsMyh2_new.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh2", feature2 = "Myh1", slot = "counts")
dev.off()	
	
	
	myh.cells <- subset(data.currentClusters, idents = c("Myonuclei","Slow Myonuclei"))
png(filename=paste(name,"FeatureScatter_Myh1vsMyh4_new.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh1", feature2 = "Myh4", slot = "counts")
dev.off()	
	
	

	
#Dotplot of Myh and Ttn gene expression in the different myonuclei (fig2c)


	data.currentClusters <- mergeData.integrated

	new.cluster.ids <- c("Myonuclei", 
		"Myonuclei",
		"Myonuclei",
		"Myonuclei",
		"FAPS-1",
		"Myonuclei",
		"FAPS-2", 
		"Newly fused myonuclei",
		"FAPS-3",
		"Slow Myonuclei",
		"MuSC",
		"Tenocyte",
		"Smooth muscular",
		"Endothelial",
		"B/T Cells",
		"Adipocyte", 
		"MTJ Myonuclei",
		"Pericyte",
		"NMJ Myonuclei")
	names(new.cluster.ids) <- levels(data.currentClusters)
	data.currentClusters <- RenameIdents(data.currentClusters, new.cluster.ids)



myh.cells <- subset(data.currentClusters, idents = c("Myonuclei", "Slow Myonuclei", "Newly fused myonuclei", "MTJ Myonuclei", "NMJ Myonuclei"))
DefaultAssay(myh.cells) <- "RNA"
png(filename=paste(name,"selectedCluster_DotPlot_Myh7-Myh2-Myh1-Myh4-Ttn_newClusters.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	DotPlot(myh.cells, features = c("Myh1","Myh2","Myh4","Myh7","Ttn"))
dev.off()


	data.currentClusters <- mergeData.integrated

	new.cluster.ids <- c("Myonuclei", 
		"Myonuclei",
		"Myonuclei",
		"Myonuclei",
		"FAPS-1",
		"Myonuclei",
		"FAPS-2", 
		"Myonuclei",
		"FAPS-3",
		"Slow Myonuclei",
		"MuSC",
		"Tenocyte",
		"Smooth muscular",
		"Endothelial",
		"B/T Cells",
		"Adipocyte", 
		"MTJ Myonuclei",
		"Pericyte",
		"NMJ Myonuclei")
	names(new.cluster.ids) <- levels(data.currentClusters)
	data.currentClusters <- RenameIdents(data.currentClusters, new.cluster.ids)






#Heatmap of genes differential expressed in the different type of Myh myonuclei (fig1d)

#~ grep('^Myh1$', rownames(myh.cells), value = FALSE)
#~ grep('^Myh2$', rownames(myh.cells), value = FALSE)
#~ grep('^Myh4$', rownames(myh.cells), value = FALSE)
#~ grep('^Myh7$', rownames(myh.cells), value = FALSE)

#~ 21509
#~ 21507
#~ 21511
#~ 18500

#~ myh.cells["Myh1", ]

length(which(myh.cells[rownames(myh.cells)[21509], ]$nCount_RNA > 30))
length(which(myh.cells[rownames(myh.cells)[21507], ]$nCount_RNA > 30))
length(which(myh.cells[rownames(myh.cells)[21511], ]$nCount_RNA > 30))
length(which(myh.cells[rownames(myh.cells)[18500], ]$nCount_RNA > 20))

myh.cells@meta.data$myh <- 'Neg'
myh.cells@meta.data$myh[which(myh.cells[rownames(myh.cells)[21509], ]$nCount_RNA > 30)] <- 'Myh1'
myh.cells@meta.data$myh[which(myh.cells[rownames(myh.cells)[21507], ]$nCount_RNA > 30)] <- 'Myh2'
myh.cells@meta.data$myh[which(myh.cells[rownames(myh.cells)[21511], ]$nCount_RNA > 30)] <- 'Myh4'
myh.cells@meta.data$myh[which(myh.cells[rownames(myh.cells)[18500], ]$nCount_RNA > 20)] <- 'Myh7'

Idents(myh.cells) <- "myh"

myh.cells.markers <- FindAllMarkers(myh.cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
myh.cells.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


myh.cells.myhClust <- subset(myh.cells, idents = c("Myh1","Myh2","Myh4","Myh7"))

Idents(myh.cells.myhClust) <- "myh"
myh.cells.myhClust <- ScaleData(object = myh.cells.myhClust, verbose = FALSE)


myh.cells.myhClust.markers <- FindAllMarkers(myh.cells.myhClust, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top50 <- myh.cells.myhClust.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
png(filename=paste(name,"_Heatmap_Myh1-2-4-7_Threshold30count_50genes.png",sep="_"),width=7 ,height=14, units="in",res = 600 ) 
DoHeatmap(myh.cells.myhClust, features = top50$gene)
dev.off()

top30 <- myh.cells.myhClust.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
png(filename=paste(name,"_Heatmap_Myh1-2-4-7_Threshold30count_30genes.png",sep="_"),width=7 ,height=14, units="in",res = 600 ) 
DoHeatmap(myh.cells.myhClust, features = top30$gene)
dev.off()

top10 <- myh.cells.myhClust.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
png(filename=paste(name,"_Heatmap_Myh1-2-4-7_Threshold30count_10genes.png",sep="_"),width=7 ,height=14, units="in",res = 600 ) 
DoHeatmap(myh.cells.myhClust, features = top10$gene)
dev.off()

top20 <- myh.cells.myhClust.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
png(filename=paste(name,"_Heatmap_Myh1-2-4-7_Threshold30count_20genes.png",sep="_"),width=7 ,height=14, units="in",res = 600 ) 
DoHeatmap(myh.cells.myhClust, features = top20$gene)
dev.off()


top100 <- myh.cells.myhClust.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)

write.table(top100,file=paste(name,"liste_top100_clusterMYH.tsv",sep="_"),sep="\t")










#Heatmap of the top 2 genes differential expressed in the different clusters (sup 1e)


	data.currentClusters <- mergeData.integrated

	new.cluster.ids <- c("Myonuclei", 
		"Myonuclei",
		"Myonuclei",
		"Myonuclei",
		"FAPS",
		"Myonuclei",
		"FAPS", 
		"Myonuclei",
		"FAPS",
		"Slow Myonuclei",
		"MuSC",
		"Tenocyte",
		"Smooth muscular",
		"Endothelial",
		"B/T Cells",
		"Adipocyte", 
		"MTJ Myonuclei",
		"Pericyte",
		"NMJ Myonuclei")
	names(new.cluster.ids) <- levels(data.currentClusters)
	data.currentClusters <- RenameIdents(data.currentClusters, new.cluster.ids)
	
	data.currentClusters <- ScaleData(object = data.currentClusters, verbose = FALSE)
	
	DefaultAssay(data.currentClusters) <- "RNA"
	
	data.currentClusters.markers <- FindAllMarkers(data.currentClusters, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	
	top2 <- data.currentClusters.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
	
	png(filename=paste(name,"Heatmap_Full-Cluster_2genes.png",sep="_"),width=7 ,height=14, units="in",res = 600 ) 
	DoHeatmap(data.currentClusters, features = top2$gene)
	dev.off()





# Coexpression of metabolic and Myh gene in myonuclei (fig5)	



	data.currentClusters <- mergeData.integrated
	new.cluster.ids <- c("Myonuclei", 
		"Myonuclei",
		"Myonuclei",
		"Myonuclei",
		"FAPS-1",
		"Myonuclei",
		"FAPS-2", 
		"Myonuclei",
		"FAPS-3",
		"Slow Myonuclei",
		"MuSC",
		"Tenocyte",
		"Smooth muscular",
		"Endothelial",
		"B/T Cells",
		"Adipocyte", 
		"MTJ Myonuclei",
		"Pericyte",
		"NMJ Myonuclei")
	names(new.cluster.ids) <- levels(data.currentClusters)
	data.currentClusters <- RenameIdents(data.currentClusters, new.cluster.ids)
	
	
	myh.cells <- subset(data.currentClusters, idents = c("Myonuclei", "Slow Myonuclei"))
png(filename=paste(name,"FeatureScatter_Myh1vsAldoa.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh1", feature2 = "Aldoa", slot = "counts")
dev.off()	
	
png(filename=paste(name,"FeatureScatter_Myh2vsAldoa.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh2", feature2 = "Aldoa", slot = "counts")
dev.off()	
	
png(filename=paste(name,"FeatureScatter_Myh4vsAldoa.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh4", feature2 = "Aldoa", slot = "counts")
dev.off()	
	
png(filename=paste(name,"FeatureScatter_Myh7vsAldoa.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh7", feature2 = "Aldoa", slot = "counts")
dev.off()	
	
	myh.cells <- subset(data.currentClusters, idents = c("Myonuclei", "Slow Myonuclei"))
png(filename=paste(name,"FeatureScatter_Myh1vsIdh2.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh1", feature2 = "Idh2", slot = "counts")
dev.off()	
	
png(filename=paste(name,"FeatureScatter_Myh2vsIdh2.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh2", feature2 = "Idh2", slot = "counts")
dev.off()	
	
png(filename=paste(name,"FeatureScatter_Myh4vsIdh2.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh4", feature2 = "Idh2", slot = "counts")
dev.off()	
	
png(filename=paste(name,"FeatureScatter_Myh7vsIdh2.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh7", feature2 = "Idh2", slot = "counts")
dev.off()	
	




sessionInfo(package = NULL)
R version 3.6.1 (2019-07-05)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] fr_FR.UTF-8/fr_FR.UTF-8/fr_FR.UTF-8/C/fr_FR.UTF-8/fr_FR.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] cowplot_1.0.0 ggplot2_3.3.2 Seurat_3.1.5  dplyr_1.0.0  

loaded via a namespace (and not attached):
 [1] httr_1.4.1         tidyr_1.1.0        jsonlite_1.6.1     viridisLite_0.3.0  splines_3.6.1      leiden_0.3.3      
 [7] assertthat_0.2.1   ggrepel_0.8.2      globals_0.12.5     pillar_1.4.4       lattice_0.20-41    glue_1.4.1        
[13] reticulate_1.16    digest_0.6.25      RColorBrewer_1.1-2 colorspace_1.4-1   htmltools_0.5.0    Matrix_1.2-18     
[19] plyr_1.8.6         pkgconfig_2.0.3    tsne_0.1-3         listenv_0.8.0      purrr_0.3.4        patchwork_1.0.1   
[25] scales_1.1.1       RANN_2.6.1         Rtsne_0.15         tibble_3.0.1       generics_0.0.2     ellipsis_0.3.1    
[31] withr_2.2.0        ROCR_1.0-11        pbapply_1.4-2      lazyeval_0.2.2     cli_2.0.2          survival_3.2-3    
[37] magrittr_1.5       crayon_1.3.4       future_1.17.0      fansi_0.4.1        nlme_3.1-148       MASS_7.3-51.6     
[43] ica_1.0-2          tools_3.6.1        fitdistrplus_1.1-1 data.table_1.12.8  lifecycle_0.2.0    stringr_1.4.0     
[49] plotly_4.9.2.1     munsell_0.5.0      cluster_2.1.0      irlba_2.3.3        compiler_3.6.1     rsvd_1.0.3        
[55] rlang_0.4.6        grid_3.6.1         ggridges_0.5.2     rstudioapi_0.11    RcppAnnoy_0.0.16   rappdirs_0.3.1    
[61] htmlwidgets_1.5.1  igraph_1.2.5       gtable_0.3.0       codetools_0.2-16   reshape2_1.4.4     R6_2.4.1          
[67] gridExtra_2.3      zoo_1.8-8          future.apply_1.5.0 uwot_0.1.8         KernSmooth_2.23-17 ape_5.4           
[73] stringi_1.4.6      parallel_3.6.1     Rcpp_1.0.4.6       vctrs_0.3.1        sctransform_0.2.1  png_0.1-7         
[79] tidyselect_1.1.0   lmtest_0.9-37     
