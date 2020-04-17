library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)


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

# Figures 

mergeData.integrated <- FindNeighbors(mergeData.integrated, dims = 1:20)
mergeData.integrated <- FindClusters(mergeData.integrated, resolution = 0.5)

# Figures

# UMAP: légendes par samples ou cluster
	DimPlot(object = mergeData.integrated, reduction = "umap", group.by = "orig.ident", label = FALSE)
	DimPlot(mergeData.integrated, reduction = "umap", label = TRUE)

# clusters Names

	new.cluster.ids <- c("Myh1+2 Cells -1", 
		"Myh4 Cells -1",
		"Myh1+2 Cells -2",
		"Myh4 Cells -2",
		"FAPS-1",
		"Myh4 Cells -3",
		"FAPS-2", 
		"Unknow myonuclei",
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

# UMAP with cluster names
png(filename=paste(name,"UMAP.png",sep="_"),width=10 ,height=10, units="in",res = 600 ) 
	DimPlot(mergeData.integrated, reduction = "umap", label = FALSE)
dev.off()

png(filename=paste(name,"UMAP_withClusterNames.png",sep="_"),width=10 ,height=10, units="in",res = 600 ) 
	DimPlot(mergeData.integrated, reduction = "umap", label = TRUE)
dev.off()

png(filename=paste(name,"UMAP_byOrig.png",sep="_"),width=10 ,height=10, units="in",res = 600 ) 
	DimPlot(mergeData.integrated, reduction = "umap", label = FALSE, group.by = "orig.ident")
dev.off()
   
# QC
	VlnPlot(mergeData.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, nrow = 4, split.by = "orig.ident")


###########################
# Differential expression###########
###########################

	DefaultAssay(mergeData.integrated) <- "RNA"

# Nouvelle colonne cluster_sample
	mergeData.integrated$cluster.sample <- paste(mergeData.integrated$orig.ident, mergeData.integrated$seurat_clusters, sep = "_")


# Marqueurs spécifiques de chaque cluster face au reste des noyaux, et commun aux différentes conditions

	Idents(mergeData.integrated) <- "seurat_clusters"

	ncluster <- length(levels(mergeData.integrated$seurat_clusters)) - 1
	for (i in 0:ncluster) {
		conservedMarkers <- FindConservedMarkers(mergeData.integrated, ident.1 = i, grouping.var = "orig.ident", verbose = FALSE)
		currentCluster.Name <- new.cluster.ids[i+1]
		nameFile <- paste(name,"FindConservedMarkers_Cluster",currentCluster.Name,".txt",sep="_")
		write.table(conservedMarkers,file=nameFile,sep="\t")
	}
png(filename=paste(name,"FeaturePlot.png",sep=""),width=10 ,height=7, units="in",res = 600 ) 
	FeaturePlot(mergeData.integrated, features = c("Myh1","Myh2","Myh4","Myh7","Ttn","Pdgfra","Pax7","Pparg","Pecam1","Ptprc","Myh11","Col11a1","Col24a1","Colq"))
dev.off()

# Travail Marqueurs Myh[X] vs Cluster

NMJ.markers <- FindMarkers(mergeData.integrated, ident.1 = "NMJ Myonuclei", ident.2 = c("Myh1+2 Cells -1","Myh4 Cells -1","Myh1+2 Cells -2","Myh4 Cells -2","Myh4 Cells -3","Myh7 Cells"), min.pct = 0.25)
Newly_fused.markers <- FindMarkers(mergeData.integrated, ident.1 = "Newly fused myonuclei", ident.2 = c("Myh1+2 Cells -1","Myh4 Cells -1","Myh1+2 Cells -2","Myh4 Cells -2","Myh4 Cells -3","Myh7 Cells"), min.pct = 0.25)

write.table(MTJ.markers,file=paste(name,"liste_clusterMTJ_vs_allMYH.tsv",sep="_"),sep="\t")
write.table(NMJ.markers,file=paste(name,"liste_clusterNMJ_vs_allMYH.tsv",sep="_"),sep="\t")
write.table(Unknow myonuclei,file=paste(name,"liste_clusterUnknow_myonuclei_vs_allMYH.tsv",sep="_"),sep="\t")

########################################################################


	data.currentClusters <- mergeData.integrated

	new.cluster.ids <- c("Myonuclei", 
		"Myonuclei",
		"Myonuclei",
		"Myonuclei",
		"FAPS-1",
		"Myonuclei",
		"FAPS-2", 
		"Unknow myonuclei",
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
dev.off()


##########################################"
########################################################################


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


##########################################"
########################################################################


	data.currentClusters <- mergeData.integrated

	new.cluster.ids <- c("Myonuclei", 
		"Myonuclei",
		"Myonuclei",
		"Myonuclei",
		"FAPS-1",
		"Myonuclei",
		"FAPS-2", 
		"Unknow myonuclei myonuclei",
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
	
	Unknow myonuclei <- FindMarkers(data.currentClusters, ident.1 = "Unknow myonuclei", ident.2 = "Myonuclei", min.pct = 0.25, only.pos = TRUE)
	Unknow myonuclei2 <- FindMarkers(data.currentClusters, ident.2 = "Unknow myonuclei", ident.1 = "Myonuclei", min.pct = 0.25, only.pos = TRUE)

	
data.currentClusters <- subset(data.currentClusters, idents = c("Unknow myonuclei","Myonuclei"))
data.currentClusters <- ScaleData(object = data.currentClusters, verbose = FALSE)
	
top30 <- head(Unknow myonuclei.markers, n=30)
top30.2 <- head(Unknow myonuclei.markers2, n=30)
png(filename=paste(name,"_Heatmap_Unknow myonuclei-vs-Myonuclei_30genes.png",sep="_"),width=7 ,height=14, units="in",res = 600 ) 
DoHeatmap(data.currentClusters, features = c(rownames(top30),rownames(top30.2)))
dev.off()


##########################################"




# Cells number for each clusters and samples

	tableData.nCells <- table(mergeData.integrated@active.ident, mergeData.integrated@meta.data$orig.ident)
	tableNOrig <- paste(name,"nCells_byCluster_bySample.txt",sep="_")
	write.table(tableData.nCells,file=tableNOrig,sep="\t")

# Travail sur les Myh[X]
	
	# Figure expression globale des myh
	myh.features <- c("Myh1","Myh2","Myh4","Myh7")
	FeaturePlot(mergeData.integrated, features = myh.features)
	RidgePlot(mergeData.integrated, features = myh.features, ncol = 2)
	
	# table expression Myh
	myh.cells <- subset(mergeData.integrated, idents = c("Myh2 Cells","Myh7 Cells","Myh1 Cells","Myh4 Cells"))
	table(myh.cells@active.ident, myh.cells@meta.data$orig.ident)
	
	# Figures comparaison d'expression des myh avec selection préalable des clusters
	myh.cells <- subset(mergeData.integrated, idents = c("Myh2 Cells","Myh7 Cells"))
png(filename=paste(name,"FeatureScatter_Myh7vsMyh2.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh7", feature2 = "Myh2", slot = "counts")
dev.off()	
	
	myh.cells <- subset(mergeData.integrated, idents = c("Myh2 Cells","Myh1 Cells"))
png(filename=paste(name,"FeatureScatter_Myh1vsMyh2.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh1", feature2 = "Myh2", slot = "counts")
dev.off()	
	
	
	myh.cells <- subset(mergeData.integrated, idents = c("Myh4 Cells -2","Myh1 Cells"))
png(filename=paste(name,"FeatureScatter_Myh1vsMyh4.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Myh1", feature2 = "Myh4", slot = "counts")
dev.off()	
	
	
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
png(filename=paste(name,"FeatureScatter_Idh2vsAldoa.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	FeatureScatter(myh.cells, feature1 = "Idh2", feature2 = "Aldoa", slot = "counts")
dev.off()	
	
	
	
#~ 	myh.cells <- subset(mergeData.integrated, features = c("Myh1", "Myh2", "Myh4", "Myh7"))


myh.cells <- subset(mergeData.integrated, idents = c("0","1","2","3","5","7","9","16","18"))
DefaultAssay(myh.cells) <- "RNA"
png(filename=paste(name,"_selectedCluster_DotPlot_Myh7-Myh2-Myh1-Myh4-Ttn.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
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



myh.cells <- subset(data.currentClusters, idents = c("Myonuclei", "Slow Myonuclei", "MTJ Myonuclei", "NMJ Myonuclei"))
DefaultAssay(myh.cells) <- "RNA"
png(filename=paste(name,"selectedCluster_DotPlot_Myh7-Myh2-Myh1-Myh4-Ttn_no-NewMyonuclei.png",sep="_"),width=7 ,height=7, units="in",res = 600 ) 
	DotPlot(myh.cells, features = c("Myh1","Myh2","Myh4","Myh7","Ttn"))
dev.off()


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
	

