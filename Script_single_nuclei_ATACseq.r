# To analyse snATAC from adult skeletal muscle I follow the vignette from Signac software :
# https://satijalab.org/signac/articles/mouse_brain_vignette.html
# https://satijalab.org/signac/articles/merging.html
# https://satijalab.org/signac/articles/motif_vignette.html
# https://satijalab.org/signac/articles/integration.html
# We merged 2 datasets called MAT and Steph from 2 replicates


library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)


# Pre-processing workflow of the 2 datasets

counts <- Read10X_h5("/Users/matthieu/Desktop/ATAC2/filtered_peak_bc_matrix.h5")

metadata <- read.csv(file = "/Users/matthieu/Desktop/ATAC2/singlecell.csv", header = TRUE, row.names = 1)


MAT <- CreateSeuratObject(counts = counts, assay = 'peaks', project = 'ATAC', min.cells = 1, meta.data = metadata)

fragment.path <- '/Users/matthieu/Desktop/ATAC2/fragments.tsv.gz'


MAT <- SetFragments(object = MAT, file = fragment.path)


counts2 <- Read10X_h5("/Users/matthieu/Desktop/ATAC2/ATAC3_steph/filtered_peak_bc_matrix.h5")


metadata2 <- read.csv(file = "/Users/matthieu/Desktop/ATAC2/ATAC3_steph/singlecell.csv", header = TRUE, row.names = 1)


Steph <- CreateSeuratObject(counts = counts2, assay = 'peaks', project = 'ATAC', min.cells = 1, meta.data = metadata2)


fragment.path2 <- '/Users/matthieu/Desktop/ATAC2/ATAC3_steph/fragments.tsv.gz'


Steph <- SetFragments(object = Steph, file = fragment.path2)



#Computing QC Metrics


MAT <- NucleosomeSignal(object = MAT)
Steph <- NucleosomeSignal(object = Steph)


MAT$pct_reads_in_peaks <- MAT$peak_region_fragments / MAT$passed_filters * 100

MAT$blacklist_ratio <- MAT$blacklist_region_fragments / MAT$peak_region_fragments

VlnPlot(object = MAT, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'), pt.size = 0.1, ncol = 4) + NoLegend()


Steph$pct_reads_in_peaks <- Steph$peak_region_fragments / Steph$passed_filters * 100


Steph$blacklist_ratio <- Steph$blacklist_region_fragments / Steph$peak_region_fragments

VlnPlot(object = Steph, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'), pt.size = 0.1, ncol = 4) + NoLegend()


MAT$nucleosome_group <- ifelse(MAT$nucleosome_signal > 10, 'NS > 10', 'NS < 10')

FragmentHistogram(object = MAT, group.by = 'nucleosome_group', region = 'chr1-1-10000000')


gene.ranges <- genes(EnsDb.Mmusculus.v79)

gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]


tss.ranges <- GRanges(seqnames = seqnames(gene.ranges), ranges = IRanges(start = start(gene.ranges), width = 2), strand = strand(gene.ranges))



seqlevelsStyle(tss.ranges) <- 'UCSC'

tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')

MAT <- TSSEnrichment(object = MAT, tss.positions = tss.ranges[1:2000])

Steph <- TSSEnrichment(object = Steph, tss.positions = tss.ranges[1:2000])

MAT$high.tss <- ifelse(MAT$TSS.enrichment > 2, 'High', 'Low')

TSSPlot(MAT, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()

Steph$high.tss <- ifelse(Steph$TSS.enrichment > 2, 'High', 'Low')

TSSPlot(Steph, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()

MAT <- subset(MAT, subset = peak_region_fragments > 3000 & peak_region_fragments < 100000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.025 & nucleosome_signal < 10 & TSS.enrichment > 2)


Steph <- subset(Steph, subset = peak_region_fragments > 3000 & peak_region_fragments < 100000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.025 & nucleosome_signal < 10 & TSS.enrichment > 2)




#Merging the 2 datasets


combined.peaks <- UnifyPeaks(object.list = list(MAT, Steph), mode = "reduce")

combined.peaks

MAT.counts <- FeatureMatrix(fragments = GetFragments(MAT), features = combined.peaks, sep = c(":", "-"), cells = colnames(MAT))

Steph.counts <- FeatureMatrix(fragments = GetFragments(Steph), features = combined.peaks, sep = c(":", "-"), cells = colnames(Steph))

MAT[['peaks']] <- CreateAssayObject(counts = MAT.counts)

Steph[['peaks']] <- CreateAssayObject(counts = Steph.counts)

MAT$dataset <- 'MAT'

Steph$dataset <- 'Steph'

combined <- merge(x = MAT, y = list(Steph), add.cell.ids = c("M", "S"))

DefaultAssay(combined) <- "peaks" 

#Normalization and linear dimensional reduction

combined <- RunTFIDF(combined)

combined <- FindTopFeatures(combined, min.cutoff = 20)

combined <- RunSVD(combined, reduction.key = 'LSI_', reduction.name = 'lsi', irlba.work = 400)

combined <- RunSVD(combined, reduction.key = 'LSI_', reduction.name = 'lsi', irlba.work = 400)

# Non-linear dimension reduction and clustering

combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')

DimPlot(combined, group.by = 'dataset', pt.size = 0.1)

saveRDS(combined, file = "/Users/matthieu/Desktop/ATAC2/merged/combined.rds")

combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)

combined <- FindClusters(object = combined, algorithm = 3, resolution = 1.2, verbose = FALSE)

DimPlot(object = combined, label = TRUE) + NoLegend()


# Create a gene activity matrix


combined <- SetFragments(combined, "/Users/matthieu/Desktop/ATAC2/fragment_merged/fragments.tsv.gz")

fragment.path <- '/Users/matthieu/Desktop/ATAC2/fragment_merged/fragments.tsv.gz'

gene.coords <- genes(EnsDb.Mmusculus.v79, filter = ~ gene_biotype == "protein_coding")

seqlevelsStyle(gene.coords) <- 'UCSC'

genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')

genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)

gene.activities <- FeatureMatrix(fragments = fragment.path, features = genebodyandpromoter.coords, cells = colnames(combined), chunk = 10)

gene.key <- genebodyandpromoter.coords$gene_name

names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)

rownames(gene.activities) <- make.unique(gene.key[rownames(gene.activities)])

gene.activities <- gene.activities[rownames(gene.activities)!="",]

combined[['RNA']] <- CreateAssayObject(counts = gene.activities)

combined <- NormalizeData(object = combined, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = median(combined$nCount_RNA))



Integrating with snRNA-seq data

library(Signac)
library(Seurat)
library(JASPAR2018)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)


Muscle_RNA <- readRDS("/Users/matthieu/Desktop/Single_rnaseq/mergeDataIntegrated.rds")

Muscle_RNA <- FindVariableFeatures(object = Muscle_RNA, nfeatures = 5000)

transfer.anchors <- FindTransferAnchors(reference = Muscle_RNA, query = combined, reduction = 'cca', dims = 1:40)

predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = Muscle_RNA$seurat_clusters, weight.reduction = combined[['lsi']], dims = 2:30)

combined <- AddMetaData(object = combined, metadata = predicted.labels)


plot1 <- DimPlot(Muscle_RNA, group.by = 'seurat_clusters', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(combined, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
plot1 + plot2



#Find differentially accessible peaks between clusters in subset of myonuclei MuSC and FAPS

combined.subset <- subset(combined, idents = c("Myh4","Myh2+1", "Myh7", "MuSC", "FAPS"))


gene.ranges <- genes(EnsDb.Mmusculus.v79)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

#Fast Myh Locus

CoveragePlot(object = combined.subset, region = "chr11:67148009-67276263", sep = c(":", "-"), annotation = gene.ranges, peaks = StringToGRanges(regions = rownames(combined.subset), sep = c(":", "-")), extend.upstream = 5000, extend.downstream = 5000, ncol = 1)

#Slow Myh Locus

CoveragePlot(object = combined.subset, region = "chr14:54928952-55003114", sep = c(":", "-"), annotation = gene.ranges, peaks = StringToGRanges(regions = rownames(combined.subset), sep = c(":", "-")), extend.upstream = 5000, extend.downstream = 5000, ncol = 1)

#Pvalb Locus

CoveragePlot(object = combined.subset, region = "chr15:78181565-78215953", sep = c(":", "-"), annotation = gene.ranges, peaks = StringToGRanges(regions = rownames(combined.subset), sep = c(":", "-")), extend.upstream = 5000, extend.downstream = 5000, ncol = 1)


#Atp2a2 Locus

CoveragePlot(object = combined.subset, region = "chr5:122433397-122525168", sep = c(":", "-"), annotation = gene.ranges, peaks = StringToGRanges(regions = rownames(combined.subset), sep = c(":", "-")), extend.upstream = 5000, extend.downstream = 5000, ncol = 1)

#Col22a1 Locus

combined.subset2 <- subset(combined, idents = c("Myh4","Myh2+1", "Myh7", "MTJ", "FAPS"))


CoveragePlot(object = combined.subset2, region = "chr15:71736186-72093835", sep = c(":", "-"), annotation = gene.ranges, peaks = StringToGRanges(regions = rownames(combined.subset2), sep = c(":", "-")), extend.upstream = 5000, extend.downstream = 5000,ncol = 1)



#Heatmap genes with more chromatin accessibility in the differents clusters


library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

DefaultAssay(combined) <- "RNA"

combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top5 <- combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(combined, features = top5$gene) + NoLegend()

top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(combined, features = top10$gene) + NoLegend()

top2 <- combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
DoHeatmap(combined, features = top2$gene) + NoLegend()


#Motif analysis with Signac and chromvar

library(Signac)
library(Seurat)
library(JASPAR2018)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
set.seed(1234)


combined <- RunChromVAR(object = combined, genome = BSgenome.Mmusculus.UCSC.mm10)
DefaultAssay(combined) <- 'chromvar'

differential.activity <- FindMarkers(object = combined, ident.1 = 'Myh4', only.pos = TRUE, test.use = 'LR', latent.vars = 'nCount_peaks')

MotifPlot(object = combined, motifs = head(rownames(differential.activity)), assay = 'peaks')




#Coemmbedding single nuclei RNAseq and single nuclei ATAC seq data with SEURAT

library(Seurat)
library(ggplot2)
library(patchwork)
peaks <- Read10X_h5("/Users/maire/Desktop/ATAC/M4/outs/filtered_peak_bc_matrix.h5")
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = "/Users/maire/Desktop/ATAC/Mus_musculus.GRCm38.93.gtf", 
    seq.levels = c(1:19, "X", "Y"), upstream = 2000, verbose = TRUE)

muscle.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10x_ATAC")
muscle.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
meta <- read.table("/Users/maire/Desktop/ATAC/M4/outs/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
    stringsAsFactors = FALSE)
meta <- meta[colnames(muscle.atac), ]
muscle.atac <- AddMetaData(muscle.atac, metadata = meta)
muscle.atac <- subset(muscle.atac, subset = nCount_ATAC > 5000)
muscle.atac$tech <- "atac"
DefaultAssay(muscle.atac) <- "ACTIVITY"
muscle.atac <- FindVariableFeatures(muscle.atac)
muscle.atac <- NormalizeData(muscle.atac)
muscle.atac <- ScaleData(muscle.atac)
DefaultAssay(muscle.atac) <- "ATAC"
VariableFeatures(muscle.atac) <- names(which(Matrix::rowSums(muscle.atac) > 100))
muscle.atac <- RunLSI(muscle.atac, n = 50, scale.max = NULL)
muscle.atac <- RunUMAP(muscle.atac, reduction = "lsi", dims = 1:50)


#merge of Soleus and quadriceps  snRNAseq 

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)

sol.data <- Read10X(data.dir = "/Users/maire/Desktop/ATAC/s/filtered_feature_bc_matrix")
sol <- CreateSeuratObject(counts = sol.data, project = "Sol")

quad.data <- Read10X(data.dir = "/Users/maire/Desktop/ATAC/q/filtered_feature_bc_matrix")
quad <- CreateSeuratObject(counts = quad.data, project = "quad")

merge.data <- merge(sol, y = quad, add.cell.ids = c("sol", "quad"), project = "SOL_quad")

merge.data[["percent.mt"]] <- PercentageFeatureSet(merge.data, pattern = "^mt-")
mergeData <- subset(merge.data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

mergeData.list <- SplitObject(object = mergeData, split.by = "orig.ident")

for (i in 1:length(x = mergeData.list)) {
    mergeData.list[[i]] <- NormalizeData(object = mergeData.list[[i]], verbose = FALSE)
    mergeData.list[[i]] <- FindVariableFeatures(object = mergeData.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}


mergeData.anchors <- FindIntegrationAnchors(object.list = mergeData.list, dims = 1:30)
mergeData.integrated <- IntegrateData(anchorset = mergeData.anchors, dims = 1:30)

DefaultAssay(mergeData.integrated) <- "integrated"


mergeData.integrated <- ScaleData(object = mergeData.integrated, verbose = FALSE)
mergeData.integrated <- RunPCA(object = mergeData.integrated, npcs = 30, verbose = FALSE)
mergeData.integrated <- RunUMAP(object = mergeData.integrated, reduction = "pca", dims = 1:30)


mergeData.integrated <- FindNeighbors(mergeData.integrated, dims = 1:20)
mergeData.integrated <- FindClusters(mergeData.integrated, resolution = 0.5)

DimPlot(mergeData.integrated, reduction = "umap", label = TRUE)

FeaturePlot(mergeData.integrated, features = c("Myh1","Myh2","Myh4","Myh7","Ttn","Pdgfra","Pax7","Pparg","Pecam1","Ptprc","Myh11","Col11a1","Col24a1","Colq"))

new.cluster.ids <- c("Fast Myh", 
		"Fast Myh",
		"FAPS",
		"Slow Myh7",
		"FAPS",
		"FAPS",
		"Tenocyte", 
		"MuSC",
		"Endothelial",
		"B/T Cells",
		"Smooth muscular ",
		"MTJ",
		"Adipocyte",
		"Pericyte",
		"NMJ")
	names(new.cluster.ids) <- levels(mergeData.integrated)
	mergeData.integrated <- RenameIdents(mergeData.integrated, new.cluster.ids)


mergeData.integrated$tech <- "rna"

transfer.anchors <- FindTransferAnchors(reference = mergeData.integrated, query = muscle.atac, features = VariableFeatures(object = mergeData.integrated), reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")



celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = mergeData.integrated$integrated_snn_res.0.5, weight.reduction = muscle.atac[["lsi"]])


muscle.atac <- AddMetaData(muscle.atac, metadata = celltype.predictions)




muscle.atac.filtered <- subset(muscle.atac, subset = prediction.score.max > 0.5)

muscle.atac.filtered$predicted.id <- factor(muscle.atac.filtered$predicted.id, levels = levels(mergeData.integrated))



genes.use <- VariableFeatures(mergeData.integrated)


refdata <- GetAssayData(mergeData.integrated, assay = "integrated", slot = "data")[genes.use, ]

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = muscle.atac[["lsi"]])

muscle.atac[["integrated"]] <- imputation

coembed <- merge(x = mergeData.integrated, y = muscle.atac)

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)


coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)


coembed <- RunUMAP(coembed, dims = 1:30)

coembed$seurat_clusters <- ifelse(!is.na(coembed$seurat_clusters), coembed$seurat_clusters, coembed$predicted.id)

p1 <- DimPlot(coembed, group.by = "tech")

p2 <- DimPlot(coembed, group.by = "seurat_clusters", label = TRUE, repel = TRUE)


p1 + p2










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
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] patchwork_1.0.1            EnsDb.Mmusculus.v79_2.99.0 ensembldb_2.10.2           AnnotationFilter_1.10.0   
 [5] GenomicFeatures_1.38.2     AnnotationDbi_1.48.0       Biobase_2.46.0             GenomicRanges_1.38.0      
 [9] GenomeInfoDb_1.22.1        IRanges_2.20.2             S4Vectors_0.24.4           BiocGenerics_0.32.0       
[13] Signac_0.2.5               cowplot_1.0.0              ggplot2_3.3.2              Seurat_3.1.5              
[17] dplyr_1.0.0               

loaded via a namespace (and not attached):
  [1] Rtsne_0.15                  colorspace_1.4-1            ellipsis_0.3.1              ggridges_0.5.2             
  [5] XVector_0.26.0              rstudioapi_0.11             leiden_0.3.3                listenv_0.8.0              
  [9] ggrepel_0.8.2               ggfittext_0.9.0             bit64_0.9-7                 fansi_0.4.1                
 [13] codetools_0.2-16            splines_3.6.1               jsonlite_1.6.1              Rsamtools_2.2.3            
 [17] ica_1.0-2                   dbplyr_1.4.4                cluster_2.1.0               png_0.1-7                  
 [21] uwot_0.1.8                  sctransform_0.2.1           compiler_3.6.1              httr_1.4.1                 
 [25] assertthat_0.2.1            Matrix_1.2-18               lazyeval_0.2.2              cli_2.0.2                  
 [29] htmltools_0.5.0             prettyunits_1.1.1           tools_3.6.1                 rsvd_1.0.3                 
 [33] igraph_1.2.5                gtable_0.3.0                glue_1.4.1                  GenomeInfoDbData_1.2.2     
 [37] RANN_2.6.1                  reshape2_1.4.4              rappdirs_0.3.1              Rcpp_1.0.4.6               
 [41] vctrs_0.3.1                 Biostrings_2.54.0           ape_5.4                     nlme_3.1-148               
 [45] rtracklayer_1.46.0          ggseqlogo_0.1               lmtest_0.9-37               stringr_1.4.0              
 [49] globals_0.12.5              lifecycle_0.2.0             irlba_2.3.3                 XML_3.99-0.3               
 [53] future_1.17.0               MASS_7.3-51.6               zlibbioc_1.32.0             zoo_1.8-8                  
 [57] scales_1.1.1                ProtGenerics_1.18.0         hms_0.5.3                   SummarizedExperiment_1.16.1
 [61] RColorBrewer_1.1-2          gggenes_0.4.0               curl_4.3                    memoise_1.1.0              
 [65] reticulate_1.16             pbapply_1.4-2               gridExtra_2.3               biomaRt_2.42.1             
 [69] stringi_1.4.6               RSQLite_2.2.0               BiocParallel_1.20.1         matrixStats_0.56.0         
 [73] rlang_0.4.6                 pkgconfig_2.0.3             bitops_1.0-6                lattice_0.20-41            
 [77] ROCR_1.0-11                 purrr_0.3.4                 GenomicAlignments_1.22.1    htmlwidgets_1.5.1          
 [81] bit_1.1-15.2                tidyselect_1.1.0            RcppAnnoy_0.0.16            plyr_1.8.6                 
 [85] magrittr_1.5                R6_2.4.1                    generics_0.0.2              DelayedArray_0.12.3        
 [89] DBI_1.1.0                   pillar_1.4.4                withr_2.2.0                 fitdistrplus_1.1-1         
 [93] survival_3.2-3              RCurl_1.98-1.2              tibble_3.0.1                future.apply_1.5.0         
 [97] tsne_0.1-3                  crayon_1.3.4                KernSmooth_2.23-17          BiocFileCache_1.10.2       
[101] plotly_4.9.2.1              progress_1.2.2              grid_3.6.1                  data.table_1.12.8          
[105] blob_1.2.1                  digest_0.6.25               tidyr_1.1.0                 openssl_1.4.1              
[109] munsell_0.5.0               viridisLite_0.3.0           askpass_1.1                

