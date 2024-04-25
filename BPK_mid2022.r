library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(monocle3)
library(RColorBrewer)
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(metap)
library(tibble)
library(AnnotationHub)
library(ensembldb)
library(DESeq2)
library(metap)
library(multtest)
library("devtools")
library(infercnv)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(scran)
library(scDblFinder)
library(harmony)
library(SeuratWrappers)
library(harmony)


bpk2022 = readRDS(file = "~/OneDrive - University of Massachusetts Boston/Ph.D. Work/Shared_SPJR/scRNA-seq/RDS_obj/bpk12.rds")
s2 = readRDS(file = "~/OneDrive - University of Massachusetts Boston/Ph.D. Work/Shared_SPJR/scRNA-seq/RDS_obj/s2.rds")

DimPlot(bpk2022, label = TRUE)
DimPlot(s2, label = TRUE)

bpk2022.sub = subset(bpk2022, idents = c("B","LP","L"))
s2.sub = subset(s2, idents = c("B","ML","LP","HSL","TPCs"))

DimPlot(s2.sub)

bpk.mid = merge(x = s2.sub, y = bpk2022.sub)

bpk.mid[["percent.mt"]] = PercentageFeatureSet(bpk.mid, pattern = "^mt-")
    head(bpk.mid[["percent.mt"]])
VlnPlot(bpk.mid, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3, pt.size = 0, group.by = "orig.ident")

plot1 = FeatureScatter(bpk.mid, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
plot1

plot2 = FeatureScatter(bpk.mid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
plot2

bpk.mid = subset(bpk.mid, subset = nFeature_RNA>200 & nFeature_RNA< 7500 & percent.mt <5)
bpk.mid = NormalizeData(bpk.mid, normalization.method = "LogNormalize", scale.factor = 1e4)

bpk.mid = FindVariableFeatures(bpk.mid, selection.method = "vst", nfeatures = 2250)
    bpk.mid.50 = head(VariableFeatures(bpk.mid),50)
    bpk.mid.plot = VariableFeaturePlot(bpk.mid)
    bpk.mid.top50.p = LabelPoints(plot = bpk.mid.plot, points = bpk.mid.50, repel= T)
all.genes = rownames(bpk.mid)
bpk.mid = ScaleData(bpk.mid, features = all.genes)


bpk.mid = RunPCA(object = bpk.mid, features = VariableFeatures(bpk.mid), genes.print = 10)
ElbowPlot(bpk.mid)

bpk.mid = FindNeighbors(bpk.mid, dims = 1:15)
bpk.mid = FindClusters(bpk.mid, resolution = 0.3)
bpk.mid = RunUMAP(bpk.mid, dims = 1:10) #range
    p1 = DimPlot(bpk.mid, reduction = "umap", split.by = "orig.ident", pt.size = .5)
    p2 = DimPlot(bpk.mid, reduction = "umap", label = TRUE)
p1+p2 

DimPlot(bpk.mid, reduction = "umap", group.by = "orig.ident", pt.size = .5)

FeaturePlot(bpk.mid, features = c("Krt14","Krt8", "Wfdc18"))


DimPlot(bpk.mid, label = T)
FeaturePlot(bpk.mid, features = c("Krt14","Lgals7","Wfdc18"))



x = FindMarkers(bpk.mid, ident.1 = 1, ident.2 = 5, min.pct = 0.25)
write.table(x, file = "~/Desktop/basal_s2vbpk2022.xls", quote = FALSE, sep="\t")

head(x)

library(tibble)
rownames_to_column(x, var = "genes")
colnames(x) = c("genes",'p_val','avg_log2FC','pct.1',"pct.2","p_val_adj")

BiocManager::install("harmony")

harmony::RunHarmony(bpk.mid)

?harmony
