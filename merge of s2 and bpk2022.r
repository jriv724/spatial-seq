#merge of s2 and bpk2022
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(monocle3)
library(ArchR)
library(RColorBrewer)
library(SeuratWrappers)
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
library("funcutils")
library(ggbeeswarm)
#BiocManager::install("DESeq2")
bpk2022 = readRDS(file = "~/OneDrive - University of Massachusetts Boston/Ph.D. Work/Shared_SPJR/scRNA-seq/RDS_obj/bpk12.rds")
s2 = readRDS(file = "~/OneDrive - University of Massachusetts Boston/Ph.D. Work/Shared_SPJR/scRNA-seq/RDS_obj/s2.rds")

bpk2022$cells = "mk2"
s2$cells = "mk1"


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
    p1 = DimPlot(bpk.mid, reduction = "umap", group.by = "orig.ident", pt.size = .5)
    p2 = DimPlot(bpk.mid, reduction = "umap", label = TRUE)
p1+p2 

tumor_rsids = subset(s2, ident = c('TPCs'))
x = rownames(tumor_rsids@meta.data)
y = (bpk.mid@meta.data[x,])
barcodes = rownames(y)

DimPlot(bpk.mid, cells.highlight = barcodes, cols.highlight = "firebrick", cols = "gray82", pt.size = .25, label = T)



bpk.mid$log10GenesPerUMI = log10(bpk.mid$nFeature_RNA) / log10(bpk.mid$nCount_RNA)
bpk.mid$mitoRatio = PercentageFeatureSet(object = bpk.mid, pattern = "^MT-")
bpk.mid
mitoRatio = bpk.mid$mitoRatio / 100

meta = bpk.mid@meta.data 
meta$cells = rownames(meta)
meta$sample = NA
meta$sample[which(str_detect(meta$cells, "^ctrl"))] = "ctrl"
meta$sample[which(str_detect(meta$cells, "^4nqo1"))] = "4nqo1"

meta = meta %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
View(meta)

bpk.mid@meta.data = meta
save(bpk.mid, file="~/s4.s2.RData")

#number of cells included in sample
meta %>% 
  ggplot(aes(x=sample, fill = sample)) + 
  geom_bar()+
  theme_classic()+
  theme(axis.text.x =  element_text(angle=45, vjust = 1, hjust = 1))+
  theme(plot.title=element_text(hjust=0.5, face="bold"))+
  ggtitle("Number of Cells")

# Visualize the number UMIs/transcripts per cell to search for outliers
meta %>% 
  	ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 600)

#visualize the distribution of genes detected per cell via histogram 
meta %>% 
  	ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 500)+
    scale_x_continuous(limits=c(300,4000))


# Visualize the distribution of genes detected per cell via boxplot
meta %>% 
  	ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  	geom_boxplot() + 
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells vs NGenes")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
meta %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
    scale_x_continuous(limits=c(100,40000))+
  	facet_wrap(~sample)

# Visualize the distribution of mitochondrial gene expression detected per cell
meta %>% 
  	ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
meta %>%
  	ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8)

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat = subset(x = bpk.mid, 
                         subset= (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (mitoRatio < 0.10))
#Within our data we will have many genes with zero counts. These genes can dramatically reduce the average expression for a cell and so we will remove them from our data. We will start by identifying which genes have a zero count in each cell:
# Extract counts
counts = GetAssayData(bpk.mid, slot = "counts")
# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero = counts > 0 
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes = Matrix::rowSums(nonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts = counts[keep_genes,]
# Reassign to filtered Seurat object
filtered = CreateSeuratObject(filtered_counts, meta.data=filtered_seurat@meta.data)

meta.filt = filtered@meta.data 
meta.filt$cells = rownames(meta.filt)
meta.filt$sample = NA
meta.filt$sample[which(str_detect(meta$cells, "^ctrl"))] = "ctrl"
meta.filt$sample[which(str_detect(meta$cells, "^4nqo1"))] = "4nqo1"

# Visualize the number UMIs/transcripts per cell to search for outliers
meta.filt %>% 
  	ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 600)



#####SCT integration analysis - this method will emphasize removal of replicating cells and control for cell cycle#####
# Download cell cycle genes for organism at https://github.com/hbc/tinyatlas/tree/master/cell_cycle. Read it in with:
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv") 
cell_cycle_genes <- read.csv(text = cc_file)
#gene names are in ENSBL IDs we will pull out actual gene names
#annotation hub
ah = AnnotationHub()
# Access the Ensembl database for organism
ahDb = query(ah,
             pattern = "Mus musculus", "EnsDb",
             ignore.case = TRUE)

id = ahDb %>%
      mcols() %>%
      rownames() %>%
      tail(n = 1)

# Download the appropriate Ensembldb database
edb = ah[[id]]
  #str(edb)
  #View(edb)

# Extract gene-level information from database
annotations = genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations = annotations %>%
        dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

#make compendium of gene names to match against ENSMBL IDs
cc.markers = dplyr::left_join(cell_cycle_genes, annotations,by=c("geneID"="gene_id"))

#S-Phase genes
sphase = cc.markers %>% 
        dplyr::filter(phase =="S") %>%
        pull("gene_name")
#G2/M genes
g2m = cc.markers %>% 
        dplyr::filter(phase =="G2/M") %>%
        pull("gene_name")


write.table(sphase, file = "/Users/joshuarivera/Desktop/sphase_genelist.xls", sep="\t", quote = FALSE,
            col.names=NA)
write.table(g2m, file = "/Users/joshuarivera/Desktop/g2m_genelist.xls", sep="\t", quote = FALSE,
            col.names=NA)

####regulating cell cycles for S2
s2.cc = CellCycleScoring(s2,
                      g2m.features = g2m,
                      s.features = sphase)
s4.cc = CellCycleScoring(s4,
                      g2m.features = g2m,
                      s.features = sphase)

s2.cc = RunPCA(s2.cc)



DimPlot(s2, reduction = "umap", group.by = "Phase")
DimPlot(s4, reduction = "umap", group.by = "Phase")
DimPlot(s4, reduction = "umap")


s2.cc.all.genes = row.names(s2)
s2.cc = ScaleData(s2, vars.to.regress = c(sphase,g2m), features = s2.cc.all.genes)
s2.cc = RunPCA(s2.cc, features = VariableFeatures(s2), genes.print=10)

s2.cc = FindNeighbors(s2.cc, dims=1:15)
s2.cc = FindClusters(s2.cc, resolution=.75)
s2.cc = RunUMAP(s2.cc, dims=1:15)
DimPlot(s2.cc, reduction = "umap", pt.size = 1, label = TRUE, label.size = 6)
DimPlot(s2.cc,
        reduction = 'umap',
        group.by="Phase")

DimPlot(s4.cc,
        reduction = 'umap',
        group.by="Phase")
#0 - luminal 
#1 - luminal
#2 - stroma
#3 - basal
#4 - 
#5
#6
#7

FeaturePlot(s2.cc, features = c("Krt8","Krt14","Wfdc18","Mki67","Lyz2","Vim","Prlr","Top2a"), ncol = 4, cols = c("gold2","red"))
FeaturePlot(s2.cc, features = c("Mki67","Top2a"), ncol = 2, cols = c("gold2","red"), label=T)
FeaturePlot(s2.cc, features = c("Krt8","Krt14","Prlr","Mfge8","Elf5","Sox9","Aldh1a3","Itgb3"), ncol = 4, cols = c("gold2","red"), label=T)
FeaturePlot(s4.cc, features = c("Krt8","Krt14","Prlr","Mfge8","Elf5","Sox9","Aldh1a3","Itgb3"), ncol = 4, cols = c("gold2","red"), label=T)

DotPlot(s2.cc, features = c("Krt8","Krt14","Wfdc18","Mki67","Lyz2","Vim","Prlr"),
        cols = c("yellow","red"))

color = scale_color_brewer(palette = "RdYlBu")


          
#######regulating cell cycle for S4-PBS
s4.cc = CellCycleScoring(s4,
                      g2m.features = g2m,
                      s.features = sphase)

s4.cc = RunPCA(s4.cc)

DimPlot(s4.cc,
        reduction = 'pca',
        group.by="Phase")

DimPlot(s4.cc, reduction = "umap", group.by = "Phase")


s4.cc.all.genes = row.names(s4)
s4.cc = ScaleData(s4, vars.to.regress = c(sphase,g2m), features = s4.cc.all.genes)
s4.cc = RunPCA(s2.cc, features = VariableFeatures(s4), genes.print=10)

s4.cc = FindNeighbors(s4.cc, dims=1:15)
s4.cc = FindClusters(s4.cc, resolution=.2)
s4.cc = RunUMAP(s4.cc, dims=1:15)
DimPlot(s4.cc, reduction = "umap", pt.size = 1, label = TRUE, label.size = 6)

sum(table(s4.cc$seurat_clusters))

FeaturePlot(s4.cc, features = c("Krt8","Krt14","Wfdc18","Mki67","Lyz2","Vim","Prlr"), ncol = 3)


########clasic seurat #######
s4.s2[["percent.mt"]] = PercentageFeatureSet(s4.s2, pattern = "^mt-")
s4.s2[["percent.mt"]]

VlnPlot(s4.s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

plot1 = FeatureScatter(s4.s2, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot1
plot2 = FeatureScatter(s4.s2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

s4.s2 = subset(s4.s2, subset = nFeature_RNA > 200 & nFeature_RNA < 8500 & percent.mt < 5)
s4.s2 = NormalizeData(s4.s2)

s4.s2 = FindVariableFeatures(s4.s2, selection.method = "vst", nfeatures=  2000)
s4.s2.top10 = head(VariableFeatures(s4.s2),20)
s4.s2.plot = VariableFeaturePlot(s4.s2)
s4.s2.top10.p = LabelPoints(plot = s4.s2.plot, points = s4.s2.top10, repel = T, ynudge = 0, xnudge=0)

all.genes = row.names(s4.s2)
s4.s2 = ScaleData(s4.s2, features = all.genes)


s4.s2 = RunPCA(object = s4.s2, features = VariableFeatures(object = s4.s2), genes.print = 10)
ElbowPlot(s4.s2)
print(s4.s2[["pca"]], dims = 1:8, nfeatures = 10)
DimPlot(s4.s2, reduction = "pca")

s4.s2@meta.data$ <- c(rep("V", ncol(stim.sparse)), rep("CTRL", ncol(ctrl.sparse)))

s4.s2 <- RunHarmony(s4.s2, "dataset")
s4.s2 <- RunUMAP(s4.s2, reduction = "harmony")

#DoHeatmap(object=s4.s2,group.by = "orig.ident")
s4.s2 = FindNeighbors(s4.s2, dims= 1:8)
s4.s2 = FindClusters(s4.s2, resolution = .3)
s4.s2 = RunUMAP(s4.s2, dims=1:8) #20

p1 = DimPlot(s4.s2, reduction = 'umap', label = T, split.by = "orig.ident")
p2 = DimPlot(s4.s2, reduction = 'umap', label = T)

p1/p2

FeaturePlot(s4.s2, features = "Ly6a")
#0-luminal
#1-basal
#2-stroma
#3-stroma
#4-luminal
#5-endothelial
#6-luminal
#7-immune
#8-stroma
#9-stroma
#10-basal
#11-immune
#12-luminal
#13-doublets
#14-doublets

s4.s2.current.cluster.ids = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)
s4.s2.new.cluster.ids =c("HSL","B","Fib","ML","LP","End","LP","Im","P","RST","Im","Im","dblt","dblt")
names(s4.s2.new.cluster.ids) = levels(s4.s2)     
s4.s2 = RenameIdents(s4.s2, s4.s2.new.cluster.ids)
DimPlot(s4.s2, reduction="umap", label=TRUE, pt.size=.75, label.size = 4)


subs4.s2 = subset(s4.s2, idents = c("LP","B","RST","HSL","ML"))
DimPlot(subs4.s2, reduction = 'umap', label=T, label.size = 4, pt.size = .75)
DimPlot(subs4.s2, reduction = 'umap', label=T, label.size = 4, pt.size = .75, group.by = "orig.ident")

s2.s4.umap = cbind("Barcode" = rownames(Embeddings(object = s4.s2, reduction = "umap")), Embeddings(object = s4.s2, reduction = "umap"))
write.table(s2.s4.umap, file="/Users/joshuarivera/Desktop/s2vs4/s2vs4_coords.csv", sep = ",", quote = F, row.names = F, col.names = T)


pos.cytokines = c("Arid5a","Stat1","Hspd1","Ncl","Hdac2","Ddt","Mif","Bsg","Xbp1","Pik3r1","Clu","Cd200","Mapk13","Gpsm3","H2-T23","Egr1","Cd14",
                "Cd74","Il17b","Gata3","Thbs1","B2m","Fermt1","Cebpb","F3","Bcl10","Park7","Hspb1","Hmgb1","Gapdh","Pycard","Hras","Irf7","Hmgb2","Cx3cl1","Cadm1")


subs4.s2 = AddModuleScore(subs4.s2, features = list(pos.cytokines), name = "pos.cytokines")

names(subs4.s2[[]])

FeaturePlot(object = subs4.s2, features = "pos.cytokines1") +
            scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
 
subs4.s2.meta = as.data.frame(subs4.s2$seurat_clusters)
subs4.s2.meta$orig.ident = as.factor(subs4.s2$orig.ident)
colnames(subs4.s2.meta) = c("cluster", "orig.ident")
subs4.s2.meta$apop = as.numeric(subs4.s2$pos.cytokines1)
subs4.s2.meta$cytokine = as.numeric(subs4.s2$pos.cytokines1)
head(subs4.s2.meta)
#rm(subs4.s2.meta)

subs4.s2.meta = as_tibble(subs4.s2.meta)

subs4.s2.meta = subs4.s2.meta.tib[subs4.s2.meta$cluster %in% "rst",] 


ggplot(subs4.s2.meta, aes(x = orig.ident, y = apop, fill =orig.ident ))+
  geom_violin(trim=FALSE)+
  stat_summary(fun=mean, geom="point", size = .5) +
  stat_summary(fun.data="mean_sdl",
               geom= "pointrange") +
  scale_fill_manual(values=c("#F37B7D", "#B9C558", "#56B4E9"))+
  geom_quasirandom()+
  theme_classic()


s4.s2$orig.ident
#########performing integration of scRNA-seq datasets######https://satijalab.org/seurat/articles/integration_introduction.html
#BPK treated vs BPK vehicle mid-point mammary gland
drug.comb = SplitObject(s4.s2, split.by = "orig.ident") 

#function to normalize and identify variable features performing variance stabilizing transformation from each dataset independently
drug.comb = lapply(X= drug.comb, FUN = function(x){
  x = NormalizeData(x)
  x = FindVariableFeatures(x, selected.method = "vst", nfeatures = 2000)
})

str(drug.comb)
#### select features that are repeatedly variable across datasets for integration 
features = SelectIntegrationFeatures(object.list = drug.comb)

#<b> Perform integration <b>##
#We then identify anchors using the FindIntegrationAnchors() function, which takes a list of Seurat objects as input, 
#and use these anchors to integrate the two datasets together with IntegrateData().
epithelial.anchors = FindIntegrationAnchors(object.list = drug.comb, anchor.features = features)

# this command creates an 'integrated' data assay
treat.combined = IntegrateData(anchorset = epithelial.anchors)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(treat.combined) = "integrated"

treat.combined = ScaleData(treat.combined, verbose = F)
treat.combined = RunPCA(treat.combined, npcs = 30, verbose = F)
treat.combined = RunUMAP(treat.combined, reduction= "pca", dims = 1:18) #15
treat.combined = FindNeighbors(treat.combined, reduction = "pca", dims = 1:30)
treat.combined = FindClusters(treat.combined, resolution = .3)


DimPlot(treat.combined, reduction='umap', split.by = "orig.ident", label = T)
 FeaturePlot(treat.combined, features= c("Prlr","Mfge8","Krt8","Krt14","Lyz2","Csn3"))


 #positive regulation of cytokine production
  #GO:0001819
pos.cytokines = c("Arid5a","Stat1","Hspd1","Ncl","Hdac2","Ddt","Mif","Bsg","Xbp1","Pik3r1","Clu","Cd200","Mapk13","Gpsm3","H2-T23","Egr1","Cd14",
                "Cd74","Il17b","Gata3","Thbs1","B2m","Fermt1","Cebpb","F3","Bcl10","Park7","Hspb1","Hmgb1","Gapdh","Pycard","Hras","Irf7","Hmgb2","Cx3cl1","Cadm1")

treat.combined[["RNA"]]@meta.features = data.frame(row.names = rownames(treat.combined[["RNA"]]))

treat.combined = AddModuleScore(treat.combined@assays$RNA, features = list(pos.cytokines), name = "pos.cytokines")

names(treat.combined[[]])

FeaturePlot(object = treat.combined, features = "pos.cytokines1") +
            scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
 


#comparing s2 v s4
t.test(x = s2.apop.clust8.df, y = s4.apop.clust8.df, alternative = "two.sided", paired = FALSE)
#Welch Two Sample t-test
#data:  s2.apop.clust8.df and s4.apop.clust8.df
#t = 3.2353, df = 78.48, p-value = 0.00178
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
 #0.02702342 0.11346201
#sample estimates:
#mean of x mean of y 
#0.1906305 0.1203878 

 
 ######
#not including LP population
treat.combined = RenameIdents(treat.combined, `0`="luminal", `1` = "basal", `2` = "stroma", 
    `3` = "luminal", `4` = "luminal", `5` = "endothelial", `6` = "luminal", `7` = "immune", `8` = "stroma", `9` = "stroma", 
    `10` = "basal", `11` = "immune", `12` = "stroma", `13` = "luminal", `14` = "stroma",`15`= "luminal/immune doublet")
#Isolating LP population#
treat.combined.lp = RenameIdents(treat.combined, "luminal"="luminal", "basal" = "basal", "stroma" = "stroma", 
    "luminal" = "LP", "luminal" = "LP", "endothelial" = "endothelial", "luminal" = "LP", "immune" = "immune","stroma" = "stroma", "stroma" = "stroma", 
    "basal" = "basal", "immune" = "immune", "stroma" = "stroma", "luminal" = "luminal", "stroma" = "stroma","luminal/immune doublet"= "luminal/immune doublet")


#pooling LP and B
treat.combined = RenameIdents(treat.combined, `0`="luminal", `1` = "basal", `2` = "stroma", 
    `3` = "E", `4` = "luminal", `5` = "endothelial", `6` = "luminal", `7` = "immune", `8` = "stroma", `9` = "stroma", 
    `10` = "E", `11` = "immune", `12` = "stroma", `13` = "luminal", `14` = "stroma",`15`= "luminal/immune doublet")

p1 = DimPlot(treat.combined, reduction = "umap", label = T)
p2 = DimPlot(treat.combined, reduction = 'umap', split.by = "orig.ident")
p1+p2

#for dissertation proposal
VlnPlot(treat.combined, idents = c(3,10,1), split.by = "orig.ident", features = "Ly6a")


FeaturePlot(treat.combined, features = c("Elf5","Cd14","Kit","Aldh1a1",
                             "Aldh1a3", "Foxc1","Foxc2","Itga2",
                             "Itgb3","Tspan8","Cd55","Prlr"), ncol= 5, label = T)

# For performing differential expression after integration, we switch back to the original data
DefaultAssay(treat.combined) = 'RNA'
basal.markers = FindConservedMarkers(treat.combined, ident.1 ="basal", grouping.var= 'orig.ident' )
  head(basal.markers)
    basal.conserved = rownames(basal.markers[1:20,])
FeaturePlot(treat.combined, features = basal.conserved, ncol = 5)



luminal.markers = FindConservedMarkers(treat.combined, ident.1 = "luminal", grouping.var = 'orig.ident')
  head(luminal.markers)  
    luminal.conserved = rownames(luminal.markers[1:20,])
FeaturePlot(treat.combined, features = luminal.conserved, ncol = 5)


####luminal DGE####
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
library(patchwork)
theme_set(theme_cowplot())


###all luminal DGE####
luminal.dge = subset(treat.combined, idents = "luminal")
  Idents(luminal.dge)= "orig.ident"
  DimPlot(luminal.dge, reduction = 'umap')
avg.lum = as.data.frame(log1p(AverageExpression(luminal.dge)$RNA))
  str(avg.lum)
avg.lum$gene = rownames(luminal.dge$RNA)

o.lum = order(avg.lum$midpoint_treated)
avg.lum.o = avg.lum$gene[o.lum]
  list(avg.lum.o[1:10])
  lum.names = avg.lum.o[1:10]
lum.genes.to.label = luminal.markers[1:10]

  
####basal dge####
bsl.dge = subset(treat.combined, idents = "basal")
  Idents(bsl.dge)= "orig.ident"
  DimPlot(bsl.dge, reduction = 'umap')
avg.bsl = as.data.frame(log1p(AverageExpression(bsl.dge)$RNA))
  str(avg.bsl)
avg.bsl$gene = rownames(bsl.dge$RNA)

bsl.genes.to.label = basal.markers[1:10]



o = order(avg.bsl$midpoint_treated)
avg.bsl.o = avg.bsl$gene[o]
  list(avg.bsl.o[1:10])
  bsl.names = avg.bsl.o[1:10]
  
tail(avg.bsl.o)
print(rev(avg.bsl.o))

###LP dge###
lp.dge = subset(treat.combined.lp, idents = "LP")
  Idents(lp.dge)= "orig.ident"
  DimPlot(lp.dge, reduction = 'umap')
avg.lp = as.data.frame(log1p(AverageExpression(lp.dge)$RNA))
  str(avg.lp)
avg.lp$gene = rownames(lp.dge$RNA)
  #avg.lp$gene[1:10]
lp.genes.to.label = avg.lp$gene[1:10]

o.lp = order(avg.lp$s2_treated_het)
avg.lp.o = avg.lp$gene[o.lp]
  list(avg.lp.o[1:10])
  lp.names = avg.lp.o[1:10]

p1 = ggplot(avg.lum, aes(s4_treated_het,midpoint_treated)) + geom_point() + ggtitle("Luminal Epithelial cells")
p1 = LabelPoints(plot = p1, points = luminal.markers[1:10,1], repel = TRUE, ynudge = 0)

p2 = ggplot(avg.bsl, aes(s4_treated_het,midpoint_treated)) + geom_point() + ggtitle("Basal Epithelial cells")
p2 = LabelPoints(plot = p2, points = basal.markers[1:10,1], repel = TRUE)

p3 = ggplot(avg.lp, aes(s4_treated_het,s2_treated_het)) + geom_point() + ggtitle("Luminal Progenitor cells")
p3 = LabelPoints(plot = p1, points = lp.names, repel = TRUE, ynudge = 0)


p1 + p2


#####heatmap attempt######
#bsl
avg.bsl %>% as_tibble()
pheatmap(avg.bsl[1:2,], 
    color = heat_colors, 
    cluster_rows = T, 
    show_rownames = F,
    annotation = meta.filt, 
    border_color = NA, 
    fontsize = 10, 
    scale = "column", 
    fontsize_row = 10, 
    height = 20)

VlnPlot(treat.combined, features = c("B2m","Cd74","Iigp1","Isg15","Pfn1","Ifitm3","Eef2"), split.by = "orig.ident" )



plots.b = VlnPlot(bsl.dge, features = c("Cd74","Ly6a","Ifitm3","Iigp1"), split.by = "orig.ident",
    pt.size = .3, combine = FALSE)
wrap_plots(plots = plots.b, ncol = 4)

VlnPlot(fvb50, features = c(c("Cd74","Ly6a","Iigp1","Sfn","Cldn3","Krt18")))
#ly6a is aka sca-1



plots.l = VlnPlot(luminal.dge, features = c("Cd74","Stmn1","","Ifitm3","Pfn1"), split.by = "orig.ident",
    pt.size = .3, combine = FALSE)

p1= FeaturePlot(treat.combined, features = "Cd74", split.by = "orig.ident", cols = c("grey","firebrick1") )

basal.luminal.sub = subset(treat.combined, idents = c("luminal", "basal"))
p1 = FeaturePlot(basal.luminal.sub, features = "Cd74", cols = c("grey","firebrick1"), split.by = "orig.ident", pt.size = 1.5)
p2 = FeaturePlot(basal.luminal.sub, features = "Ly6a", cols = c("grey","firebrick1"), split.by = "orig.ident", pt.size = 1)
p3 = FeaturePlot(basal.luminal.sub, features = "Iigp1", cols = c("grey","firebrick1"), split.by = "orig.ident", pt.size = 1)
p4 = FeaturePlot(basal.luminal.sub, features = "Ifitm3", cols = c("grey","firebrick1"), split.by = "orig.ident", pt.size = 1)
p5 = FeaturePlot(basal.luminal.sub, features = "Pfn1", cols = c("grey","firebrick1"), split.by = "orig.ident", pt.size = 1)
p6 = FeaturePlot(basal.luminal.sub, features = "Ltf", cols = c("grey","firebrick1"), split.by = "orig.ident", pt.size = 1)
FeaturePlot(basal.luminal.sub, features = "Pfn1", cols = c("grey","firebrick1"))

(p1+p2)/(p3+p4)

p1 = FeaturePlot(s2, features = "Cd74", cols = c("azure3","red"))
p2 = FeaturePlot(s2, features = "Ly6a", cols = c("azure3","red"))
p3 = FeaturePlot(s2, features = "Iigp1", cols = c("azure3","red"))
p4 = FeaturePlot(s2, features = "Ifitm3", cols = c("azure3","red"))
FeaturePlot(s2, features = "Pfn1", cols = c("azure3","red"))
p1

(p1/p2)



e.dge = subset(treat.combined, idents = "E")
DimPlot(e.dge, reduction = 'umap')
p1 = VlnPlot(e.dge, features = "Cd74", split.by = "orig.ident")
p2 = VlnPlot(e.dge, features = "Ly6a", split.by = "orig.ident")
p3 = VlnPlot(e.dge, features = "Iigp1", split.by = "orig.ident")
p4 = FeaturePlot(basal.luminal.sub, features = "Ifitm3", cols = c("grey","firebrick1"), split.by = "orig.ident", pt.size = 1)
p5= FeaturePlot(basal.luminal.sub, features = "Pfn1", cols = c("grey","firebrick1"), split.by = "orig.ident", pt.size = 1)

########subsset immune population and DGE#######
Im = c("Ptprc")
My = c("Cd74","Lyz2")
DC = c("Traf1","Cd209a","Napsa","Flt3")
Mac = c("Csf1r","Fcgr3","Adgre1","Ms4a7")
Ma  = c("Mrc1","Cd209f","Cd163")
Mb  = c("Mmp12","Mmp13","Spic")
Ly  = c("Cd3g","Cd3d","Cd3e")
NK  = c("Gzma","Ncr1","Itgae")
Tcd8= c("Cd8a","Cd8b1")
Tcd4= c("Cd4")
B   = c("Cd79a","Cd79b")

FeaturePlot(treat.combined, features = Im, label = TRUE)
FeaturePlot(treat.combined, features = My, label = TRUE)
FeaturePlot(treat.combined, features = DC, label = TRUE)
FeaturePlot(treat.combined, features = Mac, label = TRUE)
FeaturePlot(treat.combined, features = Ly, label = TRUE)
FeaturePlot(treat.combined, features = NK, label = TRUE)
FeaturePlot(treat.combined, features = Tcd8, label = TRUE)
FeaturePlot(treat.combined, features = Tcd4, label = TRUE)
FeaturePlot(treat.combined, features = B, label = TRUE)

#note: main difference between PBS and 4NQO treated mammary tissue is the presense of lymphoid derived cells (CD3+ cells)
DefaultAssay(treat.combined) = 'RNA'
im.markers = FindConservedMarkers(treat.combined, ident.1 ="immune", grouping.var= 'orig.ident' )
  head(im.markers)
    im.conserved = rownames(im.markers[1:20,])
FeaturePlot(treat.combined, features = im.conserved[1:5], ncol = 5)


im.dge = subset(treat.combined, idents = "immune")
  Idents(im.dge)= "orig.ident"
  DimPlot(im.dge, reduction = 'umap')
avg.im = as.data.frame(log1p(AverageExpression(im.dge)$RNA))
  str(avg.im)
avg.im$gene = rownames(im.dge$RNA)

im.genes.to.label = im.markers[1:10]

o = order(avg.im$midpoint_treated)
avg.im.o = avg.im$gene[o]
  list(avg.im.o[1:10])
  im.names = avg.im.o[1:10]
  
tail(avg.im.o)
print(rev(avg.im.o))


p1 = ggplot(avg.im, aes(s4_treated_het,midpoint_treated)) + geom_point() + ggtitle("Immune cells")
p1 = LabelPoints(plot = p1, points = luminal.markers[1:10,1], repel = TRUE, ynudge = 0)

