library(R.matlab)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(Matrix)
library(DSAVE)

fig____path = "Z:/projects/Single-cell modeling/figures/"
fig____path_loc = "D:/SingleCellModeling/figures/"

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

setwd("C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling")


###############################
# Load the dataset
###############################


UMICounts = readRDS("D:/SingleCellModeling/LC3/rawUMISparse.rds")


cs = colSums(UMICounts)
plot(density(cs))

cellAnnotTmp = read_tsv("D:/SingleCellModeling/LC3/GSE131907_Lung_Cancer_cell_annotation.txt")#other computer
unique(cellAnnotTmp$Cell_type )#10
unique(cellAnnotTmp$Cell_subtype )#50
unique(cellAnnotTmp$Sample )

#synchronize the rows in cellAnnot with the cols in UMICounts
library(qdapTools)

lookupTable = tibble(x=cellAnnotTmp$Index,y=1:ncol(UMICounts))
newIndices = lookup(colnames(UMICounts), lookupTable)
sum(is.na(newIndices))#ok, no NA
cellAnnot = cellAnnotTmp[newIndices,]

sum(is.na(cellAnnot$Cell_subtype))

dim(cellAnnot)#208506      7
cellAnnot[is.na(cellAnnot$Cell_subtype),]#30,065 NAs, so not a small number - we will not use these cells I guess

#there are many clusters with enough cells, but we need to split it up in at least normal and cancer
#split in a few such categories
sampleIds = unique(cellAnnot$Sample)

#lets focus on the early stage patients, tumor and healthy lung tissue
#these are the ones labeled "LUNG_T*" and "LUNG_N*"
tumorSamp = sampleIds[startsWith(sampleIds, "LUNG_T")]
normalSamp = sampleIds[startsWith(sampleIds, "LUNG_N")]

#########################
#Run DSAVE on the data and create supplementary fig
#########################

set.seed(1)
ids = cellAnnot$Index[!is.na(cellAnnot$Cell_subtype) & (cellAnnot$Cell_subtype == "NK") & cellAnnot$Sample %in% normalSamp ]
length(ids)#4752
d = UMICounts[,colnames(UMICounts) %in% ids]
dim(d)#29634 4752
varNK = DSAVEGetTotalVariationPoolSize(d,upperBound = 50, lowerBound = 5e-1)

ids = cellAnnot$Index[!is.na(cellAnnot$Cell_subtype) & (cellAnnot$Cell_subtype == "AT2") & cellAnnot$Sample %in% normalSamp]
d = UMICounts[,colnames(UMICounts) %in% ids]
dim(d)#29634 2020
varAT2 = DSAVEGetTotalVariationPoolSize(d,upperBound = 50, lowerBound = 5e-1)
#DSAVEPlotTotalVariation(list(varNK), list("NK cells"), 4)

dsBulk = DSAVE::bulkTotalVar1vs1[[4]] [[2]]

#transform the data for ggplot
x = c(0,2200, varNK$poolSizes, varAT2$poolSizes )
y = c(dsBulk, dsBulk, varNK$Rs, varAT2$Rs)
labels = c("Bulk T", "N: NK", "N: AT2")
group = factor(c(1,1,rep(2, length(varNK$Rs)),rep(3, length(varAT2$Rs))),1:3, labels)


df = tibble(x=x, y=y, Dataset=group)

#col = c(rgb(0,0,0),gg_color_hue(6),rgb(0.5,0.5,0.5))
color_palette = c('#000000','#888888','#FF0000','#00BB00','#DCDC55','#6B97EC','#BC976B', '#0000FF', '#444444', '#CC00CC') 

cutDf = df#[df$y >= 0.6,]


pSup1 = ggplot(cutDf, aes(x = x, y = y, color=Dataset, linetype=Dataset, size=Dataset)) +
  geom_line() +
  scale_linetype_manual(values = c(2,1,1), labels = labels) +
  scale_color_manual(values = color_palette[c(1,4,5)], labels = labels) +
  scale_size_manual(values = rep(1.3,3), labels = labels) +
  #ylim(1,1.0015) +
  #ggplot2::labs(y=expression("DSAVE score"), x=expression(Log[2]*"(number of cells)"), title="DSAVE") +
  ggplot2::labs(y=expression("DSAVE score"), x="Number of cells", title="DSAVE evaluation") +
  ggplot2::theme_bw() +#+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  ggplot2::theme(legend.title = element_blank(),legend.position="right", legend.text=element_text(color='black',size=14)) + guides(colour = guide_legend(nrow = 4), size = guide_legend(nrow = 4), linetype = guide_legend(nrow = 4)) +
  ggplot2::theme(text = element_text(size=14), axis.text.x = element_text(color='black', size=14), axis.text.y = element_text(color='black', size=14))
pSup1

#conclusion: aim for cluster with at least 1500-2000 cells

figurePath = "D:/SingleCellModeling/figures/"

ggsave(
  paste0(figurePath, "Fig 4 Sup 1.png"),
  plot = pSup1, device = "png",
  width = 5, height = 5, dpi = 300)



###############################
# Analyze with Seurat
###############################


dim(UMICounts)

UMICountsTumor = UMICounts[, (cellAnnot$Sample %in% tumorSamp) & (!is.na(cellAnnot$Cell_subtype)) & (cellAnnot$Cell_type != "Undetermined")]
dim(UMICountsTumor)#29634 38489
annotTumor2 = cellAnnot[(cellAnnot$Sample %in% tumorSamp) & (!is.na(cellAnnot$Cell_subtype)) & (cellAnnot$Cell_type != "Undetermined"),]
dim(annotTumor2)#38489     7

#cellAnnot$Index[1:20]
#colnames(UMICounts)[1:20]#looks good


library(Seurat)
#test to run through the date 3 data and see if it clusters according to the clusters specified
#memory.limit(size = 40000)#otherwise we run out of memory

#now filter rare types to reduce the legend
subtypesTumor = unique(annotTumor2$Cell_subtype)
cellsPerSubtype = rep(NA, length(subtypesTumor))
for (i in 1:length(subtypesTumor)) {
  cellsPerSubtype[i] = sum(annotTumor2$Cell_subtype == subtypesTumor[i])
}

tibble(subtypesTumor[cellsPerSubtype > 1000], cellsPerSubtype[cellsPerSubtype > 1000]) #we use 1600 as a limit
cellSubtypesTumor3 = subtypesTumor[cellsPerSubtype > 1600]

UMICountsTumor3 = UMICountsTumor[,annotTumor2$Cell_subtype %in% cellSubtypesTumor3] 
annotTumor3 = annotTumor2[annotTumor2$Cell_subtype %in% cellSubtypesTumor3,]

set.seed(1)
seurObj = CreateSeuratObject(counts = UMICountsTumor3, project = "MCOR3", min.cells = 0, min.features = 0)
seurObj[["authors_cell_type"]] = annotTumor3$Cell_type
seurObj[["authors_cell_subtype"]] = annotTumor3$Cell_subtype
seurObj[["authors_sample"]] = annotTumor3$Sample
seurObj <- NormalizeData(seurObj, normalization.method = "LogNormalize", scale.factor = 10000)
seurObj <- FindVariableFeatures(seurObj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurObj)
seurObj <- ScaleData(seurObj, features = all.genes)
seurObj <- RunPCA(seurObj, features = VariableFeatures(object = seurObj))
seurObj <- FindNeighbors(seurObj, dims = 1:15)
seurObj <- FindClusters(seurObj, resolution = 0.5)
seurObj <- RunUMAP(seurObj, dims = 1:15)
#pA = DimPlot(seurObj)
#pA

#DimPlot(seurObj, group.by = "authors_cell_type")
DimPlot(seurObj, group.by = "authors_cell_subtype")

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=16)
cl1 = color_list[c(1,2,4,5,7,9,10,12,13,15)]
cl2 = color_list[c(3,6,8,11,14,16)]

p4A = DimPlot(seurObj, group.by = "authors_cell_subtype", cols=cl1)
p4A = p4A + labs(x="UMAP x", y="UMAP y", title="Tumor")
p4A = p4A + theme(legend.position = "bottom", legend.text=element_text(size=14), text = element_text(size=14),plot.title=element_text(size=16.8, family="Arial", face="plain"))
p4A = p4A + guides(colour = guide_legend(nrow=5,override.aes = list(size=3)))
p4A
figurePath = "D:/SingleCellModeling/figures/"

ggsave(
  paste0(figurePath, "Fig 4A.png"),
  plot = p4A, device = "png",
  width = 4.5, height = 5.92, dpi = 300)

p4C = DimPlot(seurObj, group.by = "authors_sample")
p4C = p4C + labs(x="UMAP x", y="UMAP y", title="Tumor, per sample")
p4C = p4C + theme(legend.position = "bottom", legend.text=element_text(size=14), text = element_text(size=14),plot.title=element_text(size=16.8, family="Arial", face="plain"))
p4C = p4C + guides(colour = guide_legend(nrow=6,override.aes = list(size=3)))
p4C
figurePath = "D:/SingleCellModeling/figures/"

ggsave( #we get some warnings about the font, but it seems to work
  paste0(figurePath, "Fig 4C_2.png"),
  plot = p4C, device = "png",
  width = 4.5, height = 6.14, dpi = 300)



#Now using metabolic genes:
##############################

metGenes = read_tsv("data/metabolicGenesHuman.txt", col_names = FALSE)[[1]]
length(metGenes)#3060, ok, same as in matlab
UMICountsTumor3Met = UMICountsTumor3[rownames(UMICountsTumor3) %in% metGenes,]
dim(UMICountsTumor3Met)#2912 27504, looks ok

set.seed(1)
seurObj = CreateSeuratObject(counts = UMICountsTumor3Met, project = "MCOR3", min.cells = 0, min.features = 0)
seurObj[["authors_cell_type"]] = annotTumor3$Cell_type
seurObj[["authors_cell_subtype"]] = annotTumor3$Cell_subtype
seurObj <- NormalizeData(seurObj, normalization.method = "LogNormalize", scale.factor = 10000)
seurObj <- FindVariableFeatures(seurObj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurObj)
seurObj <- ScaleData(seurObj, features = all.genes)
seurObj <- RunPCA(seurObj, features = VariableFeatures(object = seurObj))
seurObj <- FindNeighbors(seurObj, dims = 1:15)
seurObj <- FindClusters(seurObj, resolution = 0.5)
seurObj <- RunUMAP(seurObj, dims = 1:15)
#pA = DimPlot(seurObj)
#pA

#DimPlot(seurObj, group.by = "authors_cell_type")
#DimPlot(seurObj, group.by = "authors_cell_subtype")

pSup2A = DimPlot(seurObj, group.by = "authors_cell_subtype", cols=cl1)
pSup2A = pSup2A + labs(x="UMAP x", y="UMAP y", title="Tumor, metabolic genes")
pSup2A = pSup2A + theme(legend.position = "bottom", legend.text=element_text(size=14), text = element_text(size=14),plot.title=element_text(size=16.8, family="Arial", face="plain"))
pSup2A = pSup2A + guides(colour = guide_legend(nrow=5,override.aes = list(size=3)))
pSup2A

ggsave(
  paste0(figurePath, "Fig 4 Sup 2A.png"),
  plot = pSup2A, device = "png",
  width = 4.5, height = 5.92, dpi = 300)

#######################################
# Also Seurat for the normal tissue samples
#######################################

UMICountsNormal = UMICounts[, (cellAnnot$Sample %in% normalSamp) & (!is.na(cellAnnot$Cell_subtype)) & (cellAnnot$Cell_type != "Undetermined")]
dim(UMICountsNormal)#29634 36381
annotNormal2 = cellAnnot[(cellAnnot$Sample %in% normalSamp) & (!is.na(cellAnnot$Cell_subtype)) & (cellAnnot$Cell_type != "Undetermined"),]
dim(annotNormal2)#36381     7, ok

#set.seed(1)
#seurObj = CreateSeuratObject(counts = UMICountsNormal, project = "MCOR3", min.cells = 0, min.features = 0)
#seurObj[["authors_cell_type"]] = annotNormal2$Cell_type
#seurObj[["authors_cell_subtype"]] = annotNormal2$Cell_subtype
#seurObj <- NormalizeData(seurObj, normalization.method = "LogNormalize", scale.factor = 10000)
#seurObj <- FindVariableFeatures(seurObj, selection.method = "vst", nfeatures = 2000)
#all.genes <- rownames(seurObj)
#seurObj <- ScaleData(seurObj, features = all.genes)
#seurObj <- RunPCA(seurObj, features = VariableFeatures(object = seurObj))
#seurObj <- FindNeighbors(seurObj, dims = 1:15)
#seurObj <- FindClusters(seurObj, resolution = 0.5)
#seurObj <- RunUMAP(seurObj, dims = 1:15)
#pA = DimPlot(seurObj)
#pA

#DimPlot(seurObj, group.by = "authors_cell_type")

#now filter rare types to reduce the legend
subtypesNormal = unique(annotNormal2$Cell_subtype)
cellsPerSubtype = rep(NA, length(subtypesNormal))
for (i in 1:length(subtypesNormal)) {
  cellsPerSubtype[i] = sum(annotNormal2$Cell_subtype == subtypesNormal[i])
}

tibble(subtypesNormal[cellsPerSubtype > 1000], cellsPerSubtype[cellsPerSubtype > 1000]) #we use 1600 as a limit
cellSubtypesNormal3 = subtypesNormal[cellsPerSubtype > 1600]#all clusters included > 2000 cells

UMICountsNormal3 = UMICountsNormal[,annotNormal2$Cell_subtype %in% cellSubtypesNormal3] 
annotNormal3 = annotNormal2[annotNormal2$Cell_subtype %in% cellSubtypesNormal3,]

set.seed(1)
seurObj = CreateSeuratObject(counts = UMICountsNormal3, project = "MCOR3", min.cells = 0, min.features = 0)
seurObj[["authors_cell_type"]] = annotNormal3$Cell_type
seurObj[["authors_cell_subtype"]] = annotNormal3$Cell_subtype
seurObj <- NormalizeData(seurObj, normalization.method = "LogNormalize", scale.factor = 10000)
seurObj <- FindVariableFeatures(seurObj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurObj)
seurObj <- ScaleData(seurObj, features = all.genes)
seurObj <- RunPCA(seurObj, features = VariableFeatures(object = seurObj))
seurObj <- FindNeighbors(seurObj, dims = 1:15)
seurObj <- FindClusters(seurObj, resolution = 0.5)
seurObj <- RunUMAP(seurObj, dims = 1:15)
#pA = DimPlot(seurObj)
#pA

#DimPlot(seurObj, group.by = "authors_cell_type")
#DimPlot(seurObj, group.by = "authors_cell_subtype")

p4B = DimPlot(seurObj, group.by = "authors_cell_subtype", cols=cl2)
p4B = p4B + labs(x="UMAP x", y="UMAP y", title="Normal tissue")
p4B = p4B + theme(legend.position = "bottom", legend.text=element_text(size=14), text = element_text(size=14),plot.title=element_text(size=16.8, family="Arial", face="plain"))
p4B = p4B + guides(colour = guide_legend(nrow=3,override.aes = list(size=3)))
p4B

ggsave(
  paste0(figurePath, "Fig 4B.png"),
  plot = p4B, device = "png",
  width = 4.5, height = 5.5, dpi = 300)


#Now using metabolic genes:
##############################

metGenes = read_tsv("data/metabolicGenesHuman.txt", col_names = FALSE)[[1]]
length(metGenes)#3060, ok, same as in matlab
UMICountsNormal3Met = UMICountsNormal3[rownames(UMICountsNormal3) %in% metGenes,]
dim(UMICountsNormal3Met)#2912 27504, looks ok

set.seed(1)
seurObj = CreateSeuratObject(counts = UMICountsNormal3Met, project = "MCOR3", min.cells = 0, min.features = 0)
seurObj[["authors_cell_type"]] = annotNormal3$Cell_type
seurObj[["authors_cell_subtype"]] = annotNormal3$Cell_subtype
seurObj <- NormalizeData(seurObj, normalization.method = "LogNormalize", scale.factor = 10000)
seurObj <- FindVariableFeatures(seurObj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurObj)
seurObj <- ScaleData(seurObj, features = all.genes)
seurObj <- RunPCA(seurObj, features = VariableFeatures(object = seurObj))
seurObj <- FindNeighbors(seurObj, dims = 1:15)
seurObj <- FindClusters(seurObj, resolution = 0.5)
seurObj <- RunUMAP(seurObj, dims = 1:15)
#pA = DimPlot(seurObj)
#pA

#DimPlot(seurObj, group.by = "authors_cell_type")
#DimPlot(seurObj, group.by = "authors_cell_subtype")

pSup2B = DimPlot(seurObj, group.by = "authors_cell_subtype", cols=cl2)
pSup2B = pSup2B + labs(x="UMAP x", y="UMAP y", title="Normal tissue, metabolic genes")
pSup2B = pSup2B + theme(legend.position = "bottom", legend.text=element_text(size=14), text = element_text(size=14),plot.title=element_text(size=16.8, family="Arial", face="plain"))
pSup2B = pSup2B + guides(colour = guide_legend(nrow=3,override.aes = list(size=3)))
pSup2B

ggsave(
  paste0(figurePath, "Fig 4 Sup 2B.png"),
  plot = pSup2B, device = "png",
  width = 4.5, height = 5.5, dpi = 300)


###########################
# Create profiles for fig 2, without bootstraps
###########################


UMICountsTumor3
annotTumor3
cellSubtypesTumor3
  
UMICountsNormal3
annotNormal3
cellSubtypesNormal3




numSamp = length(cellSubtypesTumor3) + length(cellSubtypesNormal3)

pooledSamp = Matrix(0, nrow = nrow(UMICounts), ncol=numSamp)
rownames(pooledSamp) = rownames(UMICounts)
colnames(pooledSamp) = c(paste0("T: ",cellSubtypesTumor3),paste0("N: ",cellSubtypesNormal3))
#first the tumor samples
for (i in 1:length(cellSubtypesTumor3)) {
  subMat = UMICountsTumor3[, (annotTumor3$Cell_subtype == cellSubtypesTumor3[i])]
  pooledSamp[,i] = rowSums(subMat)
}
for (i in 1:length(cellSubtypesNormal3)) {
  subMat = UMICountsNormal3[, (annotNormal3$Cell_subtype == cellSubtypesNormal3[i])]
  pooledSamp[,length(cellSubtypesTumor3) + i] = rowSums(subMat)
}
pooledSamp
tibb = as_tibble(pooledSamp) %>% add_column(gene = rownames(pooledSamp), .before = 1)
tibb

write_tsv(tibb, file="data/LC3/ModelDataLC3.txt")

#######################################
#Investigate the UMIs per cell distribution per cell subtype
#######################################

UMIsPerCellPerCluster = vector(mode = "list", length = ncol(pooledSamp))

#first the tumor samples
for (i in 1:length(cellSubtypesTumor3)) {
  subMat = UMICountsTumor3[, (annotTumor3$Cell_subtype == cellSubtypesTumor3[i])]
  UMIsPerCellPerCluster[[i]] = as.numeric(colSums(subMat))
}

#Then the normal samples
for (i in 1:length(cellSubtypesNormal3)) {
  subMat = UMICountsNormal3[, (annotNormal3$Cell_subtype == cellSubtypesNormal3[i])]
  UMIsPerCellPerCluster[[length(cellSubtypesTumor3) + i]] = as.numeric(colSums(subMat))
}

x = NULL
group = NULL

for (i in 1:length(UMIsPerCellPerCluster)) {
  x = c(x,UMIsPerCellPerCluster[[i]])
  group = c(group, rep(i,length(UMIsPerCellPerCluster[[i]])))
}
df = tibble(x = x, Cluster = factor(group, 1:length(UMIsPerCellPerCluster), colnames(pooledSamp)))

df2 = df[df$x < 10000 ,]

color_palette = c('#000000','#888888','#FF0000','#00BB00','#DCDC55','#6B97EC','#BC976B', '#0000FF', '#444444', '#CC00CC') 

labels = colnames(pooledSamp)

pSup3 = ggplot(df2, aes(x = x, color = Cluster, linetype = Cluster, size=Cluster)) +
  #geom_density(aes(colour=Cluster),show.legend=FALSE, size=1.5) +
  stat_density(aes(x=x, colour=Cluster,), geom="line", position="identity") +
  scale_color_manual(values = color_palette[c(1:10,1:6)], labels = labels) +
  scale_linetype_manual(values = c(rep("solid",10),rep("dashed",6)), labels=labels) +
  scale_size_manual(values = rep(1.3,16), labels=labels) +
  labs(x="UMIs per cell", y="Density", title="UMIs per cell per cluster") + theme_bw() + 
  ggplot2::theme_bw() +
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  ggplot2::theme(legend.title = element_blank(),legend.position="right", legend.text=element_text(color='black',size=14)) + #guides(colour = guide_legend(nrow = 4), size = guide_legend(nrow = 4), linetype = guide_legend(nrow = 4)) +
  ggplot2::theme(text = element_text(size=14), axis.text.x = element_text(color='black', size=14), axis.text.y = element_text(color='black', size=14))

pSup3

meanUMIsPerCell = rep(NA, length(UMIsPerCellPerCluster))
for (i in 1:length(meanUMIsPerCell)) {
  meanUMIsPerCell[i] = mean(UMIsPerCellPerCluster[[i]])
}
tibble(c(cellSubtypesTumor3, cellSubtypesNormal3), meanUMIsPerCell)


ggsave(
  paste0(fig____path_loc, "LC_UMIsPerCellDens.png"),
  plot = pSup3, device = "png",
  width = 8, height = 5, dpi = 300)


##################################
# Generate bootstraps for the selected cell subtypes
##################################



#first the tumor samples
for (i in 1:length(cellSubtypesTumor3)) {
  print(paste0(i," of ", length(cellSubtypesTumor3)))
  subMat = UMICountsTumor3[, (annotTumor3$Cell_subtype == cellSubtypesTumor3[i])]
  pooledBootstraps = Matrix(0, nrow = nrow(subMat), ncol=100)
  rownames(pooledBootstraps) = rownames(subMat)
  for (j in 1:100) {
    #bootstrap,i.e. sample the same number of samples with replacement 
    sel2 = sample(ncol(subMat),ncol(subMat), replace = TRUE) 
    #length(sel1[sel2])#Test: 1809, ok
    pooledBootstraps[,j] = rowSums(subMat[,sel2])
  }
  tibbBstr = as_tibble(pooledBootstraps) %>% add_column(gene = rownames(pooledBootstraps), .before = 1)
  sampleName = cellSubtypesTumor3[i];
  sampleName = gsub("/", "_", sampleName)
  fn = paste0("data/LC3/Run5/Bootstrap_T_", sampleName,".txt");
  write_tsv(tibbBstr, file=fn)
}

#Then the normal samples
for (i in 1:length(cellSubtypesNormal3)) {
  print(paste0(i," of ", length(cellSubtypesNormal3)))
  subMat = UMICountsNormal3[, (annotNormal3$Cell_subtype == cellSubtypesNormal3[i])]
  pooledBootstraps = Matrix(0, nrow = nrow(subMat), ncol=100)
  rownames(pooledBootstraps) = rownames(subMat)
  for (j in 1:100) {
    #bootstrap,i.e. sample the same number of samples with replacement 
    sel2 = sample(ncol(subMat),ncol(subMat), replace = TRUE) 
    #length(sel1[sel2])#Test: 1809, ok
    pooledBootstraps[,j] = rowSums(subMat[,sel2])
  }
  tibbBstr = as_tibble(pooledBootstraps) %>% add_column(gene = rownames(pooledBootstraps), .before = 1)
  sampleName = cellSubtypesNormal3[i];
  sampleName = gsub("/", "_", sampleName)
  fn = paste0("data/LC3/Run5/Bootstrap_N_", sampleName,".txt");
  write_tsv(tibbBstr, file=fn)
}
