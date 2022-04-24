library(Seurat)
library(BUSpaRse)
library(tidyverse)
library(textTinyR)
library(qdapTools)
#install.packages("qdapTools")

#install.packages("anndata")
#anndata::install_anndata()
library(anndata)
library(tidyverse)

setwd("C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling")


#read the data in backed mode, it will not fit in memory otherwise
mcor3Counts = read_h5ad("D:/SingleCellModeling/MCOR3/CountMatrices/10xv3_gene_tsne_new.h5ad", backed = 'r')
genesTmp = mcor3Counts$var_names
spl = str_split(genesTmp, "_", n = 2)
genes = genesTmp #allocate
for (i in 1:length(spl)) {
  genes[i] = spl[[i]][1]
}

metGenes = read_tsv("data/metabolicGenesMouse.txt", col_names = FALSE)[[1]]

selMetGenes = genes %in% metGenes
sum(selMetGenes) #2784
X = mcor3Counts$layers$get('X')


allData = Matrix::t(X)
allMeta = mcor3Counts$obs

#figure out which cell types to include
cellTypes = unique(allMeta$cluster_label)
cellsPerCluster = rep(NA, length(cellTypes))
for (i in 1:length(cellTypes)) {
  cellsPerCluster[i] = sum(allMeta$cluster_label == cellTypes[i] )
}
sel = cellsPerCluster >= 450
sum(sel)

cellTypesFilt = cellTypes[sel]
allLargeData = allData[,allMeta$cluster_label %in% cellTypesFilt]
rownames(allLargeData) = genes
allLargeMeta = allMeta[allMeta$cluster_label %in% cellTypesFilt,]





#####################################################################
#The original paper reports large batch effects between batches run at different time points. 
#Therefore, for the UMAP plots, let's focus on one
#So, this is a different way to cut the matrix, use all genes here
#####################################################################
unique(mcor3Counts$obs$batch)
dates = unique(mcor3Counts$obs$date)
sum(mcor3Counts$obs$date == dates[1])#19212
sum(mcor3Counts$obs$date == dates[2])#25308
sum(mcor3Counts$obs$date == dates[3])#26845

date3Data = X[mcor3Counts$obs$date == dates[3],]
date3Data = Matrix::t(date3Data)
date3Meta = mcor3Counts$obs[mcor3Counts$obs$date == dates[3],]
dim(date3Meta)

rownames(date3Data) = genes
#saveRDS(date3Data, file="D:/SingleCellModeling/MCOR3/CountMatrices/Date3CountMatrix.rds")
#saveRDS(date3Meta, file="D:/SingleCellModeling/MCOR3/CountMatrices/Date3MetaData.rds")

#date3Data = readRDS(file="D:/SingleCellModeling/MCOR3/CountMatrices/Date3CountMatrix.rds")
#date3Meta = readRDS(file="D:/SingleCellModeling/MCOR3/CountMatrices/Date3MetaData.rds")

#aggregate the rows with the same gene name:
library(Matrix.utils)
date3Data2 = aggregate.Matrix(date3Data, rownames(date3Data), fun='sum')

############################################
#Generate UMAP figures from the date3 data:
############################################
memory.limit(size = 100000)

#Only include the selected clusters:
sel2 = date3Meta$cluster_label %in% cellTypes[sel]
date3Data2Filt = date3Data2[,sel2]

library(Seurat)
#test to run through the date 3 data and see if it clusters according to the clusters specified
gc()
set.seed(1)
seurObj = CreateSeuratObject(counts = date3Data2Filt, project = "MCOR3", min.cells = 0, min.features = 0)
seurObj[["authors_cluster"]] = date3Meta$cluster_label[sel2]
seurObj[["authors_class"]] = date3Meta$class_label[sel2]
seurObj[["authors_subclass"]] = date3Meta$subclass_label[sel2]
seurObj <- NormalizeData(seurObj, normalization.method = "LogNormalize", scale.factor = 10000)
seurObj <- FindVariableFeatures(seurObj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurObj)
seurObj <- ScaleData(seurObj, features = all.genes)
seurObj <- RunPCA(seurObj, features = VariableFeatures(object = seurObj))
seurObj <- FindNeighbors(seurObj, dims = 1:15)
seurObj <- FindClusters(seurObj, resolution = 0.5)
seurObj <- RunUMAP(seurObj, dims = 1:15)

pA1 = DimPlot(seurObj, group.by = "authors_cluster")
pA1 = pA1 + labs(x="UMAP x", y="UMAP y", title="All genes")
pA1 = pA1 + theme(legend.position = "none", legend.text=element_text(size=14), text = element_text(size=14),plot.title=element_text(size=16.8, family="Arial", face="plain"))
pA1


figurePath = "D:/SingleCellModeling/figures/"


ggsave(
  paste0(figurePath, "MCOR3UMAPAllGenes.png"),
  plot = pA1, device = "png",
  width = 4.5, height = 4.5, dpi = 300)


rm(seurObj)
gc()


#pClass = DimPlot(seurObj, group.by = "authors_class")
#pSubClAll = DimPlot(seurObj, group.by = "authors_subclass") + ggtitle("All genes, authors' subclass")
#pFeat = FeaturePlot(seurObj, c("Naa30","Man1a"))

#pClass

##########################
#Now with metabolic genes
##########################

metGenes = read_tsv("data/metabolicGenesMouse.txt", col_names = FALSE)[[1]]
selMetGenes = rownames(date3Data2Filt) %in% metGenes
sum(selMetGenes) #2781

date3Data2FiltMet = date3Data2Filt[selMetGenes,]
set.seed(1)
seurObjMet = CreateSeuratObject(counts = date3Data2FiltMet, project = "MCOR3", min.cells = 0, min.features = 0)
seurObjMet[["authors_cluster"]] = date3Meta$cluster_label[sel2]
seurObjMet[["authors_class"]] = date3Meta$class_label[sel2]
seurObjMet[["authors_subclass"]] = date3Meta$subclass_label[sel2]
seurObjMet <- NormalizeData(seurObjMet, normalization.method = "LogNormalize", scale.factor = 10000)
seurObjMet <- FindVariableFeatures(seurObjMet, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurObjMet)
seurObjMet <- ScaleData(seurObjMet, features = all.genes)
seurObjMet <- RunPCA(seurObjMet, features = VariableFeatures(object = seurObjMet))
seurObjMet <- FindNeighbors(seurObjMet, dims = 1:15)
seurObjMet <- FindClusters(seurObjMet, resolution = 0.5)
seurObjMet <- RunUMAP(seurObjMet, dims = 1:15)
#pA = DimPlot(seurObjMet)
#pA

pA2 = DimPlot(seurObjMet, group.by = "authors_cluster")
pA2 = pA2 + labs(x="UMAP x", y="UMAP y", title="Metabolic genes")
pA2 = pA2 + theme(legend.position = "right", legend.text=element_text(size=14), text = element_text(size=14),plot.title=element_text(size=16.8, family="Arial", face="plain"))
pA2


figurePath = "D:/SingleCellModeling/figures/"


ggsave(
  paste0(figurePath, "MCOR3UMAPMetGenes.png"),
  plot = pA2, device = "png",
  width = 6.75, height = 4.5, dpi = 300)

rm(seurObj)
rm(seurObjMet)

##################################
# Generate bootstrap data for the 16 populations - use all cells
##################################


unique(allMeta$class_label[allMeta$cluster_label %in% cellTypesFilt]) #all low quality and non-neuronal gone, good

#Now create bootstraps, 100 per sample

allLargeData = as(allLargeData, "dgCMatrix")#much faster for some reason


for (i in 1:length(cellTypesFilt)) {
  print(paste0(i," of ", length(cellTypesFilt)))
  subMat = allLargeData[, allLargeMeta$cluster_label == cellTypesFilt[i]]
  dim(subMat)
  pooledBootstraps = Matrix(0, nrow = nrow(subMat), ncol=100)
  rownames(pooledBootstraps) = rownames(subMat)
  for (j in 1:100) {
    # bootstrap,i.e. sample the same number of samples with replacement 
    sel3 = sample(ncol(subMat),ncol(subMat), replace = TRUE) 
    #length(sel1[sel2])#Test: 1809, ok
    pooledBootstraps[,j] = rowSums(subMat[,sel3])
  }
  tibbBstr = as_tibble(pooledBootstraps) %>% add_column(gene = rownames(pooledBootstraps), .before = 1)
  sampleName = ctFilt[i];
  sampleName = gsub("/", "_", sampleName)
  fn = paste0("data/MCOR3/Run5/Bootstrap_", sampleName,".txt");
  write_tsv(tibbBstr, file=fn)
  gc()
}

#############################################################
#also save the profiles as they are, without any bootstrapping
#############################################################
ctExpr = Matrix(0, nrow = nrow(UMICountsFilt), ncol=length(ctFilt))
colnames(ctExpr) = ctFilt
rownames(ctExpr) = rownames(UMICountsFilt)
for (i in 1:length(ctFilt)) {
  print(paste0(i," of ", length(ctFilt)))
  subMat = UMICountsFilt[, metaFilt$cluster_label == ctFilt[i]]
  dim(subMat)
  ctExpr[,i] = rowSums(subMat)
}

tibbBstr = as_tibble(ctExpr) %>% add_column(gene = rownames(ctExpr), .before = 1)
fn = "data/MCOR3/pooledCT.txt";
write_tsv(tibbBstr, file=fn)



#######################################
#Investigate the UMIs per cell distribution per cell subtype
#######################################

allLargeData = as(allLargeData, "dgCMatrix")#much faster for some reason

UMIsPerCellPerCluster = vector(mode = "list", length = length(cellTypesFilt))

for (i in 1:length(cellTypesFilt)) {
  subMat = allLargeData[, allLargeMeta$cluster_label == cellTypesFilt[i]]
  UMIsPerCellPerCluster[[i]] = as.numeric(colSums(subMat))
}


x = NULL
group = NULL

for (i in 1:length(UMIsPerCellPerCluster)) {
  x = c(x,UMIsPerCellPerCluster[[i]])
  group = c(group, rep(i,length(UMIsPerCellPerCluster[[i]])))
}
df = tibble(x = x, Cluster = factor(group, 1:length(UMIsPerCellPerCluster), cellTypesFilt))

df2 = df[df$x < 1000000 ,]

color_palette = c('#000000','#888888','#FF0000','#00BB00','#DCDC55','#6B97EC','#BC976B', '#0000FF', '#444444', '#CC00CC') 

labels = cellTypesFilt

pSup3 = ggplot(df2, aes(x = x, color = Cluster, linetype = Cluster, size=Cluster)) +
  #geom_density(aes(colour=Cluster),show.legend=FALSE, size=1.5) +
  stat_density(aes(x=x, colour=Cluster,), geom="line", position="identity") +
  scale_color_manual(values = color_palette[c(1:6,1:6,1:5)], labels = labels) +
  scale_linetype_manual(values = c(rep("solid",6),rep("dashed",6), rep("dotted",5)), labels=labels) +
  scale_size_manual(values = rep(1.3,17), labels=labels) +
  labs(x="UMIs per cell", y="Density", title="") + theme_bw() + 
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

figurePath = "D:/SingleCellModeling/figures/"

ggsave(
  paste0(figurePath, "MCOR3_UMIsPerCellDens.png"),
  plot = pSup3, device = "png",
  width = 12, height = 8, dpi = 300)



#########################
#Run DSAVE on the data and create supplementary fig
#########################

set.seed(1)

#sum(rownames(allLargeMeta) != colnames(allLargeData))#0, ok, so they are in order

library(DSAVE)

d = allLargeData[, allLargeMeta$cluster_label == cellTypesFilt[14]]
dim(d)#[1] 24575   498
varVipChat = DSAVEGetTotalVariationPoolSize(d,upperBound = 50, lowerBound = 5e-1, poolSizes = c(100,150,200,245))

d = allLargeData[, allLargeMeta$cluster_label == cellTypesFilt[11]]
dim(d)#[1] 24575   1295
varL6NPTrh1 = DSAVEGetTotalVariationPoolSize(d,upperBound = 50, lowerBound = 5e-1, poolSizes = c(100,150,200,400,640))

d = allLargeData[, allLargeMeta$cluster_label == cellTypesFilt[9]]
dim(d)#[1] 24575   4272
varL5ITS100b = DSAVEGetTotalVariationPoolSize(d,upperBound = 50, lowerBound = 5e-1, poolSizes = c(100,150,200,400,800,1600,2100))


dsBulk = DSAVE::bulkTotalVar1vs1[[4]] [[2]]

#transform the data for ggplot
x = c(0,2150, varVipChat$poolSizes, varL6NPTrh1$poolSizes, varL5ITS100b$poolSizes)
y = c(dsBulk, dsBulk, varVipChat$Rs, varL6NPTrh1$Rs, varL5ITS100b$Rs)
labels = c("Bulk T", "Vip Chat", "L6 NP Trh_1", "L5 IT S100b")
group = factor(c(1,1,rep(2, length(varVipChat$Rs)),rep(3, length(varL6NPTrh1$Rs)), rep(4, length(varL5ITS100b$Rs))),1:4, labels)


df = tibble(x=x, y=y, Dataset=group)

#col = c(rgb(0,0,0),gg_color_hue(6),rgb(0.5,0.5,0.5))
color_palette = c('#000000','#888888','#FF0000','#00BB00','#DCDC55','#6B97EC','#BC976B', '#0000FF', '#444444', '#CC00CC') 

cutDf = df#[df$y >= 0.6,]


pSup1 = ggplot(cutDf, aes(x = x, y = y, color=Dataset, linetype=Dataset, size=Dataset)) +
  geom_line() +
  scale_linetype_manual(values = c(2,1,1,1), labels = labels) +
  scale_color_manual(values = color_palette[c(1,2,5,3)], labels = labels) +
  scale_size_manual(values = rep(1.3,4), labels = labels) +
  #ylim(1,1.0015) +
  #ggplot2::labs(y=expression("DSAVE score"), x=expression(Log[2]*"(number of cells)"), title="DSAVE") +
  ggplot2::labs(y=expression("DSAVE score"), x="Number of cells", title="") +
  ggplot2::theme_bw() +#+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  ggplot2::theme(legend.title = element_blank(),legend.position="right", legend.text=element_text(color='black',size=14)) + guides(colour = guide_legend(nrow = 4), size = guide_legend(nrow = 4), linetype = guide_legend(nrow = 4)) +
  ggplot2::theme(text = element_text(size=14), axis.text.x = element_text(color='black', size=14), axis.text.y = element_text(color='black', size=14))
pSup1

#conclusion: aim for cluster with at least 1500-2000 cells

figurePath = "D:/SingleCellModeling/figures/"

ggsave(
  paste0(figurePath, "Fig 3 Sup DSAVE.png"),
  plot = pSup1, device = "png",
  width = 8, height = 5, dpi = 300)












