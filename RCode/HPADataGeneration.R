#Generates DSAVE data for the clusters, decides which are large enough, and 
#saves the DSAVE data, and generates bootstrap data for the selected clusters
library(tidyverse)
memory.limit(size = 100000)

fig__dir = "D:/SingleCellModeling/MetAtlas/Figures/"

library(Matrix.utils) #for transpose of sparse matrix

sparseMat = readRDS("D:/SingleCellModeling/MetAtlas/SparseMat.RDS")
metaData = readRDS("D:/SingleCellModeling/MetAtlas/metaData.RDS")

tissClust = metaData[,c(1,3)]
tissClust2 = unique(tissClust)
unique(tissClust2$Tissue)#26 tissues

nm = rep(NA,dim(tissClust2)[1]) #444 in total
for (i in 1:length(nm))  {
  print(i)
  nm[i] = sum(tissClust$Tissue == tissClust2$Tissue[i] & tissClust$Cluster == tissClust2$Cluster[i])
}

nm #not that bad


#sum(nm > 1500) #125
sum(nm > 400) #297
#run DSAVE on these and see how they fare, and estimate a minimal number of cells per cluster from that
library(DSAVE)
tissClust3 = tissClust2[nm>400,]
dsaveRes = vector(mode = "list", length = dim(tissClust3)[1])
for (i in 1:(dim(tissClust3)[1])) {
  print(i)
  dSub = sparseMat[,metaData$Tissue == tissClust3$Tissue[i] & metaData$Cluster == tissClust3$Cluster[i]]
  res = DSAVEGetTotalVariationPoolSize(dSub, upperBound = 50, lowerBound = 5e-1, poolSizes = c(100,200,500,750,1000,1500,2000,2500,3000,4000,5000))
  dsaveRes[[i]] = res
}

dsBulk = DSAVE::bulkTotalVar1vs1[[4]] [[2]]

saveRDS(dsaveRes, "D:/SingleCellModeling/MetAtlas/DSAVEData.RDS")
saveRDS(tissClust3, "D:/SingleCellModeling/MetAtlas/DSAVETissues.RDS")




#Set the filter manually from the figure
#We see that we need about 1,000 cells for adipose tissue on average in the figure
#translate this to UMIs per cluster, and assume it is the same for all

dAdi = sparseMat[,metaData$Tissue == "Adipose"]
dim(dAdi)#20090 83536
UMIsPerAdiCell = colSums(dAdi)
reqUMIsPerClust = mean(UMIsPerAdiCell)*1000#5009333

#Now, go through each cluster and
# 1) Check if it has enough UMIs
# 2) Check if it has
tissClust = metaData[,c(1,3)]
tissClust2 = unique(tissClust) #444 long
unique(tissClust2$Tissue)#26 tissues

gc()

UMIsPerClust = rep(NA,dim(tissClust2)[1])
UMIsPerCell = colSums(sparseMat)
for (i in 1:length(UMIsPerClust))  {
  print(i)
  sel = tissClust$Tissue == tissClust2$Tissue[i] & tissClust$Cluster == tissClust2$Cluster[i]
  UMIsPerClust[i] = sum(UMIsPerCell[sel])
}
sum(UMIsPerClust >= reqUMIsPerClust)
sel = UMIsPerClust >= reqUMIsPerClust
tissClustToModel = tissClust2[sel]
saveRDS(tissClustToModel, "D:/SingleCellModeling/MetAtlas/tissClustToModel.RDS")

#now generate the bootstrap data
set.seed(1)

for (i in 1:nrow(tissClustToModel)) {
  print(paste0(i," of ", nrow(tissClustToModel)))
  sel = tissClust$Tissue == tissClustToModel$Tissue[i] & tissClust$Cluster == tissClustToModel$Cluster[i]
  subMat = sparseMat[, sel]
  pooledBootstraps = Matrix(0, nrow = nrow(subMat), ncol=100)
  rownames(pooledBootstraps) = rownames(subMat)
  for (j in 1:100) {
    #bootstrap,i.e. sample the same number of samples with replacement 
    sel2 = sample(ncol(subMat),ncol(subMat), replace = TRUE) 
    pooledBootstraps[,j] = rowSums(subMat[,sel2])
  }
  tibbBstr = as_tibble(pooledBootstraps) %>% add_column(gene = rownames(pooledBootstraps), .before = 1)
  sampleName = paste0(tissClustToModel$Tissue[i],'_',tissClustToModel$Cluster[i]);
  fn = paste0("D:/SingleCellModeling/MetAtlas/Bootstraps/Bootstrap_", sampleName,".txt");
  write_tsv(tibbBstr, file=fn)
}

###############################
#Generate mapping to cell types
###############################

#Now the cell types
d = read_tsv("D:/SingleCellModeling/HPA/rna_single_cell_type_tissue.tsv")
dSub = unique(d[,3:5]) #444, matches perfectly with the counts file, so this most likely worked
#hmm, the tissues are named differently... annoying...
ctTissues = unique(dSub$Tissue)
mdTissues =unique(metaData$Tissue )
#so, we replace the tissues with the names in the metadata
print(tibble(ctTissues, mdTissues), n=100) #manually checked them, they come in the same order
convTab = tibble(x = ctTissues, y = mdTissues)
#remove "c-" from the cluster
print(tibble(substring(dSub$Cluster, 3, 100), dSub$Cluster),n=500)#manually checked them, ok
library(qdapTools)
d2 = tibble(tissueCT = dSub$Tissue, cluster = substring(dSub$Cluster, 3, 100), celltype = dSub$`Cell type` , tissueMD = lookup(dSub$Tissue, convTab))
print(d2,n=500) #Looks good. 
write_tsv(d2, "D:/SingleCellModeling/MetAtlas/tissue_ct.txt")
