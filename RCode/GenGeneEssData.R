
library(tidyverse)
library(ggplot2)
library(qdapTools)
library("ggpubr")
DepMapPath = "F:/DepMap/"
dataFolder = "C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/data/"


d = read_csv(paste0(DepMapPath, "CCLE_expression_full.csv"))
dim(d)#1377 52055

#transpose
ds = as_tibble(cbind(gene = names(d)[-1], t(as.matrix(d[,-1]))))
colnames(ds)[-1] = d[[1]]

#make sure the columns are numeric instead of chr
ds2 = ds
for (i in 2:ncol(ds)) {
  ds2[[i]] = as.numeric(ds[[i]])

}

#convert the gene names

pattern = ".*\\(([A-Z0-9]*)\\)"
newGenes = str_match(ds$gene, pattern)
ds2$gene = newGenes[,2]
ds2

#fix the ERCC genes
ds2$gene[is.na(newGenes[,2])] = ds$gene[is.na(newGenes[,2])];

#now convert the data. It is currently as log2(TPM + 1)
dsTPM = ds2
dsTPM[,-1] = 2^ds2[,-1] - 1
colSums(dsTPM[,-1])#very minor roundoff differences, ok

write_tsv(dsTPM, paste0(dataFolder, "DepMap_tpm_ens.txt"))

