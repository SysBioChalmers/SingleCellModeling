#devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
#library(loomR)

#BiocManager::install("rhdf5")

#lfile <- connect(filename = "D:/SingleCellModeling/L4/tissue-stability-human-lung-10XV2.loom", mode = "r+")
library(rhdf5)
library(tidyverse)

setwd("C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling")


h5h = H5Fopen("D:/SingleCellModeling/L4/tissue-stability-human-lung-10XV2.loom")
#& gives handle, $ loads object

plot(density(h5h$col_attrs$n_molecules[h5h$col_attrs$n_molecules < 5000]))
#reasonable to filter at 1000
sel = h5h$col_attrs$n_molecules >= 1000
mm = h5h&'matrix'

inds = which(sel)
length(inds)#82649

#have to read the data in chunks not to run out of memory

memory.limit(size = 70000)
m = mm[inds[1:20000],]
rowSums(m)

m1 = as(m, "sparseMatrix")
rm(m)
m = mm[inds[20001:35000],]
m2 = as(m, "sparseMatrix")
rm(m)

m = mm[inds[35001:50000],]
m3 = as(m, "sparseMatrix")
rm(m)

m = mm[inds[50001:65000],]
m4 = as(m, "sparseMatrix")
rm(m)

m = mm[inds[65001:length(inds)],]
m5 = as(m, "sparseMatrix")
rm(m)

mtot = cbind(Matrix::t(m1),Matrix::t(m2),Matrix::t(m3),Matrix::t(m4),Matrix::t(m5))
dim(mtot)#58347 82649

genes = h5h$row_attrs$Gene
length(genes)
dim(mtot)#looks ok
rownames(mtot) = genes;
cellIds = h5h$col_attrs$CellID

colnames(mtot) = cellIds[sel]

#saveRDS(mtot, file="D:/SingleCellModeling/L4/mtot_lung.RDS");

library(textTinyR)
mtotPooled = sparse_Sums(mtot, rowSums = TRUE)
sum(mtotPooled)

lungd = tibble(gene=as.character(rownames(mtot)), pool=mtotPooled)
print(lungd$gene,)

write_tsv(lungd,"data/L4/L4_all_pooled_lung.txt")

#########################
#now the spleen data
#########################

h5h = H5Fopen("D:/SingleCellModeling/L4/tissue-stability-human-spleen-10XV2.loom")
#& gives handle, $ loads object

plot(density(h5h$col_attrs$n_molecules[h5h$col_attrs$n_molecules < 5000]))
#reasonable to filter at 1000
sel = h5h$col_attrs$n_molecules >= 1000
mm = h5h&'matrix'

inds = which(sel)
length(inds)#105660

memory.limit(size = 70000)
m = mm[inds[1:20000],]

m1 = as(m, "sparseMatrix")
rm(m)
m = mm[inds[20001:35000],]
m2 = as(m, "sparseMatrix")
rm(m)

m = mm[inds[35001:50000],]
m3 = as(m, "sparseMatrix")
rm(m)

m = mm[inds[50001:65000],]
m4 = as(m, "sparseMatrix")
rm(m)

m = mm[inds[65001:80000],]
m5 = as(m, "sparseMatrix")
rm(m)

m = mm[inds[80001:95000],]
m6 = as(m, "sparseMatrix")
rm(m)

m = mm[inds[95001:length(inds)],]
m7 = as(m, "sparseMatrix")
rm(m)

mtot = cbind(Matrix::t(m1),Matrix::t(m2),Matrix::t(m3),Matrix::t(m4),Matrix::t(m5),Matrix::t(m6),Matrix::t(m7))
dim(mtot)#58347 105660

genes = h5h$row_attrs$Gene
length(genes)
dim(mtot)#looks ok
rownames(mtot) = genes;
cellIds = h5h$col_attrs$CellID

colnames(mtot) = cellIds[sel]

#saveRDS(mtot, file="D:/SingleCellModeling/L4/mtot_spleen.RDS");

library(textTinyR)
mtotPooled = sparse_Sums(mtot, rowSums = TRUE)
sum(mtotPooled)

spleend = tibble(gene=as.character(rownames(mtot)), pool=mtotPooled)
print(spleend$gene,)

write_tsv(spleend,"data/L4/L4_all_pooled_spleen.txt")



