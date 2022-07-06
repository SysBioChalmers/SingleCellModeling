#This code extracts the data from a huge tsv file and saves it as a sparse matrix
library(tidyverse)
memory.limit(size = 100000)

fig__dir = "D:/SingleCellModeling/MetAtlas/Figures/"

#use the data.table package, it is better for reading large datasets
#install.packages("data.table")
library(data.table)
Sys.setenv("VROOM_CONNECTION_SIZE" = 500000) #needed to be able to read long lines

d = fread("D:/SingleCellModeling/HPA/rna_single_cell_read_count.tsv")
dim(d) # 566108  20093
library(Matrix.utils) #for transpose of sparse matrix

#need to convert to sparse in chunks - will run out of memory otherwise...
gc()
sparseMat1 = Matrix(as.matrix(d[1:100000,c(-1,-2,-3)]), sparse=TRUE)
gc()
sparseMat2 = Matrix(as.matrix(d[100001:200000,c(-1,-2,-3)]), sparse=TRUE)
gc()
sparseMat3 = Matrix(as.matrix(d[200001:300000,c(-1,-2,-3)]), sparse=TRUE)
gc()
sparseMat4 = Matrix(as.matrix(d[300001:400000,c(-1,-2,-3)]), sparse=TRUE)
gc()
sparseMat5 = Matrix(as.matrix(d[400001:500000,c(-1,-2,-3)]), sparse=TRUE)
gc()
sparseMat6 = Matrix(as.matrix(d[500001:(dim(d)[1]),c(-1,-2,-3)]), sparse=TRUE)
gc()

metaData = d[,1:3]
rm(d)
gc()

#merge into 1 matrix
sparseMat12 = rbind2(sparseMat1,sparseMat2)
rm(sparseMat1)
rm(sparseMat2)
gc()
sparseMat45 = rbind2(sparseMat4,sparseMat5)
rm(sparseMat4)
rm(sparseMat5)
gc()
sparseMat123 = rbind2(sparseMat12,sparseMat3)
rm(sparseMat12)
rm(sparseMat3)
gc()
sparseMat12345 = rbind2(sparseMat123,sparseMat45)
rm(sparseMat123)
rm(sparseMat45)
gc()
sparseMat123456 = rbind2(sparseMat12345,sparseMat6)
rm(sparseMat12345)
rm(sparseMat6)
gc()
dim(sparseMat123456)
sparseMat = sparseMat123456
rm(sparseMat123456)


sparseMat2 = t(sparseMat)
dim(sparseMat2)
sparseMat = sparseMat2
colnames(sparseMat) = metaData$Cell
saveRDS(sparseMat, "D:/SingleCellModeling/MetAtlas/SparseMat.RDS")
saveRDS(metaData, "D:/SingleCellModeling/MetAtlas/metaData.RDS")

#Some tests to verify that the matrix is correct (The numbers were looked up manually in the text file)
#test one number from each chunk
sparseMat[7,2] == 47 #ok
sparseMat[3,119510] == 4 #ok
sparseMat[7,207622] == 2 #ok
sparseMat[7,312933] == 2 #ok
sparseMat[1,406090] == 11 #ok
sparseMat[7,521545] == 3 #ok
#all looks ok!

rm(sparseMat)
rm(metaData)
gc()

