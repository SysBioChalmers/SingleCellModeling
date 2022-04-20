library(R.matlab)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(Matrix)

#This code expects you to have ~48 GB of RAM
setwd("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Single-cell Modeling")

UMICounts = readRDS("D:/SingleCellModeling/LC3/GSE131907_Lung_Cancer_raw_UMI_matrix.rds") #other computer

UMICounts2part1 = as(UMICounts[,1:50000], "Matrix")#strangely, this becomes sparse automatically?
UMICounts2part2 = as(UMICounts[,50001:100000], "Matrix")
UMICounts2part3 = as(UMICounts[,100001:150000], "Matrix")
UMICounts2part4 = as(UMICounts[,150001:180000], "Matrix")
UMICounts2part5 = as(UMICounts[,180001:(dim(UMICounts)[2])], "Matrix")

rm(UMICounts)
UMICounts = cbind(UMICounts2part1,UMICounts2part2,UMICounts2part3,UMICounts2part4,UMICounts2part5)
dim(UMICounts)#29634 208506, same size as before
rm(UMICounts2part1)
rm(UMICounts2part2)
rm(UMICounts2part3)
rm(UMICounts2part4)
rm(UMICounts2part5)

saveRDS(UMICounts, "D:/SingleCellModeling/LC3/rawUMISparse.rds")
