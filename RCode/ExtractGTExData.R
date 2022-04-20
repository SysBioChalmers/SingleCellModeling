library(tidyverse)
setwd("C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling")

metaDat = read_tsv("F:/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

tissues = unique(metaDat$SMTSD)#55 types
length(metaDat$SMTSD)#22951

#grab 5 of each



RNASeq = read_tsv("F:/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", skip=2)
dim(RNASeq)#56200 17384, so, not all of the filtered are present

idsInRNASeq = colnames(RNASeq)[c(-1,-2)]


metaDatFilt = metaDat[metaDat$SAMPID %in% idsInRNASeq,]
dim(metaDatFilt)# 17382    63


#get number of samples per tissue
sampPerTissue = rep(NA, length(tissues))
for (i in 1:length(tissues)) {
  sel = metaDatFilt$SMTSD == tissues[i]
  sampPerTissue[i] = sum(sel)
}
sampPerTissue

tissuesFilt = tissues[sampPerTissue >= 5]
#select the 5 first of each type
ids = rep(NA, 5*length(tissuesFilt))
for (i in 1:length(tissuesFilt)) {
  sel = metaDatFilt$SMTSD == tissuesFilt[i]
#  print(sum(sel))
  ids[(1:5) + (i-1)*5] = metaDatFilt$SAMPID[which(sel)[1:5]];
}

library(qdapTools)
lookupTable = tibble(x = colnames(RNASeq), y = 1:length(colnames(RNASeq)))
indices = lookup(ids, lookupTable)


subMatrix = RNASeq[,c(1,2,indices)]
dim(subMatrix)#looks good
#add one column with tissue

tiss = rep(tissuesFilt, 1, each=5)

#remove version from gene
genes = substring(subMatrix[[1]], 1, str_length("ENSG00000186056"))
subMatrix[[1]] = genes

#test
idTest = colnames(subMatrix)[c(-1,-2)]
#make a loop to make it super simple and reliable
tissTest = rep("", length(idTest))
for (i in 1:length(idTest)) {
  tissTest[i] = metaDat$SMTSD[metaDat$SAMPID == idTest[i]]
}
all(tissTest == tiss) #ok

write_tsv(subMatrix, 'data/gtexIndSamp.txt')
write_tsv(tibble(tiss), 'data/gtexIndSampTissues.txt', col_names = FALSE)

#Now apply different normalization strategies than TPM, and do it on all data together

#first TMM
#we need the counts, load them:
counts = read_tsv("F:/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip=2)
dim(counts)#56200 17384, exactly the same as tpm, good
#get the count data for the ids selected above
selCounts = counts[,c(1,2,indices)]
dim(selCounts) #267, looks good
#convert to pseudo-counts
total_counts = colSums(selCounts[,c(-1,-2)]);
tpmSums = colSums(subMatrix[,c(-1,-2)])#not exactly 1 M

pseudoCounts = selCounts;
pseudoCounts[,c(-1,-2)] = t(t(subMatrix[,c(-1,-2)]) * total_counts / tpmSums)
colSums(pseudoCounts[,c(-1,-2)])[1]#48320215
total_counts[1]#48320215 , looks ok

#only keep the gene symbols to make it possible to merge with the other datasets
pseudoCounts = pseudoCounts[,-1]
colnames(pseudoCounts)[1] = "gene"
pseudoCounts[[1]]
pseudoCounts[[1]][duplicated(pseudoCounts[[1]])] #most of these genes are not of any interest, let's just merge all rows that have the same name

merged1 = pseudoCounts %>% group_by(gene) %>% summarize(across(everything(), sum))
dim(merged1)#54592   266
colSums(merged1[,-1])[1]#48320215, looks good

#Now merge with the rest.
L4Lung = read_tsv("data/L4/L4_all_pooled_lung.txt")
colnames(L4Lung)[2] = "L4Lung"
L4Spleen = read_tsv("data/L4/L4_all_pooled_spleen.txt")
colnames(L4Spleen)[2] = "L4Spleen"
LC3Data = read_tsv("data/LC3/ModelDataLC3.txt")

L4Lung2 = L4Lung %>% group_by(gene) %>% summarize(across(everything(), sum))
dim(L4Lung)
dim(L4Lung2)#removes a bunch
L4Spleen2 = L4Spleen %>% group_by(gene) %>% summarize(across(everything(), sum))
LC3Data2 = LC3Data %>% group_by(gene) %>% summarize(across(everything(), sum))
dim(LC3Data)
dim(LC3Data2)#does nothing


merged2 = inner_join(merged1, L4Lung2, by="gene")
dim(merged2)# 36718   267
merged3 = inner_join(merged2, L4Spleen2, by="gene")
dim(merged3)# 36718   268
merged4 = inner_join(merged3, LC3Data2, by="gene")
dim(merged4)# 20759   276

#Does TMM normalization using edgeR - mostly copied from the cell type profiles paper
TMMNorm <- function(ds) {
  matr = as.matrix(ds[,-1])
  #using TMM from edgeR:
  normFactors <- edgeR::calcNormFactors(matr,NULL,"TMM")
  
  libSizes = colSums(matr)
  effectiveLibSizes = libSizes * normFactors
  
  #I need to transpose the data matrix back and forth to get the
  #row wise division to work...
  matr = t(t(matr) * (10^6/ effectiveLibSizes))
  
  #There are some issues with roundoff or similar making the mean of the sum of all 
  #genes not to be exactly 10^6. Fix this:
  sumPerSamp = colSums(matr)
  meanSampSum = mean(sumPerSamp)
  
  matr = matr * 10^6 / meanSampSum
  
  ds[,-1] = matr;
  
  return (ds)
}



#now apply TMM
tmmData = TMMNorm(merged4)

#Test
#cs = colSums(tmmData[,-1])
#mean(cs)#ok
cs#pretty diverse scaling, interesting

#assemble tissues
tissAll = c(tiss, "L4 Lung", "L4 Spleen", rep("LC3",dim(LC3Data)[2]-1))
length(tissAll)#283
dim(tmmData)#284, ok, gene as well
#now save

write_tsv(tmmData, 'data/gtexIndSampTMM.txt')
write_tsv(tibble(tissAll), 'data/gtexIndSampTissuesTMM.txt', col_names = FALSE)

#also test to do quantile normalization
quantileData = merged4
#TPM first
cs = colSums(quantileData[,-1])
for(i in 1:length(cs)) {
  quantileData[,i+1] = quantileData[,i+1]/cs[i] * 10^6
}
colSums(quantileData[,-1])#looks ok

#now do the quantile
quantileData[,-1] = preprocessCore::normalize.quantiles(
  as.matrix(quantileData[,-1]))

colSums(quantileData[,-1])#hmm, not exactly 10^6, but very close

#plot(density(quantileData[[2]][quantileData[[2]] < 200] ))
#plot(density(quantileData[[2]][quantileData[[2]] < 200] ))#looks the same, ok

#save
write_tsv(quantileData, 'data/gtexIndSampQuantile.txt')
#the ids are the same
write_tsv(tibble(tissAll), 'data/gtexIndSampTissuesQuantile.txt', col_names = FALSE)


#also TMM-normalize the HCA CB per sample data
hcaPerSample = read_tsv('data/pseudoBulkModelDataCPM.txt') #here we assume that they have 1 M counts, which may be a bit of a stretch - I doubt it matters much though
hcaPerSampleTMM = TMMNorm(hcaPerSample)
mean(colSums(hcaPerSampleTMM[,-1]))#ok
colnames(hcaPerSampleTMM)[1] = "gene"

write_tsv(hcaPerSampleTMM, 'data/pseudoBulkModelDataTMM.txt')


