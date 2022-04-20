library(R.matlab)
library(tidyverse)
library(ggplot2)
library(ggpubr)


fig____path = "Z:/projects/Single-cell modeling/figures/"

setwd("C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#####################################
#first run DSAVE
#####################################
#devtools::install_github("SysBioChalmers/DSAVE-R")
library(DSAVE)
library(Seurat)
set.seed(1)
d = Read10X(data.dir = "data/exportToR/datasets/hcat", gene.column=1)
varHcat = DSAVEGetTotalVariationPoolSize(d,upperBound = 50, lowerBound = 5e-1)
d = Read10X(data.dir = "data/exportToR/datasets/hcab", gene.column=1)
varHcab = DSAVEGetTotalVariationPoolSize(d,upperBound = 50, lowerBound = 5e-1)
d = Read10X(data.dir = "data/exportToR/datasets/t68k", gene.column=1)
varT68k = DSAVEGetTotalVariationPoolSize(d,upperBound = 50, lowerBound = 5e-1)
d = Read10X(data.dir = "data/exportToR/datasets/lccm", gene.column=1)
varLccm = DSAVEGetTotalVariationPoolSize(d,upperBound = 50, lowerBound = 5e-1)
d = Read10X(data.dir = "data/exportToR/datasets/lcct", gene.column=1)
varLcct = DSAVEGetTotalVariationPoolSize(d,upperBound = 50, lowerBound = 5e-1)
d = Read10X(data.dir = "data/exportToR/datasets/tcd8", gene.column=1)
varTcd8 = DSAVEGetTotalVariationPoolSize(d,upperBound = 50, lowerBound = 5e-1)

#load the melt matrix
d = as.matrix(read_csv("data/exportToR/datasets/melt/matrix.txt", col_names = FALSE))
genes = read_csv("data/exportToR/datasets/melt/genes.txt", col_names = FALSE)
rownames(d) = genes[[1]]
varMelt = DSAVEGetTotalVariationPoolSize(d,upperBound = 50, lowerBound = 5e-1)

dsBulk = DSAVE::bulkTotalVar1vs1[[4]] [[2]]

#transform the data for ggplot
x = log2(c(0,10000, varMelt$poolSizes, varHcat$poolSizes, varHcab$poolSizes, varT68k$poolSizes, varLccm$poolSizes, varLcct$poolSizes, varTcd8$poolSizes))
y = c(dsBulk, dsBulk, varMelt$Rs, varHcat$Rs, varHcab$Rs, varT68k$Rs, varLccm$Rs, varLcct$Rs, varTcd8$Rs)
labels = c("Bulk T", "Mel T (Sm2)", "HCA CB T", "HCA CB B", "PBMC68k T", "LC T", "LC Mac", "TCD8 T")
group = factor(c(1,1,rep(2, length(varMelt$Rs)),rep(3, length(varHcat$Rs)),rep(4, length(varHcab$Rs)),rep(5, length(varT68k$Rs)),rep(6, length(varLccm$Rs)),rep(7, length(varLcct$Rs)),rep(8, length(varTcd8$Rs))),1:8, labels)


df = tibble(x=x, y=y, Dataset=group)

col = c(rgb(0,0,0),gg_color_hue(6),rgb(0.5,0.5,0.5))

cutDf = df#[df$y >= 0.6,]


pB = ggplot(cutDf, aes(x = x, y = y, color=Dataset, linetype=Dataset, size=Dataset)) +
  geom_line() +
  scale_linetype_manual(values = c(2,1,1,1,1,1,1,1), labels = labels) +
  scale_color_manual(values = col[c(1,1,2,3,4,5,6,7)], labels = labels) +
  scale_size_manual(values = rep(1,8), labels = labels) +
  #ylim(1,1.0015) +
  ggplot2::labs(y=expression("DSAVE score"), x=expression(Log[2]*"(number of cells)"), title="DSAVE") +
  ggplot2::theme_bw() +#+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  ggplot2::theme(legend.title = element_blank(),legend.position="bottom", legend.text=element_text(color='black',size=14)) + guides(colour = guide_legend(nrow = 4), size = guide_legend(nrow = 4), linetype = guide_legend(nrow = 4)) +
  ggplot2::theme(text = element_text(size=14), axis.text.x = element_text(color='black', size=14), axis.text.y = element_text(color='black', size=14))
pB


#################################
# Model distance using tINIT 3
#################################

d = readMat("data/PoolSizeModelRes.mat")
hcatData = as.numeric(d$d[,,1]$hcat)
t68kData = as.numeric(d$d[,,1]$t68k)
bulkVal = as.numeric(d$d[,,1]$bulk)
#also get the x values
d = readMat("data/PSPoolSizes.mat")
hcatX = d$d[,,1]$poolSizesHcat
t68kX = d$d[,,1]$poolSizesT68k

x = log2(c(hcatX,t68kX))
x = c(min(x),max(x),x)
y = c(bulkVal,bulkVal,hcatData,t68kData)
labels = c( "Bulk T", "HCA CB T", "PBMC68k T")
datas = factor(c(rep(0,2),rep(1,length(hcatX)), rep(2,length(t68kX))), c(0,1,2), labels)
col = c(rgb(0,0,0),gg_color_hue(6),rgb(0.5,0.5,0.5))

df = tibble(x=x, y=y, Dataset = datas)
pA = ggplot(df, aes(x = x, y = y, color=Dataset, linetype=Dataset, size=Dataset)) +
  geom_line() +
  scale_linetype_manual(values = c(2,1,1), labels = labels) +
  scale_color_manual(values = col[c(1,2,4)], labels = labels) +
  scale_size_manual(values = rep(1,2,3), labels = labels) +
  #ylim(1,1.0015) +
  ggplot2::labs(y=expression("Jaccard index"), x=expression(Log[2]*"(number of cells)"), title="tINIT 3") +
  ggplot2::theme_bw() +#+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  ggplot2::theme(legend.title = element_blank(),legend.position="bottom", legend.text=element_text(color='black', size=14)) + guides(colour = guide_legend(nrow = 2), size = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  ggplot2::theme(text = element_text(size=14), axis.text.x = element_text(color='black', size=14), axis.text.y = element_text(color='black', size=14))
pA




##################
# Contaminated samples

poolCont = readMat("data/PSresPoolSizeVsReactionScoresJaccardCont.mat")

X10x = poolCont$resPoolSizeVsReactionScoresJaccardCont[,,1]$X
Y10x = poolCont$resPoolSizeVsReactionScoresJaccardCont[,,1]$Y

#transform the data for ggplot
x = log2(c(as.vector(X10x[,1]), as.vector(X10x[,2]), as.vector(X10x[,3]), as.vector(X10x[,4]), as.vector(X10x[,5])))
y = c(as.vector(Y10x[,1]), as.vector(Y10x[,2]), as.vector(Y10x[,3]), as.vector(Y10x[,4]), as.vector(Y10x[,5]))
labels = c("0%", "2%", "5%", "10%", "20%")
group = factor(rep(1:5, 1, each = nrow(X10x)), 1:5, labels)


df = tibble(x=x, y=y, Contamination=group)

col = c(rgb(0,0,0),gg_color_hue(6),rgb(0.5,0.5,0.5))

cutDf = df[df$y >= 0.8,]


pC = ggplot(cutDf, aes(x = x, y = y, color=Contamination, linetype=Contamination, size=Contamination)) +
  geom_line() +
  scale_linetype_manual(values = c(1,1,1,1,1), labels = labels) +
  scale_color_manual(values = col[c(1,2,3,4,5)], labels = labels) +
  scale_size_manual(values = rep(1,5), labels = labels) +
  #ylim(1,1.0015) +
  ggplot2::labs(y=expression("Jaccard index"), x=expression(Log[2]*"(number of cells)"), title="Contamination") +
  ggplot2::theme_bw() + #+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  ggplot2::theme(legend.title = element_blank(),legend.position="bottom", legend.text=element_text(color='black', size=14)) + guides(colour = guide_legend(nrow = 3), size = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3)) +
  ggplot2::theme(text = element_text(size=14), axis.text.x = element_text(color='black', size=14), axis.text.y = element_text(color='black', size=14))

pC


#save pA, pB and pC as separate figures - build that figure separately

ggsave(
  paste0(fig____path, "Fig2A.png"),
  plot = pA,
  width = 4, height = 4.03, dpi = 300)

ggsave(
  paste0(fig____path, "Fig2B.png"),
  plot = pB,
  width = 4, height = 4.5, dpi = 300)

ggsave(
  paste0(fig____path, "Fig2C.png"),
  plot = pC,
  width = 4, height = 4.26, dpi = 300)

###################################################
#Variation across samples supplementary
###################################################

d = readMat("data/PoolSizeModelRes.mat")
bulk = as.numeric(d$d[,,1]$bulkAll)
pseudoBulkTPM = as.numeric(d$d[,,1]$pseudoBulkAll)
pseudoBulkTMM = as.numeric(d$d[,,1]$pseudoBulkTMMAll)
bulkVal = as.numeric(d$d[,,1]$bulk)


y = c(bulk, pseudoBulkTPM, pseudoBulkTMM)
names = c("Bulk TMM", "Pseudo-bulk CPM", "Pseudo-bulk TMM")
group = factor(rep(c(1,2,3),1,each=28), 1:3, as.character(names)[1:3])  
df = tibble(y=y, group=group)

color_palette <- c('#B5D39B','#E7B56C','#6B97BC')  # light green, light yellow, light blue
pX = ggplot(df, aes(x = group, y = y, fill = group)) +
  geom_violin(trim=F, show.legend=F, scale='count') +
  scale_fill_manual(values=color_palette) +
  theme_classic() + 
  ylab("Jaccard index") +
  xlab("") +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,
                                   color='black', size=12),
        axis.text.y = element_text(color='black', size=12),
        axis.line.x = element_blank())
pX

ggsave(
  paste0(fig____path, "Fig2SupSampleVar.png"),
  plot = pX,
  width = 4.5, height = 4.5, dpi = 300)




###################################################
#Structural comparisons
###################################################
structuralTPM = readMat("data/structCompTPM.mat")
acrossGTExTPM = as.numeric(structuralTPM$d[,,1]$acrossGTEx)
withinGTExTPM = as.numeric(structuralTPM$d[,,1]$withinGTEx)
withinScLungTPM = as.numeric(structuralTPM$d[,,1]$withinScLung)
scVsGTExLungTPM = as.numeric(structuralTPM$d[,,1]$scVsGTExLung)
scVsGTExAnyTPM = as.numeric(structuralTPM$d[,,1]$scVsGTExAny);

#first violin plot for TPM
##########################

y = c(acrossGTExTPM, withinGTExTPM, withinScLungTPM, scVsGTExLungTPM, scVsGTExAnyTPM)
comparison = factor(c(rep(0,length(acrossGTExTPM)),rep(1,length(withinGTExTPM)),rep(2,length(withinScLungTPM)),rep(3,length(scVsGTExLungTPM)),rep(4,length(scVsGTExAnyTPM))),0:4,c("Across bulk tissues", "Within bulk tissues", "Within LC3 immune", "Immune across tech", "Across tissues and tech"))

df = tibble(y=y, Comparison=comparison)

color_palette <- c('#B5D39B','#E7B56C','#6B97BC','#BC976B','#BC556B')  # light green, light yellow, light blue
#  pdf('GEA_compare_MCC.pdf', width=3.7, height=3.8)
pI = ggplot(df, aes(x = Comparison, y = y, fill = Comparison)) +
  geom_violin(trim=F, show.legend=F, scale='width') +
  scale_fill_manual(values=color_palette) +
  theme_classic() +
  ggtitle("Structural variation") +
  ylab("Jaccard index") +
  xlab('') +
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=30, hjust=1, vjust=1,
                                   color='black', size=14),
        axis.text.y = element_text(color='black', size=14),
        axis.line.x = element_blank(),
        plot.margin = margin(0,0,0,1.3, "cm")) #+
pI

ggsave(
  paste0(fig____path, "Fig2E.png"),
  plot = pI,
  width = 4, height = 4, dpi = 300)


#Then comparison of normalizations
structuralTMM = readMat("data/structCompTMM.mat")
acrossGTExTMM = as.numeric(structuralTMM$d[,,1]$acrossGTEx)
withinGTExTMM = as.numeric(structuralTMM$d[,,1]$withinGTEx)
withinScLungTMM = as.numeric(structuralTMM$d[,,1]$withinScLung)
scVsGTExLungTMM = as.numeric(structuralTMM$d[,,1]$scVsGTExLung)
scVsGTExAnyTMM = as.numeric(structuralTMM$d[,,1]$scVsGTExAny);

structuralQ = readMat("data/structCompQuantile.mat")
acrossGTExQ = as.numeric(structuralQ$d[,,1]$acrossGTEx)
withinGTExQ = as.numeric(structuralQ$d[,,1]$withinGTEx)
withinScLungQ = as.numeric(structuralQ$d[,,1]$withinScLung)
scVsGTExLungQ = as.numeric(structuralQ$d[,,1]$scVsGTExLung)
scVsGTExAnyQ = as.numeric(structuralQ$d[,,1]$scVsGTExAny);

y = c(acrossGTExTPM, withinGTExTPM, withinScLungTPM, scVsGTExLungTPM, scVsGTExAnyTPM,acrossGTExTMM, withinGTExTMM, withinScLungTMM, scVsGTExLungTMM, scVsGTExAnyTMM,acrossGTExQ, withinGTExQ, withinScLungQ, scVsGTExLungQ, scVsGTExAnyQ)
comparison = factor(rep(c(rep(0,length(acrossGTExTPM)),rep(1,length(withinGTExTPM)),rep(2,length(withinScLungTPM)),rep(3,length(scVsGTExLungTPM)),rep(4,length(scVsGTExAnyTPM))),3),0:4,c("Across bulk tissues", "Within bulk tissues", "Within LC3 immune", "Immune across tech", "Across tissues and tech"))
ltot = length(acrossGTExTPM) + length(withinGTExTPM) + length(withinScLungTPM) + length(scVsGTExLungTPM) + length(scVsGTExAnyTPM)
normalization = factor(c(rep(0,ltot),rep(1,ltot),rep(2,ltot)),0:2,c("TPM", "TMM", "Quantile"))

df = tibble(y=y, Comparison=comparison, Normalization = normalization)

color_palette <- c('#B5D39B','#E7B56C','#6B97BC')  # light green, light yellow, light blue

pJ = ggplot(data=df, aes(x=Comparison, y=y, fill=Normalization)) +
  geom_boxplot(outlier.shape = NA) +
  labs( y="Jaccard index", x="", title = "Effect of normalization") +
  scale_fill_manual(values=color_palette) +
  theme_classic() +
  coord_cartesian(ylim = c(0.805, 1)) +
  #ylim(0.82,1) +
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=30, hjust=1, vjust=1,
                                   color='black', size=14),
        axis.text.y = element_text(color='black', size=14),
        axis.line.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        plot.margin = margin(0,0,0,1.3, "cm")) #+
pJ

ggsave(
  paste0(fig____path, "Fig2F.png"),
  plot = pJ,
  width = 4.4, height = 5, dpi = 300)






###################################################
#t-sne Figs
###################################################

noNormData = readMat("data/allTSNE.mat")
x = as.numeric(noNormData$d[,,1]$tsneX)
y = as.numeric(noNormData$d[,,1]$tsneY)

#use the TMM tissues file, it is the same
tissues = read_tsv('data/gtexIndSampTissuesTMM.txt', col_names = FALSE)[[1]]


lungImmune = c(1,3,5,6,7,8,9,10,11,12,13,14,16)
lungEpithelial = c(2,4,15)


selLung = tissues == "Lung"
selBrain = startsWith(tissues,"Brain")
selSpleen = tissues == "Spleen"
selBlood = tissues == "Whole Blood"
selL4Lung = tissues == "L4 Lung"
selL4Spleen = tissues == "L4 Spleen"
selscLC = tissues == "LC3"
selLC3Immune = which(selscLC)[lungImmune]
selLC3Epi = which(selscLC)[lungEpithelial]

tissType = rep(4,length(tissues))
tissType[selLung] = 0
tissType[selBrain] = 1
tissType[selSpleen] = 2
tissType[selBlood] = 3
tissType[selL4Lung] = 5
tissType[selL4Spleen] = 6
#tissType[selscLC] = 6
tissType[selLC3Immune] = 7
tissType[selLC3Epi] = 8
labels = c("GTEx lung","GTEx brain","GTEx spleen","GTEx blood","GTEx other","L4 sc lung","L4 sc spleen","LC3 sc immune", "LC3 sc epithelial")
tt = factor(tissType, 0:8, labels)
df = tibble(x=x, y=y, tissue = tt)

#color_palette = c('#B5D39B','#E7B56C','#6B97BC','#BC976B','#BC556B')  
color_palette = c('#202020','#DC556B','#6B97EC','#BC976B','#55BC6B')  


pF = ggplot(df, aes(x = x, y = y, color=tissue, shape=tissue)) +
  geom_point(size=2, stroke = 2) +
  #scale_linetype_manual(values = c(1,1,1,1,1), labels = labels) +
  scale_color_manual(values = color_palette[c(1,2,3,4,5,1,3,1,1)], labels = labels) +
  scale_shape_manual(values = c(19,19,19,19,19,0,0,2,3), labels=labels) +
  #scale_size_manual(values = rep(20,7), labels = labels) +
  #ylim(1,1.0015) +
  ggplot2::labs(y=expression("t-SNE y"), x="t-SNE x", title="Structural comparison") +
  ggplot2::theme_bw() + #+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  ggplot2::theme(legend.title = element_blank()) +
  ggplot2::theme(legend.title = element_blank(),legend.position="bottom", legend.text=element_text(size=14)) + guides(colour = guide_legend(nrow = 3), size = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3)) +
  ggplot2::theme(text = element_text(size=14),
                 axis.text.x = element_text(color='black', size=14),
                 axis.text.y = element_text(color='black', size=14))
pF

ggsave(
  paste0(fig____path, "Fig2D.png"),
  plot = pF,
  width = 8, height = 9, dpi = 300)

#ggsave(
#  paste0(fig____path, "Fig2D.svg"),
#  plot = pF,
#  width = 8, height = 9, dpi = 300)


##############################
# TMM
##############################

tmmData = readMat("data/tmmTSNE.mat")
x = as.numeric(tmmData$d[,,1]$tsneX)
y = as.numeric(tmmData$d[,,1]$tsneY)
tissues = read_tsv('data/gtexIndSampTissuesTMM.txt', col_names = FALSE)[[1]]

lungImmune = c(1,3,5,6,7,8,9,10,11,12,13,14,16)
lungEpithelial = c(2,4,15)


selLung = tissues == "Lung"
selBrain = startsWith(tissues,"Brain")
selSpleen = tissues == "Spleen"
selBlood = tissues == "Whole Blood"
selL4Lung = tissues == "L4 Lung"
selL4Spleen = tissues == "L4 Spleen"
selscLC = tissues == "LC3"
selLC3Immune = which(selscLC)[lungImmune]
selLC3Epi = which(selscLC)[lungEpithelial]

tissType = rep(4,length(tissues))
tissType[selLung] = 0
tissType[selBrain] = 1
tissType[selSpleen] = 2
tissType[selBlood] = 3
tissType[selL4Lung] = 5
tissType[selL4Spleen] = 6
#tissType[selscLC] = 6
tissType[selLC3Immune] = 7
tissType[selLC3Epi] = 8
labels = c("GTEx lung","GTEx brain","GTEx spleen","GTEx blood","GTEx other","L4 sc lung","L4 sc spleen","LC3 sc immune", "LC3 sc epithelial")
tt = factor(tissType, 0:8, labels)
df = tibble(x=x, y=y, tissue = tt)

#color_palette = c('#B5D39B','#E7B56C','#6B97BC','#BC976B','#BC556B')  
color_palette = c('#202020','#DC556B','#6B97EC','#BC976B','#55BC6B')  


pG = ggplot(df, aes(x = x, y = y, color=tissue, shape=tissue)) +
  geom_point(size=2, stroke = 2) +
  #scale_linetype_manual(values = c(1,1,1,1,1), labels = labels) +
  scale_color_manual(values = color_palette[c(1,2,3,4,5,1,3,1,1)], labels = labels) +
  scale_shape_manual(values = c(19,19,19,19,19,0,0,2,3), labels=labels) +
  #scale_size_manual(values = rep(20,7), labels = labels) +
  #ylim(1,1.0015) +
  ggplot2::labs(y=expression("t-SNE y"), x="t-SNE x", title="") +
  ggplot2::theme_bw() + #+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  ggplot2::theme(legend.title = element_blank()) +
  ggplot2::theme(legend.title = element_blank(),legend.position="bottom", legend.text=element_text(size=14)) + guides(colour = guide_legend(nrow = 3), size = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3)) +
  ggplot2::theme(text = element_text(size=14),
                 axis.text.x = element_text(color='black', size=14),
                 axis.text.y = element_text(color='black', size=14))
pG

ggsave(
  paste0(fig____path, "Fig2 Sup TMM.png"),
  plot = pG,
  width = 8, height = 9, dpi = 300)

##############################
# Quantile
##############################

qData = readMat("data/quantileTSNE.mat")
x = as.numeric(qData$d[,,1]$tsneX)
y = as.numeric(qData$d[,,1]$tsneY)
tissues = read_tsv('data/gtexIndSampTissuesQuantile.txt', col_names = FALSE)[[1]]

lungImmune = c(1,3,5,6,7,8,9,10,11,12,13,14,16)
lungEpithelial = c(2,4,15)


selLung = tissues == "Lung"
selBrain = startsWith(tissues,"Brain")
selSpleen = tissues == "Spleen"
selBlood = tissues == "Whole Blood"
selL4Lung = tissues == "L4 Lung"
selL4Spleen = tissues == "L4 Spleen"
selscLC = tissues == "LC3"
selLC3Immune = which(selscLC)[lungImmune]
selLC3Epi = which(selscLC)[lungEpithelial]

tissType = rep(4,length(tissues))
tissType[selLung] = 0
tissType[selBrain] = 1
tissType[selSpleen] = 2
tissType[selBlood] = 3
tissType[selL4Lung] = 5
tissType[selL4Spleen] = 6
#tissType[selscLC] = 6
tissType[selLC3Immune] = 7
tissType[selLC3Epi] = 8
labels = c("GTEx lung","GTEx brain","GTEx spleen","GTEx blood","GTEx other","L4 sc lung","L4 sc spleen","LC3 sc immune", "LC3 sc epithelial")
tt = factor(tissType, 0:8, labels)
df = tibble(x=x, y=y, tissue = tt)

#color_palette = c('#B5D39B','#E7B56C','#6B97BC','#BC976B','#BC556B')  
color_palette = c('#202020','#DC556B','#6B97EC','#BC976B','#55BC6B')  


pH = ggplot(df, aes(x = x, y = y, color=tissue, shape=tissue)) +
  geom_point(size=2, stroke = 2) +
  #scale_linetype_manual(values = c(1,1,1,1,1), labels = labels) +
  scale_color_manual(values = color_palette[c(1,2,3,4,5,1,3,1,1)], labels = labels) +
  scale_shape_manual(values = c(19,19,19,19,19,0,0,2,3), labels=labels) +
  #scale_size_manual(values = rep(20,7), labels = labels) +
  #ylim(1,1.0015) +
  ggplot2::labs(y=expression("t-SNE y"), x="t-SNE x", title="") +
  ggplot2::theme_bw() + #+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  ggplot2::theme(legend.title = element_blank()) +
  ggplot2::theme(legend.title = element_blank(),legend.position="bottom", legend.text=element_text(size=14)) + guides(colour = guide_legend(nrow = 3), size = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3)) +
  ggplot2::theme(text = element_text(size=14),
                 axis.text.x = element_text(color='black', size=14),
                 axis.text.y = element_text(color='black', size=14))
pH

ggsave(
  paste0(fig____path, "Fig2 Sup Quantile.png"),
  plot = pH,
  width = 8, height = 9, dpi = 300)


############################
#Statistical test for on/off
############################

dat <- data.frame(
  "ct_1" = c(1, 99),
  "ct_2" = c(99, 1),
  row.names = c("On", "Off"),
  stringsAsFactors = FALSE
)
dat
test <- fisher.test(dat)
test#p-value < 2.2e-16

dat <- data.frame(
  "ct_1" = c(0, 8),
  "ct_2" = c(8, 0),
  row.names = c("On", "Off"),
  stringsAsFactors = FALSE
)
dat
test <- fisher.test(dat)
test#p-value 0.0001554, not too bad


dat <- data.frame(
  "ct_1" = c(0, 6),
  "ct_2" = c(6, 0),
  row.names = c("On", "Off"),
  stringsAsFactors = FALSE
)
dat
test <- fisher.test(dat)
test


















