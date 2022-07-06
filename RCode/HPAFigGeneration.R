library(tidyverse)
memory.limit(size = 100000)

fig__dir = "D:/SingleCellModeling/MetAtlas/Figures/"

library(Matrix.utils) #for transpose of sparse matrix

dsBulk = DSAVE::bulkTotalVar1vs1[[4]] [[2]]

dsaveRes = readRDS("D:/SingleCellModeling/MetAtlas/DSAVEData.RDS")

#transform the data for ggplot
x = NULL
y = NULL
tiss = NULL
group = NULL
for (i in 1:length(dsaveRes)) {
  x = c(x,dsaveRes[[i]]$poolSizes)  
  y = c(y, dsaveRes[[i]]$Rs)
  tiss = c(tiss, rep(tissClust3$Tissue[i], length(dsaveRes[[i]]$Rs)))
  group = c(group, rep(i,length(dsaveRes[[i]]$Rs)))
}

#labels = c('Bulk', 1:(dim(tissClust3)[1]))
#labels = group
#grp = factor(group, 1:max(group), labels)
grp = as.factor(group)
tss = as.factor(tiss)

df = tibble(x=x, y=y, CT=grp, tiss=tss)

#col = c(rgb(0,0,0),gg_color_hue(6),rgb(0.5,0.5,0.5))
#color_palette = c('#000000','#888888','#FF0000','#00BB00','#DCDC55','#6B97EC','#BC976B', '#0000FF', '#444444', '#CC00CC') 

cutDf = df[df$x <= 3000,]


p1 = ggplot(cutDf, aes(x = x, y = y, group = CT, color=tiss)) +
  geom_hline(yintercept=dsBulk, linetype="dashed", color = "black")+
  geom_line() +
#  scale_linetype_manual(values = c(2,1,1,1), labels = labels) +
#  scale_color_manual(values = color_palette[c(1,2,5,3)], labels = labels) +
#  scale_size_manual(values = rep(1.3,4), labels = labels) +
  #ylim(1,1.0015) +
  #ggplot2::labs(y=expression("DSAVE score"), x=expression(Log[2]*"(number of cells)"), title="DSAVE") +
  ggplot2::labs(y=expression("DSAVE score"), x="Number of cells", title="") +
  ggplot2::theme_bw() +#+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  #ggplot2::theme(legend.title = element_blank(),legend.position="right", legend.text=element_text(color='black',size=14)) + guides(colour = guide_legend(nrow = 20), linetype = guide_legend(nrow = 4)) +
  ggplot2::theme(legend.title = element_blank(),legend.position = "none",legend.text=element_blank()) +
  ggplot2::theme(text = element_text(size=14), axis.text.x = element_text(color='black', size=14), axis.text.y = element_text(color='black', size=14)) +
  facet_wrap(~tiss, nrow=4)
p1

ggsave(
  paste0(fig__dir, "DSAVERes.png"),
  plot = p1,
  width = 15, height = 10, dpi = 300)


#for simplicity, filter on clusters having at least the number of molecules as the average of 1000 cells
#in the adipose cluster



#Set the filter manually from the figure

#Adipose:
#there is one outlier, but it will not affect the mean that much
dAdi = d[d$Tissue == "Adipose",]
dim(dAdi)
#adiClust = sort(unique(d$Cluster))
UMIsPerAdiCell = rowSums(dAdi[,c(-1,-2,-3)])
#We see that we need about 1,000 cells for adipose tissue on average in the figure
#translate this to UMIs per cluster, and assume it is the same for all
reqUMIsPerClust = mean(UMIsPerAdiCell)*1000#5009333

tissClust = d[,c(1,3)]
tissClust2 = unique(d[,c(1,3)])
unique(tissClust2$Tissue)#26 tissues

gc()

sparseMat = Matrix(as.matrix(d[,c(-1,-2,-3)]), sparse=TRUE)

UMIsPerClust = rep(NA,dim(tissClust2)[1])
UMIsPerCell = rowSums(d[,c(-1,-2,-3)]) #takes a few minutes to run
for (i in 1:length(UMIsPerClust))  {
  print(i)
  sel = tissClust$Tissue == tissClust2$Tissue[i] & tissClust$Cluster == tissClust2$Cluster[i]
  UMIsPerClust[i] = sum(UMIsPerCell[sel])
}

nm #not that bad


