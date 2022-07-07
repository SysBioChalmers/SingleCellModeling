#Plots the DSAVE score for all tissues
library(tidyverse)
library(DSAVE)

fig__dir = "D:/SingleCellModeling/MetAtlas/Figures/"

library(Matrix.utils) #for transpose of sparse matrix

dsBulk = DSAVE::bulkTotalVar1vs1[[4]] [[2]]

dsaveRes = readRDS("D:/SingleCellModeling/MetAtlas/DSAVEData.RDS")
tissClust3 = readRDS("D:/SingleCellModeling/MetAtlas/DSAVETissues.RDS")

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

grp = as.factor(group)
tss = as.factor(tiss)

df = tibble(x=x, y=y, CT=grp, tiss=tss)

cutDf = df[df$x <= 3000,]

p1 = ggplot(cutDf, aes(x = x, y = y, group = CT, color=tiss)) +
  geom_hline(yintercept=dsBulk, linetype="dashed", color = "black")+
  geom_line() +
  scale_x_continuous(breaks = c(1,1000,2000,3000), labels = c("0", "1k", "2k", "3k")) +
  ggplot2::labs(y=expression("DSAVE score"), x="Number of cells", title="") +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  ggplot2::theme(legend.title = element_blank(),legend.position = "none",legend.text=element_blank()) +
  ggplot2::theme(text = element_text(size=14), axis.text.x = element_text(color='black', size=14), axis.text.y = element_text(color='black', size=14)) +
  facet_wrap(~tiss, nrow=4)
p1

ggsave(
  paste0(fig__dir, "DSAVEResMetAtlas.png"),
  plot = p1,
  width = 15, height = 10, dpi = 300)
#install.packages("svglite")
ggsave(
  paste0(fig__dir, "DSAVEResMetAtlas.svg"),
  plot = p1,
  width = 15, height = 10, dpi = 300)
