library(R.matlab)
library(tidyverse)
library(ggplot2)
library(ggpubr)


fig____path = "Z:/projects/Single-cell modeling/figures/"


setwd("C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling")

###################################################
#New tINIT
###################################################

newtINITData = readMat("data/tsneGTExNewtINIT.mat")
x = as.numeric(newtINITData$d[,,1]$tsneX)
y = as.numeric(newtINITData$d[,,1]$tsneY)
tissues = as.character(unlist(newtINITData$d[,,1]$tissues))


selLung = tissues == "Lung"
selBrain = startsWith(tissues,"Brain")
selSpleen = tissues == "Spleen"
selBlood = tissues == "Whole Blood"
selLiver = tissues == "Liver"
selKidney = startsWith(tissues,"Kidney")
selColon = startsWith(tissues,"Colon")


tissType = rep(7,length(tissues))
tissType[selLung] = 0
tissType[selBrain] = 1
tissType[selSpleen] = 2
tissType[selBlood] = 3
tissType[selLiver] = 4
tissType[selKidney] = 5
tissType[selColon] = 6
labels = c("Lung","Brain","Spleen","Blood","Liver","Kidney","Colon","Other")
tt = factor(tissType, 0:7, labels)
df = tibble(x=x, y=y, tissue = tt)

#color_palette = c('#B5D39B','#E7B56C','#6B97BC','#BC976B','#BC556B')  
color_palette = c('#202020','#DC556B','#6B97EC','#BC976B','#55BC6B')  


pNew = ggplot(df, aes(x = x, y = y, color=tissue, shape=tissue)) +
  geom_point(size=2, stroke = 2) +
  #scale_linetype_manual(values = c(1,1,1,1,1), labels = labels) +
  scale_color_manual(values = color_palette[c(1,2,3,4,1,2,3,4)], labels = labels) +
  scale_shape_manual(values = c(19,19,19,19,0,0,0,0), labels=labels) +
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
pNew

ggsave(
  paste0(fig____path, "Fig1SupNewtINIT.png"),
  plot = pNew,
  width = 8, height = 9, dpi = 300)

ggsave(
  paste0(fig____path, "Fig1SupNewtINIT.svg"),
  plot = pNew,
  width = 8, height = 9, dpi = 300)

###################################################
#Old tINIT
###################################################

oldtINITData = readMat("data/tsneGTExOldtINIT.mat")
x = as.numeric(oldtINITData$d[,,1]$tsneX)
y = as.numeric(oldtINITData$d[,,1]$tsneY)
tissues = as.character(unlist(oldtINITData$d[,,1]$tissues))


selLung = tissues == "Lung"
selBrain = startsWith(tissues,"Brain")
selSpleen = tissues == "Spleen"
selBlood = tissues == "Whole Blood"
selLiver = tissues == "Liver"
selKidney = startsWith(tissues,"Kidney")
selColon = startsWith(tissues,"Colon")


tissType = rep(7,length(tissues))
tissType[selLung] = 0
tissType[selBrain] = 1
tissType[selSpleen] = 2
tissType[selBlood] = 3
tissType[selLiver] = 4
tissType[selKidney] = 5
tissType[selColon] = 6
labels = c("Lung","Brain","Spleen","Blood","Liver","Kidney","Colon","Other")
tt = factor(tissType, 0:7, labels)
df = tibble(x=x, y=y, tissue = tt)

#color_palette = c('#B5D39B','#E7B56C','#6B97BC','#BC976B','#BC556B')  
color_palette = c('#202020','#DC556B','#6B97EC','#BC976B','#55BC6B')  


pOld = ggplot(df, aes(x = x, y = y, color=tissue, shape=tissue)) +
  geom_point(size=2, stroke = 2) +
  #scale_linetype_manual(values = c(1,1,1,1,1), labels = labels) +
  scale_color_manual(values = color_palette[c(1,2,3,4,1,2,3,4)], labels = labels) +
  scale_shape_manual(values = c(19,19,19,19,0,0,0,0), labels=labels) +
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
pOld

ggsave(
  paste0(fig____path, "Fig1SupOldtINIT.png"),
  plot = pOld,
  width = 8, height = 9, dpi = 300)

ggsave(
  paste0(fig____path, "Fig1SupOldtINIT.svg"),
  plot = pOld,
  width = 8, height = 9, dpi = 300)

