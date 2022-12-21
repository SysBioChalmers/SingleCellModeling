library(ggplot2)
library(scales) # for muted function
library(tidyverse)

figPath = "Z:/projects/Single-cell modeling/figures/"

setwd("C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling")

library(R.matlab)

##################################
# Structural comparison
##################################

maTsne = readMat("data/MetAtlas/metAtlasTsne.mat")$d[,,1] 

cellTypesOrig = unlist(maTsne$celltype)
tissues = unlist(maTsne$tissue)

#tsne plot
x = as.numeric(maTsne$tsneX)
y = as.numeric(maTsne$tsneY)


#make the cell types more general, into families such as immune etc.

print(tibble(ct = cellTypesOrig, tissue = tissues), n=500)
cellTypes = cellTypesOrig
cellTypes[cellTypes == "plasma cells"] = "lymphocytes"
cellTypes[cellTypes == "t-cells"] = "lymphocytes"
cellTypes[cellTypes == "macrophages"] = "other immune cells"
cellTypes[cellTypes == "dendritic cells"] = "other immune cells"
cellTypes[cellTypes == "mixed immune cells"] = "other immune cells"
cellTypes[cellTypes == "undifferentiated cells"] = "organ-spec. cells"
cellTypes[cellTypes == "oligodendrocytes"] = "organ-spec. cells"
cellTypes[cellTypes == "breast glandular cells"] = "organ-spec. cells"
cellTypes[cellTypes == "breast myoepithelial cells"] = "organ-spec. cells"
cellTypes[cellTypes == "basal respiratory cells"] = "organ-spec. cells"
cellTypes[cellTypes == "respiratory ciliated cells"] = "organ-spec. cells"
cellTypes[cellTypes == "club cells"] = "organ-spec. cells"
cellTypes[cellTypes == "glandular and luminal cells"] = "organ-spec. cells"
cellTypes[cellTypes == "endometrial stromal cells"] = "organ-spec. cells"
cellTypes[cellTypes == "endometrial ciliated cells"] = "organ-spec. cells"
cellTypes[cellTypes == "squamous epithelial cells"] = "organ-spec. cells"
cellTypes[cellTypes == "proximal tubular cells"] = "organ-spec. cells"
cellTypes[cellTypes == "hepatocytes"] = "organ-spec. cells"
cellTypes[cellTypes == "granulosa cells"] = "organ-spec. cells"
cellTypes[cellTypes == "theca cells"] = "organ-spec. cells"
cellTypes[cellTypes == "ductal cells"] = "organ-spec. cells"
cellTypes[cellTypes == "pancreatic endocrine cells"] = "organ-spec. cells"
cellTypes[cellTypes == "mixed cell types"] = "organ-spec. cells"
cellTypes[cellTypes == "exocrine glandular cells"] = "organ-spec. cells"
cellTypes[cellTypes == "cytotrophoblasts"] = "organ-spec. cells"
cellTypes[cellTypes == "extravillous trophoblasts"] = "organ-spec. cells"
cellTypes[cellTypes == "syncytiotrophoblasts"] = "organ-spec. cells"
cellTypes[cellTypes == "hofbauer cells"] = "organ-spec. cells"
cellTypes[cellTypes == "basal prostatic cells"] = "organ-spec. cells"
cellTypes[cellTypes == "urothelial cells"] = "organ-spec. cells"
cellTypes[cellTypes == "prostatic glandular cells"] = "organ-spec. cells"
cellTypes[cellTypes == "paneth cells"] = "organ-spec. cells"
cellTypes[cellTypes == "intestinal goblet cells"] = "organ-spec. cells"
cellTypes[cellTypes == "langerhans cells"] = "organ-spec. cells"
cellTypes[cellTypes == "basal keratinocytes"] = "organ-spec. cells"
cellTypes[cellTypes == "suprabasal keratinocytes"] = "organ-spec. cells"
cellTypes[cellTypes == "proximal enterocytes"] = "organ-spec. cells"
cellTypes[cellTypes == "late spermatids"] = "organ-spec. cells"
cellTypes[cellTypes == "early spermatids"] = "organ-spec. cells"
cellTypes[cellTypes == "adipocytes"] = "organ-spec. cells"



unique(cellTypes)

#Create a cell type factor that is ordered differently
ct = rep(NA, length(cellTypes))
ct[cellTypes == "organ-spec. cells"] = 1
ct[cellTypes == "lymphocytes"] = 2
ct[cellTypes == "other immune cells"] = 3
ct[cellTypes == "endothelial cells"] = 4
ct[cellTypes == "fibroblasts"] = 5
ct[cellTypes == "skeletal myocytes"] = 6
ct[cellTypes == "cardiomyocytes"] = 7
ct[cellTypes == "smooth muscle cells"] = 8
ct[cellTypes == "inhibitory neurons"] = 9
ct[cellTypes == "excitatory neurons"] = 10

#ctNames = c("organ-spec. cells", "lymphocytes", "other immune cells", "endothelial cells", "fibroblasts", "skeletal myocytes", 
#            "cardiomyocytes", "smooth muscle cells", "inhibitory neurons", "excitatory neurons")
ctNames = c("organ-spec.", "lymphocytes", "other immune", "endothelial", "fibroblasts", "skel. myocytes", 
            "cardiomyocytes", "smooth muscle", "inh. neurons", "exc. neurons")

df = tibble(x=x, y=y, tissue = as.factor(tissues), celltype = factor(ct, 1:10, ctNames))
#df = tibble(x=x, y=y, tissue = as.factor(tissues), celltype = as.factor(cellTypes))

#color_palette = c('#000000','#888888','#FF0000','#00BB00','#DC556B','#6B97EC','#BC976B', '#0000FF')  # light green, light yellow, light blue

#labels = unlist(S_5$sampleIds)

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

cl1 = ggplotColours(n=19)
cl1[1] = "#666666"
cl2 = c("#999999", ggplotColours(n=length(unique(cellTypes))-1))



pA = ggplot(df, aes(x = x, y = y, color=tissue, shape=tissue)) +
  geom_point(size=2, stroke = 2) +
  scale_color_manual(values = cl1) +
  scale_shape_manual(values = c(1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1)) +
  labs(y=expression("t-SNE Y"), x="t-SNE X", title="") +
  theme_bw() + 
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  theme(legend.title = element_blank(),legend.position="bottom", legend.text=element_text(size=14)) +
  theme(text = element_text(size=14),
                 axis.text.x = element_text(color='black', size=14),
                 axis.text.y = element_text(color='black', size=14)) +
  guides(colour = guide_legend(nrow=7,override.aes = list(size=3)))

pA

pB = ggplot(df, aes(x = x, y = y, color=celltype, shape=celltype)) +
  geom_point(size=2, stroke = 2) +
  scale_color_manual(values = cl2) +
  scale_shape_manual(values = c(1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6)) +
  ggplot2::labs(y=expression("t-SNE Y"), x="t-SNE X", title="") +
  ggplot2::theme_bw() + 
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  ggplot2::theme(legend.title = element_blank(),legend.position="bottom", legend.text=element_text(size=14)) + #guides(colour = guide_legend(nrow = 4), size = guide_legend(nrow = 4), linetype = guide_legend(nrow = 4)) +
  ggplot2::theme(text = element_text(size=14),
                 axis.text.x = element_text(color='black', size=14),
                 axis.text.y = element_text(color='black', size=14)) +
  guides(colour = guide_legend(nrow=4,override.aes = list(size=3)))
pB



ggsave(
  paste0(figPath, "Fig5A.png"),
  plot = pA,
  width = 6, height = 8, dpi = 1200) 

ggsave(
  paste0(figPath, "Fig5A.eps"),
  plot = pA,
  width = 6, height = 8, dpi = 300) 

ggsave(
  paste0(figPath, "Fig5B.png"),
  plot = pB,
  width = 6, height = 7.3, dpi = 1200) 

ggsave(
  paste0(figPath, "Fig5B.eps"),
  plot = pB,
  width = 6, height = 7.3, dpi = 300) 


