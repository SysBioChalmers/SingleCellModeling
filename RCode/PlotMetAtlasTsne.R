library(ggplot2)
library(scales) # for muted function

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
  paste0(fig____path, "Fig5A.png"),
  plot = pA,
  width = 6, height = 8, dpi = 300) 

ggsave(
  paste0(fig____path, "Fig5A.svg"),
  plot = pA,
  width = 6, height = 8, dpi = 300) 

ggsave(
  paste0(fig____path, "Fig5B.png"),
  plot = pB,
  width = 6, height = 7.3, dpi = 300) 

ggsave(
  paste0(fig____path, "Fig5B.svg"),
  plot = pB,
  width = 6, height = 7.3, dpi = 300) 





####################################
# Tasks comparison
####################################

#Set 5
#######################

T_5 = readMat("data/LC3/LC3_tasks_run_5.mat")$t.5[,,1] 
T_5
tasks5 = unlist(T_5$tasks)
cellTypes5 = unlist(T_5$sampleIds)

cellTypes5 = str_replace(cellTypes5, fixed("N_"),"N:")
cellTypes5 = str_replace(cellTypes5, fixed("T_"),"T:")

data5 = T_5$data
onOff5 = data5;
onOff5[,] = NA
onOff5[data5 <=1] = 0;
onOff5[data5 <=10 & data5 >1] = 1;
onOff5[data5 <=50 & data5 >10] = 2;
onOff5[data5 <=89 & data5 >50] = 3;
onOff5[data5 <=98 & data5 >89] = 4;
onOff5[data5 >=99] = 5;

#Filter the tasks to only show tasks that have a clear on or off

taskSel = (rowSums(onOff5 == 0) > 0) & (rowSums(onOff5 == 5) > 0)
sum(taskSel)#14, good number to present


#filt = onOff1[,1] != onOff1[,2]
onOff5Filt = onOff5[taskSel,]

vals = rep(NA,nrow(onOff5Filt)*ncol(onOff5Filt))
for (i in 1:ncol(onOff5Filt)) {
  vals[(1:nrow(onOff5Filt))+(nrow(onOff5Filt)*(i-1))] = onOff5Filt[,i];
}

shortenedTasks5 = str_replace(tasks5, fixed(" (physiological substrates, physiological excretion)"),"")
shortenedTasks5 = str_replace(shortenedTasks5, fixed(" (minimal substrates, physiological excretion)"),"")
shortenedTasks5 = str_replace(shortenedTasks5, fixed(" (minimal substrates, minimalexcretion)"),"")
shortenedTasks5 = str_replace(shortenedTasks5, fixed(" (minimal substrates with AA, minimal excretion)"),"")
shortenedTasks5 = str_replace(shortenedTasks5, fixed(" de novo"),"")

valsFactor = factor(vals, c(0,1,2,3,4,5), c("< 2%", "2-10%", "11-50%", "51-89%", "90-98%", "> 98%"))
samples = factor(rep(1:16,1,each=nrow(onOff5Filt)), 1:length(cellTypes5), cellTypes5)
tasks = factor(rep(1:nrow(onOff5Filt),ncol(onOff5Filt)), 1:nrow(onOff5Filt), shortenedTasks5[taskSel])

df = tibble(vals = valsFactor,samples,tasks)

p4D = ggplot(df, aes(samples, tasks)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = vals), color="dark gray") + # background colours are mapped according to the value column
  scale_fill_manual(values=c("white", '#cccccc','#aaaaaa','#777777','#444444', "black")) +
  ggtitle("Task analysis") +
  theme(legend.title=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        text = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=90, hjust = 1,vjust=0.5,size = 14, color="black"),
        axis.text.y = element_text(size = 14, color="black"))
  
p4D


ggsave(
  paste0(fig____path, "Fig4D.png"),
  plot = p4D,
  width = 8.7, height = 5.4, dpi = 300)

ggsave(
  paste0(fig____path, "Fig4D.svg"),
  plot = p4D,
  width = 8.7, height = 5.4, dpi = 300)


