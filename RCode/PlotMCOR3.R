library(ggplot2)
library(scales) # for muted function

figPath = "Z:/projects/Single-cell modeling/figures/"

setwd("C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling")

library(R.matlab)

##################################
# Structural comparison
##################################

#Set 5
#######################

S_5 = readMat("data/MCOR3/MCOR3_rxns_run_5.mat")$s.5[,,1] 
S_5
dim(S_5$data) # 10376 rxns

data5 = S_5$data
onOff5 = data5;
onOff5[,] = NA
onOff5[data5 <=1] = 0;
onOff5[data5 <=10 & data5 >1] = 1;
onOff5[data5 <=50 & data5 >10] = 2;
onOff5[data5 <=89 & data5 >50] = 3;
onOff5[data5 <=98 & data5 >89] = 4;
onOff5[data5 >=99] = 5;

sel = (rowSums(onOff5 == 0) > 0) & (rowSums(onOff5 == 5) > 0)
sum(sel) #387, maybe a bit too much to present
length(sel)
onOff5[sel,]

#test a tsne plot
x = as.numeric(S_5$tsneX)
y = as.numeric(S_5$tsneY)

res = prcomp(t(data5))
x = res$x[,1]
y = res$x[,3]
df = tibble(x=x, y=y, tissue = as.factor(unlist(S_5$sampleIds)))

color_palette = c('#000000','#888888','#FF0000','#00BB00','#DC556B','#6B97EC','#BC976B', '#0000FF')  # light green, light yellow, light blue

labels = unlist(S_5$sampleIds)

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

palette2 <- ggplotColours(n=17) #should be the same as in Seurat


pC = ggplot(df, aes(x = x, y = y, color=tissue, shape=tissue)) +
  geom_point(size=2, stroke = 2) +
  scale_color_manual(values = palette2, labels = labels) + #using the same colors as in Seurat
  scale_shape_manual(values =               c(1,1,1,1,1,2,2,3,3,3,1,1,1,2,4,4,5), labels=labels) +
  ggplot2::labs(y=expression("PC 3"), x="PC 1", title="Structural comparison") +
  ggplot2::theme_bw() + #+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  ggplot2::theme(legend.title = element_blank()) +
  ggplot2::theme(legend.title = element_blank(),legend.position="right", legend.text=element_text(size=14)) + #guides(colour = guide_legend(nrow = 4), size = guide_legend(nrow = 4), linetype = guide_legend(nrow = 4)) +
  ggplot2::theme(text = element_text(size=14),
                 axis.text.x = element_text(color='black', size=14),
                 axis.text.y = element_text(color='black', size=14))
pC

ggsave(
  paste0(figPath, "Fig3C.png"),
  plot = pC,
  width = 6, height = 5, dpi = 300)


#################################################################
# Supplementary to show that there is no clear grouping on layer
#################################################################

x = res$x[,1]
y = res$x[,2]
df = tibble(x=x, y=y, tissue = as.factor(unlist(S_5$sampleIds)))

pSup1A = ggplot(df, aes(x = x, y = y, color=tissue, shape=tissue)) +
  geom_point(size=2, stroke = 2) +
  scale_color_manual(values = palette2, labels = labels) + #using the same colors as in Seurat
  scale_shape_manual(values = c(1,rep(2,6),rep(3,7),rep(4,3)), labels=labels) +
  ggplot2::labs(y=expression("PC 2"), x="PC 1", title="") +
  ggplot2::theme_bw() + #+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  ggplot2::theme(legend.title = element_blank()) +
  ggplot2::theme(legend.title = element_blank(),legend.position="right", legend.text=element_text(size=14)) + #guides(colour = guide_legend(nrow = 4), size = guide_legend(nrow = 4), linetype = guide_legend(nrow = 4)) +
  ggplot2::theme(text = element_text(size=14),
                 axis.text.x = element_text(color='black', size=14),
                 axis.text.y = element_text(color='black', size=14))
pSup1A

x = res$x[,3]
y = res$x[,4]
df = tibble(x=x, y=y, tissue = as.factor(unlist(S_5$sampleIds)))

pSup1B = ggplot(df, aes(x = x, y = y, color=tissue, shape=tissue)) +
  geom_point(size=2, stroke = 2) +
  scale_color_manual(values = palette2, labels = labels) + #using the same colors as in Seurat
  scale_shape_manual(values = c(1,rep(2,6),rep(3,7),rep(4,3)), labels=labels) +
  ggplot2::labs(y=expression("PC 4"), x="PC 3", title="") +
  ggplot2::theme_bw() + #+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())
  ggplot2::theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  ggplot2::theme(legend.title = element_blank()) +
  ggplot2::theme(legend.title = element_blank(),legend.position="right", legend.text=element_text(size=14)) + #guides(colour = guide_legend(nrow = 4), size = guide_legend(nrow = 4), linetype = guide_legend(nrow = 4)) +
  ggplot2::theme(text = element_text(size=14),
                 axis.text.x = element_text(color='black', size=14),
                 axis.text.y = element_text(color='black', size=14))
pSup1B

figSup1 = ggarrange(pSup1A,pSup1B, nrow=1, ncol=2, labels=c("A","B"), font.label = list(size = 24))

ggsave(
  paste0(figPath, "FigSup3_1.png"),
  plot = figSup1,
  width = 12, height = 5, dpi = 300)

####################################
# Tasks comparison
####################################

#Set 5
#######################

T_5 = readMat("data/MCOR3/MCOR3_tasks_run_5.mat")$t.5[,,1] 
T_5
tasks5 = unlist(T_5$tasks)
cellTypes5 = unlist(T_5$sampleIds)
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
sum(taskSel)#13, good number to present


onOff5Filt = onOff5[taskSel,]

vals = rep(NA,nrow(onOff5Filt)*ncol(onOff5Filt))
for (i in 1:ncol(onOff5Filt)) {
  vals[(1:nrow(onOff5Filt))+(nrow(onOff5Filt)*(i-1))] = onOff5Filt[,i];
}

shortenedTasks5 = str_replace(tasks5, fixed(" (physiological substrates, physiological excretion)"),"")
shortenedTasks5 = str_replace(shortenedTasks5, fixed(" (minimal substrates, physiological excretion)"),"")
shortenedTasks5 = str_replace(shortenedTasks5, fixed(" (minimal substrates, minimalexcretion)"),"")
shortenedTasks5 = str_replace(shortenedTasks5, fixed(" (minimal substrates with AA, minimal excretion)"),"")
shortenedTasks5 = str_replace(shortenedTasks5, fixed(" de novo synthesis"),"")

valsFactor = factor(vals, c(0,1,2,3,4,5), c("< 2%", "2-10%", "11-50%", "51-89%", "90-98%", "> 98%"))
samples = factor(rep(1:17,1,each=nrow(onOff5Filt)), 1:17, cellTypes5)
tasks = factor(rep(1:nrow(onOff5Filt),ncol(onOff5Filt)), 1:nrow(onOff5Filt), shortenedTasks5[taskSel])

df = tibble(vals = valsFactor,samples,tasks)

p3D = ggplot(df, aes(samples, tasks)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = vals), color="dark gray") + # background colours are mapped according to the value column
  scale_fill_manual(values=c("white", '#cccccc','#aaaaaa','#777777','#222222', "black")) +
  ggtitle("Task analysis") +
  theme(legend.title=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        text = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=90, hjust = 1,vjust=0.5,size = 14, color="black"),
        axis.text.y = element_text(size = 14, color="black"))
  
p3D



ggsave(
  paste0(figPath, "Fig3D.png"),
  plot = p3D,
  width = 7, height = 5, dpi = 300)


