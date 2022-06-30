library(ggplot2)

figPath = "Z:/projects/Single-cell modeling/figures/"

setwd("C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling")

  ############################################
  # Fig 1B - Gene essentiality
  ############################################

  # load data
  fns = c("data/geneEss_old_tINIT.txt",
          "data/geneEss_newalg.txt")
  
  names = c("tINIT","ftINIT")
  
  gea_res = NULL
  for (i in 1:length(fns)) {
    x = read.delim(file = fns[i], sep='\t', stringsAsFactors=F)
    x$model = names[i]
    gea_res = rbind(gea_res, x)
  }
  gea_res$model = factor(gea_res$model, as.character(names)[1:2])  # to enforce the model order

  
  color_palette <- c('#B5D39B','#6B97BC')  # light green, light blue

  p1B = ggplot(gea_res, aes(x = model, y = MCC, fill = model)) +
    geom_violin(trim=F, show.legend=F, scale='count') +
    scale_fill_manual(values=color_palette) +
    theme_classic() + 
    ylab('MCC') +
    xlab('') +
    theme(text = element_text(size=14),
          axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,
                                     color='black', size=14),
          axis.text.y = element_text(color='black', size=14),
          axis.line.x = element_blank()) +
    ylim(c(0.1,0.38)) # + #manipulate these numbers to include all data
    #ylim(c(0,0.5)) # +
  p1B
  

  ggsave(
    paste0(figPath, "FigGeneEss.png"),
    plot = p1B,
    width = 2.5, height = 3, dpi = 300)
  
  ggsave(
    paste0(figPath, "FigGeneEss.eps"),
    plot = p1B,
    width = 2.5, height = 3, dpi = 300)
  
  
########################################
# Fig. 1C - Execution times
########################################
  
  library(R.matlab)
  
  tINITTimes = as.numeric(readMat("data/old_tinit_exec_times.mat")$old.tinit.exec.times)
  tINIT2Times = as.numeric(readMat("data/new_tinit_exec_times.mat")$new.tinit.exec.times)
  
  names = c("tINIT","ftINIT")
  
  # include the next two lines to compare biomass vs. all-task methods
  # fnames <- c(fnames, paste('BMonly_results_', modelnames, '_DepMap_06.txt', sep=''))
  # modelnames <- c(modelnames, modelnames)
  
  times = c(tINITTimes,tINIT2Times)
  method = factor(c(rep(0,length(tINITTimes)), rep(1,length(tINIT2Times))),c(0,1), names)
  
  df = tibble(times = times, method=method)
  

  color_palette = c('#B5D39B','#6B97BC')  # light green, light blue
  
  p1C = ggplot(df, aes(x = method, y = times, fill = method)) +
    geom_boxplot() +
    scale_fill_manual(values=color_palette) +
    #scale_y_continuous(trans = 'log10', breaks = c(seq(10,90,by=10),seq(100,900,by=100),seq(1000,3000,by=1000)), labels = c("10",rep("",8),"100", rep("",8), "1000", rep("",2)), limits = c(10,3500)) +
    scale_y_continuous(trans = 'log10', limits = c(10,3500)) +
    theme_classic() + 
    ylab('Execution time (s)') +
    xlab('') +
    theme(text = element_text(size=14),
          axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,
                                     color='black', size=14),
          axis.text.y = element_text(color='black', size=14),
          axis.line.x = element_blank(), legend.title = element_blank(), legend.position = "none") +
    geom_hline(yintercept=0, color='black')
  p1C #the warning here is ok, it is there because not all combinations of x and y exist
  
  
  ggsave(
    paste0(figPath, "FigExecTimes.png"),
    plot = p1C,
    width = 2.5, height = 3, dpi = 300)
  
  ggsave(
    paste0(figPath, "FigExecTimes.eps"),
    plot = p1C,
    width = 2.5, height = 3, dpi = 300)
  
  
    