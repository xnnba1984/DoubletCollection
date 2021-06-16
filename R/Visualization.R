#' Boxplot
#'
#' Plot a boxplot that shows the measurement distributions for different methods on multiple datasets.
#' Methods are sorted based on average performance.
#' @param result A dataframe with approriate format for visualization.
#' @param measurement A character to show the title of boxplot.
#' @param save Whether saving the plot to a file.
#' @param path The local directory to save the plot.
#'
#' @return No return value.
#' @export
#' @import ggplot2
#'
#' @examples
#' Plot_Boxplot(result = result.auprc, measurement = 'AUPRC', save=T, name = 'AUPRC.png', path = getwd())
#' Plot_Boxplot(result = result.auroc, measurement = 'AUROC', save=T, name = 'AUROC.png', path = getwd())
#'
Plot_Boxplot<- function(result, measurement, save=T, path=getwd(), name){

  result.ave <- aggregate(result$value, list(result$method), mean, na.rm = TRUE)
  colnames(result.ave) <- c('method','ave_value')
  result.ave <- result.ave[order(result.ave$ave_value),]
  result$method <- factor(result$method, levels = result.ave$method)

  print(
  ggplot(result, aes(x=method, y=value, fill=method)) + geom_boxplot() + theme_bw() +
    labs(x=NULL, y=NULL, title=measurement) + theme(text = element_text(size=15)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=22),
          axis.text.y = element_text(size=22),
          plot.title = element_text(hjust = 0.5,size=25))+
    theme(axis.ticks=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=3))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position = "none")
  )
  if(save){
    cat("Save plots to ", path, '\n')
    ggsave(name, path = path)
  }
}


#' Lineplot
#'
#' Plot a lineplot across different parameters of scRNA-seq data for multiple doublet-detection methods.
#' @param result A dataframe with approriate format for visualization.
#' @param type A character to show the parameter of scRNA-seq data.
#' @param measurement A character to show the measurement of lineplot.
#'
#' @return No return values.
#' @export
#' @import ggplot2
#'
#' @examples
#'
#' result.auprc <- ListToDataframe(auprc.list.all, type = 'lineplot')
#' result.auroc <- ListToDataframe(auroc.list.all, type = 'lineplot')
#' Plot_Lineplot(result = result.auprc, type = 'Heterogeneity', measurement = 'AUPRC')
#' Plot_Lineplot(result = result.auroc, type = 'Doublet Rate', measurement = 'AUROC')
#'
Plot_Lineplot<- function(result, type, measurement, save=T, path=getwd(), name){

  print(
  ggplot(result, aes(x=as.numeric(parameter), y=value, col=method, group=method)) + geom_line(size=1) + geom_point(size=2) + theme_bw()+
    labs(x=NULL, y=measurement, title=type) +
    theme(plot.title = element_text(hjust = 0.5, size=25), axis.text.x = element_text(size=22),
          axis.text.y = element_text(size=22), text=element_text(size=20),
          panel.border = element_rect(colour = "black", fill=NA, size=3)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank())
  )
  if(save){
    cat("Save plots to ", path, '\n')
    ggsave(name, path = path)
  }
}


#' Barplot
#'
#' Plot Barplot for visualizing DE gene analysis.
#' @param result A data frame with appropriate format for visualization.
#' @param measurement A title of visualization.
#'
#' @return No return values.
#' @export
#' @import ggplot2
#'
#' @examples
#' Plot_Barplot(table.DE.all[table.DE.all$measurement=='precision',], 'Precision')
#' Plot_Barplot(table.DE.all[table.DE.all$measurement=='recall',], 'Recall')
#' Plot_Barplot(table.DE.all[table.DE.all$measurement=='tnr',], 'TNR')
#'
Plot_Barplot <- function(result, measurement, save=T, path=getwd(), name){

  print(
  ggplot(result, aes(fill=method, y=value, x = DE_method)) +
    geom_bar(position=position_dodge(width=0.8), stat="identity", width = .7) + theme_bw() +
    labs(x=NULL, y=NULL, fill='', title=measurement) +
    theme(plot.title = element_text(hjust = 0.5, size=25),
          axis.text.y = element_text(size=22),
          axis.text.x = element_text(size=22),
          panel.border = element_rect(colour = "black", fill=NA, size=3)) +
    coord_flip(ylim=c(min(result$value*.995),max(result$value*1.005))) +
    theme(text = element_text(size=20)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  )
  if(save){
    cat("Save plots to ", path, '\n')
    ggsave(name, path = path)
  }
}


#' The temporally DE genes version of Barplot
#'
#' Plot barplot for visualizing temporally DE genes analysis.
#' @param result A data frame with approriate format for visualization.
#' @param title A title of the visualization.
#'
#' @return No return values.
#' @export
#' @import ggplot2
#'
#' @examples
#' table.DE.temp <- ListToDataframe(l = de.temp.result.all, type = 'barplot')
#' Plot_Barplot_temp(result = table.DE.temp, title = 'Temporally DE Genes')
#'
Plot_Barplot_temp <- function(result, title, save=T, path=getwd(), name){
  result$measurement <- factor(result$measurement, levels = rev(c('precision','recall','tnr')))
  print(
  ggplot(result, aes(fill=method, y=value, x = measurement)) +
    geom_bar(position=position_dodge(width=0.8), stat="identity", width = .7) + theme_bw() +
    labs(x=NULL, y=NULL, fill='', title=title) +
    theme(plot.title = element_text(hjust = 0.5, size=25),
          axis.text.y = element_text(size=22),
          axis.text.x = element_text(size=22),
          panel.border = element_rect(colour = "black", fill=NA, size=3)) +
    coord_flip(ylim=c(min(result$value*.995),max(result$value*1.005))) +
    theme(text = element_text(size=20)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  )
  if(save){
    cat("Save plots to ", path, '\n')
    ggsave(name, path = path)
  }
}

#' Heatmap
#'
#' Plot Heatmap for visualizing cell clustering result.
#' @param result A dataframe with approriate format for visualization.
#' @param cluster A number of true clusters in the dataset (numeric).
#'
#' @return No return values.
#' @export
#' @import ggplot2
#' @import ggthemes
#'
#' @examples
#' table.cluster <- ListToDataframe(l = result.cluster.all.rate, type = 'heatmap'); table.cluster
#' Plot_Heatmap(result = table.cluster, cluster = 4)
#' Plot_Heatmap(result = table.cluster, cluster = 6)
#' Plot_Heatmap(result = table.cluster, cluster = 8)
#'
Plot_Heatmap<- function(result, cluster, save=T, path=getwd(), name){
  result <- result[result$cluster_correct==cluster,]
  result$correct <- result$Louvain_clustering==result$cluster_correct
  result.correct <- aggregate(result$correct, list(result$method), mean, na.rm = TRUE)
  colnames(result.correct) <- c('method', 'rate_correct')
  result.correct <- result.correct[order(result.correct$rate_correct),]
  result$method <- factor(result$method, levels = result.correct$method)
  result$Louvain_clustering <- ifelse(result$Louvain_clustering != cluster, cluster+1, result$Louvain_clustering)
  result$removal_rate <- as.numeric(result$removal_rate)

  title <- paste(as.character(cluster), 'Cell Types', sep = ' ')

  print(
  ggplot(result, aes(x=method, y=removal_rate, fill=Louvain_clustering)) + geom_tile(color="white", size=1) +
    labs(x=NULL, y='% of Removed Droplets \n(Detected Doublets)', title= title) +
    theme_tufte(base_family="Helvetica") +
    theme(axis.ticks=element_blank()) +
    theme(axis.text.x = element_text(margin = margin(l = 0), angle = 45, hjust = .9,size=22),
          plot.title = element_text(hjust = 0.5, size=25,vjust = -0.1),
          axis.text.y = element_text(size=22),
          text=element_text(size=20),
          legend.title = element_blank()) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       breaks = seq(from = min(result$removal_rate), to = max(result$removal_rate), by = 0.04)) +
    scale_fill_gradient(low="#fee391", high="#cc4c02",limits = c(cluster, cluster+1),
                        breaks = seq(from=cluster, to=cluster+1, by = 1), guide = "legend",
                        labels = c("Correct", "Incorrect"))
  )
  if(save){
    cat("Save plots to ", path, '\n')
    ggsave(name, path = path)
  }
}


#' Lineplots for distributed computing
#'
#' Plot lineplots of detection accuracy across different batch numbers and methods.
#' @param result A data frame with approriate format to visualize.
#' @param data Name of the dataset to visualize.
#' @param measurement "AUPRC" or "AUROC".
#'
#' @return No return values.
#' @export
#' @import ggplot2
#'
#' @examples
#' Plot_Lineplot_Distributed(result = table.batch, data = 'pbmc-ch', measurement = 'AUPRC')
#'
Plot_Lineplot_Distributed<- function(result, data, measurement, save=T, path=getwd(), name){

  print(
  ggplot(result, aes(x=as.numeric(batch), y=value, col=method, group=method)) + geom_line(size=1) + geom_point(size=2) + theme_bw() +
    labs(x='Batch Number', y=measurement, title=data) +
    theme(plot.title = element_text(hjust = 0.5, size=25), axis.text.x = element_text(size=22),
          axis.text.y = element_text(size=22), text=element_text(size=20),
          panel.border = element_rect(colour = "black", fill=NA, size=3)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank())
  )
  if(save){
    cat("Save plots to ", path, '\n')
    ggsave(name, path = path)
  }
}

