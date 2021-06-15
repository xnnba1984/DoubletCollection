#' Read rds datasets
#'
#' Read all rds files under current working directory or a user-specified directory.
#' @param path A string of a user-specified directory.
#'
#' @return A list of data matrix and doublet annotations.
#' @export
#'
#' @examples
#' data.list <- ReadData(path = ".../real_datasets")
#'
ReadData <- function(path = NULL){
if(is.null(path)){
  path <- getwd()
}
print(path)
file.names <- list.files(path)
data.list <- list()
data.list$count <- list()
data.list$label <- list()
for(file in file.names){
  if(grepl('.rds', file, fixed = T)){
    cat('Read ', file, '\n', sep = '')
    name <- strsplit(file, ".", fixed = T)[[1]][1]
    data <- readRDS(paste(path, '/', file, sep = ''))
    data.list$count[[name]] <- data[[1]]
    data.list$label[[name]] <- data[[2]]
  }
}
return(data.list)
}


#' Transform a list to a dataframe
#'
#' Transform a list to a dataframe for visualization.
#' @param l A list of AUPRCs or AUROCs of different doublet-detection methods on multiple datasets.
#' @param type A charactor of "boxplot", "lineplot", "heatmap", "barplot", or "distributed" that refer to different format.
#'
#' @return A dataframe with approriate format for visualization.
#' @export
#'
#' @examples
#' result.auprc <- ListToDataframe(auprc.list.all, 'boxplot')
#' result.auroc <- ListToDataframe(auroc.list.all, 'boxplot')
#'
ListToDataframe <- function(l, type){
  if(type=='boxplot'){
    result <- data.frame(dataset=rep(names(l), each = length(l[[1]])),
                         method=rep(names(l[[1]]), length(l)),
                         value=unlist(l), row.names = NULL)
  }
  if(type=='lineplot'){
    result <- data.frame(parameter=rep(names(l), each = length(l[[1]])),
                         method=rep(names(l[[1]]), length(l)),
                         value=unlist(l), row.names = NULL)
  }
  if(type=='heatmap'){
    for(rate in names(l)){
      l1 <- l[[rate]]
      for(data in names(l1)){
        l2 <- l1[[data]]
        for(method in names(l2)){
          cluster.num <- nlevels(l2[[method]])
          l[[rate]][[data]][[method]] <- cluster.num
        }
      }
    }
    result <- data.frame(removal_rate=rep(names(l), each = length(unlist(l[[1]]))),
                         cluster_correct = rep(names(l[[1]]), each = length(unlist(l[[1]][1]))),
                         method = rep(names(l[[1]][[1]])),
                         Louvain_clustering=unlist(l), row.names = NULL)
  }
  if(type=='barplot'){
    result <- data.frame(method=rep(names(l), each = length(l[[1]])),
                         measurement=rep(names(l[[1]]), length(l)),
                         value=unlist(l), row.names = NULL); result
  }
  if(type=='distributed'){
    result <- data.frame(batch=rep(names(l), each = length(l[[1]])),
                         method=rep(names(l[[1]]), length(l)),
                         value=unlist(l), row.names = NULL); result
  }
  return(result)
}


#' Split data into different batches
#'
#' Split data into different batches.
#' @param count A scRNA-seq data matrix.
#' @param label A vector of doublet annotations.
#' @param batch A batch number.
#'
#' @return A list of randomly split datasets
#' @export
#'
#' @examples
#'
SplitData <- function(count, label, batch){
  count.list <- list()
  label.list <- list()
  index <- 1:dim(count)[2]
  index <- sample(x=index, size = length(index), replace = F)
  index.batch <- split(index, cut(seq_along(index), batch, labels = F))
  for(i in 1:length(index.batch)){
    #i=1
    count.batch <- count[,index.batch[[i]]]; dim(count.batch)
    count.list[[as.character(i)]] <- count.batch
    label.batch <- label[index.batch[[i]]]
    label.list[[as.character(i)]] <- label.batch
  }
  data.split <- list(count=count.list, label=label.list)
  return(data.split)
}










