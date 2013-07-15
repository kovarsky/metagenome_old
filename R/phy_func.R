phy_func <- function(setnum){
  library(reshape)
  library(phytools)
  library(rjson)
  library(ggplot2)
  library(phangorn)
  #При помощи новой функции, которая не выкидывает элементы с нулевыми расстояниями
  #получаем матрицы расстояний
  cat("-----set: ",setnum,'------\n')
  PATH <- '/data1/bio/runs-kovarskiy/metagenome_data/'
  filename <- paste(PATH,setnum,'fastasamples.list',sep='/')
  data <- read.table(file=filename,stringsAsFactors=F)$V1
  old_data <- load_distances2(setnum,treshold=0.6)
  old_data <- old_data[[3]][data,data]
  print(dim(old_data))
  if (!is.null(dim(old_data))){
    if (dim(old_data)[1] >= 15){
      #Строим дерево по матрице расстояний
      up_tree <- upgma(old_data)
      splits <- as.matrix(as.splits(up_tree))
      edges <- up_tree$edge[order(up_tree$edge[,1]),]
      #Получаем список пар разбиений (названия листьев)
      res <- list()
      #matchings <- as.matching(up_tree, labels=TRUE)
      for (i in 1:up_tree$Nnode){
        temp <- list()
        temp[[1]] <- subset(colnames(splits),splits[edges[(i*2-1),2],]==1)
        temp[[2]] <- subset(colnames(splits),splits[edges[(i*2),2],]==1)
        #if (length(intersect(temp[[1]], temp[[2]])) != 0){
        #  temp[[2]] <- subset(colnames(splits),splits[up_tree$edge[i,2],]==0)
        #}
        #temp[[1]] <- colnames(splits)[which(splits[up_tree$edge[i,1],]==1)]
        #temp[[2]] <- colnames(splits)[which(splits[up_tree$edge[i,2],]==1)]
        res[[i]] <- temp
      }
      #Сохраняем список в файл
      sink(paste(PATH,setnum,"newtree.json",sep='/'))
      cat(toJSON(res))
      sink()
    }
  }
}