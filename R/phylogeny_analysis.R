
PATH <- '/data1/bio/runs-kovarskiy/metagenome_data/'
setnum <- 2
nsamples <- dim(data)[1]
#data <- matrix(c("r","a","y","g","g","a","c","-","c","t","c","g","a","a","t","g","g","a","t","-","c","t","c","a","a","a","t","-","g","a","c","c","c","t","?","g"),
#                dimnames = list(c("t1", "t2", "t3"),NULL), nrow=3, byrow=TRUE)


contrast <- matrix(data = c(1,0,0,0,0,
                           0,1,0,0,0,
                           0,0,1,0,0,
                           0,0,0,1,0,
                           0,0,0,0,1,
                           0,0,0,0,0),
                          ncol = 5, byrow = TRUE)
dimnames(contrast) <- list(c("A","C","G","T","-","?"), c("a", "c", "g", "t", "-"))
phydata <- phyDat(data, type="USER", contrast=contrast)
#phydata <- phyDat(data, type="USER", levels=c("A","C","G","T","-"),ambiguity = c("?"))

nj_tree <- NJ(phydata)
dm1 <- dist.dna(phydata)
dm2 <- dist.hamming(phydata)
nj_tree <- NJ(dm2)
plot(nj_tree)
up_tree <- upgma(dm2)
plot(up_tree)

###Old_matrix
library(reshape)
library(phytools)
library(rjson)
library(ggplot2)
library(phangorn)
#При помощи новой функции, которая не выкидывает элементы с нулевыми расстояниями
#получаем матрицы расстояний
setnum <- 0
filename <- paste(PATH,setnum,'fastasamples.list',sep='/')
data <- read.table(file=filename,stringsAsFactors=F)$V1
old_data <- load_distances2(setnum,treshold=0.6)
old_data <- old_data[[3]][data,data]
#Получение таблицы (списка вершин и ребер)
temp <- old_data
temp[upper.tri(x=temp,diag=FALSE)] <- NA
nett <- melt(temp)
nett <- nett[which(!is.na(nett$x))]
#Строим дерево по матрице расстояний
up_tree <- upgma(old_data)
splits <- as.matrix(as.splits(up_tree))
edges <- up_tree$edge[order(up_tree$edge[,1]),]
#Получаем список пар разбиений (названия листьев)
res <- list()
#matchings <- as.matching(up_tree, labels=TRUE)
for (i in 1:up_tree$Nnode){
  print(i)
  temp <- list()
  temp[[1]] <- subset(colnames(splits),splits[edges[(i*2-1),2],]==1)
  temp[[2]] <- subset(colnames(splits),splits[edges[(i*2),2],]==1)
  #if (length(intersect(temp[[1]], temp[[2]])) != 0){
  #  temp[[2]] <- subset(colnames(splits),splits[up_tree$edge[i,2],]==0)
  #}
  print(intersect(temp[[1]], temp[[2]]))
  #temp[[1]] <- colnames(splits)[which(splits[up_tree$edge[i,1],]==1)]
  #temp[[2]] <- colnames(splits)[which(splits[up_tree$edge[i,2],]==1)]
  res[[i]] <- temp
}
#Сохраняем список в файл
sink(paste(PATH,setnum,"tree.json",sep='/'))
cat(toJSON(res))
sink()

#!!!test
temp <- list(c(1), list(1,2))
sink(paste(PATH,setnum,"temp.json",sep='/'))
cat(toJSON(temp))
sink()


#Фнукция грузит матрицу расстояний, не выкидывая нулевые
load_distances2 <- function(setnum, suffix='alt_4', treshold,all=TRUE,group='USA',drop_outlier=FALSE){
  cross_cov_name <- paste('cross_cov_',suffix,'.out', sep='')
  cross_cov <- load_matrix(setnum, cross_cov_name, header = T)
  cross_cov <- cross_cov + t(cross_cov)
  cross_diff_name <- paste('cross_diffs_',suffix,'.out', sep='')
  cross_diff <- load_matrix(setnum, cross_diff_name, header = T)
  cross_diff <- cross_diff + t(cross_diff)
  temp <- c(as.character(data.cov_distrs$names[which(data.cov_distrs$frac > treshold & data.cov_distrs$setnums == setnum)]),'REF')
  if (!all){
    group_names <- data.sample_names2$names[data.sample_names2$country==group]
    temp <- temp[temp%in%group_names]
  }
  #str(temp)
  good_samples <- which(colnames(cross_cov)%in%temp)
  #print("SRS014923"%in%temp)
  #print(good_samples)
  short_cross_cov <- cross_cov[good_samples,good_samples]
  #short_cross_cov <- drop_less(cross_cov, treshold)
  #short_cross_cov <- cross_cov
  short_cross_diff <- cross_diff[colnames(short_cross_cov),colnames(short_cross_cov)]
  distance <- short_cross_diff / short_cross_cov
  distance <- drop_NAN(distance)
  #distance <- drop_null(distance)
  if (drop_outlier){
    distance <- drop_outliers(distance)
  }
  
  #distance <- drop_less(distance)
  return(list(short_cross_cov, short_cross_diff,distance))
}
#tree <- rtree(nsamples)
#fitter <- pml(tree, phydata)