###Old_matrix
library(reshape)
library(phytools)
library(rjson)
library(ggplot2)
library(phangorn)
library(igraph)
PATH <- '/data1/bio/runs-kovarskiy/metagenome_data/'
setnum <- 40
frac <- 0.5
posnums <- 100
#Получаем все матрицы с фиксированной долей позиций:
distances <-  load_distances_bootstrapped(setnum,frac,posnums)
trees <- lapply(distances, upgma)
class(trees) <- "multiPhylo"
netx <- consensusNet(trees)
maxnum <- max(netx$edge)
leafnum <- length(netx$tip.label)
matches <- data.frame(old=c(1:maxnum), new=c(netx$tip.label, (leafnum+1):maxnum))
netframe <- as.data.frame(apply(netx$edge,2,function(x){matches$new[match(x,matches$old)]}))
write.table(netframe, file=paste('/data1/bio/runs-kovarskiy/metagenome_data/', setnum,'/network.tsv',sep=''), sep='\t',col.names=FALSE,quote=FALSE,row.names=FALSE)

filename <- paste(PATH,setnum,'fastasamples.list',sep='/')
data <- read.table(file=filename,stringsAsFactors=F)$V1
old_data1 <- load_distances2(setnum,treshold=0.6)
old_data1 <- old_data[[3]][data,data]
#Строим дерево по матрице расстояний
up_tree <- upgma(old_data)
splits <- as.matrix(as.splits(up_tree))
edges <- up_tree$edge[order(up_tree$edge[,1]),]
#Получаем список пар разбиений (названия листьев)
res <- list()




#Фнукция грузит матрицу расстояний, не выкидывая нулевые
load_distances_bootstrapped <- function(setnum, frac, times=100){
  PATH <- '/data1/bio/runs-kovarskiy/metagenome_data/'
  namefile <- paste(PATH,setnum,'newnames.list',sep='/')
  newnames <- read.table(file=namefile,stringsAsFactors=F)$V1
  #np.savetxt('%s/sampled_cross_cov_%.1f_%d.out' % (path, frac, time), cross_cov, delimiter='\t',fmt='%.1f')
  i <- 0
  distances <- list()
  while(i < times){
    cross_cov_name <- paste('sampled_cross_cov_',frac,'_',i,'.out', sep='')
    cross_cov <- load_matrix(setnum, cross_cov_name, header = F)
    cross_cov <- cross_cov + t(cross_cov)
    cross_diff_name <- paste('sampled_cross_diffs_',frac,'_',i,'.out', sep='')
    cross_diff <- load_matrix(setnum, cross_diff_name, header = F)
    cross_diff <- cross_diff + t(cross_diff)
    dimnames(cross_cov) <- list(newnames,newnames)
    dimnames(cross_diff) <- list(newnames,newnames)
    #temp <- c(as.character(data.cov_distrs$names[which(data.cov_distrs$frac > 0 & data.cov_distrs$setnums == setnum)]),'REF')
    #good_samples <- which(colnames(cross_cov)%in%temp)
    #short_cross_cov <- cross_cov[good_samples,good_samples]
    #short_cross_cov <- drop_less(cross_cov, treshold)
    #short_cross_cov <- cross_cov
    #short_cross_diff <- cross_diff[colnames(short_cross_cov),colnames(short_cross_cov)]
    distance <- cross_diff / cross_cov
    #distance <- drop_NAN(distance)
    #distance <- drop_null(distance)
    #distance <- drop_less(distance)
    i <- i + 1
    distances[[i]] <- distance
  }
  return(distances)
}


load_subsetted_matrix <- function(setnum, filename, names,header=F){
  path <- paste(PATH, setnum, '/', filename, sep='')
  namepath <- paste(PATH,setnum,'/names.list',sep='')
  data <- read.csv(path, sep='\t', header=F)
  namedata <- read.csv(namepath, sep='\t', header=F)
  head_names <- c(as.character(namedata[,1]), 'REF')
  if (header){
    #head_names <- lapply(data[1,],function(x){return(str_replace(x,'\\#?\\s',''))})
    #head_names[which(is.na(head_names))] <- 'ref'
    #data <- as.matrix(data[-1,])
    #data <- apply(data,2,as.numeric)
    colnames(data) <- head_names
    rownames(data) <- head_names
  }
  return(data)
}