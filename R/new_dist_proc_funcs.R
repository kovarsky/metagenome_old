PATH <- '/data1/bio/runs-kovarskiy/metagenome_data/'
#PATH <- '/data1/bio/runs-kovarskiy/onkokids/snpdata/'
library (stringr)
library (ggplot2)
library(reshape)
library(gplots)
library(RColorBrewer)
library(MASS)
library(fpc)

cols.gentleman <- function(ncol=500) {
  library(RColorBrewer)
  hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(ncol)
  return(rev(hmcol))
}


temp_preproc_abund <- function(setnum){
  temp <- mds.data.final_3[which(mds.data.final_3$setnum==setnum),]
  temp$country <- as.character(temp$country)
  temp <- temp[which(temp$country %in% c('CHN','EUR','USA')),]
  temp$country <- as.factor(temp$country)
  strainname <- as.character(strain_meta.data$short_descript[which(strain_meta.data$set==setnum)])
  temp <- temp[,c("short_names","country","clusters")]
  temp$clusters <- as.integer(temp$clusters)
  names(temp) <- c("names","country",strainname)
  return(temp)
}


load_matrix <- function(setnum, filename, header=F){
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

load_distances <- function(setnum, suffix='final_2', treshold,all=TRUE,group='USA',drop_outlier=TRUE){
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
  #print(good_samples)
  short_cross_cov <- cross_cov[good_samples,good_samples]
  #short_cross_cov <- drop_less(cross_cov, treshold)
  #short_cross_cov <- cross_cov
  short_cross_diff <- cross_diff[colnames(short_cross_cov),colnames(short_cross_cov)]
  distance <- short_cross_diff / short_cross_cov
  distance <- drop_NAN(distance)
  distance <- drop_null(distance)
  if (drop_outlier){
    distance <- drop_outliers(distance)
  }
  
  #distance <- drop_less(distance)
  return(list(short_cross_cov, short_cross_diff,distance))
}

get_mean_dists <- function(setnums, suffix='final_2', treshold){
  res <- data.frame(setnum=c(), mean=c(), median=c(), sd=c())
  for (set in setnums){
    temp <- load_distances(set, suffix, treshold)
    v <- as.vector(temp[[3]])
    data <- data.frame(setnum=set, mean=mean(v), median=median(v), sd=sd(v))
    res <- rbind(res,data)
  }
  return(res)
}

make_data_frame_mds <- function(data){
  mdsdata <- isoMDS(data)
  mdsdata <- as.data.frame(mdsdata$points)
  mdsdata$names <- rownames(mdsdata)
  res <- merge(mdsdata, data.sample_names2, by.x=c("names"), by.y=c("names"), all.x=T,all.y=F)
  res$country <- as.character(res$country)
  res$country[which(is.na(res$country))] <- 'RUS'
  res$country[which(res$names == 'REF')] <- 'REF'
  res$country <- as.factor(res$country)
  clusters <- pamk(data,krange=1:4,diss=TRUE)
  avg.widths <- clusters$pamobject$silinfo$clus.avg.widths
  nc <- clusters$nc
  print(nc)
  print(avg.widths)
  clusters <- as.data.frame(clusters$pamobject$clustering)
  names(clusters) <- c('clusters')
  clusters$names <- row.names(clusters)
  if (nc > 1){
    if (all(avg.widths < 0.2)){
      clusters$clusters <- 1
    }
  }
  clusters$clusters <- as.factor(clusters$clusters)
  res <- merge(res,clusters,by.x=c('names'),by.y=c('names'),all.x=T,all.y=T)
  return(res)
}

load_all_sets <- function(setnums, flag='final_2',treshold,all=TRUE,group='USA',drop_outlier=TRUE){
  len <- length(setnums)
  for (i in 1:len){
    print('==============================')
    print(setnums[i])
    data <- load_distances(setnums[i],flag,treshold,all,group,drop_outlier)
    if (dim(data[[3]])[1] >= 10){
      temp <- make_data_frame_mds(data[[3]])
      temp$setnum <- setnums[i]
      if (i==1){
        res <- temp
      }
      else{
        res <- rbind(res,temp)
      }
    }
  }
  return(res)
}

get_stat_info <- function(dlist){
  dlist <- lapply(dlist,function(x){return(x[1:322,1:322])})
  res <- data.frame(cross_cov=as.vector(dlist[[1]]),cross_diff=as.vector(dlist[[2]]),dist=as.vector(dlist[[3]]))
  return(res)
}

########################################################################
## Dropping Functions

drop_NAN <- function(dist_m){
  which_NA <- function(vector){
    return(length(which(is.na(vector))))
  }
  dist_m <- as.matrix(as.dist(dist_m))
  while(length(which(is.na(dist_m)))!=0){
    na_count_vector <- apply(dist_m,1,which_NA)
    max_nums <- which(na_count_vector==max(na_count_vector))
    del_num <- sample(length(max_nums),1)
    del_ind <- max_nums[del_num]
    dist_m <- dist_m[-del_ind,-del_ind]
  }
  return(dist_m)
}

drop_less <- function(dist_m, treshold=10^4){
  which_null <- function(vector){
    return(length(which(vector < treshold)))
  }  
  while((length(which(dist_m < treshold))>dim(dist_m)[1])&(dim(dist_m)[1] >= 10)){
    ###   Why len? Not which_null?'''
    count <- apply(dist_m,2,which_null)
    max_ <- max(count[count>1])
    max_indices <- which(count==max_)
    ind <- max_indices[sample(length(max_indices),1,1)]
    dist_m <- dist_m[-ind,-ind]
  }
  return(dist_m)
}


drop_outliers <- function(dist_m){
  print(dim(dist_m))
  mean_val <- mean(dist_m)
  print(mean_val)
  dev <- sd(as.vector(dist_m))
  tresh <- mean_val + 2 * dev
  which_more <- function(vector){
    return(length(which(vector>tresh)))
  }
  out_count_vector <- apply(dist_m,1,which_more)
  len <- length(out_count_vector)
  frac = 0.7
  len_tresh <- len * frac
  #print(len_tresh)
  del_ind <- which(out_count_vector > len_tresh)
  #print(del_ind)
  if (length(del_ind)>0){
    dist_m <- dist_m[-del_ind,-del_ind]
  }
  
  return(dist_m)
}

drop_null <- function(dist_m){
  which_null <- function(vector){
    return(length(which(vector == 0)))
  } 
  if (!all(dist_m==0)){
    while(length(which(dist_m==0))>dim(dist_m)[1]){
      ###   Why len? Not which_null?'''
      count <- apply(dist_m,2,which_null)
      max_ <- max(count[count>1])
      max_indices <- which(count==max_)
      ind <- max_indices[sample(length(max_indices),1,1)]
      dist_m <- dist_m[-ind,-ind]
    }
  }
  return(dist_m)
}

