#PATH <- '/data1/bio/runs-kovarskiy/onkokids/snpdata/'
PATH <- '/data1/bio/runs-kovarskiy/onkokids_fewbacs/'
load_cov_distrs <- function(setnum){
  path <- paste(PATH,setnum,sep='')
  files <- list.files(path, setnum, pattern = '\\.hist',full.names=T)
  res <- data.frame(cover=c(), count=c(), names=c())
  for (file in files){
    fname <- substr(basename(file),1,nchar(basename(file))-5)
    temp <- try(read.csv(file,sep='\t',header=F))
    #temp <- read.csv(file,sep='\t',header=F)
    if (class(temp) != "try-error"){
      names(temp) <- c("cover", "count")
      summ <- sum(temp$count)
      total <- data.onko.strains$length[data.onko.strains$set == setnum]
      delta <- total - summ
      add_info <- data.frame(cover=c(0),count=c(delta))
      temp <- rbind(temp, add_info)
      temp$names <- fname
      res <- rbind(res,temp)
    }
    else{
      total <- data.onko.strains$length[data.onko.strains$set == setnum]
      temp <- data.frame(cover=c(0),count=c(total))
      temp$names <- fname
      res <- rbind(res,temp)
    }
  }
  return(res)
}

cov_distrs_proc <- function(setnums){
  i <- 1
  for (setnum in setnums){
    fname <- paste(PATH, setnum, 'pileups.short',sep='/')
    print(fname)
    namelist <- as.character(read.csv(fname,sep='\t',header=F)[,1])
    data <- load_cov_distrs(setnum)
    data$nts <- as.numeric(data$cover * data$count)
    #data <- data[which(data$cover!=0),]
    data$count[which(data$cover==0)] <- 0
    res <- aggregate(data[,c("count", "nts")],by=list(data$names),sum)
    length <- data.onko.strains$length[data.onko.strains$set == setnum]
    res$redund <- res$nts / length
    res$uncov <- length - res$count
    res$length <- length
    res$setnum <- setnum
    res$frac <- res$count/length
    res$mean <- res$nts / res$count
    names(res) <- c("names", "cover", "nts", "redund", "uncov", "length", "setnums","frac", "mean")
    res2 <- res
    res2$names <- factor(res2$names,levels=namelist)
    res2 <- res2[order(res2$names),]
    res2$mean[which(is.nan(res2$mean))] <- 0.0
    outname <- paste(PATH,setnum,'means.list',sep='/')
    write.table(res2$mean, outname, sep='\t',row.names=FALSE,quote=FALSE,col.names=FALSE)
    if (i == 1){
      rres <- res
    }
    else{
      rres <- rbind(rres, res)
    }
    i <- i + 1
  }
  return(rres)
}

load_matrix <- function(setnum, filename, header=F){
  path <- paste(PATH, setnum, '/', filename, sep='')
  namepath <- paste(PATH,setnum,'/pileups.short',sep='')
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

load_distances <- function(setnum, suffix='alt_2', treshold){
  cross_cov_name <- paste('cross_cov_',suffix,'.out', sep='')
  cross_cov <- load_matrix(setnum, cross_cov_name, header = T)
  cross_cov <- cross_cov + t(cross_cov)
  cross_diff_name <- paste('cross_diffs_',suffix,'.out', sep='')
  cross_diff <- load_matrix(setnum, cross_diff_name, header = T)
  cross_diff <- cross_diff + t(cross_diff)
  temp <- c(as.character(data.onko.coverages$names[which(data.onko.coverages$frac > treshold & data.onko.coverages$setnums == setnum)]),'REF')
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
  #distance <- drop_outliers(distance)
  #distance <- drop_less(distance)
  return(list(short_cross_cov, short_cross_diff,distance))
}

get_mean_dists <- function(setnums, suffix='alt_2', treshold){
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
  res <- merge(mdsdata, data.onko.meta, by.x=c("names"), by.y=c("names"), all.x=T,all.y=F)
  res$patient <- as.character(res$patient)
  res$patient[which(is.na(res$patient))] <- 'non_patient'
  res$patient[which(res$names == 'REF')] <- 'REF'
  res$patient <- as.factor(res$patient)
#   clusters <- pamk(data,krange=1:4,diss=TRUE)
#   avg.widths <- clusters$pamobject$silinfo$clus.avg.widths
#   nc <- clusters$nc
#   print(nc)
#   print(avg.widths)
#   clusters <- as.data.frame(clusters$pamobject$clustering)
#   names(clusters) <- c('clusters')
#   clusters$names <- row.names(clusters)
#   if (nc > 1){
#     if (all(avg.widths < 0.2)){
#       clusters$clusters <- 1
#     }
#   }
#   clusters$clusters <- as.factor(clusters$clusters)
#   res <- merge(res,clusters,by.x=c('names'),by.y=c('names'),all.x=T,all.y=T)
  return(res)
}

load_all_sets <- function(setnums, flag='alt_2',treshold){
  len <- length(setnums)
  res <- 0
  for (i in 1:len){
    print('==============================')
    print(setnums[i])
    data <- load_distances(setnums[i],flag,treshold)
    if (dim(data[[3]])[1] >= 4){
      temp <- make_data_frame_mds(data[[3]])
      temp$setnum <- setnums[i]
      if (class(res)=="numeric"){
        res <- temp
      }
      else{
        temp <- try(rbind(res,temp))
        if (class(temp) != "try-error"){
          res <- temp
        }
      }
    }
    else{
      print('!')
    }
  }
  return(res)
}
