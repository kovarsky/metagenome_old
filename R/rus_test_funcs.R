PATH <- '/data1/bio/runs-kovarskiy/metagenome_data/rus_test/'
load_cov_distrs <- function(){
  path <- PATH
  files <- list.files(path, setnum, pattern = '\\.hist',full.names=T)
  res <- data.frame(cover=c(), count=c(), names=c())
  for (file in files){
    fname <- substr(basename(file),1,nchar(basename(file))-5)
    temp <- try(read.csv(file,sep='\t',header=F))
    #temp <- read.csv(file,sep='\t',header=F)
    if (class(temp) != "try-error"){
      names(temp) <- c("cover", "count")
      summ <- sum(temp$count)
      total <- data.rus.strains$length[1]
      delta <- total - summ
      add_info <- data.frame(cover=c(0),count=c(delta))
      temp <- rbind(temp, add_info)
      temp$names <- fname
      res <- rbind(res,temp)
    }
    else{
      total <- data.rus.strains$length[1]
      temp <- data.frame(cover=c(0),count=c(total))
      temp$names <- fname
      res <- rbind(res,temp)
    }
  }
  return(res)
}

cov_distrs_proc <- function(){
    fname <- paste(PATH, 'pileups.short',sep='/')
    print(fname)
    namelist <- as.character(read.csv(fname,sep='\t',header=F)[,1])
    data <- load_cov_distrs()
    data$nts <- as.numeric(data$cover * data$count)
    #data <- data[which(data$cover!=0),]
    data$count[which(data$cover==0)] <- 0
    res <- aggregate(data[,c("count", "nts")],by=list(data$names),sum)
    length <- data.rus.strains$length[1]
    res$redund <- res$nts / length
    res$uncov <- length - res$count
    res$length <- length
    res$frac <- res$count/length
    res$mean <- res$nts / res$count
    names(res) <- c("names", "cover", "nts", "redund", "uncov", "length", "frac", "mean")
    res2 <- res
    res2$names <- factor(res2$names,levels=namelist)
    res2 <- res2[order(res2$names),]
    res2$mean[which(is.nan(res2$mean))] <- 0.0
    outname <- paste(PATH,'means.list',sep='/')
    write.table(res2$mean, outname, sep='\t',row.names=FALSE,quote=FALSE,col.names=FALSE)
  return(res)
}

load_matrix <- function(filename, header=F){
  path <- paste(PATH,'/', filename, sep='')
  namepath <- paste(PATH,'/pileups.short',sep='')
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

load_distances <- function(suffix='alt_2', treshold){
  cross_cov_name <- paste('cross_cov_',suffix,'.out', sep='')
  cross_cov <- load_matrix(cross_cov_name, header = T)
  cross_cov <- cross_cov + t(cross_cov)
  cross_diff_name <- paste('cross_diffs_',suffix,'.out', sep='')
  cross_diff <- load_matrix(cross_diff_name, header = T)
  cross_diff <- cross_diff + t(cross_diff)
  temp <- c(as.character(data.rus.coverages$names[which(data.rus.coverages$frac > treshold)]),'REF')
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
  distance <- drop_outliers(distance)
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
  res <- merge(mdsdata, data.rus.meta, by.x=c("names"), by.y=c("names"), all.x=T,all.y=F)
  str(res)
  res$sample <- as.character(res$sample)
  #res$protocol[which(is.na(res$protocol))] <- 'non_protocol'
  res$sample[which(res$names == 'REF')] <- 'REF'
  res$sample <- as.factor(res$sample)
  res$sequen <- as.character(res$sequen)
  #res$protocol[which(is.na(res$protocol))] <- 'non_protocol'
  res$sequen[which(res$names == 'REF')] <- 'NONE'
  res$sequen <- as.factor(res$sequen)
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
    

