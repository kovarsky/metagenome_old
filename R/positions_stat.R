stat_table <- function(dir,prefix){
  file <- paste(dir,prefix,'varpos_variants.txt',sep='')
  data <- read.csv(file,header=F,sep='\t')
  names(data) <- c('chrom','pos','A','T','G','C')
  norm_func <- function(v){
    return(sqrt(sum(v*v)))
  }
  data$norm <- apply(data[,3:6],1,norm_func)
  data$sum <- apply(data[,3:6],1,sum)
  #data$diff <- data$sum - data$max
  #data$frac <- data$diff/data$sum
  data$flag <- prefix
  print (prefix)
  print ('--done!--')
  return(data)
}

filter_data <- function(data,covtresh=0.5,fractresh=0.1){
  maxcov <- max(data$sum)
  data <- data[-which(is.na(data$frac)),]
  data <- data[which(data$frac > fractresh),]
  data <- data[which((data$sum/maxcov) > covtresh),]
  return(data)
}

merge_data <- function(data_list){
  len <- length(data_list)
  if (len < 2){
    return(data_list[[1]]) 
  }
  else{
    data <- data_list[[1]]
    for (i in 2:len){
      data <- rbind(data,data_list[[i]])
    }
    return(data)
  } 
}

read_all_tables <- function(dir,prefixes){
  data_list <- list()
  for (prefix in prefixes){
    data_list[[prefix]] <- stat_table(dir,prefix)
  }
  return(data_list)
}

pos_processing <- function(data_list,covtresh,fractresh){
  len <- length(data_list)
  for (i in 1:len){
    data_list[[i]] <- filter_data(data_list[[i]],covtresh,fractresh)
  }
  res <- merge_data(data_list)
}

comm_pos <- function(data1,data2){

  label1 = data1$flag[1]
  label2 = data2$flag[1]
  len <- dim(data1)[1]
  res <- data.frame(res=rep(0, len))
  for (i in 3:6){
    buff <- data1[,i] * data2[,i]
    #print(buff)
    res$res <- res$res + buff
    #res <- cbind(res,data1[,i])
    #res <- cbind(res,data2[,i])
    #print(res)
  }

  norm <- data1$norm * data2$norm
  res$res <- sqrt(res$res/norm)
  #res$norm <- norm
  res$pos <- data1$pos
  res[,label1] <- data1$sum
  res[,label2] <- data2$sum
  res$chrom <- data1$chrom
  ind <- which(res$res < 0.5)
  print(length(ind))
  res <- res[ind,]
  ind <- which(res[,label1]>10 & res[,label2]>10)
  print(length(ind))
  return(res[ind,c("chrom","pos")])
}

whole_comm_pos <- function(data_list, name = '/data1/bio/runs-kovarskiy/additional_pileup_files/temp_pos.list'){
  len <- length(data_list)
  combs <- combn(len,2)
  num <- dim(combs)[2]
  for (i in 1:num){
    buff <- comm_pos(data_list[[combs[1,i]]], data_list[[combs[2,i]]])
    if (i > 1){
      res <- merge(res,buff,by.x=c('chrom','pos'),by.y=c('chrom','pos'),all.x=T,all.y=T)
    }
    else{
      res <- buff
    }
  }
  write.table(x=res,file=name,sep='\t',col.names=F,row.names=F,quote=F)
  return(res)
}

read_new_tables <- function(dir,prefix)
{
  filename <- paste(dir,prefix,'newvarpos_variants.txt',sep='/')
  data <- read.csv(filename,sep='\t',header=F)
  names(data) <- c('chroms','pos','A','T','G','C','sum','diff','max_alt')
  covfilename <- paste(dir,prefix,'common_coverage_summary.list',sep='/')
  covdata <- read.csv(covfilename,sep='\t',header=F)
  names(covdata) <- c('sum','times')
  return(list(data,covdata))
}

read_cov_table <- function(dir,prefix)
{
  covfilename <- paste(dir,prefix,'cover_summ.list',sep='/')
  covdata <- read.csv(covfilename,sep='\t',header=F)
  names(covdata) <- c('sum','times')
  return(covdata)
}

plot_cov_stat <- function(strainnums){
  len <- length(strainnums)
  j = 0
  for (i in strainnums){
    print(i)
    setnum <- as.character(strain_meta.data$set[i])
    prefix <- strain_meta.data$str_subset[i]
    strainname <- paste(strain_meta.data$short_descript[i],as.character(i),sep='_')
    dir <- paste(CONSTDIR,setnum,sep='')
    covdata <- read_cov_table(dir,prefix)
    covdata$strain_num <- i
    covdata$strainname <- strainname
    if (j == 0){
      res <- covdata
    }
    else{
      res <- rbind(res,covdata)
    }
    j <- j + 1
  }
  return(res)
}

get_cdf <- function(datalist, max_alt_flag=F){
  #maxes <- apply(data[,c('A','T','G','C')],1,max)
  #data$diff <- data$sum - maxes
  data <- datalist[[1]]
  covdata <- datalist[[2]]
  data$ind <- 1
  sum_aggr <- aggregate(data$ind, by=list(data$sum),sum)
  names(sum_aggr) <- c('sum','total_num')
  if (max_alt_flag){
    sum_diff_aggr <- aggregate(data$ind,by=list(data$max_alt,data$sum),sum)
  }
  else{
    sum_diff_aggr <- aggregate(data$ind,by=list(data$diff,data$sum),sum)
  }
  names(sum_diff_aggr) <- c('diff','sum','num')
  merged_data <- merge(sum_diff_aggr,covdata,by.x=c('sum'),by.y=c('sum'),all.x=T,all.y=F)
  merged_data$freq <- merged_data$num/merged_data$times
  merged_data$variab <- merged_data$diff/merged_data$sum
  merged_data <- merged_data[which(merged_data$diff != 0),]
  return(merged_data)
  #ecdf_data$freq <- ecdf_data$num/ecdf_data$total_num
  #return(data)
}

aggregation <- function(data, tresholds){
  agg <- function(tresh){
    temp <- data[data$variab > tresh,]
    res <- aggregate(temp$freq,by=list(temp$sum),sum)
    freqname<-paste('freq',as.character(tresh),sep='')
    names(res) <- c('sum',freqname)
    return(res)
  }
  len <- length(tresholds)
  for (i in 1:len){
    tresh <-tresholds[i]
    print(i)
    if (i==1){
      res <- agg(tresh)
    }
    else{
      res <- merge(res,agg(tresh),by.x=c('sum'),by.y=c('sum'),all.x=F,all.y=F)
    }
  }
  return(res)
}

aggregation_by_diff <- function(data, tresholds){
  agg <- function(tresh){
    temp <- data[data$diff > tresh,]
    res <- aggregate(temp$freq,by=list(temp$sum),sum)
    numname<-paste('more',as.character(tresh),sep='')
    names(res) <- c('sum',numname)
    return(res)
  }
  len <- length(tresholds)
  for (i in 1:len){
    tresh <-tresholds[i]
    print(i)
    if (i==1){
      res <- agg(tresh)
    }
    else{
      res <- merge(res,agg(tresh),by.x=c('sum'),by.y=c('sum'),all.x=F,all.y=F)
    }
  }
  return(res)
}

CONSTDIR <- '/data1/bio/runs-kovarskiy/metagenome_data/'

plot_position_stat <- function(strainnums, by.diffs = F, max_alt_flag=F){
  len <- length(strainnums)
  j = 0
  for (i in strainnums){
    print(i)
    num <- which(strain_meta.data$set == i)
    setnum <- as.character(strain_meta.data$set[num])
    #prefix <- strain_meta.data$str_subset[num]
    strainname <- paste(strain_meta.data$short_descript[num],num,sep='_')
    #dir <- paste(CONSTDIR,setnum,sep='')
    datalist <- read_new_tables(CONSTDIR,setnum)
    merged_data <- get_cdf(datalist,max_alt_flag)
    if (by.diffs){
      agg <- aggregation_by_diff(merged_data,c(0,1,2,5))
    }
    else{
      agg <- aggregation(merged_data,c(0.0,0.01,0.02,0.05))
    }
    agg$strain_num <- i
    agg$strainname <- strainname
    if (j == 0){
      res <- agg
    }
    else{
      res <- rbind(res,agg)
    }
    j <- j + 1
  }
  return(res)
}

cover_stat <- function(setnum){
  path <- '/data1/bio/runs-kovarskiy/metagenome_data/'
  statfilename <- paste(path, setnum, sep='','/total_sample_coverages.tsv')
  statdata <- read.csv(statfilename, sep=' ', head=F)
  namesfile <- paste('/data1/bio/runs-kovarskiy/metagenome_data',setnum,'pileups.list',sep='/')
  snames <- read.csv(namesfile, sep='', header=F)
  snames <- snames[,1]
  names(statdata) <- c('names','covered_positions','sum_of_cover_depth')
  statdata$names <- snames
  statdata <- merge(x=statdata,y=data.sample_names2,by.x=c('names'),by.y=c('names'),all.x=T,all.y=F)
  return(statdata)
}

cover_stats <- function(setnums){
  len <- length(setnums)
  for (i in 1:len){
    data <- try(cover_stat(setnums[i]))
    if (typeof(data)!='character'){
      summn <- sum(as.numeric(data$sum_of_cover_depth))
      summc <- sum(as.numeric(data$covered_positions))
      temp <- data.frame(setnum=setnums[i],summn=summn,summc=summc)
      print(temp)
      if (i == 1){
        res <- temp
      }
      else{
        res <- rbind(res,temp)
      }
    }
  }
  return(res)
}

finstat_plot <- function(position_stat_data){
  strain_num_ind <- which(names(position_stat_data)=='strain_num')
  temp <- melt(position_stat_data[,-strain_num_ind],id.vars=c('sum','strainname'))
  p <- ggplot(temp,aes(x=sum,y=value,group=variable,fill=variable))+geom_density(stat="identity",alpha=0.5,size=0.3)+facet_wrap(~strainname)+scale_y_continuous(limits=c(0,0.5))
  return(p)
}

snp_rate_stat <- function(setnums, treshs){
  len <- length(setnums)
  for (i in 1:len){
    setnum <- setnums[[i]]
    snp_rate <- snp_rate_load(setnum, treshs)
      if (i==1){
        res <- snp_rate
      }
      else{
        res <- rbind(res, snp_rate)
      }
    }
  return(res)
}

snp_rate_load <- function(setnum, treshs=c(5)){
  print(setnum)
  snp_stat_file <- paste('/data1/bio/runs-kovarskiy/metagenome_data',setnum,'new_snp_stat2.tsv',sep='/')
  namesfile <- paste('/data1/bio/runs-kovarskiy/metagenome_data',setnum,'pileups.list',sep='/')
  snames <- read.csv(namesfile, sep='', header=F)
  snames <- snames[,1]
  #snames <- rep(snames, rep(length(treshs),length(snames)))
  new_snp_stat <- read.csv(snp_stat_file, sep='', header=F)
  names(new_snp_stat) <- c('names','varpos_cov','snps')
  new_snp_stat$treshs <- treshs
  cov_file <- paste('/data1/bio/runs-kovarskiy/metagenome_data',setnum,'total_sample_coverages.tsv',sep='/')
  cov_stat <- read.csv(cov_file, sep='', header=F)
  names(cov_stat) <- c('tcov','tnts')
  
  cov_stat$names <- snames
  res <- merge(new_snp_stat, cov_stat, by.x=c('names'), by.y=c('names'), all.x=F, all.y=F)
  res$frac <- res$snps/res$tcov
  res <- merge(res, data.sample_names2, by.x=c('names'), by.y=c('names'),all.x=F,all.y=F)
  res$setnum <- setnum
  return(res)
}

snp_rate_plot <- function(setnums,treshold, covtreshold=10^6){
  cumul_cov <- get_cumul_cov(setnums)
  res <-  snp_rate_stat(setnums,c(0.005, 0.01, 0.02, 0.05))
  rate_stat <- res
  res <- res[which(res$treshs==treshold),]
  res <- res[which(res$tcov>covtreshold),]
  agg <- aggregate(res$frac,by=list(res$setnum),mean)
  res$set <- res$setnum
  res <- res[,-11]
  res <- merge(res,cumul_cov,by.x=c("set"),by.y=c("set"),all.x=T,all.y=T)
  #names(agg) <- c("set",'snpfrac')
  #agg <- merge(agg,cumul_cov,by.x=c("set"),by.y=c("set"),all.x=T,all.y=T)
  return(list(res, rate_stat))
}

get_cumul_cov <- function(setnums){
  res <- data.frame(set=c(), possum=c(), ntsum=c(), samplnum=c())
  for (setnum in setnums){
    fname <- paste('/data1/bio/runs-kovarskiy/metagenome_data',setnum,'total_sample_coverages.tsv',sep='/')
    data <- read.csv(fname, sep=' ', header=F)
      possum <- sum(as.numeric(data[,2]),na.rm=T)
      ntsum <- sum(as.numeric(data[,3]),na.rm=T)
      if (is.na(ntsum)){
        print(data)
      }
      samplnum <- length(which(!is.na(data[,3])))
      temp <- data.frame(set=setnum, possum=possum, ntsum=ntsum, samplnum=samplnum)
      res <- rbind(res, temp)
    }
  res <- merge(res,strain_meta.data[,c('set','length','short_descript')],by.x=c("set"),by.y=c("set"),all.x=F,all.y=F)
  return(res)
}


#plot: ggplot(res,aes(x=flag,ymin=pos,ymax=pos+100,size=50,color=flag)) + geom_linerange() + coord_flip() + scale_size(range = c(1, 50)) + scale_alpha(range=c(0.1,0.5))