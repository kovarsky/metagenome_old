PATH <- '/data1/bio/runs-kovarskiy/metagenome_data/'
library(reshape)

diversity_load <- function(setnum, treshs=c(0,2,5,10), flag=''){
  
  diversity_name <- paste(PATH,setnum,'/diversity',flag,'.tsv',sep='')
  data.diversity <- read.csv(diversity_name, sep='\t', header=F)
  tresholds <- treshs
  len <- length(tresholds)
  res <- data.frame(names=c(),cover=c(),diversity=c(),tresh=c())
  for (i in 1:len){
    temp <- data.diversity[,c(1,1+i,1+i+len)]
    names(temp) <- c('names','cover','diversity')
    temp$tresh <- tresholds[i]
    res <- rbind(res,temp)
  }
  res$diversity <- res$diversity / res$cover
  
  res.cumulative <- res[res$names=='cumulative',]
  res <- merge(res,data.sample_names2,by.x=c("names"),by.y=c("names"),all.x=F, all.y=F)
  
  return(list(res, res.cumulative))
}

cross_hetero_load <- function(setnum, tresh){
  library(reshape)
  diversity_name <- paste(PATH,setnum,'/diversity.tsv',sep='')
  data.diversity <- read.csv(diversity_name, sep='\t', header=F)
  tresholds <- c(0,2,5,10)
  len <- length(tresholds)
  samples <- data.diversity[,1]
  samples <- samples[which(samples!='cumulative')]
  
  crosscov_name <- paste(PATH,setnum,'/cross_cov_',tresh,'.out',sep='')
  data.crosscov <- as.matrix(read.csv(crosscov_name, sep='\t', header=F))
  data.crosscov <- data.crosscov + t(data.crosscov)
  dimnames(data.crosscov) <- list(samples,samples)
  crosshetero_name <- paste(PATH,setnum,'/cross_hetero_',tresh,'.out',sep='')
  data.crosshetero <- as.matrix(read.csv(crosshetero_name, sep='\t', header=F))
  data.crosshetero <- data.crosshetero + t(data.crosshetero)
  dimnames(data.crosshetero) <- list(samples, samples)
  data.res <- data.crosshetero / data.crosscov
  
  melted.data.res <- melt(data.res)
  names(melted.data.res) <- c('samples1','samples2','heterogenity')
  melted.data.crosscov <- melt(data.crosscov)
  names(melted.data.crosscov) <- c('samples1','samples2','crosscov')
  melted.res <- merge(melted.data.res, melted.data.crosscov, by.x = c('samples1', 'samples2'), by.y = c('samples1', 'samples2'),all.x=T,all.y=T)
  
  return(list(data.res, melted.res))
}


get_err_occurence <- function(setnum,prefix){
  filename <- paste(PATH,setnum,'/counts',prefix,'.tsv',sep='')
  print('...reading')
  data <- read.csv(filename, sep='\t', header=F)
  print('Done!')
  len <- dim(data)[2]
  data <- data[,c(len-1,len)]
  data$diff <- data[,2] - data[,1]
  data <- data[data$diff<=20,]
  means <- aggregate(data$diff,by=list(data[,1]),FUN=mean)
  counts <- count(data[,1])
  sds <- aggregate(data$diff,by=list(data[,1]),FUN=sd)
  res <- cbind(means,sds[,2],counts[,2])
  names(res) <- c('cover','mean','sd','count')
  res$flag <- prefix
  return(res)
}

heterogen_spectrum_load <- function(setnum, prefix='', breaks=(seq(0.05,0.5,0.05))){
  filename <- paste(PATH,setnum,'/alt_spectra',prefix,'.out',sep='')
  samplefilename <- paste(PATH,setnum,'/samplenames.list',sep='')
  data <- read.csv(filename, sep='\t', header=F)
  len <- length(breaks)
  maxlen <- dim(data)[2]
  res <- data[,c(1:len, maxlen)]
  res[,len] <- res[,len] + data[,c(len + 1)]
  names(res) <- c(breaks, 'total')
  samplenames <- read.csv(samplefilename, sep='\t', header=F)
  res$names <- samplenames[,1]
  res <- melt(res, id=c("names",'total'))
  names(res) <- c("names","total","treshold","altnum")
  res$treshold <- as.numeric(as.character(res$treshold))
  res$setnum <- setnum
  res <- merge(res, data.sample_names2, by.x=c("names"),by.y=c("names"),all.x=T,all.y=T)
  na_ind <- which(is.na(res$short_names))
  res$country <- as.character(res$country)
  res$country[na_ind] <- 'RUS'
  res$short_names[na_ind] <- as.character(res$names[na_ind])
  res$groups[na_ind] <- substr(res$names[na_ind],1,3)
  res$frac <- res$altnum / res$total
  return(res)
}

heterogen_spectra_load <- function(setnums, prefix='', breaks=(seq(0.05,0.5,0.05)),tresh=5000){
  len <- length(setnums)
  for (i in 1:len){
    if (i == 1){
      res <- heterogen_spectrum_load(setnums[i], prefix=prefix, breaks=breaks)
    }
    else{
      temp <- heterogen_spectrum_load(setnums[i], prefix=prefix, breaks=breaks)
      res <- rbind(res, temp)
    }
  }
  return(res)
}

strain_aggregation <- function(spectra, total_treshold=5000, alt_treshold=0.1){
  spectra <- spectra[which((spectra$treshold>=alt_treshold) & (spectra$total>=total_treshold)),]
  res <- aggregate(spectra$frac,by=list(spectra$names,spectra$setnum,spectra$total,spectra$country),FUN=sum)
  names(res) <- c("names","setnum","total","country","frac")
  res_by_set <- aggregate(res$frac,by=list(res$setnum),FUN=mean)
  names(res_by_set) <- c("setnum","mean_frac")
  return(list(res_by_set, res))
}