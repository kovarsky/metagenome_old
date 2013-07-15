PATH <- '/data1/bio/runs-kovarskiy/metagenome_data/'
#PATH <- '/data1/bio/runs-kovarskiy/onkokids/snpdata'
compute_snp_rates <- function(setnums){
  res <- data.frame(names=c(),cover=c(),count=c(),snpcount=c(),frac=c(),setnum=c())
  temp <- data.cov_distrs[,c("names", "redund","setnums","mean","nts")]
  temp$len_frac <- data.cov_distrs$frac
  temp$tcovd <- data.cov_distrs$cover
  for (setnum in setnums){
    print(setnum)
    snp_distr <- load_snp_distrs(setnum)
    cov_distr <- load_cov_distrs(setnum)
    merged_data <- merge(cov_distr[cov_distr$names %in% snp_distr$names,],snp_distr,by.x=c("names","cover"),by.y=c("names","cover"),all.x=T,all.y=T)
    merged_data <- merged_data[-which(merged_data$cover==0),]
    merged_data$snpcount[which(is.na(merged_data$snpcount))] <- 0
    merged_data$frac <- merged_data$snpcount/merged_data$count
    merged_data$setnum <- setnum
    merged_data <- merge(temp[temp$setnums==setnum,],merged_data,by.x=c("names"),by.y=c("names"),all.x=T,all.y=T)
    res <- rbind(res,merged_data)
  }
  #res$frac[which(is.na(res$frac))] <- 0
  return(res)
}

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
      total <- strain_meta.data$length[strain_meta.data$set == setnum]
      delta <- total - summ
      print(delta)
      add_info <- data.frame(cover=c(0),count=c(delta))
      temp <- rbind(temp, add_info)
      temp$names <- fname
      res <- rbind(res,temp)
    }
    else{
      total <- strain_meta.data$length[strain_meta.data$set == setnum]
      temp <- data.frame(cover=c(0),count=c(total))
      temp$names <- fname
      res <- rbind(res,temp)
    }
  }
  return(res)
}

load_hetero_distrs <- function(setnum){
  path <- paste(PATH,setnum,sep='')
  files <- list.files(path, setnum, pattern = '\\.hethist2',full.names=T)
  res <- data.frame(cover=c(), hetero=c(), hetcount=c(), names=c())
  for (file in files){
    fname <- substr(basename(file),1,nchar(basename(file))-9)
    temp <- try(read.csv(file,sep='\t',header=F))
    #temp <- read.csv(file,sep='\t',header=F)
    if (class(temp) != "try-error"){
      names(temp) <- c("cover", "hetero", "hetcount")
      #summ <- sum(temp$count)
      temp$names <- fname
      res <- rbind(res,temp)
    }
    else{
      #total <- strain_meta.data$length[strain_meta.data$set == setnum]
      #temp <- data.frame(cover=c(0),count=c(total))
      #temp$names <- fname
      #res <- rbind(res,temp)
    }
  }
  return(res)
}

load_hetero_spectras <- function(covdata){
  len <- dim(covdata)[2]
  res <- data.frame(eq=c(),max=c(),total=c(),names=c(),redund=c())
  for (i in 1:len){
    print(i)
    setnum <- covdata$setnums[i]
    fname <- covdata$names[i]
    redund <- covdata$redund[i]
    path <- paste(PATH,setnum,'/',sep='/')
    full_name <- paste(path,fname,'.het',sep='')
    temp <- read.csv(full_name,sep='\t',header=F)
    names(temp) <- c("eq","max","total")
    temp$names <- fname
    temp$redund <- redund
    res <- rbind(res,temp)
  }
  res$frac <- res$max/res$total
  inds <- which(res$eq == res$max)
  res$frac[inds] <- 1 - res$frac[inds]
  res$tfrac <- res$total/res$redund
  res$bins <- cut(res$frac,breaks=c(seq(0,1,0.1)))
  res$tbins <- cut(res$tfrac,breaks=c(seq(0,4,0.1)))
  return(res)
}

load_snp_counts <- function(setnums,suffix='test__2_1',treshold=1){
  res <- data.frame(cover=c(),freq=c(),setnum=c(),count=c())
  for (setnum in setnums){
    print(setnum)
    path <- paste(PATH,setnum,'/snp_counts_',suffix,sep='')
    data <- read.csv(path,sep='\t',header=F)
    data <- data[,c(6,7)]
    names(data) <- c("summ","snp")
    data <- data[data$snp>treshold,]
    data <- count(data$summ)
    names(data) <- c("cover","freq")
    path2 <- paste(PATH,setnum,'/poscount_',suffix,sep='')
    data2 <- read.csv(path2,sep='\t',header=F)
    names(data2) <- c("cover","count")
    m <- merge(data,data2,by.x=c('cover'),by.y=c('cover'),all.x=T,all.y=F)
    m$setnum <- setnum
    res <- rbind(res,m)
  }
  res$frac <- res$freq/res$count
  return(res)
}

load_snp_read_counts <- function(setnums,suffix='test__2_1',treshold=1){
  res <- data.frame(cover=c(),freq=c(),setnum=c(),count=c())
  for (setnum in setnums){
    print(setnum)
    path <- paste(PATH,setnum,'/snp_counts_',suffix,sep='')
    data <- read.csv(path,sep='\t',header=F)
    path3 <- paste(PATH,setnum,'/mutual_freqs_',suffix,sep='')
    temp <- read.csv(path3,sep='\t',header=F)
    temp <- temp[,2] + temp[,3] + temp[,4] + temp[,5]
    data <- data[,c(6,7)]
    names(data) <- c("summ","snp")
    data$summ <- temp
    data <- data[data$snp>treshold,]
    data <- count(data$summ)
    names(data) <- c("cover","freq")
    path2 <- paste(PATH,setnum,'/read_counts_',suffix,sep='')
    data2 <- read.csv(path2,sep='\t',header=F)
    names(data2) <- c("cover","count")
    m <- merge(data,data2,by.x=c('cover'),by.y=c('cover'),all.x=T,all.y=F)
    m$setnum <- setnum
    res <- rbind(res,m)
  }
  res$frac <- res$freq/res$count
  return(res)
}

load_snp_distrs <- function(setnum){
  path <- paste(PATH,setnum,sep='')
  files <- list.files(path, setnum, pattern = '\\.refsnphist',full.names=T)
  res <- data.frame(cover=c(), count=c(), names=c())
  for (file in files){
    fname <- substr(basename(file),1,nchar(basename(file))-11)
    temp <- try(read.csv(file,sep='\t',header=F))
    #temp <- read.csv(file,sep='\t',header=F)
    if (class(temp) != "try-error"){
      names(temp) <- c("cover", "snpcount")
      summ <- sum(temp$snpcount)
      total <- strain_meta.data$length[strain_meta.data$set == setnum]
      delta <- total - summ
      add_info <- data.frame(cover=c(0),snpcount=c(delta))
      temp <- rbind(temp, add_info)
      temp$names <- fname
      res <- rbind(res,temp)
    }
    else{
      total <- strain_meta.data$length[strain_meta.data$set == setnum]
      add_info <- data.frame(cover=c(0),snpcount=c(total))
      add_info$names <- fname
      res <- rbind(res,add_info)
    }
  }
  return(res)
}



aggregate_cov_distrs <- function(setnums){
  i <- 1
  for (setnum in setnums){
    print(setnum)
    data <- load_cov_distrs(setnum)
    data$nts <- as.numeric(data$cover * data$count)
    #data <- data[which(data$cover!=0),]
    data$count[which(data$cover==0)] <- 0
    res <- aggregate(data[,c("count", "nts")],by=list(data$names),sum)
    length <- strain_meta.data$length[strain_meta.data$set == setnum]
    res$redund <- res$nts / length
    res$uncov <- length - res$count
    res$length <- length
    res$setnum <- setnum
    res$frac <- res$count/length
    res$mean <- res$nts / res$count
    names(res) <- c("names", "cover", "nts", "redund", "uncov", "length", "setnums","frac", "mean")
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
    print(dim(res))
    length <- strain_meta.data$length[strain_meta.data$set == setnum]
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

cum_func <- function(setnum, param_num){
  total_len <- strain_meta.data$length[which(strain_meta.data$set == args$setnum)]
  cov_distr <- load_cov_distrs(setnum)
  fname <- paste(PATH,setnum,'/stat_test_res',param_num,'.csv',sep='')
  params <- read.csv(fname,sep='\t')
  cov_distr <- cov_distr[cov_distr$names%in%params$fname,]
  func <- function(x, param_num, args, len){
    res <- apply(sapply(1:param_num, function(i){args[1,i+param_num] * len * dpois(x, args[1,i])}), 1, sum)
    return(res)
  }
  entry_num <- length(params$fname)
  for (i in 1:entry_num){
    fname <- as.character(params$fname[i])
    max_arg <- max(cov_distr$cover[cov_distr$names==fname])
    args <- params[params$fname==fname,]
    counts <- func(0:max_arg, param_num, args, total_len)
    temp <- data.frame(cover=c(0:max_arg),count=counts,names=fname)
    if (i == 1){
      res <- temp
    }
    else {
      res <- rbind(res, temp)
    }
  }
  return(list(res, cov_distr))
}

get_distr_parameters <- function(setnum, param_num,zero=T,frac=3,pattern="\\.hist$"){
  path <- paste(PATH,setnum,sep='')
  files <- list.files(path, setnum, pattern = pattern,full.names=T)
  res <- data.frame(cover=c(), count=c(), names=c())
  init_params <- c(0, 10, 11, 12)
  weights <- c(0.2, 0.21, 0.22, 0.23)
  res_weights <- matrix(0,nrow=0,ncol=param_num)
  res_params <- matrix(0,nrow=0,ncol=param_num)
  fnames <- data.frame(fname=c(),frac=c())
  for (file in files){
    fname <- substr(basename(file),1,nchar(basename(file))-5)
    print(fname)
    temp <- try(read.csv(file,sep='\t',header=F))
    if (class(temp)!='try-error'){
      names(temp) <- c("cover", "count")
      summ <- sum(temp$count)
      total <- strain_meta.data$length[strain_meta.data$set == setnum]
      delta <- total - summ
      add_info <- data.frame(cover=c(0),count=c(delta))
      temp <- rbind(temp, add_info)
      total_frac <- sum(temp$cover * as.numeric(temp$count)) / total
      if (total_frac >= frac){
        #       flag <-  TRUE
        #       while(flag){
        #         
        #       }
        params <- try(parameter_estimation(init_params[1:param_num],weights[1:param_num],temp$cover,temp$count,total,zero=zero))
        if (class(params) != "try-error"){
          res_weights <- rbind(res_weights, params[[2]])
          res_params <- rbind(res_params, params[[1]])  
          fnames_temp <- data.frame(fname=fname,frac=total_frac)
          fnames <- rbind(fnames, fnames_temp)
        }
      }
      #temp$names <- fname
      #res <- rbind(res,temp)
    }
  }
  return(list(res_params, res_weights, fnames))
}

general_test <- function(setnums,zero=T){
  for (param_num in 2:2){
    i <- 1
    for (setnum in setnums){
      params <- get_distr_parameters(setnum, param_num,zero=zero)
      temp_table <- cbind(as.data.frame(params[[1]]),as.data.frame(params[[2]]),params[[3]])
      temp_table$setnum <- setnum
      
      outname_temp <- paste(PATH,'/', setnum, '/stat_test_res',param_num,'.csv',sep='')
      write.table(temp_table, outname_temp, sep='\t')
      if (i == 1){
        res <- temp_table
      }
      else{
        res <- rbind(res, temp_table)
      }
      i <- i + 1
    }
    outname <- paste(PATH,'stat_test_res',param_num,'.csv',sep='')
    write.table(res, outname, sep='\t')
  }
}

tail_poisson <- function(x, lambda){
  hlim <- 1000

  tail_p <- function(x){
    sum(sapply(c((x+1):hlim),function(i){dpois(i, lambda) / i}))
  }
  res <- sapply(x,tail_p)
  return(res/sum(res))
}

tail_nbinom <- function(x,params){
  hlim <- 1000

  tail_p <- function(x){
    sum(sapply(c((x+1):hlim),function(i){dnbinom(i, params[1],params[2]) / i}))
  }
  res <- sapply(x,tail_p)
  return(res/sum(res))
}


loglik <- function(lambdas, x, y, weights){
  len <- length(weights)
  summs <- sapply(1:len,function(arg){weights[arg] * dpois(x, lambdas[arg])})
  summs <- apply(summs,1,sum)
  good_inds <- which(summs!=0)
  res <- sum (log(summs[good_inds]) * y[good_inds] )
  return(-res)
}

loglik2 <- function(lambdas, x, y, weights){
  summs <- weights * dnbinom(x, lambdas[1], lambdas[2])
  good_inds <- which(summs!=0)
  res <- sum (log(summs[good_inds]) * y[good_inds] )
  return(-res)
}

get_weights <- function(distrs, x, y, weights, L, params){
  len <- length(distrs)
  M <- sapply(1:len, function(ind) {return(weights[ind] * distrs[[ind]](x, params[[ind]]))})
  summ <- apply(M,1,sum)
  #print('----')
  #print(M)
  #print(weights)
  #print(params)
  good_inds <- which(summ!=0)
  gammas <- M[good_inds,] / summ[good_inds]
  #print(good_inds)
  weights <- apply(gammas, 2, function(gamma){return(sum(gamma * y[good_inds]))}) / L
  return(weights)
}

parameter_estimation_v2 <- function(distr_num, x, y, L,zero=T){
  distrs <- c(dpois, tail_poisson, dpois, dpois)[1:distr_num]
  params <- c(0, 10, 11, 12)[1:distr_num]
  weights <- rep(1/distr_num, distr_num)
  print(weights)
  y <- y/L
  for (i in 1:20){
    weights <- get_weights(distrs, x, y, weights, L, params)
    delta_y <- (y - weights[1] * distrs[[1]](x, params[1]) - weights[2] * distrs[[2]](x, params[2]))
    out <- optim(par = params[3:distr_num], fn=loglik,x=x,y=delta_y, weights=weights[3:distr_num], method = "L-BFGS-B",lower=c(0,0))
    print(out$par)
    params[3:distr_num] <- out$par
    params[2] <- params[3]
    # if (zero){
    #   params[1] <- 0
    # }
    print('-----------')
    print(params)
    print(weights)
  }
  return(list(params, weights))
}

parameter_estimation_v3 <- function(distr_num, x, y, L,zero=T){
  dnbinom2 <- function(x,params){
    return(dnbinom(x,params[1],params[2]))
  }
  distrs <- c(dpois, tail_nbinom, dnbinom2, dpois)[1:distr_num]
  params <- list(0, c(11,1/2), c(11,1/2), 12)[1:distr_num]
  weights <- rep(1/distr_num, distr_num)
  print(weights)
  y <- y/L
  for (i in 1:20){
    weights <- get_weights(distrs, x, y, weights, L, params)
    delta_y <- (y - weights[1] * distrs[[1]](x, params[[1]]) - weights[2] * distrs[[2]](x, params[[2]]))
    out <- optim(par = params[[3]], fn=loglik2,x=x,y=delta_y, weights=weights[3], method = "L-BFGS-B",lower=c(0,0.01),upper=c(200,0.9))
    print(out$par)
    params[[3]] <- out$par
    params[[2]] <- params[[3]]
    # if (zero){
    #   params[1] <- 0
    # }
    print('-----------')
    print(params)
    print(weights)
  }
  return(list(params, weights))
}



parameter_estimation <- function(params, weights, x, y, L,zero=F){
  len <- length(weights)
  get_distrs <- function(params){
    res <- list()
    for (p in params){
      res <- c(res, dpois, recursive=F)
    }
    return(res)
  }
  distrs <- get_distrs(params)
  for (i in 1:20){
    weights <- get_weights(distrs, x, y, weights, L, params)
    out <- optim(par = params, fn=loglik,x=x,y=y, weights=weights, method = "L-BFGS-B",lower=c(0,0))
    params <- out$par
    if (zero){
      params[1] <- 0
    }
    print('-----------')
    print(params)
    print(weights)
    distrs <- get_distrs(params)
  }
  return(list(params, weights,distrs))
}

