library(ggplot2)
library(reshape)
preprocess_dist_m <- function(dist_m){
  data <- melt(dist_m)
  names(data) <- c('names1','names2','distance')
  data <- data[which(data$names1!=data$names2),]
  return(data)
}

add_new_group <- function(data, group_frame){
  A1 <- group_frame
  names(A1) <- paste(names(A1),1,sep='')
  A2 <- group_frame
  names(A2) <- paste(names(A2),2,sep='')
  dm <- merge(data,A1,by.x=c('names1'),by.y=c('names1'),all.x=T,all.y=F)
  dm <- merge(dm,A2,by.x=c('names2'),by.y=c('names2'),all.x=T,all.y=F)
  dm <- dm[dm$names1!=dm$names2,]
  return(dm)
}

quantile_partition <- function(group_frame){
#   drop_ind <- which(is.na(group_frame[,2]))
#   if (length(drop_ind)>0){
#     group_frame <- group_frame[-drop_ind,]
#   }
  quantiles <- quantile(group_frame[,2])
  quantiles <- cut(group_frame[,2],quantiles,labels=c('1stQ','2ndQ','3rdQ','4thQ'))
  group_frame$quantiles <- quantiles
  return(group_frame)
}

get_distance_histograms <- function(dm){
  agg <- means_compute(dm)
  pp <- ggplot(dm,aes(x=distance,color=country1))+geom_histogram(binwidth=0.001)+facet_grid(country1 ~ country2)+geom_vline(aes(xintercept = x), data = agg)
  return(pp)
}

wilcox_test <- function(dm){
  countries <- levels(dm$country1)
  #print(countries)
  p_vals <- c()
  groups <- c()
  for (country1 in countries){
    for (country2 in countries){
      if (country2 != country1){
        groups <- c(groups,paste(country1,'_vs_',country2,sep=''))
        sub <- dm[dm$country1==country1,]
        p_val <- wilcox.test(sub$distance[sub$country2==country1],sub$distance[sub$country2==country2],alternative='less')
        p_vals <- c(p_vals,p_val$p.value)
      }
    }
  }
  res <- data.frame(groups=groups,p_values=p_vals)
  return(res)
}

distance_histogram <- function(setnum,suffix='final_2', treshold =5000){
  #distance_m_list <- get_distance_m_list(setnum)
  distance_m_list <- load_distances(setnum, suffix = suffix, treshold = treshold)
  distance_m <- distance_m_list[[3]]
  dm <- melt(distance_m)
  names(dm) <- c('names1','names2','distance')  
  A1 <- data.frame(names1=data.sample_names2$names,country1=data.sample_names2$country)
  A2 <- data.frame(names2=data.sample_names2$names,country2=data.sample_names2$country)
  dm <- merge(dm,A1,by.x=c('names1'),by.y=c('names1'),all.x=F,all.y=F)
  dm <- merge(dm,A2,by.x=c('names2'),by.y=c('names2'),all.x=F,all.y=F)
  dm <- dm[dm$names1!=dm$names2,]
  dm$country1 <- as.factor(as.character(dm$country1))
  print(dim(dm))
  p_vals <- wilcox_test(dm)
  agg <- means_compute(dm)
  dm$country <- NA
  inds <- which(dm$country1==dm$country2)
  dm$country[inds] <- as.character(dm$country1[inds])
#   pal <-  brewer.pal(4,"Set1")
#   pal <- pal[c(1,2,4)]
#   names(pal) <-levels(statdata$country)
  #pp <- ggplot(dm,aes(x=distance,color=country1))+geom_histogram(binwidth=0.001)+facet_grid(country1 ~ country2)+geom_vline(aes(xintercept = x), data = agg)
  ##Hot to mark mean line?
  ##Also need to make colors to be uniform
  return(list(dm,p_vals,agg))
}

get_p_vals <- function(setnums, suffix='final_2', treshold=5000){
  #len <- length(setnums)
  res <- data.frame(groups=c(),p_values=c(),set=c())
  for (setnum in setnums){
    print(setnum)
    temp <- try(distance_histogram(setnum=setnum,suffix=suffix,treshold=treshold)[[2]])
    if (class(temp) != "try-error"){
      temp$set <- setnum
      res <- rbind(res, temp)
    }
  }
  return(res)
}

cross_naming <- function(names1,names2){
  len <- length(names1)
  v <- rep(x='',times=len)
  f1 <- function(el){
    return(paste(el,'_vs_',el,sep=''))
  }
  f2 <- function(el){
    return(paste(el,'_vs_others',sep=''))
  }
  v[which(names1==names2)] <- sapply(names1[which(names1==names2)],f1)
  v[which(names1!=names2)] <- sapply(names1[which(names1!=names2)],f2)
  return(v)
}
#dir=NA,contignum=NA,mean_CHN,mean_EUR,mean_RUS,mean_USA,CHN_p,EUR_p,RUS_p,USA_P
p_vals.data <- data.frame(matrix(ncol=5, nrow= 0))
means.data <- data.frame(matrix(ncol=18, nrow= 0))

means_compute <- function(dm){
  agg <- aggregate(dm[,c('distance')],by=list(dm[,4],dm[,5]),FUN=mean)
  names(agg) <- c(names(dm[,c(4:5)]),'x')
  #delete repeated els
  return(agg)
}

ref_mean_dist <- function(setnums, prefix, treshold,selector='REF'){
  i <- 0
  for (setnum in setnums){
    data <- load_distances(setnum, prefix, treshold)
    dist_ref <- data.frame(names=names(data[[3]][selector,]), dist=data[[3]][selector,])
    dist_cov <- data.frame(names=names(data[[1]][selector,]), cov=data[[1]][,selector],diff=data[[2]][,selector])
    dist_ref <- merge(dist_ref, dist_cov,by.x=c("names"),by.y=c("names"),all.x=F,all.y=F)
    dist_ref <- merge(dist_ref, data.sample_names2,by.x=c("names"),by.y=c("names"),all.x=T,all.y=F)
    dist_ref$setnum <- setnum
    if (i==0){
      res <- dist_ref
    }
    else{
      res <- rbind(res, dist_ref)
    }
    i <- i + 1
  }
  res <- res[which(res$names != 'REF'),]
  means <- aggregate(res$dist, by=list(res$country,res$setnum),mean)
  medians <- aggregate(res$dist, by=list(res$country,res$setnum),median)
  sds <- aggregate(res$dist, by=list(res$country,res$setnum),sd)
  agr_res <- data.frame(country=means[,1], setnum=means[,2], mean=means[,3], median=medians[,3], sd=sds[,3])
  return(list(res,agr_res))
}


