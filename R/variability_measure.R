#Измерение средних расстояний штаммов
#Важно исключать реплики
#То же касается теста Уилкоксона

treshold <- 0.6
setnums <- c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,47,48,49,56,57,58,59,62,68,73,74,75,82,85,87,89,96,98,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122)
res <- data.frame(setnum=c(),mean_dist=c(),nums=c())
for (setnum in setnums){
  dist_m <- load_distances(setnum,'alt_4',0.6,group='CHN',all=TRUE)[[3]]
  if (dim(dist_m)[1]>=10){
    temp <- data.frame(setnum=c(setnum),mean_dist=c(mean(dist_m[lower.tri(dist_m,diag=F)])),nums=c(dim(dist_m)[1]))
    res <- rbind(res,temp)
  }
}

res2 <- aggregate(data.cov_distrs$nts,by=list(data.cov_distrs$setnums),sum)
names(res2) <- c("setnum","total_nts")
m <- merge(res,res2,by.x=c("setnum"),by.y=c("setnum"),all.x=T,all.y=F)
m$set <- m$setnum
m <- merge(m,strain_meta.data,by.x=c("set"),by.y=c("set"),all.x=T,all.y=F)

ggplot(m,aes(x=length,y=mean_dist))+geom_line()+geom_point(data=m[m$setnum==2,],aes(x=length,y=mean_dist),color="red",size=5)
ggplot(m,aes(x=total_nts,y=mean_dist))+geom_line()+geom_point(data=m[m$setnum==2,],aes(x=total_nts,y=mean_dist),color="red",size=5)
###
library(gplots)
distance_m_list <- load_distances(0, suffix = 'alt_4', treshold = 0.3)[[3]]
heatmap.2(distance_m_list)


###


distance_histogram <- function(setnum,suffix='final_2', treshold =5000){
  #distance_m_list <- get_distance_m_list(setnum)
  distance_m_list <- load_distances(setnum, suffix = suffix, treshold = treshold)
  distance_m <- distance_m_list[[3]]
  nnames <- rownames(distance_m)
  subdata <- data.sample_names2[which(data.sample_names2$names %in% nnames),]
  good_names <- as.character(subdata$names[which(!duplicated(subdata$subj_id) | is.na(subdata$subj_id))])
  distance_m <- distance_m[good_names,good_names]
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