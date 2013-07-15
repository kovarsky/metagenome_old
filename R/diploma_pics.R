#saturation curves:
ggplot(data.cov_distrs,aes(x=redund,y=frac))+geom_point()+facet_wrap(~setnums,ncol=7,scales="free")


#partly:
nums <- c( 41,32,36,32,6,9,110,47,2)
nums <- c(6)
p <- ggplot(data.cov_distrs[data.cov_distrs$setnums %in% nums,],aes(x=redund,y=frac))+geom_point()+facet_wrap(~setnums,ncol=3,scales="free")
p <- p + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(size=15),axis.text.y  = element_text(size=15),axis.title.y = element_text(size=20))+xlab(label='С, среднее покрытие')
#p <- p + ylab(label=paste(expression('F'[0]),', покрытая часть генома',sep=''))
p + ylab(label='F, покрытая часть генома')

#partly with colors
nums <- c(0,2,21,28,34)
m <- merge(data.cov_distrs[data.cov_distrs$setnums %in% nums,],data.sample_names2,by.x=c("names"),by.y=c("names"),all.x=T,all.y=F)
p <- ggplot(m,aes(x=redund,y=frac,color=country))+geom_point()+facet_wrap(~setnums,ncol=3,scales="free")
p <- p + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(size=15),axis.text.y  = element_text(size=15),axis.title.y = element_text(size=20))+xlab(label='С, среднее покрытие')
p + scale_x_continuous(lim=c(0,50))

#Проверка отношенией двух определений среднего покрытия
subdata <- data.cov_distrs[data.cov_distrs$setnums %in% nums,]
str(subdata)
subdata$nfrac <- subdata$redund/subdata$mean
ggplot(subdata,aes(x=redund,y=nfrac))+geom_point()+facet_wrap(~setnums,ncol=3,scales="free")


#err occurence
errs1 <- get_err_occurence(0, '_eur')
errs2 <- get_err_occurence(0, '_usa')
errs3 <- get_err_occurence(0, '_chn')
out2 <- rbind(errs1,errs2,errs3)
p <- ggplot(out2,aes(x=cover,y=mean,color=flag)) + geom_point()+scale_x_continuous(lim=c(50,500))+scale_y_continuous(lim=c(0,2.5))+geom_smooth(method=lm,color='black',aes(group=flag))
p <- p + scale_color_manual(values=brewer.pal(3, "Set1"),name="Страны",
                      labels=c("Китай", "Дания", "США"))
p <- p + theme(axis.title.x = element_text(size=15),axis.text.x  = element_text(size=15),axis.text.y  = element_text(size=15),axis.title.y = element_text(size=15))+xlab(label='n, покрытие')
p <- p + ylab(label='M(n), ожидаемое кол-во ошибок')
p <- p + theme(legend.text = element_text(size = 13))
p <- p + theme(legend.title = element_text(size = 13,face="bold"))


p#mds
temp_mds2 <- load_all_sets(setnums=c(21,28),flag='alt_3',treshold=0.1,group='CHN',all=TRUE)
ggplot(temp_mds2,aes(x=V1,y=V2,color=country))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")
temp_mds_4 <- load_all_sets(setnums=c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,47,48,49,56,57,58,59,62,68,73,74,75,82,85,87,89,96,98,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122),flag='alt_4',treshold=0.3,group='CHN',all=TRUE)
data.p_vals.alt_4<- get_p_vals(setnums=c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,47,48,49,56,57,58,59,62,68,73,74,75,82,85,87,89,96,98,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122),'alt_4',0.3)
temp_mds_3 <- load_all_sets(setnums=c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,46,47,48,49,53,54,55,56,57,58,59,62,68,69,73,74,75,76,82,85,87,89,96,98,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122),flag='alt_3',treshold=0.3,group='CHN',all=TRUE)
data.p_vals.alt_4<- get_p_vals(setnums=c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,46,47,48,49,53,54,55,56,57,58,59,62,68,69,73,74,75,76,82,85,87,89,96,98,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122),'alt_4',0.3)
ggplot(temp_mds_4,aes(x=V1,y=V2,color=clusters))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")
ggplot(temp_mds_4,aes(x=V1,y=V2,color=country))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")
ggplot(temp_mds_4,aes(x=V1,y=V2,color=state))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")
ggplot(temp_mds_4,aes(x=V1,y=V2,color=gender))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")

diabete_sep_4 <- load_all_sets(setnums=c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,47,48,49,56,57,58,59,62,68,73,74,75,82,85,87,89,96,98,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122),flag='alt_4',treshold=0.25,group='CHN',all=FALSE)
ggplot(diabete_sep_4,aes(x=V1,y=V2,color=state))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")
gender_sep_4 <- load_all_sets(setnums=c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,47,48,49,56,57,58,59,62,68,73,74,75,82,85,87,89,96,98,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122),flag='alt_4',treshold=0.25,group='USA',all=FALSE)
ggplot(gender_sep_4,aes(x=V1,y=V2,color=gender))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")

ggplot(temp_mds_3,aes(x=V1,y=V2,color=clusters))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")
ggplot(temp_mds_3,aes(x=V1,y=V2,color=country))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")


sub.mds_4 <- temp_mds_4[temp_mds_4$setnum %in% c(6),]

ggplot(sub.mds_4,aes(x=V1,y=V2,color=country))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=4,scales="free")
ggplot(sub.mds_4,aes(x=V1,y=V2,color=clusters))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=4,scales="free")
agg.p_vals <- aggregate(data.p_vals.alt_4$p_values,by=list(data.p_vals.alt_4$set),sum)

###6 set
sub.mds_4 <- temp_mds_4[temp_mds_4$setnum %in% c(6),]
sub.mds_4$category <- '1'
sub.mds_4$category[sub.mds_4$V1 < 0] <- '2'
p2 <- ggplot(sub.mds_4,aes(x=V1,y=V2,color=category))+geom_point(size=4)#+scale_color_brewer(palette="Set1")#+facet_wrap(~setnum,ncol=4,scales="free")
p2 <- p2 + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(size=15),axis.text.y  = element_text(size=15),axis.title.y = element_text(size=20))
p2 <- p2 + theme(legend.text = element_text(size = 13),legend.title=element_text(size = 13))+scale_color_manual(values=brewer.pal(2, "Set1"),name="Кластеры:",labels=c("1", "2"))
temp2 <- sub.mds_4[,c("names","category")]
nums <- c(6)
temp <- data.cov_distrs[data.cov_distrs$setnums %in% nums,]
temp <- merge(temp,temp2,by.x=c("names"),by.y=c("names"),all.x=T,all.y=F)
p <- ggplot(temp,aes(x=redund,y=frac,color=category))+geom_point()+geom_line()+scale_color_brewer(palette="Set1")#+facet_wrap(~setnums,ncol=3,scales="free")
p <- p + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(size=15),axis.text.y  = element_text(size=15),axis.title.y = element_text(size=20))+xlab(label='С, среднее покрытие')
#p <- p + ylab(label=paste(expression('F'[0]),', покрытая часть генома',sep=''))
p <- p + ylab(label='F, покрытая часть генома')+ theme(legend.position="none")
multiplot(p,p2,cols=2)


ggplot(sub.mds_4,aes(x=V1,y=V2,color=clusters))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=4,scales="free")
agg.p_vals <- aggregate(data.p_vals.alt_4$p_values,by=list(data.p_vals.alt_4$set),sum)



sub.mds_4 <-temp_mds_4[temp_mds_4$setnum %in% c(121,40,43,89),]
sub.mds_4 <-temp_mds_4[temp_mds_4$setnum %in% c(121,40,43),]
sub.mds_4 <-temp_mds_4[temp_mds_4$setnum %in% c(121,43,89),]
sub.mds_4$set <- sub.mds_4$setnum
m <- merge(strain_meta.data, sub.mds_4,by.x=c("set"),by.y=c("set"),all.x=F,all.y=T)
p <- ggplot(m,aes(x=V1,y=V2,color=country))+geom_point(size=3)+facet_wrap(~short_descript,ncol=3,scales="free")
p <- p + scale_color_manual(values=brewer.pal(4, "Set1"),name="Страны",labels=c("Китай", "Дания", "Референс","США"))
p <- p + theme(legend.text = element_text(size = 15),legend.title=element_text(size = 15),strip.text.x = element_text(size=15))
p <- p +  theme(axis.text = element_text(size = 14,color="black"),strip.text.x = element_text(size=18))#,axis.text.x  = element_text(color="black", size=14))
p
m <- merge(strain_meta.data, sub.mds_4,by.x=c("set"),by.y=c("set"),all.x=F,all.y=T)
#countries:

p <- ggplot(m,aes(x=V1,y=V2,color=country))+geom_point(size=2.5)+facet_wrap(~setnum,ncol=4,scales="free")
p <- p + scale_color_manual(values=brewer.pal(5, "Set1"),name="Страны",labels=c("Китай", "Дания", "Референс","США","Россия"))
p <- p + theme(legend.text = element_text(size = 13),legend.title=element_text(size = 13),strip.text.x = element_text(size=12))
p +  theme(axis.text = element_text(size = 10),strip.text.x = element_text(size=12))
#clusters:

p <- ggplot(m,aes(x=V1,y=V2,color=clusters))+geom_point(size=2.5)+facet_wrap(~setnum,ncol=4,scales="free")
p <- p + scale_color_manual(values=brewer.pal(5, "Set1"),name="Кластеры",labels=c("1", "2", "3","4","5"))
p <- p + theme(legend.text = element_text(size = 13),legend.title=element_text(size = 13),strip.text.x = element_text(size=12))
p +  theme(axis.text = element_text(size = 10),strip.text.x = element_text(size=12))
 #principal components
M <- load_distances(19,"alt_3",0.3)
pca <- prcomp(M[[3]])
tt <- princomp(M[[3]])
unclass(loadings(tt))[2]
PCbiplot2(pca,x='PC1',y='PC2')


#p-vals of segregation
data.p_vals.alt_3.04 <- get_p_vals(setnums=c(0,1,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,36,37,38,39,40,41,42,43,44,46,47,49,53,54,55,56,57,58,59,62,68,69,73,74,75,76,82,85,87,89,96,98,100,102,104,106,107,108,110,112,113,114,115,116,117,118,119,120,121,122),'alt_4',0.3)
data.p_vals.alt_3.03 <- get_p_vals(setnums=c(115),'alt_3',0.5)
da
#ref_mean_dists
#res1 <- ref_mean_dist(c(0,1,2,6,7,9,10,34,48),'final_3',0.3)
res1 <- ref_mean_dist(c(34),'final_3',0.3)
p1 <- ggplot(res1[[1]],aes(x=cov,y=dist))+geom_point()+facet_wrap(~setnum,scales="free")
p1 <- p1 + theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(size=15,color="black"),axis.text.y  = element_text(size=15,color="black"),axis.title.y = element_text(size=18))+xlab(label= 'N(млн. нуклеотидов), количество покрытых позиций,\n удовлетворяющих условию (1)')
p1 <- p1 + ylab(label='f, количество SNP на единицу длины генома\n(частота полиморфизмов)')
p1 <- p1 + theme(legend.text = element_text(size = 13))
p1 <- p1 + theme(legend.title = element_text(size = 13,face="bold"))
p1 + scale_x_continuous(labels=fff)
#res2 <- ref_mean_dist(c(0,1,2,6,7,9,10,34,48),'alt_3',0.3)
res2 <- ref_mean_dist(c(34),'alt_3',0.3)
p2 <- ggplot(res2[[1]],aes(x=cov,y=dist))+geom_point()+facet_wrap(~setnum,scales="free")
p2 <- p2 + theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(size=15,color="black"),axis.text.y  = element_text(size=15,color="black"),axis.title.y = element_text(size=18))+xlab(label= 'N(млн. нуклеотидов), количество покрытых позиций,\n удовлетворяющих условию (2)')
p2 <- p2 + ylab(label='f, количество SNP на единицу длины генома\n(частота полиморфизмов)')
p2 <- p2 + theme(legend.text = element_text(size = 13))
p2 <- p2 + theme(legend.title = element_text(size = 13,face="bold"))
p2 + scale_x_continuous(labels=fff)

#распределение снипов в конкретном образце в зависимости от покрытие
temp <- data.snp_characteristics[data.snp_characteristics$setnum==2,]
temp <- temp[temp$len_frac>0.3,]
temp$quarts <- cut(temp$cover / temp$mean, breaks=seq(0,2,0.5))
temp$quartsfrac <- temp$cover / temp$mean
temp$flag <- 'less'
temp$flag[temp$redund>=4.522294] <- 'more'

ggplot(temp[which(!is.na(temp$quarts)),], aes(x=quarts,y=frac,fill=flag))+geom_boxplot() + scale_y_continuous(lim=c(0,0.05))
ggplot(temp[which(!is.na(temp$quarts)),], aes(x=quartsfrac,y=frac,color=flag))+geom_point() + scale_y_continuous(lim=c(0,0.05))
ggplot(temp[temp$snpcount>20,],aes(x=cover,y=frac,color=names))+geom_line()+scale_x_continuous(lim=c(0,100))

#
p <- ggplot(data.cov_distrs[data.cov_distrs$setnums==2 & data.cov_distrs$frac>0.3,],aes(x=redund))+geom_histogram()
p
median(data.cov_distrs$redund[data.cov_distrs$setnums==2 & data.cov_distrs$frac>0.3])
#Labeller

labeller <- function(var, value){
  value <- as.character(value)
  if (var=="tbins") { 
    value <- substitute(c[i] %in% B, list(B = value))
  }
  return(value)
}



#Пики гетерогенности
temp <- spec[spec$names=="SRS011239",]
temp$tbins <- cut(temp$tfrac,breaks=c(seq(0,2,0.5)))
p_het <- ggplot(temp[which(!is.na(temp$tbins)),],aes(x=frac))+geom_histogram(binwidth=0.01)+scale_x_continuous(lim=c(0.05,0.95))+facet_wrap(~tbins,scales="free")
p_het <- p_het + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(size=15),axis.text.y  = element_text(size=15),axis.title.y = element_text(size=20))
p_het <- p_het + ylab(label='Количество позиций')
p_het <- p_het + xlab(label=expression(frac(n[i]^major - n[i]^alt,c[i]) ~ ', разность частот поддержек вариантов')) + theme(strip.text.x = element_text(size=15))
p_het# + ggtitle(expression('title'^2',asd'))

#abundance_correlation
#112
aa <- temp_preproc_abund(17)
rr <- process_abund(abund.data.matrix,aa)
pp <- abund_corr_plots(rr)
pp

temp_preproc_abund <- function(setnum){
  temp <- temp_mds_3[which(temp_mds_3$setnum==setnum),]
  temp$country <- as.character(temp$country)
  temp <- temp[which(temp$country %in% c('CHN','EUR','USA')),]
  temp$country <- as.factor(temp$country)
  strainname <- as.character(strain_meta.data$short_descript[which(strain_meta.data$set==setnum)])
  temp <- temp[,c("short_names","country","clusters")]
  temp$clusters <- as.integer(temp$clusters)
  names(temp) <- c("names","country",strainname)
  return(temp)
}

process_abund <- function(datm, cluster_data){
  cluster_data <- aa
  samples <- cluster_data$names
  datm <- abund.data.matrix
  n <- dim(cluster_data)[2]
  cluster_data <- subset(cluster_data,select=c(3:n))
  strains <- names(cluster_data)
  inds <- which(colnames(datm) %in% strains)
  #residue <- datm[samples,-inds]
  residue <- datm[samples,]
  datm <- subset(datm[samples,],select=c(inds))
  nstrains <- length(strains)
  nsamples <- length(samples)
  nameslist <- list()
  for (i in 1:nstrains){
    print(i)
    clustnum <- max(cluster_data[,i],na.rm=T)
    tempdata <- matrix(NA,nrow=nsamples,ncol=clustnum)
    clustnames <- paste(strains[i],c(1:clustnum),sep='.')
    nameslist[[i]] <- c(strains[i],clustnames)
    colnames(tempdata) <- clustnames
    for (j in 1:nsamples){
      if (is.na(cluster_data[j,strains[i]])){
        tempdata[j,] <- rep(NA,clustnum)
      }
      else{
        tempdata[j,cluster_data[j,strains[i]]] <- datm[samples[j],strains[i]]
      }
    }
    if (i == 1){
      res <- tempdata
    }
    else{
      res <- cbind(res,tempdata)
    }
  }
  res <- cbind(res,residue)
  treshold <- 0
  res <- res[,-which(apply(res,2,sum) <= treshold)]
  return(list(res,nameslist))
}



abund_corr_plots <- function(proc_abund_res){
  library(reshape)
  library(ggplot2)
  res <- proc_abund_res[[1]]
  nameslist <- proc_abund_res[[2]]
  #
  #res[which(is.na(res))] <- 0
  corrm <- cor(res,method="spearman",use="pairwise.complete.obs")
  #corrm <- cor(res,method="pearson")
  len <- length(nameslist)
  for (i in 1:len){
    print(i)
    strainnames <- nameslist[[i]]
    temp <- corrm[strainnames,]
    inds <- which(apply(temp,2,function(x){return(any(abs(x)>0.5))}))
    temp <- temp[,inds]
    plottemp <- melt(temp)
    names(plottemp) <- c("strain1","strain2","corr_coeff")
    maxes <- apply(temp[-1,],2,function(i) {return(max(abs(i)))})
    #print(temp[,which(maxes>0.4 & maxes != 1.0)])
    p <- ggplot(plottemp, aes(x=strain2,y=strain1)) + geom_tile(aes(fill = corr_coeff,colour = "white")) + scale_fill_gradient2(low="red",mid="black",high="green",midpoint=0)#gradientn(values=seq(-1, 1, length = 3),colours = c("#FF0000","#000000","#000000","#00FF00"))
    plot(p)
    #print('--')
    #print(p)
    #print(p, vp=viewport(layout.pos.row=1,layout.pos.col=1:9))
    
  }
  return(corrm)
}


#abundance_histogram
#нужно нормальзовать!
strainnames <- as.character(strain_meta.data$short_descript)
shortnames <- as.character(data.sample_names2$short_names)
temp <- abund.data.matrix[c(shortnames),c(strainnames)]
temp2 <- melt(temp)
res <- aggregate(temp2$value,by=list(temp2$X2),sum)
levs <- as.character(res$Group.1)[order(res$x,decreasing=T)]
res$Group.1 <- factor(x=as.character(res$Group.1),levels=levs)