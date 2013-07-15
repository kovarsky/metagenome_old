#with correction
data.snp_characteristics <- compute_snp_rates(c(2))
temp <- data.snp_characteristics#[data.snp_characteristics$setnum==2,]
inds<-which(temp$cover>2)
temp3 <- temp[inds,]
#with correction
temp3 <- temp3[temp3$mean<=temp3$cover,]
temp2 <- aggregate(temp3[,c("count","snpcount")],by=list(temp3$names,temp3$setnum,temp3$len_frac,temp3$redund),sum)
names(temp2) <- c("names","setnum","len_frac","redund", "count","snpcount")
temp2$frac <- temp2$snpcount/temp2$count
temp2 <- merge(temp2,data.sample_names2,by.x=c("names"),by.y=c("names"),all.x=T,all.y=T)
ggplot(temp2[temp2$len_frac>0.2 & temp2$setnum %in% c(0,1,2,3,6,7,8,12,9,24,18,34,43,31,32),],aes(x=len_frac,y=frac))+geom_point(size=2)+facet_wrap(~setnum,ncol=5,scales="free")+scale_color_brewer(palette="Set1")
f_data <- temp2[temp2$len_frac>0.2,]

temp$frac <- temp$snpcount/temp$count
inds<-which(temp$cover>2)
temp3 <- temp[inds,]
bb <- sort(unique(temp3$redund),decreasing=T)
goods <- bb[c(1,10,20,30,40,50,60)]
ggplot(temp3[temp3$redund%in%goods & temp3$setnum %in% c(0,1,2,3,6,7,8,12,9,24,18,34,43,31,32),],aes(x=cover,y=frac,color=names))+geom_line()+scale_x_continuous(lim=c(0,100))

#without correction
temp3 <- temp[inds,]
temp2 <- aggregate(temp3[,c("count","snpcount")],by=list(temp3$names,temp3$setnum,temp3$len_frac,temp3$redund),sum)
names(temp2) <- c("names","setnum","len_frac","redund", "count","snpcount")
temp2$frac <- temp2$snpcount/temp2$count
temp2 <- merge(temp2,data.sample_names2,by.x=c("names"),by.y=c("names"),all.x=T,all.y=T)
ggplot(temp2[temp2$len_frac>0.2,],aes(x=len_frac,y=frac))+geom_point(size=4)+facet_wrap(~setnum,ncol=8,scales="free")+scale_color_brewer(palette="Set1")
f_data <- temp2[temp2$len_frac>0.2,]




#Кластеризация
data <- load_distances(1,"alt_2",0.2,all=F,group="CHN",drop_outlier=TRUE)
data.mds <- make_data_frame_mds(data[[3]])
ggplot(data.mds,aes(x=V1,y=V2,color=country))+ geom_point() + scale_color_brewer(palette="Set1")

source('distribution_funcs.R')
source('new_dist_proc_funcs.R')
#alt_4
temp_mds <- load_all_sets(setnums=c(0,1,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,36,37,38,39,40,41,42,43,44,47,56,57,58,59,62,68,73,74,75,82,85,87,89,96,98,100,102,103,104,106,107,108,109,110,112,113,114,115,116,117,118,119,120,121,122),flag='alt_4',treshold=0.3,group='CHN',all=TRUE)
data.p_vals.alt_4.03 <- get_p_vals(setnums=c(0,1,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,36,37,38,39,40,41,42,43,44,47,56,57,58,59,62,68,73,74,75,82,85,87,89,96,98,100,102,103,104,106,107,108,109,110,112,113,114,115,116,117,118,119,120,121,122),'alt_4',0.3)
data.mds.alt_4.03 <- temp
temp <- data.mds.alt_3.01
temp <- load_all_sets(setnums=c(0,1,2),flag='alt_4',treshold=0.3)
ggplot(temp_mds,aes(x=V1,y=V2,color=clusters))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")
ggplot(temp_mds,aes(x=V1,y=V2,color=country))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")
###
M <- load_distances(19,"alt_3",0.3)
pca <- prcomp(M[[3]])
tt <- princomp(M[[3]])
unclass(loadings(tt))[2]
PCbiplot2(pca,x='PC1',y='PC2')
###

temp <- load_all_sets(setnums=c(9,19,32,40,41,42),flag='alt_3',treshold=0.1,group='CHN',all=TRUE)
p1 <- ggplot(temp,aes(x=V1,y=V2,color=country))+geom_point(size=1.5)+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")
p1

#CHINESE FOOD
temp <- load_all_sets(setnums=c(0,1,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,36,37,38,39,40,41,42,43,44,47,56,57,58,59,62,68,73,74,75,82,85,87,89,96,98,100,102,103,104,106,107,108,109,110,112,113,114,115,116,117,118,119,120,121,122),flag='alt_4',treshold=0.3,group='CHN',all=FALSE)
ggplot(temp,aes(x=V1,y=V2,color=country))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")

temp <- load_all_sets(setnums=c(0,1,4,6),flag='alt_4',treshold=0.3,group='CHN',all=FALSE)
#USAAAA

temp <- load_all_sets(setnums=c(0,1,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,36,37,38,39,40,41,42,43,44,47,56,57,58,59,62,68,73,74,75,82,85,87,89,96,98,100,102,103,104,106,107,108,109,110,112,113,114,115,116,117,118,119,120,121,122),flag='alt_4',treshold=0.2,group='USA',all=FALSE)
ggplot(temp,aes(x=V1,y=V2,color=gender))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")
ggplot(temp,aes(x=V1,y=V2,color=factor(subj_id)))+geom_point()+facet_wrap(~setnum,ncol=8,scales="free")






#alt_3
temp <- load_all_sets(setnums=c(0,1,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,36,37,38,39,40,41,42,43,44,46,47,49,53,54,55,56,57,58,59,62,68,69,73,74,75,76,82,85,87,89,96,98,100,102,104,106,107,108,110,112,113,114,115,116,117,118,119,120,121,122),flag='alt_3',treshold=0.3)
data.p_vals.alt_3.03 <- get_p_vals(setnums=c(0,1,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,36,37,38,39,40,41,42,43,44,46,47,49,53,54,55,56,57,58,59,62,68,69,73,74,75,76,82,85,87,89,96,98,100,102,104,106,107,108,110,112,113,114,115,116,117,118,119,120,121,122),'alt_3',0.3)
data.mds.alt_3.03 <- temp

ggplot(temp,aes(x=V1,y=V2,color=clusters))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=9,scales="free")
ggplot(temp,aes(x=V1,y=V2,color=country))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=9,scales="free")



temp2 <- load_all_sets(setnums=c(1,6,13,14,15,16,17,18,20,21,22,23,24,25,28,29,30,32,36,37,38,39,42,43,44,49,55,56,57,58,68,69,73,74,75,76,82,85,87,89,98,100,102,104),flag='final_2',treshold=0.1)
ggplot(temp2,aes(x=V1,y=V2,color=clusters))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")

#Характеристика снипов 12-го штамма
data.snp_characteristics <- compute_snp_rates(c(12))
temp <- data.snp_characteristics#[data.snp_characteristics$setnum==2,]
temp2 <- temp[temp$setnums==12,]
temp2$flag <- 'less'
temp2$flag[which(temp2$mean>temp2$cover)] <- 'more'
temp3 <- temp2[which(temp2$snpcount>10),]
temp3 <- temp3[temp3$cover<=100,]

ggplot(temp3,aes(x=cover,y=frac,color=flag)) + geom_point(alpha=0.3)+scale_y_continuous(lim=c(0,0.1))+scale_x_continuous(lim=c(0,100))+scale_color_brewer(palette='Set1')

v1 <- dpois(x=1:100,lambda=0.05)
poiss_d <- data.frame(frac=v1,cover=1:100)

errs1 <- get_err_occurence(12, '_eur')
errs2 <- get_err_occurence(12, '_usa')
errs3 <- get_err_occurence(12, '_chn')
out <- rbind(errs1,errs2,errs3)

temp <- data.cov_distrs[which(data.cov_distrs$setnums %in% c(2)),]
p2 <- ggplot(temp,aes(x=redund,y=frac))+geom_point()+facet_wrap(~setnums,ncol=6,scales="free")
ggplot(temp,aes(x=redund,y=frac))+geom_point()+facet_wrap(~setnums,ncol=6,scales="free")

a <- read.csv('/data1/bio/runs-kovarskiy/metagenome_data/2/normalized_mutual_profile_test__2_1',sep='\t',header=F)
ggplot(temp,aes(x=redund,y=frac))+geom_point()+facet_wrap(~setnums,ncol=8,scales="free")



abc <- load_snp_read_counts(c(0,1,2),treshold=5)
ggplot(abc,aes(x=cover,y=frac))+geom_point(alpha=0.8)+scale_y_continuous(lim=c(0,0.2))+facet_grid(~setnum,scale="free")


a <- load_snp_counts(c(0,1,2,3,6,7,8,9,10,11,12,13,14,15,16,17),treshold=5)


