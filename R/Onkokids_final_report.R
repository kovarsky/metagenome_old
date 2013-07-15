source('onko_funcs.R')
library(MASS)
library(ggplot2)
PATH <- '/data1/bio/runs-kovarskiy/onkokids/snpdata/'
data.onko.strains <- read.csv('/data1/bio/runs-kovarskiy/onkokids/snpdata/strain_metadata.tsv',sep='\t',header=F)
names(data.onko.strains) <- c("descript","length","set","short_descript")
data.onko.cov_profiles <- load_cov_distrs(c(1:17))
data.onko.coverages <- cov_distrs_proc(c(1:17))
data.onko.meta <- read.csv('/data1/bio/runs-kovarskiy/onkokids/onkokids/data/onko_extra_meta.txt',sep='\t')[,c(1:2)]
names(data.onko.meta) <- c('names','patient')
data.onko.meta$short_names <- data.onko.meta$patient

#
data.onko.coverages.old <- data.onko.coverages
data.onko.strains.old <- data.onko.strains
#

temp <- load_all_sets(c(1:4),flag='alt_2',treshold=0.1)
data.onko.strains$setnum <- data.onko.strains$set

temp <- merge(temp,data.onko.strains[,c("setnum","short_descript")],by.x=c("setnum"),by.y=c("setnum"),all.x=T,all.y=F)
ggplot(temp,aes(x=V1,y=V2,color=patient,label=names))+geom_point(size=2,color="black")+geom_text()+facet_wrap(~setnum,scales="free")+scale_color_brewer(palette="Dark2")+ggtitle("MDS for multiple strains")
ggplot(temp,aes(x=V1,y=V2,color=patient,label=names))+geom_point(size=2,color="black",alpha=0.6)+geom_text(size=3.5,aes(color=patient))+facet_wrap(~short_descript,scales="free")+scale_color_brewer(palette="Dark2",name="Subjects")+ggtitle("MDS for multiple strains")

ggplot(temp[temp$setnum==17,],aes(x=V1,y=V2,color=patient,label=names))+geom_point(size=2,aes(color=patient))+geom_text(size=3,aes(y=V2+0.0002),color="black")+facet_wrap(~short_descript,scales="free")+scale_color_brewer(palette="Dark2",name="Subjects")+ggtitle("MDS map for Streptococcus thermophilus")

ggplot(temp[temp$setnum==7,],aes(x=V1,y=V2,color=patient,label=names))+geom_point(size=2,aes(color=patient))+geom_text(size=3,aes(y=V2+0.0002),color="black")+facet_wrap(~short_descript,scales="free")+scale_color_brewer(palette="Dark2",name="Subjects")+ggtitle("MDS map for Enterococcus faecium DO")

ggplot(temp[temp$setnum==11,],aes(x=V1,y=V2,color=patient,label=names))+geom_point(size=2,aes(color=patient))+geom_text(size=3,aes(y=V2+0.0002),color="black")+facet_wrap(~short_descript,scales="free")+scale_color_brewer(palette="Dark2",name="Subjects")+ggtitle("MDS map for Klebsiella")


##################################################################

source('onko_funcs.R')
library(MASS)
library(ggplot2)
PATH <- '/data1/bio/runs-kovarskiy/onkokids_fewbacs/'
data.onko.strains <- read.csv('/data1/bio/runs-kovarskiy/onkokids_fewbacs/strain_metadata.tsv',sep='\t',header=F)
names(data.onko.strains) <- c("descript","length","set","short_descript")
data.onko.cov_profiles <- load_cov_distrs(c(1:4))
data.onko.coverages <- cov_distrs_proc(c(1:4))
data.onko.meta <- read.csv('/data1/bio/runs-kovarskiy/onkokids/onkokids/data/onko_extra_meta.txt',sep='\t')[,c(1:2)]
names(data.onko.meta) <- c('names','patient')
data.onko.meta$short_names <- data.onko.meta$patient


temp <- data.onko.coverages.old[data.onko.coverages.old$set %in% c(7,10,11,17),]
m <- data.frame(old=c(7,10,11,17),new=c(1,2,3,4))
temp$setnums<-m$new[match(temp$setnums,m$old)]
temp$dataset <- 'old'
data.onko.coverages$dataset <- 'new'
data.onko.new_and_old <- rbind(temp[which(temp$names %in% data.onko.coverages$names),],data.onko.coverages)
#data.onko.new_and_old <- merge(temp[,c("names","setnums","frac","redund","nts")],data.onko.coverages[,c("names","setnums","frac","redund","nts")],by.x=c("names","setnums"),by.y=c("names","setnums"),all.x=F,all.y=T)
#names(data.onko.new_and_old) <- c("names","setnums","frac.old","redund.old","nts.old","frac.new","redund.new","nts.new")
data.onko.strains$setnums <- data.onko.strains$set
data.onko.new_and_old <- merge(data.onko.new_and_old,data.onko.strains[,c("setnums","short_descript")],by.x=c("setnums"),by.y=c("setnums"))
ggplot(data.onko.new_and_old, aes(x=names, y=redund,fill=dataset)) + geom_histogram(binwidth=.5, position="dodge",stat="identity")+facet_wrap(~short_descript)+scale_y_log10()+scale_fill_brewer(palette="Set1")
ggplot(data.onko.new_and_old, aes(x=names, y=frac,fill=dataset)) + geom_histogram(binwidth=.5, position="dodge",stat="identity")+facet_wrap(~short_descript)+scale_fill_brewer(palette="Set1")
ggplot(data.onko.new_and_old,aes(y=frac,x=redund,color=dataset))+geom_point()+geom_line(group=names)+facet_wrap(~short_descript)+scale_color_brewer(palette="Set1")


PATH <- '/data1/bio/runs-kovarskiy/onkokids_fewbacs/'
temp <- load_all_sets(c(1:4),flag='alt_3',treshold=0.1)
ggplot(temp,aes(x=V1,y=V2,color=patient,label=names))+geom_point(size=3,aes(color=patient))+geom_text(size=3)+facet_wrap(~setnum,scales="free")+scale_color_brewer(palette="Dark2")+ggtitle("MDS for multiple strains")
data.onko.strains$setnum <- data.onko.strains$setnums
temp2 <- merge(temp,data.onko.strains,by.x=c("setnum"),by.y=c("setnum"),all.x=T,all.y=T)
ggplot(temp2,aes(x=V1,y=V2,color=patient,label=names))+geom_point(size=3,aes(color=patient))+geom_text(size=3)+facet_wrap(~descript,scales="free")+scale_color_brewer(palette="Dark2")+ggtitle("MDS for multiple strains")
#ggplot(temp,aes(x=V1,y=V2,color=patient,label=names))+geom_point(size=2,color="black",alpha=0.6)+geom_text(size=3.5,aes(color=patient))+facet_wrap(~short_descript,scales="free")+scale_color_brewer(palette="Dark2",name="Subjects")+ggtitle("MDS for multiple strains")

ggplot(temp[temp$setnum==1,],aes(x=V1,y=V2,color=patient,label=names))+geom_point(size=2,aes(color=patient))+geom_text(size=3,aes(y=V2+0.0002),color="black")+facet_wrap(~short_descript,scales="free")+scale_color_brewer(palette="Dark2",name="Subjects")+ggtitle("MDS map for Streptococcus thermophilus")

ggplot(temp[temp$setnum==7,],aes(x=V1,y=V2,color=patient,label=names))+geom_point(size=2,aes(color=patient))+geom_text(size=3,aes(y=V2+0.0002),color="black")+facet_wrap(~short_descript,scales="free")+scale_color_brewer(palette="Dark2",name="Subjects")+ggtitle("MDS map for Enterococcus faecium DO")

ggplot(temp[temp$setnum==11,],aes(x=V1,y=V2,color=patient,label=names))+geom_point(size=2,aes(color=patient))+geom_text(size=3,aes(y=V2+0.0002),color="black")+facet_wrap(~short_descript,scales="free")+scale_color_brewer(palette="Dark2",name="Subjects")+ggtitle("MDS map for Klebsiella")
