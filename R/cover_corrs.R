
library(gplots)
PATH <- '/data1/bio/runs-kovarskiy/metagenome_data/'
source('distribution_funcs.R')
source('cover_corr_funcs.R')
setnum <- 2
win_size <- 100
redund_treshold <- 5

cols_3 <- brewer.pal(3, "Set1")
cols_2 <- c("red","blue")
cols_many <- cols.gentleman(10,"Spectral")
cols_many <- brewer.pal(10, "Spectral")

fname <- paste(PATH,setnum,'/normalized_cov_prof_windows_',win_size,'.csv',sep='')
samples_fname <- paste(PATH,setnum,'/goods.list',sep='')
good_samples <- read.csv(samples_fname, sep='\t', header=F)[,1]
data <- read.csv(fname,sep=',',header=F)
names(data) <- good_samples
good_samples <- as.character(data.cov_distrs$names[which(data.cov_distrs$redund > redund_treshold & data.cov_distrs$setnums == setnum)])
cut_data <- data[,c(good_samples)]
countries <- data.sample_names2$country[match(good_samples,data.sample_names2$names)]
gender <- data.sample_names2$gender[match(good_samples,data.sample_names2$names)]
subjs <- as.factor(data.sample_names2$subj_id_new[match(good_samples,data.sample_names2$names)])
ind <- data.sample_names2$gender[match(good_samples,data.sample_names2$names)]
cols_country <- cols_3[countries]
cols_subjs <- cols_many[subjs]
cols_gender <- cols_2[gender]
dinds <- which(apply(cut_data,1,function(x){return(length(which(x>0)))})>100)

heatmap.2(cor(cut_data[dinds,],method="pearson"),col=cols.gentleman(500), trace='none', scale="none",ColSideColors=c(cols_subjs),RowSideColors=c(cols_country))
#heatmap.2(cor(cut_data,method="pearson"),col=cols.gentleman(500), trace='none', scale="none",ColSideColors=c(cols_subjs),RowSideColors=c(cols_gender))
ggplot(cut_data,aes(x=cut_data[,1],y=cut_data[,2])) + geom_point(alpha=0.1)
#multiplot(p1, p1, p1, p1, cols=2)

new_stat <- cov_filter(data=cut_data,treshold=0)
ggplot(new_stat,aes(x=count,y=mean))+geom_point()

inds <- which(apply(cut_data,1,function(x){return(all(x>0))}))
heatmap.2(cor(cut_data[inds,],method="pearson"),col=cols.gentleman(500), trace='none', scale="none",RowSideColors=my_cols)

p1 <- heatmap.2(cor(t(cut_data),method="pearson"),trace='none', scale="none",RowSideColors=my_cols)

mat <- t(as.matrix(cut_data))/data.cov_distrs$redund[which(data.cov_distrs$redund > redund_treshold & data.cov_distrs$setnums == setnum)]
mat <- t(mat)
colnames(mat) <- good_samples
cut_data2 <- as.data.frame(mat)
ggplot(cut_data2,aes(x=cut_data2[,1],y=cut_data2[,2])) + geom_point(alpha=0.1)
ggplot(cut_data2,aes(x=cut_data2[,1],y=cut_data2[,2])) + geom_point(alpha=0.1,size=4)+scale_x_continuous(lim=c(0,5))+scale_y_continuous(lim=c(0,5))

v1 <- apply(cut_data,1,mean)
v2 <- apply(cut_data,1,sd)
res <- data.frame(means=v1,sd=v2)

cut_data2 <- cut_data
cut_data2[cut_data2<=0.1] <- 0
heatmap.2(cor(cut_data2,method="spearman"))

fname <- paste(PATH,setnum,'/cov_prof_windows_',win_size,'.csv',sep='')
data <- read.csv(fname,sep=',',header=F)


# good color palette for heatmap
cols.gentleman <- function(ncol=500,palette="RdBu") {
  library(RColorBrewer)
  hmcol <- colorRampPalette(brewer.pal(10, palette))(ncol)
  return(rev(hmcol))
}

heatmap.2(g_onko, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),
          col=cols.gentleman(500) trace='none', scale="none" , margin=c(10,10),
          RowSideColors=my_cols )



###Heatmaps

library(gplots)
PATH <- '/data1/bio/runs-kovarskiy/metagenome_data/'
source('distribution_funcs.R')
source('cover_corr_funcs.R')
setnum <- 2
win_size <- 100
redund_treshold <- 5

cols_3 <- brewer.pal(3, "Set1")
cols_2 <- c("red","blue")
cols_many <- cols.gentleman(10,"Spectral")
cols_many <- brewer.pal(10, "Spectral")

#fname <- paste(PATH,setnum,'/normalized_cov_prof_windows_',win_size,'.csv',sep='')
#samples_fname <- paste(PATH,setnum,'/goods.list',sep='')
#good_samples <- read.csv(samples_fname, sep='\t', header=F)[,1]
#data <- read.csv(fname,sep=',',header=F)
#names(data) <- good_samples
#good_samples <- as.character(data.cov_distrs$names[which(data.cov_distrs$redund > redund_treshold & data.cov_distrs$setnums == setnum)])
#cut_data <- data[,c(good_samples)]
setnum <- 0
cut_data <- load_distances(setnum, suffix = 'alt_4', treshold = 0.3)[[3]]
good_samples <- rownames(cut_data)
countries <- data.sample_names2$country[match(good_samples,data.sample_names2$names)]
gender <- data.sample_names2$gender[match(good_samples,data.sample_names2$names)]
subjs <- as.factor(data.sample_names2$subj_id_new[match(good_samples,data.sample_names2$names)])
ind <- data.sample_names2$gender[match(good_samples,data.sample_names2$names)]
cols_country <- cols_3[countries]
cols_subjs <- cols_many[subjs]
cols_gender <- cols_2[gender]
#dinds <- which(apply(cut_data,1,function(x){return(length(which(x>0)))})>100)
heatmap.2(cut_data,col=cols.gentleman(500), trace='none', scale="none",ColSideColors=c(cols_subjs),RowSideColors=c(cols_country))