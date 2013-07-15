PATH <- '/data1/bio/runs-kovarskiy/metagenome_data/'
source('distribution_funcs.R')
setnum <- 3
cov_histograms <- load_cov_distrs(setnum)
a <- load_hetero_distrs(setnum)
b <- load_cov_distrs(setnum)
m <- merge(a,b,c("cover","names"),c("cover","names"),all.x=T,all.y=T)
m$hetero[which(is.na(m$hetero))] <- 0
m$hetcount[which(is.na(m$hetcount))] <- 0
m$mean_het <- m$hetero/m$count
m$het_site_rate <- m$hetcount/m$count
temp <- data.cov_distrs[data.cov_distrs$setnums==setnum & data.cov_distrs$redund>=100,]
tempnames <- temp$names
res <- m[m$names%in%tempnames,]
res <- merge(res,temp[,c("names","redund")],by.x=c("names"),by.y=c("names"),all.x=T,all.y=F)
res$ncover <- res$cover/res$redund
res$bins <- cut(res$ncover,breaks=c(seq(0,3,0.1)))

agg_res <- aggregate(x=res[,c("hetero","hetcount","count")],by=list(res$bins,res$names),sum)
names(agg_res) <- c("intervals","names","hetero","hetcount","count")
agg_res$nhetero <- agg_res$hetero/agg_res$count
agg_res$nhetcount <- agg_res$hetcount/agg_res$count
ggplot(agg_res,aes(x=intervals,y=nhetero,color=names))+geom_point()
ggplot(agg_res,aes(x=intervals,y=nhetcount,color=names))+geom_point()
###
# cov_hists <- cov_histograms[cov_histograms$names%in%temp$names,]
# cov_hists <- merge(cov_hists,temp[,c("names","redund")],by.x=c("names"),by.y=c("names"),all.x=T,all.y=F)
# cov_hists$tfrac <- cov_hists$cover/cov_hists$redund
# cov_hists$cover_bins <- cut(cov_hists$tfrac,breaks=c(seq(0,4,0.1)))
# cov_hists <- aggregate(cov_hists$count,by=list(cov_hists$names,cov_hists$cover_bins),sum)
# names(cov_hists) <- c("names","cover_bins","count")
spec <- load_hetero_spectras(covdata=temp)
spec$buff <-1
agg_spec <- aggregate(spec$buff,by=list(spec$names,spec$bins,spec$tbins),sum)
names(agg_spec) <- c("names","frac_bins","cover_bins","hetcount")
# agg_spec <- merge(cov_hists,agg_spec,by.x=c("names","cover_bins"),by.y=c("names","cover_bins"),all.x=T,all.y=T)
# agg_spec$hetfrac <- agg_spec$hetcount/agg_spec$count
ggplot(agg_spec[agg_spec$frac_bins!="(0,0.1]",],aes(x=frac_bins,y=hetfrac,color=cover_bins))+facet_wrap(~names,scales="free")+geom_point()

ggplot(spec[spec$total>50,],aes(x=frac))+geom_histogram(binwidth=0.01)+scale_x_continuous(lim=c(0.05,0.95))+facet_wrap(~names)
###
temp <- spec[spec$names=="SRS011239",]
temp$tbins <- cut(temp$tfrac,breaks=c(seq(0,2,0.5)))
ggplot(temp[which(!is.na(temp$tbins)),],aes(x=frac))+geom_histogram(binwidth=0.01)+scale_x_continuous(lim=c(0.05,0.95))+facet_wrap(~tbins,scales="free")

ggplot(agg_res[agg_res$names=="SRS011239",],aes(x=intervals,y=nhetcount,color=names))+geom_point()
hetname <- paste(PATH,setnum,"SRS011239.het",sep='/')
data <- read.csv(hetname,sep='\t',header=F)
names(data) <- c('eq','max','total')
data$max_frac <- data$max/data$total
inds <- which(data$eq==data$max)
data$max_frac[inds] <- 1 - data$max_frac[inds]
data$bins <- cut(data$max_frac,breaks=c(seq(0,1,0.1)))
ggplot(data[data$total>100,],aes(x=max_frac))+geom_histogram(binwidth=0.01)+scale_x_continuous(lim=c(0.1,0.9))
ggplot(data,aes(x=total,color=bins))+geom_density()
###


ggplot(m[m$names%in%tempnames,],aes(x=cover,y=mean_het,color=names))+stat_smooth()+scale_x_continuous(lim=c(5,200))+scale_y_continuous(lim=c(0,0.02))
ggplot(m[m$names%in%tempnames,],aes(x=cover,y=mean_het,color=names))+geom_line()+scale_x_continuous(lim=c(5,400))+scale_y_continuous(lim=c(0,0.02))
ggplot(m[m$names%in%tempnames,],aes(x=cover,y=het_site_rate,color=names))+geom_line()+scale_x_continuous(lim=c(5,400))+scale_y_continuous(lim=c(0,0.02))
ggplot(m[m$names%in%tempnames,],aes(x=cover,y=count,color=names))+geom_line()+scale_x_continuous(lim=c(5,600))+scale_y_continuous(lim=c(0,100000))
#total:
data.cov_distrs[data.cov_distrs$setnums==setnum & data.cov_distrs$redund>=200,]

fname <- paste(PATH,setnum,'sumheterogenityalt',sep='/')
total_het <- read.csv(fname,sep='\t',header=F)
names(total_het) <- c('cover','heterogenity','count')
total_het$diversity <- total_het$heterogenity / total_het$count
ggplot(total_het,aes(x=cover,y=diversity))+geom_line()
ggplot(total_het,aes(x=cover,y=diversity))+geom_point(alpha=0.3)+stat_smooth()+scale_x_continuous(lim=c(0,5000))+scale_y_continuous(lim=c(0,0.05))