library(RGraphics) # support of the "R graphics" book, on CRAN
library(gridExtra) 
g1 <- tableGrob(head(iris))
string <- "
This famous (Fisher's or Anderson's) iris data set gives the
measurements in centimeters of the variables sepal length and width
and petal length and width, respectively, for 50 flowers from each of
3 species of iris. The species are Iris setosa, versicolor, and
virginica.
"
g2 <- splitTextGrob(string)
#"Sepal.Length" "Sepal.Width"  "Petal.Length" "Petal.Width"  "Species"
g3 <- qplot(Sepal.Length,  Petal.Length, data=iris, colour=Species)
grid.arrange(g1, g2, g3, ncol=1, main="The iris data")





strain_plots <- function(setnum){
  descript <- paste(strain_meta.data$descript[strain_meta.data$set==setnum])
  sampletable <- count(as.character(keke$country))
  names(sampletable) <- c("country", "num")
  g0 <- tableGrob(sampletable)
  grid.newpage()
  h <- grobHeight(g0)
  w <- grobWidth(g0)
  g0_title <- textGrob("Количество образцов", y=unit(0.5,"npc") + 0.5*h, vjust=0, gp=gpar(fontsize=16))
  g0 <- gTree(children=gList(g0, g0_title))
  g1 <- ggplot(temp_mds_4[temp_mds_4$setnum==setnum,],aes(x=V1,y=V2,color=clusters))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")+ggtitle(label="MDS-проекция матрицы расстояний\n, цветом обозначены страны")
  g2 <- ggplot(temp_mds_4[temp_mds_4$setnum==setnum,],aes(x=V1,y=V2,color=country))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")+ggtitle(label="MDS-проекция матрицы расстояний\n, цветом обозначены детектируемые кластеры")
  g3 <- ggplot(data.cov_distrs[data.cov_distrs$setnums==setnum,],aes(x=redund,y=frac))+geom_point()+facet_wrap(~setnums,ncol=3,scales="free")
  g3 <- g3 + theme(axis.title.x = element_text(size=12),axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12),axis.title.y = element_text(size=12))+xlab(label='С, среднее покрытие')
  #p <- p + ylab(label=paste(expression('F'[0]),', покрытая часть генома',sep=''))
  g3 <- g3 + ylab(label='F, покрытая часть генома')+ggtitle("Кривые насыщения покрытия")
  g4 <- ggplot(data.var_measure,aes(x=length,y=mean_dist,color=factor(tresh)))+
    geom_line()+geom_point(data=data.var_measure[data.var_measure$setnum==setnum,],aes(x=length,y=mean_dist),color="red",size=3)+theme(legend.text = element_text(size = 10),legend.title=element_text(size = 10))+
    scale_color_manual(values=brewer.pal(2, "Paired"),name="Порог на покрытую\n часть генома,\nпо которым\n вычисляется\nизменчивость:",labels=c("0.3", "0.6"))+
    ylab(label='Среднее расстояние между образцами\nдля данного вида')+
    xlab(label='Размер генома')+ggtitle(label="Изменчивость бактериальных видов,\nкрасной точкой показан рассматриваемый вид")
  #--------------------------------
  res <- ref_mean_dist(c(setnum),'alt_4',0.3)
  g5 <- ggplot(res[[1]],aes(x=cov,y=dist))+geom_point()+facet_wrap(~setnum,scales="free")
  g5 <- g5 + theme(axis.title.x = element_text(size=12),axis.text.x  = element_text(size=12,color="black"),axis.text.y  = element_text(size=12,color="black"),axis.title.y = element_text(size=12))+xlab(label= 'N(млн. нуклеотидов), количество покрытых позиций,\n удовлетворяющих условию')
  g5 <- g5 + ylab(label='f, количество SNP на единицу длины генома\n(частота полиморфизмов)')
  g5 <- g5 + theme(legend.text = element_text(size = 12))
  g5 <- g5 + theme(legend.title = element_text(size = 12,face="bold"))
  g5 <- g5 + scale_x_continuous(labels=fff)
  g6 <- tableGrob(data.p_vals.alt_4[data.p_vals.alt_4$set == setnum, c("groups","p_values")])
  temp <- data.p_vals.alt_4.duals[data.p_vals.alt_4.duals$set==setnum,]
  res <- aggregate(temp$p_values,by=list(temp$macro),FUN=sum)
  res <- res[res$x<=0.01,]$Group.1
  if (length(res)>0){
    diffstring <- paste("Значимые отличия в группах:", res)
  }
  else{
    diffstring <- "Значимые отличия в группах не найдены" 
  }
  grid.newpage()
  h <- grobHeight(g6)
  w <- grobWidth(g6)
  g6_footnote <- textGrob(diffstring, x=unit(0.5,"npc") - 0.5*w, y=unit(0.5,"npc") - 0.5*h, vjust=1, hjust=0,gp=gpar( fontface="italic"))
  g6 <- gTree(children=gList(g6, g6_footnote))
  png("/data1/bio/runs-kovarskiy/metagenome_data/1.png",width=8.27,height=8.27,units="in",res=600)
  grid.arrange(g2,g1,g0,g3, g4, g5, g6, ncol=2, main=descript)
  dev.off()
  #grid.arrange(grob1, g7, ncol=2, main="Akkermansia")
  
}


data.p_vals.alt_4.duals <- get_dual(data.p_vals.alt_4)





#temp_mds_4 <- load_all_sets(setnums=c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,47,48,49,56,57,58,59,62,68,73,74,75,82,85,87,89,96,98,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122),flag='alt_4',treshold=0.3,group='CHN',all=TRUE)
#data.p_vals.alt_4<- get_p_vals(setnums=c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,47,48,49,56,57,58,59,62,68,73,74,75,82,85,87,89,96,98,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122),'alt_4',0.3)
#temp_mds_3 <- load_all_sets(setnums=c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,46,47,48,49,53,54,55,56,57,58,59,62,68,69,73,74,75,76,82,85,87,89,96,98,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122),flag='alt_3',treshold=0.3,group='CHN',all=TRUE)
#data.p_vals.alt_4<- get_p_vals(setnums=c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,46,47,48,49,53,54,55,56,57,58,59,62,68,69,73,74,75,76,82,85,87,89,96,98,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122),'alt_4',0.3)


#ggplot(temp_mds_4,aes(x=V1,y=V2,color=state))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")
#ggplot(keke,aes(x=V1,y=V2,color=country))+geom_point()+scale_color_brewer(palette="Set1")+facet_wrap(~setnum,ncol=8,scales="free")

# 
# treshold <- 0.6
# setnums <- c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,47,48,49,56,57,58,59,62,68,73,74,75,82,85,87,89,96,98,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122)
# res <- data.frame(setnum=c(),mean_dist=c(),nums=c())
# for (setnum in setnums){
#   dist_m <- load_distances(setnum,'alt_4',treshold,group='CHN',all=TRUE)[[3]]
#   if (dim(dist_m)[1]>=10){
#     temp <- data.frame(setnum=c(setnum),mean_dist=c(mean(dist_m[lower.tri(dist_m,diag=F)])),nums=c(dim(dist_m)[1]))
#     res <- rbind(res,temp)
#   }
# }
# 
# res2 <- aggregate(data.cov_distrs$nts,by=list(data.cov_distrs$setnums),sum)
# names(res2) <- c("setnum","total_nts")
# m <- merge(res,res2,by.x=c("setnum"),by.y=c("setnum"),all.x=T,all.y=F)
# m$set <- m$setnum
# m <- merge(m,strain_meta.data,by.x=c("set"),by.y=c("set"),all.x=T,all.y=F)
# 
# #ggplot(m,aes(x=length,y=mean_dist))+geom_line()+geom_point(data=m[m$setnum==2,],aes(x=length,y=mean_dist),color="red",size=5)
# #ggplot(m,aes(x=total_nts,y=mean_dist))+geom_line()+geom_point(data=m[m$setnum==2,],aes(x=total_nts,y=mean_dist),color="red",size=5)
# ###
# 
# m$tresh <- treshold
# data.var_measure <- rbind(m,data.var_measure)
