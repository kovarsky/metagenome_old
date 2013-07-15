json_open <- function(filename){
  print('json_open...')
  require('rjson')
  file_as_text <- readLines(filename)
  json_data <- fromJSON(file_as_text)
  #str(json_data)
  return(json_data)
}
list_to_matrix <- function(list){
  xdim <- length(list[[1]])
  ydim <- length(list)
  res <- array(unlist(list), dim = c(xdim, xdim, ydim))
  return(res)
}
give_names <- function(datastruct, names){
  rownames(datastruct) <- names
  colnames(datastruct) <- names
  return(datastruct)
}

data_proc <- function(list, names,summ=F,contignum=1){
  print('data_proc...')
  data <- list_to_matrix(list)
  if (summ){
  res <- apply(data,c(1,2),sum)
  }
  else{
  res <- data[,,contignum]
  }
  res <- give_names(res, names)
  return(res)
}


json_distance <- function(filename, treshnum,summ=F,contignum=1){
  print('json_distance...')
  # !!! It returns matrices from JSON with mututal covered variable positions
  # and difference in these position
  #function that allows to obtain whole genome distances
  #i.e. it generalizes contig information
  json_data <- json_open(filename)
  names <- json_data[[1]]
  data.heterogen <- data_proc(json_data[[2]], names,summ,contignum)
  #data.heterogen <- c()
  data.mut_cover <- data_proc(json_data[[3]][[treshnum]], names,summ,contignum)
  data.diff <- data_proc(json_data[[4]][[treshnum]], names,summ,contignum)
  return(list(names=names, heterogen=data.heterogen, mut_cover=data.mut_cover, diff=data.diff))
}

stat_routine <- function(json_data,names,lims,contignum){
  #print(contignum)
  print(lims)
  print('stat_routine...')
  cover <- array(unlist(json_data), lims)
  #cover <- apply(cover,c(1,3),summ)
  cover <- as.data.frame(cover[,contignum,])
  print(dim(cover))
  print(length(names))
  cover$samples <- names
  return(cover)
}

get_data_cover <- function(setnum){
  dirr <- '/data1/bio/runs-kovarskiy/metagenome_data'
  dist_filename <- paste(dirr, setnum, 'distance.json',sep='/')
  stat_filename <- paste(dirr, setnum, 'stat.json',sep='/')
  data.distance <- json_distance(dist_filename,summ=F,treshnum=1,contignum=1)
  data.stats <- json_stat(stat_filename, data.distance$names)
  data.cover <- data.frame(cover=data.stats$cover[,1],samples = data.stats$cover$samples,stringsAsFactors=F)
  return(data.cover)
}

get_all_densities <- function(setnums){
  len <- length(setnums)
  for (i in 1:len){
    setnum <- setnums[i]
    temp.cover <- get_data_cover(setnum)
    temp.density <- density(temp.cover$cover)
    temp.density <- data.frame(x=temp.density$x, y=temp.density$y)
    temp.density$set <- setnum
    temp.cover$set <- setnum
    temp.cover <- merge(temp.cover,strain_meta.data[,c(-1)],by.x=c("set"),by.y=c("set"),all.x=T,all.y=F)
    if (i==1){
      res.density <- temp.density
      res.cover <- temp.cover
    }
    else{
      res.cover <- rbind(res.cover,temp.cover)
      res.density <- rbind(res.density, temp.density)
    }
  }
  return(list(res.cover,res.density))
}

json_stat <- function(filename, names){
  print('json_stat...')
  print(filename)
  ###Pretty much similar with json_distance function
  # It just returns vectors with with coverage and self-heterogenity
  json_data <- json_open(filename)
  #treshnum <- length(json_data[[1]])
  contignum <- length(json_data[[1]][[1]])
  samplnum <- length(json_data[[1]][[1]][[1]])
  lims <- c(samplnum,contignum,1)
  
  cover <- stat_routine(json_data[[1]][[1]],names,lims)
  lims[3] <- 1
  heterogen <- stat_routine(json_data[[2]],names,lims)
  return(list(cover=cover, heterogen=heterogen))
}

drop_by_coverage <- function(data.cover,distance_m,treshold=10000){
  max_cov <- max(data.cover$cover)
#   require(pastecs)
#   ddensity <- density(data.cover$cover)
#   ts_cover <- ts(ddensity$y)
#   t_points <- turnpoints(ts_cover)
#   tresh_p <- which(t_points$pits)[1]
#   treshold <- ddensity$x[tresh_p]
  #print(data.cover$cover)
  #treshold <- 1/4
  treshold <- max_cov/4
  print(treshold)
#   if (length(treshold)!=0 & !is.na(treshold)){
#     data.cover <- data.cover[data.cover$cover>treshold,]
#     distance_m <- distance_m[data.cover$samples,data.cover$samples]
#   }
  data.cover <- data.cover[data.cover$cover>treshold,]
  distance_m <- distance_m[data.cover$samples,data.cover$samples]
  return(list(distance_m,data.cover))
}

drop_by_abundance <- function(abund_info, distance_m){
  abnames <- abund_info$names[abund_info$abund!=0]
  abnames <- which(colnames(distance_m)%in%abnames)
  distance_m <- distance_m[abnames,abnames]
  abund_info <- abund_info[abund_info$names%in%abnames,]
  return(distance_m)
}

distances_histograms <- function(dir,maxnum){
  for (i in 1:maxnum){
    distance_histogram(dir,i,T)
  }
}



against_hist <- function(dm,country,contignum,dir){
  print(country)
  dm$temp<-country
  dm$temp[dm$country2!=country]<-'DIFF'
  plot <- ggplot(dm[dm$country1=='CHN',],aes(x=distance,color=temp))+geom_histogram(binwidth=0.001)+facet_grid(temp ~ .)
  img_path <- paste(dir,country,'_against_hist',contignum,'.png',sep='')
  ggsave(img_path,plot,width=9,height=7)
}


json_general <- function(dist_filename, stat_filename, treshnum,MDS=T,summ=F,contignum=1,strain_meta_ind=0,drop_by_abund=F){
  ###The main function. In essence, it is that function that makes dist_m using data from json_distance
  # I.e. from segmental data_matrices it constructs dist_m structure.
  #ALso it filters by coverage and remove NANs.
  print('json_general')
  require('MASS')
  print(dist_filename)
  data.distance <- json_distance(dist_filename, treshnum,summ,contignum)
  distance_m <- data.distance$diff/data.distance$mut_cover
  data.stats <- json_stat(stat_filename, data.distance$names)
  data.cover <- data.frame(cover=data.stats$cover[,contignum],samples = data.stats$cover$samples,stringsAsFactors=F)
  heterovector <- data.stats$heterogen[,contignum]/data.stats$cover[,contignum]
  data.heterogen <- data.frame(hetero=heterovector,samples = data.stats$cover$samples,stringsAsFactors=F)
  ###Stupid way below
  #distance_m[is.nan(distance_m)] <- 0.01
  #distance_m[distance_m==0] <- 0.01
  ###
  L <- drop_by_coverage(data.cover,distance_m)
  distance_m <- L[[1]]
  data.cover <- L[[2]]
  short_strain_name <- as.character(strain_meta.data$short_descript[strain_meta_ind])
  if (drop_by_abund){
    abund_info <- data.frame('short_names'=rownames(abund.data.matrix),'abund'=as.vector(abund.data.matrix[,short_strain_name]))
    abund_info <- merge(x=abund_info,y=data.sample_names2,by.x=c('short_names'),by.y=c('short_names'),all.x=T,all.y=T)
    distance_m <- drop_by_abundance(abund_info, distance_m)
  }
  distance_m <- drop_NAN(distance_m)
  
  if (MDS){
    distance_m <- as.dist(distance_m)
    distance <- isoMDS(distance_m)
    data.points <- as.data.frame(distance$points)
    data.points$names <- row.names(data.points)
    
    data.stats <- json_stat(stat_filename, data.distance$names)
    data.points$coverage <- data.stats$cover[,treshnum]
    data.points$heterogen <- data.stats$heterogen[,1]/data.stats$cover[,1]
    return(data.points)
  }
  else{
    return(list(distance_m,data.cover,data.heterogen))
  }
}

stat_obtain <- function(dist_filename,stat_filename,treshnum=1,summ=F,contignum=1){
  data.distance <- json_distance(dist_filename, treshnum,summ,contignum)
  data.stats <- json_stat(stat_filename, data.distance$names)
  return(data.stats)
}

json_general_hetero <- function(dist_filename, stat_filename, treshnum,MDS=T,summ=F,contignum=1,strain_meta_ind=0,drop_by_abund=F){
  print('json_general_hetero')
  require('MASS')
  data.distance <- json_distance(dist_filename,1,summ,contignum)
  distance_m <- data.distance$heterogen/data.distance$mut_cover
  data.stats <- json_stat(stat_filename, data.distance$names)
  data.cover <- data.frame(cover=data.stats$cover[,contignum],samples = data.stats$cover$samples,stringsAsFactors=F)
  heterovector <- data.stats$heterogen[,contignum]/data.stats$cover[,contignum]
  data.heterogen <- data.frame(hetero=heterovector,samples = data.stats$cover$samples,stringsAsFactors=F)
  
  ###Stupid way below
  #distance_m[is.nan(distance_m)] <- 0.01
  #distance_m[distance_m==0] <- 0.01
  ###
  L <- drop_by_coverage(data.cover,distance_m)
  distance_m <- L[[1]]
  data.cover <- L[[2]]
  short_strain_name <- as.character(strain_meta.data$short_descript[strain_meta_ind])
  if (drop_by_abund){
    abund_info <- data.frame('short_names'=rownames(abund.data.matrix),'abund'=as.vector(abund.data.matrix[,short_strain_name]))
    abund_info <- merge(x=abund_info,y=data.sample_names2,by.x=c('short_names'),by.y=c('short_names'),all.x=T,all.y=T)
    distance_m <- drop_by_abundance(abund_info, distance_m)
  }
  distance_m <- drop_NAN(distance_m)
  
  if (MDS){
    distance_m <- as.dist(distance_m)
    distance <- isoMDS(distance_m)
    data.points <- as.data.frame(distance$points)
    data.points$names <- row.names(data.points)
    
    data.stats <- json_stat(stat_filename, data.distance$names)
    data.points$coverage <- data.stats$cover[,treshnum]
    data.points$heterogen <- data.stats$heterogen[,1]/data.stats$cover[,1]
    return(data.points)
  }
  else{
    return(list(distance_m,data.cover,data.heterogen))
  }
}

###Rewrite####
distances_obtain <- function(path,json_dir_list,treshnum){
  STAT = '/stats.json'
  DIST = '/distance.json'
  distances_list <- list()
  i <- 1
  for (dir in json_dir_list){
    stat_filename <- paste(path,dir,STAT,sep='')
    dist_filename <- paste(path,dir,DIST,sep='')
    distances <- list()
    for (tresh in 1:treshnum){
      print(dist_filename)
      print(tresh)
      distance.data <- json_distance(dist_filename, tresh)
      distance_m <- distance.data$diff/distance.data$mut_cover
      distance_m <- drop_NAN(distance_m)
      distances[[tresh]] <- distance_m
    }
    distances_list[[i]] <- distances
    i <- i + 1
  }
  return(distances_list)
}

drop_NAN <- function(dist_m){
  which_NA <- function(vector){
    return(length(which(is.na(vector))))
  }
  dist_m <- as.matrix(as.dist(dist_m))
  while(length(which(is.na(dist_m)))!=0){
    na_count_vector <- apply(dist_m,1,which_NA)
    max_nums <- which(na_count_vector==max(na_count_vector))
    del_num <- sample(length(max_nums),1)
    del_ind <- max_nums[del_num]
    dist_m <- dist_m[-del_ind,-del_ind]
  }
  return(dist_m)
}


drop_outliers <- function(dist_m){
  print(dim(dist_m))
  mean_val <- mean(dist_m)
  print(mean_val)
  dev <- sd(as.vector(dist_m))
  tresh <- mean_val + 2 * dev
  which_more <- function(vector){
    return(length(which(vector>tresh)))
  }
  out_count_vector <- apply(dist_m,1,which_more)
  len <- length(out_count_vector)
  frac = 0.7
  len_tresh <- len * frac
  #print(len_tresh)
  del_ind <- which(out_count_vector > len_tresh)
  #print(del_ind)
  if (length(del_ind)>0){
    dist_m <- dist_m[-del_ind,-del_ind]
  }
  
  return(dist_m)
}

drop_null <- function(dist_m){
  which_null <- function(vector){
    return(length(which(vector == 0)))
  }  
  while(length(which(dist_m==0))>dim(dist_m)[1]){
    ###   Why len? Not which_null?'''
    count <- apply(dist_m,2,which_null)
    max_ <- max(count[count>1])
    max_indices <- which(count==max_)
    ind <- max_indices[sample(length(max_indices),1,1)]
    dist_m <- dist_m[-ind,-ind]
  }
  return(dist_m)
}

distances_comm <- function(distances_list){
  #distances_list <- unlist(distances_list,recursive=F)
  sample_num <- length(distances_list)
  corr_m <- matrix(0,sample_num,sample_num)
  for (i in 1:sample_num){
    for (j in 1:i){
      print(c(i,j))
      corr_m[i,j] <- dist_comm(distances_list[[i]],distances_list[[j]])
      corr_m[j,i] <- corr_m[i,j]
    }
  }
  return(corr_m)
}

mutality_comm <- function(distances)

dist_comm <- function(d_m1,d_m2){
  require("vegan")
  comm_names <- intersect(row.names(d_m1),row.names(d_m2))
  ind1 <- which(row.names(d_m1) %in% comm_names)
  ind2 <- which(row.names(d_m2) %in% comm_names)
  d_m1 <- d_m1[ind1,ind1]
  d_m2 <- d_m2[ind2,ind2]
  m_size <-dim(d_m1)[1]
  if (is.null(m_size)){
    m_size <- 0
  }
  if (m_size > 10){
    M <- mantel(d_m1,d_m2,method = "spearman")$statistic
  }
  else{
    print(m_size)
    M <- 0.0
  }
  return(M)
}

el_comm <- function(d_m1,d_m2){
  comm_names <- intersect(row.names(d_m1),row.names(d_m2))
  ind1 <- which(row.names(d_m1) %in% comm_names)
  ind2 <- which(row.names(d_m2) %in% comm_names)
  d_m1 <- d_m1[ind1,ind1]
  d_m2 <- d_m2[ind2,ind2]
  m_size <-dim(d_m1)[1]
  if (is.null(m_size)){
    m_size <- 0
  }
  print(m_size)
  return(m_size)
}

els_comm <- function(distances_list){
  #distances_list <- unlist(distances_list,recursive=F)
  sample_num <- length(distances_list)
  count_m <- matrix(0,sample_num,sample_num)
  for (i in 1:sample_num){
    for (j in 1:i){
      print(c(i,j))
      count_m[i,j] <- el_comm(distances_list[[i]],distances_list[[j]])
      print(count_m[i,j])
      count_m[j,i] <- count_m[i,j]
    }
  }
  return(count_m)
}


sample_naming<-function(strainlist,tresholds){
  namelist <- list()
  naming <- function(strain,treshold){
    return(paste(strain,as.character(treshold),sep='_'))
  }
  for (strain in strainlist){
    for (treshold in tresholds){
      namelist <- c(namelist,naming(strain,treshold))
    }
  }
  return(unlist(namelist))
}

json_map <- function(dist_filename,stat_filename,tresh,MDS=T,summ=F,contignum=1,hetero=F,strain_meta_ind=0,drop_by_abund=F){
  ###Pretty useless function. It chooses accroding flags make hetero dist_m or simple variant.
  if (hetero){
   M <- json_general_hetero(dist_filename,stat_filename,MDS,summ,contignum) 
  }
  else{
    M <- json_general(dist_filename,stat_filename,tresh,MDS,summ,contignum,strain_meta_ind,drop_by_abund)
  }
  if (MDS){
    M$class<- str_extract(M$names, "[PB]")
    M$nums <- str_extract(M$names, "[0-9]+")
    return(ggplot(M,aes(x=V1,y=V2,color=class,label=nums)) + geom_point() + geom_text(color=grey,alpha=0.4))
  }
  else{
    return(M)
  }
}

drop_group <- function(dist_m,group='RUS'){
  print(dim(dist_m))
  names <- data.sample_names2$names[data.sample_names2$country == group]
  dist_names <- rownames(dist_m)
  ind <- which(dist_names %in% names)
  if (length(ind)!=0){
    dist_m <- dist_m[-ind,-ind]
  }
  
  return(dist_m)
}

get_points <- function(dist_m, group,group_drop,clustnum=3){
  library(fpc)
  ###This function that firstly filter data by groups, nonzero distances
  # and remover outliers. Then using dist_m make MDS and clasterization
  #!!! It returns coordinates of points
  dist_m <- drop_null(dist_m)
  if (group_drop){
    dist_m <- drop_group(dist_m,group='RUS')
  }
  dist_m <- drop_outliers(dist_m)
  print(dim(dist_m))
  MDS_data <- isoMDS(dist_m)
  points <- as.data.frame(MDS_data$points)
  points$names <- row.names(points)
  points <- merge(points,data.sample_names2,by.x=c('names'),by.y=c('names'),all.x=T,all.y=F)
  dimm <- dim(dist_m)[1]
  cat('Number of samples in dist_m',dimm)
  clusters <- pamk(dist_m,krange=1:4,diss=TRUE)
  avg.widths <- clusters$pamobject$silinfo$clus.avg.widths
  nc <- clusters$nc
  print(nc)
  clusters <- as.data.frame(clusters$pamobject$clustering)
  names(clusters) <- c('clusters')
  clusters$names <- row.names(clusters)
  if (nc > 1){
    if (all(avg.widths < 0.2)){
      clusters$clusters <- 1
    }
  }
  clusters$clusters <- as.factor(clusters$clusters)
  points <- merge(points,clusters,by.x=c('names'),by.y=c('names'),all.x=T,all.y=T)
  return(points) 
}

get_clusters <-function(rootdir,folder,maxcontnum){
  l <- list()
  for (i in 1:maxcontnum){
    dist_path <- paste(rootdir,folder,'distance.json',sep='/')
    stat_path <- paste(rootdir,folder,'stat.json',sep='/')
    M <- json_map(dist_path,stat_path,tresh=1,MDS=F,summ=F,contignum=i,hetero=F)
    M <- M[[1]]
    dist_m <- drop_null(M)
    dist_m <- drop_group(dist_m,group='RUS')
  }
  return(kmeans(dist_m,3))
}

point_plot <- function(contignum=1,group='RUS',group_drop=F,hetero=T,setnum,clustnum=3,only_table=F,drop_by_abund=F,flag=''){
  DIR <- '/data1/bio/runs-kovarskiy/metagenome_data/'
  ind <- which(strain_meta.data$set == setnum)
  rootdir <- paste(DIR,strain_meta.data$set[ind],sep='')
  djson <- paste('distance',flag,'.json',sep='')
  dist_filename <- paste(rootdir,djson,sep='/')
  stat_filename <- paste(rootdir,'stat.json',sep='/')
  library(ggplot2)
  library(RColorBrewer)
  ###Using data from json_map[] and get_points. It merges coordinates with metadata:
  # coverage, heterogenity. Then it produces plot.
  #!!! It returns frame with coordinates and meta and plot.
  pal <-  brewer.pal(4,"Set1")
  names(pal) <-levels(data.sample_names2$country)
  M <- json_map(dist_filename,stat_filename,tresh=1,MDS=F,summ=F,contignum,hetero,strain_meta_ind=ind,drop_by_abund)
  data.cover <- M[[2]]
  data.heterogen <- M[[3]]
  names(data.cover) <- c('cover','names')
  names(data.heterogen) <- c('hetero','names')
  M <- M[[1]]
  #M2 <- M + 0.001
  P <- get_points(M,group,group_drop,clustnum)
  P <- merge(P,data.cover,by.x=c('names'),by.y=c('names'),all.x=T,all.y=F)
  P <- merge(P,data.heterogen,by.x=c('names'),by.y=c('names'),all.x=T,all.y=F)
  if (only_table){
    return(P)
  }
  #P2 <- get_points(M2)
  #inds <- which(P$names=='SRS024388')
  #P <- P[-inds,]
  else{
    name <- plot_name(contignum,setnum,dim(P)[1])  
    plot_1 <- ggplot(P,aes(x=V1,y=V2,color=country,size=cover))+geom_point()+scale_color_manual(name = "country",values = pal)+scale_size(range=c(1,3))+ggtitle(label=name)#+geom_text(aes(label=names))
    #plot_1 <- ggplot(P,aes(x=V1,y=V2,color=clusters,size=cover))+geom_point()+scale_size(range=c(1,3))+ggtitle(label=name)
    return(list(plot_1, P))
  }
  
  #return(length(which(!is.na(P[,2]))))
}

plot_name <- function(contignum,setnum,samplnum){
  ind <- which(strain_meta.data$set == setnum) 
  #cont <- strain_meta.data$cont[ind]
  short_descript <- strain_meta.data$short_descript[ind]
  length <- paste('Length = ',as.character(strain_meta.data$length[ind]),sep='')
  samplnum <- paste('Num of samples = ',as.character(samplnum),sep='')
  name <- paste(short_descript,length,samplnum,sep='\n')
  return(name)
}


multiplot <- function(rootdir,folder,maxcontnum,setnum,subsetnum){
  l <- list()
  K <- c()
  for (i in 1:maxcontnum){
    dist_path <- paste(rootdir,folder,'distance.json',sep='/')
    stat_path <- paste(rootdir,folder,'stat.json',sep='/')
    my_plot <- point_plot(dist_path,stat_path,i,hetero=F,group_drop=T,setnum=setnum,subsetnum=subsetnum)
    l <- c(l,my_plot)
    K <- c(K,dim(my_plot[[2]])[1])
    str_num <- as.character(i)
    img_path <- paste(rootdir,'/',folder,'/plot_',folder,str_num,'.png',sep='')
    ggsave(img_path,my_plot[[1]],width=9,height=7)
  }
  return(list(l,K))
}

good_proc <- function(){
  folders <- as.character(c(1:15))
  KK <- c()
  for (folder in folders){
    subsetnum <- as.numeric(folder)
    K <- multiplot('/data1/bio/runs-kovarskiy/additional_pileup_files5',folder,1,5,subsetnum)
    KK <- c(KK,K[[2]])
  }
  return(KK)
}

macro_f <- function(max_num){
  for (i in 0:max_num){
    num <- paste('0',as.character(i),sep='')
    multiplot('/data1/bio/runs-kovarskiy/additional_pileup_files',num,6,1,i)
  }
}


get_distance_m_list <- function(setnums,drop_by_abund=F){
  DIR <- '/data1/bio/runs-kovarskiy/metagenome_data/'
  distances_list <- list()
  len <- length(setnums)
  for (i in 1:len){
    setnum <- setnums[i]
    ind <- which(strain_meta.data$set == setnum)
    rootdir <- paste(DIR,strain_meta.data$set[ind],sep='')
    dist_filename <- paste(rootdir,'distance.json',sep='/')
    stat_filename <- paste(rootdir,'stat.json',sep='/')
    M <- json_map(dist_filename,stat_filename,tresh=1,MDS=F,summ=F,1,hetero=F,drop_by_abund=drop_by_abund)
    dist_m <- M[[1]]
    dist_m <- drop_group(dist_m, group = 'RUS')
    dist_m <- drop_outliers(dist_m)
#       clusters <- kmeans(dist_m,4)
#       for (clust_num in 1:4){
#         file <- paste(dir_name,'cont_',j,'_clust_',clust_num,'_varpos_variants.txt',sep='')
#         print(file)
#         write.table(names(clusters$cluster[clusters$cluster==clust_num]),file=file,quote=F,row.names=F,col.names=F)
#       }
    distances_list[[length(distances_list)+1]] <- dist_m
      #dist_m <- drop_null(dist_m)  
  }
  return(distances_list)
}

get_all_distance_m <- function(metaind){
  DIR <- '/data1/bio/runs-kovarskiy/metagenome_data'
  j <- 0
  for (i in metaind)
  {
    print(i)
    rootdir <- paste(DIR,as.character(strain_meta.data$set[i]),sep='/')
    dist_path <- paste(rootdir,'distance.json',sep='/')
    stat_path <- paste(rootdir,'stat.json',sep='/')
    M <- json_map(dist_path,stat_path,tresh=1,MDS=F,summ=F,contignum=1,hetero=F)
    M <- M[[1]]
    dist_m <- drop_null(M)
    dist_m <- drop_group(dist_m,group='RUS')
    dist_m[lower.tri(dist_m, diag=TRUE)] <- NA
#     dx <- as.data.frame(as.table(mtx))
#     dx <- subset(dx, !is.na(Freq))
    data.dist <- melt(dist_m)
    names(data.dist) <- c('names1','names2','val')
    bad_ind <- which(is.na(data.dist$val))
    data.dist <- data.dist[-bad_ind,]
    data.dist$descript <- paste(strain_meta.data$short_descript[i],as.character(strain_meta.data$set[i]),sep='')
    data.dist$sdescript <- strain_meta.data$short_descript[i]
    data.dist$length <- strain_meta.data$length[i]
    #data.dist$frac <- strain_meta.data$frac[i]
    if (j == 0){
      res <- data.dist
    }
    else{
      res <- rbind(res,data.dist)
    }
    j <- j + 1
  }
  return(res)
}

cols.gentleman <- function(ncol=500) {
  library(RColorBrewer)
  hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(ncol)
  return(rev(hmcol))
}

get_clusters_new <-function(setnum,clustnum=1,flag=''){
  library(fpc)
  DIR <- '/data1/bio/runs-kovarskiy/metagenome_data/'
  print (setnum)
  rootdir <- paste(DIR,as.character(setnum),sep='/')
  jsonname <- paste('distance',flag,'.json',sep='')
  dist_path <- paste(rootdir,jsonname,sep='/')
  stat_path <- paste(rootdir,'stat.json',sep='/')
  M <- json_map(dist_path,stat_path,tresh=1,MDS=F,summ=F,contignum=1,hetero=F,drop_by_abund=T,strain_meta_ind=which(strain_meta.data$set==setnum))
  M <- M[[1]]
  dist_m <- drop_null(M)
  dist_m <- drop_group(dist_m,group='RUS')
  dimm <- dim(dist_m)[1]
  cat('Number of samples in dist_m',dimm)
  clusters <- pamk(dist_m,krange=1:4,diss=TRUE)
  res <- as.data.frame(clusters$pamobject$clustering)
  #names(clusters) <- c('clusters')
  res$names <- row.names(res)
  return(clusters)
}

get_cluster_data <- function(samplenums,fixed_clustnum=TRUE,clustnum=3)
{   
  data <- data.sample_names2[,c('names','country')]
  for (i in samplenums)
  {
    print(i)
    if (fixed_clustnum==FALSE){
      clustnum <- strain_meta.data$clustnum[i]
    }
    setnum <- i
    clust_data <- get_clusters_new(setnum, clustnum)
    nc <- clust_data$nc
    avg.widths <- clust_data$pamobject$silinfo$clus.avg.widths
    print(clust_data$nc)
    print(clust_data$pamobject$silinfo$clus.avg.widths)
    clust_data <- clust_data$pamobject$clustering
    
    if (nc > 1){
      if (all(avg.widths < 0.1)){
        #clust_data <- data.frame(temp=1,names=names(clust_data))
      }
      else{
        clust_data <- data.frame(temp=clust_data,names=names(clust_data))
        names(clust_data) <- c(as.character(strain_meta.data$short_descript[which(strain_meta.data$set==i)]),'names')
        data <- merge(data,clust_data,by.x=c('names'),by.y=c('names'),all.x=T,all.y=F)
      }
    }
    else{
      #clust_data <- data.frame(temp=clust_data,names=names(clust_data))
    }
  }
  #nms <- str_replace(data$names,'\\_\\d+','')
  #nms <- str_replace(nms,'merged','')
  #data$names <- nms
  return(data)
}

# Spearman distance
distSpear<-function (x, ...) 
{
  result <- 1-abs(cor(t(x), method='spearman',use="pairwise.complete.obs" )) #, use='pairwise')
  result <- drop_NAN(result)
  return(as.dist(result))
}


get_cluster_stat <- function(data){
  library(ggplot2)
  len_NA <- function(v)
  {
    return(length(which(is.na(v))))
  }
  maxlen <- dim(data)[2]
  data <- data[-which(apply(data,1,len_NA)>0.9 * maxlen),]
  data.matrix <- as.matrix(data[,-c(1,2)])
  row.names(data.matrix) <- data$names
  #hr <- hclust(distSpear(data), method="ward")
  dist_m <- distSpear(t(data.matrix))
  hc <- hclust(dist_m, method="ward")
  data <- melt(data)
  names(data) <- c("samples","country","strain","clust")
  data <- data[which(data$strain %in% hc$label),]
  data$strain <- factor(data$strain,levels=hc$label[hc$order])
  p <- ggplot(data, aes(samples, strain)) + geom_tile(aes(fill = factor(clust),colour = "white")) + scale_fill_brewer(palette='Set1')+facet_grid(.~country,scales="free")
  return(p)
}

get_strain_clusters <- function(data){
  #Obtain correlation between structure of clustering
  library(ggplot2)
  len_NA <- function(v)
  {
    return(length(which(is.na(v))))
  }
  maxlen <- dim(data)[2]
  inds <- which(apply(data,1,len_NA)>0.9 * maxlen)
  if (length(inds) > 0){
    data <- data[-inds,]
  }
  data.matrix <- as.matrix(data[,-c(1,2)])
  row.names(data.matrix) <- data$names
  #hr <- hclust(distSpear(data), method="ward")
  dist_m <- distSpear(t(data.matrix))
  hc <- hclust(dist_m, method="ward")
  dist_m <- as.matrix(dist_m)
  dist_m <- melt(dist_m)
  names(dist_m) <- c('strain1','strain2','distance')
  dist_m$strain1 <- factor(dist_m$strain1,levels=hc$label[hc$order])
  dist_m$strain2 <- factor(dist_m$strain2,levels=hc$label[hc$order])
  p <- ggplot(dist_m, aes(strain1, strain2)) + geom_tile(aes(fill = distance,colour = "white")) + scale_fill_gradient(low='red',high='white')
  return(p)
  
}


multiplot_all <- function(samplenums)
{  
  # data <- data.sample_names2[data.sample_names2$country != 'RUS', c('names','country')]
  # num_of_runs <- dim(strain_meta.data)[1]
  plotlist <- list()
  j <- 1
  for (i in samplenums)
  {
    print(i)
    my_plot <- point_plot(1,hetero=F,group_drop=T,setnum=strain_meta.data$set[i],clustnum=3)
    plotlist[[j]] <- my_plot[[1]] + theme(legend.position="none",plot.title = element_text(size = rel(0.9)),axis.title.x = element_blank(),axis.title.y = element_blank()) 
    j <- j + 1
    #l <- c(l,my_plot)
    #K <- c(K,dim(my_plot[[2]])[1])
    #str_num <- as.character(strain_meta.data$contnum[i])
    img_path <- paste('/data1/bio/runs-kovarskiy/metagenome_data/imgs/','plot_all_',as.character(strain_meta.data$set[i]),'.png',sep='')
    #img_path <- paste(rootdir,'/',strain_meta.data$str_subset[i],'/plot_',strain_meta.data$str_subset[i],str_num,'.png',sep='')
    ggsave(img_path,my_plot[[1]],width=9,height=7)
  }
  return(plotlist)
}

get_all_mds_frames <- function(setnums){
  len <- length(setnums)
  for (i in 1:len){
    setnum <- setnums[i]
    temp <- point_plot(1,hetero=F,group_drop=T,setnum=setnum,clustnum=3,only_table=T,drop_by_abund=F)
    temp$strain <- strain_meta.data$short_descript[strain_meta.data$set==setnum]
    if (i == 1){
      res <- temp
    }
    else{
      res <- rbind(res, temp)
    }
  }
  return(res)
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

preprocess_abund <- function(data)
{
  #data <- drop_NA(data)
  datm <- as.matrix(data[,-1])
  datm <- 100 * datm / rowSums(datm)
  rownames(datm) <- data[,1]
  # norm to 100%
  #datm <- datm[,which(apply(datm,2,max) >= min_tax_perc)] # drop low count taxons
  datm <- datm[data.sample_names2$short_names,]
  return(datm)
}

process_abund <- function(datm, cluster_data){
  samples <- cluster_data$names
  #cluster_data <- as.data.frame(as.numeric(cluster_data[,-c(1,2)]))
  strains <- names(cluster_data)[-c(1,2)]
  cluster_data <- as.data.frame(cluster_data[,-c(1,2)])
  names(cluster_data) <- strains
  inds <- which(colnames(datm) %in% strains)
  #residue <- datm[strsamples,-inds]
  residue <- datm[samples,]
  #datm <- datm[samples,inds]
  nstrains <- length(strains)
  nsamples <- length(samples)
  nameslist <- list()
  #
  #cluster_data$buff <- 2
  #str(cluster_data)
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
        print(samples[j])
        print(strains[i])
        print(abund.data.matrix[samples[j],strains[i]])
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
    plottemp <- melt(temp)
    names(plottemp) <- c("strain1","strain2","corr_coeff")
    maxes <- apply(temp[-1,],2,function(i) {return(max(abs(i)))})
    #print(temp[,which(maxes>0.4 & maxes != 1.0)])
    p <- ggplot(plottemp, aes(x=strain2,y=strain1)) + geom_tile(aes(fill = corr_coeff,colour = "white")) + scale_fill_gradientn(values=seq(-1, 1, length = 4),colours = c("#FF0000","#000000","#000000","#000000","#00FF00"))
    plot(p)
    #print('--')
    #print(p)
    #print(p, vp=viewport(layout.pos.row=1,layout.pos.col=1:9))
    
  }
  return(corrm)
}

group_comparison <- function(distance_m, meta_info){
  levelnum <- length(levels(meta_info$spec))
  dimnames <- levels(meta_info$spec)
  comp_m <- matrix(nrow=levelnum,ncol=levelnum,dimnames=list(dimnames,dimnames))
  res <- data.frame(group1=c(),group2=c(),mean=c(),dist=c())
  for (i in 1:levelnum){
    for (j in 1:levelnum){
      names1 <- meta_info$names[which(meta_info$spec == dimnames[i])]
      names2 <- meta_info$names[which(meta_info$spec == dimnames[j])]
      sub_m <- distance_m[names1, names2]
      v <- as.vector(sub_m)
      temp <- data.frame(group1=dimnames[i], group2=dimnames[j], mean=mean(sub_m), dist=v)
      comp_m[dimnames[i],dimnames[j]] <- mean(sub_m)
      res <- rbind(res, temp)
    }
  }
  return(list(comp_m, res))
}
#nms <- str_replace(dat$names,'\\_\\d+','')
#nms <- str_replace(nms,'merged','')