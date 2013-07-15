load_gene_cover <- function(flag,setnum){
  path <- paste('/data1/bio/runs-kovarskiy/metagenome_data/',setnum,'/',sep='')
  sample_names <- read.csv(paste(path,flag,'_names.list',sep=''),header=F,sep='\t',stringsAsFactors=FALSE)$V1
  means <- data.cov_distrs[data.cov_distrs$setnums==setnum & as.character(data.cov_distrs$names) %in% sample_names,]
  means$names <- as.character(means$names)
  len <- length(means$names)
  data <- read.csv(paste(path,flag,'_gene_covers.tsv',sep=''),header=F,sep='\t')
  names(data) <- c('short_id',sample_names)
  res <- data.frame(short_id=data$short_id)
  for (i in 1:len){
    #print(means$redund[i])
    if (means$redund[i] >= 5){
      res <- cbind(res,data[,c(means$names[i])] / means$mean[i])
    }
      
    #data[,c(means$names[i])] <- 
  }
  return(res)
}


get_gene_info <- function(flag, setnum){
  data <- load_gene_cover(flag,34)
  temp <- data[,-1]
  nonzeros <- apply(temp,1,function(x){length(which(x!=0))})
  means <- apply(temp,1,function(x){sum(x)})/nonzeros
  res <- data.frame(short_id=data[,1],nonzeros,means)
  names(res) <- c("short_id",paste("nonzero",flag,sep="."),paste("means",flag,sep="."))
  return(res)
}

res1 <- get_gene_info("chn",34)
res2 <- get_gene_info("eur",34)
res3 <- get_gene_info("usa",34)
res <- cbind(res1,res2[,-1],res3[,-1])
total_res <- merge(res,anno.data[,c("short_id","RefSeq.Locus.Tag","Locus.Tag","Length","Product")],by=c("short_id"),all=T)

conserved_ids <- read.csv("/data1/bio/runs-kovarskiy/metagenome_data/34/conserved_genes_ids.list",header=F)$V1
total_res$type <- "other"
total_res$type[total_res$short_id %in% conserved_ids] <- "conservative"

ggplot(total_res,aes(x=means.chn,y=means.usa))+geom_point(alpha=0.3)
ggplot(total_res,aes(x=nonzero.chn,y=nonzero.usa))+geom_point(alpha=0.3,size=5)
ggplot(total_res,aes(x=type,y=means.usa))+geom_boxplot()


