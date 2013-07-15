setnum <- 40
filename <- paste(PATH,setnum,'fastasamples.list',sep='/')
data <- read.table(file=filename,stringsAsFactors=F)$V1
old_data <- load_distances2(setnum,treshold=0.6)
old_data <- old_data[[3]][data,data]
main_tree <- upgma(old_data)

boot_dists <- data.frame(frac=c(),dist=c(),sd=c())
fracs <- c(0.2,0.5,0.8,0.9)
for (frac in fracs){
  print(frac)
  distances <-  load_distances_bootstrapped(setnum,frac,posnums)
  trees <- lapply(distances, upgma)
  class(trees) <- "multiPhylo"
  treedists <- unlist(lapply(trees, function(x){RF.dist(main_tree,x)}))
  res <- data.frame(frac=frac,dist=mean(treedists),sd=sd(treedists))
  boot_dists <- rbind(boot_dists, res)
  }
ggplot(boot_dists,aes(x=frac,y=dist))+geom_line()+geom_errorbar(aes(ymin=dist-sd, ymax=dist+sd))



distances <-  load_distances_bootstrapped(setnum,0.5,100)
trees <- lapply(distances, upgma)
class(trees) <- "multiPhylo"
main_tree$node.label <- as.character(prop.clades(main_tree,trees))
plot.phylo(main_tree,show.node.label=T)


#plotBS(main_tree, trees, type="unrooted", bs.col="red", bs.adj=NULL)
consistancies <- read.csv(paste(PATH,setnum,'/consistancies.list',sep=''),sep='\t',header=F)$V1
main_tree$node.label <- as.character(consistancies)
plot.phylo(main_tree,show.node.label=T)