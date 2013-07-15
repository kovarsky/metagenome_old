cov_filter <- function(data, treshold){
  counts <- apply(data,1,function(x){length(which(x>treshold))})
  means <- apply(data,1,function(x){mean(x[which(x>treshold)])})
  medians <- apply(data,1,function(x){median(x[which(x>treshold)])})
  sds <- apply(data,1,function(x){sd(x[which(x>treshold)])})
  res <- data.frame(count=counts,mean=means,median=medians,sd=sds)
  return(res)
}