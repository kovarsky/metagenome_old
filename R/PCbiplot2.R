PCbiplot2 <- function(PC, x="PC1", y="PC2",title="_") {
  # PC being a prcomp object
  groups <- substr(row.names(PC$x),1,3)
  #groups <- data.frame(name)
  data <- data.frame(obsnames=row.names(PC$x), PC$x,groups=groups)
  data$names <- data$obsnames
  data <- merge(data,data.sample_names2,by.x=c("names"),by.y=c("names"),all.x=T,all.y=F)
  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.6, size=3, aes(label=obsnames))
  plot <- plot + geom_point(aes(color=country),size=3) + scale_color_brewer(palette='Set1')
  plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  #plot <- plot  + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1,hjust=1, alpha =1/2, color="red")
  #plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
  plot <- plot + ggtitle(title) + theme(plot.title = element_text(lineheight=.8, face="bold"))
  return(plot)
}