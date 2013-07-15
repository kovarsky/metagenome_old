


get_dual <- function(data){
  data$macro <- 'CHN_vs_EUR'
  data$macro[data$groups %in% c('EUR_vs_USA','USA_vs_EUR')] <- 'EUR_vs_USA'
  data$macro[data$groups %in% c('CHN_vs_USA','USA_vs_CHN')] <- 'CHN_vs_USA'
  return(data)
}

data <- get_dual(data.p_vals.alt_4)

get_any <- function(data, treshold){
  temp <- aggregate(data$p_values,by=list(data$set, data$macro),FUN=function(x){return(any(x<treshold))})
  return(temp)
}

grouped_data <- get_any(data,0.001)
temp.chn_vs_eur <- grouped_data$Group.1[which(grouped_data$Group.2=="CHN_vs_EUR" & grouped_data$x)]
temp.chn_vs_usa <- grouped_data$Group.1[which(grouped_data$Group.2=="CHN_vs_USA" & grouped_data$x)]
temp.eur_vs_usa <- grouped_data$Group.1[which(grouped_data$Group.2=="EUR_vs_USA" & grouped_data$x)]

write(temp.chn_vs_usa,file="chn_vs_usa",sep='\n')
write(temp.eur_vs_usa,file="eur_vs_usa",sep='\n')
write(temp.chn_vs_eur,file="chn_vs_eur",sep='\n')



get_all <- function(data, treshold){
  temp <- aggregate(data$p_values,by=list(data$set, data$macro),FUN=function(x){return(all(x<treshold))})
  return(temp)
}

grouped_data2 <- get_all(data,0.001)
#grouped_data2 <- grouped_data2$x
temp.chn_vs_eur <- grouped_data2$Group.1[which(grouped_data2$Group.2=="CHN_vs_EUR" & grouped_data2$x)]
temp.chn_vs_usa <- grouped_data2$Group.1[which(grouped_data2$Group.2=="CHN_vs_USA" & grouped_data2$x)]
temp.eur_vs_usa <- grouped_data2$Group.1[which(grouped_data2$Group.2=="EUR_vs_USA" & grouped_data2$x)]

write(temp.chn_vs_usa,file="chn_vs_usa",sep='\n')
write(temp.eur_vs_usa,file="eur_vs_usa",sep='\n')
write(temp.chn_vs_eur,file="chn_vs_eur",sep='\n')


