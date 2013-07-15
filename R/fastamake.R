setnums <- c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,47,48,49,56,57,58,59,62,68,73,74,75,82,85,87,89,96,98,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122)
for (setnum in setnums){
  print('------')
  print(setnum)
  data <- load_distances(setnum, suffix = 'alt_4', treshold = 0.6)
  if (length(dim(data[[1]])[1]) > 0){
    if (dim(data[[1]])[1] > 20){
      m_data <- drop_less(data[[1]], treshold=10^6)
      path <- paste('/data1/bio/runs-kovarskiy/metagenome_data',setnum,'fastasamples.list',sep='/')
      #print(path)
      #write( rownames(m_data), path)
      print(dim(m_data))
    }
  }
}


treshold <- 0.6
max_num <- 10
covers_and_countries <- merge(data.cov_distrs,data.sample_names2[,c(1,2,4)],by.x=c("names"),by.y=c("names"),all.x=T,all.y=F)
covers_and_countries <- covers_and_countries[covers_and_countries$frac > treshold,]

for (setnum in setnums){
  print(setnum)
  temp <- covers_and_countries[covers_and_countries$setnums == setnum,]
  temp <- temp[order(temp$frac),]
  eur <- as.character(temp$names[temp$country=='EUR'])[1:max_num]
  chn <- as.character(temp$names[temp$country=='CHN'])[1:max_num]
  usa <- as.character(temp$names[temp$country=='USA'])[1:max_num]
  res <- c(eur,chn,usa)
  print(length(res[which(!is.na(res))]))
  #data <- load_distances(setnum, suffix = 'alt_4', treshold = 0.6)
  #if (length(dim(data[[1]])[1]) > 0){
    #if (dim(data[[1]])[1] > 20){
     # m_data <- drop_less(data[[1]], treshold=10^6)
  path <- paste('/data1/bio/runs-kovarskiy/metagenome_data',setnum,'fastasamples.list',sep='/')
  print(path)
  #write(res[which(!is.na(res))], path,ncolumns=1)
}
