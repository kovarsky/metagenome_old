source('rus_test_funcs.R')
PATH <- '/data1/bio/runs-kovarskiy/metagenome_data/rus_test/'
data.rus.strains <- read.csv('/data1/bio/runs-kovarskiy/metagenome_data/rus_test/strain_info',sep='\t',header=F)
names(data.rus.strains) <- c("descript","length","set","short_descript")
data.rus.cov_profiles <- load_cov_distrs()
data.rus.meta <- read.csv('/data1/bio/runs-kovarskiy/metagenome_data/rus_test/rus_test.meta.tsv',sep='\t')
data.rus.coverages <- cov_distrs_proc()


data1 <- load_distances('alt_2',0.05)
temp1 <- make_data_frame_mds(data1[[3]])
ggplot(temp1,aes(x=V1,y=V2,color=sample,label=short_names,shape=sequen))+geom_point(size=2)+geom_text()



data2 <- load_distances('alt_3',0.1)
temp2 <- make_data_frame_mds(data2[[3]])
ggplot(temp2,aes(x=V1,y=V2,color=sample,label=short_names,shape=sequen))+geom_point(size=2)+geom_text()
