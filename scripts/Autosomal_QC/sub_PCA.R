###################################
### PCA plotting
### date: 28-03-2019
### version: 0.01
### authors: EL - RAG
###################################
### New
###################################

### packages

library(tidyverse)
library(optparse)
library(data.table)
library(gridExtra)
library(readxl)
library(viridis)

###test data
## opt<-list()
## opt$out<-"/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC_14mil/plots"
## opt$input<-"/Users/e.a.lopera.maya/Tasks/testdata"
## opt$input<-"/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC_14mil/6_PCA"
## '/apps/data/1000G/populationInfo/20130606_sample_info_pop_superpop.xlsx'
## opt$file<-"second_PCA.eigenvec"
## opt$file<-"PCA_1000G.eigenvec"
## opt$res<-"2nd."
## opt$res<-"1st."
## G_pops<-read_excel('~/Tasks/testdata/20130606_sample_info_pop_superpop.xlsx')

#opt$file<-"PCA_project1.eigenvec"
#########################################################################################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input path to perform QC", metavar="character"),
  
  make_option(c("-f", "--file"), type="character", default="PCA_1000G.eigenvec", 
              help="Output path to save report", metavar="character"),
  
  make_option(c("-r", "--res"), type="character", default="1st.", 
              help="Output path to save report", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="plots", 
              help="Output path to save report", metavar="character")
);

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied as input path", call.=FALSE)
}
######################################Main#######################################
###processing for 1000G

#loading data
PCA.file <- file.path(opt$input, opt$file)
G_eigenvec_table<-read.table(PCA.file, header = F)

####take reference file from the cluster
G_pops<-read_excel('/apps/data/1000G/populationInfo/20130606_sample_info_pop_superpop.xlsx')

##preparing data
colnames(G_eigenvec_table)<-c('FID','Sample ID',paste("PC",1:20))
G_eigenvec_table$`Sample ID`<-as.character(G_eigenvec_table$`Sample ID`)
G_eigenvec_table$FID<-as.character(G_eigenvec_table$FID)
G_eigenvec<-left_join(G_eigenvec_table,G_pops, by='Sample ID')
G_eigenvec[G_eigenvec$FID  %like% "gonl-",c("Population", "superpop")]<-"GoNL"
G_eigenvec[is.na(G_eigenvec$Population),c("Population", "superpop")] <- "Cohort"    
levels(G_eigenvec$Population)<-c(levels(G_eigenvec$Population),"Cohort")


#### preparing data for second PCA Plot ####
G_eigenvec$eurplot<-ifelse((G_eigenvec$superpop=="EUR") | (G_eigenvec$Population== "Cohort")| 
                             (G_eigenvec$Population== "GoNL") ,"eurplot","notplot")
pc2filt<-G_eigenvec[G_eigenvec$eurplot=="eurplot" & G_eigenvec$Population != "Cohort",]

##PCA European population definitions
multiplier<-ifelse(opt$res=="2nd.",2,3)

pc1min<-min(pc2filt[,3])-multiplier*abs(sd(pc2filt[,3]))
pc1max<-max(pc2filt[,3])+multiplier*abs(sd(pc2filt[,3]))
pc2min<-min(pc2filt[,4])-multiplier*abs(sd(pc2filt[,4]))
pc2max<-max(pc2filt[,4])+multiplier*abs(sd(pc2filt[,4]))

##PCA plot scale limits
xlimmin<-min(G_eigenvec[,3])
xlimmax<-max(G_eigenvec[,3])
ylimmin<-min(G_eigenvec[,4])
ylimmax<-max(G_eigenvec[,4])

##make samples files
G_eigenvec$sampletype<-ifelse((G_eigenvec$`PC 1`< pc1max) &
                                (G_eigenvec$`PC 1`> pc1min) &
                                (G_eigenvec$`PC 2`< pc2max) &
                                (G_eigenvec$`PC 2`> pc2min)&
                                (G_eigenvec$eurplot=="eurplot"),"eur","NE" )

eur.samples<-G_eigenvec[G_eigenvec$sampletype=="eur",c(1,2)]
sample.list<-G_eigenvec[G_eigenvec$sampletype=="eur" &
                          (G_eigenvec$Population=="Cohort"),c(1,2)]
excl.list<-G_eigenvec[G_eigenvec$sampletype=="NE" &
                        (G_eigenvec$Population=="Cohort"),c(1,2)]


eur.sample.file <- file.path(opt$input,paste0(opt$res, "eur.samples"))
write.table(eur.samples,eur.sample.file, row.names = F, quote=F)

sample.list.file <- file.path(opt$input,paste0(opt$res, "incl.samples"))
write.table(sample.list,sample.list.file, row.names = F, quote=F)

excl.sample.file <- file.path(opt$input,paste0(opt$res, "excl.samples"))
write.table(excl.list,excl.sample.file, row.names = F, quote=F)

##aesthetics

included<-nrow(sample.list)

ELcolor<-c("chartreuse3",  "plum3", "brown", "orange" , "blue","cyan1", 
           "darkgray","forestgreen","blueviolet","pink", colorspace::rainbow_hcl(12),colorspace::diverge_hcl(4),viridis(7))
n<-names(ELcolor)<-c("EUR","AFR", "Cohort","GoNL","SAS", "EAS","AMR")
pops<-unique(G_eigenvec$Population)
pat<-paste(n,collapse="|")
names(ELcolor)<-c(n,pops[grepl(pat,pops)==FALSE])


##plotting

if (opt$res=="2nd.") {pops.plot<- ggplot(G_eigenvec,aes(x=G_eigenvec$`PC 1`, y=G_eigenvec$`PC 2`))+
  geom_point(aes(colour=Population, shape=Population),alpha=0.4,size=2)+
  theme_bw()+
  scale_colour_manual(values=ELcolor)+
  guides(colour = guide_legend(override.aes = list(alpha=1)))+
  scale_shape_manual(values=c(16,17,rep(16,5)))+
  xlab("PC 1")+
  ylab("PC 2")+
  ggtitle("1000G Populations")+
  theme(text=element_text(size=10, family = 'Helvetica'),plot.title = element_text(hjust = 0.5))
} else {
pops.plot<- ggplot(G_eigenvec,aes(x=G_eigenvec$`PC 1`, y=G_eigenvec$`PC 2`))+
  geom_point(aes(colour=Population, shape=Population),alpha=0.6,size=2)+
  theme_bw()+
  scale_colour_manual(values=ELcolor) +
  guides(colour = guide_legend(override.aes = list(alpha=1)))+
  scale_shape_manual(values=c(rep(16,8),17,rep(16,19)))+
  xlab("PC 1")+
  ylab("PC 2")+
  ggtitle("1000G Populations")+
  theme(text=element_text(size=10, family = 'Helvetica'),plot.title = element_text(hjust = 0.5))
}

if (opt$res=="2nd.") superpops.plot<-"No Superpopulations" else {
  
  superpops.plot<-ggplot(G_eigenvec,aes(x=G_eigenvec$`PC 1`, y=G_eigenvec$`PC 2`))+
    geom_point(aes(colour=superpop, shape=superpop),alpha=0.6, size=2)+
    theme_bw()+
    scale_colour_manual(values=ELcolor)+ 
    guides(colour = guide_legend(override.aes = list(alpha=1)))+
    scale_shape_manual(values=c(rep(16,2),17,rep(16,4)))+ 
    xlab("PC 1")+ylab("PC 2")+ 
    ggtitle("1000G Superpopulations")+ 
    theme(text=element_text(size=10, family = 'Helvetica'), plot.title = element_text(hjust = 0.5))
}


eurpops.plot<-ggplot(G_eigenvec[G_eigenvec$eurplot=="eurplot",],aes(x=`PC 1`, y=`PC 2`))+
  coord_cartesian(xlim = c(xlimmin,xlimmax), ylim = c(ylimmin,ylimmax))+
  geom_rect(aes(xmin = pc1min, ymin = pc2min,
                xmax= pc1max,ymax = pc2max ),fill="khaki1", alpha=0.01, linetype=2)+
  geom_point(aes(colour=Population, shape=Population),alpha=0.6,size=2)+
  guides(colour = guide_legend(override.aes = list(alpha=1)))+
  theme_bw()+
  scale_colour_manual(values=ELcolor)+
  scale_shape_manual(values=c(16,17,rep(16,5)))+
  xlab("PC 1")+ylab("PC 2")+
  ggtitle("1000G European populations")+
  theme(text=element_text(size=10, family = 'Helvetica'),plot.title = element_text(hjust = 0.5))

### saving out plots. 
combine.plot.title <- paste0(opt$res," PCA analysis with 1000G", "\n", opt$name, " ", date())
PCA.plot.file <- file.path(opt$out,paste0( "08_PCA.",opt$res,"plot.tiff"))
precancestry<-ifelse(opt$res=="2nd.","Non-Finnish European","European samples")
plot.conclude<-paste0(included,"\n",precancestry)
#plotting
tiff(PCA.plot.file,  width = 2000, height = 4000, 
     units = "px", res = 300, compression = "lzw")

if (opt$res=="2nd.") {
  grid.arrange(pops.plot, 
               eurpops.plot,
               top=combine.plot.title,
               bottom=plot.conclude,
               ncol=1)
  } else {
  grid.arrange(superpops.plot,
               pops.plot, 
               eurpops.plot,
               top=combine.plot.title,
               bottom=plot.conclude,
               ncol=1)
}
  
  
dev.off()






