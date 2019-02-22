###################################
### Heterozygosity filtering and report for autosomal QC
### date: 12-12-2018
### version: 0.01
### authors: EL
###################################

## example run 

#RScript Het_autosomeQC.R -i ***/autosome_QC
#  -o repout


#########################################################################################################
#########################################################################################################
### Functions and libraries
#########################################################################################################
#########################################################################################################

###packages
library(ggplot2)
library(optparse)
library(gridExtra)
library(dplyr)

##Arguments
#########################################################################################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input path ", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="Output path to save report", metavar="character")
); 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied as input path", call.=FALSE)
}

output <- as.character(opt$out)
#opt<-list()
#opt$input<-"C:/VirtualM/Analysis2_QC"
#########################################################################################################
### Main
#########################################################################################################
##get files
het.file <- file.path(opt$input,"autosomal.het")
roh.file<-file.path(opt$input,"autosomal.hom.indiv")
cr.file<-file.path(opt$input,"CR.samples")
##Organize data
hetdata<-read.table(het.file, header = T)
hetdata$OHET<-(hetdata$N.NM.- hetdata$O.HOM.)
rohdata<-read.table(roh.file, header=T)
crdata<-read.table(cr.file,header=F)
colnames(crdata)<-c("IID", "FID","F_MISS")
crdata<-crdata%>% group_by(IID,FID)%>%summarise(F_miss=mean(F_MISS))

##unify data in a  single database
hetdata<-right_join(hetdata,rohdata,by=c("IID","FID"))
hetdata<-inner_join(hetdata,crdata,by=c("IID","FID"))
##make filters to exclude and include data
Excluded<-hetdata[abs(hetdata$OHET-mean(hetdata$OHET))>4*sd(hetdata$OHET),c(1,2)]
Included<-hetdata[abs(hetdata$OHET-mean(hetdata$OHET))<=4*sd(hetdata$OHET),c(1,2)]
##make files for included/excluded individuals
excl.file <- file.path(opt$input,"Excluded.het")
incl.file <- file.path(opt$input,"Included.het")
write.table(Excluded,excl.file, quote=F,row.names = F)
write.table(Included,incl.file, quote=F,row.names = F)

###plotting
excluded<-nrow(Excluded)
plot.conclude<-paste0(excluded," ","Samples excluded")
CR.het<-ggplot(hetdata, aes(x=log(F_miss),y=OHET))+
             geom_point(alpha=0.4, size=2)+
             ggtitle("Heterozygosity distribution by sample")+
             theme_bw()+
             xlab("log(missing)")+ylab("Heterozygosity")+
             geom_hline(yintercept=mean(hetdata$OHET)-4*sd(hetdata$OHET),color="red")+
             geom_hline(yintercept=mean(hetdata$OHET)+4*sd(hetdata$OHET),color="red")+
             theme(text=element_text(size=10, family="Helvetica"))

ROH.het<-ggplot(hetdata, aes(x=KB,y=OHET))+
         geom_point(alpha=0.4, size=2)+
         ggtitle(paste0("Heterozygosity vs long runs","\n", "of homozygosity (ROH)"))+
         theme_bw()+
         xlab("Total length of ROH (kb)")+ylab("Heterozygosity")+
         geom_hline(yintercept=mean(hetdata$OHET)-4*sd(hetdata$OHET),color="red")+
         geom_hline(yintercept=mean(hetdata$OHET)+4*sd(hetdata$OHET),color="red")+
         theme(text=element_text(size=10, family="Helvetica"))

###print plots
hetero.density.file <- file.path(output, "04_heterozygosity.tiff")
tiff(hetero.density.file,  
     width = 2500, height = 1500, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(CR.het, ROH.het, ncol=2, bottom=plot.conclude)
dev.off()





