###################################
### Heterozygosity filtering and report for autosomal QC
### date: 12-12-2018
### version: 0.01
### authors: EL
###################################

## example run 

#RScript Het_autosomeQC.R -i ***/autosome_QC
#  -o repout
# opt<-list()
# opt$input<-"/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC/4_Het" 
# opt$input<-"/groups/umcg-aad/tmp04/umcg-elopera/test/second_QC_iteration/4_Het" 
# opt$output<-"/groups/umcg-aad/tmp04/umcg-elopera/test/second_QC_iteration/4_Het" 

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
library(ggrepel)
library(data.table)

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
if (is.null(opt$out)){opt$out<-opt$input}


##opt<-list()
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
crdata$IID<-as.character(crdata$IID)
crdata$FID<-as.character(crdata$FID)
crdata<-crdata%>% group_by(IID,FID)%>%summarise(F_miss=mean(F_MISS))

##unify data in a  single database
hetdata$IID<-as.character(hetdata$IID)
hetdata$FID<-as.character(hetdata$FID)
rohdata$IID<-as.character(rohdata$IID)
rohdata$FID<-as.character(rohdata$FID)
hetdata<-right_join(hetdata,rohdata,by=c("IID","FID"))
hetdata<-inner_join(hetdata,crdata,by=c("IID","FID"))
##make filters to exclude and include data



M<-lm(OHET~KB,data=hetdata)
hetdata$outlier <- ifelse(hetdata$OHET - (coef(M)[1] + coef(M)[2]*hetdata$KB) 
                          < -766.7772|
                            hetdata$OHET>39263.88,
                          "Excluded" ,"Included")


Excluded<-hetdata[hetdata$outlier=="Excluded",c(1,2)]
Included<-hetdata[hetdata$outlier=="Included",c(1,2)]
##make files for included/excluded individuals
excl.file <- file.path(opt$input,"Excluded.het")
incl.file <- file.path(opt$input,"Included.het")
write.table(Excluded,excl.file, quote=F,row.names = F)
write.table(Included,incl.file, quote=F,row.names = F)

###plotting
excluded<-nrow(Excluded)
hetdata$excl.label<-
  ifelse(hetdata$outlier=="Excluded",sapply(hetdata$IID , FUN= function(x){substring(x, first = nchar(x)-6, last = nchar(x))}),NA)


plot.conclude<-paste0(excluded," ","Samples excluded")
combine.plot.title <- paste0(opt$res," Heterozygosity", "\n", " ", date())

CR.het<-ggplot(hetdata, aes(x=log(F_miss),y=OHET))+
             geom_point(aes(color=outlier),alpha=0.4, size=2)+
             ggtitle("Heterozygosity distribution by sample")+
             scale_color_manual(values  = c("blue","black"))+
  guides(colour = guide_legend(override.aes = list(alpha=1)))+
             theme_bw()+
             xlab("log(missing)")+ylab("Heterozygosity")+
             geom_hline(yintercept=39263.88,color="red")+
             theme(text=element_text(size=10, family="Helvetica"))


ROH.het<-ggplot(hetdata, aes(x=KB,y=OHET))+
  geom_point(aes(color=outlier),alpha=0.4, size=2)+
  ggtitle(paste0("Heterozygosity vs ROH"))+
  theme_bw()+
  geom_abline(slope = M$coefficients[2], 
              intercept = M$coefficients[1]-766.7772,colour="red",lwd=1.00)+
  scale_color_manual(values  = c("blue","black"))+
  guides(colour = guide_legend(override.aes = list(alpha=1)))+
  xlab("Total length of ROH (kb)")+ylab("Heterozygosity")+
  geom_hline(yintercept=39263.88,color="red",lwd=1.00)+
  #geom_hline(yintercept=mean(hetdata$OHET)-4*sd(hetdata$OHET),color="darkgreen",lwd=1.00)+
  theme(text=element_text(size=10, family="Helvetica"))

if(excluded<20){
  ROH.het<- ROH.het +  geom_text_repel(aes(label=excl.label),size=2.5)
  }

###print plots
hetero.density.file <- file.path(output, "04.Het.scatter.tiff")
tiff(hetero.density.file,  width = 3000, height = 1500, 
    units = "px", res = 300, compression = "lzw")
grid.arrange(CR.het, ROH.het, ncol=2,top= combine.plot.title,bottom=plot.conclude)
dev.off()


