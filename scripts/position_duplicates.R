###################################
### Duplicate markesr removal by Call rate
### date: 15-02-2019
### version: 0.01
### authors: EL
###################################

## example run 

#RScript Het_autosomeQC.R -i ***/Autosome_QC
#  -o repout

#########################################################################################################
#########################################################################################################
### Functions and libraries
#########################################################################################################
#########################################################################################################

###packages


library(ggplot2)
library(data.table)
library(dplyr)
library(optparse)

###test

#opt<-list()
#opt$input<-"C:/VirtualM/Analysis2_QC"
#opt$out<-"C:/VirtualM/Analysis2_QC"

##Arguments
#########################################################################################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input path ", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Output path to save report", metavar="character")
); 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied as input path", call.=FALSE)
}

if (is.null(opt$out)){opt$out<-opt$input}

################### Main ###########################################################

autosomal<-c(1:22,"XY")

for (chr in autosomal)
{
 

##read chromosome files  
dup.file <- file.path(opt$input,paste0("chr_",chr,".bim"))
frq.file<-file.path(opt$input,paste0("chr_",chr,".lmiss"))

#create dataframes
bimfile<-fread(dup.file,header=F)
colnames(bimfile)<-c("CHR", "SNP", "CM","POS","A1", "A2")
frqfile<-fread(frq.file, header=T)
#add call rate to the markers data
bimfile<-bimfile%>%mutate(NAME=paste(CHR,":",POS,"_",A1,"_",A2))
datfile<-left_join(bimfile, frqfile,by="SNP")
#make Exclude list with all duplicates with less call rate
datfile<-datfile%>%group_by(NAME)%>%mutate(keepCR=min(F_MISS))
excl.less.duplicates<-datfile[datfile$F_MISS>datfile$keepCR,2]

#make list of duplicates with the same call rate
undup<-datfile %>%filter(!SNP %in% excl.less.duplicates$SNP)
excl.cr.duplicates<-undup[duplicated(undup$NAME),2]
#make list of all the duplicates
excl.duplicates<-c(excl.less.duplicates[[1]],excl.cr.duplicates[[1]])

excl.file <- file.path(opt$out,paste0("chr_",chr,".excl.duplicates"))
##output files
write.table(excl.duplicates,excl.file, quote=F,row.names = F,col.names = F)

}
