###################################
### Duplicate samples removal by Call rate
### date: 15-02-2019
### version: 0.01
### authors: EL - RAG
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
## local test
## opt<-list()
## opt$input<-"~/Tasks/testdata/"
## opt$out<-"~/Tasks/testdata/"
## cluster test
## opt$wd<-"/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC/5_Relatedness/proc"
## opt$ref<-"/groups/umcg-aad/tmp04/umcg-elopera/concordanceCheck/ugli.final.pairing.concordanceCheck.txt"
## chr=1
##Arguments

codextract <- function(x) {
  x <- as.character(x)
  return(substr(x,nchar(x)-9,nchar(x)))
}

#########################################################################################################
option_list = list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="Input path ", metavar="character"),
  
  make_option(c("-r", "--ref"), type="character", default=NULL, 
              help="Output path to save report", metavar="character")
); 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$wd)){
  print_help(opt_parser)
  stop("At least one argument must be supplied as input path", call.=FALSE)
}

if (is.null(opt$out)){opt$out<-opt$input}

################### Main ###########################################################


##read chromosome files  
dup.file<-fread(opt$ref,data.table=F,header=T)
frq_file<-file.path(opt$wd,"full_autosomal_rel.temp.imiss")
frq.file<-fread(frq_file,data.table=F,header=T)


## name columns on dataframe to match

colnames(dup.file)<-c("PSEUDOID", "Sample_ID", "clean_ID","GWASID","GONL_ID")
frq.file$clean_ID<-sapply(frq.file$IID,FUN=codextract)

##JOIN THE DATAFRAMES
datfile<-left_join(dup.file, frq.file,by="clean_ID")

#make Exclude list with all duplicates with less call rate
datfile<-datfile%>%group_by(PSEUDOID)%>%mutate(keepCR=min(F_MISS))%>%as.data.frame()
excl.less.duplicates<-datfile[which(!is.na(datfile$F_MISS) &datfile$F_MISS>datfile$keepCR),c("IID","FID")]

#make list of duplicates with the same call rate
undup<-datfile %>%filter(!IID %in% excl.less.duplicates[,"IID"])
excl.s.duplicates<-undup[duplicated(undup$PSEUDOID),c("IID","FID")]
#make list of all the duplicates
excl.duplicates<-rbind(excl.less.duplicates,excl.s.duplicates)

excl.file <- file.path(opt$wd, "intended.duplicates")
##output files
write.table(excl.duplicates,excl.file, quote=F,row.names = F,col.names = F)


