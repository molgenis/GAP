###################################
### Duplicate samples removal by Call rate
### date: 15-02-2019
### version: 2
### authors: EL - RAG
###################################
### new
### 04-04-2022
### removed handlers for specific sample names

## example run 
## the duplicated sample with the lesser missing rate will be kept
### make sure you use a headless duplicated sample files & that your duplicated samples are named as [samplename]_1, [samplename]_2... etc
#RScript Het_autosomeQC.R -w [working dir] -r [duplicated sample file]

#########################################################################################################
#########################################################################################################
### Functions and libraries
#########################################################################################################
#########################################################################################################

###packages
library(data.table)
library(dplyr)
library(optparse)

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
##read files
##read chromosome files  
dup.file<-fread(opt$ref,data.table=F,header=F)
frq_file<-file.path(opt$wd,"full_autosomal_rel.temp.imiss")
frq.data<-fread(frq_file,data.table=F,header=T)

## create clean (non-duplicated) IDs
frq.data$clean_ID<-gsub(frq.data$IID,pattern="_[0-9]",replacement = "")
dup.file$clean_ID<-gsub(dup.file$V1,pattern="_[0-9]",replacement = "")

##JOIN THE DATAFRAMES
datfile<-left_join( dup.file,frq.data ,by="clean_ID")

#make Exclude list with all duplicates with less call rate
datfile<-datfile%>%group_by(clean_ID)%>%mutate(keepCR=min(F_MISS))%>%as.data.frame()
excl.less.duplicates<-datfile[which(!is.na(datfile$F_MISS) &datfile$F_MISS>datfile$keepCR),c("IID","FID")]

#make list of duplicates with the same call rate
undup<-datfile %>%filter(!IID %in% excl.less.duplicates[,"IID"])
excl.s.duplicates<-undup[duplicated(undup$PSEUDOID),c("IID","FID")]
#make list of all the duplicates
excl.duplicates<-rbind(excl.less.duplicates,excl.s.duplicates)

excl.file <- file.path(opt$wd, "intended.duplicates")
##output files
write.table(excl.duplicates,excl.file, quote=F,row.names = F,col.names = F)


