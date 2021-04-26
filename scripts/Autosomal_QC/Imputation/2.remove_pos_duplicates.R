

library(optparse)
library(data.table)
library(stringr)
library(tidyverse)

#
#opt<-list()
#opt$input<-"/groups/umcg-ugli/tmp04/projects/merged_general_QC/8_second_QC_iteration/10_pre_imputation/X_pre_impute"

#########################################################################################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="file for duplicated", metavar="character"),
  
  make_option(c("-c", "--chr"), type="character", default="1", 
              help="chromosome", metavar="character")
);
opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#########################################################################################################


### list position duplicates
chromosomes<-opt$chr
wkdir<-file.path(opt$input)
extractvector<-c()
for (chrfile in chromosomes) {
  
  extractfile<-fread(paste0(wkdir,"/chr_",chrfile,".bim"),data.table = F,header = F)
  dup.pos<-extractfile[which(duplicated(extractfile$V4)),"V4"]
  dup.vars<-extractfile[extractfile$V4 %in% dup.pos, "V2"]
  extractvector<-c(extractvector,dup.vars)
}
write.table(extractvector,paste0(wkdir,"/dupvars.remove"),quote = F,row.names = F,col.names = F)



