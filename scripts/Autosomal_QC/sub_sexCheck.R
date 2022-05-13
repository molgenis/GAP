###################################
### Sex check evaluation
### date: 09-04-2019
### version: 0.01
### authors: EL - RAG
###################################
## New
## 04-04-2022
## removed all the code to report concordance and duplicates by plate
## removed all code dependent on UGLI sample name structure
## 04-07-2019
## changed  lines assigning a value to a match, they produced error when match was not found added an ifelse statement instead
## added while loop to make as many plots as needed by every 50 plates
###################################
## New
## 09-04-2019
## Fixed changing colors for tile.plot
###################################

library(tidyverse)
library(data.table)
library(optparse)
library(gridExtra)
library(viridis)
library(RColorBrewer)

#########################################################################################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input path to perform QC", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="genoQC_Report", 
              help="Output path to save report", metavar="character"),
  
  make_option(c("-p", "--phenotypes"), type="character", default=NULL, 
              help="pedigree file with sex info", metavar="character"),
  
  make_option(c("-d", "--dup"), type="character", default=NULL, 
              help="path to the duplicates file", metavar="character")
); 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied as input path", call.=FALSE)
}

output <- as.character(opt$out)
#output <- file.path(output, "plots")
dir.create(paste0(output,"/sex_check/"), recursive = T)
output<-paste0(output,"/sex_check/")

#### These two constants define the size of the plots given the number of plates in each batch.  
tile.plot.height.factor <- 150
bar.plot.width.factor <- 123
#########################################################################################################
### Main
#########################################################################################################
## Read phenotypes, we are currently using the info from the spreadsheets.
phenos <- fread(opt$phenotypes, data.table = FALSE,header=F)
plink.sex <- fread(opt$input, data.table=FALSE)
##harmonize names
phenos$Sample_ID<- gsub(phenos$V2,pattern="_[0-9]",replacement = "")
plink.sex$IID<-gsub(plink.sex$IID,pattern="_[0-9]",replacement = "")
## see concordance
plink.sex$pheno.sex <- phenos$V5[match(plink.sex$IID, phenos$Sample_ID)]
plink.sex$plink.sex <- ifelse(plink.sex$SNPSEX == 0, NA, ifelse(plink.sex$SNPSEX == 1, "M", "F"))
plink.sex$sex.concordance <- plink.sex$SNPSEX == plink.sex$pheno.sex
##calculate results for summary report
nonna<-length(which(!is.na(plink.sex$sex.concordance)))
conc<-length(which(plink.sex$sex.concordance==T))
### flag samples
plink.sex$sex.concordance[which(plink.sex$sex.concordance == "FALSE")] <- "Non concordant"
plink.sex$sex.concordance[which(is.na(plink.sex$pheno.sex))] <- "No data in pedigree file"
plink.sex$sex.concordance[which(plink.sex$sex.concordance == "TRUE")] <- "OK"
plink.sex$sex.concordance[which(plink.sex$SNPSEX==0)] <- "Failed genetic imputation"
report<-table(plink.sex$sex.concordance)
report[5]<-c(nonna)
names(report)[5]<-"samples with sex information from pedigree file"
report<-as.matrix(report)
### write reports
write.table(report,paste0(output,"sex_concordance.rep"),sep='\t',quote = F,col.names = F)
write.table(plink.sex,paste0(output,"all.samples.concordance.txt"),sep='\t',quote = F,row.names = F)


