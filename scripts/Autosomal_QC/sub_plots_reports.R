###################################----Init----####################
# QC report for genotyping data
# date: 04-03-2019
# version: 0.8 (unnoficcial)
# authors: EL - RAG
##_________________________________
### ----Change log----
##_________________________________
##19-08-2019
## added option to dismiss duplicated snps if there are no more
##20-05-2019
## removed relatedness plot function and section
## amplified CR barplot
## changed names of plots to the QCmanual format
## added folder 0_pre for initial files

##04-03-2019
##excluded PCA analysis to an individual Rscript
##added ref-AF plots with outliers report
##added function to plot ref-AF plots
##added MAF plots after PCA Analysis

##21-02-2019
## added duplicate SNPs comment
## corrected fread bugg that caused to misread columns on snps files
## changed input so " >95%" call rate is now just "high"
## smoothed chr density plots
## adapted reading for parsimonius bash script 
## 11-02-2019
## Added PCA plots

##05-02-2019
## Corrected refMAF bug that caused the pipeline to stop
## added percentajes tothe labels of the barplots
## harmonized chromosome nomenclature on chromosome legends
## Added proper format signatures to relatedness plots
## Corrected header reading in snp and sample file

##11-01-2019
## New arguents - project name, sample info sheet 
## New plot showing starting samples and plate numbers 

## Added relatedness analysis figures. 

## Updated plot MAF correlation to properly asses reference and alternative alleles
## included MAF correlation with goNL 

## changed geom_density -> stat_density using position = "identity" and geom = "line 
## typo on plot title.

## 20180125
## Usage full reference MAF table instead of calling each ref dataset independently
## Changed density calculation of scatter plot - > geom_hex()
## Added MAF comparison with EXaC and gnomAD reference cohort 

##_______________________________
####----example run----  
##_______________________________

#RScript genotypeQC.R -i ***/plink_autosome_QC \
#  -o 
#  -r /groups/umcg-wijmenga/tmp04/umcg-raguirre/pln_ugli/eurMAF_1000g/ 


###__________________________________________________________________________##
########### ------Functions and libraries------#############################
##___________________________________________________________________________
##___________________________________________________________________________

library(tidyverse)
#library(ggsci)
library(data.table)
library(grid)
library(optparse)
#library(ggridges)
library(gridExtra)
library(MASS)
library(viridis)
library(readxl)

### Check ref and alt alleles are the same. 
check.ref.alt.alleles <- function(query.alt, query.ref, ref.alt, ref.ref){
  res <-c()
  if(sum(is.na(c(query.alt, query.ref, ref.alt, ref.ref))) >= 1){
    res <- c(NA)
  }
  else if(query.alt == ref.alt & query.ref == ref.ref) {
    res <- "same"
  } else if(query.alt != ref.alt & query.ref != ref.ref & query.alt == ref.ref & query.ref == ref.alt){
    res <- "flip"
  } else if(sum(c(query.alt, query.ref) %in% c("A", "T")) == 2 | sum(c(query.alt, query.ref) %in% c("G", "C")) == 2){
    res <- "ambiguous" 
  } else if(nchar(paste0(query.alt, query.ref)) > 2){
    res <- "more than 2 alleles"
  }else{
    res <- "no match"
  }
  return(res)
}

###MAF outlier plot###
maf_outlier_plot <- function(cohort.bim.dat, 
                             ref.panel.corrected.AF, 
                             cohort.name){
  
  cohort.bim.dat.pro <- cohort.bim.dat[,c("chr", "bp","cohort.maf", ref.panel.corrected.AF)]
  colnames(cohort.bim.dat.pro) <- c("chr", "bp","cohort.maf", "maf.corrected")
  
  ## remove all NA's fron the reference AF
  cohort.bim.dat.pro <- cohort.bim.dat.pro[complete.cases(cohort.bim.dat.pro),]
  ## selecting all points outside of main group. 
  cohort.bim.dat.pro$lm.residuals <-  lm(cohort.bim.dat.pro$cohort.maf ~ 
                                           cohort.bim.dat.pro$maf.corrected)$residuals
  
  cohort.bim.dat.pro$outlier_threshold <- ifelse(abs(cohort.bim.dat.pro$lm.residuals) > 
                                                   sd(cohort.bim.dat.pro$lm.residuals)*4,
                                                 "Outlier", "OK")
  
  outlier.perc <- format(round(
    length(which(cohort.bim.dat.pro$outlier_threshold == "Outlier")) / 
      sum(!is.na(cohort.bim.dat.pro$outlier_threshold))*100, 
    digits = 2))
  
  hla.cohort.bim.dat <- cohort.bim.dat.pro[which(cohort.bim.dat.pro$chr == 6 & 
                                                   cohort.bim.dat.pro$bp > 28477797 & cohort.bim.dat.pro$bp < 35000000),]
  hla.outlier.perc <- format(round(
    length(which(hla.cohort.bim.dat$outlier_threshold == "Outlier")) / 
      sum(!is.na(cohort.bim.dat.pro$outlier_threshold))*100, 
    digits = 4))
  
  
  main.title <- paste0("Outlier = 4xSD(residuals(cohort.maf~ref.af))", "\n",
                       "Total number of markers=", 
                       format(sum(!is.na(cohort.bim.dat.pro$outlier_threshold)), scientific = T, big.mark=","), "\n",
                       "Outliers=", outlier.perc, "%\n",
                       "Outliers in HLA locus=", hla.outlier.perc, "%")
  
  corrected.final.v2 <- ggplot(cohort.bim.dat.pro[which(!is.na(cohort.bim.dat.pro$outlier_threshold)),], 
                               aes(y=cohort.maf , x= maf.corrected))+
    geom_point(size=0.75, alpha= 0.5, aes(color=outlier_threshold))+
    scale_color_brewer(palette = "Dark2", name="")+
    geom_abline(color="black", slope=1, intercept = sd(cohort.bim.dat.pro$lm.residuals)*4, lwd=1.25, linetype=3)+
    geom_abline(color="black", slope=1, intercept = -sd(cohort.bim.dat.pro$lm.residuals)*4, lwd=1.25, linetype=3)+
    theme_classic()+
    ylab("Cohort MAF")+
    xlab(paste0(cohort.name, " AF"))+
    ggtitle(main.title)+
    theme(text=element_text(size=10, family = "Helvetica") ,legend.position = "bottom", legend.text = element_text(size = 10))
  
  return(corrected.final.v2)
}




##Arguments
### Set up arguments for testing 
# opt <- list()
# opt$input <- "/groups/umcg-aad/tmp04/umcg-elopera/Autosome_QC_B_part8"
# opt$out <- "/groups/umcg-aad/tmp04/umcg-elopera/Autosome_QC_B_part8/plots"
# opt$name <- "test01"
# opt$sampleinfo <- "/groups/umcg-aad/tmp04/projects/IBD_part/run01/jobs/IBD_part.csv"
# opt$ref <- "/groups/umcg-wijmenga/tmp04/umcg-raguirre/pln_ugli/af.ref.data.txt"

###############################----Parsing arguments Options-----#############################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input path to perform QC", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="genoQC_Report", 
              help="Output path to save report", metavar="character"),
  
  make_option(c("-n", "--name"), type="character", default="genoQC_Report", 
              help="Name of the project", metavar="character"), 
  
  make_option(c("-s", "--sampleinfo"), type="character", default=NULL, 
              help="path to sample info spreadsheet", metavar="character"), 
  
  make_option(c("-r", "--ref"), type="character", default=NULL, 
              help="Path to the reference Allelic Frequency of SNPs", metavar="character")
); 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied as input path", call.=FALSE)
}
#
output <- as.character(opt$out)
#_________________________________________________________________________________________________________
############################################ ------Main------#############################################
#_________________________________________________________________________________________________________

#####--Plate information sheet---####### 
# Date and project name 
# Table of all plates used in the analysis 
# Table showing number of starting samples 
# Cohort within freeze? 

# testing
 #opt <- list(); opt$sampleinfo <- "/groups/umcg-aad/tmp04/projects/IBD_part/run01/jobs/IBD_part.csv"

if(is.null(opt$sampleinfo) == FALSE){
  sample.info <- fread(as.character(opt$sampleinfo), data.table = FALSE)
  plates.in.analysis <- as.character(unique(sample.info$Sample_Plate))
  
  if(length(plates.in.analysis) %% 2 == 1){
    plates.in.analysis <-   c(plates.in.analysis, "")
  }
  
  plates.in.analysis.matrix <- as.data.frame(matrix(plates.in.analysis, ncol= 4))
  #plates.in.analysis.matrix <- as.data.frame( matrix(paste0("test", sample(1:200, size = 30, replace = FALSE)),nrow= 4))
  plate.table <- tableGrob(plates.in.analysis.matrix, ttheme_minimal(base_size= 8, base_family= "Helvetica"),
                           rows=NULL, cols = NULL)
  
  cover.plot.file <- file.path(output, "01_coverPlot.tiff")
  cover.plot.title <- paste0(opt$name, "\n", date(), "\n", "List of plates used: ")
  
  tiff(cover.plot.file, width = 2000, height = 2500, 
       units = "px", res = 300, compression = "lzw")
  
  grid.arrange(plate.table, top= cover.plot.title)
  
  dev.off()
}

############# Barplots of call rates########

## loading samples and SNPs
samples_all.file <- file.path(opt$input,"0_pre/full.ind")
if (file.exists(samples_all.file)){pre="0_pre/"} else{pre=""}
samples_all.file <- file.path(opt$input,pre,"full.ind")
samples_all <- fread(samples_all.file, data.table = FALSE, header=F)

snps_all.file <- file.path(opt$input,pre,"full.snps")
snps_all <- fread(snps_all.file, sep=" ", header=F)

samples_80.file <- file.path(opt$input, paste0("1_CR80/","incl80.samples"))
samples_80 <- fread(samples_80.file, data.table = FALSE, header=F)

snps_80.file <- file.path(opt$input, paste0("1_CR80/","incl80.vars"))
snps_80 <- fread(snps_80.file,sep=" ", header=F)

samples_high.file <- file.path(opt$input, paste0("2_CR_high/", "inclhigh.samples"))
samples_high<-fread(samples_high.file, data.table = FALSE, header=F)
colnames(samples_high)<-c("IID")
samples_high$IID<-as.character(samples_high$IID)

snps_high.file <- file.path(opt$input, paste0("2_CR_high/","inclhigh.vars"))
snps_high <- fread(snps_high.file, data.table = FALSE,sep=" ", header=F)

dups.file<- file.path(opt$input,pre,"extr.dups")
if(file.info(dups.file)$size==0|file.exists(dups.file)==FALSE ){
  dups<-data.frame(ID=as.numeric(character()))
  
}  else {dups<-fread(dups.file, data.table = FALSE,sep=" ", header=F)}
   

###group data to make the barplots
samples <- c("All"=nrow(samples_all),"CR>80"=nrow(samples_80), "CR>99"=nrow(samples_high))
samples_lab<-c(paste0(samples[1],"\n","(100%)"),
               paste0(samples[2],"\n","(",round(samples[2]*100/samples[1],1),"%",")") ,
               paste0(samples[3],"\n","(",round(samples[3]*100/samples[1],1),"%",")"))
snps <- c("All"=nrow(snps_all),"CR>80"=nrow(snps_80), "CR>99"=nrow(snps_high) )
snps_lab<-c(paste0(snps[1],"\n","(100%)"),
               paste0(snps[2],"\n","(",round(snps[2]*100/snps[1],1),"%",")") ,
               paste0(snps[3],"\n","(",round(snps[3]*100/snps[1],1),"%",")"))



incl <- data.frame(samples, snps,"group"= names(snps))

#variables to calculate the ylim of the barplots
z_samples <- mean(c(samples["All"]-samples["CR>80"], samples["CR>80"]-samples["CR>99"]))+1
z_snps <- mean(c(snps["All"]-snps["CR>80"], snps["CR>80"]-snps["CR>99"]))+1

##barplot for samples
barplot.samples <- ggplot(incl, aes(group,samples))+
  geom_bar(stat = "identity", aes(colour=group, fill=group), position = position_dodge(width = 0.5))+
  coord_cartesian(ylim = c(round(samples["CR>99"]-z_samples-4,digits=0),
                           round(samples["All"]+2,digits=0)))+
  ylab("Number of SNPs")+
  ylab("Number of Samples")+
  theme_classic()+
  geom_text(aes(label=samples_lab), vjust=0, size=2.5)+
  theme(text=element_text(size=10, family = 'Helvetica'), legend.position = "none")

##barplot for snps
barplot.snps <- ggplot(incl, aes(group,snps)) +
  geom_bar(stat = "identity", aes(colour=group, fill=group), 
           position = position_dodge(width = 0.5))+
  coord_cartesian(ylim = c(round(snps["CR>99"]-z_snps-400,digits=0),
                           round(snps["All"]*1.0025,digits=0)))+
  ylab("Number of markers")+
  theme_classic()+
  geom_text(aes(label=snps_lab), vjust=0,size=2.5)+
  theme(text=element_text(size=10, family = 'Helvetica'), legend.position = "none")


samples.snps.barplot.file <- file.path(output, "02.CR.barplot.tiff")
samples.snps.barplot.title <- paste0("Number of samples and markers at different call rate (CR) thresholds", 
                                     "\n",
                                     "(from autosomal and pseudo-autosomal chromosomes)", 
                                     "\n", 
                                     opt$name, " ", date())
dups.comment<-paste0("Y axis does not start from zero","\n", 
                     "The total number of markers excludes ", nrow(dups),
                     " markers that were found to be duplicated by position")

tiff(samples.snps.barplot.file,  
     width = 1500, height = 1500, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(barplot.samples, barplot.snps, 
             ncol=2, top=textGrob(samples.snps.barplot.title, gp=gpar(fontsize=9,font=8)),
             bottom=textGrob(dups.comment, gp=gpar(fontsize=6,font=8)))
dev.off()


############### Final CR distribution#######
##look for filepath of CR80
#create list of filenames
cr.files.path <- file.path(opt$input, "1_CR80/")
cr.files <- list.files(cr.files.path, pattern = "\\.2\\.imiss$", full.names = TRUE)

#create list of files
cr.dat.list <- lapply(cr.files, fread, data.table=FALSE)
names(cr.dat.list) <- gsub(basename(cr.files), replacement = '', pattern = '.imiss')
#select columns to be used
cr.dat.list <- lapply(cr.dat.list, function(x){x[,c("IID", "N_MISS", "N_GENO","F_MISS")]})
#add chr column to data list
cr.dat.list<- lapply(names(cr.dat.list), 
                     function(chr){newdf <- cbind(cr.dat.list[[chr]], CHR=chr) ; newdf})
##set right type of variables
cr.dat.list <- lapply(cr.dat.list, function(df){mutate_at(df, .vars = c("IID","CHR"), as.character)})
cr.dat.list <- lapply(cr.dat.list, function(df){mutate_at(df, .vars = c("N_MISS", "N_GENO","F_MISS"), as.numeric)})
##merge chr files into one database (per chromosome database)
cr.dat <- bind_rows(cr.dat.list)
#extract the samples below high CR
cr.dat <-semi_join(cr.dat,samples_high,by="IID")
##get the genomic CR data base
cr.datG<-cr.dat %>% group_by(IID) %>%
  summarise(N_MISST=sum(N_MISS),
            N_GENOT=sum(N_GENO),
            Call_rate=((N_GENOT-N_MISST)/N_GENOT))

##Plot CR distribution
# define levels of CHR to plot CHRs in order. 
cr.dat$CHR <- gsub(cr.dat$CHR, pattern = "chr_", replacement ="" )
cr.dat$CHR <- gsub(cr.dat$CHR, pattern = "\\..*", replacement ="" )
cr.dat$CHR <- factor(cr.dat$CHR, levels=c(1:22, "XY"))

chr.labels <- c(1:22,25)
names(chr.labels) <-  c(1:22,'XY')

#by chromosome 
cr.by.chr.dens <- ggplot(cr.dat, aes(x=1-F_MISS))+
  stat_density(aes(colour= CHR),geom= "line", position= "identity",adjust=5)+
  facet_wrap(~CHR, scales = 'free')+
  scale_color_manual(labels=names(chr.labels), values=rainbow(23))+
  xlab("Call rate")+
  theme_bw()+
  theme(text = element_text(size=10, family='Helvetica'), legend.position = "none")
##by chromosome merged
cr.merged.chrs.dens <- ggplot(cr.dat, aes(x=1-F_MISS))+
  stat_density(aes(colour= CHR), geom= "line", position= "identity",adjust=5)+
  scale_color_manual(labels=names(chr.labels), values=rainbow(23))+
  xlab("Call rate")+
  theme_bw()+
  theme(text = element_text(size=10, family='Helvetica'))
#by full autosomes+XY
cr.overall <- ggplot(cr.datG,aes(x=Call_rate))+
  stat_density(geom= "line", position= "identity")+
  xlab("Call rate")+
  theme_bw()+
  theme(text=element_text(size=10,family='Helvetica'))


cr.plots.file <- file.path(output, "02.CR.density.a.tiff")
cr.plot.tilte <- paste0("Distribution of call rate (CR) from autosomal and pseudo-autosomal chromosomes", "\n", opt$name, " ", date())
tiff(cr.plots.file,  
     width = 3000, height = 1300, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(cr.merged.chrs.dens, cr.overall, nrow=1, top= cr.plot.tilte)
dev.off()


cr.plots.by.chr.file <- file.path(output, "02.CR.density.b.tiff")
cr.plot.grid.tilte <- paste0("Distribution of call rate from autosomal and pseudo-autosomal chromosomes", 
                             "\n", opt$name, " ", date())
tiff(cr.plots.by.chr.file,  
     width = 2800, height = 3200, 
     units = "px", res = 300, compression = "lzw")
plot(cr.by.chr.dens+ggtitle(cr.plot.grid.tilte))
dev.off()

#____________________________________________________________________
############### MAF and HW distributions across all chromosomes#######
#____________________________________________________________________
### outputhigh
maf.hw.path <- file.path(opt$input, "/3_MAF_HWE/")
#maf.hw.path <- "/groups/umcg-wijmenga/tmp04/umcg-raguirre/pln_ugli/genotype"
#maf.hw.path <- "/Users/raulaguirre/PLN_UGLI/Data/testGenotypeData"

maf.files <- list.files(maf.hw.path, pattern = ".frq$", full.names = TRUE)
hw.files <- list.files(maf.hw.path, pattern = ".hwe$", full.names = TRUE)

## Distribution of MAF across all SNPs with sufficient call rate 
maf.dat.list <- lapply(maf.files, fread, data.table=FALSE)
maf.dat.list <- lapply(maf.dat.list, function(x){x[,c("CHR", "SNP", "MAF")]})
maf.dat.list <- lapply(maf.dat.list, function(df){mutate_at(df, .vars = c("CHR"), as.character)})
maf.dat.list <- lapply(maf.dat.list, function(df){mutate_at(df, .vars = c("MAF"), as.numeric)})
maf.dat <- bind_rows(maf.dat.list)
rm(maf.dat.list)

## Distribution of HW across all SNPs with sufficient call rate 
hw.dat.list <- lapply(hw.files, fread, data.table=FALSE)
hw.dat.list <- lapply(hw.dat.list, function(x){x[,c("CHR", "SNP", "P")]})
hw.dat.list <- lapply(hw.dat.list, function(df){mutate_at(df, .vars = c("CHR"), as.character)})
hw.dat.list <- lapply(hw.dat.list, function(df){mutate_at(df, .vars = c("P"), as.numeric)})
hw.dat <- bind_rows(hw.dat.list)
rm(hw.dat.list)

## Merge HW and MAF info 
hw.dat$MAF <- maf.dat$MAF[match(hw.dat$SNP ,maf.dat$SNP)] 

## Define levels in CHR to plot them in numerical order
hw.dat$CHR <- factor(hw.dat$CHR, levels= c(1:22, 25))
maf.dat$CHR <- factor(maf.dat$CHR, levels= c(1:22, 25))



maf.dist.plot.chr <- ggplot(maf.dat, aes(x=MAF))+
  stat_density(aes(color= CHR), position= "identity", geom= "line")+
  ggtitle("MAF distribution per chromosome")+
  scale_color_manual(labels=names(chr.labels), values=rainbow(23))+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

maf.dist.plot.all <- ggplot(maf.dat, aes(x=MAF))+
  stat_density(color="black",  position= "identity", geom= "line")+
  ggtitle("MAF distribution")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))


hw.dist.plot.chr <- ggplot(hw.dat, aes(x=-log10(P)))+
  stat_density(aes(color= CHR), position= "identity", geom= "line")+
  scale_color_manual(labels=names(chr.labels), values=rainbow(23))+
  ggtitle("HW pVal distribution  \n per chromosome")+
  xlab("-log10(HW-P)")+
  xlim(c(0, 20))+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

hw.maf.dist.plot.chr <- ggplot(hw.dat[which(hw.dat$MAF > 0.01),], aes(x=-log10(P)))+
  stat_density(aes(color= CHR), position= "identity", geom= "line")+
  geom_vline(xintercept = 6)+
  ggtitle("HW pVal distribution of \n markers with a MAF > 0.01")+
  scale_color_manual(labels=names(chr.labels), values=rainbow(23))+
  xlim(c(0, 20))+
  xlab("-log10(HW-P)")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

hw.maf.dist.plot.all <- ggplot(hw.dat[which(hw.dat$MAF > 0.01),], aes(x=-log10(P)))+
  stat_density(position= "identity", geom= "line")+
  geom_vline(xintercept = 6)+
  xlim(c(0, 20))+
  ggtitle("HW pVal distribution of \n markers with a MAF > 0.01")+
  xlab("-log10(HW-P)")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))


# layout matrix for grid.arrange
lay <- rbind(c(1,1,1,2,2,2),
             c(3,3,4,4,5,5)
)
#concluding information
excluded<-fread(paste0(maf.hw.path, "excl_HW.snps"),data.table = FALSE,sep=" ", header=F)


### printing out the maf.hw distributions
maf.hw.title <- paste0("MAF and HW distributions (from autosomal and pseudo-autosomal chromosomes)","\n", 
                       opt$name, " ", date())
maf.hw.conclude<-paste0(nrow(excluded)," markers excluded ", snps[3]-nrow(excluded), " remaining.")
maf.hw.tiff.file <- file.path(output, "03.MAFHW.density.tiff")
tiff(maf.hw.tiff.file, width = 3000, height = 2000, units = "px", res = 300, compression = "lzw")
grid.arrange(maf.dist.plot.chr, maf.dist.plot.all, 
             hw.dist.plot.chr, hw.maf.dist.plot.chr, 
             hw.maf.dist.plot.all, 
             layout_matrix = lay, top=maf.hw.title, 
             bottom=textGrob(maf.hw.conclude, gp=gpar(fontsize=6,font=8)))
dev.off()


#_____________________________________________________________________________________________
############### Concordance of cohort MAF with Allelic frequencies in different reference cohorts########
#_________________________________________________________________________________________________________
###
# Read reference AF from pre-processed file

ref.maf.dat <- fread(opt$ref ,data.table = FALSE)

###
# Read bim files from cohort 
# opt <-list(); opt$input <- "/groups/umcg-aad/tmp04/umcg-elopera/plink_autosome_QC"
cohort.bim.path <- file.path(opt$input, "2_CR_high")
cohort.bim.files <- list.files(cohort.bim.path, pattern = ".bim$", full.names = TRUE)
cohort.bim.list <- lapply(cohort.bim.files, fread, data.table=FALSE)
cohort.bim.dat <- bind_rows(cohort.bim.list)
rm(cohort.bim.list)
colnames(cohort.bim.dat) <- c("chr", "snp", "morgan", "bp", "A1", "A2")
# generate chr:pos IDs
cohort.bim.dat$id <- paste0(cohort.bim.dat[,1],":",cohort.bim.dat[,4])

#################################
#### add the MAF to the bim file. 
cohort.bim.dat$cohort.maf <- maf.dat$MAF[match(as.character(cohort.bim.dat$snp), as.character(maf.dat$SNP))]
#### Remove all SNPS which have a MAF of 0 in the cohort
zero.maf.index <- which(cohort.bim.dat$cohort.maf == 0)
if (length(zero.maf.index) >= 1){
  cohort.bim.dat <- cohort.bim.dat[-zero.maf.index,]
}

one.maf.index <- which(cohort.bim.dat$cohort.maf == 1)
if (length(one.maf.index) >= 1){
  cohort.bim.dat <- cohort.bim.dat[-one.maf.index,]
}

#################################
#### add the HW pvalue - remove all SNPS with a pValue <= 0.000005
cohort.bim.dat$cohort.hw.p <- hw.dat$P[match(as.character(cohort.bim.dat$snp), as.character(hw.dat$SNP))]

hw.index <- which(cohort.bim.dat$cohort.hw.p <= 0.000005)
if (length(hw.index) >= 1){
  cohort.bim.dat <- cohort.bim.dat[-hw.index,]
}


################################ 1000G
snp.index <- match(cohort.bim.dat$id, ref.maf.dat$id)

### I would no longer remove markers that are not present in 1000G
#cohort.bim.dat <- cohort.bim.dat[which(cohort.bim.dat$id %in% ref.maf.dat$id),]

cohort.bim.dat$ref.alt.allele <- ref.maf.dat$ref.alt.allele[snp.index]
cohort.bim.dat$ref.ref.allele <- ref.maf.dat$ref.ref.allele[snp.index]
cohort.bim.dat$ref.maf <- ref.maf.dat$ref.maf[snp.index]

## check that the major and minor allele are the same. 
# check allelic concordance 
allele.matrix <- cbind(query.alt= cohort.bim.dat$A1, 
                       query.ref= cohort.bim.dat$A2, 
                       ref.alt= cohort.bim.dat$ref.alt.allele, 
                       ref.ref= cohort.bim.dat$ref.ref.allele)

cohort.bim.dat$allele.check <- apply(allele.matrix, 1, FUN = function(x){check.ref.alt.alleles(query.alt= x[1], query.ref= x[2],  ref.alt= x[3], ref.ref= x[4])})

cohort.bim.dat$ref.maf.corrected <- cohort.bim.dat$ref.maf
cohort.bim.dat$ref.maf.corrected[which(cohort.bim.dat$allele.check == "flip")] <- (1-as.numeric(cohort.bim.dat$ref.maf[which(cohort.bim.dat$allele.check == "flip")]))
cohort.bim.dat$ref.maf.corrected <- as.numeric(cohort.bim.dat$ref.maf.corrected)
cohort.bim.dat$cohort.maf <- as.numeric(cohort.bim.dat$cohort.maf)


################################ goNL 
# opt$goNLref <- "/groups/umcg-wijmenga/tmp04/umcg-raguirre/pln_ugli/gonl_MAF/goNL_MAF.edR.INFO"

cohort.bim.dat$goNl.maf <- ref.maf.dat$goNl.maf[snp.index]
cohort.bim.dat$goNl.alt <- ref.maf.dat$goNl.alt[snp.index]
cohort.bim.dat$goNl.ref <- ref.maf.dat$goNl.ref[snp.index]

goNl.allele.matrix <- cbind(query.alt=cohort.bim.dat$A1, query.ref=cohort.bim.dat$A2, 
                            ref.alt=cohort.bim.dat$goNl.alt, ref.ref=cohort.bim.dat$goNl.ref)

cohort.bim.dat$goNl.allele.check <- apply(goNl.allele.matrix, 1, 
                                          FUN = function(x){check.ref.alt.alleles(
                                            query.alt= x[1], query.ref= x[2],  ref.alt= x[3], ref.ref= x[4])
                                          })

cohort.bim.dat$goNl.maf.corrected <- cohort.bim.dat$goNl.maf
cohort.bim.dat$goNl.maf.corrected[which(cohort.bim.dat$goNl.allele.check == "flip")] <- (1-as.numeric(cohort.bim.dat$goNl.maf[which(cohort.bim.dat$goNl.allele.check == "flip")]))
cohort.bim.dat$goNl.maf.corrected <- as.numeric(cohort.bim.dat$goNl.maf.corrected)


################################ EXaC 


cohort.bim.dat$exac.maf <- ref.maf.dat$exac.af[snp.index]
cohort.bim.dat$exac.maf.nfe <- ref.maf.dat$exac.af.NFE[snp.index]
cohort.bim.dat$exac.alt <- ref.maf.dat$exac.alt[snp.index]
cohort.bim.dat$exac.ref <- ref.maf.dat$exac.ref[snp.index]

exac.allele.matrix <- cbind(query.alt=cohort.bim.dat$A1, query.ref=cohort.bim.dat$A2, 
                            ref.alt=cohort.bim.dat$exac.alt, ref.ref=cohort.bim.dat$exac.ref)

cohort.bim.dat$exac.allele.check <- apply(exac.allele.matrix, 1, 
                                          FUN = function(x){check.ref.alt.alleles(
                                            query.alt= x[1], query.ref= x[2],  ref.alt= x[3], ref.ref= x[4])
                                          })

cohort.bim.dat$exac.maf.corrected <- cohort.bim.dat$exac.maf
cohort.bim.dat$exac.maf.corrected[which(cohort.bim.dat$exac.allele.check == "flip")] <- (1-as.numeric(cohort.bim.dat$exac.maf[which(cohort.bim.dat$exac.allele.check == "flip")]))
cohort.bim.dat$exac.maf.corrected <- as.numeric(cohort.bim.dat$exac.maf.corrected)

cohort.bim.dat$exac.maf.corrected.nfe <- cohort.bim.dat$exac.maf.nfe
cohort.bim.dat$exac.maf.corrected.nfe[which(cohort.bim.dat$exac.allele.check == "flip")] <- (1-as.numeric(cohort.bim.dat$exac.maf.corrected.nfe[which(cohort.bim.dat$exac.allele.check == "flip")]))
cohort.bim.dat$exac.maf.corrected.nfe <- as.numeric(cohort.bim.dat$exac.maf.corrected.nfe)

################################ gnomad

cohort.bim.dat$gnomad.maf <- ref.maf.dat$gnomad.af[snp.index]
cohort.bim.dat$gnomad.maf.nfe <- ref.maf.dat$gnomad.af.NFE[snp.index]
cohort.bim.dat$gnomad.alt <- ref.maf.dat$gnomad.alt[snp.index]
cohort.bim.dat$gnomad.ref <- ref.maf.dat$gnomad.ref[snp.index]

gnomad.allele.matrix <- cbind(query.alt=cohort.bim.dat$A1, query.ref=cohort.bim.dat$A2, 
                              ref.alt=cohort.bim.dat$gnomad.alt, ref.ref=cohort.bim.dat$gnomad.ref)

cohort.bim.dat$gnomad.allele.check <- apply(gnomad.allele.matrix, 1, 
                                            FUN = function(x){check.ref.alt.alleles(
                                              query.alt= x[1], query.ref= x[2],  ref.alt= x[3], ref.ref= x[4])
                                            })

cohort.bim.dat$gnomad.maf.corrected <- cohort.bim.dat$gnomad.maf
cohort.bim.dat$gnomad.maf.corrected[which(cohort.bim.dat$gnomad.allele.check == "flip")] <- (1-as.numeric(cohort.bim.dat$gnomad.maf[which(cohort.bim.dat$gnomad.allele.check == "flip")]))
cohort.bim.dat$gnomad.maf.corrected <- as.numeric(cohort.bim.dat$gnomad.maf.corrected)


cohort.bim.dat$gnomad.maf.corrected.nfe <- cohort.bim.dat$gnomad.maf.nfe
cohort.bim.dat$gnomad.maf.corrected.nfe[which(cohort.bim.dat$gnomad.allele.check == "flip")] <- (1-as.numeric(cohort.bim.dat$gnomad.maf.corrected.nfe[which(cohort.bim.dat$gnomad.allele.check == "flip")]))
cohort.bim.dat$gnomad.maf.corrected.nfe <- as.numeric(cohort.bim.dat$gnomad.maf.corrected.nfe)

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#plots
#1000G
maf.ref.cohort.hex <- ggplot(cohort.bim.dat[which(cohort.bim.dat$allele.check== "same"| cohort.bim.dat$allele.check== "flip"),], 
                             aes(x=ref.maf.corrected , y = cohort.maf))+
  ggtitle("AF from EUR 1000g vs cohort MAF")+
  ylab("MAF at cohort")+
  xlab("AF at EUR 1000G")+
  geom_hex(alpha= 1, bins= 150)+
  scale_fill_viridis(breaks= c(1,10,100,1000,10000), limits=c(1, 10000), 
                     labels=c(1, expression(10^1), expression(10^2), expression(10^3), expression(10^4)),
                     name="Counts", trans="log10")+
  theme_classic()+
  guides(fill=guide_colorbar(barheight = 0.5, barwidth = 5))+
  theme(text=element_text(size=10, family = "Helvetica") ,legend.position = "bottom", legend.text = element_text(size = 10))


### goNL plots
gonl.maf.ref.cohort.hex <- ggplot(cohort.bim.dat[which(cohort.bim.dat$goNl.allele.check== "same"| cohort.bim.dat$goNl.allele.check== "flip"),], 
                                  aes(x=goNl.maf.corrected , y = cohort.maf))+
  ggtitle("AF from GoNL vs cohort MAF")+
  ylab("MAF at cohort")+
  xlab("AF at GoNL")+
  geom_hex(alpha= 1, bins= 150)+
  scale_fill_viridis(breaks= c(1,10,100,1000,10000), limits=c(1, 10000), 
                     labels=c(1, expression(10^1), expression(10^2), expression(10^3), expression(10^4)),
                     name="Counts", trans="log10")+
  theme_classic()+
  guides(fill=guide_colorbar(barheight = 0.5, barwidth = 5))+
  theme(text=element_text(size=10, family = "Helvetica") ,legend.position = "bottom", legend.text = element_text(size = 10))



### exac ALL 
exac.maf.ref.cohort.hex <- ggplot(cohort.bim.dat[which(cohort.bim.dat$exac.allele.check== "same"| cohort.bim.dat$exac.allele.check== "flip"),], 
                                  aes(x=exac.maf.corrected , y = cohort.maf))+
  ggtitle("AF from EXaC vs cohort MAF")+
  ylab("MAF at cohort")+
  xlab("AF at EXaC")+
  geom_hex(alpha= 1, bins= 150)+
  scale_fill_viridis(breaks= c(1,10,100,1000,10000), limits=c(1, 10000), 
                     labels=c(1, expression(10^1), expression(10^2), expression(10^3), expression(10^4)),
                     name="Counts", trans="log10")+
  theme_classic()+
  guides(fill=guide_colorbar(barheight = 0.5, barwidth = 5))+
  theme(text=element_text(size=10, family = "Helvetica") ,legend.position = "bottom", legend.text = element_text(size = 10))

### exac NFE
exac.maf.NFE.ref.cohort.hex <- ggplot(cohort.bim.dat[which(cohort.bim.dat$exac.allele.check== "same"| cohort.bim.dat$exac.allele.check== "flip"),], 
                                      aes(x=exac.maf.corrected.nfe , y = cohort.maf))+
  ggtitle("AF from EXaC vs cohort MAF")+
  ylab("MAF at cohort")+
  xlab("AF at EXaC (NFE popoulation)")+
  geom_hex(alpha= 1, bins= 150)+
  scale_fill_viridis(breaks= c(1,10,100,1000,10000), limits=c(1, 10000), 
                     labels=c(1, expression(10^1), expression(10^2), expression(10^3), expression(10^4)),
                     name="Counts", trans="log10")+
  theme_classic()+
  guides(fill=guide_colorbar(barheight = 0.5, barwidth = 5))+
  theme(text=element_text(size=10, family = "Helvetica") ,legend.position = "bottom", legend.text = element_text(size = 10))


### gnomad 
gnomad.ref.cohort.hex <- ggplot(cohort.bim.dat[which(cohort.bim.dat$gnomad.allele.check== "same"| cohort.bim.dat$gnomad.allele.check== "flip"),], 
                                aes(x=gnomad.maf.corrected , y = cohort.maf))+
  ggtitle("AF from Gnomad vs cohort MAF")+
  ylab("MAF at cohort")+
  xlab("AF at gnomAD")+
  geom_hex(alpha= 1, bins= 150)+
  scale_fill_viridis(breaks= c(1,10,100,1000,10000), limits=c(1, 10000), 
                     labels=c(1, expression(10^1), expression(10^2), expression(10^3), expression(10^4)),
                     name="Counts", trans="log10")+
  theme_classic()+
  guides(fill=guide_colorbar(barheight = 0.5, barwidth = 5))+
  theme(text=element_text(size=10, family = "Helvetica") ,legend.position = "bottom", legend.text = element_text(size = 10))

### gnomad NFE
gnomad.NFE.ref.cohort.hex <- ggplot(cohort.bim.dat[which(cohort.bim.dat$gnomad.allele.check== "same"| cohort.bim.dat$gnomad.allele.check== "flip"),], 
                                    aes(x=gnomad.maf.corrected.nfe , y = cohort.maf))+
  ggtitle("AF from gnomAD vs cohort MAF")+
  ylab("MAF at cohort")+
  xlab("AF at gnomAD (NFE popoulation) ")+
  geom_hex(alpha= 1, bins= 150)+
  scale_fill_viridis(breaks= c(1,10,100,1000,10000), limits=c(1, 10000), 
                     labels=c(1, expression(10^1), expression(10^2), expression(10^3), expression(10^4)),
                     name="Counts", trans="log10")+
  theme_classic()+
  guides(fill=guide_colorbar(barheight = 0.5, barwidth = 5))+
  theme(text=element_text(size=10, family = "Helvetica") ,legend.position = "bottom", legend.text = element_text(size = 10))


### saving out plots. 
combine.plot.title <- paste0("MAF comparison with reference sets AF", "\n", opt$name, " ", date())
maf.ref.cohort.file <- file.path(output, "05.snp_ref.scatter.tiff")

#plotting
tiff(maf.ref.cohort.file,  width = 2000, height = 4000, units = "px", res = 300, compression = "lzw")
grid.arrange(top=combine.plot.title,
             maf.ref.cohort.hex,
             gonl.maf.ref.cohort.hex,
             exac.maf.ref.cohort.hex,
             exac.maf.NFE.ref.cohort.hex,
             gnomad.ref.cohort.hex,
             gnomad.NFE.ref.cohort.hex,
             ncol=2)
dev.off()

# /groups/umcg-wijmenga/tmp04/umcg-raguirre/pln_ugli/af.ref.data.txt
maf_outlier_plot_1000g <- maf_outlier_plot(cohort.bim.dat, ref.panel.corrected.AF = "ref.maf.corrected", cohort.name="1000G")
maf_outlier_plot_goNl <- maf_outlier_plot(cohort.bim.dat, ref.panel.corrected.AF = "goNl.maf.corrected", cohort.name="goNl")
maf_outlier_plot_exac <- maf_outlier_plot(cohort.bim.dat, ref.panel.corrected.AF = "exac.maf.corrected", cohort.name="ExaC")
maf_outlier_plot_exac_nfe <- maf_outlier_plot(cohort.bim.dat, ref.panel.corrected.AF = "exac.maf.corrected.nfe", cohort.name="ExaC-NFE")
maf_outlier_plot_gnomad <- maf_outlier_plot(cohort.bim.dat, ref.panel.corrected.AF = "gnomad.maf.corrected", cohort.name="gNomad")
maf_outlier_plot_gnomad_nfe <- maf_outlier_plot(cohort.bim.dat, ref.panel.corrected.AF = "gnomad.maf.corrected.nfe", cohort.name="gNomad-NFE")

maf.ref.cohort.outlier.file <- file.path(output, "05.out_ref.scatter.tiff") 
tiff(maf.ref.cohort.outlier.file,  width = 2000, height = 4000, units = "px", res = 300, compression = "lzw")
grid.arrange(top=combine.plot.title, 
             maf_outlier_plot_1000g,
             maf_outlier_plot_goNl,
             maf_outlier_plot_exac,
             maf_outlier_plot_exac_nfe,
             maf_outlier_plot_gnomad,
             maf_outlier_plot_gnomad_nfe,
             ncol=2)
dev.off()

####___________________________________________
############### MAF after PCA########

maf.path <- file.path(opt$input, "6_PCA/")
if (file.exists(paste0(maf.path,"chr_1.frq"))){ 

maf.files <- list.files(maf.path, pattern = ".frq$", full.names = TRUE)

## Distribution of MAF across all SNPs with sufficient call rate 
maf.dat.list <- lapply(maf.files, fread, data.table=FALSE)
maf.dat.list <- lapply(maf.dat.list, function(x){x[,c("CHR", "SNP", "MAF")]})
maf.dat.list <- lapply(maf.dat.list, function(df){mutate_at(df, .vars = c("CHR"), as.character)})
maf.dat.list <- lapply(maf.dat.list, function(df){mutate_at(df, .vars = c("MAF"), as.numeric)})
maf.dat <- bind_rows(maf.dat.list)
rm(maf.dat.list)


## Define levels in CHR to plot them in numerical order
maf.dat$CHR <- factor(maf.dat$CHR, levels= c(1:22, 25))

maf.dist.plot.chr <- ggplot(maf.dat, aes(x=MAF))+
  stat_density(aes(color= CHR), position= "identity", geom= "line")+
  ggtitle("MAF distribution per chromosome")+
  scale_color_manual(labels=names(chr.labels), values=rainbow(23))+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

maf.dist.plot.all <- ggplot(maf.dat, aes(x=MAF))+
  stat_density(color="black",  position= "identity", geom= "line")+
  ggtitle("MAF distribution")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

### printing out the maf.hw distributions
maf.title <- paste0("MAF distribution After PCA","\n", 
                    opt$name, " ", date())
maf.final.tiff.file <- file.path(opt$out, "10_final.maf.distributions.tiff")
tiff(maf.final.tiff.file, width = 2500, height = 1500, units = "px", res = 300, compression = "lzw")
grid.arrange(maf.dist.plot.chr, maf.dist.plot.all, 
             ncol=2, top=maf.title)
dev.off()
}
############################END################################
cat("\n[INFO]\t Finished plotting QC report")

####Done
