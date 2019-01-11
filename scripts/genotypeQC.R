###################################
### QC report for genotyping data
### date: 11-01-2019
### version: 0.01
### authors: EL - RAG
###################################
### New
###################################

## New arguents - project name, sample info sheet 
## New plot showing starting samples and plate numbers 

## Added relatedness analysis figures. 

## Updated plot MAF correlation to properly asses reference and alternative alleles
## included MAF correlation with goNL 

## changed geom_density -> stat_density using position = "identity" and geom = "line 
## typo on plot title.
######################################################################
## example run  

#RScript genotypeQC.R -i ***/plink_autosome_QC \
#  -o 
#  -r /groups/umcg-wijmenga/tmp04/umcg-raguirre/pln_ugli/eurMAF_1000g/ 


#########################################################################################################
#########################################################################################################
### Functions and libraries
#########################################################################################################
#########################################################################################################

library(tidyverse)
#library(ggsci)
library(data.table)
library(optparse)
#library(ggridges)
library(gridExtra)
library(MASS)
library(viridis)
# taken from https://slowkow.com/notes/ggplot2-color-by-density/
# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

## Check ref and alt alleles are the same. 
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


##Arguments
## Set up arguments for testing 
# opt <- list()
# opt$input <- "/groups/umcg-aad/tmp04/umcg-elopera/plink_autosome_QC"
# opt$out <- "/groups/umcg-wijmenga/tmp04/umcg-raguirre/pln_ugli/qcPlots_test"
# opt$name <- "test01"
# opt$sampleinfo <- "/groups/umcg-aad/tmp04/projects/IBD_part/run01/jobs/IBD_part.csv"
# opt$refMaf <- "/groups/umcg-wijmenga/tmp04/umcg-raguirre/pln_ugli/eurMAF_1000g"
# opt$goNLrefMaf <- "/groups/umcg-wijmenga/tmp04/umcg-raguirre/pln_ugli/gonl_MAF/goNL_MAF.edR.INFO"

#########################################################################################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input path to perform QC", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="genoQC_Report", 
              help="Output path to save report", metavar="character"),
  
  make_option(c("-n", "--name"), type="character", default="genoQC_Report", 
              help="Name of the project", metavar="character"), 
  
  make_option(c("-s", "--sampleinfo"), type="character", default=NULL, 
              help="path to sample info spreadsheet", metavar="character"), 
  
  make_option(c("-r", "--refMaf"), type="character", default=NULL, 
              help="Path to the reference MAF of SNPs", metavar="character"),

    make_option(c("-g", "--goNLrefMaf"), type="character", default=NULL, 
              help="Path to the reference MAF of SNPs", metavar="character")
); 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied as input path", call.=FALSE)
}
#
output <- as.character(opt$out)


#########################################################################################################
### Main
#########################################################################################################

############### First plot - 
# Date and project name 
# Table of all plates used in the analysis 
# Table showing number of starting samples 
# Cohort within freeze? 

# testing
# opt <- list(); opt$sampleinfo <- "/groups/umcg-aad/tmp04/projects/IBD_part/run01/jobs/IBD_part.csv"

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


############### Samples and SNPs passing call rate QC 

## loading samples and SNPs
samples_all.file <- file.path(opt$input,"full.ind")
samples_all <- fread(samples_all.file, data.table = FALSE)

snps_all.file <- file.path(opt$input,"full.snps")
snps_all <- fread(snps_all.file)


samples_80.file <- file.path(opt$input, paste0("1_CR80/","incl80.samples"))
samples_80 <- fread(samples_80.file, data.table = FALSE)

snps_80.file <- file.path(opt$input, paste0("1_CR80/","incl80.vars"))
snps_80 <- fread(snps_80.file)

samples_95.file <- file.path(opt$input, paste0("2_CR95/", "incl95.samples"))
samples_95<-fread(samples_95.file, data.table = FALSE)
colnames(samples_95)<-c("IID")
samples_95$IID<-as.character(samples_95$IID)

snps_95.file <- file.path(opt$input, paste0("2_CR95/","incl95.vars"))
snps_95 <- fread(snps_95.file, data.table = FALSE)

################################################
###group data to make the barplots
samples <- c("All"=nrow(samples_all),"CR>80"=nrow(samples_80), "CR>95"=nrow(samples_95))
snps <- c("All"=nrow(snps_all),"CR>80"=nrow(snps_80), "CR>95"=nrow(snps_95) )

incl <- data.frame(samples, snps,"group"= names(snps))

#variables to calculate the ylim of the barplots
z_samples <- mean(samples["All"]-samples["CR>80"], samples["CR>80"]-samples["CR>95"])
z_snps <- mean(snps["All"]-snps["CR>80"], snps["CR>80"]-snps["CR>95"])

##barplot for samples
barplot.samples <- ggplot(incl, aes(group,samples))+
  geom_bar(stat = "identity", aes(colour=group, fill=group), position = position_dodge(width = 0.5))+
  coord_cartesian(ylim = c(1.002*samples["All"]-3.5*z_samples,1.002*samples["All"] ))+ 
  ylab("Number of Samples")+
  theme_bw()+
  geom_text(aes(label=samples), vjust=0)+
  theme(text=element_text(size=10, family = 'Helvetica'), legend.position = "none")

##barplot for snps
barplot.snps <- ggplot(incl, aes(group,snps)) +
  geom_bar(stat = "identity", aes(colour=group, fill=group), 
           position = position_dodge(width = 0.5))+
  coord_cartesian(ylim = c(1.002*snps["All"]-3.5*z_snps,1.002*snps["All"]))+
  ylab("Number of SNPs")+
  theme_bw()+
  geom_text(aes(label=snps), vjust=0)+
  theme(text=element_text(size=10, family = 'Helvetica'), legend.position = "none")


samples.snps.barplot.file <- file.path(output, "02_samples.snps.barplot.tiff")
samples.snps.barplot.title <- paste0("Number of samples and genotypes at multiple call rate thresholds", "\n",
                                     "(from autosomal and pseudo-autosomal chromosomes)", "\n", 
                                     opt$name, " ", date())

tiff(samples.snps.barplot.file,  
     width = 1500, height = 1000, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(barplot.samples, barplot.snps, 
             ncol=2, top=samples.snps.barplot.title)
dev.off()


############### FInal CR distribution
###look for filepath of CR80####
#create list of filenames
cr.files.path <- file.path(opt$input, "1_CR80/")
cr.files <- list.files(cr.files.path, pattern = ".imiss", full.names = TRUE)

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
#extract the samples below 95 CR
cr.dat <-semi_join(cr.dat,samples_95,by="IID")
##get the genomic CR data base
cr.datG<-cr.dat %>% group_by(IID) %>%
  summarise(N_MISST=sum(N_MISS),
            N_GENOT=sum(N_GENO),
            Call_rate=((N_GENOT-N_MISST)/N_GENOT))

##Plot CR distribution
# define levels of CHR to plot CHRs in order. 
cr.dat$CHR <- gsub(cr.dat$CHR, pattern = "chr_", replacement ="" )
cr.dat$CHR <- factor(cr.dat$CHR, levels=c(1:22, "XY"))

#by chromosome 
cr.by.chr.dens <- ggplot(cr.dat, aes(x=1-F_MISS))+
  stat_density(aes(colour= CHR),geom= "line", position= "identity")+
  facet_wrap(~CHR, scales = 'free')+
  xlab("Call rate")+
  theme_bw()+
  theme(text = element_text(size=10, family='Helvetica'), legend.position = "none")
##by chromosome merged
cr.merged.chrs.dens <- ggplot(cr.dat, aes(x=1-F_MISS))+
  stat_density(aes(colour= CHR), geom= "line", position= "identity")+
  xlab("Call rate")+
  theme_bw()+
  theme(text = element_text(size=10, family='Helvetica'))
#by full autosomes+XY
cr.overall <- ggplot(cr.datG,aes(x=Call_rate))+
  stat_density(geom= "line", position= "identity")+
  xlab("Call rate")+
  theme_bw()+
  theme(text=element_text(size=10,family='Helvetica'))


cr.plots.file <- file.path(output, "03.1_callRate.plots.tiff")
cr.plot.tilte <- paste0("Distribution of call rate from autosomal and pseudo-autosomal chromosomes", "\n", opt$name, " ", date())
tiff(cr.plots.file,  
     width = 3000, height = 1300, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(cr.merged.chrs.dens, cr.overall, nrow=1, top= cr.plot.tilte)
dev.off()


cr.plots.by.chr.file <- file.path(output, "03.2_callRate.by.chr.plots.tiff")
cr.plot.grid.tilte <- paste0("Distribution of call rate from autosomal and pseudo-autosomal chromosomes", 
                             "\n", opt$name, " ", date())
tiff(cr.plots.by.chr.file,  
     width = 2800, height = 3200, 
     units = "px", res = 300, compression = "lzw")
plot(cr.by.chr.dens+ggtitle(cr.plot.grid.tilte))
dev.off()

###########################################################################
############### MAF and HW distributions across all chromosomes 
### output95
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
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

maf.dist.plot.all <- ggplot(maf.dat, aes(x=MAF))+
  stat_density(color="black",  position= "identity", geom= "line")+
  ggtitle("MAF distribution")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))


hw.dist.plot.chr <- ggplot(hw.dat, aes(x=-log10(P)))+
  stat_density(aes(color= CHR), position= "identity", geom= "line")+
  ggtitle("HW pVal distribution  \n per chromosome")+
  xlab("-log10(HW-P)")+
  xlim(c(0, 20))+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

hw.maf.dist.plot.chr <- ggplot(hw.dat[which(hw.dat$MAF > 0.01),], aes(x=-log10(P)))+
  stat_density(aes(color= CHR), position= "identity", geom= "line")+
  geom_vline(xintercept = 4)+
  ggtitle("HW pVal distribution of \n SNPs with a MAF > 0.01")+
  xlim(c(0, 20))+
  xlab("-log10(HW-P)")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

hw.maf.dist.plot.all <- ggplot(hw.dat[which(hw.dat$MAF > 0.01),], aes(x=-log10(P)))+
  stat_density(position= "identity", geom= "line")+
  geom_vline(xintercept = 4)+
  xlim(c(0, 20))+
  ggtitle("HW pVal distribution of \n SNPs with a MAF > 0.01")+
  xlab("-log10(HW-P)")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))


# layout matrix for grid.arrange
lay <- rbind(c(1,1,1,2,2,2),
             c(3,3,4,4,5,5)
)

### printing out the maf.hw distributions
maf.hw.title <- paste0("MAF and HW distributions (from autosomal and pseudo-autosomal chromosomes)","\n", opt$name, " ", date())
maf.hw.tiff.file <- file.path(output, "04_maf.hw.distributions.tiff")
tiff(maf.hw.tiff.file, width = 3000, height = 2000, units = "px", res = 300, compression = "lzw")
grid.arrange(maf.dist.plot.chr, maf.dist.plot.all, 
             hw.dist.plot.chr, hw.maf.dist.plot.chr, 
             hw.maf.dist.plot.all, 
             layout_matrix = lay, top=maf.hw.title)
dev.off()



############################################################
############### Correlation of MAF with 1000G EUR populations
############### Needs to be corrected for alleles 

# Read bim files from cohort 
# opt <-list(); opt$input <- "/groups/umcg-aad/tmp04/umcg-elopera/plink_autosome_QC"
cohort.bim.path <- file.path(opt$input, "2_CR95")
cohort.bim.files <- list.files(cohort.bim.path, pattern = ".bim$", full.names = TRUE)
cohort.bim.list <- lapply(cohort.bim.files, fread, data.table=FALSE)
cohort.bim.dat <- bind_rows(cohort.bim.list)
rm(cohort.bim.list)
colnames(cohort.bim.dat) <- c("chr", "snp", "morgan", "bp", "A1", "A2")
# generate chr:pos IDs
cohort.bim.dat$id <- paste0(cohort.bim.dat[,1],":",cohort.bim.dat[,4])

################################
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

############################## 
#### add the HW pvalue - remove all SNPS with a pValue <= 0.000005
cohort.bim.dat$cohort.hw.p <- hw.dat$P[match(as.character(cohort.bim.dat$snp), as.character(hw.dat$SNP))]

hw.index <- which(cohort.bim.dat$cohort.hw.p <= 0.000005)
if (length(hw.index) >= 1){
  cohort.bim.dat <- cohort.bim.dat[-hw.index,]
}


################################ 1000G
# Read 1000G info 
#ref.maf.path <- "/groups/umcg-wijmenga/tmp04/umcg-raguirre/pln_ugli/eurMAF_1000g/"
ref.maf.path <- opt$refMaf
ref.maf.files <- list.files(ref.maf.path, pattern = "_EURMAF_1000g.INFO$", full.names = TRUE)
ref.maf.list <- lapply(ref.maf.files, fread, data.table=FALSE)
ref.maf.dat <- bind_rows(ref.maf.list)
rm(ref.maf.list)

###eliminate the SNPs with more than two alleles????
ref.maf.dat<-ref.maf.dat[!grepl("\\,",ref.maf.dat[,5]),]
ref.maf.dat[,5]<-as.numeric(ref.maf.dat[,5])
##

# generate chr:pos IDs
ref.maf.dat$id <- paste0(ref.maf.dat$CHROM, ":",ref.maf.dat$POS)

### I would no longer remove markers that are not present in 1000G
#cohort.bim.dat <- cohort.bim.dat[which(cohort.bim.dat$id %in% ref.maf.dat$id),]

cohort.bim.dat$ref.alt.allele <- ref.maf.dat$ALT[match(cohort.bim.dat$id, ref.maf.dat$id)]
cohort.bim.dat$ref.ref.allele <- ref.maf.dat$REF[match(cohort.bim.dat$id, ref.maf.dat$id)]
cohort.bim.dat$ref.maf <- ref.maf.dat$EUR_AF[match(cohort.bim.dat$id, ref.maf.dat$id)]

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
# Read goNL info 
# opt$goNLrefMaf <- "/groups/umcg-wijmenga/tmp04/umcg-raguirre/pln_ugli/gonl_MAF/goNL_MAF.edR.INFO"

gonl.maf <- fread(opt$goNLrefMaf, data.table = FALSE)
gonl.maf$id <- paste0(gonl.maf$CHROM, ":", gonl.maf$POS)

cohort.bim.dat$goNl.maf <- gonl.maf$MAF[match(cohort.bim.dat$id, gonl.maf$id)]
cohort.bim.dat$goNl.alt <- gonl.maf$ALT[match(cohort.bim.dat$id, gonl.maf$id)]
cohort.bim.dat$goNl.ref <- gonl.maf$REF[match(cohort.bim.dat$id, gonl.maf$id)]

goNl.allele.matrix <- cbind(query.alt=cohort.bim.dat$A1, query.ref=cohort.bim.dat$A2, 
                       ref.alt=cohort.bim.dat$goNl.alt, ref.ref=cohort.bim.dat$goNl.ref)

cohort.bim.dat$goNl.allele.check <- apply(goNl.allele.matrix, 1, 
                                          FUN = function(x){check.ref.alt.alleles(
                                              query.alt= x[1], query.ref= x[2],  ref.alt= x[3], ref.ref= x[4])
                                            })

cohort.bim.dat$goNl.maf.corrected <- cohort.bim.dat$goNl.maf
cohort.bim.dat$goNl.maf.corrected[which(cohort.bim.dat$goNl.allele.check == "flip")] <- (1-as.numeric(cohort.bim.dat$goNl.maf[which(cohort.bim.dat$goNl.allele.check == "flip")]))
cohort.bim.dat$goNl.maf.corrected <- as.numeric(cohort.bim.dat$goNl.maf.corrected)


#plots
maf.ref.cohort.scatter <- ggplot(cohort.bim.dat[which(cohort.bim.dat$allele.check== "same"| cohort.bim.dat$allele.check== "flip"),], 
                                 aes(x=ref.maf.corrected , y = cohort.maf))+
  geom_point(size=0.5,alpha= 0.6, color="blue")+
  ggtitle("MAF from EUR 1000g vs cohort MAF")+
  ylab("MAF at cohort")+
  xlab("MAF at EUR 1000G")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

# calculate density for scatter plot. 
cohort.bim.dat.filter <- cohort.bim.dat[which(cohort.bim.dat$allele.check== "same" | cohort.bim.dat$allele.check== "flip"),]
cohort.bim.dat.filter <- cohort.bim.dat.filter[!is.na(cohort.bim.dat.filter$ref.maf.corrected),]
cohort.bim.dat.filter$density <- get_density(cohort.bim.dat.filter$ref.maf.corrected, cohort.bim.dat.filter$cohort.maf)
maf.ref.cohort.scatter.dens <- ggplot(cohort.bim.dat.filter, aes(x=ref.maf.corrected , y = cohort.maf))+
  geom_point(size=0.5, alpha= 0.6, aes(color=density))+
  scale_color_viridis()+
  ggtitle("MAF from EUR 1000g vs cohort MAF")+
  ylab("MAF at cohort")+
  xlab("MAF at EUR 1000G")+
  guides(colour=guide_colourbar(barheight = 0.5, barwidth = 4, title = "Density"))+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"), legend.position = "bottom")

### goNL plots
gonl.maf.ref.cohort.scatter <- ggplot(cohort.bim.dat[which(cohort.bim.dat$goNl.allele.check== "same"| cohort.bim.dat$goNl.allele.check== "flip"),], 
                                 aes(x=goNl.maf.corrected , y = cohort.maf))+
  geom_point(size=0.5, alpha= 0.6, color="orange")+
  ggtitle("MAF from GoNL vs cohort MAF")+
  ylab("MAF at cohort")+
  xlab("MAF at GoNL")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

gonl.cohort.bim.dat.filter <- cohort.bim.dat[which(cohort.bim.dat$goNl.allele.check == "same" | cohort.bim.dat$goNl.allele.check== "flip"),]
gonl.cohort.bim.dat.filter <- gonl.cohort.bim.dat.filter[!is.na(gonl.cohort.bim.dat.filter$goNl.maf.corrected),]
gonl.cohort.bim.dat.filter$density <- get_density(gonl.cohort.bim.dat.filter$goNl.maf.corrected, gonl.cohort.bim.dat.filter$cohort.maf)
gonl.maf.ref.cohort.scatter.dens <- ggplot(gonl.cohort.bim.dat.filter, aes(x=goNl.maf.corrected , y = cohort.maf))+
  geom_point(size=0.5, alpha= 0.6, aes(color=density))+
  scale_color_viridis()+
  ggtitle("MAF from GoNL vs cohort MAF")+
  ylab("MAF at cohort")+
  xlab("MAF at GoNL")+
  guides(colour=guide_colourbar(barheight = 0.5, barwidth = 4, title = "Density"))+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"), legend.position = "bottom")


### saving out plots. 
combine.plot.title <- paste0("MAF correlation with reference genomes (1000G and GoNL)", "\n", opt$name, " ", date())
maf.ref.cohort.file <- file.path(output, "05_maf.ref.cohort.scatter.tiff")

#plotting
tiff(maf.ref.cohort.file,  width = 2000, height = 2000, units = "px", res = 300, compression = "lzw")
grid.arrange(maf.ref.cohort.scatter, gonl.maf.ref.cohort.scatter, 
             maf.ref.cohort.scatter.dens, gonl.maf.ref.cohort.scatter.dens,
             top=combine.plot.title,
             ncol=2)
dev.off()
####Done



############################################################
############### Plot of relatedness

related.file <- file.path(opt$input, paste0("5_Relatedness/","autosomal_rel.genome"))
reldata1<-read.table(related.file, header = T)


##options for  plots Z1 vs Z0
z1vz0<- ggplot(reldata1,aes(x=Z0,y=Z1))+
  stat_sum(aes(size = factor(..n..)), geom = "point")+
  scale_size_discrete(range = c(2, 8))+
  theme_bw()+
  ggtitle("Sample relatedness")+
  theme(text=element_text(size=10,family="Helvetica"))+
  coord_cartesian(xlim = c(0, 1), ylim=c(0,0.6))


##pi_hat vs expected plot
obsvexp<- ggplot(reldata1,aes(x=PI_HAT,y=EZ))+
  geom_point()+
  theme_bw()+
  ggtitle("Relatednes coefficient")+
  theme(text=element_text(size=10,family="Helvetica"))+
  coord_cartesian(xlim = c(0, 1), ylim=c(0,1))+
  xlab("Observed R.C.")+
  ylab("Expected R.C.")

relatedness.file <- file.path(output, "relatedness.tiff")
tiff(relatedness.file,  
     width = 3000, height = 1500, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(z1vz0, obsvexp, ncol=2)
dev.off()

############################################################
cat("\n[INFO]\t Finished plotting QC report")
