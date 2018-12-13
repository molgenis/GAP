###################################
### QC report for genotyping data
### date: 30-11-2018
### version: 0.01
### authors: EL - RAG
###################################

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


##Arguments
#########################################################################################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input path to perform QC", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="Outout path to save report", metavar="character"), 
  
  make_option(c("-r", "--refMaf"), type="character", default="out.txt", 
              help="Path to the reference ", metavar="character")
); 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied as input path", call.=FALSE)
}


output <- as.character(opt$out)

#########################################################################################################
### Main
#########################################################################################################


############### Samples and SNPs passing call rate QC 

## loading samples and SNPs
samples_all.file <- file.path(opt$input,"full.ind")
samples_all <- read.table(samples_all.file)

snps_all.file <- file.path(opt$input,"full.snps")
snps_all <- read.table(snps_all.file)


samples_80.file <- file.path(opt$input, paste0("1_CR80/","incl80.samples"))
samples_80 <- read.table(samples_80.file)

snps_80.file <- file.path(opt$input, paste0("1_CR80/","incl80.vars"))
snps_80 <- read.table(snps_80.file)

samples_95.file <- file.path(opt$input, paste0("2_CR95/", "incl95.samples"))
samples_95<-read.table(samples_95.file)
colnames(samples_95)<-c("IID")
samples_95$IID<-as.character(samples_95$IID)

snps_95.file <- file.path(opt$input, paste0("2_CR95/","incl95.vars"))
snps_95 <- read.table(snps_95.file)

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
  theme_bw()+
  geom_text(aes(label=samples), vjust=0)+
  theme(text=element_text(size=10, family = 'Helvetica'))

##barplot for snps
barplot.snps <- ggplot(incl, aes(group,snps)) +
  geom_bar(stat = "identity", aes(colour=group, fill=group), 
           position = position_dodge(width = 0.5))+
  coord_cartesian(ylim = c(1.002*snps["All"]-3.5*z_snps,1.002*snps["All"]))+
  theme_bw()+
  geom_text(aes(label=snps), vjust=0)+
  theme(text=element_text(size=10, family = 'Helvetica'))


samples.snps.barplot.file <- file.path(output, "samples.snps.barplot.tiff")
tiff(samples.snps.barplot.file,  
     width = 2000, height = 1500, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(barplot.samples, barplot.snps, ncol=2)
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
            Call_rate=((N_GENOT-N_MISST)/N_GENOT) )
##Plot CR distribution
#by chromosome
cr.by.chr.dens <- ggplot(cr.dat, aes(x=1-F_MISS))+
  geom_density(aes(colour= CHR))+
  facet_wrap(~CHR, scales = 'free')+
  theme_bw()+
  theme(text = element_text(size=10, family='Helvetica'))
##by chromosome merged
cr.merged.chrs.dens <- ggplot(cr.dat, aes(x=1-F_MISS))+
  geom_density(aes(colour= CHR))+
  theme_bw()+
  theme(text = element_text(size=10, family='Helvetica'))
#by full autosomes+XY
cr.overall <- ggplot(cr.datG,aes(x=Call_rate))+
  geom_density()+
  theme_bw()+
  theme(text=element_text(size=10,family='Helvetica'))


cr.plots.file <- file.path(output, "callRate.plots.tiff")
tiff(cr.plots.file,  
     width = 2000, height = 3000, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(cr.merged.chrs.dens, cr.overall, ncol=1)
dev.off()


cr.plots.by.chr.file <- file.path(output, "callRate.by.chr.plots.tiff")
tiff(cr.plots.by.chr.file,  
     width = 2800, height = 3200, 
     units = "px", res = 300, compression = "lzw")
plot(cr.by.chr.dens)
dev.off()

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

## 
hw.dat$MAF <- maf.dat$MAF[match(hw.dat$SNP ,maf.dat$SNP)] 


maf.dist.plot.chr <- ggplot(maf.dat, aes(x=MAF))+
  stat_density(geom="line", aes(color= CHR))+
  ggtitle("MAF distribution per chromosome")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

maf.dist.plot.all <- ggplot(maf.dat, aes(x=MAF))+
  stat_density(geom="line", color="black")+
  ggtitle("MAF distribution")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))


hw.dist.plot.chr <- ggplot(hw.dat[which(hw.dat$P > 0.000005),], aes(x=-log10(P)))+
  stat_density(geom="line", aes(color= CHR))+
  ggtitle("HW pVal distribution  \n per chromosome")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

hw.maf.dist.plot.chr <- ggplot(hw.dat[which(hw.dat$P > 0.000005 & hw.dat$MAF > 0.1),], aes(x=-log10(P)))+
  stat_density(geom="line", aes(color= CHR))+
  ggtitle("HW pVal distribution of \n SNPs with a MAF > 0.1")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

hw.maf.dist.plot.all <- ggplot(hw.dat[which(hw.dat$P > 0.000005 & hw.dat$MAF > 0.1),], aes(x=-log10(P)))+
  stat_density(geom="line")+
  ggtitle("HW pVal distribution of \n SNPs with a MAF > 0.1")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))


# layout matrix for grid.arrange
lay <- rbind(c(1,1,1,2,2,2),
             c(3,3,4,4,5,5)
)

### printing out the maf.hw distributions

maf.hw.tiff.file <- file.path(output, "maf.hw.distributions.tiff")
tiff(maf.hw.tiff.file, width = 3000, height = 2000, units = "px", res = 300, compression = "lzw")
grid.arrange(maf.dist.plot.chr, maf.dist.plot.all, 
             hw.dist.plot.chr, hw.maf.dist.plot.chr, hw.maf.dist.plot.all, 
             layout_matrix = lay)
dev.off()



############################################################
############### Correlation of MAF with 1000G EUR populations

# Read bim files from cohort 
#cohort.bim.path <- "/groups/umcg-wijmenga/tmp04/umcg-raguirre/pln_ugli/CR95"
cohort.bim.path <- file.path(opt$input, "2_CR95")
cohort.bim.files <- list.files(cohort.bim.path, pattern = ".bim$", full.names = TRUE)
cohort.bim.list <- lapply(cohort.bim.files, fread, data.table=FALSE)
cohort.bim.dat <- bind_rows(cohort.bim.list)
rm(cohort.bim.list)
# generate chr:pos IDs
cohort.bim.dat$id <- paste0(cohort.bim.dat[,1],":",cohort.bim.dat[,4])

# Read 1000g info 
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

cohort.bim.dat <- cohort.bim.dat[which(cohort.bim.dat$id %in% ref.maf.dat$id),]

cohort.bim.dat$ref.ref.allele <- ref.maf.dat$REF[match(cohort.bim.dat$id, ref.maf.dat$id)]
cohort.bim.dat$cohor.maf <- maf.dat$MAF[match(cohort.bim.dat[,2],maf.dat$SNP)]
cohort.bim.dat$ref.maf <- ref.maf.dat$EUR_AF[match(cohort.bim.dat$id, ref.maf.dat$id)]
## check that the major and minor allele are the same. 
# check allelic concordance 
swap.index <- which(cohort.bim.dat$ref.ref.allele != cohort.bim.dat[,6])
cohort.bim.dat$ref.maf.corrected <- cohort.bim.dat$ref.maf
cohort.bim.dat$ref.maf.corrected[swap.index] <- (1-as.numeric(cohort.bim.dat$ref.maf[swap.index]))
cohort.bim.dat$ref.maf.corrected <- as.numeric(cohort.bim.dat$ref.maf.corrected)
cohort.bim.dat$cohor.maf <- as.numeric(cohort.bim.dat$cohor.maf)

#plots
maf.ref.cohort.scatter <- ggplot(cohort.bim.dat, aes(x=ref.maf.corrected , y = cohor.maf))+
  geom_jitter(alpha= 0.6)+
  ggtitle("MAF from EUR 1000g vs cohort maf")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

# calculate density for scatter plot. 
cohort.bim.dat$density <- get_density(cohort.bim.dat$ref.maf.corrected, cohort.bim.dat$cohor.maf)
maf.ref.cohort.scatter.dens <- ggplot(cohort.bim.dat, aes(x=ref.maf.corrected , y = cohor.maf))+
  geom_point(alpha= 0.6, aes(color=density))+
  scale_color_viridis()+
  ggtitle("MAF from EUR 1000g vs cohort maf")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))


maf.ref.cohort.file <- file.path(output, "maf.ref.cohort.scatter.tiff")
tiff(maf.ref.cohort.file,  width = 1400, height = 2000, units = "px", res = 300, compression = "lzw")
grid.arrange(maf.ref.cohort.scatter, maf.ref.cohort.scatter.dens, ncol=1)
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



cat("\n[INFO]\t Finished plotting QC report")

