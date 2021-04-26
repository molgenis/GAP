###################################
### Genotype concordance GSA vs GWAS chip 
### date: 11-07-2019
### version: 0.01
### authors: EL - RAG
###################################
### New
# Added extra plots for MAF 
# Optimized concordance check 
# Defined filtered plink files for GWAS and GoNL 
# ### date: 11-07-2019
# Removed poor quality variants from the GoNL WGS file, which improves the concordance results. 
###################################

library(data.table)
library(tidyverse)
library(optparse)
library(gridExtra)
library(snpStats)

##cluster test
#opt<-list()
#opt$out<-"/groups/umcg-aad/tmp04/umcg-elopera/testdata/famcheck/"

#########################################################################################################
option_list = list(
  make_option(c("-p", "--plink"), type="character", default="/groups/umcg-aad/tmp04/umcg-elopera/ugli_blood_gsa/pairing.dat", 
              help="Path to plink file, it assumes .bim .fam and .bed", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="./plots/", 
              help="Output path to save report", metavar="character"),
  make_option(c("-l", "--llref"), type="character", default="/groups/umcg-aad/tmp04/umcg-elopera/concordanceCheck/gwas_og_plink_chr/filtered_merged_LL_GWAS", 
              help="reference file for Lifelines gwas data", metavar="character"),
  make_option(c("-g", "--gonlref"), type="character", default="/groups/umcg-aad/tmp04/umcg-elopera/concordanceCheck/gonl_og_plink_chr/filtered_merged_goNL_OG", 
              help="reference file for GoNL data", metavar="character"),
  make_option(c("-d", "--datapairing"), type="character", default="/groups/umcg-aad/tmp04/umcg-elopera/concordanceCheck/ugli.final.pairing.concordanceCheck.txt", 
              help="reference for pairing gonl and lifelines with UGLI data", metavar="character")
); 
opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


output <- opt$out
# output <- "/groups/umcg-aad/tmp04/umcg-elopera/concordanceCheck/Test1"
dir.create(path = output, recursive = T, showWarnings = T )

#########################################################################################################
### Reference and pre-filtered plink files of reference data (LifeLines GWAS array and GoNL WGS)
### Pairing file which links back the IDs from UGLI to LifeLines GWAS array and GoNL WGS
#########################################################################################################

gwas.plink.file <- opt$llref
gonl.plink.file <- opt$gonlref
pairing.file <- opt$datapairing
#########################################################################################################
### Main
#########################################################################################################


## Process merged input plink to include only samples overlapping with GoNL and GWAS LL 
ugli.file <- as.character(opt$plink)
#ugli.file <- "/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC_err/5_Relatedness/proc/full_data"
ugli.fam <- fread(paste0(ugli.file, ".fam"), data.table = FALSE)
pairing <- fread(pairing.file, data.table = FALSE)

pairing$ugli.fam <- NA 
pairing$ugli.fam[which(pairing$ugli.id.clean %in% ugli.fam$V1)] <- pairing$ugli.id.clean[which(pairing$ugli.id.clean %in% ugli.fam$V1)]
pairing$ugli.fam[which(pairing$SAMPLE_IDENTIFIER_Rotterdam %in% ugli.fam$V1)] <- pairing$ugli.id.clean[which(pairing$SAMPLE_IDENTIFIER_Rotterdam %in% ugli.fam$V1)]
pairing$ugli.fam[which(is.na(pairing$ugli.overlaps))] <- NA

pairing.ugli.fam.keep.file <- file.path(output, "pairing.ugli.fam.keep.txt")
write.table(data.frame(t1= pairing$ugli.fam,
                       t2=pairing$ugli.fam),
                      file = pairing.ugli.fam.keep.file, 
                      sep="\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

ugli.concordance.file <- paste0(output,"/ugli.concordanceCheck")


plink.system.call <- paste0("ml plink \\\n
                            plink \\
                            --bfile ",ugli.file," \\
                            --allow-no-sex \\
                            --make-bed \\
                            --keep ", pairing.ugli.fam.keep.file," \\
                            --out ", ugli.concordance.file, "\n")

system(plink.system.call)

###########################################################################
###########################################################################
##### Concordance check in GoNL WGS

gonl.plink <- snpStats::read.plink(bed = gonl.plink.file)
ugli.plink <- snpStats::read.plink(bed = ugli.concordance.file)


## get correct index of samples. First from ugli. then subset again for which ones are actually in the fam file. 
ugli.sample.index <- pairing$ugli.fam[which(!is.na(pairing$gonl.id) & !is.na(pairing$ugli.fam))]

## double check actual samples loaded in the plink file ugli.plink$genotypes@.Data
ugli.sample.index <- ugli.sample.index[which(ugli.sample.index %in% rownames(ugli.plink$genotypes@.Data))]
gonl.sample.index <- pairing$gonl.id[match(ugli.sample.index, pairing$ugli.fam)]



snp.index <- intersect(gonl.plink$map$snp.name, ugli.plink$map$snp.name)

gonl.geno <- gonl.plink$genotypes@.Data[gonl.sample.index,snp.index]
gonl.stats <- snpStats::col.summary(new("SnpMatrix", gonl.geno))
gonl.geno <- as.matrix(gonl.geno)
gonl.geno <- apply(gonl.geno, 2, as.numeric)
rownames(gonl.geno) <- ugli.sample.index

#ugli.sample.index.short <- ugli.sample.index[ugli.sample.index %in% rownames(ugli.plink$genotypes@.Data)]
ugli.geno <- ugli.plink$genotypes@.Data[ugli.sample.index,snp.index]
ugli.stats <- snpStats::col.summary(new("SnpMatrix", ugli.geno))
ugli.geno <- as.matrix(ugli.geno)
ugli.geno <- apply(ugli.geno, 2, as.numeric)
rownames(ugli.geno) <- ugli.sample.index



ugli.snp.index <- match(snp.index, ugli.plink$map$snp.name)
gonl.snp.index <- match(snp.index, gonl.plink$map$snp.name)


allelles.ok.snps <- (ugli.plink$map$allele.1[ugli.snp.index] == gonl.plink$map$allele.1[gonl.snp.index] &
                     ugli.plink$map$allele.2[ugli.snp.index] == gonl.plink$map$allele.2[gonl.snp.index])

swap.index <-  ugli.plink$map$allele.1[ugli.snp.index][!allelles.ok.snps] == gonl.plink$map$allele.2[gonl.snp.index][!allelles.ok.snps] & 
               ugli.plink$map$allele.2[ugli.snp.index][!allelles.ok.snps] == gonl.plink$map$allele.1[gonl.snp.index][!allelles.ok.snps]

#### remove all the swaps that could not be confirmed. 

conflicted.snps <- ugli.plink$map$snp.name[ugli.snp.index][!allelles.ok.snps][which(!swap.index)]

swap.index[which(swap.index==FALSE)] <- TRUE

if(length(swap.index) >1 & all(swap.index)) {
  swap.geno.mat <- gonl.geno[,which(!allelles.ok.snps)]
  swap.geno.mat <- ifelse(swap.geno.mat == 01, 03, ifelse(gonl.geno[,which(!allelles.ok.snps)] == 03, 01, 02))
  gonl.geno[,which(!allelles.ok.snps)] <- swap.geno.mat
} 

ugli.geno <- ugli.geno[,(colnames(ugli.geno) %in% conflicted.snps ==FALSE)]
gonl.geno <- gonl.geno[,(colnames(gonl.geno) %in% conflicted.snps ==FALSE)]


concordance.mat.gonl <- ugli.geno == gonl.geno
  
n.snps <- length(intersect(colnames(gonl.geno), colnames(ugli.geno)))
sample.concordance.gonl <- data.frame(sample.gonl.id=gonl.sample.index,
                                 sample.ugli.id=ugli.sample.index,
                                 concordance= apply(concordance.mat.gonl, 1, function(x){(sum(x)/n.snps) * 100}))


n.samples <- length(gonl.sample.index)
snp.concordance.gonl <- data.frame(snp=colnames(ugli.geno), 
                              concordance= apply(concordance.mat.gonl, 2, function(x){(sum(x)/n.samples) * 100})
                              )

snp.concordance.gonl$maf.ugli <- ugli.stats$MAF[match(snp.concordance.gonl$snp, rownames(ugli.stats))]
snp.concordance.gonl$maf.gonl <- gonl.stats$MAF[match(snp.concordance.gonl$snp, rownames(gonl.stats))]

## adding also the call rate from gonl as variants from the raw file haven't been included. 
snp.concordance.gonl$call.rate.gonl <- gonl.stats$Call.rate[match(snp.concordance.gonl$snp, rownames(gonl.stats))]
snp.concordance.gonl$HW.z.gonl <- gonl.stats$z.HWE[match(snp.concordance.gonl$snp, rownames(gonl.stats))]

snp.concordance.gonl$delta.maf <- snp.concordance.gonl$maf.ugli-snp.concordance.gonl$maf.gonl

## remove all variants for which gonl reports a call rate lower than 0.99
snp.concordance.gonl <- snp.concordance.gonl[which(snp.concordance.gonl$call.rate.gonl > 0.99),]
snp.concordance.gonl <- snp.concordance.gonl[which(!is.na(snp.concordance.gonl$HW.z.gonl)),]
snp.concordance.gonl <- snp.concordance.gonl[which(snp.concordance.gonl$maf.gonl != 0),]
#rm(gonl.plink, n.samples, n.snps, concordance.mat.gonl, ugli.geno, gonl.geno)
###########################################################################
###########################################################################
##### Concordance check in LifeLines GWAS array

gwas.plink <- snpStats::read.plink(bed = gwas.plink.file)
ugli.plink <- snpStats::read.plink(bed = ugli.concordance.file)


pairing <- pairing[which(pairing$gwas.ll.id %in% rownames(gwas.plink$genotypes@.Data) &
                         pairing$ugli.fam %in% rownames(ugli.plink$genotypes@.Data)),]

ugli.sample.index <- pairing$ugli.fam[which(!is.na(pairing$ugli.fam) & !is.na(pairing$gwas.ll.id))]
gwas.sample.index <- pairing$gwas.ll.id[which(!is.na(pairing$ugli.fam) &  !is.na(pairing$gwas.ll.id))]

#ugli.sample.index <- intersect(ugli.sample.index, rownames(ugli.plink$genotypes@.Data))
#gwas.sample.index <- intersect(gwas.sample.index, rownames(gwas.plink$genotypes@.Data))

snp.index <- intersect(gwas.plink$map$snp.name, ugli.plink$map$snp.name)


gwas.geno <- gwas.plink$genotypes@.Data[gwas.sample.index,snp.index]
gwas.stats <- snpStats::col.summary(new("SnpMatrix", gwas.geno))
gwas.geno <- as.matrix(gwas.geno)
gwas.geno <- apply(gwas.geno, 2, as.numeric)
rownames(gwas.geno) <- gwas.sample.index

ugli.geno <- ugli.plink$genotypes@.Data[ugli.sample.index,snp.index]
ugli.stats <- snpStats::col.summary(new("SnpMatrix", ugli.geno))
ugli.geno <- as.matrix(ugli.geno)
ugli.geno <- apply(ugli.geno, 2, as.numeric)
rownames(ugli.geno) <- gwas.sample.index

###########################################################################
## check minor and major allele are in the same order. 

ugli.snp.index <- match(snp.index, ugli.plink$map$snp.name)
gwas.snp.index <- match(snp.index, gwas.plink$map$snp.name)

allelles.ok.snps <- (ugli.plink$map$allele.1[ugli.snp.index] == gwas.plink$map$allele.1[gwas.snp.index] &
                       ugli.plink$map$allele.2[ugli.snp.index] == gwas.plink$map$allele.2[gwas.snp.index])

swap.index <-  ugli.plink$map$allele.1[ugli.snp.index][!allelles.ok.snps] == gwas.plink$map$allele.2[gwas.snp.index][!allelles.ok.snps] & 
  ugli.plink$map$allele.2[ugli.snp.index][!allelles.ok.snps] == gwas.plink$map$allele.1[gwas.snp.index][!allelles.ok.snps]


if(length(swap.index) >= 1 & all(swap.index)) {
  
  swap.geno.mat <- gwas.geno[,which(!allelles.ok.snps)]
  swap.geno.mat <- ifelse(swap.geno.mat == 01, 03, ifelse(gwas.geno[,which(!allelles.ok.snps)] == 03, 01, 02))
  gwas.geno[,which(!allelles.ok.snps)] <-swap.geno.mat
  
}

concordance.mat.gwas <- ugli.geno == gwas.geno

n.snps <- length(snp.index)
sample.concordance.gwas <- data.frame(sample.gwas.id=gwas.sample.index,
                                 sample.ugli.id=ugli.sample.index,
                                 concordance= apply(concordance.mat.gwas, 1, function(x){(sum(x)/n.snps) * 100})
                                 )

n.samples <- length(gwas.sample.index)
snp.concordance.gwas <- data.frame(snp=snp.index, 
                              concordance= apply(concordance.mat.gwas, 2, function(x){(sum(x)/n.samples) * 100}),
                              swap= allelles.ok.snps)


snp.concordance.gwas$maf.ugli <- ugli.stats$MAF[match(snp.concordance.gwas$snp, rownames(ugli.stats))]
snp.concordance.gwas$maf.gwas <- gwas.stats$MAF[match(snp.concordance.gwas$snp, rownames(gwas.stats))]
snp.concordance.gwas$delta.maf <- snp.concordance.gwas$maf.ugli-snp.concordance.gwas$maf.gwas

###########################################################################
###########################################################################
##### Plotting
output.plots <- paste0(output, "/concordancePlots/")
dir.create(output.plots, recursive = TRUE)

################################
######## SNP concordance plot

gwas.snp.concordance.plot <- ggplot(snp.concordance.gwas, aes(x=as.numeric(concordance), stat(count)))+
  geom_density(bw=1.5, geom="line", fill= "#5DBCD2", alpha= 0.6)+
  #xlim(c(95,100))+
  #ylim(c(0,50000))+
  xlim(c(50,100))+
  ggtitle(paste0("Concordance in GWAS array \n n=", nrow(snp.concordance.gwas)))+
  ylab("Number of SNPs")+
  xlab("Concordance (in %)")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 11))

gwas.snp.concordance.gonl <- ggplot(snp.concordance.gonl, aes(x=as.numeric(concordance), stat(count)))+
  geom_density(bw=1.5, geom="line", fill="#FF4F00", alpha= 0.6)+
  #xlim(c(95,100))+
  #ylim(c(0,50000))+
  xlim(c(50,100))+
  ggtitle(paste0("Concordance in GoNL WGS \n n=", nrow(snp.concordance.gonl)))+
  ylab("Number of SNPs")+
  xlab("Concordance (in %)")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 11))


snp.concordance.plot.file <- file.path(output.plots, paste0("snp.concordance.plot.tiff"))
tiff(snp.concordance.plot.file, width = 1800, height = 800, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(gwas.snp.concordance.plot, gwas.snp.concordance.gonl, nrow=1)
dev.off()

################################
######## Sample concordance plot

sample.concordance.plot.gwas <- ggplot(sample.concordance.gwas, aes(x=concordance,  stat(count)))+
  geom_density(bw=1, geom="line", fill= "#5DBCD2", alpha= 0.6)+
  xlim(c(75,100))+
  ggtitle(paste0("Concordance in GWAS array \n n=", nrow(sample.concordance.gwas)))+
  ylab("Number of samples")+
  xlab("Concordance (in %)")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 11))

sample.concordance.plot.gonl <- ggplot(sample.concordance.gonl, aes(x=concordance,  stat(count)))+
  geom_density(bw=1, geom="line", fill="#FF4F00", alpha= 0.6)+
  xlim(c(75,100))+
  ggtitle(paste0("Concordance in GoNL WGS \n n=", nrow(sample.concordance.gonl)))+
  ylab("Number of samples")+
  xlab("Concordance (in %)")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 11))

sample.concordance.plot.file <- file.path(output.plots, paste0("sample.concordance.plot.tiff"))
tiff(filename = sample.concordance.plot.file, width = 1800, height = 800, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(sample.concordance.plot.gwas, sample.concordance.plot.gonl, nrow=1)
dev.off()



################################
######## MAF plots 

maf.plot.gwas <- ggplot(snp.concordance.gwas, aes(x=maf.gwas, y=maf.ugli))+
  geom_point(alpha= 0.5, size=0.5, color= "#5DBCD2")+
  ylab("MAF in UGLI (GSEA array)")+
  xlab("MAF in LL (GWAS array)")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica",size = 10))

log.maf.plot.gwas <- ggplot(snp.concordance.gwas, aes(x=log(maf.gwas), y=log(maf.ugli)))+
  geom_point(alpha= 0.5, size=0.5, color= "#5DBCD2")+
  ylab("log(MAF) in UGLI (GSEA array)")+
  xlab("log(MAF) in LL (GWAS array)")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 10))

delta.maf.concordance.plot.gwas <- ggplot(snp.concordance.gwas, aes(x=concordance, y=abs(delta.maf)))+
  geom_point(alpha= 0.5, size=0.5, color= "#5DBCD2")+
  geom_rug(size=0.1, alpha=0.8, color="#5DBCD2")+
  ylab("\u0394 abs(MAF)")+
  xlab("Concordance (in %)")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 10))

delta.hex.maf.concordance.plot.gwas <- ggplot(snp.concordance.gwas, aes(x=concordance, y= abs(delta.maf)))+
  geom_hex(alpha= 0.8, bins= 20)+
  scale_fill_continuous(low = "gray", high = "#5DBCD2", trans = "log10")+
  guides(guides(fill = guide_colourbar(barwidth = 0.5, barheight = 10, ticks = FALSE)))+
  ylab("\u0394 abs(MAF)")+
  xlab("Concordance (in %)")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 10))


hist.concordance.plot.gwas <- ggplot(snp.concordance.gwas, aes(x=concordance, y= stat(count)))+
  geom_histogram(alpha= 0.7, fill= "#5DBCD2", binwidth = 10)+
  xlab("Concordance (in %)")+
  ylab("Number of SNPs")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 10))

boxplot.concordance.gwas <- ggplot(snp.concordance.gwas, aes(x=maf.ugli , y= concordance))+
  geom_boxplot(fill= "#5DBCD2", aes(group = cut_width(maf.ugli, 0.10)))+
  ylab("Concordance (in %)")+
  xlab("MAF in UGLI")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 10))




maf.plot.gonl <- ggplot(snp.concordance.gonl, aes(x=maf.gonl, y=maf.ugli))+
  geom_point(alpha= 0.5, size=0.5, color= "#FF4F00")+
  ylab("MAF in UGLI (GSEA array)")+
  xlab("MAF in GoNL")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica",size = 10))


log.maf.plot.gonl <- ggplot(snp.concordance.gonl, aes(x=log(maf.gonl), y=log(maf.ugli)))+
  geom_point(alpha= 0.5, size=0.5, color= "#FF4F00")+
  ylab("log(MAF) in UGLI (GSEA array)")+
  xlab("log(MAF) in GoNL")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 10))

delta.maf.concordance.plot.gonl <- ggplot(snp.concordance.gonl, aes(x=concordance, y= abs(delta.maf)))+
  geom_point(alpha= 0.5, size=0.5, color= "#FF4F00")+
  geom_rug(size=0.1, alpha=0.8, color="#FF4F00")+   
  ylab("\u0394 abs(MAF)")+
  xlab("Concordance (in %)")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 10))

delta.hex.maf.concordance.plot.gonl <- ggplot(snp.concordance.gonl, aes(x=concordance, y= abs(delta.maf)))+
  geom_hex(alpha= 0.8, bins= 20)+
  scale_fill_continuous(low = "gray", high = "#FF4F00", trans = "log10")+
  guides(guides(fill = guide_colourbar(barwidth = 0.5, barheight = 10, ticks = FALSE)))+
  ylab("\u0394 abs(MAF)")+
  xlab("Concordance (in %)")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 10))

hist.concordance.plot.gonl <- ggplot(snp.concordance.gonl, aes(x=concordance, y= stat(count)))+
  geom_histogram(alpha= 0.7, fill= "#FF4F00", binwidth = 10)+
  xlab("Concordance (in %)")+
  ylab("Number of SNPs")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 10))

boxplot.concordance.gonl <- ggplot(snp.concordance.gonl, aes(x=maf.ugli , y= concordance))+
  geom_boxplot(fill= "#FF4F00", aes(group = cut_width(maf.ugli, 0.10)))+
  ylab("Concordance (in %)")+
  xlab("MAF in UGLI")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 10))



maf.concordance.plot.file <- file.path(output.plots, paste0("maf.concordance.plot.tiff"))
tiff(filename = maf.concordance.plot.file, width = 4200, height = 1600, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(maf.plot.gwas, log.maf.plot.gwas, delta.maf.concordance.plot.gwas, delta.hex.maf.concordance.plot.gwas, hist.concordance.plot.gwas, boxplot.concordance.gwas,
             maf.plot.gonl, log.maf.plot.gonl, delta.maf.concordance.plot.gonl, delta.hex.maf.concordance.plot.gonl, hist.concordance.plot.gonl, boxplot.concordance.gonl,
             nrow=2)
dev.off()


snp.concordance.gwas$gonl.concordance <- snp.concordance.gonl[rownames(snp.concordance.gwas),"concordance"]
combined.plot <- ggplot(snp.concordance.gwas, aes(x=concordance , y= gonl.concordance))+
  geom_point(alpha= 0.5, size=0.5, aes(color= maf.ugli))+
  guides(guides(color = guide_colourbar(barwidth = 0.5, barheight = 10, ticks = FALSE)))+
  xlab("Concordance (in %) GWAS")+
  ylab("Concordance (in %) GoNL")+
  theme_classic()+
  theme(text = element_text(family = "Helvetica", size = 10))


combined.concordance.plot.file <- file.path(output.plots, paste0("combined.concordance.plot.tiff"))
tiff(filename = combined.concordance.plot.file, width = 900, height = 800, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(combined.plot, nrow=1)
dev.off()


### write out reports
gonl.sample.report.file <- file.path(output, paste0("gonl.sample.report.csv"))
write.csv(x =sample.concordance.gonl , gonl.sample.report.file, quote = FALSE)

gwas.sample.report.file <- file.path(output, paste0("gwas.sample.report.csv"))
write.csv(x =sample.concordance.gwas , gwas.sample.report.file, quote = FALSE)

gonl.snp.report.file <- file.path(output, paste0("gonl.snp.report.csv"))
write.csv(x =snp.concordance.gonl , gonl.snp.report.file, quote = FALSE)
 
gwas.snp.report.file <- file.path(output, paste0("gwas.snp.report.csv"))
write.csv(x =snp.concordance.gwas , gwas.snp.report.file, quote = FALSE)


