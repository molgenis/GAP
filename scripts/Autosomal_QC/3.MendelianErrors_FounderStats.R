###################################
### Medelian erros and Stats on founders.
### date: 12-08-2019
### version: 0.01
### authors: EL - RAG
###################################

library(data.table)
library(tidyverse)
library(optparse)
library(gridExtra)
library(grid)
library(scales)

##cluster test
#opt<-list()
#opt$out<-"/groups/umcg-aad/tmp04/umcg-elopera/testdata/famcheck/"

## functions

HW.MAF.dist.plot  <- function(maf.dat, hw.dat, out, name=""){
  
  maf.dist.plot.chr <- ggplot(maf.dat, aes(x=MAF))+
    stat_density(aes(color= CHR), position= "identity", geom= "line")+
    ggtitle("MAF distribution per chromosome")+
    scale_color_manual(labels=(chr.labels), values=rainbow(23))+
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
    ggtitle("HW pVal distribution of \n SNPs with a MAF > 0.01")+
    scale_color_manual(labels=names(chr.labels), values=rainbow(23))+
    xlim(c(0, 20))+
    xlab("-log10(HW-P)")+
    theme_bw()+
    theme(text=element_text(size=10, family="Helvetica"))
  
  hw.maf.dist.plot.all <- ggplot(hw.dat[which(hw.dat$MAF > 0.01),], aes(x=-log10(P)))+
    stat_density(position= "identity", geom= "line")+
    geom_vline(xintercept = 6)+
    xlim(c(0, 20))+
    ggtitle("HW pVal distribution of \n SNPs with a MAF > 0.01")+
    xlab("-log10(HW-P)")+
    theme_bw()+
    theme(text=element_text(size=10, family="Helvetica"))
  
  
  # layout matrix for grid.arrange
  lay <- rbind(c(1,1,1,1,2,2),
               c(3,3,4,4,5,5)
  )
  
  
  ### printing out the maf.hw distributions
  maf.hw.title <- paste0("MAF and HW distributions from Founders","\n", name," ", date())
  maf.hw.conclude<-paste0(sum(hw.dat$P<0.0001)," markers excluded ", (nrow(hw.dat)/3), " remaining.")
  maf.hw.tiff.file <- file.path(out, paste0(name, "MAFHW.density.tiff"))
  tiff(maf.hw.tiff.file, width = 3000, height = 2100, units = "px", res = 300, compression = "lzw")
  grid.arrange(maf.dist.plot.chr, maf.dist.plot.all, 
               hw.dist.plot.chr, hw.maf.dist.plot.chr, 
               hw.maf.dist.plot.all, 
               layout_matrix = lay, top=maf.hw.title, 
               bottom=textGrob(maf.hw.conclude, gp=gpar(fontsize=9,font=8)))
  dev.off()
}


#########################################################################################################
option_list = list(
  make_option(c("-p", "--plink"), type="character", default=NULL, 
              help="Path to plink files index, it assumes a bed, bim and fam file with the same file name", metavar="character"),
  make_option(c("-x", "--xplink"), type="character", default=NULL, 
              help="Path to plink files index, it assumes a bed, bim and fam file with the same file name", metavar="character"),
  make_option(c("-m", "--mendel"), type="numeric", default=2, 
              help="Thrteshold for % of parent-sibling pairs to exclude SNPs", metavar="character"),
  make_option(c("-p", "--pairing"), type="character", default=2, 
              help=, metavar="character"),
  make_option(c("-o", "--out"), type="character", default="./famCheck_genotypeQC", 
              help="Output path to save report", metavar="character")
); 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

## testing
# opt <- list()
# opt$plink <- "/groups/umcg-ugli/tmp04/projects/merged_general_QC/second_QC_iteration/QCed_autosomes"
# opt$xplink <- "/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC/second_QC_iteration/X_QC/chr_X"
# opt$out <- "/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC/second_QC_iteration/"
# opt$mendel <- 1
# opt$pairing <- "/groups/umcg-ugli/tmp04/projects/merged_general_QC/5_Relatedness/family_corrections/output_files/pairing_corrected.dat"
# opt$famly <- "/groups/umcg-ugli/tmp04/projects/merged_general_QC/5_Relatedness/family_corrections/output_files/LifeLines_Corrected_family_info.dat"


mendelError.perc.threshold <- opt$mendel
out.mendel <- file.path(opt$out, "7_MendelErrors")
dir.create(out.mendel)
out.founder <- file.path(opt$out, "8_FounderStats")
dir.create(out.founder)

##### Run plink command to merge all chromosomes

all.chr.list <- list.files(opt$plink, pattern = ".bim", full.names = T)
all.chr.list <- gsub(all.chr.list, pattern = ".bim$", replacement = "")

all.chr.list.file <- file.path(out.mendel, "plink.all.chr.list")
write.table(as.data.frame(all.chr.list), file= all.chr.list.file, quote = FALSE, row.names = FALSE, col.names = FALSE)

all.chr.Xchr.list <- c(all.chr.list, opt$xplink)
all.chr.Xchr.list <- gsub(all.chr.Xchr.list, pattern = ".bim$", replacement = "")
all.chr.Xchr.list.file <- file.path(out.mendel, "plink.all.chr.Xchr.list")
write.table(as.data.frame(all.chr.Xchr.list), file= all.chr.Xchr.list.file, quote = FALSE, row.names = FALSE, col.names = FALSE)

###merge all chromosomes
merged.plink <- file.path(out.mendel, "merged")
plink.merge.chr.call <- paste0(
  "ml plink \n",
  "plink --merge-list ", all.chr.list.file, " --out ", merged.plink, " \n"
) 
system(plink.merge.chr.call)

merged.Xchr.plink <- file.path(out.mendel, "merged.Xchr")
plink.merge.chr.Xchr.call <- paste0(
  "ml plink \n",
  "plink --merge-list ", all.chr.Xchr.list.file, " --out ", merged.Xchr.plink, " \n"
) 
system(plink.merge.chr.Xchr.call)


##### Load the family info to the fam file from merged.plink
pairing <- fread(as.character(opt$pairing), data.table = FALSE)

fam.file <- paste0(merged.Xchr.plink, ".fam")
fam <- fread(fam.file, data.table=FALSE)

fam$id.clean <- gsub(fam$V1, pattern = "^OV365|^OV0365", replacement = "")
pairing$id.clean <- gsub(pairing$SAMPLE_IDENTIFIER_Rotterdam, pattern = "^OV365|^OV0365", replacement = "")
fam$PSEUDOIDEXT <- pairing$PSEUDOIDEXT[match(fam$id.clean, pairing$id.clean)]

# only pseudoID
family.info <- fread(opt$famly, data.table = FALSE)
#family.info <- read.delim(opt$famly, header = T)

new.fam <- data.frame(FID= family.info$FAM_ID[match(fam$PSEUDOIDEXT, family.info$PSEUDOIDEXT)],
                      IID= fam$PSEUDOIDEXT,
                      father= family.info$FATHER_PSEUDOID[match(fam$PSEUDOIDEXT, family.info$PSEUDOIDEXT)],
                      mother= family.info$MOTHER_PSEUDOID[match(fam$PSEUDOIDEXT, family.info$PSEUDOIDEXT)],
                      Sex= pairing$Gender[match(fam$PSEUDOIDEXT, pairing$PSEUDOIDEXT)]
)

new.fam$Sex <- ifelse(new.fam$Sex == "F", 2, 1)
#new.fam$Sex[which(is.na(new.fam$Sex))] <- 0
new.fam$genotype <- 1
new.fam$genotype[which(new.fam$IID %in% pairing$PSEUDOIDEXT)] <- 2


autosome.fam.file <- paste0(merged.plink, ".fam")
autosome.fam <- fread(autosome.fam.file, data.table=FALSE)

autosome.fam$id.clean <- gsub(autosome.fam$V1, pattern = "^OV365|^OV0365", replacement = "")
autosome.fam$PSEUDOIDEXT <- pairing$PSEUDOIDEXT[match(autosome.fam$id.clean, pairing$id.clean)]
new.autosome.fam <- data.frame(FID= family.info$FAM_ID[match(fam$PSEUDOIDEXT, family.info$PSEUDOIDEXT)],
                               IID= fam$PSEUDOIDEXT,
                               father= family.info$FATHER_PSEUDOID[match(fam$PSEUDOIDEXT, family.info$PSEUDOIDEXT)],
                               mother= family.info$MOTHER_PSEUDOID[match(fam$PSEUDOIDEXT, family.info$PSEUDOIDEXT)],
                               Sex= pairing$Gender[match(fam$PSEUDOIDEXT, pairing$PSEUDOIDEXT)]
)

new.autosome.fam$Sex <- ifelse(new.autosome.fam$Sex == "F", 2, 1)
new.autosome.fam$Sex[which(is.na(new.autosome.fam$Sex))] <- 0
new.autosome.fam$genotype <- 1
new.autosome.fam$genotype[which(new.autosome.fam$IID %in% pairing$PSEUDOIDEXT)] <- 2

#########################
##### write out the modified fam files with updated families from the famliy.info file. 
## Read pairing info in order to mark the samples that we have genotyped.
write.table(new.fam, file=fam.file, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(new.autosome.fam, file=autosome.fam.file, 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)



##### Extract only females for X chr founder analysis. 
##### With phenotypes from the pairing file generate a plink file for the X chromosome that would include on females. 
females.in.fam.file <- file.path(out.mendel, "females.in.fam.txt")
females.in.fam <- new.fam[which(new.fam$Sex == 2), c("FID","IID")]
write.table(females.in.fam, file=females.in.fam.file, row.names = FALSE, quote = FALSE, col.names = FALSE, sep="\t")

xchr.founders.plink.file <- file.path(out.founder, "xchr.founders")
xchr.founder.plink.call <- paste0(
  "ml plink \n",
  "plink --bfile ", merged.Xchr.plink, " \\",
  "--chr X \\",
  "--make-bed \\",
  "--keep ", females.in.fam.file, " \\",
  "--no-pheno"," \\",
  "--out ", xchr.founders.plink.file, " \n"
)
system(xchr.founder.plink.call)


####################################################################################################
##### Run plink command to get mendelian erros. 
plink.mendel.call <- paste0(
  "ml plink \n",
  "plink ", 
  "--bfile ",  merged.Xchr.plink,  " ", 
  "--no-pheno ",
  "--mendel ", 
  "--out ", out.mendel , "/cohort.autosome.Xchr",
  "\n"
)
system(plink.mendel.call)

## Run PEDSTATS to count the total number of parent siblings within our cohort.
## To do so the .fam file needs to be modified to be used and input file for PEDSTATS
## and extra column needs to be added to consider smaples which have been genotyped. 


#made.up.parents <- unique(c(grep("^S", fam$V4),
#                                     grep("^S", fam$V3)))
#no.parents <- which(fam$V3 == 0 | fam$V4 == 0)
#fam$genotype[unique(made.up.parents, no.parents)] <- 1

## write.out modified fam. 
pedstat.fam <- family.info[,c("FAM_ID", "PSEUDOIDEXT", "FATHER_PSEUDOID", "MOTHER_PSEUDOID", "GENDER1M2F")]
pedstat.fam$genotype <- ifelse(pedstat.fam$PSEUDOIDEXT %in% new.fam$IID, 2, 1)


##### double check pedigree for missing parents IDs

pedstat.fam.file <- file.path(out.mendel, "/pedstat.fam.txt")
write.table(pedstat.fam, file=pedstat.fam.file, sep="\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
ped.dat.file <- file.path(out.mendel, "/pedstat.fam.dat")
write.table(rbind(c("A","genotype")), file=ped.dat.file, sep="\t", row.names = FALSE, quote = FALSE, col.names = FALSE)


## system call for PED 
out.mendel.ped <- file.path(out.mendel, "/res.pedstats.txt")
ped.tool.path <- "/home/umcg-raguirregamboa/pedstats-0.6.12/pedstats"
ped.mendel.call <- paste0(
  "cd ",  out.mendel, " \n",
  ped.tool.path , " ", 
  "-p ",  pedstat.fam.file,  " ", 
  "-d ",  ped.dat.file,  " ", 
  "--pairs ", 
  "--affectedFor:genotype ", ">" , out.mendel.ped,
  "\n" 
)
system(ped.mendel.call)

## load results from PED and get the total number of parent sibling relationships within the pedigree 

ped.start.line <- system(paste0("grep -n ", "'^PAIR\\sSTATISTICS' ", out.mendel.ped, "\n"), intern = T)
ped.start.line <- as.numeric(gsub(ped.start.line, pattern = ":.*", replacement = ""))+1

ped.end.line <- system(paste0("grep -n ", "'^Pair\\sCounts\\sby\\sAffection\\sStatus:' ", out.mendel.ped, "\n"), intern = T)
ped.end.line <- as.numeric(gsub(ped.end.line, pattern = ":.*", replacement = ""))-2

## read only the selected lines to load the PEDStat results. 
ped.lines <- c(ped.start.line:ped.end.line)
ped.info <- lapply(ped.lines, function(i) scan(file=out.mendel.ped,  nlines = 1, skip = i, what="a")) 
ped.info <- ped.info[-which(sapply(ped.info, length ) == 0)]
ped.info.bind <- do.call(rbind, ped.info)

#number of parent sibling relationships
n.parent.child <- as.numeric(ped.info.bind[grep(ped.info.bind[,1], pattern = "Parent-Child"),2])


## plot results, distribution of mendial errors across SNPs. 

mendel.snp.file <- file.path(out.mendel , paste0("/cohort",".lmendel"))
mendel.snp <- fread(mendel.snp.file, data.table = FALSE)

snp.mendel.labels <- c("0", "1-5", "6-10", "11-59", "51-100", paste0("101-", max(mendel.snp$N)))
mendel.snp$binned.N <- cut(mendel.snp$N, breaks = c(-Inf, 1, 5, 10, 50, 100, max(mendel.snp$N)), labels = snp.mendel.labels)

mendel.snp.pData <- as.data.frame(table(mendel.snp$binned.N))
mendel.snp.pData$Percentage <- (mendel.snp.pData$Freq/sum(mendel.snp.pData$Freq))*100

mendel.snp$Percentage.PSp <- (mendel.snp$N/n.parent.child)



mendel.snp.plot <- ggplot(mendel.snp.pData, aes(x=Var1, y=Freq))+
  geom_bar(stat="identity", position="dodge")+
  geom_text(aes(label=Freq, y=Freq+5), position = position_dodge(0.9), vjust = 0)+
  ggtitle("Distribution of mendelian errors per SNP")+
  ylab("Number of SNPs")+
  xlab("Number of errors")+
  theme_classic()+
  theme(text = element_text(size=9, family='Helvetica'))

mendel.snp.perc.plot <- ggplot(mendel.snp.pData, aes(x=Var1, y=Percentage))+
  geom_bar(stat="identity", position="dodge")+
  geom_text(aes(label=paste0(round(Percentage, digits=2), "%"), y=Percentage+1.5), position = position_dodge(0.9), vjust = 0)+
  ggtitle("Distribution of mendelian errors per SNP")+
  ylab("Percentage of SNPs")+
  xlab("Number of errors")+
  theme_classic()+
  theme(text = element_text(size=9, family='Helvetica'))

n.snps.out <- sum(mendel.snp$Percentage.PSp*100 >= mendelError.perc.threshold , na.rm = T)
mendel.snp.perc.PSp.plot <- ggplot(mendel.snp, aes(x=Percentage.PSp))+
  # since scale_x_continuous(labels = percent) mendelError.perc.threshold needs to be /100 
  geom_vline(xintercept = mendelError.perc.threshold/100, color="red")+ 
  geom_histogram(stat="count", position="dodge", binwidth = 0.01)+
  ggtitle(label = paste0("% of PS (", format(n.parent.child, big.mark=",",scientific=FALSE) ,") that are affected per SNP\n",
                         format(n.snps.out, big.mark=",",scientific=FALSE), " SNPs do not pass the ", mendelError.perc.threshold, "% threshold", "\n", 
                         "(Autosomes and X chr)"))+
  scale_x_continuous(labels = percent)+
  ylab("Density")+
  xlab("PS pairs")+
  theme_classic()+
  theme(text = element_text(size=9, family='Helvetica'))

####################
mendel.plot.file <- file.path(out.mendel, "mendel.snp.plot.tiff")
tiff(mendel.plot.file,  
     width = 1200, height = 900*3, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(mendel.snp.plot, 
             mendel.snp.perc.plot, 
             mendel.snp.perc.PSp.plot, 
             nrow=3)
dev.off()

mendel.out.snps <- file.path(out.mendel, "mendel.snp.out.txt")
write.csv(mendel.snp[which(mendel.snp$Percentage.PSp*100 >= mendelError.perc.threshold),], file=mendel.out.snps)

########################################


mendel.sample.file <- file.path(out.mendel, paste0("/cohort",".imendel"))
mendel.sample <- fread(mendel.sample.file, data.table = FALSE)

sample.mendel.labels <- c("0", "1-100", "101-1,000", "1001-10,000", "10,001-100,000", paste0("10,0001-", format(max(mendel.sample$N))))
mendel.sample$binned.N <- cut(mendel.sample$N, breaks = c(-Inf, 1, 100, 1000, 10000, 100000, max(mendel.sample$N)), labels = sample.mendel.labels)

mendel.sample.pData <- as.data.frame(table(mendel.sample$binned.N ))
mendel.sample.pData$Percentage <- (mendel.sample.pData$Freq/sum(mendel.sample.pData$Freq))*100

mendel.sample.plot <- ggplot(mendel.sample.pData, aes(x=Var1, y=Freq))+
  geom_bar(stat="identity", position="dodge")+
  geom_text(aes(label=Freq, y=Freq+5), position = position_dodge(0.9), vjust = 0)+
  ggtitle("Distribution of mendelian errors per Sample")+
  ylab("Number of Samples")+
  xlab("Number of errors")+
  theme_classic()+
  theme(text = element_text(size=10, family='Helvetica'), axis.text.x = element_text(angle = 45, hjust = 1))

mendel.sample.perc.plot <- ggplot(mendel.sample.pData, aes(x=Var1, y=Percentage))+
  geom_bar(stat="identity", position="dodge")+
  geom_text(aes(label=paste0(round(Percentage, digits=2), "%"), y=Percentage+5), position = position_dodge(0.9), vjust = 0)+
  ggtitle("Distribution of mendelian errors per Sample")+
  ylab("Percentage of Samples")+
  xlab("Number of errors")+
  theme_classic()+
  theme(text = element_text(size=10, family='Helvetica'), axis.text.x = element_text(angle = 45, hjust = 1))



mendel.plot.file <- file.path(out.mendel, "mendel.sample.plot.tiff")
tiff(mendel.plot.file,  
     width = 1200, height = 1200*2, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(mendel.sample.plot, mendel.sample.perc.plot, nrow=2)
dev.off()


############################################################
############################################################
##### Run plink command to get stats within founders using 

## autosomes chromosomes only
plink.founderStats.call <- paste0(
  "ml plink \n",
  "plink ", 
  "--bfile ", merged.plink,  " ", 
  "--freq ", 
  "--hardy ",
  "--filter-founders ",
  "--out ", out.founder, "/FounderStats \n"
)

system(plink.founderStats.call)


plink.founderStats.Xchr.call <- paste0(
  "ml plink \n",
  "plink ", 
  "--bfile ", xchr.founders.plink.file,  " ", 
  "--freq ", 
  "--hardy ",
  "--filter-founders ",
  "--out ", out.founder, "/Xchr.FounderStats \n"
)

system(plink.founderStats.Xchr.call)

## x chromosome only, including only females. 

maf.file <- file.path(out.founder, "FounderStats.frq")
maf.dat <- fread(maf.file, data.table = FALSE)
maf.dat$CHR <-as.factor(maf.dat$CHR)

hw.file <- file.path(out.founder, "FounderStats.hwe")
hw.dat <- fread(hw.file, data.table = FALSE)
hw.dat$CHR <-as.factor(hw.dat$CHR)
hw.dat$MAF <- maf.dat$MAF[match(hw.dat$SNP,maf.dat$SNP)]


xchr.maf.file <- file.path(out.founder, "Xchr.FounderStats.frq")
xchr.maf.dat <- fread(xchr.maf.file, data.table = FALSE)
xchr.maf.dat$CHR <-as.factor(xchr.maf.dat$CHR)

xchr.hw.file <- file.path(out.founder, "Xchr.FounderStats.hwe")
xchr.hw.dat <- fread(xchr.hw.file, data.table = FALSE)
xchr.hw.dat$CHR <-as.factor(xchr.hw.dat$CHR)
xchr.hw.dat$MAF <- xchr.maf.dat$MAF[match(xchr.hw.dat$SNP,xchr.maf.dat$SNP)]


HW.MAF.dist.plot(maf.dat= maf.dat, hw.dat=hw.dat, out=out.founder, name="Autosome")

xchr.maf.dist.plot.all <- ggplot(xchr.maf.dat, aes(x=MAF))+
  stat_density(color="black",  position= "identity", geom= "line")+
  ggtitle("MAF distribution")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

xchr.hw.maf.dist.plot.all <- ggplot(xchr.hw.dat[which(xchr.hw.dat$MAF > 0.01),], aes(x=-log10(P)))+
  stat_density(position= "identity", geom= "line")+
  geom_vline(xintercept = 6)+
  xlim(c(0, 20))+
  ggtitle("HW pVal distribution of \n SNPs with a MAF > 0.01")+
  xlab("-log10(HW-P)")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

### printing out the maf.hw distributions
xchr.maf.hw.title <- paste0("MAF and HW distributions from Founders","\n", "X_chr"," ", date())
xchr.maf.hw.conclude <-paste0(sum(xchr.hw.dat$P<0.0001)," markers excluded ", (nrow(xchr.hw.dat)), " remaining.")
xchr.maf.hw.tiff.file <- file.path(out.founder, paste0("X_chr", "MAFHW.density.tiff"))
tiff(xchr.maf.hw.tiff.file, width = 2000, height = 1500, units = "px", res = 300, compression = "lzw")
grid.arrange(xchr.maf.dist.plot.all, xchr.hw.maf.dist.plot.all, 
             top=xchr.maf.hw.title, nrow=1,
             bottom=textGrob(xchr.maf.hw.conclude, gp=gpar(fontsize=9,font=8)))
dev.off()



final.list.excluded.snps <- unique(c(xchr.hw.dat$SNP[which(xchr.hw.dat$P<0.0001)],
                                     hw.dat$SNP[which(hw.dat$P<0.0001)],
                                     mendel.snp$SNP[which(mendel.snp$Percentage.PSp >= 0.01)]))

final.list.excluded.snps.file <- file.path(out.founder, "founder.mendelErrors.final.list.excluded.snps")
write.table(final.list.excluded.snps, file=final.list.excluded.snps.file,col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)


#####

autoQCed="/groups/umcg-ugli/tmp04/projects/merged_general_QC/second_QC_iteration/QCed_autosomes"
out_final="/groups/umcg-ugli/tmp04/projects/merged_general_QC/second_QC_iteration/final_QCed_autosomes_X"

mkdir ${out_final}

for chr in {1..22} "XY"
do

### create plink files and call_rate stats for individuals and SNPs
plink \
--bfile ${autoQCed}/chr_${chr} \
--exclude "/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC/second_QC_iteration//8_FounderStats/founder.mendelErrors.final.list.excluded.snps" \
--make-bed \
--out ${out_final}/chr_${chr}

done


