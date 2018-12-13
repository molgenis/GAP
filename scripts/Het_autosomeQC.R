###################################
### Heterozygosity filtering and report for autosomal QC
### date: 12-12-2018
### version: 0.01
### authors: EL
###################################

## example run 

#RScript Het_autosomeQC.R -i ***/autosome_QC
#  -o repout


#########################################################################################################
#########################################################################################################
### Functions and libraries
#########################################################################################################
#########################################################################################################

###packages
library(ggplot2)
library(optparse)
library(gridExtra)

##Arguments
#########################################################################################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input path ", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="Output path to save report", metavar="character")
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

het.file <- file.path(opt$input,"autosomal.het")

hetdata<-read.table(het.file, header = T)
hetdata$OHET<-(hetdata$N.NM.- hetdata$O.HOM.)

Excluded<-hetdata[abs(hetdata$OHET-mean(hetdata$OHET))>4*sd(hetdata$OHET),c(1,2)]
Included<-hetdata[abs(hetdata$OHET-mean(hetdata$OHET))<=4*sd(hetdata$OHET),c(1,2)]


excl.file <- file.path(opt$input,"Excluded.het")
incl.file <- file.path(opt$input,"Included.het")

write.table(Excluded,excl.file, quote=F,row.names = F)
write.table(Included,incl.file, quote=F,row.names = F)

density.het<-ggplot(hetdata, aes(x=OHET))+
             geom_density()+
             ggtitle("Heterozigosity distribution by sample")+
             theme_bw()+
             xlab("Heterozygosity")+
             geom_vline(xintercept=mean(hetdata$OHET)-4*sd(hetdata$OHET))+
             geom_vline(xintercept=mean(hetdata$OHET)+4*sd(hetdata$OHET))+
             theme(text=element_text(size=10, family="Helvetica"))

hetero.density.file <- file.path(output, "hetero.density.tiff")
tiff(hetero.density.file,  
     width = 1500, height = 1500, 
     units = "px", res = 300, compression = "lzw")
grid.arrange(density.het, ncol=1)
dev.off()
