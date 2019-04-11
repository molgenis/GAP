###################################
### Sex check evaluation
### date: 03-04-2019
### version: 0.01
### authors: EL - RAG
###################################
### New
###################################

library(tidyverse)
library(data.table)
library(optparse)
library(gridExtra)
library(viridis)
library(RColorBrewer)

# opt <- list()
# opt$input <- "/groups/umcg-aad/tmp04/umcg-elopera/X_QC_RD_part10/imputed.sexcheck"
# opt$phenotypes <- "/groups/umcg-ugli/tmp04/projects/UGLI_RD_part10/run01/results/UGLI_RD_part10.csv"
# opt$output <- "/groups/umcg-wijmenga/tmp04/umcg-raguirre/pln_ugli/genotypeQC_test/x_chr_part10"

#########################################################################################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input path to perform QC", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="genoQC_Report", 
              help="Output path to save report", metavar="character"),
  
  make_option(c("-p", "--phenotypes"), type="character", default=NULL, 
              help="path to sample info spreadsheet", metavar="character")
); 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied as input path", call.=FALSE)
}

output <- as.character(opt$out)
#output <- file.path(output, "plots")
#dir.create(output, recursive = T)

#########################################################################################################
### Main
#########################################################################################################

## Read phenotypes, we are currently using the info from the spreadsheets.

phenos <- fread(opt$phenotypes, data.table = FALSE)

plink.sex <- fread(opt$input, data.table=FALSE)
plink.sex$pheno.sex <- phenos$Gender[match(plink.sex$IID, phenos$Sample_ID)]
plink.sex$plink.sex <- ifelse(plink.sex$SNPSEX == 0, NA, ifelse(plink.sex$SNPSEX == 1, "M", "F"))

plink.sex$sex.concordance <- plink.sex$pheno.sex == plink.sex$plink.sex

## Evaluate sex concordance by plate
plink.sex$plate.id <- gsub(plink.sex$IID, pattern = "OV365", replacement = "")
plink.sex$plates <- gsub(plink.sex$IID, pattern = "OV365", replacement = "")
plink.sex$plates <- gsub(plink.sex$plates, pattern = "OV0365", replacement = "")
plink.sex$plates <- gsub(plink.sex$plates, pattern = "_.*", replacement = "")

plink.sex$plates.col <- gsub(plink.sex$IID, pattern = ".*_", replacement = "")
plink.sex$plates.col <- gsub(plink.sex$plates.col, pattern = "[A-Z]", replacement = "")

plink.sex$plates.row <- gsub(plink.sex$IID, pattern = ".*_", replacement = "")
plink.sex$plates.row <- gsub(plink.sex$plates.row, pattern = "[0-9][0-9]", replacement = "")


all.plates <- unique(plink.sex$plates)
all.plates.concordance <- lapply(all.plates, function(x){sum(plink.sex[which(plink.sex$plates == x), "sex.concordance"], na.rm = T)})
names(all.plates.concordance) <- all.plates

melt.all.plates.concordance <- melt(all.plates.concordance)
melt.all.plates.concordance$n.samples.plate <- table(plink.sex$plates)[all.plates]

melt.all.plates.concordance$percentage <- melt.all.plates.concordance$value/melt.all.plates.concordance$n.samples.plate

melt.all.plates.concordance$n.samples.wrong <- melt.all.plates.concordance$n.samples.plate - melt.all.plates.concordance$value

conc.group.labels <- rev(c("100%", "99-95%", "95-90%", "90-80%", "80-60%", "60-0%"))
melt.all.plates.concordance$grouped <- cut(melt.all.plates.concordance$percentage, right = F,
                                           breaks = c(Inf,1 ,0.95, 0.90, 0.80, 0.60, 0), 
                                           labels = conc.group.labels)

melt.all.plates.concordance$grouped <- factor(melt.all.plates.concordance$grouped, levels= conc.group.labels)
bar.colors <- brewer.pal(n = 11, name = "Spectral")[c(1,3,5,7,9,11)]
names(bar.colors) <- conc.group.labels

concordance.per.plate.plot <- ggplot(melt.all.plates.concordance, aes(x=L1, y=percentage))+
                                geom_bar(stat = "identity", aes(fill= grouped))+
                                scale_fill_manual(name="", values = bar.colors)+
                                theme_classic()+
                                ggtitle("Concordance between reported sex and genotyped imputed sex")+
                                ylab("Concordance")+
                                xlab("Plates")+
                                theme(legend.position = "bottom",
                                      axis.text.x = element_text(angle = 45, hjust = 1),
                                      text = element_text(family = "Helvetica", size = 10)
                                )

bar.plot.width.factor <- 123
concordance.per.plate.plot.file <- file.path(output, paste0("concordance.per.plate.plot.tiff"))
### Save plot and table with information. 
tiff(concordance.per.plate.plot.file, width = bar.plot.width.factor*length(all.plates), height = 900, 
     units = "px", res = 300, compression = "lzw")
plot(concordance.per.plate.plot)
dev.off()


colnames(melt.all.plates.concordance)[1:2] <- c("Plate", "n.concordant.samples")
all.plates.concordance.file <-  file.path(output, paste0("all.plates.concordance.txt"))
write.table(melt.all.plates.concordance, file= all.plates.concordance.file, 
            sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE )



################################################################################
################################################################################
### Plot each plate present in the batch 
### a plot to visualize their position on a plate will be made

#Create a blank dataframe with all possible samples given by all.plates 
row.numbers <- c(paste0(c("0"),c(1:9)), c(10:12))
plate.df <- data.frame(Plate= rep(all.plates, each=96), 
                       row= rep(LETTERS[c(1:8)], each=length(c(1:12)), times=length(all.plates)),
                       col= rep(row.numbers, times= length(all.plates)*length(LETTERS[c(1:8)]))
                       )
plate.df$plate.id <- paste0(plate.df$Plate, "_",plate.df$row,  plate.df$col)

plate.df$col <- factor(plate.df$col, levels = rev(levels(plate.df$col)))

plink.sex$sex.concordance[which(is.na(plink.sex$sex.concordance))] <- "Failed genetic imputation"
plate.df$Status <- "No data in batch"
plate.df$Status[match(plink.sex$plate.id,plate.df$plate.id)] <- plink.sex$sex.concordance
plate.df$Status[which(plate.df$Status == "FALSE")] <- "Non concordant"
plate.df$Status[which(plate.df$Status == "TRUE")] <- "OK"


tile.plot.height.factor <- 150

tile.colors <- viridis_pal(option = "A")(4)
names(tile.colors) <- levels(as.factor(plate.df$Status))

tile.plot <- ggplot(plate.df, aes(x=row, y= col))+
  geom_tile(aes(fill= Status, width=0.9, height=0.9),color="gray")+
  facet_wrap(~Plate, ncol= 5, scales = "free")+
  scale_fill_manual(values = tile.colors)+
  theme_bw()+
  ggtitle("Sex concordance across plates")+
  ylab("")+
  xlab("")+
  theme(legend.position = "bottom", strip.background = element_blank(),
    text = element_text(family = "Helvetica", size = 9))

plot.file <- file.path(output, paste0("plate_sex_check.tiff"))

tiff(plot.file, width = 2000, height = tile.plot.height.factor*length(all.plates), 
     units = "px", res = 300, compression = "lzw")
plot(tile.plot)
dev.off()


plate.df$imputed.sex <- NA
plate.df$imputed.sex[match(plink.sex$plate.id,plate.df$plate.id)] <- plink.sex$plink.sex
plate.df$imputed.sex.F <- NA
plate.df$imputed.sex.F[match(plink.sex$plate.id,plate.df$plate.id)] <- plink.sex$F
plate.df$reported.sex <- NA 
plate.df$reported.sex[match(plink.sex$plate.id,plate.df$plate.id)] <- plink.sex$pheno.sex

all.samples.concordance.file <-  file.path(output, paste0("all.samples.concordance.txt"))
write.table(plate.df, file= all.samples.concordance.file, 
            sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE )


