###################################
### Sex check evaluation
### date: 09-04-2019
### version: 0.01
### authors: EL - RAG
###################################
## New
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

# opt <- list()
# opt$input <- "/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC/X_QC/0_pre/imputed.sexcheck"
# opt$phenotypes <- "/groups/umcg-aad/tmp04/umcg-elopera/ugli_blood_gsa/pairing.dat"
# opt$output <- "/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC/plots"
# opt$dup<-"/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC/5_Relatedness/proc2/equal.samples"

##function to substract the unique part of the sample ID
codextract <- function(x) {
  x <- as.character(x)
  return(substr(x,nchar(x)-9,nchar(x)))
}

#########################################################################################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input path to perform QC", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="genoQC_Report", 
              help="Output path to save report", metavar="character"),
  
  make_option(c("-p", "--phenotypes"), type="character", default=NULL, 
              help="path to sample info spreadsheet", metavar="character"),
  
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

phenos <- fread(opt$phenotypes, data.table = FALSE)
plink.sex <- fread(opt$input, data.table=FALSE)

##harmonize names
colnames(phenos)<-c("PSEUDOIDEXT","Sample_ID", "Age", "BIRTHYEAR", "Gender")
phenos$Sample_ID<-sapply(phenos$Sample_ID,FUN=codextract)
plink.sex$IID<-sapply(plink.sex$IID,FUN=codextract)


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

#### while loop to split the data in chunks of 50 plates and plot by 50 plates (or less)

limit<-0
n=0
while(limit<length(all.plates)){
  limit<-ifelse(50*(n+1)>length(all.plates),length(all.plates),50*(n+1))
  
  concordance.per.plate.plot <- ggplot(melt.all.plates.concordance[c(((50*n)+1):limit),], aes(x=L1, y=percentage))+
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
  
  concordance.per.plate.plot.file <- file.path(output, paste0("06.sexcheck.barplot_",n,".tiff"))
  ### Save plot and table with information. 
  tiff(concordance.per.plate.plot.file, width = bar.plot.width.factor*length(all.plates[((50*n)+1):limit]), height = 900, 
       units = "px", res = 300, compression = "lzw")
  plot(concordance.per.plate.plot)
  dev.off()
  
  
  ################################################################################
  ### Plot each plate present in the batch 
  ### a plot to visualize their position on a plate will be made
  
  #Create a blank dataframe with all possible samples given by all.plates 
  row.numbers <- c(paste0(c("0"),c(1:9)), c(10:12))
  plate.df <- data.frame(Plate= rep(all.plates[((50*n)+1):limit], each=96), 
                         row= rep(LETTERS[c(1:8)], each=length(c(1:12)), times=length(all.plates[((50*n)+1):limit])),
                         col= rep(row.numbers, times= length(all.plates[((50*n)+1):limit])*length(LETTERS[c(1:8)]))
  )
  plate.df$plate.id <- paste0(plate.df$Plate, "_",plate.df$row,  plate.df$col)
  
  plate.df$col <- factor(plate.df$col, levels = rev(levels(plate.df$col)))
  
  plink.sex$sex.concordance[which(is.na(plink.sex$sex.concordance))] <- "Failed genetic imputation"
  plate.df$Status <-ifelse(plate.df$plate.id %in% plink.sex$plate.id,
                           plink.sex$sex.concordance[match(plate.df$plate.id,plink.sex$plate.id)],"No data in batch or QC-removed" )
  plate.df$Status[which(plate.df$Status == "FALSE")] <- "Non concordant"
  plate.df$Status[which(plate.df$Status == "TRUE")] <- "OK"
  
  
  tile.colors <- c(viridis_pal(option = "A")(4)[1],"beige",viridis_pal(option = "A")(4)[3],"#99FFCC")
  names(tile.colors) <- c("Failed genetic imputation", "No data in batch or QC-removed", "Non concordant", "OK")
  
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
  
  plot.file <- file.path(output, paste0("06.sexcheck.plates_",n,".tiff"))
  
  tiff(plot.file, width = 2000, height = tile.plot.height.factor*length(all.plates[((50*n)+1):limit]), 
       units = "px", res = 300, compression = "lzw")
  plot(tile.plot)
  dev.off()
  
  
  plate.df$imputed.sex <-ifelse(plate.df$plate.id %in% plink.sex$plate.id,
                                plink.sex$plink.sex[match(plate.df$plate.id,plink.sex$plate.id)],NA )
  plate.df$imputed.sex.F <-ifelse(plate.df$plate.id %in% plink.sex$plate.id,
                                plink.sex$F[match(plate.df$plate.id,plink.sex$plate.id)],NA )
  plate.df$reported.sex <-ifelse(plate.df$plate.id %in% plink.sex$plate.id,
                                plink.sex$pheno.sex[match(plate.df$plate.id,plink.sex$plate.id)],NA )

  all.samples.concordance.file <-  file.path(output, paste0("all.samples.concordance_",n,".txt"))
  write.table(plate.df, file= all.samples.concordance.file, 
              sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE )
  
  ##read duplicate data and plot if existing
  if (file.exists(opt$dup)) {
    dups <- fread(opt$dup, data.table = FALSE)
    
    ##column to mark as duplicate
    dups$ID1<-sapply(dups[,1],FUN=codextract)
    dups$ID2<-sapply(dups[,2],FUN=codextract)
    plate.df$markdup<-ifelse((plate.df$plate.id %in% dups$ID2|plate.df$plate.id %in% dups$ID1),1,NA)
    ##column to ascertain the duplicate pair
    dups$groupdup <- 1:nrow(dups)
    dup.pairs <- melt(dups, id.vars = "groupdup", measure.vars=c("ID1", "ID2"))
    plate.df$exact.dup<-dup.pairs$groupdup[match(plate.df$plate.id,dup.pairs$value)]
    
    ## sex check with amrked duplicates
    tile.plot.marked <- ggplot(plate.df, aes(x=row, y= col))+
      geom_tile(aes(fill= Status,color=as.factor(markdup), width=0.8, height=0.8),size=0.8)+
      facet_wrap(~Plate, ncol= 5, scales = "free")+
      scale_fill_manual(values = tile.colors)+
      theme_bw()+
      ggtitle("Sex concordance and duplicate plates")+
      ylab("")+
      xlab("")+
      guides(color=F)+
      scale_color_manual(values="black",na.value="gray")+
      theme(legend.position = "bottom", strip.background = element_blank(),
            text = element_text(family = "Helvetica", size = 9))
    
    plot.file <- file.path(output, paste0("sex_dup_marked_",n,".tiff"))
    
    tiff(plot.file, width = 2000, height = tile.plot.height.factor*length(all.plates[((50*n)+1):limit]), 
         units = "px", res = 300, compression = "lzw")
    plot(tile.plot.marked)
    dev.off()
    
    ## sex check with oaired duplicates
    paired.colors <- rep(c("#a50026","#d73027","#f46d43","#fdae61", "#fee090","#762a83","#1b7837","#7fbc41","#74add1","#4575b4","#313695"),20)
    
    tile.plot.paired <- ggplot(plate.df, aes(x=row, y= col))+
      geom_tile(aes(fill= Status,color=as.factor(exact.dup), width=0.8, height=0.8),size=0.8,na.rm = T)+
      facet_wrap(~Plate, ncol= 5, scales = "free")+
      scale_fill_manual(values = tile.colors)+
      theme_bw()+
      ggtitle("Sex concordance and paired duplicated samples")+
      ylab("")+
      xlab("")+
      scale_color_manual(values=paired.colors,na.value="gray")+
      guides(color=F)+
      theme(legend.position = "bottom", strip.background = element_blank(),
            text = element_text(family = "Helvetica", size = 9))
    
    plot.file <- file.path(output, paste0("sex_dup_paired_",n,".tiff"))
    
    tiff(plot.file, width = 2000, height = tile.plot.height.factor*length(all.plates[((50*n)+1):limit]), 
         units = "px", res = 300, compression = "lzw")
    plot(tile.plot.paired)
    dev.off()
  }
  n=n+1
}


colnames(melt.all.plates.concordance)[1:2] <- c("Plate", "n.concordant.samples")
all.plates.concordance.file <-  file.path(output, paste0("all.plates.concordance.txt"))
write.table(melt.all.plates.concordance, file= all.plates.concordance.file, 
            sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE )
