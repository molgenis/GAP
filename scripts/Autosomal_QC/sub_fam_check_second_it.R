###################################
### Family relationship check
### date: 09-10-2019
### version: 6.0
### authors: EL - RAG
###################################
## New
## 09-10-2019
## included function to make dummy and map file within the code
## included lines to change the names of columns to the usage in the functions. (must me reviewed if reference files change from pairing to sample-info)
## 13-09-2019
## added line to create output directory before the first plot
## added call to function to make dummy files if it doesnt exist already
## added --relations options to add a file with the list of FamilyIDs with relations that are aoutside KING understanding. namely, known twins added here, can change the error flag, from the .kin output, for plotting
## filter for at least "two individuals" for the plot now refers to "two individuals in the pedigree" instead of two individuals with genetic ionformation
## added MZ/duplicates to the first degree family filter before plotting
## unrelated family errors are also plotted now
## added plotting for the people not present in the pedigree file with a dummy family name
## 30-07-2019
## added command to attach dummy pedigree for individuals with genetic information but not pedigree
## removed function to plot pedigrees with ggplot
## added unrelated flags from the family errors to the Cranefoot plot for new_found
###################################
# 16-04-2018
# Added file.exists checks 
# changed king input parameters to include -fam, -bed and -bim files to allow different names and not have to duplicate files. 
###################################

############# Packages ############
library(data.table)
library(tidyverse)
library(optparse)
library(gridExtra)
library(reshape)

##cluster test
## opt<-list()
## opt$code<-"/groups/umcg-aad/tmp04/umcg-elopera/ugli_blood_gsa/pairing.dat"
## opt$info<-"/groups/umcg-aad/tmp04/umcg-elopera/ugli_blood_gsa/LifeLines_families_info_withspouse.dat"
## opt$plink<-"/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC/5_Relatedness/proc/full_data.no.dup"
## opt$out<-"/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC/testdata/test_5/plots/"
## opt$workdir<-"/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC/5_Relatedness/proc2/"
## opt$king<-"/groups/umcg-aad/tmp04/umcg-elopera/tools/KING/king"
## opt$dummy<-"/groups/umcg-aad/tmp04/umcg-elopera/ugli_blood_gsa/corrected_v7_dummy_allped" 
## opt$crane<-"/groups/umcg-aad/tmp04/umcg-elopera/tools/Cranefoot/example/cranefoot"
## opt$makeped<- TRUE

######### Functions ###############

#######function to substract the unique part of the sample ID
codextract <- function(x) {
  x <- as.character(x)
  return(substr(x,nchar(x)-9,nchar(x)))
}
################
## Make dummy famfile from the pedigree
make_dummy_files<-function(pedfile,plink_path,dummy_file){
  ######################
  ##read files
  info.table <- fread(pedfile, data.table=F,sep = "\t")
  snps<-tail(data.table::fread(paste0(plink_path,".bim"), data.table=F))
  ##make ped file
  complete.ped<-info.table[,c(2,1,3,4,5)]### from pedigree file
  complete.ped$V6<-"-9"
  num.snps<- nrow(snps)
  for (snp in 1:(2*num.snps)){
    complete.ped[paste0("C",snp)]<-0
  }
  write.table(snps,paste0(dummy_file,".map"),quote = F,row.names = F,col.names = F)
  write.table(complete.ped,paste0(dummy_file,".ped"),quote = F,row.names = F,col.names = F)
}

################################################################################
### Pedigree function with Cranefoot
pedigree_crane<-function(ls,
                         pedfile,king,
                         king_0,unexpected=TRUE,
                         Fs.grade=FALSE,
                         fam_batch,crane.path){
  
  if(Fs.grade==TRUE){rel_vector=c("PO","FS")} else {rel_vector=king_0$InfType}
  
  ####generate single family (or 2 families) file
  if(unexpected==TRUE) {
    ###processing from kin0 information
    IDS<-unique(as.vector(rbind(king_0[which((king_0$FID1==ls|king_0$FID2==ls) & 
                                               king_0$InfType %in% rel_vector),"ID2"], 
                                king_0[which((king_0$FID1==ls|king_0$FID2==ls) & 
                                               king_0$InfType %in% rel_vector),"ID1"] )))
    add_FID<-as.vector(fam_batch[fam_batch$PSEUDOIDEXT %in% IDS,"FAM_ID"])
    family.ped<-pedfile[which(pedfile$FAM_ID==ls|pedfile$FAM_ID %in% add_FID),c(1,2,3,4,5,7) ] 
  } else {
    ###processing from kin information
    family.ped<-pedfile[pedfile$FAM_ID==ls ,c(1,2,3,4,5,7) ]
  }
  names(family.ped)<-c("IID","FAM_ID","FATHER_PSEUDOID","MOTHER_PSEUDOID","GENDER1M2F","AGE")
  
  ##retrieve parents to make them individuals
  extra.IDs <- unique(c(family.ped[,"FATHER_PSEUDOID"],family.ped[,"MOTHER_PSEUDOID"]))
  `%!in%` = Negate(`%in%`)## create negation of %in% function
  extra.IDs <- extra.IDs[ which(extra.IDs %!in% family.ped$PSEUDOIDEXT)]   # remove extra IDs (parents) which are already in the final.fam$PSEUDOIDEXT
  extra.IDs.sex <- ifelse(extra.IDs %in% family.ped$FATHER_PSEUDOID, 1, 2) #define sex of new inidivudals (parents)
  #fill data for parents
  FAM_ID<-ifelse(extra.IDs.sex==1,family.ped$FAM_ID[match(extra.IDs,family.ped$FATHER_PSEUDOID)],
                 family.ped$FAM_ID[match(extra.IDs,family.ped$MOTHER_PSEUDOID)])
  AGE<-ifelse(extra.IDs %in% family.ped$IID, family.ped$AGE[match(extra.IDs,family.ped$IID)]," ")
  
  extra.fam <- data.frame(extra.IDs,
                          FAM_ID,
                          rep(0, length(extra.IDs)),
                          rep(0, length(extra.IDs)),
                          extra.IDs.sex,
                          AGE)
  #unite parents with offspring in a unique database
  colnames(extra.fam) <- colnames(family.ped) 
  final.fam <- rbind(family.ped, extra.fam)
  final.fam<-final.fam[final.fam$IID!=0,]
  
  ##Generate phenotypes files (sex, genetic info)
  #final.fam <-family.ped
  final.fam$GENDER1M2F<-ifelse(final.fam$GENDER1M2F==1,"M","F")## gender info
  final.fam$Genetic_info<-ifelse(final.fam$IID %in% fam_batch$PSEUDOIDEXT,"550077","999999")## genetic info colored "purple"
  
  #### include relations from kin or kin0
  if(unexpected==TRUE) {
    ###for kin 0
    relations<-king_0[which(king_0$InfType!="UN" & 
                              (king_0$FID2 %in% unique(final.fam$FAM_ID) & king_0$FID1 %in% unique(final.fam$FAM_ID)) &
                              king_0$InfType %in% rel_vector), ]  
    
    relations2<-king[which( king$InfType=="UN" 
                            &  king$Error==1 & ((king$ID1 %in% unique(final.fam$IID)| king$ID2 %in% unique(final.fam$IID)) )), c(1,2,1,3,4,7:15) ]
    
    colnames(relations2)<-colnames(relations)
    relations<-rbind(relations,relations2)%>%sort_df(vars="InfType")
    
    if (nrow(relations)==0){gen.fam<-"incomplete family information"}
    else {
      relations$Familial_errors<- 1:nrow(relations)
      gen.fam<-gather(relations, key="mode", "IID", "ID1", "ID2")## group each relation
      ###stablish error filters
      err.row.index<-c(seq(1:nrow(gen.fam)))
      err.col.index<-c("FID1","FID2","Familial_errors","IID","message","InfType")
      event<-"new_found."
    }
  } else {
    relations<-king[which( (king$ID1 %in% unique(final.fam$IID)| king$ID2 %in% unique(final.fam$IID)) ), ]
    relations$Familial_errors<- 1:nrow(relations)
    gen.fam<-gather(relations, key="mode", "IID", "ID1", "ID2")## group each relation
    ###stablish error filters
    err.row.index<-which(gen.fam$Error.tag=="Error")
    err.col.index<-c("FID","Familial_errors","IID","Error.tag","message","InfType")
    event<-"error."
  }
  if (gen.fam=="incomplete family information") {return(paste0("Family ",ls," incomplete"))} ####incomplete will be reported when there is no pedigree onformation  for genetic samples
  else {
    
    gen.fam$GENDER1M2F<-family.ped$GENDER1M2F[match(gen.fam$IID,family.ped$IID)]
    ##add  text with the information of relationships
    gen.fam$message<-ifelse(gen.fam$mode=="ID1",
                            paste0(gen.fam$InfType," with ",gen.fam[which(gen.fam$Familial_errors==gen.fam$Familial_errors & gen.fam$mode=="ID2"), "IID"]),
                            paste0(gen.fam$InfType," with ",gen.fam[which(gen.fam$Familial_errors==gen.fam$Familial_errors & gen.fam$mode=="ID1"), "IID"]))
    
    ###separate familial error in a different table
    Gen_errors<-gen.fam[err.row.index,err.col.index]
    
    ####make text columns for information of relationships
    n.text.col<-max(table(Gen_errors$IID))
    fillrels<-function(x) {
      k<-c("IID"=x,
           Gen_errors[Gen_errors$IID==x ,"message"], rep("",n.text.col-length(Gen_errors[Gen_errors$IID==x,"message"])))
      return(k)
    }
    mesdf<-lapply(unique(Gen_errors$IID), FUN= fillrels )
    mesdf<-data.frame(matrix(unlist(mesdf), nrow=length(mesdf), byrow=T))
    names(mesdf)<-c("IID",paste0("Genetic_relationship_",seq(1:n.text.col))) 
    mesdf$IID<-as.character(mesdf$IID)
    Gen_errors$IID<-as.character(Gen_errors$IID)
    Gen_errors<-right_join(Gen_errors,mesdf,by="IID")
    arrowID<-unique(as.vector(Gen_errors[Gen_errors$InfType!="UN","IID"]))
    Gen_errors[Gen_errors$IID %!in% arrowID,"Familial_errors"]<-" "
    ###create directory to save the result plots
    dir.create(paste0(opt$out,"Crane_fam_scripts"),recursive = TRUE, showWarnings = FALSE)
    crane.dir<-paste0(opt$out,"Crane_fam_scripts/")
    
    ## write necessary input files for cranefoot
    write.table(final.fam,paste0(opt$workdir,"family.ped"),quote = F,row.names = F,col.names = T,sep = "\t")
    write.table(Gen_errors,paste0(opt$workdir,"fam.error"),quote = F,row.names = F,col.names = T,sep = "\t")
    write.table(final.fam,paste0(opt$workdir,"phenotype_1"),quote = F,row.names = F,col.names = T,sep = "\t")
    
    ## create config file for cranefoot
    config.file<-data.frame(
      cbind(
        
        c("PedigreeFile","PedigreeName","NameVariable","FatherVariable",
          "MotherVariable","GenderVariable","ColorVariable","TextVariable","TextVariable","ArrowVariable",
          rep("TextVariable",n.text.col)),
        
        c(paste0(opt$workdir,"family.ped"),paste0(crane.dir,"07.",event ,ls,".ped"),
          "IID","FATHER_PSEUDOID","MOTHER_PSEUDOID","GENDER1M2F","Genetic_info","IID","AGE",
          "Familial_errors", paste0("Genetic_relationship_",seq(1:n.text.col)) ),
        
        c("","","","","",paste0(opt$workdir,"phenotype_1"),paste0(opt$workdir,"phenotype_1"),
          paste0(opt$workdir,"phenotype_1"),paste0(opt$workdir,"phenotype_1"),paste0(opt$workdir,"fam.error"), 
          rep(paste0(opt$workdir,"fam.error"),n.text.col) )
        
      )
    )
    write.table(config.file,paste0(opt$workdir,"CFG"),quote = F,row.names = F,col.names = F,sep = "\t")
    
    ## create system call for cranefoot
    crane.system.call <- paste0(c(opt$crane)," ",paste0(opt$workdir,"CFG"))
    system(crane.system.call)
    return(paste0("Family ",ls," done"))
  }
}

################################################################################

#########################################################################################################
option_list = list(
  make_option(c("-p", "--plink"), type="character", default=NULL, 
              help="Path to plink files index, it assumes a bed, 
              bim and fam file with the same file name", 
              metavar="character"),
  
  make_option(c("-c", "--code"), type="character", default=NULL, 
              help="Path to pairing ID's file", metavar="character"),
  
  make_option(c("-i", "--info"), type="character", default=NULL, 
              help="Phenotype and pedigree information file", metavar="character"),
  
  make_option(c("-k", "--king"), type="character", default=NULL, 
              help="Path to excecutable king", metavar="character"),
  
  make_option(c("-d", "--dummy"), type="character", default=NULL, 
              help="Path to dummyfiles index", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="./famCheck_genotypeQC", 
              help="Output path to save report", metavar="character"),
  
  make_option(c("-C", "--cranepath"), type="character", 
              default="/groups/umcg-aad/tmp04/umcg-elopera/tools/Cranefoot/example/cranefoot", 
              help="path to cranefoot executable file", metavar="character"),
  
  make_option(c("-M", "--makeped"), type="logical", 
              default=FALSE, 
              help="TRUE if familial data is complete to create", metavar="character"),
  
  make_option(c("-r", "--relations"), type="character", 
              default=NULL, 
              help="path to list of FAM-IDs with  known relations (twins, etc)", 
              metavar="character"),
  
  make_option(c("-w", "--workdir"), type="character", 
              default="opt$out", help="processing directory", metavar="character")
); 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
#########################################################################################################
### Main
#########################################################################################################

# routine check if files exists 
if(file.exists(opt$code) == FALSE){
  stop(paste0("[ERROR]\t Pairing file does not exist, input file:\n", opt$code, "\n"))
} else if (all(file.exists(paste0(opt$plink, c(".bim", ".fam", ".bed")))) == FALSE){
  stop(paste0("[ERROR]\t At least one of the plink files (.bim .fam .bed) does not exist:\n", opt$plink, "\n"))
}else if (file.exists(opt$info) == FALSE){
  stop(paste0("[ERROR]\t Info/Phenotype files does not exist, input file:\n", opt$info, "\n"))
} else if(file.exists(opt$king) == FALSE){
  stop(paste0("[ERROR]\t Path to KING does not exist, input file:\n", opt$king, "\n"))
} else{
  cat("[INFO]\t All files from arguments exists \n")
}

cat("[INFO]\t Reading input files")
pairing.table <- fread(opt$code, data.table=F)
info.table <- fread(opt$info, data.table=F)
fam.table <- fread(file = paste0(opt$plink,".fam"))

####chequeo de IDs para generalizar

# Matching sample IDs from ".fam" file and "opt$code"  file. 
fam.table$sust <- sapply(fam.table$V1, FUN=codextract)  
pairing.table$sust <- sapply(pairing.table[,2], FUN=codextract) 
fam.table$PSEUDOIDEXT <- as.character(pairing.table[,1]
                                      [match(fam.table$sust, pairing.table$sust)])
# merge pedigree file with the modified .fam file. 
colnames(info.table)<-c ("PSEUDOIDEXT","FAM_ID" ,"FATHER_PSEUDOID","MOTHER_PSEUDOID","GENDER1M2F","PARTNEREXT")
fam.table <- left_join(fam.table, info.table, by="PSEUDOIDEXT")

#create dir for outpue
dir.create(path = opt$workdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path = opt$out, recursive = TRUE, showWarnings = FALSE)

## If there are samples in the plink for which a "PSEUDOIDEXT" is not assigned. Then we add complete this info in the fam file
n.na.pseudoID <- sum(is.na(fam.table$PSEUDOIDEXT))
if(n.na.pseudoID>=1){
  cat("[WARNING] a total of", n.na.pseudoID, 
      "samples from the .fam file are not present in the pairing.table, fake pedigree info will be introduce to assess this/these sample/s")
  
  na.pseudoID.index <- which(is.na(fam.table$PSEUDOIDEXT))
  fam.table[na.pseudoID.index,c("FATHER_PSEUDOID", "MOTHER_PSEUDOID", "GENDER1M2F", "PARTNEREXT")] <- matrix(0, ncol= 4, nrow= n.na.pseudoID)
  fam.table[na.pseudoID.index,c("PSEUDOIDEXT", "FAM_ID")] <- cbind(fam.table[na.pseudoID.index,"V1"],fam.table[na.pseudoID.index,"V1"])
}

## If there are persons with a PSEUDOIDEXT but no pedrigree info.
n.na.pedigree <- sum(is.na(fam.table$FAM_ID))
if(n.na.pedigree>=1){
  cat("[WARNING] a total of", n.na.pedigree, 
      "samples from the .fam file do not have pedigree information present in the pairing file, dummy pedigree information will be created for them")
  
  na.pedigree.index <- which(is.na(fam.table$FAM_ID))
  
  fam.table[na.pedigree.index,"FAM_ID"] <- fam.table[na.pedigree.index,"PSEUDOIDEXT"]
  fam.table[na.pedigree.index,c("FATHER_PSEUDOID", "MOTHER_PSEUDOID", "GENDER1M2F", "PARTNEREXT")] <- matrix(0, ncol= 4, nrow= n.na.pedigree)
  
  sub_dummy_ped<-data.frame(cbind(fam.table[na.pedigree.index,c("PSEUDOIDEXT",
                                                          "FAM_ID",
                                                          "FATHER_PSEUDOID",
                                                          "MOTHER_PSEUDOID", 
                                                          "GENDER1M2F", 
                                                          "PARTNEREXT")]))
  
  info.table<-as.data.frame(rbind(info.table,sub_dummy_ped))
                                  
}

#cols.for.new.fam <- c(9,8,10:12,6) ### order -> "FAM_ID", "PSEUDOIDEXT", "FATHER_PSEUDOID", "MOTHER_PSEUDOID", "GENDER1M2F", "V6"
new.fam <- fam.table[,c("FAM_ID", "PSEUDOIDEXT", "FATHER_PSEUDOID", "MOTHER_PSEUDOID", "GENDER1M2F", "V6")]

## If there are any duplicated "PSEUDOIDEXT" we change the name to avoid confusions and later report them
n.duplicated.PSEUDOIDEXT <- sum(duplicated(new.fam$PSEUDOIDEXT))
if(n.duplicated.PSEUDOIDEXT >= 1){
  cat("[WARNING] a total of", n.duplicated.PSEUDOIDEXT, 
      "samples are duplicated, a *_dup* tag will be added to the sample ID")
  
  dup.table <- table(new.fam$PSEUDOIDEXT[which(duplicated(new.fam$PSEUDOIDEXT))])
  for(i.dup.id in names(dup.table)){
    i.dup.index <- which(new.fam$PSEUDOIDEXT == i.dup.id)
    
    n.dup.sufix.vector <- c("", paste0("_", 1:c(dup.table[i.dup.id])))
    new.fam$PSEUDOIDEXT[i.dup.index] <- paste0(unique(new.fam$PSEUDOIDEXT[i.dup.index]), n.dup.sufix.vector)
  }
}
## keep the pairing of the newly named duplicated samples in a different dataframe for later
pair_dup<-data.frame(cbind(new.fam$PSEUDOIDEXT,fam.table$V2,fam.table$PSEUDOIDEXT))
colnames(pair_dup)<-c("dup_PSEUDOIDEXT","Plate_pos","PSEUDOIDEXT")

#### create a new folder for outpur  -> copy input files for king -> change working directory to run king -> lounch king through system()
new.fam.file <- file.path(opt$workdir,"batchinfo.fam")
write.table(new.fam, file=new.fam.file, row.names = FALSE, col.names = FALSE, quote = FALSE )
setwd(dir = opt$workdir)

## if we have the complete data then we merge the whole pedigree information with it, else we just check family in the batch
if (opt$makeped==TRUE  ) {
  ###in case dummy file is not constructed already, build it from the pedigree file
  if (file.exists(paste0(opt$dummy,".map"))==F){ 
    make_dummy_files(pedfile=opt$info,plink_path=opt$plink,dummy_file=opt$dummy)
    
  }
  ###merge new.fam with the whole pedigree 
  plink.system.call <- paste0("ml plink;",
                              " plink ",
                              " --bed ", paste0(opt$plink,".bed"),
                              " --fam ", paste0(opt$workdir,"batchinfo.fam"),
                              " --bim ", paste0(opt$plink,".bim"),
                              " --merge ", paste0(opt$dummy,".ped")," ", paste0(opt$dummy,".map"),
                              " --merge-mode 4 ",
                              " --make-bed",
                              " --out fullped")
  
  system(plink.system.call)
  #king system call based on http://people.virginia.edu/~wc9c/KING/manual.html#INPUT
  king.system.call <- paste0(c(opt$king),
                             " -b ", "fullped.bed",
                             " --related ", 
                             " --degree 2",
                             " --prefix famCheck_genotypeQC")
  
  system(king.system.call)
  
} else{
  
  #king system call based on http://people.virginia.edu/~wc9c/KING/manual.html#INPUT
  king.system.call <- paste0(c(opt$king),
                             " -b ", paste0(opt$plink,".bed"),
                             " --fam ", paste0(opt$workdir,"batchinfo.fam"),
                             " --bim ", paste0(opt$plink,".bim"),
                             " -- related", 
                             " --degree 2",
                             " --prefix famCheck_genotypeQC")
  system(king.system.call)
  
  cat('[INFO] KING is done processing\n\n')
}


##########
#Read king output. 
#####

# .king0 file contains the results for all possible pairs for which there is no genetec relationship implicated in the pedigree information.
king_0.file <- file.path(opt$workdir,"famCheck_genotypeQC.kin0")
king_0 <- fread(king_0.file, data.table = FALSE)

# .king file contains the results for all pairs with same family ID.
king.file <- file.path(opt$workdir,"famCheck_genotypeQC.kin")
king <- fread(king.file, data.table = FALSE)

### correct the errors caused by twins or other known relationships (as there is no way to indicate twins or actual fourth degrees from original pedirgree files)
if (is.null(opt$relations)==FALSE){
  ls<-fread(opt$relations,data.table=F)
  ls<-ls$V1
  for (i in 1:length(ls)){
    king[king$FID==ls[i] & king$Error==1,"Error"]<-0
  }
}

##########
#Prepare data for plotting
#####

# Number of erros across family relationships. 
king$Error.tag <- cut(king$Error, breaks=c(-Inf, 0.49 ,0.51,1), 
                      labels= c("Ok", "Warning", "Error"), right = T)

## Evaluate all errors from king output. 
king.error <- king[which(king$Error.tag == "Error"),]

## king.error can be further filtered to only include FS, PO and UN
king.error.filtered <- king.error[which(king.error$InfType %in% c("FS", "PO", "UN")),]
UN.error.king <- king.error.filtered[king.error.filtered$InfType == "UN",]

#for(i.un.error in 1:nrow(UN.error.king)){
#  info.table[info.table$FAM_ID %in% UN.error.king$FID[i.un.error],]

  
#}

##########
#Plotting
#####

n.errors.families.bar <- ggplot(king, aes(x=Error.tag))+
  geom_bar(stat = "count")+
  geom_text(aes(label=..count..),stat='count',position=position_dodge(1), size=4,vjust=0)+
  xlab("")+
  ggtitle("Number of concordant and non concordant genetic relationships")+
  theme_classic()+
  theme(text = element_text(size=10, family = "Helvetica"))

king.error$InfType <- factor(king.error$InfType, levels= c("Dup/MZ" , "PO", "FS", "2nd", "3rd","4th","UN"))
n.errors.fam.relation <- ggplot(king.error, aes(x= InfType))+
  geom_bar(stat = "count")+
  geom_text(aes(label=..count..),stat='count',position=position_dodge(1), size=4, vjust=0)+
  xlab("Relationship inferred by genetics")+
  ggtitle("Types of family relationships corrected from pedigree")+
  theme_classic()+
  theme(text = element_text(size=10, family = "Helvetica"))

infered.expected.hex <- ggplot(king, aes(x=Kinship, y=Z0))+
  geom_hex()+
  scale_fill_viridis_c()+
  theme_classic()+
  theme(text = element_text(size=10, family = "Helvetica"))

king_0$InfType <- factor(king_0$InfType, levels= c("Dup/MZ" , "PO", "FS", "2nd", "3rd"))
non.family.inferences.plot <- ggplot(king_0, aes(x=InfType))+
  geom_bar(stat = "count")+
  geom_text(aes(label=..count.. ),stat='count',position=position_dodge(1.5), size=4,vjust=0)+
  xlab("Relationship inferred by genetics")+
  ggtitle("Number non annotated genetic relationships")+
  theme_classic()+
  theme(text = element_text(size=10, family = "Helvetica"))

plot.file <- file.path(opt$out, paste0("07.famCheck_plots.tiff"))
tiff(plot.file,  width = 3000, height = 3500, units = "px", res = 300, compression = "lzw")
grid.arrange(n.errors.families.bar, n.errors.fam.relation,
             infered.expected.hex, non.family.inferences.plot,
             nrow=2)
dev.off()

###create duplicate samples only report
duplicate_samples<-rbind(king_0[king_0$InfType=="Dup/MZ", c("ID1","ID2")],king[king$InfType=="Dup/MZ", c("ID1","ID2")])
if (nrow(duplicate_samples)!=0){
  duplicate_samples$ID1 <- as.character(pair_dup$Plate_pos
                                        [match(duplicate_samples$ID1, pair_dup$dup_PSEUDOIDEXT)])
  duplicate_samples$ID2 <- as.character(pair_dup$Plate_pos
                                        [match(duplicate_samples$ID2, pair_dup$dup_PSEUDOIDEXT)])
  write.table(duplicate_samples,paste0(opt$workdir,"equal.samples"),quote = F,row.names = F)
}

####

### if the questionary information is compete we will want to make the pedigrees check, otherwise, just looking
### duplicates should be enough

if (opt$makeped==TRUE) {
  #nombrar columnas
  colnames(pairing.table)[c(1,3)]<-c("PSEUDOIDEXT",  "Age")
  info.table$Age<-pairing.table$Age[match(info.table$PSEUDOIDEXT,pairing.table$PSEUDOIDEXT)]
  #info.table$Birth_year<-pairing.table$BIRTHYEAR[match(info.table$PSEUDOIDEXT,pairing.table$PSEUDOIDEXT)]
  info.table[which(is.na(info.table$Age)),"Age"]<-" "
  ############
  ###Filter list of families for plotting according to: 1. number of individuals with genetic information more than 2, 1, first degree
  names(new.fam)<-c("FAM_ID","PSEUDOIDEXT","FATHER_PSEUDOID","MOTHER_PSEUDOID","GENDER1M2F","V6")
  
  ####for kin0: list of families with more than 3 memebers and new FIRST GRADE relationships
  fstdeg <- king_0[king_0$InfType=="FS" | king_0$InfType=="PO"| king_0$InfType %like% "MZ", ]
  fam_list0<-c()
  for (i in 1:nrow(fstdeg)) {
    size1<-nrow(info.table[which(info.table$FAM_ID==fstdeg$FID1[i] |info.table$PSEUDOIDEXT==fstdeg$ID1[i]), ])
    size2<-nrow(info.table[which(info.table$FAM_ID==fstdeg$FID2[i] |info.table$PSEUDOIDEXT==fstdeg$ID2[i]), ])
    if (size1>2) {fam_list0<-c(fam_list0,fstdeg$FID1[i])} 
    else { if (size2>2) {fam_list0<-c(fam_list0,fstdeg$FID2[i])}}
  }
  fam_list0<-unique(fam_list0)
  
  ####for kin: list of families with more than 3 memebers and reported FIRST GRADE, and unrelated, relationships with ERROR
  fstdeg2 <- king[which(king$Error==1 & (king$InfType=="FS" | king$InfType=="PO"| king$InfType %like% "MZ"| king$InfType=="UN")), ]
  fam_list<-c()
  for (i in 1:nrow(fstdeg2)) {
    size1<-nrow(info.table[which(info.table$FAM_ID==fstdeg2$FID[i]),])
    if (size1>2 ) {fam_list<-c(fam_list,fstdeg2$FID[i])} 
  }
  fam_list<-unique(fam_list)  
    
  ####make padigrees with cranefoot
  ##list
  all_error.list<-lapply(fam_list,FUN= pedigree_crane,pedfile=info.table,king=king,king_0=king_0,unexpected=F,
                         Fs.grade=F,fam_batch=new.fam,crane.path=opt$crane)
  
  all_new_found.list<-lapply(fam_list0,FUN= pedigree_crane,pedfile=info.table,king=king,king_0=king_0,
                             unexpected=T,Fs.grade=F,fam_batch=new.fam,crane.path=opt$crane)
  
  write.table(unlist(all_error.list),paste0(opt$out,"Crane_all_error.FID.list"),quote = F,row.names = F,col.names = F)
  write.table(unlist(all_new_found.list),paste0(opt$out,"Crane_all_newfound.FID.list"),quote = F,row.names = F,col.names = F)
  #########################################################################################################
}

####
##done





