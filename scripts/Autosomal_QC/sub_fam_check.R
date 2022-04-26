###################################
### Family relationship check
### date: 09-04-2019
### version: 0.01
### authors: EL - RAG
###################################
### New
# 04-04-2022
# removed plate specific reports
# 16-04-2018
# Added file.exists checks 
# changed king input parameters to include -fam, -bed and -bim files to allow different names and not have to duplicate files. 
###################################

library(data.table)
library(tidyverse)
library(optparse)
library(gridExtra)
library(reshape)


#########################################################################################################
option_list = list(
  make_option(c("-p", "--plink"), type="character", default=NULL, 
              help="Path to plink files index, it assumes a bed, bim and fam file with the same file name", metavar="character"),
  
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
  
  make_option(c("-w", "--workdir"), type="character", 
              default="opt$out", help="processing directory", metavar="character")
); 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)



#########################################################################################################
### Main
#########################################################################################################

# routine check if files exists 
if (all(file.exists(paste0(opt$plink, c(".bim", ".fam", ".bed")))) == FALSE){
  stop(paste0("[ERROR]\t At least one of the plink files (.bim .fam .bed) does not exist:\n", opt$plink, "\n"))
}else if (file.exists(opt$info) == FALSE){
  stop(paste0("[ERROR]\t Info/Phenotype files does not exist, input file:\n", opt$info, "\n"))
} else if(file.exists(opt$king) == FALSE){
  stop(paste0("[ERROR]\t Path to KING does not exist, input file:\n", opt$king, "\n"))
} else{
  cat("[INFO]\t All files from arguments exists \n")
}

cat("[INFO]\t Reading input files")
info.table <- fread(opt$info, data.table=F)
fam.table <- fread(file = paste0(opt$plink,".fam"))

# Matching sample IDs from ".fam" file info  file. Remove duplicate indicator
fam.table$sust <- gsub(fam.table$V2,pattern="_[0-9]",replacement = "")
info.table$sust<-gsub(info.table$V2,pattern="_[0-9]",replacement = "")
# merge infor table with .fam file. 
fam.table <- left_join(fam.table, info.table, by="sust")

##organize file
fam.table<-fam.table[,c(8,7,10:13,1)]
names(fam.table)<-c("FAM_ID","IID","FATHER_PSEUDOID","MOTHER_PSEUDOID","GENDER1M2F","PARTNEREXT","IIDx")
fam.table<-data.frame(fam.table)
fam.table$PARTNEREXT<-sapply(fam.table$PARTNEREXT,function(x){ifelse(x==-9,0,x)})

## If there are samples in the plink for which a "V2" is not assigned. Then we add complete this info in the fam file
n.na.pseudoID <- sum(is.na(fam.table$FAM_ID))
if(n.na.pseudoID>=1){
  cat("[WARNING] a total of", n.na.pseudoID, 
      "samples from the .fam file are not present in the pairing.table, fake pedigree info will be introduce to assess this/these sample/s")
  
  na.pseudoID.index <- which(is.na(fam.table$FAM_ID))
  fam.table[na.pseudoID.index,c("FATHER_PSEUDOID", "MOTHER_PSEUDOID", "GENDER1M2F", "PARTNEREXT")] <- data.frame(matrix(0, ncol= 4, nrow= n.na.pseudoID))
  fam.table[na.pseudoID.index,c("IID", "FAM_ID")] <- cbind(fam.table[na.pseudoID.index,"IIDx"],fam.table[na.pseudoID.index,"IIDx"])
}

## If there are persons with an IID but no pedrigree info.
n.na.pedigree <- sum(is.na(fam.table$FAM_ID))
if(n.na.pedigree>=1){
  cat("[WARNING] a total of", n.na.pedigree, 
      "samples from the .fam file do not have pedigree information present in the pairing file")
  
  na.pedigree.index <- which(is.na(fam.table$FAM_ID))
  
  fam.table[na.pedigree.index,"FAM_ID"] <- fam.table[na.pedigree.index,"V2"]
  fam.table[na.pedigree.index,c("FATHER_PSEUDOID", "MOTHER_PSEUDOID", "GENDER1M2F", "PARTNEREXT")] <- matrix(0, ncol= 4, nrow= n.na.pedigree)
}

#cols.for.new.fam <- c(9,8,10:12,6) ### order -> "FAM_ID", "V2", "FATHER_PSEUDOID", "MOTHER_PSEUDOID", "GENDER1M2F", "V6"
new.fam <- cbind(fam.table[,c(1:5)],"V6"=rep(-9,nrow(fam.table)))

#### create a new folder for outpur  -> copy input files for king -> change working directory to run king -> lounch king through system()
new.fam.file <- file.path(opt$workdir,"batchinfo.fam")
write.table(new.fam, file=new.fam.file, row.names = FALSE, col.names = FALSE, quote = FALSE )
setwd(dir = opt$workdir)

## if we have the complete data the we merge the whole questionary information with it, else we just check family in the batch

  #king system call based on http://people.virginia.edu/~wc9c/KING/manual.html#INPUT
  king.system.call <- paste0(c(opt$king),
                             " -b ", paste0(opt$plink,".bed"),
                             " --fam ", paste0(opt$workdir,"/batchinfo.fam"),
                             " --bim ", paste0(opt$plink,".bim"),
                             " --related", 
                             " --degree 2",
                             " --prefix famCheck_genotypeQC")
  system(king.system.call)
  

##########
#Read king output. 
#####

# .king0 file contains the results for all possible pairs for which there is no genetec relationship implicated in the pedigree information.
king_0.file <- file.path(opt$workdir,"famCheck_genotypeQC.kin0")
king_0 <- fread(king_0.file, data.table = FALSE)

# .king file contains the results for all pairs with same family ID.
king.file <- file.path(opt$workdir,"famCheck_genotypeQC.kin")
king <- fread(king.file, data.table = FALSE)

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

#table(as.data.frame(info.table[info.table$FAM_ID %in% unique(UN.error.king$FID),])$FAM_ID)

for(i.un.error in 1:nrow(UN.error.king)){
  info.table[info.table$FAM_ID %in% UN.error.king$FID[i.un.error],]
  
}

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
  write.table(duplicate_samples,paste0(opt$workdir,"/equal.samples"),quote = F,row.names = F)

####

### if the questionary information is compete we will want to make the pedigrees check, otherwise, just looking
### duplicates should be enough
 if (opt$makeped==TRUE) {
  info.table$Age<-" "
  if(file.exists(opt$pairing) == TRUE){
  pairing.table <- fread(opt$pairing, data.table=F)
  info.table$Age<-pairing.table$age[match(info.table$IID,pairing.table$IID)]
  info.table$Birth_year<-pairing.table$BIRTHYEAR[match(info.table$IID,pairing.table$IID)]
  info.table[which(is.na(info.table$Age)),"Age"]<-" "
  }
 
  ############
  ###prepare list of families for plotting
  names(new.fam)<-c("FAM_ID","IID","FATHER_PSEUDOID","MOTHER_PSEUDOID","GENDER1M2F","V6")

  ####for kin0: list of families with more than 3 memebers and new FIRST GRADE relationships
  fstdeg <- king_0[king_0$InfType=="FS" | king_0$InfType=="PO", ]
  fam_list0<-c()
  for (i in 1:nrow(fstdeg)) {
    size1<-nrow(new.fam[which(new.fam$FAM_ID==fstdeg$FID1[i] |new.fam$IID==fstdeg$ID1[i]), ])
    size2<-nrow(new.fam[which(new.fam$FAM_ID==fstdeg$FID2[i] |new.fam$IID==fstdeg$ID2[i]), ])
    if (size1>2) {fam_list0<-c(fam_list0,fstdeg$FID1[i])} 
    else { if (size2>2) {fam_list0<-c(fam_list0,fstdeg$FID2[i])}}
  }
  fam_list0<-unique(fam_list0)
  
  ####for kin: list of families with more than 3 memebers and reported FIRST GRADE relationships with ERROR
  fstdeg2 <- king[which(king$Error==1 & (king$InfType=="FS" | king$InfType=="PO")), ]
  
  fam_list<-c()
  for (i in 1:nrow(fstdeg2)) {
    size1<-nrow(new.fam[which(new.fam$FAM_ID==fstdeg2$FID[i]),])
    if (size1>2 ) {fam_list<-c(fam_list,fstdeg2$FID[i])} 
  }
  fam_list<-unique(fam_list)
  
  
  #########################################################################################################
  ### Pedigree function with Cranefoot
  
  pedigree_crane<-function(ls,pedfile,king,king_0,unexpected=TRUE,Fs.grade=F,fam_batch,crane.path){
    
    ########requirements
    #ls: vector containing the families to be viewed in a pedigree
    #pedfile: .famfile with column names= "FAM_ID","V2","FATHER_PSEUDOID","MOTHER_PSEUDOID","GENDER1M2F","PARTNEREXT"
    #king:kin file with expected genetic relationships
    #king_0; kin0 file with unexpected genetic relationships
    #unexpected: TRUE to visualize unexpected (kin0) genetic relationships, FALSE to evaluate pedigree-known (.kin) relationships
    #Fs.grade: TRUE to evaluate only first grade relationships, FAlSE to see all relationshipsin the families listed in ls
    #fam_batch: the fam file from the batch
    #crane.path: path to the cranefoot executable
    
    if(Fs.grade==T){rel_vector=c("PO","FS")} else {rel_vector=king_0$InfType}
    
    ####generate ethe single family (or 2 families) file
    if(unexpected==T) {
      ###processing from kin0 information
      IDS<-rbind(king_0[which((king_0$FID1==ls|king_0$FID2==ls) & 
                                king_0$InfType %in% rel_vector),"ID2"], 
                 king_0[which((king_0$FID1==ls|king_0$FID2==ls) & 
                                king_0$InfType %in% rel_vector),"ID1"] )
      family.ped<-pedfile[which(pedfile$FAM_ID==ls|pedfile$V2 %in% IDS),c(1,2,3,4,5,7) ] 
    } else {
      ###processing from kin information
      family.ped<-pedfile[pedfile$FAM_ID==ls ,c(1,2,3,4,5,7) ]
    }
    names(family.ped)<-c("IID","FAM_ID","FATHER_PSEUDOID","MOTHER_PSEUDOID","GENDER1M2F","AGE")
    
    
    ##retrieve parents to make them individuals
    extra.IDs <- unique(c(family.ped[,"FATHER_PSEUDOID"],family.ped[,"MOTHER_PSEUDOID"]))
    `%!in%` = Negate(`%in%`)## create negation of %in% function
    extra.IDs <- extra.IDs[ which(extra.IDs %!in% family.ped$V2)]   # remove extra IDs (parents) which are already in the final.fam$V2
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
    final.fam$Genetic_info<-ifelse(final.fam$IID %in% fam_batch$V2,"550077","999999")## genetic info colored "purple"
    
    #### include relations from kin
    
    if(unexpected==T) {
      ###for kin 0
      relations<-king_0[which(king_0$InfType!="UN" & 
                                (king_0$FID2 %in% unique(final.fam$FAM_ID) & king_0$FID1 %in% unique(final.fam$FAM_ID))
                              & king_0$InfType %in% rel_vector), ]
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
    if (gen.fam=="incomplete family information") {return(paste0("Family ",ls," incomplete"))} 
    ####incomplete will be reported when there is no pedigree information  for genetic samples
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
      
      ###crate directory to save the result plots
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
            "MotherVariable","GenderVariable","ColorVariable","TextVariable","TextVariable","ArrowVariable",rep("TextVariable",n.text.col)),
          
          c(paste0(opt$workdir,"family.ped"),paste0(crane.dir,"07.",event ,ls,".ped"),
            "IID","FATHER_PSEUDOID","MOTHER_PSEUDOID","GENDER1M2F","Genetic_info","IID","AGE",
            "Familial_errors", paste0("Genetic_relationship_",seq(1:n.text.col)) ),
          
          c("","","","","",paste0(opt$workdir,"phenotype_1"),paste0(opt$workdir,"phenotype_1"),paste0(opt$workdir,"phenotype_1"),
            paste0(opt$workdir,"phenotype_1"),paste0(opt$workdir,"fam.error"), rep(paste0(opt$workdir,"fam.error"),n.text.col) )
          
        )
      )
      write.table(config.file,paste0(opt$workdir,"CFG"),quote = F,row.names = F,col.names = F,sep = "\t")
      
      ## create system call for cranefoot
      crane.system.call <- paste0(c(opt$crane)," ",paste0(opt$workdir,"CFG"))
      system(crane.system.call)
      return(paste0("Family ",ls," done"))
    }
  }
  
  ###########################################################################################
  
  ####make padigrees with cranefoot
  ##list
  all_error.list<-lapply(fam_list,FUN= pedigree_crane,pedfile=info.table,king=king,king_0=king_0,unexpected=F,
                         Fs.grade=F,fam_batch=new.fam,crane.path=opt$crane)
  
  all_new_found.list<-lapply(fam_list0,FUN= pedigree_crane,pedfile=info.table,king=king,king_0=king_0,
                             unexpected=T,Fs.grade=F,fam_batch=new.fam,crane.path=opt$crane)
  
  write.table(unlist(all_error.list),paste0(opt$out,"Crane_all_error.FID.list"),quote = F,row.names = F,col.names = F)
  write.table(unlist(all_new_found.list),paste0(opt$out,"Crane_all_newfound.FID.list"),quote = F,row.names = F,col.names = F)
}
  
