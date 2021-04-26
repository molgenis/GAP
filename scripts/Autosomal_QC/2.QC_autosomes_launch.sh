#!/bin/bash
#SBATCH --time=15:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=150G

module load plink
module load RPlus


##inpit and putput variables
InputDir="/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC"
GeneralQCDir="/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC"


### Reference files
MAFref="/groups/umcg-wijmenga/tmp04/umcg-raguirre/pln_ugli/af.ref.data.txt" ## reference for external marker concordance
samplesheet="/groups/umcg-ugli/tmp04/projects/UGLI${UGLIbatch}/run01/results/UGLI${UGLIbatch}.csv" ## Identificator list for samples

ref1000G="/apps/data/1000G/phase3/1000G_all" ## reference for PCA from 1000genomes
gonlref="/groups/umcg-aad/tmp04/umcg-elopera/thgref/gonl_SNV_INDELs_ab" ## reference for PCA from GoNL
commonSNPs="/groups/umcg-aad/tmp04/umcg-elopera/thgref/gonl_1KGenomes_GSA_common_SNPs" ## PCA reference of common snps from GoNL and 1000G, prunned according to doi:10.1038/ejhg.2014.19

pairing_ID="/groups/umcg-aad/tmp04/umcg-elopera/ugli_blood_gsa/pairing.dat" ## pairing file for the IDs of family and samples
pedigree_ref="/groups/umcg-aad/tmp04/umcg-elopera/ugli_blood_gsa/LifeLines_Corrected_family_info.dat" ## pedigree information
dummy_pedigree="/groups/umcg-aad/tmp04/umcg-elopera/ugli_blood_gsa/corrected_v2_dummy_allped" ## dummy plink file with information of all the pedigree
king_tool="/groups/umcg-aad/tmp04/umcg-elopera/tools/KING/king" ## exact location of the KING executable
cranefoot_tool="/groups/umcg-aad/tmp04/umcg-elopera/tools/Cranefoot/example/cranefoot" ## exact location of the Cranefoot executable

all_pairing_file="/groups/umcg-aad/tmp04/umcg-elopera/concordanceCheck/ugli.final.pairing.concordanceCheck.txt" ## pairing file for external samples
merged_conc1="/groups/umcg-aad/tmp04/umcg-elopera/concordanceCheck/gwas_og_plink_chr/filtered_merged_LL_GWAS" ## merged file 1, for the previous gwas study
merged_conc2="/groups/umcg-aad/tmp04/umcg-elopera/concordanceCheck/gonl_og_plink_chr/filtered_merged_goNL_OG" ## merged file 2, for the goNL study

###create working directories
mkdir -p "${GeneralQCDir}/"
mkdir -p "${GeneralQCDir}/0_pre"
mkdir -p "${GeneralQCDir}/1_CR80/"
mkdir -p "${GeneralQCDir}/2_CR_high/"
mkdir -p "${GeneralQCDir}/3_MAF_HWE/"
mkdir -p "${GeneralQCDir}/plots/"
mkdir -p "${GeneralQCDir}/4_Het/"
mkdir -p "${GeneralQCDir}/4_Het/proc/"
mkdir -p "${GeneralQCDir}/5_Relatedness/"
mkdir -p "${GeneralQCDir}/5_Relatedness/proc/"
mkdir -p "${GeneralQCDir}/6_PCA"
mkdir -p "${GeneralQCDir}/6_PCA/proc2"
mkdir -p "${GeneralQCDir}/6_PCA/proc"
mkdir -p "${GeneralQCDir}/X_QC"
mkdir -p "${GeneralQCDir}/Y_QC"
mkdir -p "${GeneralQCDir}/MT_QC"
 
for chr in {1..22} "XY" "X" "MT" "Y"
do
sed -i 's|'_[0-9]_'|'_'|g' ${InputDir}/chr_${chr}.fam ###turn repeated samples into the original sample name DNA###_[A-Z]#
done

mv  ${InputDir}/chr_Y.* ${GeneralQCDir}/Y_QC/
mv  ${InputDir}/chr_MT.* ${GeneralQCDir}/MT_QC/
mv  ${InputDir}/chr_X.* ${GeneralQCDir}/X_QC/
mv  ${InputDir}/chr_* ${GeneralQCDir}/0_pre/

 for chr in {1..22} "XY"
 do

  ### create plink files and call_rate stats for individuals and SNPs
   plink --bfile ${GeneralQCDir}/0_pre/chr_${chr} \
         --missing \
         --out ${GeneralQCDir}/0_pre/chr_${chr}
 
  ##create list of samples to exclude on the criteria callrate<=80
  sed -i 's|'_[0-9]_'|'_'|g' ${GeneralQCDir}/0_pre/chr_${chr}.fam ###turn repeated samples into the original sample name DNA###_[A-Z]#
   awk '$6>0.20 {print $1, $2}' ${GeneralQCDir}/0_pre/chr_${chr}.imiss > ${GeneralQCDir}/1_CR80/chr_${chr}.extr80_sam.temp
done
 
## create files with duplicate SNPs to exlude 
Rscript /groups/umcg-aad/tmp04/umcg-elopera/position_duplicates.R -i ${GeneralQCDir}/0_pre

##creates list of individuals and dup snps excluded for all the autosomes
cat ${GeneralQCDir}/1_CR80/chr_*.extr80_sam.temp |sort -u > ${GeneralQCDir}/1_CR80/extr80.samples
rm ${GeneralQCDir}/1_CR80/chr_*.extr80_sam.temp
cat ${GeneralQCDir}/0_pre/chr_*.excl.duplicates > ${GeneralQCDir}/0_pre/extr.dups
rm ${GeneralQCDir}/0_pre/*.excl.duplicates
rm ${GeneralQCDir}/0_pre/*nosex*

for chr in {1..22} "XY"
 do
 # exclude individuals with callrate<=80 (creates excluded individuals_file) creates data set with individual_ callrate>80
   plink --bfile ${GeneralQCDir}/0_pre/chr_${chr}  \
         --make-bed \
         --remove ${GeneralQCDir}/1_CR80/extr80.samples \
         --exclude ${GeneralQCDir}/0_pre/extr.dups \
         --allow-no-sex \
         --out ${GeneralQCDir}/1_CR80/chr_${chr}
 
 #calculates callrate stats from the previously filtered datafile
   plink  --bfile ${GeneralQCDir}/1_CR80/chr_${chr} \
          --missing \
          --allow-no-sex \
          --out ${GeneralQCDir}/1_CR80/chr_${chr}

 awk '$5>0.20 {print $2}' ${GeneralQCDir}/1_CR80/chr_${chr}.lmiss > ${GeneralQCDir}/1_CR80/chr_${chr}.extr80_var.temp
 done
 
##file of snps to exclude
cat ${GeneralQCDir}/1_CR80/chr_*.extr80_var.temp > ${GeneralQCDir}/1_CR80/extr80.vars
rm ${GeneralQCDir}/1_CR80/*.temp ##remove chrosome files for excluded samples and markers
rm ${GeneralQCDir}/1_CR80/*nosex* 
 for chr in {1..22} "XY"
 do
   # exclude  SNPs with CR<80
   plink --bfile ${GeneralQCDir}/1_CR80/chr_${chr}  \
         --make-bed \
         --exclude ${GeneralQCDir}/1_CR80/extr80.vars\
         --out ${GeneralQCDir}/1_CR80/chr_${chr}.2

   #calculates callrate stats from the previously filtered datafile
   plink  --bfile ${GeneralQCDir}/1_CR80/chr_${chr}.2 \
          --missing \
          --out ${GeneralQCDir}/1_CR80/chr_${chr}.2

   ##create list of SNPs snd samples to exclude on the criteria callrate<=high
   awk '$6>0.01 {print $1, $2}' ${GeneralQCDir}/1_CR80/chr_${chr}.2.imiss > ${GeneralQCDir}/2_CR_high/chr_${chr}.extrhigh_sam.temp
   ##information for the heterozygosity analysis
   awk '$6<0.01 {print $1, $2,$6}' ${GeneralQCDir}/1_CR80/chr_${chr}.imiss > ${GeneralQCDir}/2_CR_high/chr_${chr}.incl_CR_sam
 done
cat ${GeneralQCDir}/2_CR_high/chr_*.extrhigh_sam.temp|sort -u > ${GeneralQCDir}/2_CR_high/extrhigh.samples


 for chr in {1..22} "XY"
 do
   ## exclude individuals with callrate<=high (creates excluded individuals_file) creates data set with  individual_ callrate>99
   plink --bfile ${GeneralQCDir}/1_CR80/chr_${chr}  \
         --make-bed \
         --remove ${GeneralQCDir}/2_CR_high/extrhigh.samples \
         --out ${GeneralQCDir}/2_CR_high/chr_${chr}
 
  #calculates callrate stats from the previously filtered datafile
   plink  --bfile ${GeneralQCDir}/2_CR_high/chr_${chr} \
          --missing \
          --out ${GeneralQCDir}/2_CR_high/chr_${chr}
          
 awk '$5>0.01 {print $2}' ${GeneralQCDir}/1_CR80/chr_${chr}.lmiss > ${GeneralQCDir}/2_CR_high/chr_${chr}.extrhigh_var.temp
 done
cat ${GeneralQCDir}/2_CR_high/chr_*.extrhigh_var.temp > ${GeneralQCDir}/2_CR_high/extrhigh.vars
rm ${GeneralQCDir}/2_CR_high/*.temp ##remove chrosome files for excluded samples and markers

 for chr in {1..22} "XY"
 do
   ## exclude SNPs callrate>99
   plink --bfile ${GeneralQCDir}/2_CR_high/chr_${chr}  \
         --make-bed \
         --exclude ${GeneralQCDir}/2_CR_high/extrhigh.vars \
         --out ${GeneralQCDir}/2_CR_high/chr_${chr}.2

 done
 
##creates merged files of included individuals and SNPs to be used of further analysis
 
cat ${GeneralQCDir}/0_pre/chr_*.fam|sort -u|awk '{print$2}' > ${GeneralQCDir}/0_pre/full.ind
cat ${GeneralQCDir}/0_pre/chr_*.bim|awk '{print$2}' > ${GeneralQCDir}/0_pre/untouched.snps
sort ${GeneralQCDir}/0_pre/extr.dups ${GeneralQCDir}/0_pre/untouched.snps|uniq -u >${GeneralQCDir}/0_pre/full.snps 

cat ${GeneralQCDir}/1_CR80/chr_*.2.fam|sort -u|awk '{print$2}' > ${GeneralQCDir}/1_CR80/incl80.samples
cat ${GeneralQCDir}/1_CR80/chr_*.2.bim|awk '{print$2}' > ${GeneralQCDir}/1_CR80/incl80.vars
 
cat ${GeneralQCDir}/2_CR_high/chr_*.2.fam|sort -u|awk '{print$2}' > ${GeneralQCDir}/2_CR_high/inclhigh.samples
cat ${GeneralQCDir}/2_CR_high/chr_*.2.bim|awk '{print$2}' > ${GeneralQCDir}/2_CR_high/inclhigh.vars
 
 
##################################################################################################
################-------------MAF and HWE filtering--------########################################

###Eliminate outlier markers from  H-WE 
for chr in {1..22} "XY"
  do
 
    ## calculate MAF and HWE 
   plink --bfile ${GeneralQCDir}/2_CR_high/chr_${chr}.2 \
         --freq \
         --hardy \
         --out ${GeneralQCDir}/3_MAF_HWE/chr_${chr}
 
  awk '$9<0.000001 {print $2}' ${GeneralQCDir}/3_MAF_HWE/chr_${chr}.hwe|sort > ${GeneralQCDir}/3_MAF_HWE/highhw_${chr}.temp
  awk '$5==0 {print $2}' ${GeneralQCDir}/3_MAF_HWE/chr_${chr}.frq|sort > ${GeneralQCDir}/3_MAF_HWE/zeroMAF_${chr}.temp
  cat ${GeneralQCDir}/3_MAF_HWE/highhw_${chr}.temp ${GeneralQCDir}/3_MAF_HWE/zeroMAF_${chr}.temp  > ${GeneralQCDir}/3_MAF_HWE/extr_${chr}hw
  rm ${GeneralQCDir}/3_MAF_HWE/*.temp
  done
cat ${GeneralQCDir}/3_MAF_HWE/extr_*hw> ${GeneralQCDir}/3_MAF_HWE/excl_HW.snps
rm  ${GeneralQCDir}/3_MAF_HWE/extr_*hw
for chr in {1..22} "XY"
 do
   ## extract markers with MAF>0.01 and WHE< 1x 10 exp(-6)
   plink --bfile ${GeneralQCDir}/2_CR_high/chr_${chr}.2 \
         --make-bed \
         --exclude ${GeneralQCDir}/3_MAF_HWE/excl_HW.snps\
         --out ${GeneralQCDir}/3_MAF_HWE/chr_${chr}
 done
cat ${GeneralQCDir}/3_MAF_HWE/chr_*.bim|awk '{print$2}' > ${GeneralQCDir}/3_MAF_HWE/incl_HW.snps

##################################################################################################
################-------------Heterozygosity for samples--------###################################

#create list to merge files
find ${GeneralQCDir}/3_MAF_HWE/ -name "*.bim" > ${GeneralQCDir}/4_Het/proc/allchr.list;
sed -i 's/.bim//g' ${GeneralQCDir}/4_Het/proc/allchr.list;
#merge all chromosomes into one single genotype ...set of files (.fam, .bim, .bed)
plink --merge-list ${GeneralQCDir}/4_Het/proc/allchr.list \
       --out ${GeneralQCDir}/4_Het/proc/full_autosomal_het.temp

##first, separate a list of SNPs to exclude by LD (--indep [SNPwindow] [shift] [LD threshold in 1/(1-r2)])
plink --bfile ${GeneralQCDir}/4_Het/proc/full_autosomal_het.temp \
      --indep 50 5 2.5 --out ${GeneralQCDir}/4_Het/proc/full_extract.temp ;
#then exclude this list from the working files
plink --bfile ${GeneralQCDir}/4_Het/proc/full_autosomal_het.temp \
      --extract ${GeneralQCDir}/4_Het/proc/full_extract.temp.prune.in \
      --make-bed --out ${GeneralQCDir}/4_Het/proc/pruned_autosomal_het.temp;
 
#perform heterozigocity
plink --het \
       --bfile ${GeneralQCDir}/4_Het/proc/pruned_autosomal_het.temp \
       --homozyg \
       --out ${GeneralQCDir}/4_Het/autosomal
#erase temp files 
rm ${GeneralQCDir}/4_Het/proc/*temp*
#retrieve Call rate information from samples
cat ${GeneralQCDir}/2_CR_high/chr_*.incl_CR_sam> ${GeneralQCDir}/4_Het/CR.samples
rm ${GeneralQCDir}/2_CR_high/chr_*.incl_CR_sam
## create file with samples to exclude (het>4sd) and heterozygosity density plot
 Rscript /groups/umcg-aad/tmp04/umcg-elopera/Het_autosomeQC.R -i ${GeneralQCDir}/4_Het \
   -o ${GeneralQCDir}/plots
 
## Create QCed files corrected by heterozygosity 
for chr in {1..22} "XY"
 do
   plink --bfile ${GeneralQCDir}/3_MAF_HWE/chr_"${chr}" \
         --make-bed \
         --remove ${GeneralQCDir}/4_Het/Excluded.het \
         --out ${GeneralQCDir}/4_Het/chr_${chr}
 done
 
 
##################################################################################################
###############---------Relatedness and identity by descent (IBD)------###########################

#create list to merge files
find ${GeneralQCDir}/4_Het/ -name "*.bim" > ${GeneralQCDir}/5_Relatedness/proc/allchr.list;
sed -i 's/.bim//g' ${GeneralQCDir}/5_Relatedness/proc/allchr.list;
#merge all chromosomes into one single genotype ...set of files (.fam, .bim, .bed)
plink --merge-list ${GeneralQCDir}/5_Relatedness/proc/allchr.list \
      --missing \
      --out ${GeneralQCDir}/5_Relatedness/proc/full_autosomal_rel.temp
#separate the HLA SNPs
awk '{if ($1 == 6 && $4 >= 28477797 && $4 <= 35000000) print $2}' \
${GeneralQCDir}/5_Relatedness/proc/full_autosomal_rel.temp.bim > ${GeneralQCDir}/5_Relatedness/proc/HLAexclude.txt

## apply filters to reduce SNP number
 plink --bfile ${GeneralQCDir}/5_Relatedness/proc/full_autosomal_rel.temp \
       --exclude ${GeneralQCDir}/5_Relatedness/proc/HLAexclude.txt \
       --make-bed \
       --out ${GeneralQCDir}/5_Relatedness/proc/full_data


## exclude intentionally duplicated samples
Rscript /groups/umcg-aad/tmp04/umcg-elopera/sample_duplicates.R -w ${GeneralQCDir}/5_Relatedness/proc/ \
        -r $all_pairing_file
        
        
## apply filters to reduce SNP number
 plink --bfile ${GeneralQCDir}/5_Relatedness/proc/full_data \
       --remove ${GeneralQCDir}/5_Relatedness/proc/intended.duplicates \
       --make-bed \
       --out ${GeneralQCDir}/5_Relatedness/proc/full_data.no.dup


### genetic family concordance
Rscript /groups/umcg-aad/tmp04/umcg-elopera/fam_check3.R -p ${GeneralQCDir}/5_Relatedness/proc/full_data.no.dup \
                              -c ${pairing_ID} \
                              -i ${pedigree_ref} \
                              -k ${king_tool} \
                              -d ${dummy_pedigree} \
                              -C ${cranefoot_tool} \
                              -w ${GeneralQCDir}/5_Relatedness/proc2/ \
                              -M TRUE \
                              -o ${GeneralQCDir}/plots/


#remove temp files
rm ${GeneralQCDir}/5_Relatedness/proc/*temp*

#######################################################################################################
############################--------X chromosome QC and sex check--------##############################

#make working directories for X chromosome
mkdir -p "${GeneralQCDir}/X_QC/0_pre"
mkdir -p "${GeneralQCDir}/X_QC/1_CR80"
mkdir -p "${GeneralQCDir}/X_QC/2_CR_high"
mkdir -p "${GeneralQCDir}/X_QC/3_MAF_HWE"

###############--------Call rate and duplicate SNP QC for X-----#########

cat ${GeneralQCDir}/4_Het/Excluded.het  ${GeneralQCDir}/2_CR_high/extrhigh.samples ${GeneralQCDir}/1_CR80/extr80.samples > \
${GeneralQCDir}/X_QC/0_pre/excludebeforeX.samples

### create plink files and call_rate stats for individuals and SNPs
plink --bfile ${GeneralQCDir}/X_QC/chr_X \
     --make-bed  \
     --remove ${GeneralQCDir}/X_QC/0_pre/excludebeforeX.samples \
     --missing \
     --out ${GeneralQCDir}/X_QC/0_pre/chr_X

## generate list of duplicated SNPs (selecting the one with best call rate).
Rscript /groups/umcg-aad/tmp04/umcg-elopera/position_duplicates.R -i ${GeneralQCDir}/X_QC/0_pre 

##creates list of individuals and dup snps excluded for the X chromosome
cat ${GeneralQCDir}/X_QC/0_pre/chr_X.excl.duplicates > ${GeneralQCDir}/X_QC/0_pre/extr.dups

plink --bfile ${GeneralQCDir}/X_QC/0_pre/chr_X  \
     --make-bed \
     --exclude ${GeneralQCDir}/X_QC/0_pre/extr.dups \
     --out ${GeneralQCDir}/X_QC/1_CR80/chr_X

awk '$5>0.20 {print $2}' ${GeneralQCDir}/X_QC/0_pre/chr_X.lmiss > ${GeneralQCDir}/X_QC/1_CR80/chr_X.extr80_var
cat ${GeneralQCDir}/X_QC/1_CR80/chr_*.extr80_var > ${GeneralQCDir}/X_QC/1_CR80/extr80.vars

# exclude SNPs with a call rate < 80%
plink --bfile ${GeneralQCDir}/X_QC/1_CR80/chr_X  \
     --make-bed \
     --exclude ${GeneralQCDir}/X_QC/1_CR80/extr80.vars\
     --out ${GeneralQCDir}/X_QC/1_CR80/chr_X.2

#calculates callrate stats from the previously filtered datafile (removed SNPS with call rate < 80%)
plink  --bfile ${GeneralQCDir}/X_QC/1_CR80/chr_X.2 \
      --missing \
      --out ${GeneralQCDir}/X_QC/1_CR80/chr_X

#### NO SAMPLES are to be removed based on X chromosome call rate. 

awk '$5>0.01 {print $2}' ${GeneralQCDir}/X_QC/1_CR80/chr_X.lmiss > ${GeneralQCDir}/X_QC/2_CR_high/extrhigh.vars

## exclude SNPs callrate>99
plink --bfile ${GeneralQCDir}/X_QC/1_CR80/chr_X.2 \
     --make-bed \
     --exclude ${GeneralQCDir}/X_QC/2_CR_high/extrhigh.vars \
     --out ${GeneralQCDir}/X_QC/2_CR_high/chr_X


#### Files with variants throughout the CR filtering steps. 
cat ${GeneralQCDir}/X_QC/0_pre/chr_*.bim|awk '{print$2}' > ${GeneralQCDir}/X_QC/0_pre/untouched.snps
sort ${GeneralQCDir}/X_QC/0_pre/extr.dups ${GeneralQCDir}/X_QC/0_pre/untouched.snps|uniq -u >${GeneralQCDir}/X_QC/0_pre/full.snps 

cat ${GeneralQCDir}/X_QC/1_CR80/chr_*.2.bim|awk '{print$2}' > ${GeneralQCDir}/X_QC/1_CR80/incl80.vars
cat ${GeneralQCDir}/X_QC/2_CR_high/chr_*.bim|awk '{print$2}' > ${GeneralQCDir}/X_QC/2_CR_high/inclhigh.vars


######################### Impute sex and sexcheck ###################################
plink --bfile  ${GeneralQCDir}/X_QC/2_CR_high/chr_X  --impute-sex --make-bed --out ${GeneralQCDir}/X_QC/0_pre/imputed

###call the Rscript for plotting sex concordance
Rscript /groups/umcg-aad/tmp04/umcg-elopera/sexCheck_genotypeQC.R -i ${GeneralQCDir}/X_QC/0_pre/imputed.sexcheck \
                              -p ${pairing_ID}  \
                              -d ${GeneralQCDir}/5_Relatedness/proc2/equal.samples \
                              -o ${GeneralQCDir}/plots


grep -E 'Non concordant|Failed'  ${GeneralQCDir}/plots/sex_check/all.samples.concordance*.txt| awk '{{print $4, $4}}' > ${GeneralQCDir}/X_QC/sex.flagged

####################Filter out SNPs based on HW and MAF#############################
## calculate MAF and HWE
#plink --bfile ${GeneralQCDir}/X_QC/2_CR_high/chr_X \
     --freq \
     --hardy \
     --out ${GeneralQCDir}/X_QC/3_MAF_HWE/chr_X

###Eliminate outlier markers from  H-WE 
## awk '$5==0 {print $2}' ${GeneralQCDir}/X_QC/3_MAF_HWE/chr_X.frq > ${GeneralQCDir}/X_QC/3_MAF_HWE/zeroMAF_X.temp
  #merge both HW and MAF SNPs
 # cat ${GeneralQCDir}/X_QC/3_MAF_HWE/highhw_X.temp ${GeneralQCDir}/X_QC/3_MAF_HWE/zeroMAF_X.temp > ${GeneralQCDir}/X_QC/3_MAF_HWE/extr_Xhw
   ## extract markers with MAF==0 and WHE< 1x 10 exp(-6)
 # plink --bfile ${GeneralQCDir}/X_QC/2_CR_high/chr_X \
         --make-bed \
         --exclude ${GeneralQCDir}/X_QC/3_MAF_HWE/extr_Xhw \
         --out ${GeneralQCDir}/X_QC/3_MAF_HWE/chr_X

         #--remove ${GeneralQCDir}/X_QC/sex.exclude \

#cat ${GeneralQCDir}/X_QC/3_MAF_HWE/chr_X.bim|awk '{print$2}' > ${GeneralQCDir}/X_QC/3_MAF_HWE/incl_HW.snps

 


##################################################################################################
###########################-----------------PCA analysis---------------###########################
### filter reference Databases by the common SNPs
 plink --bfile ${ref1000G} \
       --extract  ${commonSNPs} \
       --make-bed \
       --out ${GeneralQCDir}/6_PCA/proc/thg

 plink --bfile ${gonlref} \
       --extract  ${commonSNPs}  \
       --make-bed \
       --out ${GeneralQCDir}/6_PCA/proc/goref
       

## filter also the cohort by common snps
 plink --bfile ${GeneralQCDir}/5_Relatedness/proc/full_data.no.dup \
       --extract ${commonSNPs}  \
       --make-bed \
       --out ${GeneralQCDir}/6_PCA/proc2/allchr_join


## make list of individuals to merge 
echo "${GeneralQCDir}/6_PCA/proc2/allchr_join" >> ${GeneralQCDir}/6_PCA/proc2/allfiles.list
echo "${GeneralQCDir}/6_PCA/proc/thg" >> ${GeneralQCDir}/6_PCA/proc2/allfiles.list
echo "${GeneralQCDir}/6_PCA/proc/goref" >> ${GeneralQCDir}/6_PCA/proc2/allfiles.list

###merge all data
 plink --merge-list ${GeneralQCDir}/6_PCA/proc2/allfiles.list \
       --merge-mode 3 \
       --out ${GeneralQCDir}/6_PCA/proc2/full_data
###remove SNPs with more tha 3 alleles if there are
 if [ -e ${GeneralQCDir}/6_PCA/proc2/full_data.missnp ];
  then
      plink --bfile ${GeneralQCDir}/6_PCA/proc/thg \
            --exclude ${GeneralQCDir}/6_PCA/proc2/full_data.missnp \
            --make-bed \
            --out ${GeneralQCDir}/6_PCA/proc2/allg.temp
   
      plink --bfile ${GeneralQCDir}/6_PCA/proc/goref \
            --exclude ${GeneralQCDir}/6_PCA/proc2/full_data.missnp \
            --make-bed \
            --out ${GeneralQCDir}/6_PCA/proc2/allgonl.temp
            
      plink --bfile ${GeneralQCDir}/6_PCA/proc2/allchr_join \
            --exclude ${GeneralQCDir}/6_PCA/proc2/full_data.missnp \
           --make-bed \
           --out ${GeneralQCDir}/6_PCA/proc2/allchr_join
              
              rm ${GeneralQCDir}/6_PCA/proc2/allfiles.list

      echo "${GeneralQCDir}/6_PCA/proc2/allchr_join" >> ${GeneralQCDir}/6_PCA/proc2/allfiles.list
      echo "${GeneralQCDir}/6_PCA/proc2/allg.temp" >> ${GeneralQCDir}/6_PCA/proc2/allfiles.list
      echo "${GeneralQCDir}/6_PCA/proc2/allgonl.temp" >> ${GeneralQCDir}/6_PCA/proc2/allfiles.list

      plink --merge-list ${GeneralQCDir}/6_PCA/proc2/allfiles.list \
            --out ${GeneralQCDir}/6_PCA/proc2/full_data
 fi

##extract the common snps again, filter by call rate and maf
 plink --bfile ${GeneralQCDir}/6_PCA/proc2/full_data \
      --extract  ${commonSNPs} \
      --maf 0.1 \
      --geno 0.01 \
      --make-bed \
      --out ${GeneralQCDir}/6_PCA/PCA_data1

###define sample cluster to project on
awk '{print $1,$2}' ${GeneralQCDir}/6_PCA/proc2/full_data.fam > ${GeneralQCDir}/6_PCA/proc2/merged_samples
awk 'BEGIN { FS=" " ;} {if ($2 ~ /DNA/) print $0 " Cohort"; else if($2 ~ /gonl/) print $0 " GoNl" ; else print $0 " 1000G" }' ${GeneralQCDir}/6_PCA/proc2/merged_samples> ${GeneralQCDir}/6_PCA/proc2/clusters

###make pca on ggonl thousand genomes and [plot PCA of cohort ]project chorot on them
plink -bfile ${GeneralQCDir}/6_PCA/PCA_data1 --within ${GeneralQCDir}/6_PCA/proc2/clusters --pca --pca-cluster-names 1000G GoNl --out ${GeneralQCDir}/6_PCA/PCA_1000G



###first PCA plot
Rscript /groups/umcg-aad/tmp04/umcg-elopera/PCA1.R -i ${GeneralQCDir}/6_PCA \
 -f "PCA_1000G.eigenvec" \
 -o ${GeneralQCDir}/plots \
 -r "1st."

### 2.nd PCA 
##extract selected samples by PCA
plink --bfile ${GeneralQCDir}/6_PCA/PCA_data1 \
      --keep ${GeneralQCDir}/6_PCA/1st.eur.samples \
      --maf 0.1 \
      --geno 0.01 \
      --make-bed \
      --out ${GeneralQCDir}/6_PCA/PCA_data2

###get the clustering file only for the new samples

awk '{print $2}' ${GeneralQCDir}/6_PCA/PCA_data2.fam> ${GeneralQCDir}/6_PCA/second_PCA.samples  
grep -f ${GeneralQCDir}/6_PCA/second_PCA.samples ${GeneralQCDir}/6_PCA/proc2/clusters > ${GeneralQCDir}/6_PCA/second_clusters

#2nd. PCA
plink -bfile ${GeneralQCDir}/6_PCA/PCA_data2 --within ${GeneralQCDir}/6_PCA/second_clusters --pca --pca-cluster-names GoNl 1000G --out ${GeneralQCDir}/6_PCA/second_PCA

#Rscript PCA1_v1.R -i ${GeneralQCDir}/6_PCA \
Rscript /groups/umcg-aad/tmp04/umcg-elopera/PCA1.R -i ${GeneralQCDir}/6_PCA \
 -f "second_PCA.eigenvec" \
 -o ${GeneralQCDir}/plots \
 -r "2nd."

cat ${GeneralQCDir}/6_PCA/1st.excl.samples ${GeneralQCDir}/6_PCA/2nd.excl.samples > ${GeneralQCDir}/6_PCA/excl_PCA.samples

for chr in {1..22} "XY"
 do
 ##extract selected samples by secondPCA
 plink --bfile ${GeneralQCDir}/4_Het/chr_${chr} \
       --keep ${GeneralQCDir}/6_PCA/2nd.incl.samples \
       --make-bed \
       --freq \
       --out ${GeneralQCDir}/6_PCA/chr_${chr}

 done

#######################################################
#################plot results###########################
Rscript /groups/umcg-aad/tmp04/umcg-elopera/genotypeQC_v9.R -i ${GeneralQCDir} \
 -o ${GeneralQCDir}/plots \
 -n "batch${UGLIbatch}" \
 -r ${MAFref} \

#########################################################
### Samples external concordance###

Rscript /groups/umcg-aad/tmp04/umcg-elopera/concordanceCheck_v1.2.R -p ${GeneralQCDir}/5_Relatedness/proc/full_data \
 -o ${GeneralQCDir}/plots/sample_concordance/ \
 -l $merged_conc1 \
 -g $merged_conc2 \
 -d $all_pairing_file


# opt$plink<-"/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC_err/5_Relatedness/proc/full_data"
# opt$out<-"/groups/umcg-aad/tmp04/umcg-elopera/merged_general_QC_err/plots/sample_concordance/"
# opt$llref<-"/groups/umcg-aad/tmp04/umcg-elopera/concordanceCheck/gwas_og_plink_chr/filtered_merged_LL_GWAS"
# opt$gonlref<- "/groups/umcg-aad/tmp04/umcg-elopera/concordanceCheck/gonl_og_plink_chr/filtered_merged_goNL_OG"
# opt$datapairing<- "/groups/umcg-aad/tmp04/umcg-elopera/concordanceCheck/ugli.final.pairing.concordanceCheck.txt"


