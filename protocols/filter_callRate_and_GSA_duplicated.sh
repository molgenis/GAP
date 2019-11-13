#MOLGENIS walltime=01:59:00 mem=2gb ppn=1

#string InputDir
#list chr
#string GeneralQCDir
#string callrateCutoff

set -e
set -u

module load plink
module load RPlus

makeTmpDir "${GeneralQCDir}"
tmpGeneralQCDir="${MC_tmpFile}"


###create working directories
mkdir -p "${GeneralQCDir}/"
mkdir -p "${GeneralQCDir}/0_pre"
mkdir -p "${GeneralQCDir}/1_CR80/"
mkdir -p "${GeneralQCDir}/2_CR_high/"


######### step 1. Autosome Call rate (first filter <80%)#########################
#################################################################################

### failsafe_1. Will stop if any chromosome plink file is missing
for chrNr in ${chr[@]}
  do
   chr_file=${InputDir}/$chrNr
   for type in ".bed" ".bim" ".fam"
    do
     if [ -f "$chr_file$type" ]; then
      :
     else
      echo "$chr_file does not exist"
      exit 1
     fi
    done
  done

### step 1a. call rate evaluation for markers
for chrNr in ${chr[@]}
 do
 ## create call_rate stats for individuals and markers
   plink --bfile ${InputDir}/0_pre/${chrNr} \
         --missing \
         --out ${GeneralQCDir}/0_pre/${chrNr}
 ## create list of samples to exclude on the criteria callrate<=80
   awk '$6>0.20 {print $1, $2}' ${GeneralQCDir}/0_pre/${chrNr}.imiss > ${GeneralQCDir}/1_CR80/${chrNr}.extr80_sam.temp
 done

### step 1b. create list of duplicate markers to remove; keeping the ones with the best call rate
Rscript  $EBROOTGAP/scripts//position_duplicates.R -i ${GeneralQCDir}/0_pre

## creates list of individuals to be excluded for all the autosomes
cat ${GeneralQCDir}/1_CR80/chr_*.extr80_sam.temp |sort -u > ${GeneralQCDir}/1_CR80/extr80.samples
## remove temproary file
rm ${GeneralQCDir}/1_CR80/chr_*.extr80_sam.temp
## create files with duplicate SNPs to exclude
cat ${GeneralQCDir}/0_pre/chr_*.excl.duplicates > ${GeneralQCDir}/0_pre/extr.dups
## remove unecessary intermediate files
rm ${GeneralQCDir}/0_pre/*.excl.duplicates
rm ${GeneralQCDir}/0_pre/*nosex*

### Step 1c. removing markers and duplicates from data while calculating missing rate for individuales
for chrNr in ${chr[@]}
 do
 ## exclude individuals with callrate<=80 (creates data set with individual callrate>80)
   plink --bfile ${GeneralQCDir}/0_pre/${chrNr}  \
         --make-bed \
         --remove ${GeneralQCDir}/1_CR80/extr80.samples \
         --exclude ${GeneralQCDir}/0_pre/extr.dups \
         --allow-no-sex \
         --out ${GeneralQCDir}/1_CR80/${chrNr}

## calculates callrate stats from the previously filtered datafile
   plink  --bfile ${GeneralQCDir}/1_CR80/${chrNr} \
          --missing \
          --allow-no-sex \
          --out ${GeneralQCDir}/1_CR80/${chrNr}
 ## creates list with markers to remove for each chromosome
 awk '$5>0.20 {print $2}' ${GeneralQCDir}/1_CR80/${chrNr}.lmiss > ${GeneralQCDir}/1_CR80/${chrNr}.extr80_var.temp
 done

## Concatenate file of snps to exclude
cat ${GeneralQCDir}/1_CR80/chr_*.extr80_var.temp > ${GeneralQCDir}/1_CR80/extr80.vars
## removes the unnecesary files
rm ${GeneralQCDir}/1_CR80/*.temp
rm ${GeneralQCDir}/1_CR80/*nosex*

### Step 1d. exclude markers with call rate <80% and re-calculates missing rate
 for chrNr in ${chr[@]}
 do
   # exclude  SNPs with CR<80
   plink --bfile ${GeneralQCDir}/1_CR80/${chrNr}  \
         --make-bed \
         --exclude ${GeneralQCDir}/1_CR80/extr80.vars\
         --out ${GeneralQCDir}/1_CR80/${chrNr}.2

   #calculates callrate stats from the previously filtered datafile
   plink  --bfile ${GeneralQCDir}/1_CR80/${chrNr}.2 \
          --missing \
          --out ${GeneralQCDir}/1_CR80/${chrNr}.2

   ##create list of samples to exclude on the criteria callrate<=high (default 99%)
   awk '$6>${callrateCutoff} {print $1, $2}' ${GeneralQCDir}/1_CR80/chr_${chr}.2.imiss > ${GeneralQCDir}/2_CR_high/${chrNr}.extrhigh_sam.temp
   ## creates information file for the heterozygosity analysis (step 4)
   awk '$6<${callrateCutoff} {print $1, $2,$6}' ${GeneralQCDir}/1_CR80/chr_${chr}.imiss > ${GeneralQCDir}/2_CR_high/${chrNr}.incl_CR_sam
 done
## create complete list of samples to remove
cat ${GeneralQCDir}/2_CR_high/chr_*.extrhigh_sam.temp|sort -u > ${GeneralQCDir}/2_CR_high/extrhigh.samples

###############################################################################################
#### Step 2. Remove individuals and markers with call rate<=high (default 99%) ################

### Step 2a. remove markers with 80% call rate and missing report for the hight call rate threshold
 for chrNr in ${chr[@]}
 do
   ## exclude individuals with callrate<=high (creates data set with  individual callrate>99)
   plink --bfile ${GeneralQCDir}/1_CR80/${chrNr}  \
         --make-bed \
         --remove ${GeneralQCDir}/2_CR_high/extrhigh.samples \
         --out ${GeneralQCDir}/2_CR_high/${chrNr}

  #calculates callrate stats from the previously filtered datafile
   plink  --bfile ${GeneralQCDir}/2_CR_high/${chrNr} \
          --missing \
          --out ${GeneralQCDir}/2_CR_high/${chrNr}
  ## create excluding list for each chromosome
 awk '$5>${callrateCutoff} {print $2}' ${GeneralQCDir}/1_CR80/${chrNr}.lmiss > ${GeneralQCDir}/2_CR_high/${chrNr}.extrhigh_var.temp
 done

## create list full list of markers to be excluded
cat ${GeneralQCDir}/2_CR_high/chr_*.extrhigh_var.temp > ${GeneralQCDir}/2_CR_high/extrhigh.vars
rm ${GeneralQCDir}/2_CR_high/*.temp ##remove unnecesary files

### Step 2b. remove markers with call rate less than the high threshold (default 99%)
 for chrNr in ${chr[@]}
 do
   ## exclude SNPs
   plink --bfile ${GeneralQCDir}/2_CR_high/${chrNr}  \
         --make-bed \
         --exclude ${GeneralQCDir}/2_CR_high/extrhigh.vars \
         --out ${GeneralQCDir}/2_CR_high/${chrNr}.2
 done
### Step 2c. make summary files for the call rate process. creates merged files of included individuals and SNPs to be used of further analysis
cat ${GeneralQCDir}/0_pre/chr_*.fam|sort -u|awk '{print$2}' > ${GeneralQCDir}/0_pre/full.ind
cat ${GeneralQCDir}/0_pre/chr_*.bim|awk '{print$2}' > ${GeneralQCDir}/0_pre/untouched.snps
sort ${GeneralQCDir}/0_pre/extr.dups ${GeneralQCDir}/0_pre/untouched.snps|uniq -u >${GeneralQCDir}/0_pre/full.snps

cat ${GeneralQCDir}/1_CR80/chr_*.2.fam|sort -u|awk '{print$2}' > ${GeneralQCDir}/1_CR80/incl80.samples
cat ${GeneralQCDir}/1_CR80/chr_*.2.bim|awk '{print$2}' > ${GeneralQCDir}/1_CR80/incl80.vars

cat ${GeneralQCDir}/2_CR_high/chr_*.2.fam|sort -u|awk '{print$2}' > ${GeneralQCDir}/2_CR_high/inclhigh.samples
cat ${GeneralQCDir}/2_CR_high/chr_*.2.bim|awk '{print$2}' > ${GeneralQCDir}/2_CR_high/inclhigh.vars

## remove chromosome data files in 0_pre
rm  ${GeneralQCDir}/0_pre/chr_*.bim
rm  ${GeneralQCDir}/0_pre/chr_*.fam
rm  ${GeneralQCDir}/0_pre/chr_*.bed
## remove chromosome data files in 1_CR80
rm  ${GeneralQCDir}/1_CR80/chr_*.bim
rm  ${GeneralQCDir}/1_CR80/chr_*.fam
rm  ${GeneralQCDir}/1_CR80/chr_*.bed

#######################################
