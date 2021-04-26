#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G

#############################################
### Merge multiple batches of GSA genotyping
### date=12-04-2019
### version 0.02 (unofficial)
### author: EL
#############################################
### New
## 26-04-2021
## commented
## added line to remove temporary directory
## changed memory required from 140G to 64G (untested)
#############################################

### Load packages
module load plink
module load RPlus

### Define variables (paths)
wkdir="/groups/umcg-aad/tmp04/umcg-elopera" ## the place used to process the data and make output
inputdir="/groups/umcg-ugli/tmp04/projects" ## directory with the folders fo each batch
tempdir="$wkdir/temp" ## a temporary directory

################################################################################
### Main
################################################################################
mkdir -p $tempdir #make temporary working directory

### list of batch identificator (not name, but the unique part of it)
find $inputdir/ -maxdepth 1  -name "*_part*" > $tempdir/batches.list
sed -i 's|'$inputdir/UGLI'||g' $tempdir/batches.list 


### convert files to plink format
for batch in $(cat $tempdir/batches.list)
 do

 if [ -e $inputdir/UGLI$batch/run03/ ];
  then
  genSampleDir="$inputdir/UGLI$batch/run03/results/gensamplefiles"
  else 
   if [ -e $inputdir/UGLI$batch/run02/ ];
    then
     genSampleDir="$inputdir/UGLI$batch/run02/results/gensamplefiles"
    else
     genSampleDir="$inputdir/UGLI$batch/run01/results/gensamplefiles" #.gen and .sample files directory for every batch
   fi
  fi
  
  if [ -e $genSampleDir ];
  then
   mkdir -p $tempdir/bfiles$batch 
   for chr in {1..22} "XY" "X" "MT" "Y"    
    do
     ### create plink files for each batch
      plink --data $genSampleDir/chr_$chr \
            --make-bed  \
            --out $tempdir/bfiles$batch/chr_$chr
      echo "$tempdir/bfiles$batch/chr_$chr" >> $tempdir/allb_$chr.list
    done
  fi
 done

#make results directory
mkdir -p $wkdir/merged_general_QC/

#merge all batches for every single chromosome
 for chr in {1..22} "XY" "X" "MT" "Y" 
    do
     plink --merge-list $tempdir/allb_$chr.list \
           --out $wkdir/merged_general_QC/chr_$chr
 done

rm -r $tempdir

echo "Done. Thank you for using me ;) . "


