#!/bin/bash
#SBATCH --job-name=ConcordanceCheck
#SBATCH --output=ConcordanceCheck.out
#SBATCH --error=ConcordanceCheck.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 15gb
#SBATCH --nodes 1
#SBATCH --open-mode=append

set -e
set -u


PARSED_OPTIONS=$(getopt -n "$0"  -o w:r:t: --long "workdir:ref-file:tooldir:"  -- "$@")

#
# Bad arguments, something has gone wrong with the getopt command.
#
if [ $? -ne 0 ]; then
        usage
        echo "FATAL: Wrong arguments."
        exit 1
fi


eval set -- "$PARSED_OPTIONS"

#
# Now goes through all the options with a case and using shift to analyse 1 argument at a time.
# $1 identifies the first argument, and when we use shift we discard the first argument, so $2 becomes $1 and goes again through the case.
#
while true; do
  case "$1" in
	-w|--workdir)
                case "$2" in
                *) workdir=$2 ; shift 2 ;;
            esac ;;

	-r|--reference)
                case "$2" in
                *) r=$2 ; shift 2 ;;
            esac ;;
    -t|--tooldir)
                case "$2" in
                *) tooldir=$2 ; shift 2 ;;
            esac ;; 
         --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
  esac
done



empty=""
#
# Check required options were provided.
#
if [[ -z "${workdir-}" ]]; then
        workdir="/groups/umcg-gdio/tmp04/umcg-kdelange/ConcordanceCheckGAP/" #"/groups/umcg-gd/tmp05/Concordance/"
fi
if [[ -z "${reference-}" ]]; then
        reference="/apps/data/1000G/phase1/human_g1k_v37_phiX.fasta"
fi
if [[ -z "${tooldir-}" ]]; then
	tooldir="/apps/software/ngs-utils/18.06.2/"
fi


echo "workdir: ${workdir}"
echo "reference: ${reference}"
echo "tooldir: ${tooldir}"

for vcffile in "${workdir}/NGS_vcf/"*".final.vcf"   # "${workdir}/ngs/"*".final.vcf"
do
    bedtype=($(grep -m 1 -o -P 'intervals=\[[^\]]*.bed\]' ${vcffile} | cut -d [ -f2 | cut -d ] -f1))
    echo "bedtype: ${bedtype}"
    beddir=($(dirname ${bedtype}))
    echo "beddir: ${beddir}"
    bedfile="${beddir}/captured.merged.bed"
    echo "${bedfile}"
    patientno=($(basename ${vcffile} .final.vcf))
    dnano=($(basename ${vcffile} | cut -d _ -f3 | cut -d A -f2))

    checkCount=$(ls "${workdir}/array_finalreport/concordance_DNA-${dnano}_"*".txt" | wc -l)
    arrayfile=""
   
   ##################
   ### maak mail met notificatie als er meer dan 1 file is met hetzelfde DNA nummer 
   ##################
   
    if [[ "${checkCount}" > 1 ]]
    then 
        echo "more than 1 file"
        messageGCC="Dear GCC helpdesk,\n\nThere is more than 1 sample with id:${arrayid}. \nPlease check if there is some think wrong with the concordance check.\nKind regards\nGCC"
        #printf '%b\n' "${messageGCC}" | mail -s "Concordance check error" 'helpdesk.gcc.groningen@gmail.com'
        exit 1
    else
        arrayfile=$(ls "${workdir}/array_finalreport/concordance_DNA-${dnano}_"*".txt")
    fi
    
    arrayid=($(basename ${arrayfile} .txt | cut -d _ -f2-5))
    arrayno=($(awk '{if ($1 ~ /:/ ){print $2}}' ${arrayfile} | uniq))
    
    sampleConcordanceFile=${arrayid}_ConcordanceFile.txt
    arrayTmpMap=${arrayid}.concordance.tmp.map
    arrayMapFile=${arrayid}.concordance.map
    familyList=${arrayid}_familyList.txt

    outputfolder=${workdir}/Output_Concordane
    #mkdir ${workdir}/Tmp_Concordance/${arrayid}
    reporttmpdir=${workdir}/Tmp_Concordance/${arrayid}
    
    echo "vcf-file:${vcffile}"
    echo "bedfile:${bedfile}"    
    echo "vcf-file-no:${patientno}"
    echo "dnano: ${dnano}"
    echo "arrayfile: ${arrayfile}"
    echo "arrayid: ${arrayid}"
    echo "arrayno: ${arrayno}"
    echo "outputfoler: ${outputfolder}"
    echo "reporttmpdir: ${reporttmpdir}"

    echo "starting"

    ######################################################################################

###Start protocol
    ### Input the final array file. Remove from the array file all the deletions and insertions. 
    ### The ${arrayid}_FinalReport_5.txt.tmp is the file with only SNPs, this one will be processed further.

    if test ! -e ${arrayfile};
    then
        echo "name, step, nSNPs, PercDbSNP, Ti/Tv_known, Ti/Tv_Novel, All_comp_het_called_het, Known_comp_het_called_het, Non-Ref_Sensitivity, Non-Ref_discrepancy, Overall_concordance" > ${reporttmpdir}/${sampleConcordanceFile}
        echo "[1] NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA" >> ${reporttmpdir}/${sampleConcordanceFile} 
    else
        #Check finalReport on "missing" alleles. Also, see if we can fix missing alleles somewhere in GenomeStudio
        awk '{ if ($3 != "-" || $4 != "-") print $0};' ${arrayfile} \
        > ${reporttmpdir}/${arrayid}_FinalReport.txt.tmp

        #Check finalreport on "D"alleles.
        awk '{ if ($3 != "D") print $0};' ${reporttmpdir}/${arrayid}_FinalReport.txt.tmp \
        > ${reporttmpdir}/${arrayid}_FinalReport_2.txt.tmp

        awk '{ if ( $4 != "D") print $0};' ${reporttmpdir}/${arrayid}_FinalReport_2.txt.tmp \
            > ${reporttmpdir}/${arrayid}_FinalReport_3.txt.tmp


        #Check finalreport on "I"alleles.
        awk '{ if ($3 != "I") print $0};' ${reporttmpdir}/${arrayid}_FinalReport_3.txt.tmp \
        > ${reporttmpdir}/${arrayid}_FinalReport_4.txt.tmp


        awk '{ if ($4 != "I") print $0};' ${reporttmpdir}/${arrayid}_FinalReport_4.txt.tmp \
        > ${reporttmpdir}/${arrayid}_FinalReport_5.txt.tmp

        #Push sample belonging to family "1" into list.txt
        echo 1 ${arrayno} > ${reporttmpdir}/${familyList}

    #########################################################################

        module load ngs-utils/18.06.2

        ##Create .fam, .lgen and .map file from sample_report.txt
        sed -e '1,10d' ${reporttmpdir}/${arrayid}_FinalReport_5.txt.tmp | awk '{print "1",$2,"0","0","0","1"}' | uniq > ${reporttmpdir}/${arrayid}.concordance.fam
        sed -e '1,10d' ${reporttmpdir}/${arrayid}_FinalReport_5.txt.tmp | awk '{print "1",$2,$1,$3,$4}' | awk -f $EBROOTNGSMINUTILS/RecodeFRToZero.awk > ${reporttmpdir}/${arrayid}.concordance.lgen
        sed -e '1,10d' ${reporttmpdir}/${arrayid}_FinalReport_5.txt.tmp | awk '{print $6,$1,"0",$7}' OFS="\t" | sort -k1n -k4n | uniq > ${reporttmpdir}/${arrayTmpMap}
        grep -P '^[123456789]' "${reporttmpdir}/${arrayTmpMap}" | sort -k1n -k4n > "${reporttmpdir}/${arrayMapFile}"
        grep -P '^[X]\s' "${reporttmpdir}/${arrayTmpMap}" | sort -k4n >> "${reporttmpdir}/${arrayMapFile}"
        grep -P '^[Y]\s' "${reporttmpdir}/${arrayTmpMap}" | sort -k4n >> "${reporttmpdir}/${arrayMapFile}" # still 630401 SNPs for patient 17935

    ######################################################################################
    
    ##Create .bed and other files (keep sample from sample_list.txt).

    ##Create .bed and other files (keep sample from sample_list.txt).

        module load PLINK/1.07-x86_64
        module list

        plink \
        --lfile ${reporttmpdir}/${arrayid}.concordance \
        --recode \
        --noweb \
        --out ${reporttmpdir}/${arrayid}.concordance \
        --keep ${reporttmpdir}/${familyList}

        module unload plink 
        module load plink/1.9-foss-2015b
        module list

    ##Create genotype VCF for sample
        plink \
        --recode vcf-iid \
        --ped ${reporttmpdir}/${arrayid}.concordance.ped \
        --map ${reporttmpdir}/${arrayMapFile} \
        --out ${reporttmpdir}/${arrayid}.concordance

    ##Rename plink.vcf to sample.vcf
        mv ${reporttmpdir}/${arrayid}.concordance.vcf ${reporttmpdir}/${arrayid}.genotypeArray.vcf ## nog 630431 snps

    ##Replace chr23 and 24 with X and Y
        perl -pi -e 's/^23/X/' ${reporttmpdir}/${arrayid}.genotypeArray.vcf
        perl -pi -e 's/^24/Y/' ${reporttmpdir}/${arrayid}.genotypeArray.vcf

    ##Create binary ped (.bed) and make tab-delimited .fasta file for all genotypes
        sed -e 's/chr//' ${reporttmpdir}/${arrayid}.genotypeArray.vcf | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2}}' \
        > ${reporttmpdir}/${arrayid}.genotypeArray.bed


    #Remove SNP`s from array which are not in a exon with the exon bedfile

        module load BEDTools/2.25.0-foss-2015b
        bedtools intersect -a ${reporttmpdir}/${arrayid}.genotypeArray.vcf -b ${bedfile} -header  >${reporttmpdir}/${arrayid}.genotypeArray.ExonFiltered.vcf ## nog 984 SNPs


    #Remove SNP's from array which are called homozygous reference
        awk '{ if ($10!= "0/0") print $0};' ${reporttmpdir}/${arrayid}.genotypeArray.ExonFiltered.vcf \
        > ${reporttmpdir}/${arrayid}.genotypeArray.ExonFiltered.HomozygousRefRemoved.vcf ## nog 1 SNP over

        sleep 3m

    #Count how much SNP's are in original VCF and how much proceed for Concordance
        wc -l ${reporttmpdir}/${arrayid}.genotypeArray.vcf > ${reporttmpdir}/${arrayid}_originalSNPs.txt
        wc -l ${reporttmpdir}/${arrayid}.genotypeArray.ExonFiltered.HomozygousRefRemoved.vcf > ${reporttmpdir}/${arrayid}_SNPswichproceedtoConcordance.txt

        #Change Array VCF to same name as NGS VCF
            awk '{OFS="\t"}{if ($0 ~ "#CHROM" ){ print $1,$2,$3,$4,$5,$6,$7,$8,$9,"'$arrayid'"} else {print $0}}' ${reporttmpdir}/${arrayid}.genotypeArray.ExonFiltered.HomozygousRefRemoved.vcf  > ${reporttmpdir}/${arrayid}.genotypeArray.ExonFiltered.HomozygousRefRemoved.FINAL.vcf

        #Making Array VCF index

            module load tabix/0.2.6-foss-2015b
        bgzip -c ${reporttmpdir}/${arrayid}.genotypeArray.ExonFiltered.HomozygousRefRemoved.FINAL.vcf > ${reporttmpdir}/${arrayid}.genotypeArray.ExonFiltered.HomozygousRefRemoved.FINAL.vcf.gz
            tabix -p vcf ${reporttmpdir}/${arrayid}.genotypeArray.ExonFiltered.HomozygousRefRemoved.FINAL.vcf.gz

	######################################################################################

	#Removing small indels from NGS VCF

        module load GATK/3.6-Java-1.8.0_74
        java -Xmx4g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
        -T SelectVariants \
        -R ${reference} \
        -V ${vcffile} \
        -o ${reporttmpdir}/${arrayid}.onlySNPs.FINAL.vcf \
        -selectType SNP

        #Change NGS VCF to same name as array ID
            awk '{OFS="\t"}{if ($0 ~ "#CHROM" ){ print $1,$2,$3,$4,$5,$6,$7,$8,$9,"'$arrayid'"} else {print $0}}' ${reporttmpdir}/${arrayid}.onlySNPs.FINAL.vcf  > ${reporttmpdir}/${arrayid}.onlySNPs.FINAL.ONE.vcf

    ######################################################################################

    ### Compare Array data with NGS vcf-File using GATK GenotypeConcordance

    ### Comparing VCF From NGS with Array VCF
        module load GATK/3.6-Java-1.8.0_74

        java -Xmx4g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
        -T GenotypeConcordance \
        -R ${reference} \
        -XL X:2699520-154931044 \
        -eval ${reporttmpdir}/${arrayid}.onlySNPs.FINAL.ONE.vcf \
        -comp ${reporttmpdir}/${arrayid}.genotypeArray.ExonFiltered.HomozygousRefRemoved.FINAL.vcf \
        -o ${reporttmpdir}/${arrayid}.GATK.VCF.Concordance.output.grp 

    ######################################################################################

    ## Making final output file withe the Overall Genotype Concordance score.
        
        finaloutputfile=${outputfolder}/${arrayid}.txt
        #awk '/Overall_Genotype_Concordance/ {if  == ${arrayid};print $1"\t"$4}' ${reporttmpdir}/${arrayid}.GATK.VCF.Concordance.output.grp >> ${finaloutputfile}
        grep -A2 Overall_Genotype_Concordance ${reporttmpdir}/${arrayid}.GATK.VCF.Concordance.output.grp | awk -v vari=${arrayid} '{if ($1 == vari){print $1,$4}}' >> ${finaloutputfile}
#     >> klaar: Als final output 1 txt bestand met ${arrayid}.txt met daar in 1 regel met ${arrayid}"\t"overall_GT_concordance waarde uit de ${arrayid}.GATK.VCF.Concordance.output.grp file.
#     dan teste met meer files uit de map van Marloes, neem ook DNA 97325 en 97310 mee, die moeten op 14% en 20% uitkomen. *marloes, hoe heb je dit bepaald? ik kom op 0% uit. Hier ook goed naar kijken.
#     Met MArlous overleggen over het aantal SNPs waarover concordence minimaal moet worden berekend. Martijn had het er over, dat dit het aantal snps moet zijn waardoor je uniek een patient kan callen. 
#     Eerst bepalen hoeveel SNPS je nodig hebt op een patient uniek te kunnen bepalen, dan kijken hoeveel snps we hebben om dat te doen. misschien daarna met Roan en Birgit een overlegje?

#     als klaar: in de NGS-utils repo klatsen.
#     Nog mer Roan overleggen over mappenstructuur, moet draaien op zinc-finger, gd/tmp05/Concordance
#     dan testen op zinc-finger, vcf files en array files van Marloes copieren naar zinc-finger
#     Dan met Darwin overleggen waar de output moet komen te staan. zinc-finger op de GD group, Jelko of Robert waar waar het moet komen te staan.

    fi
done