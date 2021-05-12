# Quality Control for GSA array on the UMCG HPC
Created by: Raul Aguirre-Gamboa, Esteban Lopera-Maya. \
Contact information: e.a.lopera.maya@umcg.nl. \

# Pipeline structure
Important note: this is a self-contained and specific pipeline to work on the UMCG HPC. This means, that at the time it is not fully automated to work in completeness in every kind of data and every kind of server, and many steps might be specifically designed for the UGLI data alone. You are welcome, however, to take any individual step and/or code snippet to addapt to your own conditions. \
\
This pipeline is designed to work in three main steps as indicated by the numeration in front of the name of the main scripts. Numbered scripts should be called by the user independently one ofter the other. All the scripts with the prefix "sub_" in the front of the name are automatically by the numbered scripts. \
A feedback loop should be done manually by the user, after processing also manually the output of step 2 (script 2.) to remove familial errors (as these cannot be removed automatically), the data resulting from this should go through step 2 a second time. Special adittional "sub_" scripts contain also the suffix "second_it". These should replace their counterparts in the second iteration of step 2. \

#  content

1.merge_genotype_batches.sh: This script has 2 functions \
- Transform oxford files into plink format \
- Search and and merge different batches into a single big batch \

2.QC_autosomes_launch.sh: Main script. It launches all quality control steps and uses most of the required and other reference files \
- Loads tools and reference files for all steps \
- Corrects possible sample ID differences \
- Performs call rate filtering and removes/renames duplicated SNPs \
- Performs MAF and HWE filtering \
- Performs heterozygozity filtering (requires change in second iteration) \
- Relatedness and identity by descent and duplicated samples identification (requires change in second iteration). Results are recommended to be analysed manually \
- Performs basic quality control for chromosome X (call rate and duplicate SNPs filter only) \
- Impute sex from genotype and sex check \
- PCA analysis merges the data with 1000Genomes and GoNL and projects PCs based on these. It flags but does not remove non-european samples \
- Plots output from all steps and also performs external concordance analyses. (requires external reference files) \
- Performs also internal concordance anlyses (requires internal reference files) \

3.MendelianErrors_FounderStats.R: This R script makes some additional steps \ 
- Take as input a merged QCed chromosomes as well as the chomosome X \
- Uses corrected familial information to create founders-only dataset \
- Calculates founder stats \
- Calculates mendelian errors \
- Calculates MAF and HWE filters for chormosome X with females-only

4. Pre-inputation steps: located in the folder <Imputation>. This is the process to prepare the data for inputation with the HRC reference panel, following the steps indicated by the Sanger imputation server (https://imputation.sanger.ac.uk/). \
- Remove insertions and deletions (script with number 1) \
- Remove duplicated snps (same position SNPS) (script with number 2) \
- transform to VCF and format for imputation, fix reference and alternative alleles in each SNP (use fixref plugin form BCFtools). \
***Important note: post-imputation quality control steps are not included here. This included only internal concordance evaluation with internal biobanks data. A posterior evaluation of the resulting allelic frequencies with HRC revealed ~5K SNPs with wrong reference allele asignation. These need to be removed form the current data imputed and can be found in the gearshift cluster at: "/groups/umcg-lifelines/prm03/releases/gsa_genotypes/v1/UGLI_imputed_outliers.list"
  
Acknowledgements:
This pipeline was elaborated with the financial aid from UMCG-HAP2017 as part of the project, as well as the Conacyt Fellowship, which supported RAG, and Colciencias Fellowship, which supported ELM.
  


