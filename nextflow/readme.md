#module
module load nextflow/21.10.4-Java-11-LTS

#
## nextflow command
#
```
nextflow run \
-profile local \ <slurm | local>
--samplesheet samplesheet.csv (Default gap samplesheet.)
main.nf (path the main nf script)
```

# optionals:
```
--gtcDir path where array slide directories are located.(See: gtcDir folder structure example )
-resume  (Resume a previously crashed run )
-with-report [file name] (Metrics about a workflow execution.)
-with-timeline [file name] (Timeline for all processes executed in your pipeline.)
```

#
## GAP Nextflow(NF) pipeline
#

The GAP NF pipeline consist of 5 steps:

1: gtcToFinalReport

```
Converts gct file the Illumina finalreport format using tool: Illumina: BeadArrayFiles

input: gtc files per sample, structured in gtcDir structure.
output: final report per arrayslide
```

2: mergeFinalReports

```
Converts the final report per arrayslide into one big final report.
```

3: finalReportToOptical

```
Converts the final report into Opticall input format.
```

4: OptiCall

```
Runs optiCall to call genotypes.

output: .calls and .probs files:

.calls format:
The output format is space-delimited with columns: rs, coordinate, allelesAB, pertubation value, call_1, call_2, call_3,.... 
The order of the calls is given by the order of the sample ids in the header line of the output file. 
The calls are encoded as 1 = AA, 2 = AB (heterozygote), 3 = BB, 4 = NN (no call).

.probs format:
The probs file gives probabilities for the samples in the header. For each sample at each SNP there are four probabilities, 
in the order: P(AA) P(BB) P(AB) P(NN). In cases where the maximum genotype probability is less than the probability threshold, 
the call will be NN but the posterior probabilities might not have P(NN) as the highest value.
( https://opticall.bitbucket.io )

```

5: OptiCallToGenSample

```
Converts Opticall output files into .gen and .sample Oxford file format.
```

#
## gap samplesheet format (comma separated csv)
#

Example:
```
Sample_ID,Sample_Plate,Sample_Name,Project,AMP_Plate,Sample_Well,SentrixBarcode_A,SentrixPosition_A,Scanner,Date_Scan,Replicate,Parent1,Parent2,Gender,pipeline,owner,MA1,MA2,MSM,FMS,PM1,RA1,PB2,PB1,LX1,LX2,EML,SML,ATM,XC3,XC4,sampleType,manifest,egt
1234,WG1234567,DNA-1,projectName,WG1234567,A09,12345678910,R09C02,,28-7-2021,,,,m,research,,WG1234567,WG1234567,WG1234567,WG1234567,WG1234567,WG1234567,WG1234567,WG1234567,,,,,,,,GAP,GSAMD-24v3-0-EA_20034606_A1.bpm,referentie_GSAMD_V3_20210115.egt
```	

Important columns to mention:

- Sample_Name: Is used the replace SentrixBarcode_A + SentrixPosition_A the correponding samplename.
- Project: This name is used for projectDirs and project file names.
- SentrixBarcode_A: is used to locate the directiory containing the correponding sample the gtc files in the gtcDir.
- SentrixPosition_A: describes the filename of the gtc file, which is used to map the sampleIds to the correct gtc file and barcode.
- gender: Gender of the sample what will be used to detect samplesswaps.
- manifest: Manifest files describe the SNP or probe content on a genotyping array.(BPM file format)
- egt: Cluster files describe the cluster positions for the Illumina genotyping array and
       are used when analyzing data to make the genotype call.( EGT file format )

#
## expected gtcDir folder structure:
#
```
├── gtcDir:
│   ├── 203273230154 ( SentrixBarcode_A )
│   │   ├── 203273230154_R08C02.gtc ( SentrixBarcode_A + SentrixPosition_A .gtc )
│   │   ├── 203273230154_R10C02.gtc
│   │   └── 203273230154_R12C02.gtc
            etc...
│   ├── 203281940081
│   │   ├── 203281940081_R02C01.gtc
│   │   ├── 203281940081_R09C01.gtc
│   │   └── 203281940081_R11C01.gtc
            etc...
```
