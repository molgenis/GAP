 # GAP 
   Short for Genotyping Array Pipeline.
   Consist of the following workflow:

```
   ⎛¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯⎞
   ⎜                    Iscan writes IDAT files to GATTACA {01,02}machines                ⎜
   ⎜                                                                                      ⎜
   ⎝______________________________________________________________________________________⎠
                                         v
                                         v
                                         v
   ⎛¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯⎞
   ⎜                    AGCT Pipeline conversion of IDAT files to GTC files,                     ⎜
   ⎜                    and takes place on GATTACA {01,02}machines.                            ⎜
   ⎝______________________________________________________________________________________⎠
                                         v
                                         v  > > > > > > > > > > NGS_Automated CopyRawDataToPRM [stores .IDAT and .GTC files on permanent storage system]
                                         v  > > > > > > > > > > NGS_Automated Start Pipeline [automatically start Pipeline when new IDAT files are converted to GTC files ]
                                         v
   ⎛¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯⎞
   ⎜                         Run GAP Pipeline, to create a CallRates report per project   ⎜
   ⎜                         and a PennCNV report and VCF file per patient.               ⎜
   ⎜                         Takes place on [zinc-finger of leucine-zipper].               ⎜
   ⎝______________________________________________________________________________________⎠
                                         v
                                         v  > > > > > > > > > > NGS_Automated CopyProjectDataToPRM [stores data produced by pipeline on permanent storage system]
                                         v
   ⎛¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯⎞
   ⎜                  DARWIN, calculate standard deviation                                ⎜
   ⎜                  and stores all data in array database.                               ⎜
   ⎝______________________________________________________________________________________⎠
```

  #### GAP Pipeline 

The GAP pipeline consist of 7 steps:
1: Create_PennCNV_Report_diagnostics
```
This step creates a file containing per SNP information about the log ratio and the B allel frequency of the specific snp.

Fileformat:
Name '\t' Chromosome '\t' Position '\t' Sample1.GType '\t' Sample1.LogRratio '\t' Sample1.B Allele Freq '\t' Sample2.GType '\t' Sample2.LogRratio '\t' Sample2.B Allele Freq

The log R ratio and B allele Frequency per SNP are used by Nexus (commercial software) to call CNV's

This PennCNV report contains the information for all the samples with the same sentrixbarcode (glaasje nummer). 
When Darwin can process per samples PennCNV reports, this step is no longer neseccary and will be replaced by Create_PennCNV_report_diagnocticsPerSentrixBarcode.
```
2: Create_PennCNV_report_diagnocticsPerSentrixBarcode
```
This step creates a file per patient, containing per SNP information about the log ratio and the B allel frequency of the specific snp.
The log R ratio and B allele Frequency per SNP are used by Nexus (commercial software) to call CNV's

file format:

SNP Name '\t' Sample ID '\t' Chr '\t' Position '\t' Log R Ratio '\t' B Allele Freq

```

3: Create_Callrate_file_PerSentrixBarcode
```
This step creates a file per sentrix barcode, containing Callrate information per sample .
The Callrate gives an idea of how many percent of the SNPs on the array are performing well.
If this number is below 0.97 we know the data is of pour quality.

Fileformat:

Sample ID '\t' Call Rate '\t' Gender.

```
4: Make_Final_PennCNV_report

```
This step combines the PennCNV reports from the Create_PennCNV_Report_diagnostics step and makes one file,
containing all PennCNV values from one project.

When Darwin can process per samples PennCNV reports, this step is no longer neseccary and will be replaced by Create_PennCNV_report_diagnocticsPerSentrixBarcode.

```

5: Make_Callrate_FinalReport
```
In this step the files created in the Create_Callrate_file_PerSentrixBarcode step are combined. 
So all callrates of one project are in one file.


```
6: CalculateSD
```
Per patient the SD is calculated. SD should be below 0.2.
And for each patient a VCF file is created. This vcf file is used for the concordance check.

```
7: PipelineFinished





