#!/usr/bin/env python

import sys
import argparse
import os
import csv
import re
from datetime import datetime
from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest, code2genotype, NormalizationTransform

delim = "\t"

#get commandline parameters
parser = argparse.ArgumentParser("Generate a final report from a directory of GTC files")
parser.add_argument("--manifest", help="BPM manifest file",required=True)
parser.add_argument("--samplesheet", help="Illumina sampleheet",required=False)
parser.add_argument("--excludeGTCFileIDs", help="txt file with the array IDs to be excluded (not need to add the '.gtc' extension)",required=False)
parser.add_argument("--gtc_directory", help="Directory containing GTC files",required=True)
parser.add_argument("--output_file", help="Location to write report",required=True)

args = parser.parse_args()

if os.path.isfile(args.output_file):
    sys.stderr.write("Output file already exists, please delete and re-run\n")
    sys.exit(-1)

try:
    manifest = BeadPoolManifest(args.manifest)
    samplesheet = args.samplesheet
except:
    sys.stderr.write("Failed to read data from manifest\n")
    sys.exit(-1)

# if samplesheet is given, read it and use it to replace
if args.samplesheet:
    counter=0
    found='false'
    # find rownumber for [Data] row
    with open(samplesheet, 'r') as fh:
        for line in fh :

           if re.match("^\[Data\].*", line):
               counter+=1
               found='true'
               break
           else:
              counter+=1
    fh.close()

    if found == 'false':
        counter=0

    # read sampleheet after [Data] row, and fill sampleDict with Barcode_positions and samplesIDs.
    with open(samplesheet, 'r') as f1:

       reader = csv.DictReader(f1.readlines()[counter:])
       headers = reader.fieldnames
       rowcount=0
       sampleDict ={}
       for row in reader:
           rowcount+=1
           sampleDict[row['SentrixBarcode_A'] + "_" + row['SentrixPosition_A'].upper()] =  row['Sample_ID']
    f1.close()

    if not rowcount == len(sampleDict.keys()):
        print ("Rowcount in samplesheet does not match number of samples in sampleDict: " + str(rowcount) + " vs " + str(len(sampleDict.keys())))
        exit(1)

elif not args.samplesheet:
    print("No samplesheet given.")
    sampleDict="NULL"


#If there is gtc files to be excluded, read from given file and create exlcuding dictionary
excludeIDArray={}
if args.excludeGTCFileIDs:
    file = open(args.excludeGTCFileIDs,"r")
    for line in file:
       #if IDs in file include the .gtc extension, create the without modifyingthe IDs, else add the extension.
        if '.gtc' in line:
         excludeIDArray[line.rstrip('\n')]=line.rstrip('\n')
        else:
         excludeIDArray[line.rstrip('\n')+'.gtc']=line.rstrip('\n')+'.gtc'

	
	
with open(args.output_file, "w") as output_handle:
    output_handle.write("[Header]\n")
    output_handle.write(delim.join(["Processing Date", datetime.now().strftime("%m/%d/%Y %I:%M %p")]) + "\n")
    output_handle.write(delim.join(["Content", os.path.basename(args.manifest)]) + "\n")
    output_handle.write(delim.join(["Num SNPs", str(len(manifest.names))]) + "\n")
    output_handle.write(delim.join(["Total SNPs", str(len(manifest.names))]) + "\n")

    samples = []
    for gtc_file in os.listdir(args.gtc_directory):
        if gtc_file.lower().endswith(".gtc"):
            if gtc_file in excludeIDArray:
                print(gtc_file + ' is excluded.')
            else:
                samples.append(gtc_file)


    output_handle.write(delim.join(["Num Samples", str(len(samples))]) + "\n")
    output_handle.write(delim.join(["Total Samples", str(len(samples))]) + "\n")
    output_handle.write("[Data]\n")
    output_handle.write(delim.join(["SNP Name", "Sample ID", "Chr", "MapInfo" ,"SNP" ,"REFSTRAND","SOURCESTRAND" ,"GType" ,"Allele1 - Top","Allele2 - Top","X_raw","Y_raw","X","Y","B Allele Freq","Log R Ratio","GT Score","Alleles - Plus", "Alleles - Forward"]) + "\n")

    for gtc_file in samples:
        sys.stderr.write("Processing " + gtc_file + "\n")
        gtc_file = os.path.join(args.gtc_directory, gtc_file)
        gtc = GenotypeCalls(gtc_file)
        genotypes = gtc.get_genotypes()
        raw_Xs = gtc.get_raw_x_intensities()
        raw_Ys = gtc.get_raw_y_intensities()
        normalized_intensities = GenotypeCalls( gtc_file ).get_normalized_intensities(manifest.normalization_lookups)
        base_calls = GenotypeCalls( gtc_file ).get_base_calls()
        logratios = GenotypeCalls( gtc_file ).get_logr_ratios()
        BAFs = GenotypeCalls( gtc_file ).get_ballele_freqs()
        genotype_scores = GenotypeCalls( gtc_file ).get_genotype_scores()
        plus_strand_genotypes = gtc.get_base_calls_plus_strand(manifest.snps, manifest.ref_strands)
        forward_strand_genotypes = gtc.get_base_calls_forward_strand(manifest.snps, manifest.source_strands)

        assert len(genotypes) == len(manifest.names)

        for (name, chrom, map_info, snp,ref_strand,source_strands, genotype, norm, Allele, raw_x, raw_y, BAF, logratio , genotype_score , ref_strand_genotype, source_strand_genotype) in zip(manifest.names, manifest.chroms, manifest.map_infos, manifest.snps ,manifest.ref_strands,manifest.source_strands,genotypes , normalized_intensities, base_calls, raw_Xs, raw_Ys, BAFs , logratios , genotype_scores , plus_strand_genotypes, forward_strand_genotypes):
            x_norm = norm[0]
            y_norm = norm[1]

            COMPLEMENT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C", "D": "D", "I": "I"}

            if ref_strand == 2 :
               #[A/G] to [T/C] where ref_strand = 2
               new_snp = '['+ COMPLEMENT_MAP[snp[1]] + '/' + COMPLEMENT_MAP[snp[3]] + ']'
            else:
               new_snp=snp
            # if samplesheet is given, replace Barcode_positions with sampleID.
            if args.samplesheet:
                if os.path.basename(gtc_file)[:-4] in sampleDict:
		    sampleName = sampleDict[os.path.basename(gtc_file)[:-4]]
                else:
                    sampleName = os.path.basename(gtc_file)[:-4]
            else:
                sampleName = os.path.basename(gtc_file)[:-4]

            output_handle.write(delim.join([name, sampleName, str(chrom) , str(map_info) , str(new_snp), str(ref_strand),str(source_strands), code2genotype[genotype] , Allele[:-1] , Allele[1:] , str(raw_x), str(raw_y) ,str(x_norm), str(y_norm),str(BAF),str(logratio),str(genotype_score),ref_strand_genotype, source_strand_genotype]) + "\n")
