#!/bin/python

from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest
from datetime import datetime
import sys, os, argparse

parser = argparse.ArgumentParser("Generate a final report from a directory of GTC files")
parser.add_argument("manifest", help="BPM manifest file")
parser.add_argument("gtc_directory", help="Directory containing GTC files")
parser.add_argument("output_directory", help="Directory where output has to be written")

args = parser.parse_args()

try:
    manifest = BeadPoolManifest(args.manifest)
except:
    sys.stderr.write("Failed to read data from manifest\n")
    sys.exit(-1)

samples = []
for gtc_file in os.listdir(args.gtc_directory):
	if gtc_file.lower().endswith(".gtc"):
		samples.append(gtc_file)


delim="\t"
index = 0

for gtc_file in samples:
	index += 1

	names = manifest.names
	output = open(os.path.join(args.output_directory, os.path.basename("concordance_" + gtc_file + ".txt")), "w")
    	output.write("[Header]\n")
    	output.write(delim.join(["Processing Date", datetime.now().strftime("%m/%d/%Y %I:%M %p")])+ "\n")
    	output.write(delim.join(["Content", os.path.basename(args.manifest)]) + "\n")
    	output.write(delim.join(["Num SNPs", str(len(names))]) + "\n")
    	output.write(delim.join(["Total SNPs", str(len(names))]) + "\n")
	output.write("[Data]\n")
    
	output.write(delim.join(["SNP Name", "Sample Index", "Allele1 - Forward", "Allele2 - Forward", "GC Score", "Chr", "Position"]) + "\n")

    	sys.stdout.write("Processing " + gtc_file + "\n")
    	gtc_file = os.path.join(args.gtc_directory, gtc_file)
	
	gtc = GenotypeCalls(gtc_file)
	genotypes = GenotypeCalls( gtc_file ).get_genotypes()
	chrom = manifest.chroms
	map_info = manifest.map_infos 
	forward_strand_genotypes = gtc.get_base_calls_forward_strand(manifest.snps, manifest.source_strands)
	gen_score=GenotypeCalls(gtc_file).get_genotype_scores()

	for (names, forward_strand_genotypes, chrom, map_info, gen_score) in zip(names, forward_strand_genotypes, chrom, map_info, gen_score):
	
		if forward_strand_genotypes[0] == "-":
			forward_strand_genotypes += "-"

		output.write(delim.join([names, str(index), forward_strand_genotypes[0], forward_strand_genotypes[1], str(round(gen_score,4)), chrom, str(map_info)]) + "\n")
	output.close()
