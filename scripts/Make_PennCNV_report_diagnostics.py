#!/bin/python

from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest, code2genotype
import sys
import os
import glob
import argparse

parser = argparse.ArgumentParser("Generate a final report from a directory of GTC files")
parser.add_argument("manifest", help="BPM manifest file")
parser.add_argument("gtc_directory", help="Directory containing GTC files")
parser.add_argument("output_directory", help="Directory where output has to be written")

args = parser.parse_args()

samples = []
for gtc_file in os.listdir(args.gtc_directory):
	if gtc_file.lower().endswith(".gtc"):
		samples.append(gtc_file)

for gtc_file in glob.glob(os.path.join(args.gtc_directory, '*.gtc')):
	# print gtc_file
	gtc_output = gtc_file + ".txt"
	gtc_output2 = os.path.basename(gtc_output)
	save_path = args.output_directory
	gtc_output3 = os.path.join(save_path, gtc_output2)
	print gtc_output3
	
	
	output = open(gtc_output3, "w")
	print "output is" + str(output)
	if gtc_file.lower().endswith(".gtc"):
		manifest = BeadPoolManifest(args.manifest)
		names = manifest.names
		sample_id = os.path.basename(gtc_file)[:-4]
		print sample_id
		genotypes = GenotypeCalls( gtc_file ).get_genotypes()
		chrom = manifest.chroms
		map_info = manifest.map_infos 
		map_info2 = str(map_info)
		logratio = GenotypeCalls( gtc_file ).get_logr_ratios()
		BAF = GenotypeCalls( gtc_file ).get_ballele_freqs()
		for (names, chrom, map_info, genotypes, logratio, BAF) in zip( names, chrom, map_info, genotypes, logratio, BAF):
			output.write(sample_id + "\t" + names + "\t" + chrom + "\t" + str(map_info) + "\t" + code2genotype[genotypes] + "\t" + str(logratio) + "\t" + str(BAF) + "\n")
	output.close()

