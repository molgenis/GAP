import sys
import argparse
import os
from datetime import datetime
from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest, code2genotype, NormalizationTransform

delim = "\t"

parser = argparse.ArgumentParser("Generate a final report from a directory of GTC files")
parser.add_argument("manifest", help="BPM manifest file")
parser.add_argument("gtc_directory", help="Directory containing GTC files")
parser.add_argument("output_file", help="Location to write report")

args = parser.parse_args()

if os.path.isfile(args.output_file):
    sys.stderr.write("Output file already exists, please delete and re-run\n")
    sys.exit(-1)

try:
    manifest = BeadPoolManifest(args.manifest)
except:
    sys.stderr.write("Failed to read data from manifest\n")
    sys.exit(-1)

with open(args.output_file, "w") as output_handle:
    output_handle.write("[Header]\n")
    output_handle.write(delim.join(["Processing Date", datetime.now().strftime("%m/%d/%Y %I:%M %p")])+ "\n")
    output_handle.write(delim.join(["Content", os.path.basename(args.manifest)]) + "\n")
    output_handle.write(delim.join(["Num SNPs", str(len(manifest.names))]) + "\n")
    output_handle.write(delim.join(["Total SNPs", str(len(manifest.names))]) + "\n")

    samples = []
    for gtc_file in os.listdir(args.gtc_directory):
        if gtc_file.lower().endswith(".gtc"):
            samples.append(gtc_file)

    output_handle.write(delim.join(["Num Samples", str(len(samples))]) + "\n")
    output_handle.write(delim.join(["Total Samples", str(len(samples))]) + "\n")

    output_handle.write("[Data]\n")
    output_handle.write(delim.join(["SNP Name", "Sample ID", "Chr", "MapInfo" ,"SNP" , "GType" ,"Allele1 - Top","Allele2 - Top","X_raw","Y_raw","X","Y","B Allele Freq","Log R Ratio","GT Score","Alleles - Plus", "Alleles - Forward"]) + "\n")
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
        for (name, chrom, map_info, snp, genotype, norm, Allele, raw_x, raw_y, BAF, logratio , genotype_score , ref_strand_genotype, source_strand_genotype) in zip(manifest.names, manifest.chroms, manifest.map_infos, manifest.snps ,genotypes , normalized_intensities, base_calls, raw_Xs, raw_Ys, BAFs , logratios , genotype_scores , plus_strand_genotypes, forward_strand_genotypes):

            x_norm = norm[0]
            y_norm = norm[1]

            output_handle.write(delim.join([name, os.path.basename(gtc_file)[:-4], str(chrom) , str(map_info) , str(snp), code2genotype[genotype] , Allele[:-1] , Allele[1:] , str(raw_x), str(raw_y) ,str(x_norm), str(y_norm),str(BAF),str(logratio),str(genotype_score),ref_strand_genotype, source_strand_genotype]) + "\n")

