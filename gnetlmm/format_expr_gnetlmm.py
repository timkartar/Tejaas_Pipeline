import sys, os
sys.path.append("../")
from iotools import readgtf
from collections import defaultdict
import pandas as pd
import re
import argparse


def parse_args():

	parser = argparse.ArgumentParser(description='format expression files and auxiliary files with \
												  headers and gene positions, and filter genes not in gencode!')

	parser.add_argument('--input',
						type=str,
						dest='infile',
						metavar='FILE',
						help='input expression file')

	parser.add_argument('--outfile',
						type=str,
						dest='outfile',
						metavar='FILE',
						help='output file prefix for GNETLMM formatted expr files')

	opts = parser.parse_args()
	return opts

if __name__ == '__main__':

	opts = parse_args()

	gtfpath = "/cbscratch/franco/datasets/gtex/gencode.v19.annotation.gtf.gz"
	gene_info = readgtf.gencode_v12(gtfpath, trim=False)

	gene_dict = defaultdict(list)
	for g in gene_info:
		gene_dict[g.ensembl_id] = [g.chrom, g.start, g.end]

	# GXFILE="/cbscratch/franco/datasets/gtex/Whole_Blood_Analysis.v6p.normalized.expression.txt.gencode_filter"
	GXFILE=opts.infile
	# OUTGX_gnet="/cbscratch/franco/datasets/gtex/gnetlmm/Whole_Blood_Analysis.v6p.normalized.expression.txt.gencode_filter.gnet"
	OUTGX_gnet=opts.outfile

	gx_df = pd.read_table(GXFILE, sep="\t", header=0, index_col=0)
	row_names = gx_df.index
	col_names = gx_df.columns

	filtered=0
	with open(OUTGX_gnet+".rows", 'w') as outstream:
		outstream.write("gene_ids gene_chrom gene_start gene_end\n")
		for f in row_names:
			if len(gene_dict[f]) > 0:
				line = "{:s} {:d} {:d} {:d}\n".format(f, gene_dict[f][0], gene_dict[f][1], gene_dict[f][2])
				outstream.write(line)
			else:
				filtered += 1
				# print(f)
				# raise ValueError("This gene is not in our gencode file!")

	print("Filtered out {:d} ensemble ids".format(filtered))			

	# export col info 
	with open(OUTGX_gnet+".cols", 'w') as outstream:
		outstream.write("fid\n")
		for f in col_names:
			outstream.write(f+"\n")

	# export expression matrix with no headers
	gx_df.to_csv(OUTGX_gnet+".matrix", header=False, index=False, doublequote=False, sep=" ")