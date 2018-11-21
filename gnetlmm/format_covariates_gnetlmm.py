import pandas as pd
import os, re
import argparse

# COVFILE = "/cbscratch/franco/datasets/gtex/PEER_correction/chr22/chr22_gtex_wholeblood_normalized_expression_PEER_covariates.txt"
# COVFILE is:
# 1st line is header with donor ids
# then, 1 covariate per line

def parse_args():

	parser = argparse.ArgumentParser(description='transpose covariate file and frop headers')

	parser.add_argument('--input',
						type=str,
						dest='infile',
						metavar='FILE',
						help='input covariates file')

	parser.add_argument('--output',
						type=str,
						dest='outfile',
						metavar='FILE',
						help='output GNETLMM formatted covariate file')

	parser.add_argument('--make-dummy',
						action='store_true',
						dest='makedummy',
						help="make a dummy file with ones x number of donors")

	parser.add_argument('--donors',
						type=int,
						dest='ndonors',
						help="Number of samples. Only for --make-dummy")

	opts = parser.parse_args()
	return opts

if __name__ == '__main__':
	opts = parse_args()
	if opts.makedummy and opts.ndonors != None:
		# cov_df = pd.read_table(opts.infile, sep="\t", header=0, index_col=0)
		# ndonors = cov_df.shape[1]
		# print("Found {:d} donors".format(ndonors))
		if opts.ndonors > 0:
			with open(opts.outfile, 'w') as outstream:
				outstream.writelines(["1.0\n"] * opts.ndonors)
		else:
			raise ValueError("Number of donors must be > 0.")
	else:
		cov_df = pd.read_table(opts.infile, sep="\t", header=0, index_col=0)
		cov_df.T.to_csv(opts.outfile + ".gnet", header=False, index=False, doublequote=False, sep="\t")