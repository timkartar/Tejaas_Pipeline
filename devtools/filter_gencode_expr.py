import sys, os
sys.path.append("../")
from iotools import readgtf
from iotools.readOxford import ReadOxford
from collections import defaultdict
import pandas as pd

GXFILE = sys.argv[1]

print("WARNING! only GTEx expressions! use gene_trim_dict for Cardiogenics!")
print("Gene Expr File: {:s}".format(GXFILE))

def filter_rows(df, genedict):
    # genedict = gene_dict
    # df = gx_df
    gx_gene_list = df.index
    common  = [genedict[x] for x in gx_gene_list]
    print("{:d} genes remained from {:d}".format(sum(common), len(gx_gene_list)))
    return gx_df[common]


gtfpath = "/cbscratch/franco/datasets/gtex/gencode.v19.annotation.gtf.gz"
# can't read this with current library
# gtfpath = "/cbscratch/franco/datasets/gtex/gencode.v28lift37.annotation.gtf.gz"
gene_info = readgtf.gencode_v12(gtfpath, trim=False)

gene_dict = defaultdict(lambda: False)
for g in gene_info:
    gene_dict[g.ensembl_id] = True
    
gene_info = readgtf.gencode_v12(gtfpath, trim=True)

gene_trim_dict = defaultdict(lambda: False)
for g in gene_info:
    gene_trim_dict[g.ensembl_id] = True 

gx_df = pd.read_table(GXFILE, sep="\t", header=0, index_col=0)
new_gx_df = filter_rows(gx_df, gene_dict)
new_gx_df.to_csv(GXFILE+".gencode_filter", doublequote=False, sep="\t")