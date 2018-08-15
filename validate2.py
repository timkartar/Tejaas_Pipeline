import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
import sys, os   
from utils import utils

outdir            = os.path.join("/cbscratch/franco/", "tejaas_output/validation")
tejaas_outdir     = os.path.join("/cbscratch/franco/", "tejaas_output/output")
tejaas_outdir_rakt= os.path.join("/cbscratch/franco/", "raktim_data/crop_tejaas_out")

# MatrixEQTL_outdir = os.path.join("/cbscratch/franco/", "tejaas_output/matrixEQTL_out")
MatrixEQTL_outdir = os.path.join("/cbscratch/franco/", "tejaas_output/matrixEQTL_out_filteredGT")
gtex_out_dir      = "gtex_results"
cardio_out_dir    = "cardio_results"


chroms_file = "chroms.txt"
with open(chroms_file) as instream:
    chr_list = [int(line.strip())for line in instream if len(line.strip()) > 0]      


# Load all results
gtex_beta   = 0.02  # 0.001 #
cardio_beta = 0.006 # 0.001 #

gtex_dir    = os.path.join(tejaas_outdir, gtex_out_dir, "beta_{:n}".format(gtex_beta))
cardio_dir  = os.path.join(tejaas_outdir, cardio_out_dir, "beta_{:n}".format(cardio_beta))
# tejaas_gtex_dict, tejaas_cardio_dict = load_tejaas_results(chr_list, gtex_dir, cardio_dir)
tejaas_gtex_dict = load_tejaas_results(chr_list, gtex_dir)
tejaas_cardio_dict = load_tejaas_results(chr_list, cardio_dir)

# gtex_dir   = os.path.join(MatrixEQTL_outdir, gtex_out_dir)
# cardio_dir = os.path.join(MatrixEQTL_outdir, cardio_out_dir) 
# mqtl_gtex_dict, mqtl_cardio_dict = load_matrixeqtl_results(chr_list, gtex_dir, cardio_dir)
gtex_file   = os.path.join(MatrixEQTL_outdir, gtex_out_dir, "gtex_MatrixEQTL_chr{:d}.transout")
cardio_file = os.path.join(MatrixEQTL_outdir, cardio_out_dir, "cardiogenics_MatrixEQTL_chr{:d}.transout") 
mqtl_gtex_dict = load_matrixeqtl_results(chr_list, gtex_file)
mqtl_cardio_dict = load_matrixeqtl_results(chr_list, cardio_file)



algorithms = ["TEJAAS","MatrixEQTL"]
n_algorithms = len(algorithms)
if n_algorithms == 1:
    subplotgrid = 111
    fsize       = (6,6)
if n_algorithms == 2:
    subplotgrid = 121
    fsize       = (12,6)
if n_algorithms == 3:
    subplotgrid = 131
    fsize       = (18,6)
if n_algorithms == 3:
    subplotgrid = 221
    fsize       = (15,15)
    
fig = plt.figure(figsize = fsize)
k = 0
m = 10000
offset = 1000
    
for algo_i in range(n_algorithms):
    algorithm  = algorithms[algo_i]
    subplotgrid += algo_i

    if algorithm == "TEJAAS":       
        x_values, toplot = utils.get_validation_curve(tejaas_gtex_dict, tejaas_cardio_dict, pval_thres=0.05, islog=True)
        x_values_random, random_toplot = utils.get_validation_curve(tejaas_gtex_dict, tejaas_cardio_dict, pval_thres=0.05, islog=True, randomize=True)
        
    if algorithm == "MatrixEQTL":
        x_values, toplot = utils.get_validation_curve(mqtl_gtex_dict, mqtl_cardio_dict, pval_thres=0.05, islog=False)
        x_values_random, random_toplot = utils.get_validation_curve(mqtl_gtex_dict, mqtl_cardio_dict, pval_thres=0.05, islog=False, randomize=True)
        
    ax = fig.add_subplot(subplotgrid)
    plt.title('Predicted on Cardiogenics, validated on GTEx')
    ax.plot(x_values[k:m],toplot[k:m], label=algorithm)
    ax.plot(x_values_random[k:m+offset], random_toplot[k:m+offset], label='Random SNPs - '+algorithm)
    ax.legend()
    ax.set_xlabel("# of SNPs")
    ax.set_ylabel("# of validated SNPs")

if not os.path.exists(outdir):
    os.makedirs(outdir)
fig_file = os.path.join(outdir, "chr{:d}_cardio_gtex_trans_comparison_{:n}_{:n}.png".format(chr_list[0], gtex_beta, cardio_beta ))
print(fig_file)
plt.savefig(fig_file)
# plt.show()