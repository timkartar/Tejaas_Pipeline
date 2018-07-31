import matplotlib.pyplot as plt
import mpl_stylesheet
plt.switch_backend('agg')
import numpy as np
import sys, os   
import utils

# mpl_stylesheet.banskt_presentation(fontfamily = 'system')

# gtex_out_dir = os.environ.get("GTEX_OUTDIRBASE") #"../gtex_results"
# cardio_out_dir = os.environ.get("CARDIO_OUTDIRBASE") #"../cardio_results"

outdir = os.path.join(os.environ.get("HOME"), "Tejaas_Pipeline/output")
gtex_out_dir = os.path.join(os.environ.get("HOME"), "Tejaas_Pipeline/output/gtex_results")
cardio_out_dir = os.path.join(os.environ.get("HOME"),"Tejaas_Pipeline/output/cardio_results")
chroms_file = "devtools/chroms.txt"


algorithms = ["TEJAAS", "MatrixEQTL"]
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
offset = 10000

with open(chroms_file) as instream:
    chr_list = [int(line.strip())for line in instream if len(line.strip()) > 0]

for algo_i in range(n_algorithms):
    algorithm  = algorithms[algo_i]
    subplotgrid += algo_i

    if algorithm == "TEJAAS":
        gtex_beta   = 0.001
        cardio_beta = 0.001

        gtex_dir    = os.path.join(gtex_out_dir, "beta_{:n}".format(gtex_beta))
        cardio_dir  = os.path.join(cardio_out_dir, "beta_{:n}".format(cardio_beta))
        gtex_dict, cardio_dict = utils.load_tejaas_results(chr_list, gtex_dir, cardio_dir)
        
        x_values, toplot = utils.get_validation_curve(gtex_dict, cardio_dict, pval_thres=0.05, islog=True)
        x_values_random, random_toplot = utils.get_validation_curve(gtex_dict, cardio_dict, pval_thres=0.05, islog=True, randomize=True)

    if algorithm == "MatrixEQTL":
        gtex_dir   = gtex_out_dir
        cardio_dir = cardio_out_dir 
        gtex_dict, cardio_dict = utils.load_matrixeqtl_results(chr_list, gtex_dir, cardio_dir)

        x_values, toplot = utils.get_validation_curve(cardio_dict, gtex_dict, pval_thres=0.05, islog=False)
        x_values_random, random_toplot = utils.get_validation_curve(cardio_dict, gtex_dict, pval_thres=0.05, islog=False, randomize=True)

    ax = fig.add_subplot(subplotgrid)
    plt.title('Predicted on Cardiogenics, validated on GTEx')
    ax.plot(x_values[k:m],toplot[k:m], label=algorithm)
    ax.plot(x_values_random[k:m+offset], random_toplot[k:m+offset], label='Random SNPs - '+algorithm)
    #plt.plot(check_shuf_x[k:m+offset], shuf_toplot[k:m+offset], label='TEJAAS on Shuffled genotype')
    # ax.tight_layout()
    ax.legend()
    ax.set_xlabel("# of SNPs")
    ax.set_ylabel("# of validated SNPs")

fig_file = os.path.join(outdir, "cardio_gtex_trans_comparison.png")
print(fig_file)
plt.savefig(fig_file)
# plt.show()