import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
import sys, os   
from utils import utils
from collections import defaultdict
import plotter

outdir            = os.path.join("/cbscratch/franco/", "tejaas_output/validation")
tejaas_outdir     = os.path.join("/cbscratch/franco/", "tejaas_output/output")
tejaas_outdir_rakt= os.path.join("/cbscratch/franco/", "raktim_data/crop_tejaas_out")

# MatrixEQTL_outdir = os.path.join("/cbscratch/franco/", "tejaas_output/matrixEQTL_out")
# MatrixEQTL_outdir = os.path.join("/cbscratch/franco/", "tejaas_output/matrixEQTL_out_filteredGT")
MatrixEQTL_outdir = os.path.join("/cbscratch/franco/", "tejaas_output/matrixEQTL_out_filteredGT_gencode")

basename = "forcecis_cismask_mpi_gencode_"

mod="_forcecis_cismask_mpi"
# mod="_forcetrans"
# mod="_forcecis"
tejaas_gtex_out_dir   = "gtex_results"+mod
tejaas_cardio_out_dir = "cardio_results"+mod
gtex_out_dir          = "gtex_results"
cardio_out_dir        = "cardio_results"
chroms_file           = "devtools/chroms.txt"

if not os.path.exists(outdir):
    os.makedirs(outdir)

with open(chroms_file) as instream:
    chr_list = [ int(line.strip()) for line in instream if len(line.strip()) > 0 ]      


##### Load Results #####
gtex_beta   = 0.02 # 0.001  # 0.02 #
cardio_beta = 0.006 # 0.001  # 0.006 #

method = "rr"
if method == "rr":
    gtex_dir    = os.path.join(tejaas_outdir, tejaas_gtex_out_dir, method,  "beta_{:n}".format(gtex_beta))
    cardio_dir  = os.path.join(tejaas_outdir, tejaas_cardio_out_dir, method, "beta_{:n}".format(cardio_beta))
    tejaas_gtex_dict = utils.load_tejaas_results(chr_list, gtex_dir)
    tejaas_cardio_dict = utils.load_tejaas_results(chr_list, cardio_dir)

if method == "jpa":
    gtex_dir    = os.path.join(tejaas_outdir, tejaas_gtex_out_dir, method,  "beta_{:n}".format(gtex_beta))
    cardio_dir  = os.path.join(tejaas_outdir, tejaas_cardio_out_dir, method, "beta_{:n}".format(cardio_beta))
    tejaas_gtex_dict = utils.load_tejaas_jpa_results(chr_list, gtex_dir)
    tejaas_cardio_dict = utils.load_tejaas_jpa_results(chr_list, cardio_dir)

gtex_file   = os.path.join(MatrixEQTL_outdir, gtex_out_dir, "gtex_MatrixEQTL_chr{:d}.transout")
cardio_file = os.path.join(MatrixEQTL_outdir, cardio_out_dir, "cardiogenics_MatrixEQTL_chr{:d}.transout") 

# gtex_file = gtex_file + ".truetrans.afterfix"
# cardio_file = cardio_file + ".truetrans.afterfix"

mqtl_gtex_dict = utils.load_matrixeqtl_results(chr_list, gtex_file)
mqtl_cardio_dict = utils.load_matrixeqtl_results(chr_list, cardio_file)


#############################

algorithms = ["TEJAAS"]

result_dicts = defaultdict(None)
result_dicts['tejaas_test']   = tejaas_gtex_dict
result_dicts['tejaas_valid'] = tejaas_cardio_dict
result_dicts['mqtl_gtex']     = mqtl_gtex_dict
result_dicts['mqtl_cardio']   = mqtl_cardio_dict

pval_thres = 0.05
options = defaultdict(lambda: None)
options['tejaas_method'] = method
options['pval_thres']    = pval_thres
options['zoom']          = False
options['zoom_percent']  = 0.1

plot_outdir = os.path.join(outdir, "plots")

if not os.path.exists(plot_outdir):
    os.makedirs(plot_outdir)

print("Creating Raktim's validation plot")
outfile = os.path.join(plot_outdir, basename + "chr{:d}_raktim_validation_{:n}_{:n}.png".format(chr_list[0], gtex_beta, cardio_beta ))
plotter.plot_raktim_validation( result_dicts, algorithms, options, outfile)

print("Creating 2D rank histogram")
outfile = os.path.join(plot_outdir, basename + "chr{:d}_2d_hist_{:n}_{:n}.png".format(chr_list[0], gtex_beta, cardio_beta ))
plotter.plot_rank_2dhistogram(result_dicts, algorithms, outfile)

print("Creating ROC validation plots")
outfile = os.path.join(plot_outdir, basename + "chr{:d}_roc_validation_{:n}_{:n}.png".format(chr_list[0], gtex_beta, cardio_beta ))
plotter.plot_roc_validation(result_dicts, algorithms, options, outfile)

print("Creating Empirical ROC validation plots")
options['replication'] = "empirical"
outfile = os.path.join(plot_outdir, basename + "chr{:d}_roc_empiricial_validation_{:n}_{:n}.png".format(chr_list[0], gtex_beta, cardio_beta ))
plotter.plot_roc_validation(result_dicts, algorithms, options, outfile)

algorithms = ["TEJAAS","MatrixEQTL"]
for algorithm in algorithms:
    if algorithm == "TEJAAS":
        utils.get_replication_sizes(tejaas_gtex_dict, tejaas_cardio_dict, algorithm)
        
    if algorithm == "MatrixEQTL":
        utils.get_replication_sizes(mqtl_gtex_dict, mqtl_cardio_dict, algorithm)