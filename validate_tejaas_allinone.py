import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
import sys, os   
from utils import utils
from collections import defaultdict
import plotter

outdir            = os.path.join("/cbscratch/franco/", "tejaas_output/validation")
tejaas_outdir     = os.path.join("/cbscratch/franco/", "tejaas_output/output")

basename = "alt_cismask_allinone_scan_jpa_"

test = "cardio"
test_basename = "cardiogenics"
test_title = "Cardio Monocytes"
validation = "cardio"
validation_basename = "cardiogenics"
validation_title = "Cardio Macrophages"

# validation = "gtex"
# validation_basename = "gtex"
# validation_title = "GTEx"

# tejaas_data1_out_dir  = test+"_results_cismask_mpi"
# tejaas_data2_out_dir  = test+"_results_cismask_mpi_macro"
# tejaas_data1_out_dir  = test + "_results_cismask_mpi"
# tejaas_data2_out_dir  = validation + "_results_cismask_mpi"
tejaas_data1_out_dir  = test+"_results_mono_us"
tejaas_data2_out_dir  = test+"_results_macro_us"



MatrixEQTL_outdir = os.path.join("/cbscratch/franco/", "tejaas_output/matrixEQTL_out_filteredGT_gencode")
meqtl_dir1  = "cardio_results_mono_us"
meqtl_dir2  = "cardio_results_macro_us"
# meqtl_dir1  = test +"_results"
# meqtl_dir2  = validation +"_results"
meqtl_file1 = os.path.join(MatrixEQTL_outdir, meqtl_dir1, test_basename+"_MatrixEQTL_chr{:d}.transout")
meqtl_file2 = os.path.join(MatrixEQTL_outdir, meqtl_dir2, validation_basename+"_MatrixEQTL_chr{:d}.transout") 

chroms_file = "devtools/chroms.txt"

if not os.path.exists(outdir):
    os.makedirs(outdir)

with open(chroms_file) as instream:
    chr_list = [ int(line.strip()) for line in instream if len(line.strip()) > 0 ]      


subplotgrid = 111
fsize       = (12,12)
fig = plt.figure(figsize = fsize)
ax = fig.add_subplot(subplotgrid)
# plt.title(' Cardio Mono test, Gtex validation')
plt.title("{:s} test, {:s} validation".format(test_title, validation_title))


algorithms = ["TEJAAS"]
algorithm = "TEJAAS"


plot_outdir = os.path.join(outdir, "plots")

pval_thres = 0.05
options = defaultdict(lambda: None)
method = "jpa-rr"
options['tejaas_method'] = "jpa"
options['pval_thres']    = pval_thres
# options['zoom']          = False
# options['zoom_percent']  = 0.1
scale = True

options['replication'] = "empirical"
if options['replication'] == "empirical":
    outfilemod = "empirical_"
else:
    outfilemod = "pval_thres_"

if not os.path.exists(plot_outdir):
    os.makedirs(plot_outdir)

outfile = os.path.join(plot_outdir, basename + "chr{:d}_roc_".format(chr_list[0]) +outfilemod+"validation_"+test+"_"+validation+".png")

random_lines_count = 0
##### Load Results #####
betas = [0.001, 0.006, 0.01, 0.02,  0.05, 0.1]
for i in range(len(betas)):
    for j in range(i, len(betas)):
        beta  = betas[i]
        beta2 = betas[j] 
        print("Processing sigma_beta: {:f} vs {:f}".format(beta, beta2))

        dir1  = os.path.join(tejaas_outdir, tejaas_data1_out_dir, method, "beta_{:n}".format(beta))
        dir2  = os.path.join(tejaas_outdir, tejaas_data2_out_dir, method, "beta_{:n}".format(beta2))

        if not (os.path.exists(dir1) and os.path.exists(dir2)):
            print(dir1)
            print(dir2)
            continue

        if options['tejaas_method'] == "rr" or options['tejaas_method'] == "jpa-rr":
            tejaas_dict1 = utils.load_tejaas_results(chr_list, dir1)
            tejaas_dict2 = utils.load_tejaas_results(chr_list, dir2)

        if options['tejaas_method'] == "jpa":
            tejaas_dict1 = utils.load_tejaas_jpa_results(chr_list, dir1)
            tejaas_dict2 = utils.load_tejaas_jpa_results(chr_list, dir2)

        #############################

        result_dicts = defaultdict(None)
        result_dicts['tejaas_test']   = tejaas_dict1
        result_dicts['tejaas_valid']  = tejaas_dict2

        x_vals, i_vals, n_vals, recall, ppv = plotter.get_xy_roc_validation(result_dicts, algorithm, options)
        if scale:
            i_vals = i_vals/max(i_vals)
        ax.plot(i_vals, recall, label=algorithm +"-"+ str(beta) +"-"+ str(beta2))

        if random_lines_count < 1: 
            rx_vals, ri_vals, rn_vals, rrecall, rppv = plotter.get_xy_roc_validation(result_dicts, algorithm, options, randomize=True)
            if scale:
                ri_vals = ri_vals/max(ri_vals)
            ax.plot(ri_vals, rrecall, label=algorithm + " random"+"-"+ str(beta) +"-"+ str(beta2))
            random_lines_count += 1 

        utils.get_replication_sizes(tejaas_dict1, tejaas_dict2, algorithm)

# Now add Matrix eqtl

algorithms = ["MatrixEQTL"]
algorithm = "MatrixEQTL"

mqtl_dict1 = utils.load_matrixeqtl_results(chr_list, meqtl_file1)
mqtl_dict2 = utils.load_matrixeqtl_results(chr_list, meqtl_file2)

result_dicts['mqtl_test']  = mqtl_dict1
result_dicts['mqtl_valid'] = mqtl_dict2

x_vals, i_vals, n_vals, recall, ppv = plotter.get_xy_roc_validation(result_dicts, algorithm, options)
rx_vals, ri_vals, rn_vals, rrecall, rppv = plotter.get_xy_roc_validation(result_dicts, algorithm, options, randomize=True)
if scale:
    i_vals = i_vals/max(i_vals)
    ri_vals = ri_vals/max(ri_vals)

ax.plot(i_vals, recall, label=algorithm)
ax.plot(ri_vals, rrecall, label=algorithm+"-random")

utils.get_replication_sizes(mqtl_dict1, mqtl_dict2, algorithm)


ax.legend()
ax.set_xlabel("# of SNPs")
ax.set_ylabel("# of validated SNPs")
print(outfile)
plt.savefig(outfile)