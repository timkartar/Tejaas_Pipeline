import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
import sys, os   
from utils import utils
from collections import defaultdict
import plotter
from validation_sets import PlotSet, sets

SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 20

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

outdir   = os.path.join("/cbscratch/franco/tejaas_output/validation")
inputdir = os.path.join("/cbscratch/franco/tejaas_output")

basename = "newnew_rr_genomewide"

algorithm = "TEJAAS"
testset       = PlotSet(sets[algorithm]["cardio_mono_cismask"])
validationset = PlotSet(sets[algorithm]["cardio_macro_cismask"])


algorithm = "MatrixEQTL"
meqtl_test       = PlotSet(sets[algorithm]["cardio_mono"])
meqtl_validation = PlotSet(sets[algorithm]["cardio_macro"])

tejaas_data1_out_dir  = os.path.join(inputdir, testset.datadir)
tejaas_data2_out_dir  = os.path.join(inputdir, validationset.datadir)

chroms_file = "devtools/chroms_all.txt"

plot_outdir = os.path.join(outdir, "plots")
if not os.path.exists(plot_outdir):
    os.makedirs(plot_outdir)

with open(chroms_file) as instream:
    chr_list = [ int(line.strip()) for line in instream if len(line.strip()) > 0 ]

subplotgrid = 111
fsize       = (12,12)


pval_thres = 0.05
options = defaultdict(lambda: None)
method = "jpa-rr"
options['tejaas_method'] = "rr"
options['pval_thres']    = pval_thres
options['zoom']          = False
options['zoom_percent']  = 0.05
options["scale"] = True
replication = "empirical"
options['replication'] = replication
betas = [0.01] #0.006, 0.01, 0.02,  0.05, 0.1]


# for chrom in chr_list:

fig         = plt.figure(figsize = fsize)
ax          = fig.add_subplot(subplotgrid)
plt.title("{:s} {:s} test, {:s} {:s} validation".format(testset.title, testset.tissue_title, validationset.title, validationset.tissue_title))

algorithm = "TEJAAS"
mychr_list = chr_list
if options['zoom']:
    outfile = os.path.join(plot_outdir, basename +"_roc_" +replication+ ".zoom." + \
                           testset.cismask +"_"+testset.name+"_"+testset.tissue+ \
                           "__"+validationset.name+"_"+validationset.tissue+".png")
else:
    outfile = os.path.join(plot_outdir, basename +"_roc_" +replication+ "." + \
                           testset.cismask +"_"+testset.name+"_"+testset.tissue+ \
                           "__"+validationset.name+"_"+validationset.tissue+".png")

random_lines_count = 0
##### Load Results #####
for i in range(len(betas)):
    for j in range(i, len(betas)):
        beta  = betas[i]
        beta2 = betas[j] 
        print("Processing sigma_beta: {:f} vs {:f}".format(beta, beta2))

        dir1  = os.path.join(tejaas_data1_out_dir, method, "beta_{:n}".format(beta))
        dir2  = os.path.join(tejaas_data2_out_dir, method, "beta_{:n}".format(beta2))

        print(dir1)
        print(dir2)
        if not (os.path.exists(dir1) and os.path.exists(dir2)):
            continue

        if options['tejaas_method'] == "rr" or options['tejaas_method'] == "jpa-rr":
            tejaas_dict1 = utils.load_tejaas_results(mychr_list, dir1)
            tejaas_dict2 = utils.load_tejaas_results(mychr_list, dir2)

        if options['tejaas_method'] == "jpa":
            tejaas_dict1 = utils.load_tejaas_jpa_results(mychr_list, dir1)
            tejaas_dict2 = utils.load_tejaas_jpa_results(mychr_list, dir2)

        #############################

        result_dicts = defaultdict(None)
        result_dicts['tejaas_test']   = tejaas_dict1
        result_dicts['tejaas_valid']  = tejaas_dict2

        x_vals, i_vals, n_vals, recall, ppv = plotter.get_xy_roc_validation(result_dicts, algorithm, options)
        ax.plot(i_vals, recall, label=algorithm +"-"+ str(beta) +"-"+ str(beta2))

        if random_lines_count < 1: 
            rx_vals, ri_vals, rn_vals, rrecall, rppv = plotter.get_xy_roc_validation(result_dicts, algorithm, options, randomize=True)
            ax.plot(ri_vals, rrecall, label=algorithm + " random"+"-"+ str(beta) +"-"+ str(beta2))
            random_lines_count += 1 

        # utils.get_replication_sizes(tejaas_dict1, tejaas_dict2, algorithm)

# Now add Matrix eqtl

algorithm = "MatrixEQTL"

meqtl_file1 = os.path.join(inputdir, meqtl_test.datadir, meqtl_test.transfile)
meqtl_file2 = os.path.join(inputdir, meqtl_validation.datadir, meqtl_validation.transfile) 

mqtl_dict1 = utils.load_matrixeqtl_results(mychr_list, meqtl_file1)
mqtl_dict2 = utils.load_matrixeqtl_results(mychr_list, meqtl_file2)

result_dicts['mqtl_test']  = mqtl_dict1
result_dicts['mqtl_valid'] = mqtl_dict2

x_vals, i_vals, n_vals, recall, ppv = plotter.get_xy_roc_validation(result_dicts, algorithm, options)
rx_vals, ri_vals, rn_vals, rrecall, rppv = plotter.get_xy_roc_validation(result_dicts, algorithm, options, randomize=True)

ax.plot(i_vals, recall, label=algorithm)
ax.plot(ri_vals, rrecall, label=algorithm+"-random")

# utils.get_replication_sizes(mqtl_dict1, mqtl_dict2, algorithm)

ax.legend()
ax.set_xlabel("# of SNPs")
ax.set_ylabel("# of validated SNPs")
print(outfile)
plt.savefig(outfile)