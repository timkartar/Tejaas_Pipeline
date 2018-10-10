import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
import sys, os   
from utils import utils
from collections import defaultdict
import plotter
from validation_sets import PlotSet, sets

outdir   = os.path.join("/cbscratch/franco/", "tejaas_output/validation")
inputdir = os.path.join("/cbscratch/franco/tejaas_output")

basename = "new_rr_withcis_zoom"

algorithm = "TEJAAS"
testset       = PlotSet(sets[algorithm]["cardio_mono"])
validationset = PlotSet(sets[algorithm]["cardio_macro"])

tejaas_data1_out_dir  = os.path.join(inputdir, testset.datadir)
tejaas_data2_out_dir  = os.path.join(inputdir, validationset.datadir)

chroms_file = "devtools/chroms.txt"

plot_outdir = os.path.join(outdir, "plots")
if not os.path.exists(plot_outdir):
    os.makedirs(plot_outdir)

with open(chroms_file) as instream:
    chr_list = [ int(line.strip()) for line in instream if len(line.strip()) > 0 ]

subplotgrid = 111
fsize       = (12,12)
fig         = plt.figure(figsize = fsize)
ax          = fig.add_subplot(subplotgrid)
plt.title("{:s} test, {:s} validation".format(testset.title, validationset.title))


pval_thres = 0.05
options = defaultdict(lambda: None)
method = "jpa-rr"
options['tejaas_method'] = "rr"
options['pval_thres']    = pval_thres
options['zoom']          = True
options['zoom_percent']  = 0.05
options["scale"] = False
replication = "empirical"
options['replication'] = replication

# chrs = "_".join([str(x) for x in chr_list])

for chrom in chr_list:
    print("Validating CHR {:d}".format(chrom))
    mychr_list = [chrom]
    outfile = os.path.join(plot_outdir, basename +"_chr{:s}_roc_".format(str(chrom)) +replication+"_"+testset.name+"_"+validationset.name+".png")

    random_lines_count = 0
    ##### Load Results #####
    betas = [0.001, 0.01] #0.006, 0.01, 0.02,  0.05, 0.1]
    for i in range(len(betas)):
        for j in range(i, len(betas)):
            beta  = betas[i]
            beta2 = betas[j] 
            print("Processing sigma_beta: {:f} vs {:f}".format(beta, beta2))

            dir1  = os.path.join(tejaas_data1_out_dir, method, "beta_{:n}".format(beta))
            dir2  = os.path.join(tejaas_data2_out_dir, method, "beta_{:n}".format(beta2))

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

            utils.get_replication_sizes(tejaas_dict1, tejaas_dict2, algorithm)

    # Now add Matrix eqtl

    algorithm = "MatrixEQTL"

    meqtl_test       = PlotSet(sets[algorithm]["cardio_mono"])
    meqtl_validation = PlotSet(sets[algorithm]["cardio_macro"])

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

    utils.get_replication_sizes(mqtl_dict1, mqtl_dict2, algorithm)

    ax.legend()
    ax.set_xlabel("# of SNPs")
    ax.set_ylabel("# of validated SNPs")
    print(outfile)
    plt.savefig(outfile)