import matplotlib.pyplot as plt
import numpy as np
import sys, os   
from utils import utils


def plot_raktim_validation( result_dicts, algorithms, options, outfile):
    
    # algorithms = ["TEJAAS","MatrixEQTL"]
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
    if n_algorithms == 4:
        subplotgrid = 221
        fsize       = (15,15)
        
    fig = plt.figure(figsize = fsize)
        
    for algo_i in range(n_algorithms):
        algorithm  = algorithms[algo_i]
        subplotgrid += algo_i

        if algorithm == "TEJAAS":
            tejaas_test_dict   = result_dicts["tejaas_test"]
            tejaas_valid_dict = result_dicts["tejaas_valid"]
            if options['tejaas_method'] == "rr" or options['tejaas_method'] == "jpa-rr":
                x_values, toplot = utils.get_validation_curve(tejaas_test_dict, tejaas_valid_dict, pval_thres=options['pval_thres'])
                x_values_random, random_toplot = utils.get_validation_curve(tejaas_test_dict, tejaas_valid_dict, pval_thres=options['pval_thres'], randomize=True)
            
            if options['tejaas_method'] == "jpa":
                x_values, toplot = utils.get_validation_curve_jpa(tejaas_test_dict, tejaas_valid_dict)
                x_values_random, random_toplot = utils.get_validation_curve_jpa(tejaas_test_dict, tejaas_valid_dict, randomize=True)
            
        if algorithm == "MatrixEQTL":
            mqtl_gtex_dict     = result_dicts["mqtl_gtex"]
            mqtl_cardio_dict   = result_dicts["mqtl_cardio"]
            x_values, toplot = utils.get_validation_curve(mqtl_gtex_dict, mqtl_cardio_dict, pval_thres=options['pval_thres'])
            x_values_random, random_toplot = utils.get_validation_curve(mqtl_gtex_dict, mqtl_cardio_dict, pval_thres=options['pval_thres'], randomize=True)
            
        ax = fig.add_subplot(subplotgrid)
        plt.title(' GTEx test, Cardio validation')
        if options['zoom']:
            zoom_limit = options['zoom_percent'] * len(x_values)
            x_values = x_values[:zoom_limit]
            to_plot = toplot[:zoom_limit]
            x_values_random = x_values_random[:zoom_limit]
            random_to_plot = random_toplot[:zoom_limit]
        ax.plot(x_values, toplot, label=algorithm)
        ax.plot(x_values_random, random_toplot, label='Random SNPs - '+algorithm)
        ax.legend()
        ax.set_xlabel("# of SNPs")
        ax.set_ylabel("# of validated SNPs")


    # fig_file = os.path.join(outdir, "chr"+str(chr_list[0])+"_trans_comparison.png")
    print(outfile)
    plt.savefig(outfile)
    # plt.show()


def plot_rank_2dhistogram(result_dicts, algorithms, outfile):

    # algorithms = ["TEJAAS","MatrixEQTL"]
    n_algorithms = len(algorithms)
    if n_algorithms == 1:
        subplotgrid = 111
        fsize       = (10,8)
    if n_algorithms == 2:
        subplotgrid = 121
        fsize       = (18,8)
    if n_algorithms == 3:
        subplotgrid = 131
        fsize       = (26,8)
    if n_algorithms == 4:
        subplotgrid = 221
        fsize       = (18,16)

    fig = plt.figure(figsize = fsize)
        
    for algo_i in range(n_algorithms):
        algorithm  = algorithms[algo_i]
        subplotgrid += algo_i

        if algorithm == "TEJAAS":
            tejaas_test_dict   = result_dicts["tejaas_test"]
            tejaas_valid_dict = result_dicts["tejaas_valid"]
            dict1 = tejaas_test_dict
            dict2 = tejaas_valid_dict

        if algorithm == "MatrixEQTL":
            mqtl_gtex_dict     = result_dicts["mqtl_gtex"]
            mqtl_cardio_dict   = result_dicts["mqtl_cardio"]
            dict1 = mqtl_gtex_dict
            dict2 = mqtl_cardio_dict

        cdict1, cdict2 = utils.get_compatible_snp_dicts(dict1, dict2)

        test_rank = utils.get_ranks(cdict1, descending=False)
        validation_rank = utils.get_ranks(cdict2, descending=False)
        
        allrsids = list(test_rank.keys())
        xpoints = [test_rank[k] for k in allrsids]
        ypoints = [validation_rank[k] for k in allrsids]

        ax = fig.add_subplot(subplotgrid)
        plt.title('Rank 2DHist - '+algorithm)
        # ax.scatter(xpoints, ypoints, alpha=0.1)
        pos = ax.hist2d(xpoints, ypoints, bins=50)
        ax.set_xlabel("Rank GTEx")
        ax.set_ylabel("Rank Cardiogenics")
        fig.colorbar(pos[3], ax = ax)

    print(outfile)
    plt.savefig(outfile)
    # plt.show()


def get_xy_roc_validation(result_dicts, algorithm, options, randomize=False):
        
    if algorithm == "TEJAAS":
        tejaas_test_dict   = result_dicts["tejaas_test"]
        tejaas_valid_dict  = result_dicts["tejaas_valid"]
        test_dict, validation_dict = utils.get_compatible_snp_dicts(tejaas_test_dict, tejaas_valid_dict)
        if options['tejaas_method'] == "rr" or options['tejaas_method'] == "jpa-rr":
            if options['replication'] == "empirical":
                x_vals, i_vals, n_vals, recall, ppv = utils.evaluate_replication_empirical(test_dict, validation_dict, randomize=randomize)
            else:
                x_vals, i_vals, n_vals, recall, ppv = utils.evaluate_replication(test_dict, validation_dict, randomize=randomize)
        if options['tejaas_method'] == "jpa":
            x_vals, i_vals, n_vals, recall, ppv = utils.evaluate_replication_jpa(test_dict, validation_dict, randomize=randomize)
    if algorithm == "MatrixEQTL":
        mqtl_gtex_dict     = result_dicts["mqtl_test"]
        mqtl_cardio_dict   = result_dicts["mqtl_valid"]
        validation_dict, test_dict = utils.get_compatible_snp_dicts(mqtl_gtex_dict, mqtl_cardio_dict)
        if options['replication'] == "empirical":
            x_vals, i_vals, n_vals, recall, ppv = utils.evaluate_replication_empirical(test_dict, validation_dict, randomize=randomize)
        else:
            x_vals, i_vals, n_vals, recall, ppv = utils.evaluate_replication(test_dict, validation_dict, randomize=randomize)
    
    return x_vals, np.array(i_vals), n_vals, recall, ppv


def plot_roc_validation(result_dicts, algorithms, options, outfile):

    # algorithms = ["TEJAAS" ,"MatrixEQTL"]
    n_algorithms = len(algorithms)
    if n_algorithms == 1:
        subplotgrid = 121
        fsize       = (12,6)
    if n_algorithms == 2:
        subplotgrid = 221
        fsize       = (12,12)
    if n_algorithms == 3:
        subplotgrid = 321
        fsize       = (12,18)
        
    fig = plt.figure(figsize = fsize)
       
    for algo_i in range(n_algorithms):
        algorithm  = algorithms[algo_i]

        if algorithm == "TEJAAS":
            tejaas_test_dict   = result_dicts["tejaas_test"]
            tejaas_valid_dict  = result_dicts["tejaas_valid"]

            test_dict, validation_dict = utils.get_compatible_snp_dicts(tejaas_test_dict, tejaas_valid_dict)
            
            if options['tejaas_method'] == "rr" or options['tejaas_method'] == "jpa-rr":
                if options['replication'] == "empirical":
                    x_vals, i_vals, n_vals, recall, ppv = utils.evaluate_replication_empirical(test_dict, validation_dict)
                    rx_vals, ri_vals, rn_vals, rrecall, rppv = utils.evaluate_replication_empirical(test_dict, validation_dict, randomize=True)
                else:
                    x_vals, i_vals, n_vals, recall, ppv = utils.evaluate_replication(test_dict, validation_dict)
                    rx_vals, ri_vals, rn_vals, rrecall, rppv = utils.evaluate_replication(test_dict, validation_dict, randomize=True)
            if options['tejaas_method'] == "jpa":
                x_vals, i_vals, n_vals, recall, ppv = utils.evaluate_replication_jpa(test_dict, validation_dict)
                rx_vals, ri_vals, rn_vals, rrecall, rppv = utils.evaluate_replication_jpa(test_dict, validation_dict, randomize=True)
            
        if algorithm == "MatrixEQTL":
            mqtl_gtex_dict     = result_dicts["mqtl_gtex"]
            mqtl_cardio_dict   = result_dicts["mqtl_cardio"]
            validation_dict, test_dict = utils.get_compatible_snp_dicts(mqtl_gtex_dict, mqtl_cardio_dict)
            if options['replication'] == "empirical":
                x_vals, i_vals, n_vals, recall, ppv = utils.evaluate_replication_empirical(test_dict, validation_dict)
                rx_vals, ri_vals, rn_vals, rrecall, rppv = utils.evaluate_replication_empirical(test_dict, validation_dict, randomize=True)   
            else:
                x_vals, i_vals, n_vals, recall, ppv = utils.evaluate_replication(test_dict, validation_dict)
                rx_vals, ri_vals, rn_vals, rrecall, rppv = utils.evaluate_replication(test_dict, validation_dict, randomize=True)       
        
        if options['zoom']:
            zoom_limit = options['zoom_percent'] * len(i_vals)
            i_vals = i_vals[:zoom_limit]
            ri_vals = ri_vals[:zoom_limit]
            recall = recall[:zoom_limit]
            rrecall = rrecall[:zoom_limit]
            ppv = ppv[:zoom_limit]
            rppv = rppv[:zoom_limit]

        ax = fig.add_subplot(subplotgrid)
        plt.title('Recall')
        ax.plot(i_vals, recall, label=algorithm)
        ax.plot(ri_vals, rrecall, label=algorithm + " random")
        ax.legend()
        ax.set_xlabel("# SNPs")
        ax.set_ylabel("Recall")

        subplotgrid += 1
        ax = fig.add_subplot(subplotgrid)
        plt.title('PPV')
        ax.plot(i_vals, ppv, label=algorithm)
        ax.plot(ri_vals, rppv, label=algorithm + " random")
        ax.set_ylim(0,1.1)
        ax.legend()
        ax.set_xlabel("# SNPs")
        ax.set_ylabel("PPV")
        
        subplotgrid += 1

    print(outfile)
    plt.savefig(outfile, dpi=300)
    # plt.show()