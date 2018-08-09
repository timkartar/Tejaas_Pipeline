import mpmath
import numpy as np
import os 
import random
import operator

# decimal places
mpmath.mp.dps = 500
def pval(x): return mpmath.log10(1 - 0.5 * (1 + mpmath.erf(x/mpmath.sqrt(2))))

def load_tejaas_results(chr_list, gtex_dir, cardio_dir):
    gtex_dict = dict()
    cardio_dict = dict()
    for chrm in chr_list:
        print( "Reading TEJAAS chr{:d}".format(chrm))
        dirc = os.path.join(gtex_dir, "chr{:d}".format(chrm))
        pths = [os.path.join(dirc,path) for path in os.listdir(dirc) if "_rr.txt" in path]
        for pth in pths:
            pvals= list()
            l = open(str(pth),"r").readlines()
            for line in l[1:]:    
                rsid  = line.split("\t")[0]
                P     = float(line.split("\t")[5])
                Q     = float(line.split("\t")[2])
                Mu    = float(line.split("\t")[3])
                Sigma = float(line.split("\t")[4])
                pvalue= np.log10(P) if P!=0 else pval((Q-Mu)/Sigma)
                gtex_dict[rsid] = pvalue              

        dirc = os.path.join(cardio_dir, "chr{:d}".format(chrm))
        pths = [os.path.join(dirc,path) for path in os.listdir(dirc) if "_rr.txt" in path]
        for pth in pths:
            pvals = list()
            l = open(str(pth),"r").readlines()
            for line in l[1:]:
                rsid  = line.split("\t")[0]
                P     = float(line.split("\t")[5])
                Q     = float(line.split("\t")[2])
                Mu    = float(line.split("\t")[3])
                Sigma = float(line.split("\t")[4])
                pvalue= np.log10(P) if P!=0 else pval((Q-Mu)/Sigma)
                cardio_dict[rsid] = pvalue
    return gtex_dict, cardio_dict

def load_tejaas_results_orig(chr_list, gtex_dir, cardio_dir):
    gtex_dict = dict()
    cardio_dict = dict()
    for chrm in chr_list:
        print( "Reading TEJAAS chr{:d}".format(chrm))
        dirc = os.path.join(gtex_dir, "chr{:d}".format(chrm))
        pths = [os.path.join(dirc,path) for path in os.listdir(dirc) if "_rr.txt" in path]
        for pth in pths:
            pvals= list()
            l = open(str(pth),"r").readlines()
            rsids = [line.split("\t")[0] for line in l[1:]]
            for line in l[1:]:    
                P     = float(line.split("\t")[5])
                Q     = float(line.split("\t")[2])
                Mu    = float(line.split("\t")[3])
                Sigma = float(line.split("\t")[4])
                pvals.append(np.log10(P) if P!=0 else pval((Q-Mu)/Sigma))
            for i in range(len(rsids)):
                gtex_dict[rsids[i]] = pvals[i]

        dirc = os.path.join(cardio_dir, "chr{:d}".format(chrm))
        pths = [os.path.join(dirc,path) for path in os.listdir(dirc) if "_rr.txt" in path]
        for pth in pths:
            pvals = list()
            l = open(str(pth),"r").readlines()
            rsids = [line.split("\t")[0] for line in l[1:]]
            for line in l[1:]:
                P     = float(line.split("\t")[5])
                Q     = float(line.split("\t")[2])
                Mu    = float(line.split("\t")[3])
                Sigma = float(line.split("\t")[4])
                pvals.append(np.log10(P) if P!=0 else pval((Q-Mu)/Sigma))
            for i in range(len(rsids)):
                cardio_dict[rsids[i]] = pvals[i]
    return gtex_dict, cardio_dict

def load_matrixeqtl_results(chr_list, gtex_dir, cardio_dir):
    gtex_matrix_dict = {}
    cardio_matrix_dict = {}
    for chrm in chr_list:
        print( "Reading matrixEQTL chr{:d}".format(chrm))
        results_file = os.path.join(gtex_dir, "gtex_MatrixEQTL_chr{:d}.transout".format(chrm))
        with open(results_file) as instream:
            _ = instream.readline()
            for line in instream:
                arr  = line.rstrip().split("\t")
                rsid = arr[0]
                FDR  = float(arr[5])
                if rsid not in gtex_matrix_dict:
                    gtex_matrix_dict[rsid] = np.log10(FDR)
        
        results_file = os.path.join(cardio_dir, "cardiogenics_MatrixEQTL_chr{:d}.transout".format(chrm))
        with open(results_file) as instream:
            _ = instream.readline()
            for line in instream:
                arr  = line.rstrip().split("\t")
                rsid = arr[0]
                FDR  = float(arr[5])
                if rsid not in cardio_matrix_dict:
                    cardio_matrix_dict[rsid] = np.log10(FDR)
    return gtex_matrix_dict, cardio_matrix_dict

def get_validation_curve(test_dict, validation_dict, pval_thres=0.05, islog=True, randomize=False):    
    n_snps = len(validation_dict.keys())
    isReversed = False
    if islog:
        isReversed = True
    validation_tuples = sorted(validation_dict.items(), key=operator.itemgetter(1),reverse=isReversed)
    validation_rsids = [validation_tuples[i][0] for i in range(n_snps)]
    validation_pvals = [validation_tuples[i][1] for i in range(n_snps)]

    # test_pvals = [test_dict[e] for e in test_dict.keys()]
    
    # print(np.min(validation_pvals), np.max(validation_pvals))
    # print(np.min(test_pvals), np.max(test_pvals))
    if randomize:
        random.shuffle(validation_rsids)
    
    toplot = []
    check_x = []
    positives = []
    if islog:
        pval_thres = np.log10(pval_thres)

    # print("log pval thres is: {:f}".format(pval_thres))
    # print("Nº of validation tuples: {:d}".format(len(validation_tuples)))
        
    i = 0 
    while(i < n_snps):
        try:
            if(test_dict[validation_rsids[i]] < pval_thres):
                positives.append(i)
        except:
            pass
        i = i + 1
        while(i < n_snps and validation_pvals[i] == validation_pvals[i-1]):
            try:
                if(test_dict[validation_rsids[i]] < pval_thres):
                    positives.append(i)
            except:
                pass
            i = i + 1
        check_x.append(i)
        toplot.append(len(positives))
    return check_x, toplot