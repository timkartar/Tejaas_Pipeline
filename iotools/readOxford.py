#!/usr/bin/env python


import sys
import os
import gzip
import collections
import numpy as np
import re
import pdb
import scipy.stats as ss

from utils.containers import SnpInfo

SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}


class Error(Exception):
    pass

class SAMPLE_NUMBER_ERROR(Error):
    pass

class GT_FREQS_NUMBER_ERROR(Error):
    pass

class ReadOxford:

    _read_samples_once = False
    _read_genotype_once = False
    _nloci = 0
    _nsample = 0

    
    # deafult is GTEx:
    #     - isdosage = True
    #     - data_columns = 6

    def __init__(self, gtfile, samplefile, startsnp=0, endsnp=1e15, isdosage=True, data_columns=6, chrom=None):
        if isinstance(chrom, int) and chrom > 0 and chrom < 23:
            self._chrom = chrom
        else:
            self._chrom = None
        self._gtfile = gtfile
        self._samplefile = samplefile
        self._startsnp = startsnp
        self._endsnp = endsnp
        self._data_columns = data_columns
        self._isdosage = isdosage
        self._read_genotypes()

        
    @property
    def nsample(self):
        self._read_samples()
        return self._nsample

    @property
    def samplenames(self):
        self._read_samples()
        return self._samplenames
    
    @property
    def nloci(self):
        return self._nloci

    @property
    def snpinfo(self):
        self._read_genotypes()
        return tuple(self._snpinfo)

    @property
    def dosage(self):
        return tuple(self._dosage)

    @property
    def gtnorm(self):
        return tuple(self._gtnorm)

    @property
    def gtcent(self):
        return tuple(self._gtcent)

    def _read_samples(self):
        if self._read_samples_once:
           return
        self._read_samples_once = True

        with open(self._samplefile, 'r') as samfile:
            # header = samfile.readline().strip().split()
            # header_types = samfile.readline().strip().split()
            sample = 0
            samplenames = list()
            for line in samfile:
                if re.search('^#', line):
                    continue
                sample += 1
                samplenames.append(line.strip().split()[0])
            
        self._nsample = sample
        self._samplenames = samplenames

    def _read_dosages(self):
        dosage = list()
        allsnps = list()
        linenum = 0
        print("started reading genotype")
        with gzip.open(self._gtfile, 'r') as filereader:
            for snpline in filereader:
                if linenum >= self._startsnp and linenum < self._endsnp:
                    self._nloci += 1
                    # print('reading '+str(self._nloci)+' snps, size of dosages is '+str(sys.getsizeof(dosage)), end='\r')
                    mline = snpline.split()
                    ngenotypes = len(mline) - self._data_columns
                    if ngenotypes != self._nsample:
                        raise ValueError('Number of samples differ from genotypes')

                    snp_dosage = np.array([float(x) for x in mline[self._data_columns:]])

                    maf = sum(snp_dosage) / 2 / len(snp_dosage)
                    try:                               ######## change to get the chrom numberfrom gtfile
                        chrom = int(mline[0])
                    except:
                        if self._chrom != None:
                            chrom = self._chrom
                        else:
                            chrom = -1
                    this_snp = SnpInfo(    chrom      = chrom,
                                           bp_pos     = int(mline[2]),
                                           varid      = mline[1].decode("utf-8"),
                                           ref_allele = mline[3].decode("utf-8"),
                                           alt_allele = mline[4].decode("utf-8"),
                                           maf        = maf)
                    allsnps.append(this_snp)
                    dosage.append(snp_dosage)
                linenum += 1
        return allsnps, np.array(dosage)

    def _read_genotype_freqs(self):
        dosage = list()
        allsnps = list()
        linenum = 0
        print("started reading genotype")
        with gzip.open(self._gtfile, 'r') as filereader:
            for snpline in filereader:
                if linenum >= self._startsnp and linenum < self._endsnp:
                    self._nloci += 1
                    # print('reading '+str(self._nloci)+' snps, size of dosages is '+str(sys.getsizeof(dosage)), end='\r')
                    mline = snpline.split()
                    ngenotypes = (len(mline) - self._data_columns) / 3
                    if float(ngenotypes).is_integer():
                        if ngenotypes != self._nsample:
                            print('Number of samples differ from genotypes')
                            raise SAMPLE_NUMBER_ERROR;
                    else:
                        print('Number of columns in genotype frequencies not divisible by 3')
                        raise GT_FREQS_ERROR;

                    gt_freqs = np.array([float(x) for x in mline[self._data_columns:]])

                    indsAA = np.arange(0,self._nsample)*3
                    indsAB = indsAA + 1
                    indsBB = indsAB + 1

                    snp_dosage = 2*gt_freqs[indsBB] + gt_freqs[indsAB] # [AA, AB, BB] := [0, 1, 2]

                    maf = sum(snp_dosage) / 2 / len(snp_dosage)
                    try:                          #######change this to read chrom number from gt file name
                        chrom = int(mline[0])
                    except:
                        if self._chrom != None:
                            chrom = self._chrom
                        else:
                            chrom = -1
                    this_snp = SnpInfo(    chrom      = chrom,
                                           bp_pos     = int(mline[2]),
                                           varid      = mline[1].decode("utf-8"),
                                           ref_allele = mline[3].decode("utf-8"),
                                           alt_allele = mline[4].decode("utf-8"),
                                           maf        = maf)

                    allsnps.append(this_snp)
                    dosage.append(snp_dosage)
                linenum += 1
        return allsnps, np.array(dosage)

    def _normalize_and_center_dosages(self):
        f = [snp.maf for snp in self._snpinfo]
        f = np.array(f).reshape(-1, 1)
        self._gtnorm = (self._dosage - (2 * f)) / np.sqrt(2 * f * (1 - f))
        self._gtcent = self._dosage - np.mean(self._dosage, axis = 1).reshape(-1, 1)

    def _read_genotypes(self):
        if self._read_genotype_once:
            return
        self._read_genotype_once = True
        self._read_samples() # otherwise, self._nsample is not set
        self._nloci = 0

        if self._isdosage:                    
            allsnps, dosage = self._read_dosages()
        else:
            allsnps, dosage = self._read_genotype_freqs()

        print("Read "+str(self._nloci)+" snps in "+str(self._nsample)+" samples.",end='\n')
        print("Finished readings snps")
        self._dosage = dosage
        self._snpinfo = allsnps

    def _filter_snps(self):
        # Predixcan style filtering of snps
        newsnps = list()
        for i, snp in enumerate(self._snpinfo):
            pos = snp.bp_pos
            refAllele = snp.ref_allele
            effectAllele = snp.alt_allele
            rsid = snp.varid
            # Skip non-single letter polymorphisms
            if len(refAllele) > 1 or len(effectAllele) > 1:
                continue
            # Skip ambiguous strands
            if SNP_COMPLEMENT[refAllele] == effectAllele:
                continue
            if rsid == '.':
                continue
            newsnps.append(snp)
        return newsnps

    def HWEcheck(self, genotype):
        bins = [0.66, 1.33]
        digitised_geno = np.digitize(genotype, bins).tolist()
        observed_freqs = np.array([0]*3)
        observed_freqs[0] = digitised_geno.count(0)
        observed_freqs[1] = digitised_geno.count(1)
        observed_freqs[2] = digitised_geno.count(2)
        n = sum(observed_freqs)
        p_A = (2*observed_freqs[0] + observed_freqs[1])/(2*n)
        p_a = (2*observed_freqs[2] + observed_freqs[1])/(2*n)
        X2 = n * ( (4*observed_freqs[0]*observed_freqs[2] - observed_freqs[1]**2)/
                   ((2*observed_freqs[0] + observed_freqs[1])*(2*observed_freqs[2] + observed_freqs[1]))
                 )**2
        pval = 1 - ss.chi2.cdf(X2, 1)
        return pval 

    def filter_snps(self, snpinfo, dosage):
        # Predixcan style filtering of snps
        newsnps = list()
        newdosage = list()
        nonSNP      = 0
        complements = 0
        unknown     = 0
        low_maf     = 0
        HWEfail     = 0
        for i, snp in enumerate(snpinfo):
            pos = snp.bp_pos
            refAllele = snp.ref_allele
            effectAllele = snp.alt_allele
            rsid = snp.varid
            maf = snp.maf

            # Skip non-single letter polymorphisms
            if len(refAllele) > 1 or len(effectAllele) > 1:
                nonSNP += 1
                continue
            # Skip ambiguous strands
            if SNP_COMPLEMENT[refAllele] == effectAllele:
                complements += 1
                continue
            # Skip unknown RSIDs
            if rsid == '.':
                unknown += 1
                continue
            # Skip low MAF
            if not (maf >= 0.10 and maf <=0.90):
                low_maf += 1
                continue
            if(self.HWEcheck(dosage[i]) < 0.000001):
                HWEfail += 1
                continue
            newsnps.append(snp)
            bins = [0.66, 1.33]
            # newdosage.append(np.digitize(dosage[i], bins))
            newdosage.append(dosage[i])
        print("{:d} / {:d} SNPs remain after filtering.".format(len(newsnps), len(snpinfo)))
        print(" Filtered {:d} non SNPs".format(nonSNP))
        print(" Filtered {:d} complement SNPs".format(complements))
        print(" Filtered {:d} unknown SNPs".format(unknown))
        print(" Filtered {:d} low maf SNPs".format(low_maf))
        print(" Filtered {:d} non HWE SNPs".format(HWEfail))
        return newsnps, np.array(newdosage)

    def write_dosages(self, outfile, format="gtex", filter_gt=True):
        if not self._read_genotype_once:
            raise ValueError("No dosages to write. Run read_genotypes first.",end='\n')
        print("Writing dosages to file "+outfile)

        dosage = self._dosage
        snpinfo = self._snpinfo
        
        if filter:
            snpinfo, dosage = self.filter_snps(snpinfo, dosage)

        if format == "predixcan":
            with open(outfile+".annot", 'w') as outstream2:
                with open(outfile, 'w') as outstream:
                    headers = "Id "+" ".join(self._samplenames)+"\n"
                    annotheader = "\t".join(["chr","pos","varID","refAllele","effectAllele","rsid"])
                    outstream.write(headers)
                    outstream2.write(annotheader+"\n")
                    for i, snp in enumerate(snpinfo):
                        variant_id = "_".join([str(snp.chrom), str(snp.bp_pos), snp.ref_allele, snp.alt_allele, "b37"])
                        annot_header = " ".join(["chr", "position", "VariantID", "RefAllele", "AlternativeAllele", "rsid", "rsid"])
                        dosage_row = " ".join([variant_id] + list(map(str, dosage[i])))
                        outstream.write(dosage_row+"\n")

                        annotline = "\t".join([str(snp.chrom), str(snp.bp_pos), variant_id, snp.ref_allele, snp.alt_allele, snp.varid])
                        outstream2.write(annotline+"\n")
            print("Done writing dosages")
        
        if format == "gtex":
            with gzip.open(outfile, 'wb') as outstream:
                for i, snp in enumerate(snpinfo):
                    dosage_line = " ".join([str(d) for d in dosage[i]])
                    snp_line = "{:d} {:s} {:d} {:s} {:s} {:f} ".format(snp.chrom, snp.varid, snp.bp_pos, snp.ref_allele, snp.alt_allele, snp.maf)
                    line = snp_line + dosage_line + "\n"
                    outstream.write(line.encode('utf-8'))
        # print("Writing finished.")

