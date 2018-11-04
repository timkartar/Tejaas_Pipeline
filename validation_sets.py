import copy
from collections import defaultdict


class PlotSet():

    def __init__(self, parameters):
        self.name      = parameters["name"]
        self.basename  = parameters["basename"]
        self.title     = parameters["title"] 
        self.tissue    = parameters["tissue"]  
        self.tissue_title = parameters["tissue title"]
        self.cismask   = parameters["cismask"] 
        self.peer      = parameters["peer"]
        self.samples   = parameters["samples"]
        self.datadir   = parameters["datadir"]
        self.transfile = parameters["transfile"]
        self.cisfile   = parameters["cisfile"] 

sets = defaultdict(lambda: None)
sets["TEJAAS"] = defaultdict(lambda: None)
sets["MatrixEQTL"] = defaultdict(lambda: None)

tejaas_outdir = "output3/"

params = defaultdict(lambda: None)
params["name"]         = "cardio"
params["basename"]     = "cardiogenics"
params["title"]        = "Cardiogenics"
params["tissue"]       = "mono"
params["tissue title"] = "Monocytes"
params["cismask"]      = "cismask"
params["peer"]         = ""
params["samples"]      = ""
params["datadir"]      = tejaas_outdir+"cardio_mono_cismask"
sets["TEJAAS"]["cardio_mono_cismask"] = copy.deepcopy(params)

params["cismask"]      = "cismask"
params["peer"]         = ""
params["samples"]      = "us"
params["datadir"]      = tejaas_outdir+"cardio_mono_cismask_us"
sets["TEJAAS"]["cardio_mono_cismask_us"] = copy.deepcopy(params)

params["cismask"]      = "nocismask"
params["peer"]         = ""
params["samples"]      = ""
params["datadir"]      = tejaas_outdir+"cardio_mono"
sets["TEJAAS"]["cardio_mono"]    = copy.deepcopy(params)


params = defaultdict(lambda: None)
params["name"]         = "cardio"
params["basename"]     = "cardiogenics"
params["title"]        = "Cardiogenics"
params["tissue"]       = "macro"
params["tissue title"] = "Macrophages"
params["cismask"]      = "cismask"
params["peer"]         = ""
params["samples"]      = ""
params["datadir"]      = tejaas_outdir+"cardio_macro_cismask"
sets["TEJAAS"]["cardio_macro_cismask"] = copy.deepcopy(params)

params["cismask"]      = "cismask"
params["peer"]         = ""
params["samples"]      = "us"
params["datadir"]      = tejaas_outdir+"cardio_macro_cismask_us"
sets["TEJAAS"]["cardio_macro_cismask_us"] = copy.deepcopy(params)

params["cismask"]      = "nocismask"
params["peer"]         = ""
params["samples"]      = ""
params["datadir"]      = tejaas_outdir+"cardio_macro"
sets["TEJAAS"]["cardio_macro"] = copy.deepcopy(params)




params = defaultdict(lambda: None)
params["name"]         = "gtex"
params["basename"]     = "gtex"
params["title"]        = "GTEx"
params["tissue"]       = "wb"
params["tissue title"] = "Whole Blood"
params["cismask"]      = "cismask"
params["peer"]         = ""
params["samples"]      = ""
params["datadir"]      = tejaas_outdir+"gtex_wb_cismask"
sets["TEJAAS"]["gtex_wb_cismask"] = copy.deepcopy(params)


params["cismask"]      = "nocismask"
params["peer"]         = ""
params["samples"]      = ""
params["datadir"]      = tejaas_outdir+"gtex_wb"
sets["TEJAAS"]["gtex_wb"] = copy.deepcopy(params)

params["tissue"]       = "hlv"
params["tissue title"] = "Heart Left Ventricle"
params["cismask"]      = "cismask"
params["peer"]         = ""
params["samples"]      = ""
params["datadir"]      = tejaas_outdir+"gtex_hlv_cismask"
sets["TEJAAS"]["gtex_hlv_cismask"] = copy.deepcopy(params)

params["cismask"]      = "nocismask"
params["peer"]         = ""
params["samples"]      = ""
params["datadir"]      = tejaas_outdir+"gtex_hlv"
sets["TEJAAS"]["gtex_hlv"] = copy.deepcopy(params)


params["tissue"]       = "ms"
params["tissue title"] = "Muscle Skeletal"
params["cismask"]      = "cismask"
params["peer"]         = ""
params["samples"]      = ""
params["datadir"]      = tejaas_outdir+"gtex_ms_cismask"
sets["TEJAAS"]["gtex_ms_cismask"] = copy.deepcopy(params)

params["cismask"]      = "nocismask"
params["peer"]         = ""
params["samples"]      = ""
params["datadir"]      = tejaas_outdir+"gtex_ms"
sets["TEJAAS"]["gtex_ms"] = copy.deepcopy(params)


matrixeqtl_outdir = "matrixEQTL_out_filteredGT_gencode/"

params = defaultdict(lambda: None)
params["name"]         = "gtex"
params["basename"]     = "gtex"
params["title"]        = "GTEx"
params["tissue"]       = "wb"
params["tissue title"] = "Whole Blood"
params["cismask"]      = "cismask"
params["peer"]         = ""
params["samples"]      = ""
params["datadir"]      = matrixeqtl_outdir+"gtex_wb"
params["transfile"]    = "gtex_MatrixEQTL_chr{:d}.transout"
params["cisfile"]      = "gtex_MatrixEQTL_chr{:d}.cisout"
sets["MatrixEQTL"]["gtex_wb"] = copy.deepcopy(params)

params["tissue"]       = "hlv"
params["tissue title"] = "Heart Left Ventricle"
params["datadir"]      = matrixeqtl_outdir+"gtex_hlv"
sets["MatrixEQTL"]["gtex_hlv"] = copy.deepcopy(params)

params["tissue"]       = "ms"
params["tissue title"] = "Muscle Skeletal"
params["datadir"]      = matrixeqtl_outdir+"gtex_ms"
sets["MatrixEQTL"]["gtex_ms"] = copy.deepcopy(params)


params = defaultdict(lambda: None)
params["name"]         = "cardio"
params["basename"]     = "cardiogenics"
params["title"]        = "Cardiogenics"
params["tissue"]       = "mono"
params["tissue title"] = "Monocytes"
params["cismask"]      = "cismask"
params["peer"]         = ""
params["samples"]      = ""
params["datadir"]      = matrixeqtl_outdir+"cardiogenics_mono"
params["transfile"]    = "cardiogenics_MatrixEQTL_chr{:d}.transout"
params["cisfile"]      = "cardiogenics_MatrixEQTL_chr{:d}.cisout"
sets["MatrixEQTL"]["cardio_mono"] = copy.deepcopy(params)

params["peer"]         = ""
params["samples"]      = "us"
params["datadir"]      = matrixeqtl_outdir+"cardiogenics_mono_us"
sets["MatrixEQTL"]["cardio_mono_us"] = copy.deepcopy(params)


params = defaultdict(lambda: None)
params["name"]         = "cardio"
params["basename"]     = "cardiogenics"
params["title"]        = "Cardiogenics"
params["tissue"]       = "macro"
params["tissue title"] = "Macrophages"
params["cismask"]      = "cismask"
params["peer"]         = ""
params["samples"]      = ""
params["datadir"]      = matrixeqtl_outdir+"cardiogenics_macro"
params["transfile"]    = "cardiogenics_MatrixEQTL_chr{:d}.transout"
params["cisfile"]      = "cardiogenics_MatrixEQTL_chr{:d}.cisout"
sets["MatrixEQTL"]["cardio_macro"] = copy.deepcopy(params)

params["peer"]         = ""
params["samples"]      = "us"
params["datadir"]      = matrixeqtl_outdir+"cardiogenics_macro_us"
sets["MatrixEQTL"]["cardio_macro_us"] = copy.deepcopy(params)
