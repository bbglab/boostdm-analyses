'''
Script to generate the tables to display the explanation radar plots of all observed mutations in EGFR across GBM and LUAD
Input, raw Source data: 

-- {gene}.{ttype}.prediction.tsv.gz: Saturation prediction of all mutations in gene == EGFR across ttype == GBM and LUAD
-- cohorts.tsv: IntOGen derived information about the datasets and their associated tumor type  
-- create_datasets/{cohort}.regression_data.tsv: Observed and randomized mutations for each cohort. 

Output, combined information : 

-- observed_mutations_with_shapleys_{gene}_{ttype}.tsv 

'''


import os
import gzip
import pickle
import numpy as np
import pandas as pd
from functools import reduce
import operator
import sys
# Import scripts and configuration files
import conf as conf
import oncotree as oncotree
oncotree = oncotree.Oncotree()
conf.config_params()

# EGFR in LUAD and GBM 

gene="EGFR"
path_saturation=os.path.join(conf.output_boostdm,"saturation","prediction")
df_stats = pd.read_csv(conf.cohorts_path,sep="\t")

for ttype in ["LUAD","GBM"]:
    path_file = os.path.join(path_saturation,f"{gene}.{ttype}.prediction.tsv.gz")
    df_data=pd.read_csv(path_file,sep="\t")
    cohorts_ttype = df_stats[df_stats["CANCER_TYPE"]==ttype]["COHORT"].values
    muts = []
    l_data =[]
    for cohort in cohorts_ttype:
        input_=os.path.join(conf.output_boostdm,"create_datasets",f"{cohort}.regression_data.tsv")
        k = pd.read_csv(input_,sep="\t")
        l_data.append(k[(k["response"]==1)&(k["gene"]==gene)][["gene","chr","pos","ref","alt"]]) # select observed mutations (response==1)
    df_observed = pd.concat(l_data)
    df_observed["pos"] = df_observed["pos"].astype(int)
    df_observed["chr"] = df_observed["chr"].astype(int)
    df_data=df_observed.merge(df_data) # include boostdm Predictions
    df_data.to_csv(f"source_data/observed_mutations_with_shapleys_{gene}_{ttype}.tsv",sep="\t")
