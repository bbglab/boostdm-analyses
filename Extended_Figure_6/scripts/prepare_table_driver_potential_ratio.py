'''
Script to generate the tables to display the explanation radar plots of all observed mutations in EGFR across GBM and LUAD
Input, raw Source data: 

-- Figure1/source_data/model_selection_information.tsv: information of selected models

Output, ratio of potentials drivers : 

-- source_data/potentials_driver_ratio_info.tsv

'''


import os
import gzip
import pickle
import numpy as np
import pandas as pd
import operator
import tqdm
import sys
import glob

# Import scripts and configuration files
import conf as conf
import oncotree as oncotree

oncotree = oncotree.Oncotree()
conf.config_params()

# get information of all models

df_info = pd.read_csv("../Figure1/source_data/model_selection_information.tsv",sep="\t")
df_info = df_info[df_info["selected"]] # only selected
df_info["CANCER_TYPE"]= df_info["ttype"]

# iterate over boostdm output and compute the potentials driver ratio
saturation=os.path.join(conf.output_boostdm,"saturation","prediction")
l=[]
for gene,ttype,role in df_info[["gene","ttype","ROLE"]].values.tolist():
    try:
        df = pd.read_csv(os.path.join(saturation,f"{gene}.{ttype}.prediction.tsv.gz"),sep="\t")
    except:
        print (gene,ttype)
        continue
    total = df.shape[0]
    drivers = df[df["boostDM_class"]].shape[0]
    l.append([gene,ttype,drivers/total])
    
df_c = pd.DataFrame(l,columns=["gene","ttype","potential_driver_ratio"])
output=df_c.merge(df_info)
output.to_csv("source_data/potentials_driver_ratio_info.tsv",sep="\t")