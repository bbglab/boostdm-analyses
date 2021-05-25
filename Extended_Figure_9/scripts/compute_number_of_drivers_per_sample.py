'''
Script to generate the table to display Extended Figure 8b) 
Input, raw Source data: 

-- discovery/mutations.tsv: total mutations observed in driver genes 
-- source_data/unique_observed_mutations_model_specificity_clinvar_oncokb.tsv.gz": for each mutation extract information about whether is a specific model (from Extended Figure 9b)
-- raw_data/dictionary_excess_dndscv.json: for each cohort from intogen, the excess according to dNdScv. 

Output, for each observed mutation, annotate the specificity of the model used for its prediction and whether is present in oncokb and clinvar : 

-- driver_per_sample_boostdm.tsv: table with number of drivers per sample according to boostdm


'''

import os
import gzip
import pickle
import numpy as np
import pandas as pd
from functools import reduce
import operator
import json
import sys
from collections import defaultdict
import glob
from tqdm import tqdm
import glob


# Import scripts and configuration files
import conf as conf
import oncotree as oncotree
oncotree = oncotree.Oncotree()
conf.config_params()


# load unique mutations and total observed mutations

df_muts_unique = pd.read_csv("source_data/unique_observed_mutations_model_specificity_clinvar_oncokb.tsv.gz",sep="\t")
df_muts_unique = df_muts_unique[(df_muts_unique["model_gene_ttype_specific"]==True) | (df_muts_unique["model_gene_other_ttype"]==True) | (df_muts_unique["selected_nottype_specific"]==True)][["gene","CANCER_TYPE","type_model"]].drop_duplicates() # only selected models
df_cohorts = pd.read_csv(conf.cohorts_path,sep="\t")
k={"NSCLC":"LUNG_CANCER"} # this tumor type has been renamed
df_cohorts["CANCER_TYPE"] = df_cohorts.apply(lambda row: row["CANCER_TYPE"] if not(row["CANCER_TYPE"] in k) else k[row["CANCER_TYPE"]],axis=1)

# load all mutations 
df_muts_total = pd.read_csv(conf.all_observed_mutations,sep="\t")
df_muts_total = df_muts_total.merge(df_cohorts[["COHORT","CANCER_TYPE"]]) # add cancer type

# define ttypes
leaves = oncotree.get_ttypes("CANCER")


# match them, only retain observed mutations with models

df_muts_total = df_muts_total.merge(df_muts_unique)

# for each model used, select driver mutations
def load_drivers_info(pairs,saturation):
    
    t=[]
    for gene,ttype,type_model in tqdm(pairs):
        
        if type_model=="gene/ttype specific":
        
            path=os.path.join(saturation,f"{gene}.{ttype}.prediction.tsv.gz")
            try:
                d=pd.read_csv(path,sep="\t") # out of memory, some predictions are pretty heavy
            except:
                continue
            d=d[d["boostDM_class"]==True][["gene","chr","pos","alt","boostDM_class"]].drop_duplicates().rename(columns={"alt":"mut"})
            d["CANCER_TYPE"] = ttype
            t.append(d)
        else:
            paths=glob.glob(os.path.join(saturation,f"{gene}.*.prediction.tsv.gz"))
            l=[]
            for path in paths:
                ttype_f = os.path.basename(path).split(".")[1]
                if ttype_f in leaves:
                    d=pd.read_csv(path,sep="\t") 
                    d=d[d["boostDM_class"]==True][["gene","chr","pos","alt","boostDM_class"]].drop_duplicates().rename(columns={"alt":"mut"})
                    l.append(d)
            if len(l) >0:
                c=pd.concat(l).groupby(["gene","chr","pos","mut"],as_index=False).agg({"boostDM_class":sum})
                c["CANCER_TYPE"]=ttype
                t.append(c)

    return pd.concat(t)
df_drivers=load_drivers_info(df_muts_unique[["gene","CANCER_TYPE","type_model"]].drop_duplicates().values.tolist(),os.path.join(conf.output_boostdm,"saturation","prediction"))


# we have the drivers, letch select observed mutations that are drivers
df_drivers.rename(columns={"alt":"mut"},inplace=True)
df_drivers["chr"] = df_drivers["chr"].astype(str)
df_muts_total["chr"] = df_muts_total["chr"].astype(str)
df_mut_drivers=df_muts_total.merge(df_drivers)

# now for each sample compute number of driver mutations according to boostdm 

d = json.load(open(conf.samples_intogen,'r'))
dict_info = {}
for key in d:
   
    dict_info[key]=d[key]
        
        
# count number of drivers per sample

counts_initial=df_mut_drivers.groupby(["CANCER_TYPE","sampleID"],as_index=False).agg({"boostDM_class":"count"})

# add samples with no driver mutation according to boostdm

l=[]
samples = counts_initial["sampleID"].unique()
for cohort in dict_info:
    try:
        ttype = df_cohorts[df_cohorts["COHORT"]==cohort]["CANCER_TYPE"].values[0]
    except:
        continue
    for sample in dict_info[cohort]:
        if not(sample in samples):
            l.append([ttype,sample,0.0])
counts_complementary=pd.DataFrame(l,columns=["CANCER_TYPE","sampleID","boostDM_class"])

# create the final dataframe

df_boostdm=pd.concat([counts_complementary,counts_initial]).drop_duplicates()


df_boostdm.to_csv("source_data/driver_per_sample_boostdm.tsv",sep="\t",index=False)