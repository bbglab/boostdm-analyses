'''
Script to generate the table to display Extended Figure 8b) 
Input, raw Source data: 

-- discovery/mutations.tsv: unique mutations observed in driver genes
-- model_selection_information.tsv": stats about the models, including specifity and selection information (from Figure1)

Output, for each observed mutation, annotate the specificity of the model used for its prediction : 

-- unique_observed_mutations_model_specificity.tsv.gz

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
# Import scripts and configuration files
import conf as conf
import oncotree as oncotree
oncotree = oncotree.Oncotree()
conf.config_params()

# prepare dataset by loading drivers and their associated tumor type

df_drivers = pd.read_csv(conf.drivers_path,sep="\t")
df_cohorts = pd.read_csv(conf.cohorts_path,sep="\t")

d=df_cohorts[["COHORT","CANCER_TYPE"]].set_index("COHORT").to_dict()["CANCER_TYPE"]
d["HARTWIG_LUNG_NON_SMALL_CELL"] = "LUNG_CANCER" # this tumor type has been renamed
k={"NSCLC":"LUNG_CANCER"} # this tumor type has been renamed
df_drivers["CANCER_TYPE"] = df_drivers.apply(lambda row: row["CANCER_TYPE"] if not(row["CANCER_TYPE"] in k) else k[row["CANCER_TYPE"]],axis=1) # re-annotated lung cnacers

# read all unique observed mutations and assign tumor type
df=pd.read_csv(conf.all_observed_mutations,sep="\t",usecols=["chr","pos","ref","mut","aachange","ttype","gene"]).drop_duplicates().rename(columns={"ttype":"CANCER_TYPE"}) # only unique mutations
df["CANCER_TYPE"] = df.apply(lambda row: row["CANCER_TYPE"] if not(row["CANCER_TYPE"] in k) else k[row["CANCER_TYPE"]],axis=1) # re-annotated lung cnacers
s=df_drivers[["CANCER_TYPE","SYMBOL"]].drop_duplicates().rename(columns={"SYMBOL":"gene"})
   
df=df.merge(s) # only drivers in tumor type
df.drop_duplicates(inplace=True) # make them unique 
print ("Number of unique mutatoins...",df.shape[0])


# include information about type of models, whether they have been used and whether are invoked

df_info = pd.read_csv("../Figure1/source_data/model_selection_information.tsv",sep="\t")

# gather information about general models


with gzip.open(os.path.join(conf.output_boostdm, 'model_selection', 'eval_data.pickle.gz'),'rb') as f:
    d_info_models=pickle.load(f)

# retrieve tumor-types for which we want a prediction
effective = []
prediction_folder = os.path.join(conf.output_boostdm,"saturation","prediction")


for fn in glob.glob(os.path.join(prediction_folder, '*.*.prediction.tsv.gz')):
    gene, ttype = tuple(os.path.basename(fn).split('.')[:2])
    effective.append((ttype, gene))
# effectively used
used = []
cases = []
general_effective = []
leaves = oncotree.get_ttypes("CANCER")
for k, v in d_info_models.items():
    if (k[0] in leaves) and (v[0] is not None):
        if k[0] != v[0]:
            used.append(v)
    elif (k in effective) and (v[0] is not None):
        general_effective.append(k)
        used.append(v)

print(f'cases that invoke general models: {len(used)}\ngeneral models invoked: {len(set(used))}')

df_info["invoked"] = df_info.apply(lambda row: (row["ttype"],row["gene"]) in used ,axis=1 )
df_info["effective"] = df_info.apply(lambda row: (row["ttype"],row["gene"]) in effective ,axis=1 )
# save the data
df_info.to_csv("source_data/model_information_including_generals.tsv.gz",sep="\t",compression="gzip",index=False)

# gene-ttype specific models 

df_info = df_info[df_info["selected"]]
df_info["CANCER_TYPE"]= df_info["ttype"]
df_info["model_gene_ttype_specific"] = True
df = df.merge(df_info[["gene","CANCER_TYPE","model_gene_ttype_specific","ROLE"]].drop_duplicates(),how="left")
df["model_gene_ttype_specific"].fillna(False,inplace=True)


# model other tumor type

df_info["model_gene_other_ttype"] = True
df = df.merge(df_info[["gene","model_gene_other_ttype"]].drop_duplicates(),how="left")
df["model_gene_other_ttype"].fillna(False,inplace=True)

# gene models

df_info = pd.read_csv("source_data/model_information_including_generals.tsv.gz",sep="\t")
df_info["general_ttype_model"]=df_info.apply(lambda row: True  if (((row["ttype"],row["gene"]) in d_info_models) and (row["ttype"] != d_info_models[(row["ttype"],row["gene"])][0])) or (len(oncotree.get_ttypes(row["ttype"]))>1) else False, axis=1)

df_info=df_info[(df_info["effective"])&(df_info["general_ttype_model"])]
df_info["selected_nottype_specific"] = True
df_info["CANCER_TYPE"]= df_info["ttype"]
df = df.merge(df_info[["gene","selected_nottype_specific","CANCER_TYPE"]].drop_duplicates(),how="left")
df["selected_nottype_specific"].fillna(False,inplace=True)

# save mutations




# now prepare the ratios for the heatmap plot

print ("Number of mutations", df.shape[0])

def set_type_model(row):
    if row["model_gene_ttype_specific"]:
        return "gene/ttype specific"
    if row["selected_nottype_specific"]:
        return "gene specific"
    if row["model_gene_other_ttype"]:
        return "gene other ttype"
    return "no-model"
df["type_model"] = df.apply(lambda row: set_type_model(row),axis=1)

df.to_csv("source_data/unique_observed_mutations_with_type_model.tsv.gz",sep="\t",compression="gzip",index=False)
# iterate over tumor types, to get the number of mutations total 
l_results = []
for ttype in df["CANCER_TYPE"].unique():
    tmp = df[df["CANCER_TYPE"]==ttype]
    x=pd.Series(tmp["type_model"].value_counts())
    x["CANCER_TYPE"] = ttype
    x["num_muts"] = tmp.shape[0]
    l_results.append(x)
    
# transform to dataframe

df_comb = pd.DataFrame(l_results)
df_comb.fillna(0,inplace=True)
df_comb["ttype/gene specific_n"] = df_comb.apply(lambda row: row["gene/ttype specific"] / row["num_muts"],axis=1)
df_comb["gene specific_n"]= df_comb.apply(lambda row: row["gene specific"] / row["num_muts"],axis=1)
df_comb["gene other ttype_n"] = df_comb.apply(lambda row: row["gene other ttype"] / row["num_muts"],axis=1)
df_comb.reset_index(inplace=True)
df_comb.set_index('CANCER_TYPE',inplace=True)
df_comb.drop(columns=["index"],inplace=True)
df_comb.sort_values("num_muts",ascending=False,inplace=True)
df_table_ratios = df_comb[["ttype/gene specific_n","gene specific_n","gene other ttype_n"]]


    

df_table_ratios.to_csv("source_data/ratio_models_by_tumor_type_model.tsv",sep="\t",index=True)

