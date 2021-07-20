'''
Script to generate the table to display Figure 1e and Extended Figure 1XXX
Input, raw Source data: 

-- discovery/mutations.tsv: unique mutations observed in driver genes
-- model_selection_information.tsv": stats about the models, including specifity and selection information

Output, combined information : 

-- observed_drivers_ratio.tsv

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
#df.drop(columns=["COHORT"],inplace=True)
s=df_drivers[["CANCER_TYPE","SYMBOL"]].drop_duplicates().rename(columns={"SYMBOL":"gene"})
   
df=df.merge(s) # only drivers in tumor type
df.drop_duplicates(inplace=True) # make them unique 
print ("Number of unique mutations...",df.shape[0])


# include information about type of model for each mutation...

df_info = pd.read_csv(os.path.join(os.environ["PATH_SOURCE_DATA"], "extended-figure1-source-data", "model_selection_information.tsv"), sep="\t")
df_info = df_info[df_info["selected"]]
df_info["CANCER_TYPE"]= df_info["ttype"]
df_info["model_gene_ttype_specific"] = True
df = df.merge(df_info[["gene","CANCER_TYPE","model_gene_ttype_specific","ROLE"]].drop_duplicates(),how="left")
df["model_gene_ttype_specific"].fillna(False,inplace=True)


# select only gene-ttype specific models 

v=df[df["model_gene_ttype_specific"]]

# load predictions

saturation=os.path.join(conf.output_boostdm,"saturation","prediction")
l=[]
for gene,ttype in v[["gene","CANCER_TYPE"]].drop_duplicates().values.tolist():
    try:
        df = pd.read_csv(os.path.join(saturation,f"{gene}.{ttype}.prediction.tsv.gz"),sep="\t").rename(columns={"alt":"mut"})
        df["chr"] =df["chr"].astype(str)
    except:
        #print (gene,ttype)
        continue
    df["CANCER_TYPE"] = ttype
    df = df[["boostDM_class","gene","chr","pos","mut","aachange","CANCER_TYPE"]].drop_duplicates().merge(v)
    l.append(df)

df_m = pd.concat(l).drop_duplicates()    

# compute ratios

ratios=df_m.groupby(["gene","CANCER_TYPE","ROLE"],as_index=False).agg({"boostDM_class":sum,"pos":"count"}).rename(columns={"boostDM_class":"predicted_drivers","pos":"unique_mutations"})
ratios["ratio_drivers"] =  ratios["predicted_drivers"] / ratios["unique_mutations"]


# save the info

ratios.to_csv(sys.argv[1],sep="\t",index=False)