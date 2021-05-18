'''
Script to generate the table to display Figure 1e and Extended Figure 1XXX
Input, raw Source data: 

-- discovery.tsv: raw output from boostdm with discovery index and other stats for each tumor-type gene pair
-- drivers.tsv: driver genes annotations from intOGen. 
-- ttype,f"{gene}.eval.pickle.gz": Evaluation stats (fscore50) of each gene-tumor type pair
-- linear_complexity.tsv: Linear complexity data for each gene-tumor type pair

Output, combined information : 

-- models_info.tsv 

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

# load discovery 

path_discovery_info = os.path.join(conf.output_boostdm,"discovery","discovery.tsv")
discovery = pd.read_csv(path_discovery_info, sep='\t')


# load drivers

drivers = pd.read_csv(conf.drivers_path, sep='\t')
discovery = discovery.merge(drivers.rename(columns={"SYMBOL":"gene","CANCER_TYPE":"ttype"})[["gene","ttype","ROLE"]],how="left")
discovery.drop_duplicates(inplace=True)

# add the f-score50 information 

def get_performance_pair(gene,ttype):
    with gzip.open(os.path.join(conf.output_boostdm,"evaluation",ttype,f"{gene}.eval.pickle.gz"),'rb') as f:
        values=pickle.load(f)["fscore50"]
        
    mean_fscore=np.nanmedian(values)
    q3,q1= np.percentile(values, [75 ,25])
    return mean_fscore,q3,q1
discovery[["mean_fscore50","up_fscore50","dn_fscore50"]] = discovery.apply(lambda row: pd.Series(get_performance_pair(row["gene"],row["ttype"])),axis=1)   
          
# add information about the model being selected 

discovery_tiers=[0, 0.2, 0.4, 0.6, 0.8]
mutation_tiers=[40, 30, 20, 10, 5]
fscore_threshold=.8
def is_cancer_specific(ttype):
    ttypes=oncotree.get_ttypes(ttype)
    if len(ttypes) == 1 and ttypes[0] == ttype:
        return True
    return False
    
def is_selected_model(row):

    fscore = row["mean_fscore50"]
    discovery = row["discovery_index"]
    n_muts = row["n_muts"]
    return is_cancer_specific(row["ttype"]) and (fscore >= fscore_threshold) and reduce(operator.or_, [(discovery >= dt) and (n_muts >= mt) for dt, mt in zip(discovery_tiers, mutation_tiers)])

discovery["selected"] = discovery.apply(lambda row: is_selected_model(row),axis=1 )

# add feature complexity

df_complexity = pd.read_csv("/workspace/projects/boostdm/feature_complexity/linear_complexity.tsv",sep="\t")
discovery = discovery.merge(df_complexity,how="left")

# save the combined table

discovery.to_csv(sys.argv[1],sep="\t",index=False)

