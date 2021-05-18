'''
Script to generate the tables to display the explanation radar plots of all observed mutations in EGFR across GBM and LUAD
Input, raw Source data: 

-- {gene}.{ttype}.prediction.tsv.gz: Saturation prediction of all mutations in gene == EGFR across all tumor types
-- create_datasets/{cohort}.regression_data.tsv: Observed and randomized mutations for each cohort. 

Output, combined information : 

-- observed_mutations_with_shapleys_{gene}_{ttype}.tsv 

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

def create_saturation_vectors(gene):
    saturation_folder=os.path.join(conf.output_boostdm, "saturation", "prediction")
    res = {}
    for fn in tqdm.tqdm(glob.glob(os.path.join(saturation_folder, f'{gene}.*.prediction.tsv.gz'))):
        gene = os.path.basename(fn).split('.')[0]
        ttype = os.path.basename(fn).split('.')[1]
        df = pd.read_csv(fn, sep='\t')
        df.drop_duplicates(inplace=True)
        
        df['Protein_position'] = df.apply(lambda row: get_position(row), axis=1)
        df = df[~df['aachange'].isnull()]
        df['AA'] = df.apply(lambda row: row["aachange"][0], axis=1)
        # at nucleotide level
        df3 = df.groupby(['pos', 'AA', 'Protein_position', 'alt', 'ENSEMBL_TRANSCRIPT'], as_index=False).agg({'boostDM_score': np.nanmax, 'boostDM_class': np.any})
        df3['gene'] = gene
        df3['cancer_type_annotations'] = ttype
        df3.sort_values(by='Protein_position', ascending=True, inplace=True)
        
        saturation_vector = np.array(list(map(int, df3['boostDM_class'].values)))
        aa_vector = np.array(df3['AA'].values)
        model = tuple(df.loc[0, ['selected_model_ttype', 'selected_model_gene']])
        if saturation_vector is not None:
            res[gene] = res.get(gene, {})
            res[gene].update({ttype: (model, saturation_vector, aa_vector)})
    return res

def get_PFAMs_per_transcript(gene,df_pfam):
    
    df_pfam_gene = df_pfam[(df_pfam["gene"] == gene)]
    df_pfam_gene = df_pfam_gene[["gene", "START", "END", "DOMAIN"]].drop_duplicates()
    df_pfam_gene = pd.merge(df_pfam_gene, df_names[["DOMAIN", "DOMAIN_NAME"]].drop_duplicates(), how="left")
    df_pfam_gene["POS"] = df_pfam_gene.apply(lambda row: row["START"] + ((row["END"] - row["START"]) // 2), axis=1)
    df_pfam_gene["SIZE"] = df_pfam_gene.apply(lambda row: row["END"] - row["START"] + 1, axis=1)
    df_pfam_gene["Color"] = "#998ec3"
    return df_pfam_gene

def get_position(row):

    try:
        v = int("".join(row["aachange"][1:-1]))
        return v
    except:
        return -1

# create saturation vectors 
gene="EGFR"
ttype_blueprints = create_saturation_vectors(gene)


# prepare a paired dictionary 
all_pred_specific = {}
for gene in tqdm.tqdm(ttype_blueprints):
    for ttype in ttype_blueprints[gene]:
        all_pred_specific[(ttype, gene)] = (True, ttype_blueprints[gene][ttype][1])
        
# prepare PFAM domains

df_pfam = pd.read_csv(conf.PFAM_files, sep="\t", names=["ENSEMBL_GENE", "ENSEMBL_TRANSCRIPT", "START", "END", "DOMAIN"])
df_names = pd.read_csv(conf.PFAM_info, sep="\t", names=["DOMAIN", "CLAN", "CLAN_NAME", "DOMAIN_NAME", "Long Name"])
transcript_gene = pd.read_csv(conf.path_coord,sep="\t",
                              usecols=[1,10], names=["gene", "ENSEMBL_TRANSCRIPT"]).drop_duplicates()
df_pfam = df_pfam.merge(transcript_gene)
df_pfam=get_PFAMs_per_transcript(gene,df_pfam)


# prepare degrons

df_degrons = pd.read_csv(conf.degrons_path,sep="\t")
df_uniprot = pd.read_csv(conf.uniprot_conv,sep="\t").rename(columns={"UniProtKB isoform ID":"Entry_Isoform"})
df_degrons = df_degrons.merge(df_uniprot)
df_degrons["SIZE"] = df_degrons.apply(lambda row: row["END"] - row["START"] + 1,axis=1)
d_transcripts=transcript_gene.set_index("gene").to_dict()["ENSEMBL_TRANSCRIPT"]
ts = d_transcripts[gene]
df_degrons=df_degrons[df_degrons["Transcript stable ID"]==ts]

# save the data

df_pfam.to_csv(f"source_data/PFAM_domains_{gene}.tsv",sep="\t",index=False)
df_degrons.to_csv(f"source_data/annotated_degrons_{gene}.tsv",sep="\t",index=False)
with open(f"source_data/blueprint_info_{gene}.pickle","wb") as fn:
    pickle.dump(all_pred_specific,fn)
