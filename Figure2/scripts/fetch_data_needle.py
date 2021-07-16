import sys
import os

import seaborn as sns
sns.set_style("whitegrid", {'axes.grid' : False})
import pandas as pd
import numpy as np

import conf
import oncotree


tree = oncotree.Oncotree()


def get_mutations():

    obs_mut = pd.read_csv(os.path.join(conf.output_boostdm, 'discovery', 'mutations.tsv'), sep='\t')
    obs_mut.rename(columns={'mut': 'alt'}, inplace=True)
    obs_mut['chr'] = obs_mut['chr'].astype(str)
    obs_mut['pos'] = obs_mut['pos'].astype(int)
    return obs_mut


def create_observed_dataset(gene, ttype, obs_mut):

    sat_pred = pd.read_csv(os.path.join(conf.output_boostdm,
                                        'saturation',
                                        'prediction',
                                        f'{gene}.{ttype}.prediction.tsv.gz'),
                           sep='\t')
    sat_pred['chr'] = sat_pred['chr'].astype(str)
    sat_pred['pos'] = sat_pred['pos'].astype(int)
    cohorts_under_ttype = list(zip(*tree.get_cohorts(ttype)))[0]
    obs_mut = obs_mut[(obs_mut['gene'] == gene) & (obs_mut['ttype']==ttype)]
    df = pd.merge(obs_mut, sat_pred, on=['chr', 'pos', 'alt', 'gene', 'aachange'], how='inner')
    return df.rename(columns={"ttype":"cancer_type"})


def get_position(row):

    try:
        v = int("".join(row["aachange"][1:-1]))
        return v
    except:
        return -1


def get_plot_data(data):

    df_stats = pd.read_csv(conf.cohorts_path, sep="\t")
    df_stats.rename(columns={"CANCER_TYPE": "cancer_type"},inplace=True)

    #data = data.merge(df_stats[["COHORT", "cancer_type"]])

    data["Protein_position"] = data.apply(lambda row: get_position(row), axis=1)
    data["AA"] = data.apply(lambda row: row["aachange"][0], axis=1)
    data['ID'] = data.apply(lambda x: '{}_{}'.format(x['pos'], x['alt']), axis=1)

    score_values = data['boostDM_score'].tolist()
    count_driver_unique = len(set(data[data['boostDM_class']]['ID']))
    count_driver = len(data[data['boostDM_class']]['ID'])
    count_total = len(data['ID'])
    count_total_unique = len(set(data['ID']))

    df_counts = data.groupby(["pos", "ENSEMBL_TRANSCRIPT", "cancer_type"],
                             as_index=False).agg({"sampleID": "count"})
    df_counts.rename(columns={"sampleID": "number_observed_muts"}, inplace=True)
    data = data.merge(df_counts, how='left')
    data = data.groupby(["pos", "AA", "Protein_position", "gene", "ENSEMBL_TRANSCRIPT",
                         "cancer_type", "number_observed_muts"], as_index=False).agg(
        {"boostDM_score": np.nanmax, "boostDM_class": np.any})

    data.rename(columns={"cancer_type": "cancer_type_annotations", }, inplace=True)

    return data, count_driver, count_driver_unique, count_total, count_total_unique, score_values
