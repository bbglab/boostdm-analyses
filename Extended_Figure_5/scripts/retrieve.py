# retrieve predictions and explanations

import os
import sys
sys.path.append('..')
import glob
import tqdm
import pickle
import gzip

import pandas as pd

import conf
import oncotree


tree = oncotree.Oncotree()

create_datasets_folder = os.path.join(conf.output_boostdm, 'create_datasets')
eval_folder = os.path.join(conf.output_boostdm, 'evaluation')
saturation_prediction_folder = os.path.join(conf.output_boostdm, 'saturation', 'prediction')
mutations_path = os.path.join(conf.output_boostdm, 'discovery', 'mutations.tsv')

ttypes = tree.get_ttypes('CANCER')

# mapping MoA colors to each (gene, ttype)
dg = pd.read_csv(conf.drivers_path, sep='\t')
d = dict(zip(zip(dg.SYMBOL, dg.CANCER_TYPE), map(lambda x: conf.dict_colors_role.get(x, '#808080'), dg.ROLE)))
d_gene = dict(zip(dg.SYMBOL, map(lambda x: conf.dict_colors_role.get(x, '#808080'), dg.ROLE)))

# mapping complexity to each (gene, ttype)
complex_df = pd.read_csv(conf.complexity_path, sep='\t')
complex_dict = dict(zip(zip(complex_df.gene, complex_df.ttype), complex_df.linear_complexity))

# model selection dict
with gzip.open(os.path.join(conf.output_boostdm, 'model_selection', 'eval_data.pickle.gz'), 'rb') as f:
    model_selection_dict = pickle.load(f)


def get_shaps(observed_mutations, gene=None, ttype=None):

    dg = observed_mutations
    if ttype is not None:
        cohorts = list(zip(*tree.get_cohorts(ttype)))[0]
        dg = dg[dg['ttype']==ttype].copy()
    if gene is not None:
        dg = dg[dg['gene'] == gene]
    dg.rename(columns={'mut': 'alt'}, inplace=True)

    data = []

    for fn in tqdm.tqdm(glob.glob(os.path.join(saturation_prediction_folder, '*.*.prediction.tsv.gz'))):

        g, tt = tuple((os.path.basename(fn)).split('.')[:2])
        cohorts = list(zip(*tree.get_cohorts(tt)))[0]
        dh = dg[dg['ttype']==tt]

        df_pred = pd.read_csv(fn, sep='\t')
        df_pred = df_pred[(df_pred['boostDM_score'] >= 0.5)]
        df_pred['moa'] = d.get((g, tt), d_gene.get((g, tt), '#808080'))
        df_pred['linear_complexity'] = df_pred.apply(lambda r: complex_dict.get((r['selected_model_gene'],
                                                                                 r['selected_model_ttype']), 0.),
                                                     axis=1)
        df_pred['ttype'] = tt
        df_pred['gene'] = g
        df_pred['aachange'].fillna(".", inplace=True)
        df_pred = df_pred[['chr', 'pos', 'alt'] +
                          conf.features +
                          ['moa', 'linear_complexity', 'aachange'] +
                          ['gene', 'ttype']]

        # merge
        dh['chr'] = dh['chr'].astype(str)
        df_pred['chr'] = df_pred['chr'].astype(str)
        df_data = pd.merge(dh, df_pred, on=['chr', 'pos', 'alt', 'gene', 'ttype', 'aachange'])

        data.append(df_data)

    df = pd.concat(data, axis=0)
    df.drop_duplicates(['chr', 'pos', 'alt', 'gene', 'ttype', 'aachange'], keep='first', inplace=True)
    return df


if __name__ == '__main__':
    pass


