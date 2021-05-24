import os
import gzip
import sys
sys.path.append('../')
import pickle
import json
import tqdm
import glob
import click
from multiprocessing import Pool

import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score

from mutrate import set_mutrate

import conf
import oncotree

#from .. import conf
#from .. import oncotree

tree = oncotree.Oncotree()

# load model selection

model_selection_pickle = os.path.join(conf.output_boostdm, 'model_selection', 'eval_data.pickle.gz')
with gzip.open(model_selection_pickle, 'rb') as f:
    model_selection_dict = pickle.load(f)

# cohorts and driver genes

df_drivers = pd.read_csv(conf.drivers_path, sep='\t')
driver_gene_ttypes = set(map(tuple, df_drivers[['SYMBOL', 'CANCER_TYPE']].drop_duplicates().values.tolist()))
df_drivers_small = df_drivers[['COHORT', 'CANCER_TYPE']].drop_duplicates()
cohort_ttype_dict = dict(zip(df_drivers_small['COHORT'], df_drivers_small['CANCER_TYPE']))
valid_cohorts = set(df_drivers['COHORT'].values)
driver_genes = df_drivers['SYMBOL'].drop_duplicates().values.tolist()

ttype_dict = {}
df_stats_cohort = pd.read_csv(conf.cohorts_path, sep='\t')
for tt, cohort in df_stats_cohort[['CANCER_TYPE', 'COHORT']].drop_duplicates().values:
    ttype_dict[tt] = ttype_dict.get(tt, []) + [cohort]


# load observed mutations

def load_observed():
    print('Loading observed mutations...')
    fn = os.path.join(conf.output_boostdm, 'discovery', 'mutations.tsv')
    df = pd.read_csv(fn, sep='\t')
    df['ttype'] = df['COHORT'].apply(lambda x: cohort_ttype_dict.get(x, None))
    return df


observed_mutations = load_observed()


# load saturation prediction

def load_saturation():
    print('Loading in silico saturation mutagenesis...')
    fn = '/workspace/projects/boostdm/benchmark_supplement/datasets/saturation_aggregate.tsv.gz'
    df = pd.read_csv(fn, sep='\t', low_memory=False)
    return df


saturation_prediction = load_saturation()


def get_nondriver_gene_ttype():

    """
    Get the collection of all non-driver gene-ttype pairs,
    with gene being a driver in at least another ttype
    """

    # mutrate is defined for a given gene-ttype if defined in some gene-cohort with cohort in ttype

    ttypes = tree.get_ttypes('CANCER')
    mutrate_gene_ttypes = set()
    for ttype in tqdm.tqdm(ttypes):
        for cohort in ttype_dict.get(ttype, []):
            try:
                for fn in glob.glob(os.path.join(conf.mutrate_folder, f'{cohort}.mutrate_output', 'norm_*.out.json')):
                    gene = os.path.basename(fn).split('_')[1].split('.')[0]
                    mutrate_gene_ttypes.add((gene, ttype))
            except FileNotFoundError as e:
                print(e)
                continue

    # get all mutrates from all possible genes with observed mutations
    # only genes that are classified as driver for some ttype

    gene_cohort = set(list(zip(observed_mutations['gene'], observed_mutations['COHORT'])))
    all_gene_ttypes = set(map(lambda x: (x[0], cohort_ttype_dict.get(x[1], None)), gene_cohort))

    # driver_genes = df_drivers['SYMBOL'].unique()
    # all_gene_ttypes = {x for x in all_gene_ttypes if (x[0] in driver_genes) and (x[1] is not None)}

    all_gene_ttypes = {x for x in all_gene_ttypes if (x[1] is not None)}

    # identify all non-driver pairs

    pairs = set(map(tuple, df_drivers[['SYMBOL', 'CANCER_TYPE']].drop_duplicates().values.tolist()))
    all_non_drivers = list((mutrate_gene_ttypes.intersection(all_gene_ttypes)).difference(pairs))
    return all_non_drivers


def put_chr(chromosome):

    c = str(chromosome)
    if not c.startswith('chr'):
        return 'chr' + c
    else:
        return c


def get_observed_mutations(gene, ttype, driver=True):

    df = observed_mutations[(observed_mutations['ttype'] == ttype) & (observed_mutations['gene'] == gene)].copy()
    df.drop_duplicates(['chr', 'pos', 'mut'], inplace=True)
    dg = df.merge(saturation_prediction[['gene', 'chr', 'pos', 'alt', 'boostDM_score']],
                  right_on=['gene', 'chr', 'pos', 'alt'],
                  left_on=['gene', 'chr', 'pos', 'mut'])
    if driver:
        dg = dg[dg['boostDM_score'] >= 0.5]
    sat_hash = dg.apply(lambda r: put_chr(r['chr']) + '.' + str(r['pos']) + '.' + str(r['alt']), axis=1)
    return set(sat_hash.values)


def create_observed_vs_non_observed(driver=False):

    res = {}

    for sat_fn in tqdm.tqdm(glob.glob(os.path.join(conf.output_boostdm,
                                                   'saturation', 'prediction',
                                                   '*.*.prediction.tsv.gz'))):

        gene, ttype = tuple(os.path.basename(sat_fn).split('.')[:2])

        print(gene, ttype)

        # get all possible mutations

        df_sat = pd.read_csv(sat_fn, sep='\t')
        if driver:
            df_sat = df_sat[df_sat['boostDM_score'] >= 0.5]
        sat_hash = df_sat.apply(lambda r: put_chr(r['chr']) + '.' + str(r['pos']) + '.' + str(r['alt']), axis=1)
        all_hashes = set(sat_hash.values)

        # get observed mutations

        observed = get_observed_mutations(gene, ttype, driver=driver)

        # get non observed mutations

        non_observed = all_hashes.difference(observed)

        # mark whether there is a specific model for gene, ttype

        ttype_model, gene_model = model_selection_dict.get((ttype, gene), (None, None))
        specific = False
        if (ttype_model == ttype) and (gene_model == gene):
            specific = True

        res[(ttype, gene)] = (specific, {'positive': observed, 'negative': non_observed})

    return res


def create_driver_vs_passenger(non_observed=False):

    res = {}

    leaves = tree.get_ttypes('CANCER')
    with gzip.open(os.path.join(conf.output_boostdm, 'model_selection', 'eval_data.pickle.gz'), 'rb') as g:
        model_dict = pickle.load(g)
    gene_ttypes = set(k for k, v in model_dict.items() if (k == v) and (k[0] in leaves))

    for ttype, gene in tqdm.tqdm(gene_ttypes):

        sat_fn = os.path.join(conf.output_boostdm, 'saturation', 'prediction', f'{gene}.{ttype}.prediction.tsv.gz')

        observed = set()
        if non_observed:
            observed = get_observed_mutations(gene, ttype, driver=False)

        # get all possible mutations
        try:
            df_sat = pd.read_csv(sat_fn, sep='\t')
        except FileNotFoundError as e:
            print(e)
            continue

        df_driver = df_sat[df_sat['boostDM_score'] >= 0.5]
        sat_driver = df_driver.apply(lambda r: put_chr(r['chr']) + '.' + str(r['pos']) + '.' + str(r['alt']), axis=1)
        driver = set(sat_driver.values)

        df_passenger = df_sat[df_sat['boostDM_score'] < 0.5]
        sat_passenger = df_passenger.apply(lambda r: put_chr(r['chr']) + '.' + str(r['pos']) + '.' + str(r['alt']), axis=1)
        passenger = set(sat_passenger.values)

        ttype_model, gene_model = model_selection_dict.get((ttype, gene), (None, None))
        specific = False
        if (ttype_model == ttype) and (gene_model == gene):
            specific = True

        res[(ttype, gene)] = (specific, {'positive': driver.difference(observed),
                                         'negative': passenger.difference(observed)})

    return res


saturation_annotation = os.path.join(conf.output_boostdm, 'saturation', 'annotation')


def create_non_drivers_func(item):

    gene, ttype = item
    try:
        sat_fn = next(iter(glob.glob(os.path.join(saturation_annotation, '*', f'{gene}.annotated.out.gz'))))
    except StopIteration:
        return dict()

    df_sat = pd.read_csv(sat_fn, sep='\t', low_memory=False)
    sat_hash = df_sat.apply(lambda r: put_chr(r['chr']) + '.' + str(r['pos']) + '.' + str(r['alt']), axis=1)
    all_hashes = set(sat_hash.values)

    # get observed mutations
    observed = get_observed_mutations(gene, ttype, driver=False)
    non_observed = all_hashes.difference(observed)

    return {(ttype, gene): (False, {'positive': observed, 'negative': non_observed})}


def create_non_drivers(cores=10):

    with open('./datasets/all_non_drivers.json', 'rt') as f:
        all_non_drivers = json.load(f)

    total = dict()
    with Pool(cores) as pool:
        for res in tqdm.tqdm(pool.imap(create_non_drivers_func, all_non_drivers), total=len(all_non_drivers)):
            total.update(res)

    return total


"""
def create_non_drivers(cores=10):

    res = {}
    with open('./datasets/all_non_drivers.json', 'rt') as f:
        all_non_drivers = json.load(f)
    saturation_annotation = os.path.join(conf.output_boostdm, 'saturation', 'annotation')
    for gene, ttype in tqdm.tqdm(all_non_drivers):

        # we just need to retrieve all genomic coordinates for each gene
        try:
            sat_fn = next(iter(glob.glob(os.path.join(saturation_annotation, '*', f'{gene}.annotated.out.gz'))))
        except:
            continue

        df_sat = pd.read_csv(sat_fn, sep='\t')
        sat_hash = df_sat.apply(lambda r: put_chr(r['chr']) + '.' + str(r['pos']) + '.' + str(r['alt']), axis=1)
        all_hashes = set(sat_hash.values)

        # get observed mutations
        observed = get_observed_mutations(gene, ttype, driver=False)
        non_observed = all_hashes.difference(observed)

        res[(ttype, gene)] = (False, {'positive': observed,
                                      'negative': non_observed})

    return res
"""


def create_csqn_type(csqn_type='nonsense'):

    res = {}

    for sat_fn in tqdm.tqdm(glob.glob(os.path.join(conf.path_saturation, f'*.*.prediction.out.gz'))):

        gene, ttype = tuple(os.path.basename(sat_fn).split('.')[:2])

        # get all possible mutations

        df_sat = pd.read_csv(sat_fn, sep='\t')

        df_csqn = df_sat[df_sat[f'csqn_type_{csqn_type}'] == 1]
        sat_csqn = df_csqn.apply(lambda r: put_chr(r['chr']) + '.' + str(r['pos']) + '.' + str(r['alt']), axis=1)
        csqn = set(sat_csqn.values)

        sat_all = df_sat.apply(lambda r: put_chr(r['chr']) + '.' + str(r['pos']) + '.' + str(r['alt']), axis=1)
        all = set(sat_all.values)

        ttype_model, gene_model = model_selection_dict.get((ttype, gene), (None, None))
        specific = False
        if (ttype_model == ttype) and (gene_model == gene):
            specific = True

        other = all.difference(csqn)
        res[(ttype, gene)] = (specific, {'positive': csqn, 'negative': other})

    return res


def create_feature(feature, threshold=0.0):

    """
    Separates between mutations by feature values
    If quantile=None, it interprets the feature as a boolean flag
    If quantile is not None, it uses it as quantile description
    """

    res = {}

    for sat_fn in tqdm.tqdm(glob.glob(os.path.join(conf.path_saturation, f'*.*.prediction.out.gz'))):

        gene, ttype = tuple(os.path.basename(sat_fn).split('.')[:2])

        df = pd.read_csv(sat_fn, sep='\t')

        # positive set
        dg = df[df[feature] >= threshold]
        sat = dg.apply(lambda r: put_chr(r['chr']) + '.' + str(r['pos']) + '.' + str(r['alt']), axis=1)
        positive = set(sat.values)

        # all mutations
        sat_all = df.apply(lambda r: put_chr(r['chr']) + '.' + str(r['pos']) + '.' + str(r['alt']), axis=1)
        all = set(sat_all.values)

        ttype_model, gene_model = model_selection_dict.get((ttype, gene), (None, None))
        specific = False
        if (ttype_model == ttype) and (gene_model == gene):
            specific = True

        negative = all.difference(positive)
        res[(ttype, gene)] = (specific, {'positive': positive, 'negative': negative})

    return res


def create_saturation_vectors(specific_saturation_folder):

    res = {}
    for fn in tqdm.tqdm(glob.glob(os.path.join(specific_saturation_folder, '*.*.prediction.tsv.gz'))):
        gene = os.path.basename(fn).split('.')[0]
        ttype = os.path.basename(fn).split('.')[1]
        if (gene, ttype) in driver_gene_ttypes:
            df = pd.read_csv(fn, sep='\t')
            df.drop_duplicates(inplace=True)
            if df.shape[0] == 0 or df[pd.isnull(df["aachange"])].shape[0] > 0:
                continue
            df["protein_position"] = df.apply(lambda row: int(row["aachange"][1:-1]), axis=1)
            df.sort_values("protein_position", inplace=True)
            saturation_vector = np.array(list(map(int, df['boostDM_class'].values)))
            aa_mutation_vector = np.array(df['aachange'].values)
            model = tuple(df.iloc[0, [-5, -4]])
            if saturation_vector is not None:
                res[gene] = res.get(gene, {})
                res[gene].update({ttype: (model, saturation_vector, aa_mutation_vector)})
    return res


def create_mutrate_dict_func(item):

    gene, ttype = item
    probs = np.zeros(96)
    burden = []
    for cohort in ttype_dict[ttype]:
        path_cohort = os.path.join(conf.mutrate_folder, cohort + ".mutrate_output", f"norm_{gene}.out.json")
        try:
            with open(path_cohort, 'rt') as g:
                dict_probs = json.load(g)
        except FileNotFoundError:
            return dict()
        for sample in dict_probs[gene].keys():
            burden.append(np.sum(dict_probs[gene][sample]))
            probs = np.vstack((probs, dict_probs[gene][sample]))
    normalized_probs = np.mean(probs, axis=0)
    burden_total = np.sum(burden)
    return {(gene, ttype): (normalized_probs, burden_total, len(burden))}


def create_mutrate_dict(gene_ttype_set, cores=10):

    total = dict()
    with Pool(cores) as pool:
        for res in tqdm.tqdm(pool.imap(create_mutrate_dict_func, gene_ttype_set), total=len(gene_ttype_set)):
            total.update(res)
    return total


"""
def create_mutrate_dict(gene_ttype_set):

    d_results = {}
    for gene, ttype in tqdm.tqdm(gene_ttype_set):
        probs = np.zeros(96)
        burden = []
        for cohort in ttype_dict[ttype]:
            path_cohort = os.path.join(conf.mutrate_folder,
                                       cohort + ".mutrate_output",
                                       f"norm_{gene}.out.json")
            try:
                dict_probs = json.load(open(path_cohort, 'r'))
            except FileNotFoundError:
                continue
            for sample in dict_probs[gene].keys():
                burden.append(np.sum(dict_probs[gene][sample]))
                probs = np.vstack((probs, dict_probs[gene][sample]))

        normalized_probs = np.mean(probs, axis=0)
        N = len(burden)
        burden_total = np.sum(burden)
        d_results[(gene, ttype)] = (normalized_probs, burden_total, N)

    return d_results
"""


def format_data(raw_data_instance, mutrate_table, specific=True):

    res = {}
    auc = {}

    for k, v in tqdm.tqdm(raw_data_instance.items()):
        if specific:
            if not v[0]:
                continue
        positive = v[1]['positive']
        negative = v[1]['negative']
        if (len(positive) > 0) and (len(negative) > 0):
            try:
                ttype, gene = k[0], k[1]
                positive_mutrates = set_mutrate(gene, ttype, v[1]['positive'], mutrate_table=mutrate_table)
                negative_mutrates = set_mutrate(gene, ttype, v[1]['negative'], mutrate_table=mutrate_table)
                res[k] = [positive_mutrates, negative_mutrates]

                y_true = [1] * len(positive_mutrates) + [0] * len(negative_mutrates)
                y_score = positive_mutrates + negative_mutrates
                auc[k] = roc_auc_score(y_true, y_score)

            except AssertionError:
                continue

    return res, auc


def reformat(raw):

    res = raw['mutrate']
    auc = raw['auc']
    return res, auc


@click.group()
def cli():
    pass


@cli.command()
@click.option('--driver', is_flag=True)
@click.option('--outfolder', type=click.Path())
def observed(driver, outfolder):
    res = create_observed_vs_non_observed(driver=driver)
    label = 'all'
    if driver:
        label = 'driver'
    foutput = os.path.join(outfolder, f'vectors_observed_{label}.pickle.gz')
    with gzip.open(foutput, 'wb') as f:
        pickle.dump(res, f)


@cli.command()
@click.option('--nonobserved', is_flag=True)
@click.option('--outfolder', type=click.Path())
def driverpassenger(nonobserved, outfolder):
    res = create_driver_vs_passenger(non_observed=nonobserved)
    label = 'all'
    if nonobserved:
        label = 'nonobserved'
    foutput = os.path.join(outfolder, f'vectors_driver_vs_passenger_{label}.pickle.gz')
    with gzip.open(foutput, 'wb') as f:
        pickle.dump(res, f)


@cli.command()
@click.option('--outfolder', type=click.Path())
def allnondrivers(outfolder):
    res = get_nondriver_gene_ttype()
    with open(os.path.join(outfolder, 'all_non_drivers.json'), 'wt') as f:
        json.dump(res, f)


@cli.command()
@click.option('--outfolder', type=click.Path())
@click.option('--cores', type=int)
def nondrivers(outfolder, cores):
    res = create_non_drivers(cores=cores)
    foutput = os.path.join(outfolder, f'vectors_observed_nondrivers_parallel.pickle.gz')
    with gzip.open(foutput, 'wb') as f:
        pickle.dump(res, f)


@cli.command()
@click.option('--csqn_type', type=str)
def csqn(csqn_type):
    res = create_csqn_type(csqn_type=csqn_type)
    foutput = os.path.join(conf.datasets_output_folder, f'vectors_{csqn_type}.pickle.gz')
    with gzip.open(foutput, 'wb') as f:
        pickle.dump(res, f)


@cli.command()
@click.option('--shap', is_flag=True)
def feature(shap):
    threshold = 0
    for feat in tqdm.tqdm(conf.shap_features):
        if not shap:
            feat = '_'.join(feat.split('_')[1:])
            threshold = conf.feat_threshold_dict.get(feat, None)
        if threshold is None:
            continue
        print('Processing...' + feat)
        try:
            res = create_feature(feat, threshold=threshold)
            foutput = os.path.join(conf.datasets_output_folder, f'vectors_{feat}.pickle.gz')
            with gzip.open(foutput, 'wb') as f:
                pickle.dump(res, f)
        except:
            print(feat)


@cli.command()
def phylop():
    for threshold in tqdm.tqdm(np.linspace(1, 10, 10)):
        try:
            res = create_feature('PhyloP', threshold=threshold)
            foutput = os.path.join(conf.datasets_output_folder, f'vectors_PhyloP_{int(threshold)}.pickle.gz')
            with gzip.open(foutput, 'wb') as f:
                pickle.dump(res, f)
        except:
            print(threshold)


@cli.command()
@click.option('--driver', is_flag=True)
@click.option('--outfolder', type=click.Path())
@click.option('--cores', type=int)
def mutrate(driver, outfolder, cores):
    if driver:
        res = create_mutrate_dict(driver_gene_ttypes, cores=cores)
        label = 'driver'
    else:
        with open('./datasets/all_non_drivers.json', 'rt') as f:
            all_non_drivers = json.load(f)
        res = create_mutrate_dict(all_non_drivers, cores=cores)
        label = 'nondriver'
    foutput = os.path.join(outfolder, f'dictionary_mutrate_{label}.pickle.gz')
    with gzip.open(foutput, 'wb') as f:
        pickle.dump(res, f)


@cli.command()
@click.option('--outfolder', type=click.Path())
def saturation(outfolder):
    specific_saturation_folder = os.path.join(conf.output_boostdm, 'saturation', 'prediction')
    res = create_saturation_vectors(specific_saturation_folder)
    foutput = os.path.join(outfolder, 'saturation_vectors_specific.pickle.gz')
    with gzip.open(foutput, 'wb') as f:
        pickle.dump(res, f)


if __name__ == '__main__':

    cli()
