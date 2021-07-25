import os
import sys
import glob
import re
import tqdm
import numpy as np
import itertools
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc, matthews_corrcoef
import json
from matplotlib.colors import to_hex

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
        
import numpy as np
from sklearn.metrics import matthews_corrcoef, roc_auc_score, log_loss


    os.environ['INTOGEN_DATASETS'] = "/workspace/projects/intogen_2017/pipeline/datasets/hg38_vep92_v20191009"
    intogen_data = os.environ['INTOGEN_DATASETS']

    os.environ['DRIVERS_PATH'] = "/workspace/projects/intogen_2017/postprocess/pipeline/20200213/drivers.tsv"
    drivers_path = os.environ['DRIVERS_PATH']

import matplotlib as mpl    

def config_params(font_size=7):

    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'
    mpl.rcParams['lines.linewidth'] = 0.7
    

#'/workspace/datasets/boostdm_runs/20200205/saturation_prediction/*.prediction.out.gz'
def gene_specific_and_role_dict(path_saturation_prediction, drivers_path):

    # dict of tumor-types modelled per gene
    gene_specific_models = {}
    for fn in glob.glob("{}/*.prediction.tsv.gz".format(path_saturation_prediction)):
        gene, ttype = tuple(os.path.basename(fn).split('.')[:2])
        gene_specific_models[gene] = gene_specific_models.get(gene, []) + [ttype]

        # dict of MoA
    df = pd.read_csv(drivers_path, sep='\t')
    role_dict = dict(zip(df.SYMBOL.values, df.ROLE.values))

    return gene_specific_models, role_dict


def get_tissue_specific_muts(df):  
    all_scores_cols = [i for i in df.columns if ('_score_' in i) & (i not in ['boostDM_score_SOLID', 
                                                                           'boostDM_score_NON_SOLID', 
                                                                           'boostDM_score_CANCER'])]
    list_good_mutations = df[all_scores_cols].dropna(how='all').index.tolist()
    df = df.loc[list_good_mutations]
    
    return df

def dropna_from_list(l):
    
    return [s for s in l if not np.isnan(s)]

def custom_argsort(l, key=None, ascending=False):
    
    sortable = list(map(key, l))
    index = np.argsort(sortable)
    if not ascending:
        return index[::-1]
    else:
        return index

def corrected_proportion_potential_drivers(x):
    
    return (sum(np.array(x) >= 0.5) + 1) / (len(x) + 2)
# automatize

def plot_contingency(confusion_summary, csqn_type, moa, keyword, role_dict, gene_specific_models, outpath):

    d = {g: v for g, v in confusion_summary.items() if role_dict.get(g, moa) == moa}

    genes = [k for k in d.keys() if not k.startswith('aggregate')]
    
    fscores = np.array([d[g]['fscore'] for g in genes])  # compute accuracy first, as it gives the sorting criterion
    
    recall = np.array([d[g]['recall'] for g in genes])
    precision = np.array([d[g]['precision'] for g in genes])
    true_boostdms = np.array([d[g]['true_boostdm'] for g in genes])
    false_boostdms = np.array([d[g]['false_boostdm'] for g in genes])
    fscores = np.nan_to_num(fscores)
    sortindex = np.argsort(fscores)[::-1]

    fscores = fscores[sortindex]
    recall = recall[sortindex]
    precision = precision[sortindex]
    true_boostdms = true_boostdms[sortindex]
    false_boostdms = false_boostdms[sortindex]

    genes = list(np.array(genes)[sortindex])  # sorted genes, it does not include summary track yet

    # summary track
    fscores_summary = d[f'aggregate_{moa}']['fscore']
    recall_summary = d[f'aggregate_{moa}']['recall']
    precision_summary = d[f'aggregate_{moa}']['precision']
    true_boostdms_summary = d[f'aggregate_{moa}']['true_boostdm']
    false_boostdms_summary = d[f'aggregate_{moa}']['false_boostdm']
    
    
    # beginning of the plot
    
    unit = 0.15
    height = 2.6
    width  = len(genes) * unit + 1.3
    
    fig, ax = plt.subplots(figsize=(width, height), nrows=5, gridspec_kw={
                           'height_ratios': [0.75, 0.75, 0.5, 0.5, 0.5]})
    
    # track 1: TRUE boostDM score

    true_boostdms = list(true_boostdms) + [true_boostdms_summary]
    
    x = list(itertools.chain.from_iterable([[0 + i] * len(l) for i, l in enumerate(true_boostdms)]))
    x += np.random.normal(0, 0.05, size=len(x))
    y = list(itertools.chain.from_iterable(true_boostdms))
    

    ax[0]=sns.boxplot(data=true_boostdms, showfliers=False, color = 'white',
                medianprops={'color':'black'}, linewidth = 0.4, ax = ax[0],) #zorder=0
                      #boxprops=dict(alpha=0.5), 
    
    ax[0].scatter(x, y, s=3, alpha=0.25,
                  color='DarkRed',  )  # scatter zorder=5
                   

    ax[0].set_xticklabels([])
    ax[0].set_xticks([])
    #ax[0].set_title(f'{keyword}: {csqn_type}: {moa}')
    ax[0].set_ylim(-0.1, 1.1)
    ax[0].set_xlim(-0.5, len(true_boostdms)+0.5 )
    ax[0].set_ylabel('boostDM\ntrue', rotation=0, labelpad=20)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['bottom'].set_visible(False)
    
    ax[0].hlines(0.5, -0.5, len(true_boostdms)+0.5, ls = '--',lw = 0.5, color = 'grey')

    
    # track 2: FALSE boostDM score
    false_boostdms = list(false_boostdms) + [false_boostdms_summary]
    
    x = list(itertools.chain.from_iterable([[0 + i] * len(l) for i, l in enumerate(false_boostdms)]))
    x += np.random.normal(0, 0.05, size=len(x))
    y = list(itertools.chain.from_iterable(false_boostdms))
    
    
    ax[1]=sns.boxplot(data=false_boostdms, showfliers=False, color = 'white',
                medianprops={'color':'black'}, linewidth = 0.4, ax = ax[1])
    
    ax[1].scatter(x, y, s=3, alpha=0.25,
                  color='grey')  # scatter

    ax[1].set_xticks([])
    ax[1].set_ylim(-0.1, 1.1)
    ax[1].set_xlim(-0.5, len(false_boostdms)+0.5 )

    ax[1].set_ylabel('boostDM\nfalse', rotation=0, labelpad=20)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)


    ax[1].hlines(0.5, -0.5,  len(false_boostdms)+0.5, ls = '--',lw = 0.5, color = 'grey')

    
    # track 3: fscores
    fscores = list(fscores) + [fscores_summary]
    
    sns.barplot(x=list(range(len(genes) + 1)), y=fscores, color='#FFBD29', ax=ax[2], alpha=0.7)


    ax[2].set_xticklabels([])
    ax[2].set_xticks([])
    ax[2].set_ylim(0, 1.1)
    ax[2].set_xlim(-0.5, len(fscores)+0.5 )

    ax[2].set_ylabel('F-score50', rotation=0, labelpad=15)
    ax[2].spines['right'].set_visible(False)
    ax[2].spines['top'].set_visible(False)
    ax[2].spines['bottom'].set_visible(False)

    for v, x in enumerate(fscores):
        ax[2].text(v-0.25, x+0.1, round(x, 2) , fontsize = 3)

    # track 4: recall
    recall = list(recall) + [recall_summary]

    sns.barplot(x=list(range(len(genes) + 1)), y=recall, color='DarkRed', ax=ax[3], alpha=0.7)
    ax[3].set_xticklabels([])
    ax[3].set_xticks([])
    ax[3].set_ylim(0, 1)
    ax[3].set_xlim(-0.5, len(recall)+0.5 )

    ax[3].set_ylabel('Recall', rotation=0, labelpad=15)
    ax[3].spines['right'].set_visible(False)
    ax[3].spines['top'].set_visible(False)
    ax[3].spines['bottom'].set_visible(False)

    #ax[3].hlines(0.5, 0.5, len(true_boostdms), ls = '--',lw = 0.5, color = 'grey')

    # track 5: precision
    
    precision = list(precision) + [precision_summary]

    sns.barplot(x=list(range(len(genes) + 1)), y=precision, color='#548F25', ax=ax[4], alpha=0.7)
    ax[4].set_ylim(0, 1)
    ax[4].set_xlim(-0.5, len(precision)+0.5)

    ax[4].set_ylabel('Precision', rotation=0, labelpad=15)
    ax[4].spines['right'].set_visible(False)
    ax[4].spines['top'].set_visible(False)    
    #ax[4].hlines(0.5, 0.5, len(true_boostdms), ls = '--',lw = 0.5, color = 'grey')

    # track 6: no. mutations
    ax[4].set_xticks(np.arange(len(genes) + 1))

    tolab = []
    for i1, i2, g in zip(true_boostdms, false_boostdms, genes + ['aggregate']):
        tolab.append("{} ({}|{})".format(g, len(i1), len(i2)))


    ax[4].set_xticklabels(tolab, rotation=90)


    plt.savefig(f'{outpath}/comparison.genewise.{keyword}.{moa}.{csqn_type}.validation.png', dpi=200, bbox_inches='tight')
    plt.savefig(f'{outpath}/comparison.genewise.{keyword}.{moa}.{csqn_type}.validation.svg')
    
    return  fscores_summary, recall_summary, precision_summary, true_boostdms_summary, false_boostdms_summary


def plot_summary_contingency( fscore, recall, precision, true_boostdms, false_boostdms, keyword, moa, outpath, csqn_type):

    # summary plot
    
    fig, ax = plt.subplots(figsize=(2.5, 1.25), nrows=1, ncols=2, 
                          gridspec_kw={
                           'width_ratios': [1, 0.5]})

    # track 1: boxplot

    x = list(itertools.chain.from_iterable([[0 + i] * len(l) for i, l in enumerate((true_boostdms, false_boostdms))]))
    x += np.random.normal(0, 0.05, size=len(x))
    cols = ['DarkRed'] * len(true_boostdms) + ['grey'] * len(false_boostdms)
    y = list(true_boostdms) + list(false_boostdms)

    alpha = 0.25
    if  "Clinvar" in keyword :
        alpha = 0.05

    ax[0].scatter(x, y, s=4, alpha=alpha, color=cols)  # scatter
        
    boxprops = dict(linestyle='-', linewidth=2, color='grey')
    medianprops = dict(linestyle='-', linewidth=2, color='black')
    capsprops = boxprops
    ax[0]=sns.boxplot(data=[true_boostdms, false_boostdms], showfliers=False, color = 'white',
            medianprops={'color':'black'}, linewidth = 0.7, ax = ax[0])
        

    ax[0].hlines(0.5, -0.5,  1.5, ls = '--',lw = 0.5, color = 'grey')
    


    if "Clinvar" in keyword:
        #ax[0].set_xticklabels(['true  n ={}'.format(len(true_boostdms)), 'false n={}'.format(len(false_boostdms))], rotation = 90)
        ax[0].set_xticklabels(["pathogenic", "benign"])
        ax[0].text(-0.25, 1.05, 'n ={}'.format(len(true_boostdms)) )
        ax[0].text( 0.75, 1.05, 'n ={}'.format(len(false_boostdms)) )

        #ax[0].set_title(f'{keyword}: {csqn_type}: {moa}')
    else:
        ax[0].set_xticklabels(["driver", "passenger"])
        #ax[0].set_title(f'{keyword}: {csqn_type}: {moa}')
        ax[0].text(-0.25, 1.05, 'n ={}'.format(len(true_boostdms)) )
        ax[0].text( 0.75, 1.05, 'n ={}'.format(len(false_boostdms)) )

    ax[0].set_ylabel('boostDM')
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)

    # track 2: accuracies

    bars = ( fscore, recall, precision)
    ax[1].bar(x=range(len(bars)), height=bars, 
              color=list(map(to_hex, ('#FFBD29', 'darkred', '#548F25'))))
    
    for ix, h in enumerate(bars):
        ax[1].text(ix-0.45, h+0.01, round(h, 2), fontsize = 5) 
    
    ax[1].set_xticks(range(len(bars)))
    ax[1].set_xticklabels(labels=['F-score50', 'Recall', 'Precision'], rotation = 90)
    ax[1].set_ylim(0, 1)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)

    plt.subplots_adjust(wspace = 0.75)

    plt.show()

    #plt.tight_layout()
    plt.savefig(f'{outpath}/comparison.summary.{keyword}.{moa}.{csqn_type}.validation.png', dpi=200, bbox_inches='tight')
    plt.savefig(f'{outpath}/comparison.summary.{keyword}.{moa}.{csqn_type}.validation.svg', dpi=200, bbox_inches='tight')

    


def get_confusion_summary(cross_comparison_dict, cross_comparison_dict_full, list_MOA):

    confusion_summary = {}
    
    for g in cross_comparison_dict[True]:
        
        confusion_summary[g] = {}

        true_boostdm  = cross_comparison_dict[True][g]
        false_boostdm = cross_comparison_dict[False][g]

        tp = sum(np.array(cross_comparison_dict[True][g]) >= 0.5)
        tn = sum(np.array(cross_comparison_dict[False][g]) < 0.5)
        fp = sum(np.array(cross_comparison_dict[False][g]) >= 0.5)
        fn = sum(np.array(cross_comparison_dict[True][g]) < 0.5)

        confusion_summary[g]['tp'] = tp
        confusion_summary[g]['tn'] = tn
        confusion_summary[g]['fp'] = fp
        confusion_summary[g]['fn'] = fn
        recall = 0
        precision = 0
        fscore = 0

        if (tp >0)|((fn>0) & (fp > 0)):
            ppv = (tp + 1) / (tp + fp + 2)
            npv = (tn + 1) / (tn + fn + 2)
            acc = (tp + tn + 1) / (tp + tn + fp + fn + 2)
            mcc = ((tp * tn) - (fp * fn)) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
            recall = tp / (tp + fn)
            precision = tp / (tp + fp)

            beta = 0.5
            fscore = (1 + beta ** 2) * precision * recall * (1 / (precision * beta ** 2 + recall))
        
        confusion_summary[g]['recall'] = recall
        confusion_summary[g]['precision'] = precision
        confusion_summary[g]['fscore'] = fscore
        confusion_summary[g]['true_boostdm'] = true_boostdm
        confusion_summary[g]['false_boostdm'] = false_boostdm

    for moa in [list_MOA]:

        confusion_summary[f'aggregate_{moa}'] = {}

        tp = 0
        tn = 0
        fp = 0
        fn = 0
        true_boostdm = []
        false_boostdm = []

        # get all genes with True values
        for g in cross_comparison_dict_full[True]:

            true_boostdm.extend(cross_comparison_dict_full[True][g])
            tp+=sum(np.array(cross_comparison_dict_full[True][g]) >= 0.5)
            fn+= sum(np.array(cross_comparison_dict_full[True][g]) < 0.5)

        # get all genes with False values
        for g in cross_comparison_dict_full[False]:

            false_boostdm.extend(cross_comparison_dict_full[False][g])
            tn+=sum(np.array(cross_comparison_dict_full[False][g]) < 0.5)
            fp+= sum(np.array(cross_comparison_dict_full[False][g]) >= 0.5)

        confusion_summary[f'aggregate_{moa}']['tp'] = tp
        confusion_summary[f'aggregate_{moa}']['tn'] = tn
        confusion_summary[f'aggregate_{moa}']['fp'] = fp
        confusion_summary[f'aggregate_{moa}']['fn'] = fn
        
        recall = 0
        precision = 0
        fscore = 0

        if (tp >0)|((fn>0) & (fp > 0)):
            ppv = (tp + 1) / (tp + fp + 2)
            npv = (tn + 1) / (tn + fn + 2)
            acc = (tp + tn + 1) / (tp + tn + fp + fn + 2)
            mcc = ((tp * tn) - (fp * fn)) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
            recall = tp / (tp + fn)
            precision = tp / (tp + fp)

            beta = 0.5
            fscore = (1 + beta ** 2) * precision * recall * (1 / (precision * beta ** 2 + recall))

        confusion_summary[f'aggregate_{moa}']['recall'] = recall
        confusion_summary[f'aggregate_{moa}']['precision'] = precision
        confusion_summary[f'aggregate_{moa}']['fscore'] = fscore
        confusion_summary[f'aggregate_{moa}']['true_boostdm'] = true_boostdm
        confusion_summary[f'aggregate_{moa}']['false_boostdm'] = false_boostdm
        
    return confusion_summary


def get_cross_comparison(enrichment_dict, term_true, term_false):
    
    cross_comparison_dict = {}
    cross_comparison_dict_full = {}


    try:
        assert(term_true in enrichment_dict)
        assert(term_false in enrichment_dict)

        # these both contain the full information
        d_true1  = enrichment_dict[term_true]
        d_false1 = enrichment_dict[term_false]

        common_genes = set(d_true1.keys()).intersection(set(d_false1.keys()))
        d_true = {k: v for k, v in d_true1.items() if k in common_genes}
        d_false = {k: v for k, v in d_false1.items() if k in common_genes}

        cross_comparison_dict = {True: d_true, False: d_false}
        cross_comparison_dict_full =  {True: d_true1, False: d_false1}

    except AssertionError as e:
        print('Failed')

        
    return cross_comparison_dict, cross_comparison_dict_full


# create a dictionary with instances 
def enrichment(df, csqn_type, score_terms):
    
    # get all values fun score
        
    enrichment_dict = {}
    
    for term in score_terms:

        enrichment_dict[term] = {}
        # compute collection of values
        
        # select all cases with the function and the term
        dg = df[df[f'csqn_type_{csqn_type}'] == 1]
        dg = dg.loc[lambda dg: dg['FUNC_SCORE'].str.contains(term)]
        if len(dg) == 0:
            continue

        listify = lambda r: (r['SYMBOL'], r.loc[[s for s in df.columns.tolist() 
                                                 if 'boostDM_score' in s]].values.tolist())
        
        values = dg.apply(listify, axis=1).tolist()

        genes = dg['SYMBOL'].unique()
        genes_dict = {}
        for g, l in values:
            genes_dict[g] = genes_dict.get(g, []) + dropna_from_list(l)

        # dump into dictionary
        enrichment_dict[term] = genes_dict
    
    return enrichment_dict

# main sender
def plot_boostdm_validation_bias(df, true_label, false_label, csqn_type, moa, key, role_dict, gene_specific, outpath):
    
    # remove synonymous and Nans in the score column
    df.dropna(subset=['FUNC_SCORE'], inplace=True)

    df["FUNC_SCORE"] = df["FUNC_SCORE"].astype(str)
    
    values_func_score = [re.split(',|/', effect) for effect in df.FUNC_SCORE.unique()]
    score_terms = set(list(itertools.chain.from_iterable(values_func_score)))

    df = df[df['csqn_type_synonymous'] == 0]  # filter out synonymous
    genes_wanted = [k for k, v in role_dict.items() if v == moa]
    df = df[df['SYMBOL'].isin(genes_wanted)]

    # boostDM saturations have only run for some genes, so we will filter them here
    all_sat = list(gene_specific.keys())

    df = df[df['SYMBOL'].isin(all_sat)]
        
    # create dictionary with all functional cases
    enrichment_dict = enrichment(df, csqn_type, score_terms)


    # cross comparison 
    cross_comparison_dict, cross_comparison_dict_full = get_cross_comparison(enrichment_dict, true_label,
                                                                             false_label) 

    
    # confusion matrix
    confusion_summary = get_confusion_summary(cross_comparison_dict, cross_comparison_dict_full, moa)
    keyword = f'{key}.{true_label}.{false_label}'
    
    fscore, recall, precision, true_boostdms, false_boostdms = plot_contingency(confusion_summary, csqn_type, moa, keyword, role_dict, gene_specific, outpath)
    
    plot_summary_contingency( fscore, recall, precision, true_boostdms, false_boostdms, keyword, moa, outpath, csqn_type)

    return  fscore, recall, precision


def comparison_main(input, outpath, true_label, false_label, csqn_type, label,  predictions):

    # if the holdout mode is on, we should 
    config_params(font_size=7)
    
    df = pd.read_csv(input, sep='\t')
    path_saturation_prediction = predictions  #"/workspace/datasets/boostdm_runs/boostdm-new-pipeline-discovery/output/saturation/prediction/"

    gene_specific, role_dict = gene_specific_and_role_dict(path_saturation_prediction, drivers_path)

    res = {}
    for moa in ['LoF', 'Act']:    

        print("Doing {}...".format(moa))
        
        acc, ppv, npv = plot_boostdm_validation_bias(df, true_label, false_label, csqn_type, moa, label, role_dict, gene_specific, outpath)
        res[moa] = (acc, ppv, npv)

        #except:
        #    print(moa, "not working")
    json.dump(res, open(outpath + "/" + label + ".benchmark.json", "wt"))

