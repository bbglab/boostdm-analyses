import os
import numpy as np
from scipy.stats import mannwhitneyu

import seaborn as sns
from matplotlib import gridspec
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

from conf import dict_colors_role, config_params


config_params()


def is_specific(gene, ttype, model_dict):
    """
    asserts whether the model employed for gene-ttype
    is specific -- has been fitted with mutations 
    strictly from the gene-ttype -- or not.

    it uses model_dict, an instance externally defined
    with the boostdm models
    """

    model = model_dict.get((ttype, gene), None)
    return (ttype, gene) == model


def plot_auc_vs_observed(res, auc, obs, pairs, model_dict, df_drivers, title='observed_vs_rest_non_drivers', specific=False, highlight_genes=None, density=False, saveplot=False):

    # items to highlight in the plot

    
    highlight_gene_ttypes_coord = {}
    highlight_gene_ttypes_sizes = {}
    x, y, s, c = [], [], [], []
    for k in pairs:
        gene = k[1]
        ttype = k[0]
        if specific and not is_specific(gene, ttype, model_dict):
            continue
        elif auc[k] is not None:
            ones, zeros = res[k]
            y.append(obs[k])
            x.append(auc[k])
            s.append(len(ones)+1)
            moa = df_drivers[df_drivers['SYMBOL'] == gene]['ROLE'].values[0]
            c.append(dict_colors_role[moa])
            if highlight_genes is not None:
                if (gene, ttype) in highlight_genes:
                    highlight_gene_ttypes_coord[(gene, ttype)] = x[-1], y[-1]
                    highlight_gene_ttypes_sizes[(gene, ttype)] = len(ones)
                    # print(f'{k[0]}: {k[1]}: auc: {y[-1]} log_fold_change {x[-1]} size: {len(ones)}')

    fig, ax = plt.subplots(figsize=(5.5, 5.5))

    gs = gridspec.GridSpec(figure=fig, ncols=2, nrows=2, width_ratios=[15,2], height_ratios=[2, 15])
    gs.update(hspace=0.0, wspace=0.00)

    ax0 = plt.subplot(gs[0]) # density top
    ax1 = plt.subplot(gs[1]) # null
    ax2 = plt.subplot(gs[2], sharex=ax0) # scatter
    ax3 = plt.subplot(gs[3], sharey=ax2) # density


    # ax0: density plot: oncogenes vs tumor suppressors

    if density:

        x_lof = [u for i, u in enumerate(x) if c[i] == dict_colors_role['LoF']]
        x_act = [u for i, u in enumerate(x) if c[i] == dict_colors_role['Act']]
        
        y_lof = [u for i, u in enumerate(y) if c[i] == dict_colors_role['LoF']]
        y_act = [u for i, u in enumerate(y) if c[i] == dict_colors_role['Act']]

        bandwidth = 0.03
        alpha = 0.25
        sns.distplot(x_lof, hist=False, kde=True, color=dict_colors_role['LoF'], 
                     kde_kws={'linewidth': 3, 'bw': bandwidth, 'alpha': alpha}, 
                     ax=ax0, vertical=False)
        sns.distplot(x_act, hist=False, kde=True, color=dict_colors_role['Act'], 
                     kde_kws={'linewidth': 3, 'bw': bandwidth, 'alpha': alpha}, 
                     ax=ax0, vertical=False)
        sns.distplot(y_lof, hist=False, kde=True, color=dict_colors_role['LoF'], 
                     kde_kws={'linewidth': 3, 'bw': bandwidth, 'alpha': alpha}, 
                     ax=ax3, vertical=True)
        sns.distplot(y_act, hist=False, kde=True, color=dict_colors_role['Act'], 
                     kde_kws={'linewidth': 3, 'bw': bandwidth, 'alpha': alpha}, 
                     ax=ax3, vertical=True)

    else:

        c = ['grey'] * len(c)

    # ax1: scatter plot

    ax2.scatter(x, y, alpha=0.25, s=s, c=c)
    ax2.vlines(0.5, -0.5, 2, color='grey', linestyles='dashed' , lw= 1.0)
    ax2.hlines(0, 0.3, 1., color='grey', linestyles='dashed', lw = 1.0)
    ax2.set_ylabel('observed to potential drivers',fontsize=12)
    ax2.set_xlabel('mutation probability bias\n (auROC)',fontsize=12)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlim(0.35, 1.0)
    ax2.set_ylim(0, 0.7)

    if highlight_genes is not None:
        for (gene, ttype), v in highlight_gene_ttypes_coord.items():
            ax2.text(v[0] + np.random.uniform(low=0, high=0.05, size=1), v[1], 
                     f'{gene} (ttype)',fontsize=10)
            ax2.scatter([v[0]], [v[1]], 
                        s=highlight_gene_ttypes_sizes[(gene, ttype)],
                        marker='o', color='white', edgecolors='black')

    ax0.axis('off')
    ax1.axis('off')
    ax3.axis('off')

    if density:
        ax0.set_title(title)
    else:
        ax1.set_title(title)

    if saveplot:
        plt.savefig(f'./plots/{title}.pdf', dpi=800,bbox_inches="tight")
        plt.savefig(f'./plots/{title}.png', dpi=800,bbox_inches="tight")

    plt.show()

    return y

def plot_auc_vs_observed_kde(res, auc, obs, pairs, model_dict, df_drivers, title='observed_vs_rest_non_drivers', specific=False, highlight_genes=None, density=False, saveplot=False):

    # items to highlight in the plot

    
    highlight_gene_ttypes_coord = {}
    highlight_gene_ttypes_sizes = {}
    # print(highlight_genes)
    x, y, s, c = [], [], [], []
    for k in pairs:
        gene = k[1]
        ttype = k[0]
        if specific and not is_specific(gene, ttype, model_dict):
            continue
        elif auc[k] is not None:
            ones, zeros = res[k]
            y.append(obs[k])
            x.append(auc[k])
            s.append(len(ones)+1)
            moa = df_drivers[df_drivers['SYMBOL'] == gene]['ROLE'].values[0]
            c.append(dict_colors_role[moa])
            if highlight_genes is not None:
                
                if k in highlight_genes:
                    # print (k)
                    highlight_gene_ttypes_coord[k] = x[-1], y[-1]
                    highlight_gene_ttypes_sizes[k] = len(ones)
                    # print(f'{k[0]}: {k[1]}: auc: {y[-1]} log_fold_change {x[-1]} size: {len(ones)}')

    fig, ax = plt.subplots(figsize=(5.5, 5.5))

    gs = gridspec.GridSpec(figure=fig, ncols=2, nrows=2, width_ratios=[15,2], height_ratios=[2, 15])
    gs.update(hspace=0.0, wspace=0.00)

    ax0 = plt.subplot(gs[0]) # density top
    ax1 = plt.subplot(gs[1]) # null
    ax2 = plt.subplot(gs[2], sharex=ax0) # scatter
    ax3 = plt.subplot(gs[3], sharey=ax2) # density


    # ax0: density plot: oncogenes vs tumor suppressors

    if density:

        x_lof = [u for i, u in enumerate(x) if c[i] == dict_colors_role['LoF']]
        x_act = [u for i, u in enumerate(x) if c[i] == dict_colors_role['Act']]
        
        y_lof = [u for i, u in enumerate(y) if c[i] == dict_colors_role['LoF']]
        y_act = [u for i, u in enumerate(y) if c[i] == dict_colors_role['Act']]

        bandwidth = 0.03
        alpha = 0.25
        sns.distplot(x_lof, hist=False, kde=True, color=dict_colors_role['LoF'], 
                     kde_kws={'linewidth': 3, 'bw': bandwidth, 'alpha': alpha}, 
                     ax=ax0, vertical=False)
        sns.distplot(x_act, hist=False, kde=True, color=dict_colors_role['Act'], 
                     kde_kws={'linewidth': 3, 'bw': bandwidth, 'alpha': alpha}, 
                     ax=ax0, vertical=False)
        sns.distplot(y_lof, hist=False, kde=True, color=dict_colors_role['LoF'], 
                     kde_kws={'linewidth': 3, 'bw': bandwidth, 'alpha': alpha}, 
                     ax=ax3, vertical=True)
        sns.distplot(y_act, hist=False, kde=True, color=dict_colors_role['Act'], 
                     kde_kws={'linewidth': 3, 'bw': bandwidth, 'alpha': alpha}, 
                     ax=ax3, vertical=True)

    else:

        c = ['grey'] * len(c)

    # ax1: scatter plot

    sns.kdeplot(x_act,y_act,  shade=True, shade_lowest=False,n_levels=5,ax=ax2,kde_kws={"alpha":0.3},cmap=sns.light_palette(color=dict_colors_role["Act"], input="hex", as_cmap=True))
    sns.kdeplot(x_lof,y_lof,  shade=True, shade_lowest=False,n_levels=5,ax=ax2,kde_kws={"alpha":0.3},cmap=sns.light_palette(color=dict_colors_role["LoF"], input="hex", as_cmap=True))
    ax2.scatter(x, y, alpha=0.8, s=s, c=c)
    
    ax2.vlines(0.5, -0.5, 2, color='grey', linestyles='dashed' , lw= 1.0)
    ax2.hlines(0, 0.3, 1., color='grey', linestyles='dashed', lw = 1.0)
    ax2.set_ylabel('observed to potential drivers',fontsize=12)
    ax2.set_xlabel('mutation probability bias\n (auROC)',fontsize=12)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlim(0.35, 1.0)
    ax2.set_ylim(0, 0.7)

    if highlight_genes is not None:
        for k, v in highlight_gene_ttypes_coord.items():
            ax2.text(v[0] + np.random.uniform(low=0, high=0.05, size=1), v[1], 
                     f'{k[1]} ({k[0]})',fontsize=10)
            ax2.scatter([v[0]], [v[1]], 
                        s=highlight_gene_ttypes_sizes[k], 
                        marker='o', color='white', edgecolors='black')

    ax0.axis('off')
    ax1.axis('off')
    ax3.axis('off')

    if density:
        ax0.set_title(title)
    else:
        ax1.set_title(title)

    if saveplot:
        plt.savefig(f'./plots/{title}.pdf', dpi=800, bbox_inches="tight")
        plt.savefig(f'./plots/{title}.png', dpi=800, bbox_inches="tight")

    plt.show()

    return y


def auc_randomize(res, N=100):
    random = {}
    for k in res:
        pos, neg = tuple(res[k])
        s = pos + neg
        y_true = [1] * len(pos) + [0] * len(neg)
        aucs = []
        for _ in range(N):
            np.random.shuffle(s)
            s_pos = s[:len(pos)+1]
            s_neg = s[len(pos)+1:]
            aucs.append(roc_auc_score(y_true, s_pos + s_neg))
        random[k] = np.array(aucs)
    return random


def plot_auc_volcano(res, auc, model_dict, df_drivers, output_folder,
                     title='title', specific=False, highlight_genes=None,
                     size=None, density=False, significance=None, saveplot=False, figsize=(10, 5),
                     xlim=(0.4, 0.7), ylim=(-0.2, 0.5)):

    # items to highlight in the plot
    highlight_gene_ttypes_coord = {}
    highlight_gene_ttypes_sizes = {}
    aucroc_significant = {}
    x, y, s, c = [], [], [], []

    if significance is not None:
        random_auc = auc_randomize(res)

    for k, mutrates in res.items():
        ttype = k[0]
        gene = k[1]
        if specific:
            if not is_specific(gene, ttype, model_dict):
                continue
        if auc[k] is not None:
            positives, negatives = tuple(mutrates)
            if len(positives) < 10:
                continue
            logfc = np.log10(np.median(positives) / np.median(negatives))
            auc_value = auc[k]
            x.append(logfc)
            y.append(auc_value)
            s.append(len(positives)+1)
            moa = df_drivers[df_drivers['SYMBOL'] == gene]['ROLE'].values[0]
            c.append(dict_colors_role[moa])
            if significance is not None:
                q_up = np.quantile(random_auc[k], 1 - (significance / 2))
                q_down = np.quantile(random_auc[k], significance / 2)
                if (y[-1] < q_down) or (y[-1] > q_up):
                    if (x[-1] < -0.1) or (x[-1] > 0.1):
                        aucroc_significant[(gene, ttype)] = logfc, auc_value
            if highlight_genes is not None:
                if (gene, ttype) in highlight_genes:
                    highlight_gene_ttypes_coord[(gene, ttype)] = logfc, auc_value
                    highlight_gene_ttypes_sizes[(gene, ttype)] = len(positives)
                    # print(f'{gene}: {ttype}: auc: {auc_value} log_fold_change {logfc} size: {len(positives)}, {len(negatives)}')

    fig, ax = plt.subplots(figsize=figsize)
    gs = gridspec.GridSpec(figure=fig, ncols=1, nrows=2,
                           width_ratios=[10], height_ratios=[2, 15])
    gs.update(hspace=0.0, wspace=0.00)

    ax1 = plt.subplot(gs[1])  # scatter
    ax0 = plt.subplot(gs[0], sharex=ax1)  # density

    # ax0: density plot: oncogenes vs tumor suppressors

    if density:
        x_lof = [u for i, u in enumerate(y) if c[i] == dict_colors_role['LoF']]
        x_act = [u for i, u in enumerate(y) if c[i] == dict_colors_role['Act']]
        alpha = 0.25
        sns.distplot(x_lof, hist=False, kde=True, color=dict_colors_role['LoF'],
                     kde_kws={'linewidth': 3, 'bw_method': 'silverman', 'alpha': alpha},
                     ax=ax0, vertical=False)
        sns.distplot(x_act, hist=False, kde=True, color=dict_colors_role['Act'], 
                     kde_kws={'linewidth': 3, 'bw_method': 'silverman', 'alpha': alpha},
                     ax=ax0, vertical=False)
    else:
        c = ['grey'] * len(c)

    if size is not None:
        s = [size] * len(s)

    ax1.scatter(y, x, alpha=0.25, s=s, c=c)
    ax1.vlines(0.5, -0.5, 2, color='grey', linestyles='dashed', lw=1.0)
    ax1.hlines(0, 0.3, 1., color='grey', linestyles='dashed', lw=1.0)
    ax1.set_ylabel('logFC', fontsize=12)
    ax1.set_xlabel('auROC', fontsize=12)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    ax1.set_xlim(xlim[0], xlim[1])
    ax1.set_ylim(ylim[0], ylim[1])

    if significance is not None:
        for (gene, ttype), v in aucroc_significant.items():
            moa = df_drivers[df_drivers['SYMBOL'] == gene]['ROLE'].values[0]
            color = dict_colors_role[moa]
            ax1.scatter([v[1]], [v[0]], s=30, marker='o', color='red')
            ax1.text(v[1] + np.random.uniform(low=0, high=0.01, size=1), v[0],
                     f'{gene} ({ttype})', fontsize=10,
                     color=color, weight='bold')

    if highlight_genes is not None:
        for (gene, ttype), v in highlight_gene_ttypes_coord.items():
            moa = df_drivers[df_drivers['SYMBOL'] == gene]['ROLE'].values[0]
            color = dict_colors_role[moa]
            ax1.text(v[1] + np.random.uniform(low=0, high=0.01, size=1), v[0],
                     f'{gene} ({ttype})', fontsize=10, color=color, weight='bold')
            custom_size = highlight_gene_ttypes_sizes[(gene, ttype)]
            custom_color = color
            if size is not None:
                custom_size = size
                custom_color = 'black'
            ax1.scatter([v[1]], [v[0]], s=custom_size, marker='o', color='white', edgecolors=custom_color)

    ax0.axis('off')

    if density:
        ax0.set_title(title)
    else:
        ax1.set_title(title)

    if saveplot:
        plt.savefig(os.path.join(output_folder, f'{title}.png'), dpi=300, bbox_inches='tight')
        # plt.savefig(os.path.join(output_folder, f'{title}.svg'), dpi=300, bbox_inches='tight')

    plt.show()


def comparison_boxplot(auc1, auc2, legend1, legend2, fname):

    fig, ax = plt.subplots(figsize=(3, 5))

    mycolor = 'grey'
    y = [list(auc1.values()), list(auc2.values())]
    bplot = ax.boxplot(y, vert=True, showfliers=False, widths=[0.5,0.5],
                       patch_artist=True, medianprops={'color': 'black'},
                       boxprops={'color': mycolor},
                       capprops={'color': mycolor},
                       whiskerprops={'color': mycolor})

    p = mannwhitneyu(list(auc1.values()), list(auc2.values()))[1]

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    ax.set_xticks([])
    ax.text(1, 1.1, legend1, horizontalalignment='center')
    ax.text(2, 1.1, legend2, horizontalalignment='center')

    ax.text(1., 0.35, f'N={len(auc1)}', horizontalalignment='center')
    ax.text(2, 0.35, f'N={len(auc2)}', horizontalalignment='center')

    ax.text(1.3, 1, f'p < {p:0.1e}')

    ax.set_ylabel('Probability Bias')

    ax.hlines(0.5, 0.5, 2.5, linestyles='--', linewidth=2, alpha=0.7, color='red')

    colors = ['grey', 'grey']
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

    for box in bplot['boxes']:
        box.set_alpha(0.7)

    output_folder = './raw_plots/'
    plt.savefig(os.path.join(output_folder, f'{fname}-boxplot.svg'))

    plt.show()
