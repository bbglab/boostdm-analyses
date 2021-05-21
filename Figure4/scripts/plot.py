import os
import numpy as np
from scipy.stats import mannwhitneyu

import seaborn as sns
from matplotlib import gridspec
import matplotlib.pyplot as plt

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


def plot_auc_volcano(res, auc, model_dict, df_drivers, output_folder,
                     title='', highlight=None, density=False,
                     figsize=(10, 5), xlim=(0.4, 0.7), ylim=(-0.2, 0.5)):

    # items to highlight in the plot
    highlight_gene_ttypes_coord = {}
    highlight_gene_ttypes_sizes = {}

    x, y, s, c = [], [], [], []

    for k, mutrates in res.items():
        ttype = k[0]
        gene = k[1]
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

            if highlight is not None:
                if (gene, ttype) in highlight:
                    highlight_gene_ttypes_coord[(gene, ttype)] = logfc, auc_value
                    highlight_gene_ttypes_sizes[(gene, ttype)] = len(positives)

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
                     kde_kws={'linewidth': 3, 'alpha': alpha},
                     ax=ax0, vertical=False)
        sns.distplot(x_act, hist=False, kde=True, color=dict_colors_role['Act'], 
                     kde_kws={'linewidth': 3, 'alpha': alpha},
                     ax=ax0, vertical=False)
    else:
        c = ['grey'] * len(c)

    ax1.scatter(y, x, alpha=0.25, s=s, c=c)
    ax1.vlines(0.5, -0.5, 2, color='grey', linestyles='dashed', lw=1.0)
    ax1.hlines(0, 0.3, 1., color='grey', linestyles='dashed', lw=1.0)
    ax1.set_ylabel('logFC', fontsize=12)
    ax1.set_xlabel('Probability Bias', fontsize=12)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    ax1.set_xlim(xlim[0], xlim[1])
    ax1.set_ylim(ylim[0], ylim[1])

    if highlight is not None:
        for (gene, ttype), v in highlight_gene_ttypes_coord.items():
            moa = df_drivers[df_drivers['SYMBOL'] == gene]['ROLE'].values[0]
            color = dict_colors_role[moa]
            ax1.text(v[1] + np.random.uniform(low=0, high=0.01, size=1), v[0],
                     f'{gene} ({ttype})', fontsize=10, color=color, weight='bold')
            custom_size = highlight_gene_ttypes_sizes[(gene, ttype)]
            custom_color = color
            ax1.scatter([v[1]], [v[0]], s=custom_size, marker='o', color='white', edgecolors=custom_color)

    ax0.axis('off')

    if density:
        ax0.set_title(title)
    else:
        ax1.set_title(title)

    plt.savefig(os.path.join(output_folder, f'volcano.png'))

    plt.show()
