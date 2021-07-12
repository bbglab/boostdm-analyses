from collections import defaultdict
import numpy as np
import pandas as pd
import scipy.stats
from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import cm

from conf import *

cmap = cm.get_cmap('tab10')
colors = [cmap(0), cmap(0), cmap(1), cmap(1), cmap(2), cmap(2), cmap(3), cmap(4), cmap(6), cmap(7)]
alphas = [1., 0.7, 1., 0.7, 1., 0.7, 1., 1., 1., 1.]

names = {'HotMaps_cat_1': '3D clusters tumor',
         'HotMaps_cat_2': '3D clusters pan-cancer',
         'CLUSTL_cat_1': 'Linear clusters tumor',
         'CLUSTL_cat_2': 'Linear clusters pan-cancer',
         'smRegions_cat_1': 'Domain enrichment tumor',
         'smRegions_cat_2': 'Domain enrichment pan-cancer',
         'PhyloP': 'Conservation',
         'PTM': 'Post-translational modification',
         'csqn_type_missense': 'Missense',
         'csqn_type_nonsense': 'Nonsense'}


"""Utils"""


def get_PFAMs_per_transcript(PFAM_files, PFAM_info, transcript):
    df_pfam = pd.read_csv(PFAM_files, sep="\t", names=["ENSEMBL_GENE", "ENSEMBL_TRANSCRIPT", "START", "END", "DOMAIN"])
    df_names = pd.read_csv(PFAM_info, sep="\t", names=["DOMAIN", "CLAN", "CLAN_NAME", "DOMAIN_NAME", "Long Name"])

    # Get domains
    df_pfam_gene = df_pfam[(df_pfam["ENSEMBL_TRANSCRIPT"] == transcript)]
    df_pfam_gene = df_pfam_gene[["ENSEMBL_TRANSCRIPT", "START", "END", "DOMAIN"]].drop_duplicates()
    df_pfam_gene = pd.merge(df_pfam_gene, df_names[["DOMAIN", "DOMAIN_NAME"]].drop_duplicates(), how="left")
    if df_pfam_gene.shape[0] > 0:
        df_pfam_gene["POS"] = df_pfam_gene.apply(lambda row: row["START"] + ((row["END"] - row["START"]) // 2), axis=1)
        df_pfam_gene["SIZE"] = df_pfam_gene.apply(lambda row: row["END"] - row["START"] + 1, axis=1)
        df_pfam_gene["Color"] = "#808080ff"  # color Pfam domain

    return df_pfam_gene


def gaussian_smooth(x, win):
    # if win is odd integer, then the mode occupies just one position

    gauss = scipy.stats.norm(0, 1)
    gaussian_weights = list(map(lambda x: gauss.pdf(x), np.linspace(-2, 2, win)))
    gaussian_weights /= sum(gaussian_weights)
    return np.convolve(x, gaussian_weights, mode='same')


def vertical_density(x, win, y_low, y_high, ax, **kwargs):
    # x: binary array

    x = gaussian_smooth(x, win)
    X = np.arange(len(x))
    Y = np.linspace(y_low, y_high, 100)
    Z = np.array([i for i in x for j in Y])
    Z = Z.reshape(len(Y), len(X), order='F')
    ax.contourf(X, Y, Z, 100, **kwargs)


"""Fetch mutations"""


def get_position(row):
    try:
        v = int("".join(row["aachange"][1:-1]))
        return v
    except Exception:
        return -1


def load_saturation(gene, ttype, base_path, shap_corrected=True):

    path_file = os.path.join(base_path, gene+"."+ttype+".prediction.tsv.gz")
    if os.path.exists(path_file):
        df = pd.read_csv(path_file, sep="\t")
        df.drop_duplicates(inplace=True)
        df["Protein_position"] = df.apply(lambda row: get_position(row), axis=1)

        # aggregate PTMs
        df["PTM"] = df.apply(lambda r: r['Phosphorylation'], axis=1)
        df["shap_PTM"] = df.apply(lambda r: r['shap_Phosphorylation'], axis=1)

        # SHAP-corrected feature values
        if shap_corrected:
            for c in df.columns:
                if c.startswith('shap_'):
                    df[c[5:]] = df[c].values * df[c[5:]].values  # multiply each feature by its corresponding SHAP

        # summarize to codon level
        df = df[~df['aachange'].isnull()]
        if len(df) > 0:
            df["AA"] = df.apply(lambda row: row["aachange"][0], axis=1)
            df3 = df.groupby(["AA", "Protein_position", "ENSEMBL_TRANSCRIPT"],
                             as_index=False).agg(
                                {
                                    "boostDM_score": np.nanmax,
                                    "boostDM_class": np.any,
                                    "HotMaps_cat_1": np.nanmax,
                                    "HotMaps_cat_2": np.nanmax,
                                    "smRegions_cat_1": np.nanmax,
                                    "smRegions_cat_2": np.nanmax,
                                    "CLUSTL_cat_1": np.nanmax,
                                    "CLUSTL_cat_2": np.nanmax,
                                    "csqn_type_missense": np.nanmean,
                                    "csqn_type_nonsense": np.nanmean,
                                    "PhyloP": np.nanmean,
                                    "PTM": np.nanmax,
                                })
            df3["gene"] = gene
            df3["cancer_type_annotations"] = ttype
            df3.sort_values(by='Protein_position', ascending=True, inplace=True)
            df3.reset_index(inplace=True)
            return df3
    print(f"file {path_file} not found...")
    return pd.DataFrame([])


"""Plot"""


def tracked_blueprint(gene, ttype, df_codon, outpath, show=False, plotlabel=None):

    wanted_df = df_codon[(df_codon['gene'] == gene) &
                         (df_codon['cancer_type_annotations'] == ttype)]

    for transcript, gene, ttype in wanted_df[["ENSEMBL_TRANSCRIPT",
                                              "gene",
                                              "cancer_type_annotations"]].drop_duplicates().values:

        # get PFAM domains and subset the mutation data
        subset_data_pfam = get_PFAMs_per_transcript(PFAM_files, PFAM_info, transcript)
        subset_data_muts = df_codon[
            (df_codon["ENSEMBL_TRANSCRIPT"] == transcript) &
            (df_codon["cancer_type_annotations"] == ttype)].sort_values(by='Protein_position',
                                                                        ascending=True)

        # define figure layout
        fig = plt.figure(figsize=(8, 3))

        # grid layout
        gs = gridspec.GridSpec(20, 13, figure=fig)

        border = -(len(names) + 2)  # bottom limit for blueprint scatter

        ax1 = plt.subplot(gs[:border, :12])
        ax2 = plt.subplot(gs[:border, 12], sharey=ax1)
        ax3 = plt.subplot(gs[border, :12], sharex=ax1)

        axes = []
        for i, track in enumerate(names):
            axes.append(plt.subplot(gs[border+i+1, :12], sharex=ax1))

        plot_codon_bands(subset_data_pfam, subset_data_muts, ax1, ax2, ax3)

        for i, track in enumerate(names):
            color = colors[i]
            alpha = alphas[i]
            axes[i].plot(subset_data_muts[track].values, color=color, alpha=alpha, lw=0.5)
            axes[i].spines['bottom'].set_visible(False)
            axes[i].spines['left'].set_linewidth(1)
            axes[i].spines['right'].set_visible(False)
            axes[i].spines['top'].set_visible(False)
            axes[i].set_yticks([])
            axes[i].set_ylabel(names[track], rotation=0, labelpad=45, fontsize=6, color=color)

        ax1.set_title(f'{gene} ({ttype})', fontsize=10)

        outfile_name = '{}/{}.{}.'.format(outpath, gene, ttype)
        plt.savefig(outfile_name+'svg', bbox_inches='tight')
        plt.savefig(outfile_name+'low_res.png', bbox_inches='tight', dpi=100)
        plt.savefig(outfile_name+'high_res.png', bbox_inches='tight', dpi=300)
        plt.savefig(outfile_name+'.png', bbox_inches='tight')
        if show:
            plt.show()
        plt.close(fig)


def plot_codon_bands(df_pfam_gene, df, ax_0, ax_1, ax_2):

    ax_0.set_ylabel("boostDM score", fontsize=8)
    ax_0.set_xticks(np.linspace(0, 1, 3))
    ax_0.spines['bottom'].set_visible(False)
    ax_0.spines['left'].set_linewidth(1)
    ax_0.spines['right'].set_visible(False)
    ax_0.spines['top'].set_visible(False)

    # set equivalent coordinates for the three possible mutations
    prot_pos = list(df.index.values + 1)
    d = df["boostDM_score"].values

    passenger_x, passenger_y, passenger_color = [], [], []
    driver_x, driver_y, driver_color = [], [], []
    dic_scores_pos = defaultdict(list)

    # for each of the positions
    for i, p in enumerate(prot_pos):

        dic_scores_pos[p].append(d[i])

        if d[i] < 0.5:
            passenger_x.append(p)
            passenger_y.append(d[i])
            passenger_color.append('#636363')
        else:
            driver_x.append(p)
            driver_y.append(d[i])
            driver_color.append('#ac0f0f')

    size = 1
    to_heat = []
    for p in dic_scores_pos.keys():
        to_heat.append(np.max(dic_scores_pos[p]))
    to_heat = np.array(to_heat)

    vertical_density(to_heat, 100, 0.5, 1., ax_0, alpha=0.4, cmap='Reds')

    ax_0.set_xlim(0, len(to_heat))

    ax_0.scatter(passenger_x, passenger_y, s=size, c=passenger_color, alpha=0.3)
    ax_0.scatter(driver_x, driver_y, s=size, c=driver_color, alpha=0.3)
    ax_0.set_xticks([])

    n, bins, patches2 = ax_1.hist(passenger_y + driver_y, bins=20, orientation='horizontal', alpha=1)
    for i in range(10):
        patches2[i].set_facecolor('#636363')
    for i in range(10, 20):
        patches2[i].set_facecolor('#ac0f0f')

    ax_1.spines['right'].set_visible(False)
    ax_1.spines['top'].set_visible(False)
    ax_1.spines['bottom'].set_visible(False)
    ax_1.spines['left'].set_visible(False)

    ax_2.set_ylim(0, 1)

    for i, r in df_pfam_gene.iterrows():
        start_base = r['START']
        size_base = r['SIZE']
        rect = patches.Rectangle(xy=(start_base, 0), width=size_base, height=5, color=r["Color"], alpha=0.5, zorder=2)
        ax_2.annotate(s=r["DOMAIN_NAME"], xy=(start_base + 1, 0.3), fontsize=5)
        ax_2.add_patch(rect)

    ax_2.set_xticks([])
    ax_2.set_yticks([])

    ax_0.set_yticks([0, 0.5, 1])
    ax_0.set_yticklabels([0, 0.5, 1])

    # remove axis
    ax_1.set_axis_off()

