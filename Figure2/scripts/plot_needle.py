import sys
sys.path.append('..')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import collections as mc
from matplotlib import gridspec
from collections import defaultdict
from scipy.stats import norm

import conf

np.random.seed(42)


def get_PFAMs_per_transcript(PFAM_files, PFAM_info, transcript):
    df_pfam = pd.read_csv(PFAM_files, sep="\t", names=["ENSEMBL_GENE", "ENSEMBL_TRANSCRIPT", "START", "END", "DOMAIN"])
    df_names = pd.read_csv(PFAM_info, sep="\t", names=["DOMAIN", "CLAN", "CLAN_NAME", "DOMAIN_NAME", "Long Name"])

    # Get domains
    df_pfam_gene = df_pfam[(df_pfam["ENSEMBL_TRANSCRIPT"] == transcript)]
    df_pfam_gene = df_pfam_gene[["ENSEMBL_TRANSCRIPT", "START", "END", "DOMAIN"]].drop_duplicates()
    df_pfam_gene = pd.merge(df_pfam_gene, df_names[["DOMAIN", "DOMAIN_NAME"]].drop_duplicates(), how="left")
    df_pfam_gene["POS"] = df_pfam_gene.apply(lambda row: row["START"] + ((row["END"] - row["START"]) // 2), axis=1)
    df_pfam_gene["SIZE"] = df_pfam_gene.apply(lambda row: row["END"] - row["START"] + 1, axis=1)
    df_pfam_gene["Color"] = "#998ec3"

    return df_pfam_gene


def get_positions_in_CDS(transcript, path_coord):
    df = pd.read_csv(path_coord, sep='\t', low_memory=False,
                     names=['gene', 'gene_symbol', 'prot', 'chr', 's', 'e', 'aa', 'cds', 'genpos',
                            'strand', 'transcript'])

    toappend = []
    strand = ''
    for i, row in df[df['transcript'] == transcript].sort_values(by='s').iterrows():
        toappend.extend([i for i in range(row['s'], row['e'] + 1)])
        strand = row['strand']
    if strand == -1:
        toappend = toappend[::-1]

    return toappend


def plot_score_distribution(values, ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([])

    n, bins, patches = ax.hist(values, bins=20, orientation='vertical')
    for i in range(10):
        patches[i].set_facecolor('#636363')
    for i in range(10, 20):
        patches[i].set_facecolor('#ac0f0f')
    ax.set_xlim(0, 1)


def plot_barplot_drivers(count_drivers, count_total, ax):

    ax.barh(0, (count_total - count_drivers) / count_total, height=0.3, color='grey', left=0, alpha=0.5)
    ax.barh(0, count_drivers / count_total, height=0.3, left=(count_total - count_drivers) / count_total, color='darkred', alpha=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(0, 1)
    ax.yaxis.set_label_position("right")
    ax.set_ylabel('unique', rotation=0, labelpad=18)


def plot_barplot_drivers_nonunique(count_driver, count_driver_total, ax):

    ax.barh(0, (count_driver_total - count_driver) / count_driver_total, height=0.3, color='grey', left=0)
    ax.barh(0, count_driver / count_driver_total, height=0.3,
            left=(count_driver_total-count_driver)/count_driver_total, color='darkred')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlim(0, 1)
    ax.yaxis.set_label_position("right")
    ax.set_ylabel('all observed', rotation=0, labelpad=28)


def plot_gene_full_nucleotide(subset_data_pfam, df, transcript, path_coord, ax0, ax1, all_possible=False):
    # remove those mutations not falling in CDS:
    df = df[df['AA'] != 'n']

    # Configure the axis
    ax0.set_title('Observed Mutations')
    ax0.set_ylabel("mutation count")

    ax0.spines['bottom'].set_visible(False)
    ax0.spines['left'].set_linewidth(1)
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    ax0.tick_params(axis='y', labelsize=6, pad=0.25, width=0.25, length=1.5)
    ax1.set_yticks([])

    # set equivalent coordinates for the three possible mutations
    set_coordinates = get_positions_in_CDS(transcript, path_coord)

    # we need to get the set of equivalent coordinates per gene
    equivalent_coordinates = {coord: i for i, coord in enumerate(set_coordinates)}
    vals_coord = list(equivalent_coordinates.values())
    axs = [ax0]
    for ax in axs:
        ax.set_xlim(np.min(vals_coord), np.max(vals_coord))
        ax.set_xticks([])

    # plot observed mutations
    pos_list = df["pos"].tolist()
    ys = df["number_observed_muts"].values

    d = df["boostDM_score"].values

    coordinates_mutations = []

    passenger_x = []
    passenger_y = []
    passenger_color = []

    driver_x = []
    driver_y = []
    driver_color = []

    # for each of the positions
    for i, p in enumerate(pos_list):
        if ys[i] > 0:

            coordinates_mutations.append([(equivalent_coordinates[p], 0), (equivalent_coordinates[p], ys[i] - 0.1)])

            if d[i] < 0.5:

                passenger_x.append(equivalent_coordinates[p])
                if all_possible:
                    passenger_y.append(d[i])
                else:
                    passenger_y.append(ys[i])
                passenger_color.append('#636363')

            else:
                driver_x.append(equivalent_coordinates[p])
                if all_possible:
                    driver_y.append(d[i])
                else:
                    driver_y.append(ys[i])
                driver_color.append('#ac0f0f')

    lc = mc.LineCollection(coordinates_mutations, colors='black', linewidths=1, alpha=0.3)
    ax0.add_collection(lc)

    size = 12
    ax0.scatter(passenger_x, passenger_y, s=size, c=passenger_color, alpha=0.7, label='passenger')
    ax0.scatter(driver_x, driver_y, s=size, c=driver_color, alpha=0.7, label='driver')

    leg = ax0.legend()
    leg.get_frame().set_linewidth(0.0)

    ax1.set_ylim(0, 1)

    for i, r in subset_data_pfam.iterrows():
        start_base = 3 * r['START']
        size_base = 3 * r['SIZE']
        rect = patches.Rectangle(xy=(start_base, 0), width=size_base, height=5, color=r["Color"], alpha=0.5, zorder=2)
        ax1.annotate(s=r["DOMAIN_NAME"], xy=(start_base + 1, 0.3), fontsize=7)
        ax1.add_patch(rect)


def plot_full_nucleotide_bands(df_pfam_gene, df, ax0, ax2, ax3):

    df.sort_values(by='Protein_position', ascending=True, inplace = True)
    ax0.set_ylabel("boostDM score")

    ax0.spines['bottom'].set_visible(False)
    ax0.spines['left'].set_linewidth(1)
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    ax0.tick_params(axis='y', labelsize=6, pad=0.25, width=0.25, length=1.5)

    # set equivalent coordinates for the three possible mutations
    set_coordinates = df['pos'].drop_duplicates().tolist()
    equivalent_coordinates = {coord: i for i, coord in enumerate(set_coordinates)}

    # plot observed mutations
    pos_list = df["pos"].tolist()
    d = df["boostDM_score"].values

    coordinates_mutations = []

    passenger_x = []
    passenger_y = []
    passenger_color = []

    driver_x = []
    driver_y = []
    driver_color = []

    dic_scores_pos = defaultdict(list)

    # for each of the positions
    for i, p in enumerate(pos_list):

        dic_scores_pos[equivalent_coordinates[p]].append(d[i])

        coordinates_mutations.append([(equivalent_coordinates[p], 0), (equivalent_coordinates[p], d[i] - 0.1)])

        if d[i] < 0.5:

            passenger_x.append(equivalent_coordinates[p])
            passenger_y.append(d[i])

            passenger_color.append('#636363')

        else:
            driver_x.append(equivalent_coordinates[p])
            driver_y.append(d[i])
            driver_color.append('#ac0f0f')

    size = 1
    to_heat = []
    for p, k in enumerate(dic_scores_pos.keys()):
        to_heat.append(np.max(dic_scores_pos[p]))
    to_heat = np.array(to_heat)

    vertical_density(to_heat, 100, 0.5, 1., ax0, alpha=0.4, cmap='Reds')

    ax0.set_xlim(0, len(to_heat))

    ax0.scatter(passenger_x, passenger_y, s=size, c=passenger_color, alpha=0.3)
    ax0.scatter(driver_x, driver_y, s=size, c=driver_color, alpha=0.3)
    ax0.set_xticks([])

    n, bins, patches2 = ax2.hist(passenger_y + driver_y, bins=20, orientation='horizontal', alpha=1)
    for i in range(10):
        patches2[i].set_facecolor('#636363')
    for i in range(10, 20):
        patches2[i].set_facecolor('#ac0f0f')

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    ax3.set_ylim(0, 1)

    for i, r in df_pfam_gene.iterrows():
        start_base = 3 * r['START']
        size_base = 3 * r['SIZE']
        rect = patches.Rectangle(xy=(start_base, 0), width=size_base, height=5, color=r["Color"], alpha=0.5, zorder=2)
        ax3.annotate(s=r["DOMAIN_NAME"], xy=(start_base + 1, 0.3), fontsize=7)
        ax3.add_patch(rect)

    ax3.set_xticks([])
    ax3.set_yticks([])

    ax3.set_xlabel('CDS base position')
    ax0.set_yticks([0, 0.5, 1])
    ax0.set_yticklabels([0, 0.5, 1])

    # remove axis
    ax2.set_axis_off()


def plot_full_nucleotide_bands_observed(df_pfam_gene, df, ax0, ax2, ax3):

    df.sort_values(by='Protein_position', ascending = True, inplace = True)
    ax0.set_ylabel("boostDM score")

    ax0.spines['bottom'].set_visible(False)
    ax0.spines['left'].set_linewidth(1)
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    ax0.tick_params(axis='y', labelsize=6, pad=0.25, width=0.25, length=1.5)

    # set equivalent coordinates for the three possible mutations
    set_coordinates = df['pos'].drop_duplicates().tolist()
    equivalent_coordinates = {coord: i for i, coord in enumerate(set_coordinates)}

    # plot observed mutations
    pos_list = df["pos"].tolist()
    d = df["boostDM_score"].values
    observed_type = df["observed"].values

    coordinates_mutations = []

    passenger_x = []
    passenger_y = []
    passenger_color = []
    marker_type_d,marker_type_p = [], []

    driver_x = []
    driver_y = []
    driver_color = []

    dic_scores_pos = defaultdict(list)

    # for each of the positions
    for i, p in enumerate(pos_list):

        dic_scores_pos[equivalent_coordinates[p]].append(d[i])

        coordinates_mutations.append([(equivalent_coordinates[p], 0), (equivalent_coordinates[p], d[i] - 0.1)])

        if d[i] < 0.5:

            passenger_x.append(equivalent_coordinates[p])
            passenger_y.append(d[i])
            passenger_color.append('#636363')
            if observed_type[i]:
                marker_type_p.append("o")
            else:
                marker_type_p.append("X")

        else:
            driver_x.append(equivalent_coordinates[p])
            driver_y.append(d[i])
            driver_color.append('#ac0f0f')
            if observed_type[i]:
                marker_type_d.append("o")
            else:
                marker_type_d.append("X")

    size = 3
    to_heat = []
    for p, k in enumerate(dic_scores_pos.keys()):
        to_heat.append(np.max(dic_scores_pos[p]))
    to_heat = np.array(to_heat)

    vertical_density(to_heat, 100, 0.5, 1., ax0, alpha=0.4, cmap='Reds')

    ax0.set_xlim(0, len(to_heat))
    for i in range(0,len(marker_type_p)):
        ax0.scatter(passenger_x[i], passenger_y[i], s=size, c=passenger_color[i], alpha=0.3, marker=marker_type_p[i],lw=0)
    
    for i in range(0, len(marker_type_d)):
        ax0.scatter(driver_x[i], driver_y[i], s=size, c=driver_color[i], alpha=0.3, marker=marker_type_d[i],lw=0)
    ax0.set_xticks([])

    n, bins, patches2 = ax2.hist(passenger_y + driver_y, bins=20, orientation='horizontal', alpha=1)
    for i in range(10):
        patches2[i].set_facecolor('#636363')
    for i in range(10, 20):
        patches2[i].set_facecolor('#ac0f0f')

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    ax3.set_ylim(0, 1)

    for i, r in df_pfam_gene.iterrows():
        start_base = 3 * r['START']
        size_base = 3 * r['SIZE']
        rect = patches.Rectangle(xy=(start_base, 0), width=size_base, height=5, color=r["Color"], alpha=0.5, zorder=2)
        ax3.annotate(s=r["DOMAIN_NAME"], xy=(start_base + 1, 0.3), fontsize=7)
        ax3.add_patch(rect)

    ax3.set_xticks([])
    ax3.set_yticks([])

    ax3.set_xlabel('CDS base position')
    ax0.set_yticks([0, 0.5, 1])
    ax0.set_yticklabels([0, 0.5, 1])

    # remove axis
    ax2.set_axis_off()


def gaussian_smooth(x, win):
    # if win is odd integer, then the mode occupies just one position

    gauss = norm(0, 1)
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


def plot_barplot_drivers_all_possible(df, ax):
    count_drivers = len(df[df['boostDM_class'] == True])
    count_total = len(df)
    ax.bar(0, 0, bottom=0, color='white')
    ax.bar(0, count_drivers / count_total, width=0.3, color='darkred', alpha=1)
    ax.bar(0, (count_total - count_drivers) / count_total,  bottom=count_drivers / count_total, width=0.3, color='grey', alpha=1)

    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_sketch_params(length=1)
    ax.set_xticks([])
    ax.set_ylim(0, 1)
    ax.yaxis.set_label_position("right")
    ax.yaxis.set_ticks_position("right")
    ax.set_ylabel('Potential drivers ratio')


def plot_observed_distribution(gene, ttype, data, count_driver, count_driver_unique, count_total, count_total_unique,
                               score_values, plotname=None):

    conf.config_params(font_size=6)

    wanted_df = data[(data['gene'] == gene) & (data['cancer_type_annotations'] == ttype)]

    for transcript, gene, ttype in wanted_df[["ENSEMBL_TRANSCRIPT", "gene", "cancer_type_annotations"]].drop_duplicates().values:

        # get PFAM domains and subset the mutation data
        subset_data_pfam = get_PFAMs_per_transcript(conf.PFAM_files, conf.PFAM_info, transcript)
        subset_data_muts = data[
            (data["ENSEMBL_TRANSCRIPT"] == transcript) & (data["cancer_type_annotations"] == ttype)]

        # define figure layout
        fig = plt.figure(figsize=(8, 2.25))
        plt.suptitle('{}_{}'.format(gene, ttype))
        gs = gridspec.GridSpec(11, 3, figure=fig)
        ax1 = plt.subplot(gs[1:-1, :2])
        ax2 = plt.subplot(gs[-1, :2], sharex=ax1)
        ax3 = plt.subplot(gs[0:4, 2])
        ax4 = plt.subplot(gs[7:9, 2])
        ax5 = plt.subplot(gs[9:, 2])

        # plot for each axes
        plot_score_distribution(score_values, ax3)
        plot_barplot_drivers(count_driver_unique, count_total_unique, ax5)
        plot_barplot_drivers_nonunique(count_driver, count_total, ax4)
        plot_gene_full_nucleotide(subset_data_pfam, subset_data_muts, transcript, conf.path_coord, ax1, ax2)

        ax3.set_xlabel('boostDM score')
        ax5.spines['bottom'].set_sketch_params(length=1)
        ax5.set_xlabel('passenger-driver ratio')
        ax2.set_xlabel('CDS base position')

        if plotname is not None:
            plt.savefig(plotname+'.svg', dpi=300, bbox_inches='tight')
            plt.savefig(plotname+'.png', dpi=300, bbox_inches='tight')
        plt.show()

