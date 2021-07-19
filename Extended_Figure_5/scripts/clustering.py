# hierarchical clustering functions

import sys
sys.path.append('..')

import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from seaborn import clustermap, color_palette
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import fcluster
import scipy.cluster.hierarchy as hierarchy

import conf


def get_moa_classes():

    df = pd.read_csv(conf.drivers_path, sep='\t')
    role_dict = dict(zip(zip(df.SYMBOL.values, df.CANCER_TYPE), df.ROLE.values))
    return role_dict


def fancy_dendrogram(*args, **kwargs):

    """
    adapted from:
    https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial
    """

    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d

    ddata = hierarchy.dendrogram(*args, **kwargs)

    return ddata


def var_explained(X, labels):

    centroid = np.mean(X, axis=0)  # sum of squares
    total_ss = sum([np.linalg.norm(centroid - X[i, :])**2 for i in range(X.shape[0])])
    n_clusters = []
    scores = []
    for i, flat_cluster in enumerate(labels):
        label_set = set(flat_cluster)
        cluster_ss = 0
        n_clusters.append(len(label_set))
        for label in label_set:
            indices = np.where(flat_cluster == label)[0]
            cluster_centroid = np.mean(np.take(X, indices, axis=0), axis=0)
            cluster_ss += sum([np.linalg.norm(cluster_centroid - X[j, :])**2 for j in indices])
        scores.append(1 - (cluster_ss / total_ss))
    return scores 


def generate_hierarchy(df, n_labels_bound=10):

    X = df.values
    Y = pdist(X, metric='euclidean')
    linkage = hierarchy.linkage(Y, method='ward')
    dist_matrix = squareform(Y)

    # plot dendrogram and display cophenetic distances

    ddgram = fancy_dendrogram(linkage,
                              truncate_mode='level',
                              labels=list(df.index.values),
                              leaf_rotation=90,
                              color_threshold=0,
                              above_threshold_color='green',
                              no_plot=True)

    # compute silhouette for a number of splits until reaching cophenetic distance ~ 0.1

    scores, classes, labels_list = [], [], []
    round_low = lambda x: x - abs(x - np.round(x, 3))
    upper_coph = round_low(ddgram['dcoord'][-1][1])

    if upper_coph > 0.1:

        for i in np.linspace(upper_coph, 0.1, num=1000):

            labels = fcluster(linkage, i, criterion='distance')
            n_labels = len(set(labels))

            if n_labels > n_labels_bound:
                break
                
            if n_labels in classes:
                continue
            else:
                classes.append(n_labels)
                labels_list.append(labels)
                scores.append(silhouette_score(dist_matrix, labels, metric='precomputed'))
    
    return X, linkage, scores, classes, labels_list


def draw_flat_cluster(X, linkage, columns, labels, output_path, mode='flat', moa_colors=None, title=None, fn=None):

    """mode: flat or moa"""

    if mode == 'flat':
        color_list = color_palette('Dark2', 20).as_hex()
        colors = [color_list[i] for i in labels]
    elif (mode == 'moa') and (moa_colors is not None):
        colors = moa_colors

    g = clustermap(X.T, col_linkage=linkage, col_colors=colors, 
                   row_cluster=False, figsize=(13, 8), yticklabels=columns,
                   vmin=-2.5, vmax=2.5, xticklabels=False, cmap="RdBu_r")
    
    g.ax_heatmap.set_ylabel('Features')
    g.ax_col_dendrogram.set_title(title)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    if fn is not None:
        plt.savefig(os.path.join(output_path, f'{fn}.clustermap.png'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(output_path, f'{fn}.clustermap.svg'), dpi=300, bbox_inches='tight')

    if title is not None:
        plt.title(title)

    plt.show()
    
    
def cluster_selection_plot(classes, silh_scores, variance_explained):
    
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(classes, silh_scores, '-o', label='mean silhouette')
    ax.plot(classes, variance_explained, '-o', label='var explained')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Number of Clusters')
    ax.set_ylabel('Score')
    plt.legend()
    plt.title('Clustering Selection')
    plt.show()


