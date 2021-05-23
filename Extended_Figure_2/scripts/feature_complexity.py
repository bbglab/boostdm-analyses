import os
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns


def shap_table(gene, ttype, prediction_folder):
    
    fn = os.path.join(prediction_folder, f'{gene}.{ttype}.prediction.80.30.tsv.gz')
    try:
        saturation_df = pd.read_csv(fn, sep='\t')
    except FileNotFoundError as e:
        return
    return saturation_df[saturation_df['boostDM_class']][[a for a in saturation_df.columns if a.startswith('shap')]]


def shap_table_geneset(geneset, prediction_folder):

    total_df = []
    for gene, ttype in geneset:
        shaps = shap_table(gene, ttype, prediction_folder)
        total_df.append(shaps)
    if len(total_df) > 0:
        table = pd.concat(total_df, axis=0)
        return table
    return


def shap_heatmap(geneset):
    
    table = shap_table_geneset(geneset)
    if table is not None:
        sns.heatmap(table, yticklabels=False, cbar_kws={'label': 'SHAP'})
        plt.show()


def shap_pca_plot(gene, ttype):
    
    table = shap_table_geneset([(gene, ttype)])
    if table is not None:
        X = table.values
    else:
        return
    
    # scaling
    scaler = StandardScaler()
    scaler.fit(X)
    X = scaler.transform(X)
    
    # low-rank representation
    model = PCA()
    model.fit_transform(X)
    y = [0.] + list(model.explained_variance_ratio_.cumsum())
    return y


def linear_complexity(gene, ttype, prediction_folder):
    """Close to zero (resp. one) means low (resp. high) complexity"""

    try:
        table = shap_table_geneset([(gene, ttype)], prediction_folder)
        if table is not None:
            X = table.values

        # scaling
        scaler = StandardScaler()
        scaler.fit(X)
        X = scaler.transform(X)

        # low-rank representation
        model = PCA()
        model.fit_transform(X)
        auc = model.explained_variance_ratio_.cumsum().mean()
        score = 1 - 2 * (auc - 0.5)
        return score
    except (FileNotFoundError, ValueError) as e:
        return None
