import scanpy as sc
import os
import anndata
import pandas as pd
import matplotlib.pyplot as plt
from scanpy.metrics.specificity.plot import marker_genes_distribution, one_v_max_genelist
from scanpy import metrics
from pathlib import Path
import pickle
import seaborn as sns

ROOTDIR = r''
DATA_PATH = r'data'

AW_markers = pd.read_csv(DATA_PATH + r'\Meyer_5loc_AW_selection_100pca_plus_markers.csv')
AW_markers = AW_markers[AW_markers['selection_ranking'] <= 2]
to_split = [x.split(',') for x in AW_markers['marker'].to_numpy()]
Markers_Meyer = []
for cts in to_split:
    Markers_Meyer += cts
Markers_Meyer = dict.fromkeys(Markers_Meyer,[])

gene_to_celltype = dict(zip(AW_markers['index'],AW_markers['marker']))
gene_to_celltype = {k: v.split(',') for k, v in gene_to_celltype.items()}

for key, values in gene_to_celltype.items() :
    for ct in values:
        Markers_Meyer[ct] = Markers_Meyer[ct] + [str(key)]
###################################

def save_obj(obj, name ):
    with open(Path(ROOTDIR).joinpath('data').joinpath(f'{name}.pkl'), 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open(Path(ROOTDIR).joinpath('data').joinpath(f'{name}.pkl'), 'rb') as f:
        return pickle.load(f)

Markers_Meyer = load_obj('Markers_Meyer')
Markers_Barbry = load_obj('Marker_Genes_HCA')
Markers_Meyer_CT_Deprez_df = pd.read_csv(r'C:\Users\ipmc\Documents\Metrics_results\data\Markers_Meyer_CT_Deprez.csv',sep=';')

def df_to_dict(gene_df):
    """
    Transform a gene marker by CT df into a dict

    Parameters
    ----------
    gene_df
    CT in columns, and genes stored in lines. Not same size for each CT hence df is not appropriate format.

    Returns
    -------

    """
    gene_dict = dict.fromkeys(gene_df.columns)
    for col in gene_df.columns:
        gene_dict[col] = list(gene_df[col].dropna())
    return gene_dict

def dict_to_df(gene_dict):
    """
    Transform a gene marker by CT dict into a df with lots of Nan

    :param gene_dict:
    :return:
    """



Markers_Meyer_CT_Deprez = df_to_dict(Markers_Meyer_CT_Deprez_df)

def compare(genes1,genes2):
    return len(set(genes1)&set(genes2))/len(genes1)

def compare_typing(dict1, dict2):
    comp_matrix = pd.DataFrame(index=dict1, columns=dict2,dtype='float')
    i = 0
    for key1, val1 in dict1.items():
        j = 0
        for key2, val2 in dict2.items():
            comp_matrix.iloc[i,j] = compare(val1,val2)
            j += 1
        i += 1
    return comp_matrix

plt.subplots(1,3)
sns.set(font_scale = 0.5)
sns.heatmap(compare_typing(Markers_Meyer, Markers_Barbry), ax=plt.subplot((131)), xticklabels=True,yticklabels=True)
plt.gca().title.set_text('Meyer vs Barbry')
sns.heatmap(compare_typing(Markers_Barbry, Markers_Meyer), ax=plt.subplot(132), xticklabels=True,yticklabels=True)
plt.gca().title.set_text('Barbry vs Meyer')
sns.heatmap(compare_typing(Markers_Barbry, Markers_Meyer_CT_Deprez), ax=plt.subplot(133), xticklabels=True,yticklabels=True)
plt.gca().title.set_text('Barbry vs Meyer Barb style')
plt.show()

for k,v in Markers_Meyer.items():
    while len(v) < 4:
        v.append(None)

pd.DataFrame(Markers_Meyer).to_csv(r'C:\Users\ipmc\Documents\Metrics_results\data\Marker_Meyer')

