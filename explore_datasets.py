import scanpy as sc
import os
import anndata
import pandas as pd
import matplotlib.pyplot as plt
from scanpy.metrics.specificity.plot import marker_genes_distribution, one_v_max_genelist
from scanpy.metrics.specificity.analysis import specificity_quality_control
from scanpy import metrics
from pathlib import Path
import pickle
import seaborn as sns
from collections import Counter


ROOTDIR = r''
DATA_PATH = r'data/'

full_data_dict_path = {'lukassen': 'lukassen20_lung_orig.processed.h5ad',
                       'madisson': 'madissoon19_lung.processed.h5ad',
                       'barbry': 'HCA_Barbry_Grch38_Raw_filter_Norm.h5ad',
                       'vieira_alv': 'vieira19_Alveoli_and_parenchyma_anonymised.processed.h5ad',
                       'vieira_bronch': 'vieira19_Bronchi_anonymised.processed.h5ad',
                       'vieira_nas': 'vieira19_Nasal_anonymised.processed.h5ad'}


def load_full_data(name):
    return anndata.read_h5ad(DATA_PATH + full_data_dict_path[name])

barbry = load_full_data('barbry')
madisson = load_full_data('madisson')
lukassen = load_full_data('lukassen')
vieira_alv = load_full_data('vieira_alv')
vieira_bronch = load_full_data('vieira_bronch')
vieira_nas = load_full_data('vieira_nas')


def load_obj(name):
    with open(Path(ROOTDIR).joinpath('data').joinpath(f'{name}.pkl'), 'rb') as f:
        return pickle.load(f)

Markers_Barbry = load_obj('Marker_Genes_HCA')

datas = [barbry, madisson, lukassen, vieira_alv, vieira_bronch, vieira_nas]

summary = pd.DataFrame(columns=('ncells', 'ngenes', 'nCellTypes'), index=('barbry','madisson','lukassen','vieira_alv','vieira_bronch','vieira_nas'))

i = 0

for data in datas :
    summary.iloc[i, 0] = data.n_obs
    summary.iloc[i, 1] = data.n_vars
    summary.iloc[i, 2] = len(Counter(data.obs['CellType']).keys())
    i += 1


vieira_bronch.var[vieira_bronch.var['highly_variable']==True]
sns.scatterplot(vieira_bronch.obsm['X_umap_hm'][:,0], vieira_bronch.obsm['X_umap_hm'][:,1])

specificity_quality_control(adata=madisson,
                            marker_genes=Markers_Barbry,
                            partition_key='CellType',
                            project_dir=r'C:\Users\ipmc\Documents\Metrics_results\deprez_markers_vs_datasets\deprez_vs_madisson',
                            plot_umap=True)
