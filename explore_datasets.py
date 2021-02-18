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

trav_count = pd.read_csv(DATA_PATH + 'krasnow_hlca_10x_UMIs.csv')
trav_meta = pd.read_csv(DATA_PATH + 'krasnow_hlca_10x_metadata.csv')

trav_count
travigliani = anndata.AnnData(X = trav_count[], obs = trav_meta, var = trav_count.columns)

full_data_dict_path = {'lukassen_par': 'lukassen20_lung_orig.processed.h5ad',  # parenchyma
                       'lukassen_AW': 'lukassen20_airway_orig.processed.h5ad',  # bronchi, epithelial
                       'madissoon': 'madissoon19_lung.processed.h5ad',
                       'barbry': 'HCA_Barbry_Grch38_Raw_filter_Norm.h5ad',
                       'vieira_alv': 'vieira19_Alveoli_and_parenchyma_anonymised.processed.h5ad',
                       'vieira_bronch': 'vieira19_Bronchi_anonymised.processed.h5ad',
                       'vieira_nas': 'vieira19_Nasal_anonymised.processed.h5ad',
                       'travigliani': 'facs_normal_lung_blood_scanpy.20200205.RC4.h5ad',
                       'trav2': 'droplet_normal_lung_blood_scanpy.20200205.RC4.h5ad'}


def load_full_data(name):
    return anndata.read_h5ad(DATA_PATH + full_data_dict_path[name])

barbry = load_full_data('barbry')
madissoon = load_full_data('madissoon')
lukassen_par = load_full_data('lukassen_par')
lukassen_AW = load_full_data('lukassen_AW')
vieira_alv = load_full_data('vieira_alv')
vieira_bronch = load_full_data('vieira_bronch')
vieira_nas = load_full_data('vieira_nas')

travigliani = load_full_data('travigliani')
trav2 = load_full_data('trav2')



def load_obj(name):
    with open(Path(ROOTDIR).joinpath('data').joinpath(f'{name}.pkl'), 'rb') as f:
        return pickle.load(f)

Markers_Barbry = load_obj('Marker_Genes_HCA')

datas = [barbry, madissoon, lukassen, vieira_alv, vieira_bronch, vieira_nas]

summary = pd.DataFrame(columns=('ncells', 'ngenes', 'nCellTypes'), index=('barbry','madisson','lukassen','vieira_alv','vieira_bronch','vieira_nas'))

i = 0

for data in datas :
    summary.iloc[i, 0] = data.n_obs
    summary.iloc[i, 1] = data.n_vars
    summary.iloc[i, 2] = len(Counter(data.obs['CellType']).keys())
    i += 1


vieira_bronch.var[vieira_bronch.var['highly_variable']==True]
sns.scatterplot(vieira_bronch.obsm['X_umap_hm'][:,0], vieira_bronch.obsm['X_umap_hm'][:,1])

specificity_quality_control(adata=madissoon,
                            marker_genes=Markers_Barbry,
                            partition_key='CellType',
                            project_dir=r'C:\Users\ipmc\Documents\Metrics_results\deprez_markers_vs_datasets\deprez_vs_madisson',
                            plot_umap=True)


# Celltypes nomenclature

celltypes = []
celltypes.append(list(barbry.obs['CellType'].cat.categories))
celltypes.append(list(madissoon.obs['CellType'].cat.categories))
celltypes.append(list(lukassen_par.obs['CellType'].cat.categories))
celltypes.append(list(lukassen_AW.obs['CellType'].cat.categories))
celltypes.append(list(vieira_alv.obs['CellType'].cat.categories))
celltypes.append(list(vieira_bronch.obs['CellType'].cat.categories))
celltypes.append(list(vieira_nas.obs['CellType'].cat.categories))

to_save = pd.DataFrame(celltypes,index=['barbry','madissoon','lukassen_par','lukassen_AW','vieira_alv','vieira_bronch','vieira_nas']).transpose()
to_save.to_csv(DATA_PATH + 'CellTypes_compare.csv')

