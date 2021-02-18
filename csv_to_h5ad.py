import scanpy as sc
import pandas as pd
import anndata
import csv

df = pd.read_csv('krasnow_hlca_10x_UMIs.csv', index_col=0)
metadata = pd.read_csv('krasnow_hlca_10x_metadata.csv', index_col=0)
genes = df.index

adata = anndata.AnnData(X=df, obs=metadata, var=genes)

adata.write_h5ad('travigliani.h5ad')
