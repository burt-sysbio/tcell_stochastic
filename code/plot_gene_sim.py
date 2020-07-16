# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 17:50:53 2020

@author: Philipp
"""


import numpy as np
import h5py
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from random import randrange
data = h5py.File('../output/scseq_sim.h5', 'r')


th1_genes = ["th1_"+str(i) for i in range(20)]
tfh_genes = ["tfh_"+str(i) for i in range(20)]
tr1_genes = ["tr1_"+str(i) for i in range(20)]

colnames = th1_genes + tfh_genes + tr1_genes

df_list = [pd.DataFrame(data[key]) for key in data.keys()]
data.close()


df = df_list[0]
# set elements smaller than zero to 0 but I have to check why they arose in the first place
df[df<0] = 0
adata = sc.AnnData(df.T)

# add var df to anndata object
var_df = pd.DataFrame({"gene_ids" : range(60)})
adata.var = var_df
adata.var_names = colnames
# create UMIs
obs = ["UMI"+str(np.random.randint(1000,9999)) for i in range(1898)]
adata.obs_names = obs

sc.pl.highest_expr_genes(adata, n_top=20, )

adata_th1 = adata[:,:20]
adata_tfh = adata[:,20:40]
adata_tr1 = adata[:,40:60]

cpc_th1 = adata_th1.X.sum(axis = 1)
cpc_tfh = adata_tfh.X.sum(axis = 1)
cpc_tr1 = adata_tr1.X.sum(axis = 1)

df_cpc = pd.DataFrame({"th1_genes" : cpc_th1, 
                       "tfh_genes" : cpc_tfh, 
                       "tr1_genes" : cpc_tr1})

df_cpc = pd.melt(df_cpc)

sns.catplot(data = df_cpc, x = "variable", y = "value", kind = "violin")
sns.clustermap(adata.X)


sc.tl.pca(adata)

sc.pl.pca_variance_ratio(adata, log=True)
sc.pl.pca(adata, color = ["th1_1", "tfh_1", "tr1_1"])


sc.pp.neighbors(adata)
sc.tl.louvain(adata, resolution=0.1)