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

filenames = ["scseq_environment_fate.h5", "scseq_predetermined_fate.h5"]

data = h5py.File('../output/'+filenames[1], 'r')
sc.settings.figdir = '../figures/'

th1_genes = ["th1_"+str(i) for i in range(20)]
tfh_genes = ["tfh_"+str(i) for i in range(20)]
tr1_genes = ["tr1_"+str(i) for i in range(20)]

colnames = th1_genes + tfh_genes + tr1_genes

df_list = [pd.DataFrame(data[key]) for key in data.keys()]
data.close()

# get single dataframe out of all simulations and transpose
df = df_list[0]
df = df.T

# split data frame into gene counts and corresponding cell fates
df_genes = df.iloc[:,:-1]
df_fates = df.iloc[:,-1]

# change fates to categorical values
def fate_conv(fate_arr):
    fate_arr[fate_arr == 0] = "Naive"
    fate_arr[fate_arr == 1] = "Precursor"
    fate_arr[fate_arr == 2] = "Th1"
    fate_arr[fate_arr == 3] = "Tfh"
    fate_arr[fate_arr == 4] = "Tr1"

fate_conv(df_fates)

# set elements smaller than zero to 0 but I have to check why they arose in the first place
df_genes[df_genes<0] = 0
adata = sc.AnnData(df_genes)

# add var df to anndata object
var_df = pd.DataFrame({"gene_ids" : range(60)})
adata.var = var_df
adata.var_names = colnames
# create UMIs
obs = ["UMI"+str(np.random.randint(1000,9999)) for i in range(len(df.index))]
adata.obs_names = obs


# add original fate to adata obs
adata.obs["fate"] = df_fates.values


# make a clustered heatmap
g = sns.clustermap(adata.X)
#g.savefig("../figures/heatmap_stoc_sim.png")


# plot pca
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, log=True)


sc.pl.pca(adata, color = ["th1_1", "tfh_1", "tr1_1"], save = "_stoc_sim.png")

# cluster cells
sc.pp.neighbors(adata)
sc.tl.leiden(adata, resolution=0.1)

sc.pl.pca(adata, color = ["leiden", "fate"], save = "_stoc_sim2.png")


sc.pp.neighbors(adata)
sc.tl.umap(adata, min_dist = 1.5, spread = 1.0)

fig, ax = plt.subplots()
sc.pl.umap(adata, color = ["fate"], ax = ax)
#fig.savefig("../figures/umap_stoc_sim.png")