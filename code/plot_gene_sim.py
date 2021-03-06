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
sns.set(context="poster", style = "ticks")
import matplotlib.pyplot as plt
from sinfo import sinfo
from random import randrange

# change fates to categorical values
def fate_conv(fate_arr):
    fate_arr[fate_arr == 0] = "Naive"
    fate_arr[fate_arr == 1] = "Precursor"
    fate_arr[fate_arr == 2] = "Th1"
    fate_arr[fate_arr == 3] = "Tfh"
    fate_arr[fate_arr == 4] = "Tr1"

filenames = ["scseq_fixed_probabilities.h5"]
sc.settings.figdir = '../figures/'
th1_genes = ["Th1_"+str(i) for i in range(20)]
tfh_genes = ["Tfh_"+str(i) for i in range(20)]
tr1_genes = ["Tr1_"+str(i) for i in range(20)]

colnames = th1_genes + tfh_genes + tr1_genes

savenames = ["fixed_probs"]

for file,name in zip(filenames, savenames):
    data = h5py.File('../output/'+file, 'r')
    
    
    df_list = [pd.DataFrame(data[key]) for key in data.keys()]
    data.close()
    
    # get single dataframe out of all simulations and transpose
    df = df_list[0]
    df = df.T
    
    # split data frame into gene counts and corresponding cell fates
    df_genes = df.iloc[:,:-1]
    df_fates = df.iloc[:,-1]
    
       
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
    g.savefig("../figures/heatmap_" + name +".png")

    # plot pca
    sc.tl.pca(adata)
    sc.pl.pca_variance_ratio(adata, log=True)

    sc.pl.pca(adata, color = ["Th1_1", "Tfh_1", "Tr1_1"], save = "_"+name+".svg",
              color_map = "Blues")
    
    # cluster cells
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=0.1)
    
    sc.pl.pca(adata, color = ["leiden", "fate"])

    sc.pp.neighbors(adata)
    sc.tl.umap(adata, min_dist = 1.5, spread = 1.0)
    

    sc.pl.umap(adata, color = ["Th1_1", "Tfh_1", "Tr1_1"], save = "_"+name+".png",
               title = ["IFNg", "IL21", "IL10"],
               color_map = "Blues")

sinfo(write_req_file = True)