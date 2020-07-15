# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 17:50:53 2020

@author: Philipp
"""


import numpy as np
import h5py
import pandas as pd
import scanpy as sc

data = h5py.File('../output/scseq_sim.h5', 'r')


th1_genes = ["th1_"+str(i) for i in range(5)]
tfh_genes = ["th1_"+str(i) for i in range(5)]
tr1_genes = ["tr1_"+str(i) for i in range(5)]

colnames = th1_genes + tfh_genes + tr1_genes

df_list = [pd.DataFrame(data[key]) for key in data.keys()]


df = df_list[0]
adata = sc.AnnData(df.T)
#adata.var = colnames
data.close()
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)