# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 22:23:24 2020

@author: Philipp
plot evolution of cell numebrs over time for stochastic simulation results stored in output folder
"""
import numpy as np
import h5py
import pandas as pd
import seaborn as sns

# load data
filenames = ["cell_numbers_environment_fate.h5", "cell_numbers_predetermined_fate.h5"]
data = [h5py.File('../output/'+filename, 'r') for filename in filenames]
# hdf5 file only has one group cell data
arr = [np.array(d["cell_data"]) for d in data]

colnames = ["time", "th1", "tfh", "tr1", "n_prec"]
df = [pd.DataFrame(a.T, columns = colnames) for a in arr]

fate_reg = ["environment", "predetermined"]
for d, name in zip(df,fate_reg):
    d["fate"] = name
    

df = pd.concat(df)
df = df.melt(id_vars = ["time", "fate"])

g = sns.relplot(data = df, x = "time", y = "value", col = "fate",
               hue = "variable", ci = "sd", kind = "line")


[d.close() for d in data]