# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 22:23:24 2020

@author: Philipp
"""


import numpy as np
import h5py

data = h5py.File('../output/model_output.h5', 'r')

arr = np.array(data["cell_data"])

import pandas as pd

colnames = ["time", "th1", "tfh", "tr1", "n_prec"]
df = pd.DataFrame(arr.T, columns = colnames)

df_tidy = df.melt(id_vars = ["time"])

import seaborn as sns

g = sns.relplot(data = df_tidy, x = "time", y = "value",
               hue = "variable", ci = "sd", kind = "line")


data.close()