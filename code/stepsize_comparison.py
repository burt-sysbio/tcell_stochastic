# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 15:36:14 2020

@author: Philipp
"""

import numpy as np
import seaborn as sns
import pandas as pd

df_1 = pd.read_csv("../output/step0001.csv")
df_2 = pd.read_csv("../output/step001.csv")
df_3 = pd.read_csv("../output/step01.csv")

df_list = [df_1, df_2, df_3]

df_list_tidy = [pd.melt(df, id_vars = ["time"]) for df in df_list]

ss = [0.0001, 0.001, 0.01]
for df, stepsize in zip(df_list_tidy, ss):
    df["stepsize"] = stepsize
    
df_all = pd.concat(df_list_tidy)
g = sns.relplot(data = df_all, x = "time", y = "value", hue = "variable", 
                kind = "line", col = "stepsize", ci = "sd")

g.savefig("../figures/stoc_model_stepsize_comparison.pdf")