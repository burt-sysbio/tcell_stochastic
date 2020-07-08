# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 09:59:14 2020

@author: Philipp
"""


import numpy as np
import pandas as pd
import seaborn as sns
from scipy.integrate import odeint
df = pd.read_csv("../output/single_cell.csv")
df["model"] = "stoc"

d = {
    "alpha" : 2,
    "beta" : 2,
    "r_div" : 1.0,
    "r_death" : 1.0,
    }


def single_cell(state, time, d):
    naive = state[0]
    naive_1 = state[1]
    th1 = state[2]
    
    dt_naive = -d["beta"]*naive
    dt_naive_1 = d["beta"]*(naive-naive_1)
    dt_th1 = d["beta"]*naive_1+ d["r_div"]*th1-d["r_death"]*th1
    
    return [dt_naive, dt_naive_1, dt_th1]

t = np.arange(0,20,0.01)
y0 = np.zeros(3)
y0[0] = 1
state = odeint(single_cell, y0, t, args = (d,))
df_ode = pd.DataFrame(state, columns = ["naive", "naive", "Th1"])
df_ode["time"] = t
df_ode = df_ode[["time", "Th1"]]

n_cells =500
df_ode.Th1 = df_ode.Th1*n_cells
df_ode["model"] = "ODE"

df = pd.concat([df, df_ode])
sns.relplot(data = df, x = "time", y = "Th1", hue = "model", kind = "line", ci = "sd")