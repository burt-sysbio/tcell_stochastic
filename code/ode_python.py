# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 09:40:30 2020

@author: Philipp
"""


import numpy as np
import seaborn as sns
import pandas as pd
from scipy.integrate import odeint

sns.set(context = "poster")

d = {
    "alpha" : 2,
    "beta" : 2,
    "deg_myc" : 0.5,
    "r_div" : 0,
    "r_div_th1" : 1,
    "EC50_myc" : 0.5,
    "r_death" : 1,
    "prob_th1" : 0.5,
    "prob_tfh" : 0.2,
    "prob_tr1" : 0.3,
    }

def pos_fb(c, EC50):
    return (c/(c+EC50))


def ode_model(state, time, d):
    naive = state[0]
    naive_1 = state[1]
    th1 = state[2]
    tfh = state[3]
    tr1 = state[4]
    
    myc = np.exp(-d["deg_myc"]*time)
    r_div_th1 = d["r_div_th1"]*pos_fb(myc, d["EC50_myc"])
    
    dt_naive = -d["beta"]*naive
    dt_naive_1 = d["beta"]*(naive-naive_1)
    dt_th1 = d["beta"]*naive_1*d["prob_th1"]+(r_div_th1-d["r_death"])*th1
    dt_tfh = d["beta"]*naive_1*d["prob_tfh"]+(r_div_th1-d["r_death"])*tfh
    dt_tr1 = d["beta"]*naive_1*d["prob_tr1"]+(r_div_th1-d["r_death"])*tr1
    
    dt_state = [dt_naive, dt_naive_1, dt_th1, dt_tfh, dt_tr1]
    return dt_state


t = np.arange(0,10,0.01)
y0 = np.zeros(5)
y0[0] = 1
state = odeint(ode_model, y0, t, args = (d,))

df = pd.DataFrame(state, columns = ["naive", "naive_1", "Th1", "Tfh", "Tr1"])
df = df[["Th1", "Tfh", "Tr1"]]
df["time"] = t

n_cells = 2000
df_tidy = df.melt(id_vars = ["time"])
df_tidy.value = df_tidy.value*n_cells
df_tidy["model"] = "ODE"

df_stocsim = pd.read_csv("../output/teststoring.csv")
df_tidy_stoc = df_stocsim.melt(id_vars = ["time"])
df_tidy_stoc["model"] = "Stoc."

df_all = pd.concat([df_tidy, df_tidy_stoc])
g = sns.relplot(data = df_all, x = "time", y = "value", hue = "variable", 
                kind = "line", style = "model", ci = "sd")

g.set(xlabel = "time pI", ylabel = "cells")
g.savefig("../figures/stochastic_ode_comparison.pdf")

