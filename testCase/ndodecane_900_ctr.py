'''
# T-K, P-atm
'''
import matplotlib
matplotlib.use('Agg')
#from __future__ import division
#from __future__ import print_function
import pandas as pd
import numpy as np
import time
import cantera as ct
import matplotlib.pyplot as plt
import sys
plt.style.use(['dark_background', 'science', 'high-vis'])
plt.rcParams['axes.labelsize'] = 8
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['lines.markersize'] = 3
plt.rcParams['figure.autolayout'] = True


import sys

import cantera as ct

T = 900
P = 60*ct.one_atm # atm
fuel = "C12H26"
oxidizer = "O2:0.21,N2:0.79"
phi = 1.0

gas = ct.Solution('C12H26_Yao_s54r269.xml')
gas.TP = T, P
gas.set_equivalence_ratio(phi=1.0, fuel=fuel, oxidizer=oxidizer)
r = ct.IdealGasConstPressureReactor(gas)

sim = ct.ReactorNet([r])
sim.verbose = True

# limit advance when temperature difference is exceeded
delta_T_max = 20.
r.set_advance_limit('temperature', delta_T_max)

dt = 1.e-7
t_end = 3e-4

time = []
temp = []
while sim.time < t_end:
    sim.advance(sim.time + dt)
    time.append(sim.time)
    temp.append(r.T)




log = pd.DataFrame(data={'time': time})
log['T'] = temp
log.to_csv("ndodecane_900_ctr.csv", index=False)
