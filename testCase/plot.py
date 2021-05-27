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
plt.style.use([ 'science', 'high-vis'])
plt.rcParams['axes.labelsize'] = 8
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['lines.markersize'] = 3
plt.rcParams['figure.autolayout'] = True


ctr = pd.read_csv('ndodecane_900_ctr.csv', sep=',')
OF = pd.read_csv('ndodecane_900_OF/chemFoam.out', sep='\t',skiprows=1,names=['time','T','p'])
ctrOF = pd.read_csv('ndodecane_900_ctrOF/chemFoam.out', sep='\t',skiprows=1,names=['time','T','p'])

plt.plot(ctr['time']*1000, ctr['T'],label='Cantera')
plt.plot(OF['time']*1000, OF['T'],label='chemFoam')
plt.plot(ctrOF['time']*1000, ctrOF['T'],label='chemFoam+Cantera')

plt.xlabel('Time [ms]')
plt.ylabel('T [K]')
plt.legend()
plt.savefig("compare.png",dpi=500,bbox_inches='tight',pad_inches=0.01)
plt.close()


print("finished!")