# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 15:45:27 2018

@author: Bulat
"""

import mne
import pactools
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import signal
import pandas as pd
sample_rate = 250
#%%
filenames = os.listdir("SavedData")
(b, a) = signal.butter(3, np.array([0.5, 50]) / (sample_rate / 2), 'bandpass')
Data = [signal.lfilter(b, a, np.loadtxt("SavedData\\"+x, comments="%", delimiter=",",usecols=(0,1,2)).T,1) for x in filenames]
Data = signal.lfilter(b, a, np.loadtxt("SavedData\\"+filenames[9], comments="%", delimiter=",",usecols=(0,1,2)).T,1)
passive = signal.lfilter(b, a, np.loadtxt("SavedData\\"+filenames[7], comments="%", delimiter=",",usecols=(0,1,2)).T,1)
poem = signal.lfilter(b, a, np.loadtxt("SavedData\\"+filenames[7], comments="%", delimiter=",",usecols=(0,1,2)).T,1)
Data = np.array_split(Data[1][100:],len(Data[1])/(250*70))
passive = np.array_split(passive[1][100:],len(passive[1])/(250*70))
poem = np.array_split(poem[1][100:],len(poem[1])/(250*70))
#%%
from pactools.utils import BandPassFilter
f = BandPassFilter(fs=100., fc=50., bandwidth=1., n_cycles=None)
f.plot()
for x in Data:
    plt.plot(Data[2])
#%%
from pactools import Comodulogram, REFERENCES
fs = 250
low_fq_range = np.linspace(1,40,18)
high_fq_range = np.linspace(1,40,18)
n_surrogates = 200
methods = [
    'ozkurt', 'canolty', 'tort', 'penny', 'vanwijk', 'duprelatour', 'colgin',
    'sigl', 'bispectrum'
]
mx = []
n_surrogates = 100
for x in filenames:
    for i in range(1,5):
        if "Subject"+str(i) in x:
            estimator = Comodulogram(fs=fs, low_fq_range=low_fq_range,
                                     low_fq_width=1,high_fq_range=high_fq_range,high_fq_width =1, method='duprelatour',
                                     progress_bar=True,n_jobs=-1,n_surrogates=n_surrogates)
            Data = signal.lfilter(b, a, np.loadtxt("SavedData\\"+x, comments="%", delimiter=",",usecols=(0,1,2)).T,1)
            Data = np.array_split(Data[1][100:],len(Data[1])/(250*70))
            q = estimator.fit(Data[0])
            p_value = 0.05
            estimator.plot(contour_method='comod_max', contour_level=p_value,
                       titles=['With a p-value on the distribution of maxima'])
            mx.append(estimator.get_maximum_pac())
            estimator.save("Subject"+str(i)+x.split("_")[-2],overwrite=True)
            plt.title("Subject"+str(i)+x.split("_")[-2])
            plt.savefig("pics\\"+str(x)+'.png')
            plt.show()

#%%
def nmi(active,passive):
    nmi = np.zeros(active.shape)
    for i,row in enumerate(active):
        for j,value in enumerate(row):
            nmi[i][j] = (active[i][j] - passive[i][j])/passive[i][j]
    return nmi

def app(x):
    if len(relaxed[relaxed.subject == x.subject]['cm'])>0:
        z = nmi(x.cm.comod_,relaxed[relaxed.subject == x.subject]['cm'].iloc[0].comod_)
        nmis[int(x.subject)][x.activity] = z
    
from pactools.comodulogram import read_comodulogram
filenames = os.listdir("cm")
activity = pd.DataFrame(columns=["subject","activity","cm"])
nmis = {}
for i in range(0,4):
    nmis[i] = {}
relaxed = pd.DataFrame(columns=["subject","cm"])
for x in filenames:
    if 'relaxed' in x:
        relaxed = relaxed.append({"subject":x[7],"cm":read_comodulogram("cm\\"+x)},ignore_index = True)
    else:
        activity = activity.append({"subject":x[7],"activity":x[8:],"cm":read_comodulogram("cm\\"+x)},ignore_index = True)
actovoty = activity.transform(app,axis=1)
#%%
for key,dic in nmis.items():
    print(key)
    for k,com in dic.items():
        print(k)
        plt.imshow(com.reshape(18,18))
        plt.show()