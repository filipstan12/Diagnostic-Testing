import io
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, lfilter, freqz
import os.path
import os
import argparse

dfnano = pd.read_csv("./Tables_Python_Processing/240606_Pproc_INMI.csv", encoding='utf-8', delimiter=';')
dfmeas = pd.read_csv('./Tables_Python_Processing/SampleTable.csv', delimiter=';')

print(dfmeas.head())


dfnano['CT'] = np.nan
dfnano['Observation'] = np.nan

#print(dfnano.head())

count = 0
for i, id in enumerate(dfnano['SampleID']):
    for j, idct in enumerate(dfmeas['SampleID']):
        if id == idct: #and dfnano.loc[i, 'NanoID'] == dfmeas.loc[j, 'NanoID']:
            print(id, idct)
            #print(dfct.loc[j, 'ID'], dfnano.loc[i, 'sample-ID'])
            dfnano.loc[i, 'CT'] = dfmeas.loc[j, 'CT']
            dfnano.loc[i, 'Observation'] = dfmeas.loc[j, 'Observation']

print(dfnano.head())

dfnano.to_csv('Python_combined.csv', encoding='utf-8')
print('Files stored')