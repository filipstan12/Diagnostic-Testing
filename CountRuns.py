import os

import pandas as pd

filelist = []
counter = 1

for dirpath, dirnames, filenames in os.walk("P:\SGP-20885\data\WP3\Samples\Validation INMI\.gnainfo", topdown=True):
    for f in filenames:
        if f.endswith('.gnainfo'):
            filelist.append(f)
            print(f)
            counter = counter + 1

print("Number of files: ", counter)
resultdat = []

for i, file in enumerate(filelist):
    splitfilename = file.split("_")
    results = splitfilename[-1]
    resultdat.append(results)

df = pd.DataFrame(resultdat, columns=['Name'])
print(df['Name'].value_counts())