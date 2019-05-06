# Author:  DINDIN Meryll
# Date:    05 May 2019
# Project: PopulationGenomics

import time
import tqdm
import joblib
import pandas as pd

from sklearn.decomposition import IncrementalPCA

individuals = joblib.load('data/chr21_patients.jb')
incremental = IncrementalPCA(n_components=3, copy=False)

for index in tqdm.tqdm(range(len(individuals)//8)):
    chr21 = pd.read_parquet('data/chr21.pq', columns=individuals[8*index:8*(index+1)]).values
    chr21[chr21 == '0|0'] = 0
    chr21[chr21 != 0] = 1
    chr21 = chr21.astype('int8').transpose()
    incremental.partial_fit(chr21)
    del chr21
    time.sleep(1)

joblib.dump(incremental, 'chr21_pca.jb')