#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2021-11-19'
__version__ = '0.0.1'

import argparse
import datetime as dt
from glob import glob
from sklearn.decomposition import PCA
import pandas as pd

def main():
    
    norm_method ='not_normalised'
    Data = pd.read_csv(f'Data3/{norm_method}/Data_Extracted_Peaks/Extracted_Peaks.csv',index_col=0)
    
    pca = PCA()
    components = pca.fit_transform(Data.T)
    d3 = pd.DataFrame(components)
    d31=pd.DataFrame()
    d31['PC1']=d3.loc[1]
    d31['PC2']=d3.loc[2]
    d31.plot.scatter('PC1','PC2')
    print('plotted')
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
    description="""
        This pipeline takes the normalised resampled ECIS data (produced by PCA_HUVEC.py) and calculates the slope.
        """
    )

    main()