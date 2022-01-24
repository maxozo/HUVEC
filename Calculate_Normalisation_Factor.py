

#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2021-11-19'
__version__ = '0.0.1'

import argparse
import datetime as dt
from glob import glob
import pandas as pd

def main():
    print('ļets do some PCA')
    type = 'peak' #or slope
    norm_method ='FC_norm'
    Data = pd.read_csv(f'/Users/mo11/work/HUVEC/Data3/{norm_method}/Data_Extracted_Peaks/Metrics_Calculations.csv',index_col=0)
    # Now Loop through eachof the experiments and calculate the scaling factor based on the control.
    
    # Calculate the scaling factor for each measurement based on controls.
    All_Control_Samples = Data.loc[Data.index[Data.index.str.contains('CONTROL')]]
    
    Control_Means = {}
    for exp1 in All_Control_Samples.index.str.split('_').str[0].unique():
        All_Controls = All_Control_Samples.loc[All_Control_Samples.index[All_Control_Samples.index.str.contains(exp1)]]
        Control_Means[exp1]=All_Controls.mean()
    Data_Control_Means = pd.DataFrame(Control_Means)
    
    for idx1 in Data_Control_Means.index:
        Data_Control_Means.loc[idx1]=(Data_Control_Means.loc[idx1]-Data_Control_Means.loc[idx1].min())/(Data_Control_Means.loc[idx1].max()-Data_Control_Means.loc[idx1].min())
    
    # Now apply the normalisation to all the samples
    
    
    All_Samples = Data.loc[Data.index[~Data.index.str.contains('CONTROL')]]
    for col1 in All_Samples.columns:
        All = All_Samples[col1]
        for sample1 in All.index:
            run_id = sample1.split('_')[0]
            # print(f'{run_id} {col1}')
            nor_fact = Data_Control_Means.loc[col1,run_id]
            All_Samples.loc[sample1,col1]=All_Samples.loc[sample1,col1]*nor_fact
    
    All_Samples.to_csv(f'/Users/mo11/work/HUVEC/Data3/{norm_method}/Data_Extracted_Peaks/Norm_Metrics_Calculations.csv')        
    
    print('Done')
    
    


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
    description="""
        This pipeline takes the normalised resampled ECIS data (produced by PCA_HUVEC.py) and calculates the slope.
        """
    )

    main()