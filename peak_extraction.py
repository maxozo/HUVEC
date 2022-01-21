#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2021-11-19'
__version__ = '0.0.1'

from glob import glob
import pandas as pd
import argparse
import numpy as np
import pandas as pd
# import streamlit as st
# import plotly.express as px
from scipy.spatial.distance import euclidean
from sklearn.metrics import r2_score
from fastdtw import fastdtw
from scipy import stats
from sklearn.decomposition import PCA
import datetime as dt
from scipy.stats import linregress
from scipy.signal import savgol_filter
from numpy import trapz
from scipy.integrate import simps

def get_injury_time(d,inp):
    slope_vals ={}
    intercept_vals={}
    # d.plot()
    # if (inp =='max'):
        # for recovery we want to smooth the spectra, whereas for the injury, no, since we want to detect sudden change.
    yhat = pd.Series(savgol_filter(d, 11, 4))
    f = pd.Series(yhat)
    f.index = d.index
    d=f
    # Start at the min point and step backwards untill the 
    d_min = d[d == d.min()].index[0]
    for i in range(d_min, 0, -1):
        print(i)
        # print(f"{d[i]} {d[i-1]}")
        injury_time_start = i
        if (all(list(d[i]>d.loc[i-20:i-1]))): 
            break





    # # d.plot()
    # for i in range(0, len(d), 1):
    #     start = i
    #     end = i+2

    #     slope_calc_window = d[start:end]
    #     ind2 = slope_calc_window.index[len(slope_calc_window.index)-1]
    #     ind1 = slope_calc_window.index[0]
    #     x1 = slope_calc_window[ind1]
    #     x2 = slope_calc_window[ind2]
    #     lr = linregress(slope_calc_window, slope_calc_window.index)
    #     # print(lr.intercept)
    #     # slope_calc = stats.linregress([x2,x1],[ind1,ind2])
    #     slope = (x2-x1)/(ind1-ind2)
    #     # print(slope)
    #     slope_vals[ind1]=[slope]
    #     intercept_vals[ind1]=[lr.intercept]
    #     # slope = slope_calc.slope

    # Slopes1 = pd.DataFrame(slope_vals).T
    # Intercepts1 = pd.DataFrame(intercept_vals).T
    # # Intercepts1.plot()
    # # Slopes1.plot()
    
    # strd = Slopes1.std()
    # injury_time_start =Slopes1[Slopes1 ==Slopes1.max()[0]].dropna().index[0]-10
   
    # Q1 = Slopes1.quantile(0.25)
    # Q3 = Slopes1.quantile(0.75)
    # IQR = Q3 - Q1
    # Highest_bound = Q3 +  1.5* IQR
    # Lowest_bond = Q1 - 1.5*IQR


    # if (inp =='min'):
    #     Slopes2 = Slopes1[(Slopes1> Highest_bound ) ].dropna()
    # else:

    # #     Slopes2 = Slopes1[(Slopes1> Highest_bound ) ].dropna()
    # #     Slopes3 = Slopes1[(Slopes1< Lowest_bond ) ].dropna()
    # #     t = list(Slopes2.index)
    # #     Slopes4 = Slopes1.drop(t)
    # #     strd = Slopes4.std()
    # #     Slopes2 = Slopes1[(Slopes1> 2*strd ) ].dropna()
    # #     injury_time_start = min(Slopes2.index)
    
    #     Slopes2 = Slopes1[(Slopes1> strd ) ].dropna()
    #     injury_time_start = min(Slopes2.index)
    
    # d.loc[:injury_time_start].plot()

    return injury_time_start        

def main():
    print('ļets do some PCA')
    df = pd.DataFrame()
    Data = pd.read_csv('/Users/mo11/work/HUVEC/Data3/Log_Data/Thrombin_Data_all_data_remapped_Resistance.csv',index_col=0)
    d2 = Data[Data['freq']==4000]
    d2 = d2.drop('freq',axis=1)
    d2.reset_index(drop=True,inplace=True)
    d2 = d2.fillna(0)

    # Here loop through the spectra and calculate the slope at each point
    # d2['E520_2'].plot()
    idx_all = pd.Series(d2.columns)
    
    All_experiments = set(idx_all.str.split('_').str[0])
    All_Experiment_Data = pd.DataFrame()
    All_calculations = {}
    # All_experiments=['12e','1e','9e']
    for exp1 in All_experiments:
        
        # All the experiments performed together for this run
        exp = idx_all[idx_all.str.contains(f"^{exp1}_")]
        all_injury_times = []
        for id1 in exp:
            # id1 =idx_all[idx_all.str.contains('E1192_1')].values[0]
            # id1 =idx_all[idx_all.str.contains('E520_2')].values[0]
            # id1='12e_E590_2'
            # id1 = '12e_E588_2'
            # id1='9e_E513_1'
            # id1 ='2e_E551_1'
            print(id1)
            if('EMPTY WEL' in id1):
                continue
            
            
            d1 = d2.loc[d2[id1]!= 0,id1]
            d1=d1[10:len(d1)-100]
            min_value = d1[d1==d1.min()].index[0]
            d1.plot()
            min1 = min_value-200
            if(min1<10):
                min1=10

            d=d1.iloc[min1:min_value+200]
            d.plot()
            injury_time_start = get_injury_time(d,'min')-5
            d_index = d.index
            d_reverse = d.iloc[::-1]
            d_reverse_reindex = d_reverse.reset_index()
            d_rev = d_reverse.reset_index(drop=True)
            arbitary_injury_time_end = get_injury_time(d_rev,'max')
            injury_time_end = int(d_reverse_reindex.iloc[arbitary_injury_time_end]['index'])
            # Now reverse the data and detect the end of injuty time

            all_injury_times.append(injury_time_start)
            
            Peak_window = d1.loc[injury_time_start:injury_time_start+200]
            Peak_window_only = d1.loc[injury_time_start:injury_time_end]
            nor_Peak_window_only = Peak_window_only - Peak_window_only.iloc[0]
            nor_Peak_window_only.plot()
            # Peak_window_only.plot()
            area = abs(trapz(nor_Peak_window_only))
            #nor_Peak_window_only.to_csv('test_trapz_method.csv')
            area2 = abs(simps(nor_Peak_window_only)) #7840
            print(area2)
            min_val = nor_Peak_window_only.min()
            min_val_time = nor_Peak_window_only[nor_Peak_window_only == min_val].index[0]
            rec_val_time = nor_Peak_window_only.index[len(nor_Peak_window_only.index)-1]
            rec_val = nor_Peak_window_only[rec_val_time]
            recovery_time =rec_val_time - min_val_time
            
            slope_val = (rec_val_time-min_val_time)/(rec_val-min_val)
            All_calculations[id1] = {'slope':slope_val, 'area trapz':area, 'area simpson':area2, 'recovery_time':recovery_time}
            print("area =", area)
            # d1.plot()
            
            print('plotted')
    
        print('Done with this experiment')
        Consensous_injury_start = pd.Series(all_injury_times).mode()[0]
        normalised_Peaks = d2.loc[Consensous_injury_start:Consensous_injury_start+200,exp]
        # normalised_Peaks.plot()
        normalised_Peaks=normalised_Peaks.reset_index(drop=True)
        All_Experiment_Data = pd.concat([All_Experiment_Data,normalised_Peaks], axis=1)
        # All_Experiment_Data.plot()

        print('plotted')
    All_calculations = pd.DataFrame(All_calculations).T
    cols1 = pd.Series(All_Experiment_Data.columns)
    All_Experiment_Data.to_csv('Data3/Data_Extracted_Peaks/Extracted_Peaks.csv')
    All_calculations.to_csv('Data3/Data_Extracted_Peaks/Metrics_Calculations.csv')
    control_samples = All_Experiment_Data[cols1[cols1.str.contains('CONTROL')]]
    control_samples.plot()
    # print(pd.DataFrame(slope_vals))
    
        

    # pca = PCA()
    # components = pca.fit_transform(d2.T)
    # d3 = pd.DataFrame(components)
    
    # d31=pd.DataFrame()
    # d31['PC1']=d3.loc[0]
    # d31['PC2']=d3.loc[1]
    # d31.plot.scatter('PC1','PC2')
    
    # x = np.array(df['E520_2'].fillna(0))
    # y = np.array(d2.fillna(0))
    # from scipy.stats import linregress
    # stats.linregress([x,x],y)

    # pd.to_datetime(df['DateTime']).map(dt.datetime.toordinal)
    # stats.linregress(df['DateTime'], df['E520_2'])
    # pd.to_datetime(df.index.map(dt.datetime.toordinal))
    # pd.to_datetime(df.index,format='hr')
    # DTW approach to sync the data
    # distance, path = fastdtw(x, y, dist=euclidean)

    # df['DateTime']=df.index
    # result = []
    # for i in range(0,len(path)):
    #     result.append([df['DateTime'].iloc[path[i][0]],
    #     df['Power'].iloc[path[i][0]],
    #     df['Voltage'].iloc[path[i][1]]])
    # df_sync = pd.DataFrame(data=result,columns=['DateTime','Power','Voltage']).dropna()
    # df_sync.plot()
    # df_sync = df_sync.drop_duplicates(subset=['DateTime'])
    # df_sync = df_sync.sort_values(by='DateTime')
    # df_sync.index = df_sync['DateTime']
    # df_sync.plot()
    # df_sync.to_csv('synchronized_dataset.csv',index=False)


    # Calculating the sliding window slope:


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
    description="""
        This pipeline takes the normalised resampled ECIS data (produced by PCA_HUVEC.py) and calculates the slope.
        """
    )

    main()
