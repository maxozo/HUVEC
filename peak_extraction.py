#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2021-11-19'
__version__ = '0.0.1'

import argparse
import datetime as dt
from glob import glob

import numpy as np
import pandas as pd
from fastdtw import fastdtw
from numpy import trapz
from scipy import stats
from scipy.integrate import simps
from scipy.signal import savgol_filter
# import streamlit as st
# import plotly.express as px
from scipy.spatial.distance import euclidean
from scipy.stats import linregress
from sklearn.decomposition import PCA
from sklearn.metrics import r2_score


type = 'peak' #or slope
treatment='Wounding' #Thrombin|Wounding|Prolif
measurement = 'Capacitance' # Resistance|Capacitance|Impedence
norm_method ='not_normalised' #FC_norm | mean_norm | min_max | min_max_controls | not_normalised


def get_injury_time(d1,inp,filter_width,window_size):
    # We only smoothen from the point of the min value, back, to avoid the random sudden peak
    d_min = d1[d1 == d1.min()].index[0]
    d2=d1[:d_min]
    try:
        f = pd.Series(savgol_filter(d2, window_size, filter_width))
    except:
        print('This must be an outlier since it has already selected a small window.')
        f =d2

    try:
        f.index = d2.index
        d=f
        injury_time_start=0
        # d1.plot()
        # Start at the min point and step backwards untill the 
        d_min = d[d == d.min()].index[0]        
        if inp=='min':
            
            
            for i in range(d_min-filter_width, 0, -1):
                injury_time_start = i
                if (all(list(d[i]>d.loc[i-20:i-1]))): 
                    break
        else:
            
            for i in range(d_min-filter_width, 0, -1):
                injury_time_start = i
                if (all(list(d[i]>d.loc[i-20:i-1]))): 
                    break
        failed=False   
    except:
        print("Can not step back in the slope since it is already in the minimum")        
        failed=True    


    return injury_time_start,failed      

def calculate_peak_metrics(Peak_window_only,injury_time_start,injury_time_end,min_value,d1):

    nor_Peak_window_only = Peak_window_only - Peak_window_only.iloc[0]
    # dnor_Peak_window_only.plot()
    # Peak_window_only.plot()
    area = abs(trapz(nor_Peak_window_only))
    #nor_Peak_window_only.to_csv('test_trapz_method.csv')
    area2 = abs(simps(nor_Peak_window_only)) #7840
    # print(area2)
    min_val = nor_Peak_window_only.min()
    min_val_time = nor_Peak_window_only[nor_Peak_window_only == min_val].index[0]
    rec_val_time = nor_Peak_window_only.index[len(nor_Peak_window_only.index)-1]
    rec_val = nor_Peak_window_only[rec_val_time]
    recovery_time =rec_val_time - min_val_time
    # slope_before_treatment
    before_treatment = d1.loc[:injury_time_start]
    # before_treatment.plot()
    slope_val_before_treatment = (injury_time_start-before_treatment.index[0])/(before_treatment[injury_time_start]-before_treatment.iloc[0])
    # slope after treatment
    after_treatment = d1.loc[injury_time_end:injury_time_end+200]
    # after_treatment.plot()
    slope_val_after_treatment = (200)/(after_treatment[after_treatment.index[-1]]-after_treatment.iloc[0])
    
    # time needed to recover to 50% injury area
    injury_val_at_50 = before_treatment[before_treatment.index[-1]]-Peak_window_only.iloc[0]+(nor_Peak_window_only.min())/2 #changed - )
    values_after_time_of_50 = nor_Peak_window_only.loc[min_value:][nor_Peak_window_only.loc[min_value:]>injury_val_at_50]
    # values_after_time_of_50.plot()
    try:
        time_to_recover_to_50 = values_after_time_of_50.index[0]-min_value
    except:
        time_to_recover_to_50 =None
    # slope of recovery
    slope_val = (rec_val_time-min_val_time)/(rec_val-min_val)
    return {'slope_of_recovery':slope_val,'slope_val_after_treatment':slope_val_after_treatment,'slope_val_before_treatment':slope_val_before_treatment, 'area trapz':area, 'area simpson':area2, 'recovery_time':recovery_time,'time_to_recover_to_50':time_to_recover_to_50}

def scale_to_0(normalised_Peaks):
    for column in normalised_Peaks:
        col1 = normalised_Peaks[column]
        
        normalised_Peaks[column]=col1-col1.reset_index(drop=True)[0]
    return normalised_Peaks

def invert(d1):
    # this flips the peaks into valeys and vice versa, needed to be able to apply exactly same mathed as for Resistance.
    for i,val1 in d1.iteritems():
        if val1>0:
            val1 = float(f"-{val1}")
        elif val1<0:
            val1 = float(f"{val1}".replace('-',''))   
        d1[i]=val1
    return d1
 
def get_start_and_end(d,filter_width,window_size):
    
    injury_time_start,failed1 = get_injury_time(d,'min',filter_width,window_size)
    injury_time_start=injury_time_start-5
    d_reverse = d.iloc[::-1]
    d_reverse_reindex = d_reverse.reset_index()
    d_rev = d_reverse.reset_index(drop=True)
    
    arbitary_injury_time_end,failed2 = get_injury_time(d_rev,'max',filter_width*2,window_size*2-1)
    injury_time_end = int(d_reverse_reindex.iloc[arbitary_injury_time_end]['index'])
    if failed1 | failed2:
        # print("One of the measurements Failed")
        injury_time_start=0
        injury_time_end=2
    return injury_time_end,injury_time_start 
    
def peak_detection(d1):
    min_value = d1[d1==d1.min()].index[0]
    min1 = min_value-200
    if(min1<10):
        min1=10
    d=d1.iloc[min1:]
    # d.plot()

    filter_width=8
    window_size=25
    # We try 3 different filtering strategies to select the right peak, if first approach fails it will do different smoothening of the curve.
    # Sometimes particularly in wounding assays in low frequencies smoothening introduces a sharp peak, since there is a rappid capacitance/impedence drop of. 
    # Hence performing peak seletion on non-smoothened graph is better.
    # d.plot()
    # injury_time_end,injury_time_start = get_start_and_end(d,15,21)
    injury_time_end,injury_time_start = get_start_and_end(d,filter_width,window_size)
    if (injury_time_end-injury_time_start)<50:
        print('Failed 1st smoothening')
        # the filter was too wide, so try without it. 
        filter_width=3
        window_size=11
        injury_time_end,injury_time_start = get_start_and_end(d,filter_width,window_size)    
        if (injury_time_end-injury_time_start)<50:
            print('Failed 2nd smoothening')
            filter_width=1
            window_size=5
            injury_time_end,injury_time_start = get_start_and_end(d,filter_width,window_size)
            if (injury_time_end-injury_time_start)<50:
                print('Failed 3d smoothening, lets try super smoothened line')
                filter_width=4
                window_size=5
                injury_time_end,injury_time_start = get_start_and_end(d,filter_width,window_size)
                if (injury_time_end-injury_time_start)<50:
                    print('All trys failed smoothened line')
        
    Peak_window_only = d1.loc[injury_time_start:injury_time_end]
    # Peak_window_only.plot()

    return injury_time_start,injury_time_end,min_value

def run_analysis(measurement):
    df = pd.DataFrame()

    Data = pd.read_csv(f'Data3/{norm_method}/{treatment}_Data_all_data_remapped_{measurement}.csv',index_col=0)
    All_Experiment_Data_freq = pd.DataFrame()
    All_calculations_freq = pd.DataFrame()
    
    for freq1 in list(set(Data['freq'])):
        # freq1=62.5
        d2 = Data[Data['freq']==freq1]
        d2 = d2.drop('freq',axis=1)
        d2.reset_index(drop=True,inplace=True)
        d2 = d2.fillna(0)

        # Here loop through the spectra and calculate the slope at each point
        # d2['E520_2'].plot()
        idx_all = pd.Series(d2.columns)
        
        All_experiments = set(idx_all.str.split('_').str[0])
        All_Experiment_Data = pd.DataFrame()
        All_calculations = {}
        Experiment_injury_times = {}
        # All_experiments=['12e','1e','9e']
        for exp1 in All_experiments:
            # exp1='2e'
            # exp1='16e'
            # exp1='5e'
            # exp1='15e'
            # All the experiments performed together for this run
            exp = idx_all[idx_all.str.contains(f"^{exp1}_")]
            all_injury_times = []
            # Detect the aproximate injury time using the EMPTY wells
            EMPTY_WELLS = d2[list(exp[exp.str.contains('EMPTY WELL')])]
            # experiment 8e doesnt have an empty wells in them.
            
            # EMPTY_WELLS.plot()
            for id1 in exp:
                # id1 =idx_all[idx_all.str.contains('E1192_1')].values[0]
                # id1 =idx_all[idx_all.str.contains('E520_2')].values[0]
                # id1='12e_E590_2'
                # id1 = '12e_E588_2'
                # id1='9e_E513_1'
                # id1 ='14e_E659_1'
                # id1='14e_E694_2'
                # id1='4e_E587_2' - does not ever recover to 50%
                
                # id1='3e_E530_1' - this fails in a capacitance analysis
                # id1='3e_E528_2'
                # id1='16e_E1192_1'
                # id1='13e_E612_2' -- this one is interesting, since it kind of has 2 peaks for recovery
                # id1='16e_E1081_1'
                # id1='17e_E1533_1'
                # id1='17e_E1546_1'
                # id1='16e_E1081_1'
                # id1='16e_E1081_1'
                # id1='15e_E901_2'
                # id1='15e_E901_2'
                # id1='15e_E918_2'
                
                print(id1)
                if('EMPTY WEL' in id1):
                    continue
                
                d1 = d2.loc[d2[id1]!= 0,id1] 
                if measurement=='Capacitance':
                    d1 = invert(d1)
                d3=d1.reset_index()
                d1_index = d3['index']
                d1=d1.reset_index(drop=True)
                d1=d1[10:len(d1)]
                # d1.plot()
                
                injury_time_start,injury_time_end,min_value = peak_detection(d1)
                injury_time_start_norm = d1_index[injury_time_start]

                all_injury_times.append(injury_time_start_norm)
                try:
                    Experiment_injury_times[injury_time_start_norm].append(id1)
                except:
                    Experiment_injury_times[injury_time_start_norm]=[id1]
                
                Peak_window_only = d1.loc[injury_time_start:injury_time_end]
                
                # print(id1)
                try:
                    measurements = calculate_peak_metrics(Peak_window_only,injury_time_start,injury_time_end,min_value,d1)
                    measurements['flagged as outlier']=False
                except:
                    measurements = {'slope_of_recovery':0,'slope_val_after_treatment':0,'slope_val_before_treatment':0, 'area trapz':0, 'area simpson':0, 'recovery_time':0,'time_to_recover_to_50':0}
                    measurements['flagged as outlier']=True
                
                
                All_calculations[id1] = measurements
                # d1.plot()
                print('plotted')
        
            print('Done with this experiment')
            Consensous_injury_start = pd.Series(all_injury_times).mode()[0]
            
            # We use the consensous peak measurements to detect if any of the measurements may be outliers. 
            all_injury_times2 = pd.DataFrame(all_injury_times)
            Q1 = all_injury_times2.quantile(0.25)
            Q3 = all_injury_times2.quantile(0.75)
            IQR = Q3 - Q1
            if(IQR[0]==0):
                # we handle the exception detection when all are the same besides 1
                IQR = Consensous_injury_start*0.1
            Highest_bound = Q3 +  1.5*IQR
            Lowest_bond = Q1 - 1.5*IQR
            
            Outliers1 = all_injury_times2.loc[(all_injury_times2>Highest_bound)[0]]
            Outliers2 = all_injury_times2.loc[(all_injury_times2<Lowest_bond)[0]]
            Outliers=pd.concat([Outliers1,Outliers2])
            # These would have had a wrong calculations performed, and may need to be repeated. Hence we repeat the calculations based on consensous peak and flag the measurement as an outlier.
            for i,outlier in Outliers.iterrows():
                idx = (outlier.values[0])
                ids = Experiment_injury_times[idx]
                
                for id1 in ids:
                    try:
                        d1 = d2.loc[d2[id1]!= 0,id1] 
                        if measurement=='Capacitance':
                            d1 = invert(d1)
                        d3=d1.reset_index()
                        d1_index = d3['index']
                        d1=d1.reset_index(drop=True)
                        # d1.plot()
                        start1=Consensous_injury_start-200
                        end1=Consensous_injury_start+400
                        if(start1<0):
                           start1=0 
                        if(end1>d1.index[-1]):
                            end1=d1.index[-1]
                        d1=d1[start1:end1]
                        # d1.plot()
                        injury_time_start,injury_time_end,min_value = peak_detection(d1)
                        injury_time_start_norm = d1_index[injury_time_start]
                        # if the new injury start time fails to be in the correct place (±10% of Consensous injury time), there is something wrong experimentally and sample may need to be repeated.
                        Tolerance_Range = Consensous_injury_start*0.1
                        
                        if (Consensous_injury_start-Tolerance_Range<=injury_time_start_norm<=Consensous_injury_start+Tolerance_Range):
                            Peak_window_only = d1.loc[injury_time_start:injury_time_end]
                            measurements = calculate_peak_metrics(Peak_window_only,injury_time_start,injury_time_end,min_value,d1)
                        else:
                            measurements = {'slope_of_recovery':0,'slope_val_after_treatment':0,'slope_val_before_treatment':0, 'area trapz':0, 'area simpson':0, 'recovery_time':0,'time_to_recover_to_50':0}

                        measurements['flagged as outlier']=True
                        All_calculations[id1] = measurements
                    except:
                        print('Peak is an outlier')
                        measurements = {'slope_of_recovery':0,'slope_val_after_treatment':0,'slope_val_before_treatment':0, 'area trapz':0, 'area simpson':0, 'recovery_time':0,'time_to_recover_to_50':0}
                        measurements['flagged as outlier']=True
                        All_calculations[id1] = measurements

            normalised_Peaks = d2.loc[Consensous_injury_start:Consensous_injury_start+400,exp]
            # normalised_Peaks.plot()
            normalised_Peaks=normalised_Peaks.reset_index(drop=True)
            scaled_to_0_normalised_Peaks = scale_to_0(normalised_Peaks)
            # scaled_to_0_normalised_Peaks.plot()
            All_Experiment_Data = pd.concat([All_Experiment_Data,scaled_to_0_normalised_Peaks], axis=1)
            # All_Experiment_Data.plot()
            print('plotted')
        All_calculations = pd.DataFrame(All_calculations).T
        cols1 = pd.Series(All_Experiment_Data.columns)
        All_calculations['freq']= freq1
        All_Experiment_Data['freq']= freq1
        # slope_before_treatment
        # slope after treatment
        # All_Experiment_Data.plot()
        All_Experiment_Data_freq = pd.concat([All_Experiment_Data_freq,All_Experiment_Data], axis=0)
        All_calculations_freq = pd.concat([All_calculations_freq,All_calculations], axis=0)
        control_samples = All_Experiment_Data[cols1[cols1.str.contains('CONTROL')]]
        # control_samples.plot()
    # print(pd.DataFrame(slope_vals))
    # All_Experiment_Data_freq = pd.DataFrame()
    # All_calculations_freq = pd.DataFrame()
    All_Experiment_Data_freq.to_csv(f'Data3/{norm_method}/Data_Extracted_Peaks/{treatment}_Extracted_Peaks_{measurement}.csv')
    All_calculations_freq.to_csv(f'Data3/{norm_method}/Data_Extracted_Peaks/{treatment}_Metrics_Calculations_{measurement}.csv')       
    
def main():
    print('ļets do data extractions')
    measurements=['Resistance','Capacitance','Impedence']
    for measurement in measurements:
        run_analysis(measurement)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
    description="""
        This pipeline takes the normalised resampled ECIS data (produced by PCA_HUVEC.py) and calculates the slope.
        """
    )

    main()
