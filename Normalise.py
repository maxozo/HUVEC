#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2021-11-19'
__version__ = '0.0.1'

from glob import glob
import pandas as pd
import argparse

def read_ECIS_Data(Wounding_File_Path,frequency=250):
    frequency=str(frequency)
    

    with open(Wounding_File_Path) as f:
        lines = f.readlines()
    Wounding_Data = pd.DataFrame(lines)
    Date = Wounding_Data[Wounding_Data.iloc[:,0].str.startswith('Date ')][0]
    Date_time = Date.iloc[0].split('Date , ')[1].replace('\n','')
    Date_time = int(Date_time.replace('-','').replace(',','').replace(' ','').replace(':',''))
    Frequency_indexes = Wounding_Data[Wounding_Data.iloc[:,0].str.startswith('Frequency')]
    File_length = len(Wounding_Data)+60
    try:
        try:
            start=Frequency_indexes[Frequency_indexes.iloc[:,0].str.contains(frequency)].index[1]+1
            try:
                end = Frequency_indexes[Frequency_indexes.index>start].index[0]-1
            except:
                print('Must be the last frequency in the file')
                end=len(Wounding_Data)
        except:
            print('only 1 freq available')
            start=Frequency_indexes[Frequency_indexes.iloc[:,0].str.contains(frequency)].index[0]+1
            end=len(Wounding_Data)
        
        

        skiplist1=list(range(0,start))
        skiplist2 = (list(range(end,File_length)))
        # Wounding_Data1 = pd.read_csv(Wounding_File_Path,sep='XXXXXXX',skiprows=skiplist1+skiplist2)
        Wounding_Data1 = pd.read_csv(Wounding_File_Path,skiprows=skiplist1+skiplist2)
        return Wounding_Data1,Date_time
    except:
        start=Frequency_indexes[Frequency_indexes.iloc[:,0].str.contains(frequency)].index[1]+1
        end = Frequency_indexes[Frequency_indexes.index>start].index[0]-1

        skiplist1=list(range(0,start))
        skiplist2 = (list(range(end,File_length)))
        # Wounding_Data1 = pd.read_csv(Wounding_File_Path,sep='XXXXXXX',skiprows=skiplist1+skiplist2)
        Wounding_Data1 = pd.read_csv(Wounding_File_Path,skiprows=skiplist1+skiplist2)
        return Wounding_Data1,Date_time

def process_data(title,Prolif_Data):
    Impedence_Data_Prolif = Prolif_Data[[col for col in Prolif_Data.columns if title in col]]
    cols = Impedence_Data_Prolif.iloc[0]
    cols = [w.replace(' ', '') for w in cols]
    Impedence_Data_Prolif.drop(0,inplace=True)
    Impedence_Data_Prolif.columns=cols
    Impedence_Data_Prolif=Impedence_Data_Prolif.astype(float)
    return Impedence_Data_Prolif

def normalise_timescale(Impedence_Data_Prolif,Info_Data,Time_Data_Prolif,experiment_count):
    all_data= pd.DataFrame()
    for col in Impedence_Data_Prolif.columns:
        col2=col.replace(' ','')
        Name = Info_Data[Info_Data.iloc[:,0].str.replace('0','').str.contains(col2)].index
        sample_replicate=1
        if Name.empty:
            Name = Info_Data[Info_Data.iloc[:,1].str.replace('0','').str.contains(col.replace(' ',''))].index
            sample_replicate=2
        # Name=f"{col2}_{Name[0]}_{i}"
        Name=f"{experiment_count}e_{Name[0]}_{sample_replicate}"
        Impedence_Data_Prolif.rename(columns={col:Name},inplace=True)
        data = pd.concat([pd.Series(0),Impedence_Data_Prolif[Name]])
        time =pd.concat([pd.Series(0),Time_Data_Prolif[col2]])

        data.index =pd.to_timedelta(time, unit='hr') #.dt.floor('s')
        resampled=data.resample('60S').pad()
        Data = pd.DataFrame(resampled)
        Data.columns=[Name]
        all_data=pd.concat([all_data,Data ], axis=1)
        
    d2 = pd.DataFrame(all_data)
    return d2

def get_all_frequencies(data_type,path):
    data_paths={}
    # Use First file in dictionary to listy all the frequencies recorded.
    experiment = glob(f'{path}/*')[0]
    Wounding_File_Path = (glob(f"{experiment}/*WOUND*"))[0]
    Thrombin_File_Path = (glob(f"{experiment}/*THROMBIN*"))[0]
    try:
        Proliferation_File_Path = (glob(f"{experiment}/*PROLIF*"))[0]
    except:
        Proliferation_File_Path = (glob(f"{experiment}/*prolif*"))[0]
    data_paths['Wounding_Data']=Wounding_File_Path
    data_paths['Thrombin_Data']=Thrombin_File_Path
    data_paths['Prolif_Data']=Proliferation_File_Path
    
    data_path = data_paths[data_type]
    with open(data_path) as f:
        lines = f.readlines()
    Wounding_Data = pd.DataFrame(lines)
    Frequency_indexes = Wounding_Data[Wounding_Data.iloc[:,0].str.startswith('Frequency')]
    all_frequencies = Frequency_indexes.values[0][0].replace('Frequency , ','').replace(', \n','').split(', ')
    return all_frequencies

def Impedence_Data_prelog_normalise(Impedence_Data_pre):
    # As Kevin Suggested this normalisation may help in teasing out the actual differences in data.
    CONTROLS = Impedence_Data_pre[list(Impedence_Data_pre.columns[Impedence_Data_pre.columns.str.contains('CONTROL')])]
    NOT_CONTROLS = Impedence_Data_pre[list(Impedence_Data_pre.columns[~Impedence_Data_pre.columns.str.contains('CONTROL')])]
    Mean_CONTROLS = CONTROLS.mean(axis=1)
    
    for column in Impedence_Data_pre:
        col1 = Impedence_Data_pre[column]
        
        # Here we can also add in the controls to get the control scaled normalisation
        col1_and_CONTROLS = Mean_CONTROLS.to_frame().copy()
        col1_and_CONTROLS['Sample'] = col1
        col1_and_CONTROLS = col1_and_CONTROLS.loc[~(col1_and_CONTROLS.sum(axis=1)==0)] # We remove the initial 0 to calculate the scaling factors, as this was arbitately added
        # Mean Normalisation
        mean_norm=(col1-col1.mean())/col1.std()
        mean_norm_with_controls=(col1-col1_and_CONTROLS.mean().mean())/col1_and_CONTROLS.std().std()
        
        # Min-Max Normalisation
        min_max_norm=(col1-col1.min())/(col1.max()-col1.min())
        min_max_norm_with_controls=(col1-col1_and_CONTROLS.min().min())/(col1_and_CONTROLS.max().max()-col1_and_CONTROLS.min().min())
        
        # Internal FC normalisation
        Norm_Value = col1.loc[~(col1==0)].iloc[0]
        log_norm = col1/Norm_Value
        # FC Against Control Mean intial value
        Norm_Value = Mean_CONTROLS.loc[~(Mean_CONTROLS==0)].iloc[0]
        log_control_norm = col1/Norm_Value
        
        # log_control_norm.plot()
        
        # Replace the Data with norm version
        Impedence_Data_pre[column]=min_max_norm_with_controls
        
    return Impedence_Data_pre

def normalise_against_Controls(Data):
    CONTROLS = Data[list(Data.columns[Data.columns.str.contains('CONTROL')])]
    NOT_CONTROLS = Data[list(Data.columns[~Data.columns.str.contains('CONTROL')])]
    Mean_CONTROLS = CONTROLS.mean(axis=1)
    Norm_Value = Mean_CONTROLS.loc[~(Mean_CONTROLS==0)][0]
    Data1 = NOT_CONTROLS['1e_E521_1']
    Data1.plot()
    Data_Normalised = Data1/Norm_Value
    Data_Normalised1 = Data1/Mean_CONTROLS
    Data_Normalised_not = Data1
    Data_Normalised_not.plot()
    Data_Normalised.plot()
    Mean_CONTROLS.plot()
    Data_Normalised1.plot()
    print('Done')

def main():
    path='/Volumes/GoogleDrive/My Drive/HUVEC/ECIS/ECIS Data'
    i=0

    data_types=['Thrombin_Data','Wounding_Data','Prolif_Data']

    # Have to loop through - 
    # 1) Experiments 
    # 2) Treatments - Wounding, Thrombin, Proliferation
    # 3) Frequencies
    # Will generate 3 files for each treatment, containing a scaled frequency profiles for each
    for data_type in data_types:
        print(data_type)
        all_frequencies = get_all_frequencies(data_type,path)
        print(all_frequencies)
        all_data_remapped_Impedence_all = pd.DataFrame()
        all_data_remapped_Resistance_all = pd.DataFrame()
        all_data_remapped_Capacitance_all = pd.DataFrame()
        all_Data_times_all=pd.DataFrame()
        for freq in all_frequencies:
            print(freq)
            data_paths={}
            all_data_remapped_Impedence = pd.DataFrame()
            all_data_remapped_Resistance = pd.DataFrame()
            all_data_remapped_Capacitance = pd.DataFrame()
            all_Data_times=pd.DataFrame()
            experiment_count = 0
            for experiment in glob(f'{path}/*'):
                experiment_count+=1
                i+=1
                print(experiment)
                Wounding_File_Path = (glob(f"{experiment}/*WOUND*"))[0]
                Thrombin_File_Path = (glob(f"{experiment}/*THROMBIN*"))[0]
                try:
                    Proliferation_File_Path = (glob(f"{experiment}/*PROLIF*"))[0]
                except:
                    Proliferation_File_Path = (glob(f"{experiment}/*prolif*"))[0]
                Info_sheet_path =  (glob(f"{experiment}/*INFORMATION*"))[0]


                data_paths['Wounding_Data']=Wounding_File_Path
                data_paths['Thrombin_Data']=Thrombin_File_Path
                data_paths['Prolif_Data']=Proliferation_File_Path
                File_Path= data_paths[data_type]
                Data,Date_time = read_ECIS_Data(File_Path,frequency=freq)

                Info_Data = pd.read_excel(Info_sheet_path,header=None,index_col=0).iloc[:,[0,1]]
                Info_Data=Info_Data.dropna(axis=0)
               
                
                Impedence_Data_pre = process_data('Imp.(ohm)',Data)
                Time_Data_pre = process_data('Time(hrs)',Data)
                Capacitance_Data_pre = process_data('Cap.(nF)',Data)
                Resistance_Data_pre = process_data('Res.(ohm)',Data)
            
                Impedence_Data = normalise_timescale(Impedence_Data_pre,Info_Data,Time_Data_pre,experiment_count)
                Capacitance_Data =normalise_timescale(Capacitance_Data_pre,Info_Data,Time_Data_pre,experiment_count)
                Resistance_Data =normalise_timescale(Resistance_Data_pre,Info_Data,Time_Data_pre,experiment_count)
                
                Impedence_Data = Impedence_Data_prelog_normalise(Impedence_Data)
                Capacitance_Data = Impedence_Data_prelog_normalise(Capacitance_Data)
                Resistance_Data = Impedence_Data_prelog_normalise(Resistance_Data)
                
                # Here do the normalisations against the control samples.
                # Impedence_Data = normalise_against_Controls(Impedence_Data)
                # Capacitance_Data = normalise_against_Controls(Capacitance_Data)
                # Resistance_Data = normalise_against_Controls(Resistance_Data)
                
                all_data_remapped_Impedence = pd.concat([all_data_remapped_Impedence, Impedence_Data], axis=1)
                all_data_remapped_Capacitance = pd.concat([all_data_remapped_Capacitance, Capacitance_Data], axis=1)
                all_data_remapped_Resistance = pd.concat([all_data_remapped_Resistance, Resistance_Data], axis=1)
                Data_covar = pd.DataFrame(Impedence_Data.columns)
                Data_covar['Date_time']=Date_time
                Data_covar['data_type']=data_type
                Data_covar.rename(columns={0:'Sample'},inplace=True)
                Data_covar['Group']='Sample'
                Data_covar.loc[Data_covar['Sample'].str.contains('CONTROL'),'Group']='CONTROL'
                Data_covar.loc[Data_covar['Sample'].str.contains('EMPTY WELL'),'Group']='EMPTY_WELL'
                Data_covar['label']=Data_covar['Sample'].str.split('_').str[1]
                all_Data_times = pd.concat([all_Data_times,Data_covar],axis=0)
                
                # double_check_if rescaling dont mess up data
                # Impedence_Data_pre['1e_E522_1'].reset_index()['1e_E522_1'].plot()
                # Impedence_Data['1e_E522_1'].reset_index()['1e_E522_1'].plot()
                # Log_Impedence_Data_pre['1e_E522_1'].reset_index()['1e_E522_1'].plot()
                print('checked')
                
            all_data_remapped_Resistance['freq']=freq
            all_data_remapped_Capacitance['freq']=freq
            all_data_remapped_Impedence['freq']=freq

            all_data_remapped_Resistance_all = pd.concat([all_data_remapped_Resistance_all,all_data_remapped_Resistance])
            all_data_remapped_Capacitance_all = pd.concat([all_data_remapped_Capacitance_all,all_data_remapped_Capacitance])
            all_data_remapped_Impedence_all = pd.concat([all_data_remapped_Impedence_all,all_data_remapped_Impedence])
            all_Data_times_all=pd.concat([all_Data_times_all,all_Data_times])

        all_Data_times_all_no_dubs = all_Data_times_all.drop_duplicates()
        all_Data_times_all_no_dubs['exp_id'] = all_Data_times_all_no_dubs.Sample.str.split('_').str[0]
        
        # all_Data_times_all_no_dubs.to_csv(f'Data3/{data_type}_metadata.tsv',sep='\t',index=False)
        all_data_remapped_Resistance_all.to_csv(f'Data3/min_max_controls/{data_type}_all_data_remapped_Resistance.csv')
        all_data_remapped_Capacitance_all.to_csv(f'Data3/min_max_controls/{data_type}_all_data_remapped_Capacitance.csv')
        all_data_remapped_Impedence_all.to_csv(f'Data3/min_max_controls/{data_type}_all_data_remapped_Impedence.csv')
        all_data_remapped_Capacitance.plot()
    print('Done')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
    description="""
        This pipeline takes the ECIS data recordings and normalises it all to the 60s intervals.
        """
    )

    main()