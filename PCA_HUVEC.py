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
        return Wounding_Data1
    except:
        start=Frequency_indexes[Frequency_indexes.iloc[:,0].str.contains(frequency)].index[1]+1
        end = Frequency_indexes[Frequency_indexes.index>start].index[0]-1

        skiplist1=list(range(0,start))
        skiplist2 = (list(range(end,File_length)))
        # Wounding_Data1 = pd.read_csv(Wounding_File_Path,sep='XXXXXXX',skiprows=skiplist1+skiplist2)
        Wounding_Data1 = pd.read_csv(Wounding_File_Path,skiprows=skiplist1+skiplist2)
        return Wounding_Data1

def process_data(title,Prolif_Data):
    Impedence_Data_Prolif = Prolif_Data[[col for col in Prolif_Data.columns if title in col]]
    cols = Impedence_Data_Prolif.iloc[0]
    cols = [w.replace(' ', '') for w in cols]
    Impedence_Data_Prolif.drop(0,inplace=True)
    Impedence_Data_Prolif.columns=cols
    Impedence_Data_Prolif=Impedence_Data_Prolif.astype(float)
    return Impedence_Data_Prolif

def normalise_timescale(Impedence_Data_Prolif,Info_Data,Time_Data_Prolif):
    all_data= pd.DataFrame()
    for col in Impedence_Data_Prolif.columns:
        col2=col.replace(' ','')
        Name = Info_Data[Info_Data.iloc[:,0].str.replace('0','').str.contains(col2)].index
        i=1
        if Name.empty:
            Name = Info_Data[Info_Data.iloc[:,1].str.replace('0','').str.contains(col.replace(' ',''))].index
            i=2
        # Name=f"{col2}_{Name[0]}_{i}"
        Name=f"{Name[0]}_{i}"
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
        for freq in all_frequencies:
            print(freq)
            data_paths={}
            all_data_remapped_Impedence = pd.DataFrame()
            all_data_remapped_Resistance = pd.DataFrame()
            all_data_remapped_Capacitance = pd.DataFrame()
            for experiment in glob(f'{path}/*'):
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
                Data = read_ECIS_Data(File_Path,frequency=freq)

                Info_Data = pd.read_excel(Info_sheet_path,header=None,index_col=0).iloc[:,[0,1]]
                Info_Data=Info_Data.dropna(axis=0)
               
                
                Impedence_Data_Prolif = process_data('Imp.(ohm)',Data)
                Time_Data_Prolif = process_data('Time(hrs)',Data)
                Capacitance_Data_Prolif = process_data('Cap.(nF)',Data)
                Resistance_Data_Prolif = process_data('Res.(ohm)',Data)
                
                Impedence_Data_Prolif =normalise_timescale(Impedence_Data_Prolif,Info_Data,Time_Data_Prolif)
                Capacitance_Data_Prolif =normalise_timescale(Capacitance_Data_Prolif,Info_Data,Time_Data_Prolif)
                Resistance_Data_Prolif =normalise_timescale(Resistance_Data_Prolif,Info_Data,Time_Data_Prolif)
                all_data_remapped_Impedence = pd.concat([all_data_remapped_Impedence, Impedence_Data_Prolif], axis=1)
                all_data_remapped_Capacitance = pd.concat([all_data_remapped_Capacitance, Capacitance_Data_Prolif], axis=1)
                all_data_remapped_Resistance = pd.concat([all_data_remapped_Resistance, Resistance_Data_Prolif], axis=1)
            
            all_data_remapped_Resistance['freq']=freq
            all_data_remapped_Capacitance['freq']=freq
            all_data_remapped_Impedence['freq']=freq

            all_data_remapped_Resistance_all = pd.concat([all_data_remapped_Resistance_all,all_data_remapped_Resistance])
            all_data_remapped_Capacitance_all = pd.concat([all_data_remapped_Capacitance_all,all_data_remapped_Capacitance])
            all_data_remapped_Impedence_all = pd.concat([all_data_remapped_Impedence_all,all_data_remapped_Impedence])


        all_data_remapped_Resistance_all.to_csv(f'Data2/{data_type}_all_data_remapped_Resistance.csv')
        all_data_remapped_Capacitance_all.to_csv(f'Data2/{data_type}_all_data_remapped_Capacitance.csv')
        all_data_remapped_Impedence_all.to_csv(f'Data2/{data_type}_all_data_remapped_Impedence.csv')
        all_data_remapped_Capacitance.plot()
    print('Done')
    # import matplotlib.pyplot as plt
    # plt.plot(all_data_remapped_Capacitance,label=all_data_remapped_Capacitance.columns)

    # plt.show()
    # from ppca import PPCA
    # ppca = PPCA()
    # ppca.fit(data=all_data_remapped, d=100, verbose=True)
    # variance_explained = ppca.var_exp
    # components = ppca.data
    # model_params = ppca.C



if __name__ == '__main__':

    parser = argparse.ArgumentParser(
    description="""
        This pipeline takes the ECIS data recordings and normalises it all to the 60s intervals.
        """
    )

    main()