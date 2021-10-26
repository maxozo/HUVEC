from glob import glob
import pandas as pd

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
        if Name.empty:
            Name = Info_Data[Info_Data.iloc[:,1].str.replace('0','').str.contains(col.replace(' ',''))].index
        Name=f"{col2}_{Name[0]}_{i}"
        Impedence_Data_Prolif.rename(columns={col:Name},inplace=True)
        data = pd.concat([pd.Series(0),Impedence_Data_Prolif[Name]])
        time =pd.concat([pd.Series(0),Time_Data_Prolif[col2]])

        data.index =pd.to_timedelta(time, unit='hr') #.dt.floor('s')
        resampled=data.resample('60S').pad()
        all_data=pd.concat([all_data, pd.DataFrame(resampled)], axis=1)
        
    d2 = pd.DataFrame(all_data)
    return d2

path='/Volumes/GoogleDrive/My Drive/HUVEC/ECIS/ECIS Data'
i=0
all_data_remapped_Impedence = pd.DataFrame()
all_data_remapped_Resistance = pd.DataFrame()
all_data_remapped_Capacitance = pd.DataFrame()

for experiment in glob(f'{path}/*'):
    i+=1
    # print(experiment)
    Wounding_File_Path = (glob(f"{experiment}/*WOUND*"))[0]
    Thrombin_File_Path = (glob(f"{experiment}/*THROMBIN*"))[0]
    try:
        Proliferation_File_Path = (glob(f"{experiment}/*PROLIF*"))[0]
    except:
        Proliferation_File_Path = (glob(f"{experiment}/*prolif*"))[0]
    Info_sheet_path =  (glob(f"{experiment}/*INFORMATION*"))[0]
    print(Wounding_File_Path)
    print(Thrombin_File_Path)
    print(Info_sheet_path)
    print(Proliferation_File_Path)
    print('------------------------')
    Wounding_Data = read_ECIS_Data(Wounding_File_Path,frequency=64000)
    Thrombin_Data = read_ECIS_Data(Thrombin_File_Path,frequency=4000)
    Prolif_Data = read_ECIS_Data(Proliferation_File_Path,frequency=4000)
    Info_Data = pd.read_excel(Info_sheet_path,header=None,index_col=0).iloc[:,[0,1]]
    Info_Data=Info_Data.dropna(axis=0)
    Prolif_Data=Thrombin_Data
    
    Impedence_Data_Prolif = process_data('Imp.(ohm)',Prolif_Data)
    Time_Data_Prolif = process_data('Time(hrs)',Prolif_Data)
    Capacitance_Data_Prolif = process_data('Cap.(nF)',Prolif_Data)
    Resistance_Data_Prolif = process_data('Res.(ohm)',Prolif_Data)
    
    Impedence_Data_Prolif =normalise_timescale(Impedence_Data_Prolif,Info_Data,Time_Data_Prolif)
    Capacitance_Data_Prolif =normalise_timescale(Capacitance_Data_Prolif,Info_Data,Time_Data_Prolif)
    Resistance_Data_Prolif =normalise_timescale(Resistance_Data_Prolif,Info_Data,Time_Data_Prolif)
    all_data_remapped_Impedence = pd.concat([all_data_remapped_Impedence, Impedence_Data_Prolif], axis=1)
    all_data_remapped_Capacitance = pd.concat([all_data_remapped_Capacitance, Capacitance_Data_Prolif], axis=1)
    all_data_remapped_Resistance = pd.concat([all_data_remapped_Resistance, Resistance_Data_Prolif], axis=1)

print(all_data_remapped_Resistance)

import matplotlib.pyplot as plt
plt.plot(all_data_remapped_Resistance)
plt.show()
from ppca import PPCA
ppca = PPCA()
ppca.fit(data=all_data_remapped, d=100, verbose=True)
variance_explained = ppca.var_exp
components = ppca.data
model_params = ppca.C