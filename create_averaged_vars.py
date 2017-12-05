def create_averaged_vars(Temp,Pressure,folder):
    import pandas as pd
    import os
    from dateutil import parser
    
    file_list=os.listdir(folder) #List the files in the directory
    
    end_search=""
    file_found=""
    
    for file in file_list:
        if file == "df.csv":
            end_search=True
        if file == "Pressure_ave.csv" or file == "Temp.ave.csv":
            file_found="True"
            
    if end_search!=True:
        if file_found=="True":
            Pressure_average=pd.read_csv(folder+"\\Pressure_ave.csv", index_col=[0])    
            Pressure_average.index=Pressure_average.index.map(parser.parse)
            Temp_average=pd.read_csv(folder+"\\Temp_ave.csv", index_col=[0])
            Temp_average.index=Temp_average.index.map(parser.parse)
        else:
            Temp_average=Temp[0:1]
            Pressure_average=pd.DataFrame()
            old_k=0
            k=0
            for i in range(0,len(Temp)-1):
                print(i)
                Temp_i_index=Temp.index[i]
                Pressure_i_index=Pressure.index[k]
                while Pressure_i_index<=Temp_i_index:
                    k=k+1
                    Pressure_i_index=Pressure.index[k]
            
                df=pd.DataFrame(Pressure[old_k:k].mean().to_dict(),index=[Temp_i_index])
                Pressure_average=Pressure_average.append(df)
                Temp_average=Temp_average.append(Temp[i:i+1])
                old_k=k
        df=pd.concat([Temp_average,Pressure_average],axis=1)
        df.to_csv(folder+"\\df.csv")
    else: #if df.csv is found, read in file
        df=pd.read_csv(folder+"\\df.csv",index_col=[0],parse_dates=True)

    return [df]

