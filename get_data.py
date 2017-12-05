##Checks if .csv files have been created from raw data files. If they exist 
##then they are returned, if they do not exist, they are created. 

def get_data(folder):
    import datetime as datetime
    import pandas as pd
    import os
    
    file_list=os.listdir(folder)
    T=False
    P=False
    B=False
    
    ##Check if .csv files have already been created
    for file in file_list:
        if file == "Temp.csv":
            Temp=pd.read_csv(folder+"\\Temp.csv",index_col=[0],parse_dates=True)
            T=True
        if file == "Pressure.csv":
            Pressure=pd.read_csv(folder+"\\Pressure.csv",index_col=[0],parse_dates=True)
            P=True
        if file == "BV.csv" :
            BV=pd.read_csv(folder+"\\BV.csv",index_col=[0],parse_dates=True)
            B=True
            
    if T == True and P== True and B==True:
        return Temp, Pressure, BV ##If all .csv files exist, return them
    else:   #if the .csv files do not exist, create them
        root=folder
        files=["Pressure","Temp","BV"]
        for k in files:
            filepath=root+"\\"+k+".txt"
            raw_data=pd.read_csv(filepath,sep='\t')
            date_columns=raw_data.columns[:7]
    #        date_times=[0 for x in range(len(date_columns))]
            form_date=[0 for x in range(len(raw_data))]
            for i in range(0,len(raw_data)):
                v_time=()
                for j in range(0,len(date_columns)):
                    if j == 6:
                        v_time=v_time+(int(raw_data[date_columns[6]][i]*1000000),)
                    else:
                        v_time=v_time+(raw_data[date_columns[j]][i],)
                j=j+1
                form_date[i]=datetime.datetime(*v_time[0:7])
            i=i+1
            
            pnames=raw_data.columns
            if k == files[0]:
                Pressure=pd.DataFrame(index=form_date) 
                for i in range(7,len(pnames)):
                    Pressure[pnames[i]]=raw_data[pnames[i]].values
                    i=i+1
            if k == files[1]:
                Temp=pd.DataFrame(index=form_date)  
                for i in range(7,len(pnames)):
                    Temp[pnames[i]]=raw_data[pnames[i]].values
                    i=i+1
        
            if k == files[2]:
                BV=pd.DataFrame(index=form_date) 
                BV[pnames[7]]=raw_data[pnames[7]].values
    
        Temp.to_csv(folder+"\\Temp.csv")
        Pressure.to_csv(folder+"\\Pressure.csv")
        BV.to_csv(folder+"\\BV.csv")
        
        return Temp, Pressure, BV

    