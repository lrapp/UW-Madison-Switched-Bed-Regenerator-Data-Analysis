def lin_int_cycle(df,BV,folder):
    
#old args    cycle_original,BV,folder
    import pandas as pd
    import os
    import dateutil.parser
    import numpy as np
    
    file_list=os.listdir(folder)
        

    if 'DFF.csv' in file_list:
        DFF=pd.read_csv(folder+"\\DFF.csv",index_col=[0],parse_dates=True)
        return DFF
    
    cols=df.columns
    
    unique_cntr=0

    DFF=pd.DataFrame()
    for i in range(0,len(BV)-2):
        
        df_i=df[df.index>=dateutil.parser.parse(str(BV.index[i])) ] #select all data with index greater than the ith valve swtiching start time 
        df_i=df_i[df_i.index<=dateutil.parser.parse(str(BV.index[i+1]))]
        
        df_iplus=df[df.index>=dateutil.parser.parse(str(BV.index[i+1])) ] #select all data with index greater than the ith+1 valve swtiching start time 
        df_iplus=df_iplus[df_iplus.index<=dateutil.parser.parse(str(BV.index[i+2]))]
        
        x_bv0=BV.index[i]
        x_bv1=BV.index[i+1]
        d0={}
        d1={}
        if i==0:
            for j in range(0,len(cols)):
                y0=df_i[cols[j]].values[0]                 #select the first value
                d0.update({cols[j]:y0})
                if cols[j]=='state':

                    y_bv1=df_i[cols[j]].values[-1]
                else:    

                    if len(df_iplus)!= 0:
                        y3=df_iplus[cols[j]].values[0]
                        y2=df_i[cols[j]].values[-1]
                        x2=df_i.index[-1]
                        x3=df_iplus.index[0]
                        
                        y_bv1=(y2*(x3-x_bv1)+y3*(x_bv1-x2))/(x3-x2) #interpolate
                    d1.update({cols[j]:y_bv1})

                
        else:
            df_iminus=df[df.index>=dateutil.parser.parse(str(BV.index[i-1])) ] #select all data with index greater than valve swtiching time 
            df_iminus=df_iminus[df_iminus.index<=dateutil.parser.parse(str(BV.index[i]))]

            for j in range(0,len(cols)):
                if cols[j]=='state':                    
                    y_bv0=df_i[cols[j]].values[-1]
                    y_bv1=df_i[cols[j]].values[-1]
                    d0.update({cols[j]:y_bv0})
                    d1.update({cols[j]:y_bv1})

                else:                   
                    y1=df_i[cols[j]].values[0]
                    y0=df_iminus[cols[j]].values[-1]
                    x0=df_iminus.index[-1]
                    x1=df_i.index[0]                   
                    
                    
                    y_bv0=(y0*(x1-x_bv0)+y1*(x_bv0-x0))/(x1-x0)
                    d0.update({cols[j]:y_bv0})
                    
                    if len(df_iplus)!= 0:
                        y3=df_iplus[cols[j]].values[0]
                        y2=df_i[cols[j]].values[-1]
                        x2=df_i.index[-1]
                        x3=df_iplus.index[0]
                        y_bv1=(y2*(x3-x_bv1)+y3*(x_bv1-x2))/(x3-x2)
                        d1.update({cols[j]:y_bv1})

          
    
        df_j0=pd.DataFrame(d0,index=[x_bv0]) #at x_bv0 index add d0 values
        df_j1=pd.DataFrame(d1,index=[x_bv1])
        df_j=df_j0.append(df_i) #append original cycle values to df_j0 and call it df_j
        df_j=df_j.append(df_j1)
        
        state_index_cntr=list(np.arange(0,len(df_j)))
        df_j['state_index_num']=state_index_cntr

        state_i=[]
        for ii in range(0,len(df_j)):
            state_i.append(BV.Position[i])
            
        df_j['state']=state_i
        unique_state_num=list(np.ones(len(df_j))*unique_cntr)
        df_j['unique_state_num']=unique_state_num
        unique_cntr=unique_cntr+1
        DFF=DFF.append(df_j)
        
    DFF.to_csv(folder+"\\DFF.csv")

    return DFF