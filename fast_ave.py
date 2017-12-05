def fast_ave(start,end,df_full_cols):
    import numpy as np
    import pandas as pd

    st1_list=[]
    st2a_list=[]
    st2b_list=[]
    st2c_list=[]
    st3_list=[]
    st4a_list=[]
    st4b_list=[]
    st4c_list=[]
    
    ave_list=[st1_list,st2a_list,st2b_list,st2c_list,st3_list,st4a_list,st4b_list,st4c_list]
    possible_states=['STATE1','STATE2a','STATE2b','STATE2c','STATE3','STATE4a','STATE4b','STATE4c']
    
    for i in range(start,end):
        num=len(df_full_cols.loc[df_full_cols['unique_state_num']==i])
        for k in range(0,len(possible_states)):
            if df_full_cols.loc[df_full_cols['unique_state_num']==i]['state'][0] == possible_states[k]:
                df_i=df_full_cols.loc[df_full_cols['unique_state_num']==i]
                time_d=[]
                for j in range(0,num):
                    time_d.append((df_i.index[j]-df_i.index[0]).total_seconds())
                df_i.insert(0,'rel_time',time_d)
                ave_list[k].append(df_i)
                

    min_list=[]
    for i in range(0,len(ave_list)):
        min_list.append(len(ave_list[i]))
        
    min_array=np.asarray(min_list)
    min_length=np.min(min_array[np.nonzero(min_array)])   
    
    cut_ave_list=[]
    for i in range(0,len(ave_list)):
        cut_ave_list.append(ave_list[i][0:min_length])

    mean_list=[]
    std_list=[]
    for i in range(0,len(ave_list)):
        if len(cut_ave_list[i])!=0:

            mean_list.append(pd.concat(cut_ave_list[i]).groupby('state_index_num').mean())
            std_list.append(pd.concat(cut_ave_list[i]).groupby('state_index_num').std())            
            k=0
            if len(mean_list[-1])!=len(cut_ave_list[i][0]['state'].tolist()):

                while  k < len(cut_ave_list[i]) and len(mean_list[-1])!=len(cut_ave_list[i][k]['state'].tolist()):
                    k=k+1

            mean_list[-1]['state']=cut_ave_list[i][k]['state'].tolist()
    #        std_list[-1]['state']=cut_ave_list[i][0]['state'].tolist()
    
    for i in range(0,len(std_list)):
        del mean_list[i]['unique_state_num']
        del std_list[i]['unique_state_num']
        new_cols=[x+'_ave' for x in mean_list[i].columns[0:-1].tolist()]+[mean_list[i].columns[-1]]
    #    new_cols_std=[x+'_std' for x in mean_list[i].columns[0:-1].tolist()]+[mean_list[i].columns[-1]]
        new_cols_std=[x+'_std' for x in mean_list[i].columns[0:-1].tolist()]
        mean_list[i].columns=new_cols
        std_list[i].columns=new_cols_std
    
    mean_list_state=[]
    for i in range(0,len(mean_list)):
        mean_list_state.append(mean_list[i]['state'][0])
    
    average_state_list=[]
    for i in range(0,len(mean_list)):
        average_state_list.append(pd.concat((mean_list[i],std_list[i]),axis=1))
    
    for i in range(1,len(average_state_list)):
        if average_state_list[i-1]['TI16_ave'][-1:].values[0] == average_state_list[i]['TI16_ave'][0:].values[0]:
            average_state_list[i]=average_state_list[i].drop(0)
            average_state_list[i].reset_index()
        
    #clean up 
    full=pd.concat(average_state_list,ignore_index=True)
    full1=full
    for i in range(1,len(full1)):
        if full1['rel_time_ave'][i]==0:
            full1=full1.drop(i)
    cols=full1.columns.tolist()
    st=cols.pop(cols.index('state'))
    cols2=cols[0:1]+[st]+cols[1:-1]
    full1=full1[cols2]
    full1=full1.reset_index()
    del full1['index']
      
    #check if any states are empty, if so do not include them in the final full dataframe
    asl2=[]
    if len(full.loc[full['state']=='STATE1'])!=0:
        st1=full.loc[full['state']=='STATE1']    
        asl2.append(st1)
    if len(full.loc[full['state']=='STATE2a'])!=0:
        st2a=full.loc[full['state']=='STATE2a']
        asl2.append(st2a)
    if len(full.loc[full['state']=='STATE2b'])!=0:        
        st2b=full.loc[full['state']=='STATE2b']   
        asl2.append(st2b)
    if len(full.loc[full['state']=='STATE2c'])!=0:        
        st2c=full.loc[full['state']=='STATE2c']    
        asl2.append(st2c)
    if len(full.loc[full['state']=='STATE3'])!=0:        
        st3=full.loc[full['state']=='STATE3']    
        asl2.append(st3)
    if len(full.loc[full['state']=='STATE4a'])!=0:        
        st4a=full.loc[full['state']=='STATE4a']    
        asl2.append(st4a)
    if len(full.loc[full['state']=='STATE4b'])!=0:        
        st4b=full.loc[full['state']=='STATE4b'] 
        asl2.append(st4b)
    if len(full.loc[full['state']=='STATE4c'])!=0:        
        st4c=full.loc[full['state']=='STATE4c'] 
        asl2.append(st4c)

        
    asl3=[]
    asl3.append(asl2[0].copy())
    for i in range(1,len(asl2)):
        asl3_i=asl2[i].copy()
        indexer=asl3_i.index.values
        new_vals=asl3_i['rel_time_ave'].values+asl2[i-1]['rel_time_ave'][-1:].values[0] 
        asl3_i.loc[indexer,'rel_time_ave']=new_vals     
        asl3.append(asl3_i)
        
    full2=pd.concat(asl3,ignore_index=True)
    cols=full2.columns.tolist()
    st=cols.pop(cols.index('state'))
    cols2=cols[0:1]+[st]+cols[1:-1]
    full2=full2[cols2]
        
    S1_ave=full.loc[full['state']=="STATE1"]
    S2a_ave=full.loc[full['state']=="STATE2a"]
    S2b_ave=full.loc[full['state']=="STATE2b"]
    S2c_ave=full.loc[full['state']=="STATE2c"]
    S3_ave=full.loc[full['state']=="STATE3"]
    S4a_ave=full.loc[full['state']=="STATE4a"]
    S4b_ave=full.loc[full['state']=="STATE4b"]
    S4c_ave=full.loc[full['state']=="STATE4c"]
  
    return [S1_ave,S2a_ave,S2b_ave,S2c_ave,S3_ave,S4a_ave,S4b_ave,S4c_ave,full2]
        
        

