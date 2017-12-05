##---- Main ----####
from class_def import props,Q,HT_results,main_char,Timer
        
with Timer():
    from get_data import get_data #function used to import raw data
    from create_averaged_vars import create_averaged_vars #function to average pressure data so there are equal number of temperature and pressure points
    from split_cycles import split_cycles
    from lin_int_cycle2 import lin_int_cycle2
    from fast_ave2 import fast_ave2
    from ave_for_fit2 import ave_for_fit2
    from calc_Q2 import calc_Q2
    from heat_transfer_coef_fn import ht_coef_fn
    import props5
    import os
    
    refresh_data=True
    refresh_sub=False
    calc = True
    
    file_date="10_9_17" #Specify file date    
#    root="C:\\Users\\lrapp\\OneDrive - UW-Madison\\Research\\Data Store\\Data\\"
    root="C:\\Users\\Logan\\OneDrive - UW-Madison\\Research\\Data Store\\Data\\"    
    
    folder=root+file_date
    file_list=os.listdir(folder)
    
    if refresh_data==True:
        print("Getting Data...")        
        with Timer():
            [Temp,Pressure,BV]=get_data2(folder)
        print("Get Data Complete")
        print("Averaging Data...")
        with Timer():
            [df]=create_averaged_vars(Temp,Pressure,folder)
        print("Averaged Data Complete")
        
    print("Interpolating Data...")        
    with Timer():
        df_lin=lin_int_cycle2(df,BV,folder)        
    print("Interpolating Complete")       

    print("Getting thermophysical properties...")
    with Timer():
        df_full_cols=split_cycles2(df_lin,BV,folder)
    
#    print("correcting dp")
#    with Timer():
#        L_top=5*0.0254
#        L_bottom=14*0.0254
#        g=9.81
#        T=40
#        dp_corrected=[]
#        for i in range(0,len(df_full_cols)):
#            Conditions=props(props5.f90wrap_tp(T+273.15,df_full_cols['PT01'][i]*6.89475729))
#            dp_top_i=(Conditions.density*g*L_top)*0.000145038
#            dp_bottom_i=(Conditions.density*g*L_bottom)*0.000145038
#            dp_corrected.append(df_full_cols['DP01'][i]-dp_top_i-dp_bottom_i)
##        
#    del df_full_cols['DP01']
#    df_full_cols['DP01']=dp_corrected
    
    print("getting averages")
    if calc == True:
        start=387   
        end=445
        with Timer():
            [S1_ave,S2a_ave,S2b_ave,S2c_ave,S3_ave,S4a_ave,S4b_ave,S4c_ave,full_cycle]=fast_ave2(start,end,df_full_cols)
        
        
        Regen_Results=HT_results(ave_for_fit2(full_cycle,1,1.2))
        
        q_dot_max_RE1=Regen_Results.q_dot_max_RE1
        q_dot_max_RE2=Regen_Results.q_dot_max_RE2

        Q_re=Q(calc_Q2(start,end,df_full_cols))
        
        MC=main_char()
        
        MC=ht_coef_fn(start,end,df_full_cols)
        
    