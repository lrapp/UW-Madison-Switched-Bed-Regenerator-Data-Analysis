def split_cycles(df,BV,folder):
    import dateutil.parser
    import pandas as pd
    import os
    import numpy as np
    import props5
    from class_def import props as props_vals
    
    def k1_out(Re):    
        import numpy as np
        a=-2.341
        b=-1.595e-5
        c=-1.128
        d=-2.307e-8
        k=a*np.exp(b*Re)+c*np.exp(d*Re)
        return -1*k
        
    def k2_in(Re):    
        import numpy as np
        a=5.221
        b=-4.216e-7
        c=-4.255
        d=-6.56e-7
        k=a*np.exp(b*Re)+c*np.exp(d*Re)
        return k
        
    def k1_in(Re):    
        import numpy as np
        a=-0.04707
        b=1.465e-6
        c=1.098
        d=2.552e-7
        k=a*np.exp(b*Re)+c*np.exp(d*Re)
        return k
        
    def k2_out(Re):    
        import numpy as np
        a=-0.8439
        b=-1.522e-5
        c=-1.174
        d=2.812e-8
        k=a*np.exp(b*Re)+c*np.exp(d*Re)
        return -1*k
    

    #Initialize Variables:
    cycle_st_time=[]
    cycle=[]
    H=pd.DataFrame()
    DFF_state=pd.DataFrame()
    
    file_list=os.listdir(folder)
#    for file in file_list:  
#        if 'cycle' in file_list:    #Check if 'cycle' folder exists in directory
#            if len(os.listdir(folder+"\\cycle"))!=0: #if folder exists and is not empty, read contents of folder into 'cycle'
#                cycle=[]
#                for i in range(0,len(os.listdir(folder+"\\cycle"))):
#                    cycle.append(pd.read_hdf(folder+"\\cycle\\cycle"+str(i)+".h5",'table'))
#                cycle_st_time=[]
#                for i in range(0,len(cycle)):
#                    cycle_st_time.append(cycle[i].index[0])
#                    
#                return [cycle , cycle_st_time]
    if "df_full_cols.h5" in file_list:
#        df_state = pd.read_csv(folder+"\\df_full_cols.h5",index_col=[0],parse_dates=True)
        df_state = pd.read_hdf(folder+"\\df_full_cols.h5",'table')
        return df_state
        
#        elif "cycle_store.h5" in file_list:
#            file_name=folder+"\\cycle_store.h5"
#            cycle_store=pd.HDFStore(file_name,mode='r')
#            cycle=[]
#            for i in range(0,len(cycle_store)):
#                cycle.append(cycle_store['cycle'+str(i)])
#            cycle_store.close()
#
#            cycle_st_time=[]
#            for i in range(0,len(cycle)):
#                cycle_st_time.append(cycle[i].index[0])
#            
#            return [cycle , cycle_st_time]
        
    else:
#        if file=="h_20.csv":
#            h_20=pd.read_csv(folder+"\\h_20.csv")
#            h_20.index=df.index
#        if file=="h_21.csv":
#            h_21=pd.read_csv(folder+"\\h_21.csv")
#            h_21.index=df.index
#        if file=="h_12.csv":
#            h_12=pd.read_csv(folder+"\\h_12.csv")
#            h_12.index=df.index
        
                
        #-------cycle[cycle #][0=Temp, 1=Pressure]['Instrument Name']
        for i in range(0,len(BV)-1): #loop through all Data
#            H=pd.DataFrame()
            H2=pd.DataFrame()
            
    ###Select One cycle of Data:
            df_state=df[df.index>=dateutil.parser.parse(str(BV.index[i])) ] #select all data with index greater than valve swtiching time 
            df_state=df_state[df_state.index<= dateutil.parser.parse(str(BV.index[i+1])) ] #select all data with index less than the next valve switching time
            
            if i==0:
                df_state=df_state[0:-1]
            else:
                df_state=df_state[1:-1]
            
            if len(df_state) == 0:
                continue #if no data in state 
            else:
                cycle_st_time.append(df_state.index[0])

                for j in range(0,len(df_state)):
#                    h={}
        
#                    state={'state' : BV['Position'][i]}
                    
    #               h,rho,spec_heat,viscosity,conductivity


                    
                    
                    props_list=['h','rho','mu']
#                    props_list_keys=[0,1,3]
                    props_list_keys=[6,2,11]
                    props_list_special=[1,1,10e5]
    
                    props_dict={'20':['TI20','PT03'],
                                '21':['TI21','PT01'],
                                '22':['TI22','PT03'],
                                '23':['TI23','PT01'],
                                '16':['TI16','PT03'],
                                '15':['TI15','PT03'],
                                '14':['TI14','PT03'],
                                '13':['TI13','PT03'],
                                '12':['TI12','PT03'],
                                '11':['TI11','PT01'],
                                '10':['TI10','PT01'],
                                '09':['TI09','PT01'],
                                '08':['TI08','PT01'],
                                '07':['TI07','PT01'],                                
                                }
                    
                    h2={}
    
                    for k in range(0,len(props_list)):
                        for z,tp in props_dict.items():
#                            h2.update({props_list[k]+z : props.f90wrap_enthalpy_tp(df_state[tp[0]][j]+273.15,df_state[tp[1]][j]*6.89475729)[props_list_keys[k]]/props_list_special[k]})
                            h2.update({props_list[k]+z : props5.f90wrap_tp(df_state[tp[0]][j]+273.15,df_state[tp[1]][j]*6.89475729)[props_list_keys[k]]/props_list_special[k]})
                            
#                    h2.update(state)
                    

                    D=0.00396875
                    Area=np.pi*(D/2)**2
        
                    k1_guess=1.1
                    k2_guess=1.1  
                    
                    mu1=h2['mu07']
                    mu2=h2['mu07']
                    rho1=h2['rho07']
                    rho2=h2['rho12']
#                    mu1=mu_07.get('mu07')
#                    mu2=mu_12.get('mu12')
#                    rho1=rho_07.get('rho07')
#                    rho2=rho_12.get('rho12')
                    DP1=abs(df_state['DP03'][j]*6894.75729)
                    DP2=abs(df_state['DP04'][j]*6894.75729)
                    vel1=np.sqrt(2*DP1/(rho1*k1_guess))
                    vel2=np.sqrt(2*DP2/(rho2*k2_guess))
                    Re1=(rho1*vel1*D)/mu1
                    Re2=(rho2*vel2*D)/mu2
                    
                    tol=1e-5
                    err1=1
                    err2=1
                    
                    #calculate mass flow through oriface plates using k-fit equations
                    if df_state['state'][0]=='STATE1':
                        k2=k2_in(Re2)
                        k1=k1_out(Re1)
                        while err2>tol:
                            k2_guess=k2
                            vel2=np.sqrt(2*DP2/(rho2*k2_guess))
                            Re2=(rho2*vel2*D)/(mu2)
                            k2=k2_in(Re2)
                            err2=abs(abs(k2_guess)-abs(k2))
                        
                        while err1>tol:
                            k1_guess=k1
                            vel1=np.sqrt(2*DP1/(rho1*k1_guess))
                            Re1=(rho1*vel1*D)/(mu1)
                            k1=k1_out(Re1)
                            err1=abs(abs(k1_guess)-abs(k1))
                        
                    if df_state['state'][0]=='STATE3':
                        k2=k2_out(Re2)
                        k1=k1_in(Re1)
                        while err2>tol:
                            k2_guess=k2
                            vel2=np.sqrt(2*DP2/(rho2*k2_guess))
                            Re2=(rho2*vel2*D)/(mu2)
                            k2=k2_out(Re2)
                            err2=abs(abs(k2_guess)-abs(k2))
                            
                        while err1>tol:
                            k1_guess=k1
                            vel1=np.sqrt(2*DP1/(rho1*k1_guess))
                            Re1=(rho1*vel1*D)/(mu1)
                            k1=k1_in(Re1)
                            err1=abs(abs(k1_guess)-abs(k1))
                    
                    m1={'mdot1' : rho1*Area*vel1}
                    m2={'mdot2' : rho2*Area*vel2}
                    
                    
                    h2.update(m1)
                    h2.update(m2)


                    H2_j=pd.DataFrame(h2,index=[df_state.index[j]])
                    H2=H2.append(H2_j)
        
                      
                df_state=pd.concat([df_state,H2],axis=1)
                DFF_state=DFF_state.append(df_state)
                


#                cycle.append(df_state)

#        cycle_edited=[]
#    
#        for i in range(0,len(cycle)):
#            if len(cycle[i])!=0:
#                cycle_edited.append(cycle[i]) #creat list of non-zero cycles
                
                
#        if not os.path.exists(folder+"\\cycle"): #check if 'cycle' folder exists in directory 
#            os.makedirs(folder+"\\cycle")
#        for i in range(0,len(cycle_edited)):
#            cycle_edited[i].to_hdf(folder+"\\cycle\\cycle"+str(i)+".h5",'table',mode='w') #write cycle information to 'cycle' folder
        
#        return [cycle_edited, cycle_st_time]
        L_top=5*0.0254
        L_bottom=14*0.0254
        g=9.81
        T=40
        dp_corrected=[]
        for i in range(0,len(DFF_state)):
            Conditions=props_vals(props5.f90wrap_tp(T+273.15,DFF_state['PT01'][i]*6.89475729))
            dp_top_i=(Conditions.density*g*L_top)*0.000145038
            dp_bottom_i=(Conditions.density*g*L_bottom)*0.000145038
            dp_corrected.append(DFF_state['DP01'][i]-dp_top_i-dp_bottom_i)
            
        del DFF_state['DP01']
        DFF_state['DP01']=dp_corrected
        
        DFF_state.to_hdf(folder+"\\df_full_cols.h5",'table')

        return DFF_state