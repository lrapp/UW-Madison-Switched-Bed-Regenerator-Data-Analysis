import FIT_HXR
from scipy import integrate

def ave_for_fit(full,target,valve_time):
#    target=1
    target_code=1
    n_hxrs=25
    
    dp_1=0
    dp_2=dp_1
    
    fluid_1='carbondioxide'
    fluid_2='carbondioxide'
    
    P_convert=6.89475729
    
    S1_time=list(full.loc[full['state']=='STATE1']['rel_time_ave'])[-1]
#    valve_time=1.2
    
    #t_in_1 is T12_ave of state1 
    t_in_1_RE2=full.loc[full['state']=='STATE1']["TI12"+"_ave"].mean()                
    t_in_1_RE2_std=full.loc[full['state']=='STATE1']["TI12"+"_std"].mean()      
   
    #t_in_2 is T16_ave of state3              
    t_in_2_RE2=full.loc[full['state']=='STATE3']["TI16"+"_ave"].mean()                  
    t_in_2_RE2_std=full.loc[full['state']=='STATE3']["TI16"+"_std"].mean()
      
    #p_in_1 is PT03_ave of state1                 
    p_in_1_RE2=full.loc[full['state']=='STATE1']["PT03"+"_ave"].mean()*P_convert              
    p_in_1_RE2_std=full.loc[full['state']=='STATE1']["PT03"+"_std"].mean()*P_convert    
    
    #p_in_2 is PT03_ave of state3
    p_in_2_RE2=full.loc[full['state']=='STATE3']["PT03"+"_ave"].mean()*P_convert       
    p_in_2_RE2_std=full.loc[full['state']=='STATE3']["PT03"+"_std"].mean()*P_convert
 
    m_dot_1_RE2=integrate.trapz((1/(S1_time+valve_time))*full.loc[full['state']=='STATE1']["FI01"+"_ave"],
                                full.loc[full['state']=='STATE1']['rel_time_ave'])

    m_dot_1_RE2_std=full.loc[full['state']=='STATE1']["FI01"+"_std"].mean()
    m_dot_2_RE2=integrate.trapz((1/(S1_time+valve_time))*full.loc[full['state']=='STATE3']["FI01"+"_ave"],
                                full.loc[full['state']=='STATE3']['rel_time_ave'])
      
    m_dot_2_RE2_std=full.loc[full['state']=='STATE3']["FI01"+"_std"].mean()   

    #t_in_1 is T12_ave of state1
    t_in_1_RE1=full.loc[full['state']=='STATE3']["TI07"+"_ave"].mean()                  
    t_in_1_RE1_std=full.loc[full['state']=='STATE3']["TI07"+"_std"].mean() 
    
    #t_in_2 is T16_ave of state3
    t_in_2_RE1=full.loc[full['state']=='STATE1']["TI11"+"_ave"].mean()              
    t_in_2_RE1_std=full.loc[full['state']=='STATE1']["TI11"+"_std"].mean()
                  
    #p_in_1 is PT03_ave of state1
    p_in_1_RE1=full.loc[full['state']=='STATE3']["PT01"+"_ave"].mean()*P_convert                
    p_in_1_RE1_std=full.loc[full['state']=='STATE3']["PT01"+"_std"].mean()*P_convert
    
    #p_in_2 is PT03_ave of state3    
    p_in_2_RE1=full.loc[full['state']=='STATE1']["PT01"+"_ave"].mean()*P_convert           
    p_in_2_RE1_std=full.loc[full['state']=='STATE1']["PT01"+"_std"].mean()*P_convert
    
    m_dot_1_RE1=integrate.trapz((1/(S1_time+valve_time))*full.loc[full['state']=='STATE3']["FI01"+"_ave"],
                                full.loc[full['state']=='STATE3']['rel_time_ave'])
    m_dot_1_RE1_std=full.loc[full['state']=='STATE3']["FI01"+"_std"].mean()     

    m_dot_2_RE1_std=full.loc[full['state']=='STATE1']["FI01"+"_std"].mean()    
    
    m_dot_2_RE1=integrate.trapz((1/(S1_time+valve_time))*full.loc[full['state']=='STATE1']["FI01"+"_ave"],
                                full.loc[full['state']=='STATE1']['rel_time_ave'])


    [epsilon_target_RE1,dt_min_RE1,ua_RE1,q_dot_RE1,q_dot_max_RE1,t_1_RE1,t_2_RE1]=FIT_HXR.counter_flow(fluid_1,fluid_2,target,target_code,n_hxrs,m_dot_1_RE1,m_dot_2_RE1,t_in_1_RE1,t_in_2_RE1,p_in_1_RE1,p_in_2_RE1,dp_1,dp_2)

    [epsilon_target_RE2,dt_min_RE2,ua_RE2,q_dot_RE2,q_dot_max_RE2,t_1_RE2,t_2_RE2]=FIT_HXR.counter_flow(fluid_1,fluid_2,target,target_code,n_hxrs,m_dot_1_RE2,m_dot_2_RE2,t_in_1_RE2,t_in_2_RE2,p_in_1_RE2,p_in_2_RE2,dp_1,dp_2)
   
    
    return [[epsilon_target_RE1,dt_min_RE1,ua_RE1,q_dot_RE1,q_dot_max_RE1,t_1_RE1,t_2_RE1],
            [epsilon_target_RE2,dt_min_RE2,ua_RE2,q_dot_RE2,q_dot_max_RE2,t_1_RE2,t_2_RE2]]
