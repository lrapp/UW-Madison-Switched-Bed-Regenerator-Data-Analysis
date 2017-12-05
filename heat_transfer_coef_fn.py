def ht_coef_fn(start,end,df_full_cols):
    import numpy as np
    import props5
    from ss_props import ss_props
    from ave_for_fit import ave_for_fit
    from calc_Q import calc_Q
    from fast_ave import fast_ave
    from class_def import props,Q,HT_results,main_char
    from scipy import integrate
    import pandas as pd
    
#    def calc_Q(start,end,df_full_cols):
#        def integrate_method(self, how='trapz', unit='s'):
#            '''Numerically integrate the time series.
#        
#            @param how: the method to use (trapz by default)
#            @return 
#        
#            Available methods:
#             * trapz - trapezoidal
#             * cumtrapz - cumulative trapezoidal
#             * simps - Simpson's rule
#             * romb - Romberger's rule
#        
#            See http://docs.scipy.org/doc/scipy/reference/integrate.html for the method details.
#            or the source code
#            https://github.com/scipy/scipy/blob/master/scipy/integrate/quadrature.py
#            '''
#            available_rules = set(['trapz', 'cumtrapz', 'simps', 'romb'])
#            if how in available_rules:
#                rule = integrate.__getattribute__(how)
#            else:
#                print('Unsupported integration rule: %s' % (how))
#                print('Expecting one of these sample-based integration rules: %s' % (str(list(available_rules))))
#                raise AttributeError
#            
#            result = rule(self.values, self.index.astype(np.int64) / 10**9)
#            #result = rule(self.values)
#            return result
#    
#    #    pd.TimeSeries.integrate = integrate_method
#        pd.Series.integrate = integrate_method
    def dp_fs(porosity,d_p,Conditions):
        Regen_dia=0.0492506 #1.94 inches
        A_fr=np.pi*(Regen_dia/2)**2 #frontal area 
        L=0.445
        vel=Conditions.mdot/(Conditions.density*A_fr)
        mu=Conditions.visc/1000000
        Re=(d_p*Conditions.density*vel)/mu
        
        Re_m=Re/(1-porosity)
        f2=(1.87*porosity**0.75)/(1-porosity)**0.26
        f1L=136/(1-porosity)**0.38
        f1T=29/((1-porosity)**1.45*porosity**2)
        q=np.exp(-(porosity**2*(1-porosity))*Re_m/12.6)
        fp_fs=(q*(f1L/Re_m)+(1-q)*(f2+f1T/Re_m))*((1-porosity)/porosity**3)
        

        dP_fs=((fp_fs*Conditions.density*vel**2*L)/d_p)*0.000145038
        
        return dP_fs
    
    def dp_KTA(porosity,d_p,Conditions):
        Regen_dia=0.0492506 #1.94 inches
        A_fr=np.pi*(Regen_dia/2)**2 #frontal area 
        L=0.445
        vel=Conditions.mdot/(Conditions.density*A_fr)
        mu=Conditions.visc/1000000
        Re=(d_p*Conditions.density*vel)/mu
        
        Re_m=Re/(1-porosity)
        fp=(160+3*Re_m**0.9)*((1-porosity)**2/(porosity**3*Re))
        

        dP_KTA=((fp*Conditions.density*vel**2*L)/d_p)*0.000145038
        
        return dP_KTA


    def regen_info(porosity,d_p,rho_matrix):
#        L=0.43815 #Length of Regenerator
        L=0.445
        d_p=3.175/1000 #sphere diameter
#        Regen_dia=0.049276 #1.94 inches
        Regen_dia=0.0492506
        A_fr=np.pi*(Regen_dia/2)**2 #frontal area 
    #    r_char=(porosity*d_p)/(4*(1-porosity))
    #    d_p=3/1000 #sphere diameter
        Regenerator_total_volume=A_fr*L #m^3
        void_volume=Regenerator_total_volume*porosity #m^3
        matrix_volume=Regenerator_total_volume-void_volume #m^3
    
        mass_matrix=matrix_volume*rho_matrix
#        A_ht=((4*(1-porosity))/d_p)*Regenerator_total_volume
        A_ht=6*(1-porosity)*(A_fr*L)/d_p
        
        return [mass_matrix,A_ht]
        
          
    def heat_transfer_coef(porosity,d_p,Conditions):
        mu=Conditions.visc/1000000
        cp=Conditions.cp
        cond=Conditions.cond
        m_dot=Conditions.mdot
        
        Regen_dia=0.0492506 #1.94 inches
        A_fr=np.pi*(Regen_dia/2)**2 #frontal area 
        
#        r_char=(porosity*d_p)/(4*(1-porosity))
        r_char=(porosity*d_p)/(6*(1-porosity))
       
        G=m_dot/(porosity*A_fr)
    
        Re=(4*G*r_char)/mu
        
        j_h=0.23*Re**(-0.3)
        
        Pr=(mu*(cp*1000))/cond
        
        h_bar=(j_h*G*cp)/Pr**(2/3) #kW/m^2-K
    
        return h_bar, Re,Pr
    
    def ht_coef_Achenbach(porosity,d_p,Conditions):
        mu=Conditions.visc/1000000
        cp=Conditions.cp
        cond=Conditions.cond
        m_dot=Conditions.mdot
        rho=Conditions.density
        
        Regen_dia=0.0492506 #1.94 inches
        A_fr=np.pi*(Regen_dia/2)**2 #frontal area 
        
        
        vel=m_dot/(rho*A_fr)
        
        Re=(rho*vel*d_p)/mu
        
        Re_m=Re/(1-porosity)
        
        Nu=((1.18*Re**(0.58))**4+((0.23*Re_m**0.75)**4))**(1/4)
        
        h_bar=(Nu*cond)/d_p
    
        return h_bar, Re_m
    
    def ht_coef_Gn(porosity,d_p,Conditions):
        mu=Conditions.visc/1000000
        cp=Conditions.cp
        cond=Conditions.cond
        m_dot=Conditions.mdot
        rho=Conditions.density
        
        Regen_dia=0.0492506 #1.94 inches
        A_fr=np.pi*(Regen_dia/2)**2 #frontal area 
        
        Vel=m_dot/(rho*A_fr)
        
        Re=(rho*Vel*d_p)/mu
        
#        Re_phi=(Vel*d_p)/((mu/rho)*porosity)
        Re_phi=Re/porosity
        
        fa=1+1.5*(1-porosity)
        
        Pr=(mu*(cp*1000))/cond
        
        Nu_lam=0.664*np.sqrt(Re_phi)*(Pr)**(1/3)
        
        Nu_turb=(0.037*Re_phi**(0.8)*Pr)/((1+2.443*Re_phi**(-0.1))*(Pr**(2/3)-1))
        
        Nu_sphere=2+np.sqrt(Nu_lam**2+Nu_turb**2)
        
        Nu=fa*Nu_sphere
        
        h_bar=(Nu*cond)/d_p
    
        return h_bar, Re_phi
    
    
        
    #start=100
    #end=150
    #    [S1_ave,S2a_ave,S2b_ave,S2c_ave,S3_ave,S4a_ave,S4b_ave,S4c_ave,full_cycle]=cycle_ave(folder,new_cycle,start,end,df)
    
    ##Get Averaged States
    [S1_ave,S2a_ave,S2b_ave,S2c_ave,S3_ave,S4a_ave,S4b_ave,S4c_ave,full_cycle]=fast_ave(start,end,df_full_cols)
    
    ##Get State Switching Times
    S1_swt=list(S1_ave['rel_time_ave'])[-1]
    S2a_swt=list(S2a_ave['rel_time_ave'])[-1]
    S2b_swt=list(S2b_ave['rel_time_ave'])[-1]    
    S3_swt=list(S3_ave['rel_time_ave'])[-1]
    S4a_swt=list(S4a_ave['rel_time_ave'])[-1]
    S4b_swt=list(S4b_ave['rel_time_ave'])[-1]
        
    ##Set Averaged Variables for calculations
    T_C_in_RE1=full_cycle.loc[full_cycle['state']=='STATE3']["TI07"+"_ave"].mean()  
    T_C_in_RE1_std=full_cycle.loc[full_cycle['state']=='STATE3']["TI07"+"_std"].mean() 
    
    P_C_RE1=full_cycle.loc[full_cycle['state']=='STATE3']["PT01"+"_ave"].mean()*6.89475729
    P_C_RE1_std=full_cycle.loc[full_cycle['state']=='STATE3']["PT01"+"_std"].mean()*6.89475729
    
    m_dot_C_RE1=integrate.trapz((1/(list(S1_ave['rel_time_ave'])[-1]))*full_cycle.loc[full_cycle['state']=='STATE3']["FI01"+"_ave"],S3_ave['rel_time_ave'])
    m_dot_C_RE1_EES=integrate.trapz((1/(list(S1_ave['rel_time_ave'])[-1]+S2a_swt+S2b_swt))*full_cycle.loc[full_cycle['state']=='STATE3']["FI01"+"_ave"],S3_ave['rel_time_ave'])
    m_dot_C_RE1_std=full_cycle.loc[full_cycle['state']=='STATE3']["FI01"+"_std"].mean()     
    
    T_H_in_RE1=full_cycle.loc[full_cycle['state']=='STATE1']["TI11"+"_ave"].mean()
    T_H_in_RE1_std=full_cycle.loc[full_cycle['state']=='STATE1']["TI11"+"_std"].mean()
    
    P_H_RE1=full_cycle.loc[full_cycle['state']=='STATE1']["PT01"+"_ave"].mean()*6.89475729
    P_H_RE1_std=full_cycle.loc[full_cycle['state']=='STATE1']["PT01"+"_std"].mean()*6.89475729
    
    m_dot_H_RE1=integrate.trapz((1/(list(S1_ave['rel_time_ave'])[-1]))*full_cycle.loc[full_cycle['state']=='STATE1']["FI01"+"_ave"],S1_ave['rel_time_ave'])
    m_dot_H_RE1_EES=integrate.trapz((1/(list(S1_ave['rel_time_ave'])[-1]+S4a_swt+S4b_swt))*full_cycle.loc[full_cycle['state']=='STATE1']["FI01"+"_ave"],S1_ave['rel_time_ave'])
    m_dot_H_RE1_std=full_cycle.loc[full_cycle['state']=='STATE1']["FI01"+"_std"].mean()    
    
    
    Max_Results_New=HT_results(ave_for_fit(full_cycle,1,1.2))
    Max_Results_Original=HT_results(ave_for_fit(full_cycle,1,0))
    
    valve_time=0
    
    Q_re=Q(calc_Q(start,end,df_full_cols))
    
    eff_C1=(Q_re.Q_C1_ave/(Q_re.half_cycle_time_ave+valve_time))/Max_Results_Original.q_dot_max_RE1
    eff_H1=(Q_re.Q_H1_ave/(Q_re.half_cycle_time_ave+valve_time))/Max_Results_Original.q_dot_max_RE1    

    eff_C1_new=Q_re.Q_C1_ave/(Max_Results_New.q_dot_max_RE1*(S1_swt+S2a_swt+S2b_swt))
    eff_H1_new=Q_re.Q_H1_ave/(Max_Results_New.q_dot_max_RE1*(S3_swt+S4a_swt+S4b_swt))    
    
    
    Actual_Results=HT_results(ave_for_fit(full_cycle,eff_C1,0))
    
    T_C_out_RE1=Actual_Results.t_1_RE1[0]
    T_H_out_RE1=Actual_Results.t_2_RE1[-1]
        
    T_C_ave=(T_C_in_RE1+T_C_out_RE1)/2
    T_H_ave=(T_H_in_RE1+T_H_out_RE1)/2
    
    #Set Cold side properties
    C_in_conditions=props(props5.f90wrap_tp(T_C_in_RE1+273.15,P_C_RE1))
    C_in_conditions.mdot=m_dot_C_RE1_EES
    
    #Set Hot side properties
    H_in_conditions=props(props5.f90wrap_tp(T_H_in_RE1+273.15,P_H_RE1))
    H_in_conditions.mdot=m_dot_H_RE1_EES
    
    #Get Cold Averaged Conditions
    C_ave_conditions=props(props5.f90wrap_tp(T_C_ave+273.15,P_C_RE1))
    C_ave_conditions.mdot=m_dot_C_RE1_EES
    
    #Get Hot Averaged Conditions
    H_ave_conditions=props(props5.f90wrap_tp(T_H_ave+273.15,P_H_RE1))
    H_ave_conditions.mdot=m_dot_H_RE1_EES
    
    #Calculate maximum enthalpies 
    h_H_out_max=props5.f90wrap_tp(T_C_in_RE1+273.15,P_H_RE1)[6] #Cold inlet with Hot side Pressure
    h_C_out_max=props5.f90wrap_tp(T_H_in_RE1+273.15,P_C_RE1)[6] #Hot inlet with Cold side Pressure
    
    #Calculate Specific Heats                                 
    C_p_H=(H_in_conditions.enth-h_H_out_max)/(H_in_conditions.temp-C_in_conditions.temp)
    C_p_C=(h_C_out_max-C_in_conditions.enth)/(H_in_conditions.temp-C_in_conditions.temp)
    
    
    C_dot_C=C_in_conditions.mdot*C_p_C
    C_dot_H=H_in_conditions.mdot*C_p_H
    
    C_dot_min=min(C_dot_C,C_dot_H)
    C_dot_max=max(C_dot_C,C_dot_H)
    
    C_R=C_dot_min/C_dot_max
    
    if len(S2c_ave)==0:
        switching_time=(S1_ave['rel_time_ave'].iloc[-1]+S2a_ave['rel_time_ave'].iloc[-1]+S2b_ave['rel_time_ave'].iloc[-1]+S3_ave['rel_time_ave'].iloc[-1]+S4a_ave['rel_time_ave'].iloc[-1]+S4b_ave['rel_time_ave'].iloc[-1])
    else:    
        switching_time=(S1_ave['rel_time_ave'].iloc[-1]+S2a_ave['rel_time_ave'].iloc[-1]+S2b_ave['rel_time_ave'].iloc[-1]+S2c_ave['rel_time_ave'].iloc[-1]+S3_ave['rel_time_ave'].iloc[-1]+S4a_ave['rel_time_ave'].iloc[-1]+S4b_ave['rel_time_ave'].iloc[-1]+S4c_ave['rel_time_ave'].iloc[-1])
    
    #Calculate matrix properties
    c_matrix_C=ss_props(C_in_conditions.temp)[0]
    rho_matrix_C=ss_props(C_in_conditions.temp)[1]
    
    c_matrix_H=ss_props(H_in_conditions.temp)[0]
    rho_matrix_H=ss_props(H_in_conditions.temp)[1]
    
    #average cold and hot properties
    c_matrix=np.average([c_matrix_C,c_matrix_H])
    rho_matrix=np.average([rho_matrix_C,rho_matrix_H])
    
    porosity=0.37
    d_p=3.175/1000 #1/8 inches = 3.175 mm
    
    mass_matrix,A_ht=regen_info(0.37,d_p,rho_matrix)
    
    
    h_bar_c,Re_c,Pr_c=heat_transfer_coef(porosity,d_p,C_ave_conditions) #kW/m^2-K
    h_bar_h,Re_h,Pr_h=heat_transfer_coef(porosity,d_p,H_ave_conditions)
   
#    L_top=5*0.0254
#    L_bottom=14*0.0254
#    g=9.81
#    HTCB_correction_conditions=props(props5.f90wrap_tp(40+273.15,H_ave_conditions.pressure))
#    HTCB_dp_top=(HTCB_correction_conditions.density*g*L_top)*0.000145038
#    HTCB_dp_bottom=(HTCB_correction_conditions.density*g*L_bottom)*0.000145038
#                   
#    CTHB_correction_conditions=props(props5.f90wrap_tp(40+273.15,C_ave_conditions.pressure))
#    CTHB_dp_top=(CTHB_correction_conditions.density*g*L_top)*0.000145038
#    CTHB_dp_bottom=(CTHB_correction_conditions.density*g*L_bottom)*0.000145038
    
    dP_HTCB_fs=dp_fs(porosity,d_p,H_ave_conditions)
    dP_CTHB_fs=dp_fs(porosity,d_p,C_ave_conditions)
    
    dP_HTCB_KTA=dp_KTA(porosity,d_p,H_ave_conditions)
    dP_CTHB_KTA=dp_KTA(porosity,d_p,C_ave_conditions)
    
    h_bar_c_achen,Re_c_achen=ht_coef_Achenbach(porosity,d_p,C_ave_conditions)
    
    h_bar_c_GN,Re_phi=ht_coef_Gn(porosity,d_p,C_ave_conditions)


    C_m=(2*mass_matrix*c_matrix*(1/switching_time))/C_dot_min
    C_m_e=(2*C_R*C_m)/(1+C_R)
    
    NTU_R=1/(C_dot_min)*((1/(h_bar_c)*A_ht)+(1/(h_bar_h)*A_ht))**(-1)
    
#    NTU_C=(h_bar_c*A_ht)/(C_conditions.mdot*C_conditions.cp)
#    NTU_H=(h_bar_h*A_ht)/(H_conditions.mdot*H_conditions.cp)
    
#    NTU_e_C=(2*C_R*NTU_C)/(1+C_R)
#    NTU_e_H=(2*C_R*NTU_H)/(1+C_R)

    
#    NTU=np.mean([NTU_C,NTU_H])
    
    NTU_R_e=2*C_R*NTU_R/(1+C_R)
    
    Re=np.mean([Re_c,Re_h])
    
    MC=main_char()
    MC.C_m_e=C_m_e
    MC.NTU_e=NTU_R_e
    MC.Re=Re
    MC.switching_time=switching_time
    MC.eff_C1=eff_C1
    MC.eff_H1=eff_H1
    MC.eff_C1_new=eff_C1_new
    MC.eff_H1_new=eff_H1_new
    MC.Q_C1_ave=Q_re.Q_C1_ave
    MC.Q_H1_ave=Q_re.Q_H1_ave
    MC.Q_C1_std=Q_re.Q_C1_std
    MC.Q_H1_std=Q_re.Q_H1_std
    MC.half_cycle_time_ave=Q_re.half_cycle_time_ave
    MC.S1_swt=S1_swt
    MC.S2a_swt=S2a_swt
    MC.S2b_swt=S2b_swt
    MC.S3_swt=S3_swt
    MC.S4a_swt=S4a_swt
    MC.S4b_swt=S4b_swt
        
    
    MC.q_dot_max_RE1=Max_Results_New.q_dot_max_RE1
    
    MC.T_H_RE1=T_H_in_RE1
    MC.T_C_RE1=T_C_in_RE1
    MC.T_H_RE1_std=T_H_in_RE1_std
    MC.T_C_RE1_std=T_C_in_RE1_std
    
    MC.P_H_RE1=P_H_RE1
    MC.P_C_RE1=P_C_RE1
    MC.P_H_RE1_std=P_H_RE1_std
    MC.P_C_RE1_std=P_C_RE1_std
    
    MC.mdot_H_RE1=m_dot_H_RE1
    MC.mdot_C_RE1=m_dot_C_RE1
    MC.mdot_H_RE1_std=m_dot_H_RE1_std
    MC.mdot_C_RE1_std=m_dot_C_RE1_std
    
    MC.mdot_H_RE1_EES=m_dot_H_RE1_EES
    MC.mdot_C_RE1_EES=m_dot_C_RE1_EES
    
    MC.dP_HTCB_fs=dP_HTCB_fs
    MC.dP_CTHB_fs=dP_CTHB_fs
    
    MC.dP_HTCB_KTA=dP_HTCB_KTA
    MC.dP_CTHB_KTA=dP_CTHB_KTA
    
    MC.Pr_c=Pr_c
    MC.Pr_h=Pr_h
    
#    return C_m_e,NTU_R_e,switching_time,Re
    return MC
