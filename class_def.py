import time as time
class Timer(object):
    
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()
    def __exit__(self, type, value, traceback):
        if self.name:
            print('[%s]' % self.name,)
        print('Elapsed: %s' % (time.time() - self.tstart))

class props:
    import numpy as np
    def __init__(self, properties = [0]*13, mdot = np.NaN):
        self.temp  =  properties[0]
        self.pressure  =  properties[1]
        self.density  =  properties[2]
        self.vol  =  properties[3]
        self.qual = properties[4]
        self.inte = properties[5]
        self.enth = properties[6]
        self.entr = properties[7]
        self.cv = properties[8]
        self.cp = properties[9]
        self.ssnd = properties[10]
        self.visc = properties[11]
        self.cond = properties[12]
        self.mdot = mdot
    def __repr__(self):
        units = ['K','kPa','kg/m3','m3/kg','','kJ/kg','kJ/kg','kJ/kg-K','kJ/kg-K','kJ/kg-K','m/s','uPa-s','W/m-K','kg/s']
        return ('temp'+'\t\t'+str(round(self.temp,2))+"\t"+units[0]+"\n"
                +'pressure'+'\t'+str(round(self.pressure,2))+"\t"+units[1]+"\n"
                +'density'+'\t\t'+str(round(self.density))+"\t\t"+units[2]+"\n"
                +'vol'+'\t\t\t'+str(round(self.vol))+"\t\t"+units[3]+"\n"     
                +'qual'+'\t\t'+str(self.qual)+"\t"+units[4]+"\n"  
                +'inte'+'\t\t'+str(round(self.inte,2))+"\t"+units[5]+"\n"                    
                +'enth'+'\t\t'+str(round(self.enth,2))+"\t"+units[6]+"\n"
                +'entr'+'\t\t'+str(round(self.cv,2))+"\t"+units[7]+"\n"   
                +'cv'+'\t\t\t'+str(round(self.cv,2))+"\t"+units[8]+"\n"   
                +'cp'+'\t\t\t'+str(round(self.cp,2))+"\t"+units[9]+"\n"                    
                +'ssnd'+'\t\t'+str(round(self.ssnd,2))+"\t"+units[10]+"\n"
                +'visc'+'\t\t'+str(round(self.visc,2))+"\t"+units[11]+"\n"
                +'cond'+'\t\t'+str(round(self.cond,2))+"\t"+units[12]+"\n"       
                +'mdot'+'\t\t'+str(round(self.mdot,2))+"\t\t"+units[13])
        
class Q:
    def __init__(self,results=[0]*11):
        self.Q_C1_ave=results[0]
        self.Q_H1_ave=results[1]
        self.Q_C2_ave=results[2]
        self.Q_H2_ave=results[3]
        self.Q_C1_std=results[4]
        self.Q_H1_std=results[5]
        self.Q_C2_std=results[6]
        self.Q_H2_std=results[7]
        self.half_cycle_time_ave=results[8]
        self.Q_C1_st2a=results[9]
        self.Q_C1_st2b=results[10]
        
    def __repr__(self):
        return ('Q_C1'+'\t'+str(round(self.Q_C1_ave,2))+" +/- "+str(round(self.Q_C1_std,2))+"\n"
               +'Q_H1'+'\t'+str(round(self.Q_H1_ave,2))+" +/- "+str(round(self.Q_H1_std,2))+"\n"
               +'Q_C2'+'\t'+str(round(self.Q_C2_ave,2))+" +/- "+str(round(self.Q_C2_std,2))+"\n"
               +'Q_H2'+'\t'+str(round(self.Q_H2_ave,2))+" +/- "+str(round(self.Q_H2_std,2))+"\n"
               +'half cycle time= '+str(round(self.half_cycle_time_ave,2)))

class HT_results:
    def __init__(self,results=[[0]*7,[0]*7]):
        self.epsilon_target_RE1=results[0][0]
        self.dt_min_RE1=results[0][1]
        self.ua_RE1=results[0][2]
        self.q_dot_RE1=results[0][3]
        self.q_dot_max_RE1=results[0][4]
        self.t_1_RE1=results[0][5]
        self.t_2_RE1=results[0][6]
        self.q_dot_max_RE2=results[1][4]
    def __repr__(self):
        return ('epsilon_target_RE1'+'\t'+str(round(self.epsilon_target_RE1,4))+"\n"
               +'dt_min_RE1'+'\t'+str(round(self.dt_min_RE1,2))+"\n"
               +'ua_RE1'+'\t'+str(round(self.ua_RE1,2))+"\n"
               +'q_dot_RE1'+'\t'+str(round(self.q_dot_RE1,4))+"\n"
               +'q_dot_max_RE1'+ '\t'+str(round(self.q_dot_max_RE1,4)))

class characteristics:
    def __init__(self,results=[0]*4):
        self.C_m_e=results[0]
        self.NTU_e=results[1]
        self.Re=results[3]
        self.switching_time=results[2]
    def __repr__(self):
        return ('C_m_e'+'\t'+str(round(self.C_m_e,2))+" +/- "+str(round(self.Q_C1_std,2))+"\n"
               +'NTU_e'+'\t'+str(round(self.Q_H1_ave,2))+" +/- "+str(round(self.Q_H1_std,2))+"\n"
               +'Re'+'\t'+str(round(self.Q_C2_ave,2))+" +/- "+str(round(self.Q_C2_std,2))+"\n"
               +'switching_time'+'\t'+str(round(self.Q_H2_ave,2))+" +/- "+str(round(self.Q_H2_std,2))+"\n")       
        
        
class main_char:
    def __init__(self,results=[0]*45):
        self.C_m_e=results[0]
        self.NTU_e=results[1]
        self.Re=results[2]
        self.switching_time=results[3]
        self.eff_C1=results[4]
        self.eff_C2=results[5]
        self.Q_C1_ave=results[6]
        self.Q_H1_ave=results[7]
        self.Q_C2_ave=results[8]
        self.Q_H2_ave=results[9]
        self.Q_C1_std=results[10]
        self.Q_H1_std=results[11]
        self.Q_C2_std=results[12]
        self.Q_H2_std=results[13]
        self.half_cycle_time_ave=results[14]
        self.T_H_RE1=results[15]
        self.T_C_RE1=results[16]
        self.T_H_RE1_std=results[17]
        self.T_H_RE1_std=results[18]
        self.P_H_RE1=results[19]
        self.P_C_RE1=results[20]
        self.P_H_RE1_std=results[21]
        self.P_H_RE1_std=results[22]
        self.mdot_H_RE1=results[23]
        self.mdot_C_RE1=results[24]
        self.mdot_H_RE1_std=results[25]
        self.mdot_C_RE1_std=results[26]
        self.q_dot_max_RE1=results[27]
        self.eff_H1=results[28]
        self.eff_C1_new=results[29]
        self.eff_H1_new=results[30]
        self.mdot_C_RE1_EES=results[31]
        self.mdot_H_RE1_EES=results[32]
        self.S1_swt=results[33]
        self.S2a_swt=results[34]
        self.S2b_swt=results[35]
        self.S3_swt=results[36]
        self.S4a_swt=results[37]
        self.S4b_swt=results[38]
        self.dP_HTCB_fs=results[39]
        self.dP_CTHB_fs=results[40]
        self.dP_HTCB_KTA=results[41]
        self.dP_CTHB_KTA=results[42]
        self.Pr_c=results[43]
        self.Pr_h=results[44]

    def __repr__(self):
        return ('C_m_e'+'\t'+str(round(self.C_m_e,2))+"\n"
               +'NTU_e'+'\t'+str(round(self.NTU_e,2))+"\n"
               +'Re'+'\t'+str(round(self.Re,2))+"\n"
               +"eff_C1"+'\t'+str(round(self.eff_C1,3))+'\n'
               +"eff_H1"+'\t'+str(round(self.eff_H1,3))+'\n'
               +"eff_C1_new"+"\t"+str(round(self.eff_C1_new,3))+'\n'
               +"eff_H1_new"+"\t"+str(round(self.eff_H1_new,3))+'\n'
               +'switching_time'+'\t'+str(round(self.switching_time,2))+"\n"
               +'Q_C1_ave'+str(round(self.Q_C1_ave,2))+" +/- "+str(round(self.Q_C1_std,2))+"\n"
               +'Q_H1_ave'+str(round(self.Q_H1_ave,2))+" +/- "+str(round(self.Q_H1_std,2))+"\n")       
                
               
