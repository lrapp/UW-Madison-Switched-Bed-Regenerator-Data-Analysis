import numpy as np
import pandas as pd
from scipy import integrate

def calc_Q(start,end,df_full_cols):

    
    def integrate_method(self, how='trapz', unit='s'):
        '''Numerically integrate the time series.
    
        @param how: the method to use (trapz by default)
        @return 
    
        Available methods:
         * trapz - trapezoidal
         * cumtrapz - cumulative trapezoidal
         * simps - Simpson's rule
         * romb - Romberger's rule
    
        See http://docs.scipy.org/doc/scipy/reference/integrate.html for the method details.
        or the source code
        https://github.com/scipy/scipy/blob/master/scipy/integrate/quadrature.py
        '''
        available_rules = set(['trapz', 'cumtrapz', 'simps', 'romb'])
        if how in available_rules:
            rule = integrate.__getattribute__(how)
        else:
            print('Unsupported integration rule: %s' % (how))
            print('Expecting one of these sample-based integration rules: %s' % (str(list(available_rules))))
            raise AttributeError
        
        result = rule(self.values, self.index.astype(np.int64) / 10**9)
        #result = rule(self.values)
        return result
    
#    pd.TimeSeries.integrate = integrate_method
    pd.Series.integrate = integrate_method
    
    
    Q_C1=[]
    Q_H1=[]
    Q_C2=[]
    Q_H2=[]
    Q_H2_ST2a=[]
    Q_C1_ST2a=[]
    Q_H2_ST2b=[]
    Q_C1_ST2b=[]    
    
    Q_H2_ST4a=[]
    Q_C1_ST4a=[]
    Q_H2_ST4b=[]
    Q_C1_ST4b=[]       
    
    Q16=0
    Q12=0
    Q11=0
    Q07=0
    
    possible_states=['STATE1','STATE2a','STATE2b','STATE2c',
                     'STATE3','STATE4a','STATE4b','STATE4c']
    
    possible_states=['STATE1','STATE3','STATE2a']
    Q_combos={'STATE1' : [[Q_H2,Q16,Q12],[Q_C1,Q11,Q07]],
              'STATE3' : [[Q_H1,Q11,Q07],[Q_C2,Q16,Q12]],
              'STATE2a' : [[Q_H2_ST2a,Q16,Q12],[Q_C1_ST2a,Q11,Q07]],
              'STATE2b' : [[Q_H2_ST2b,Q16,Q12],[Q_C1_ST2b,Q11,Q07]],
              'STATE4a' : [[Q_H2_ST4a,Q16,Q12],[Q_C1_ST4a,Q11,Q07]],
              'STATE4b' : [[Q_H2_ST4b,Q16,Q12],[Q_C1_ST4b,Q11,Q07]]}
    
    half_cycle_time=[]
    for i in range(start,end):
        df_i=df_full_cols.loc[df_full_cols['unique_state_num']==i]
        Q16=(integrate.trapz(df_i['h16']*df_i['FI01'],df_i.index.astype(np.int64)/10**9))
        Q12=(integrate.trapz(df_i['h12']*df_i['FI01'],df_i.index.astype(np.int64)/10**9))
        Q11=(integrate.trapz(df_i['h11']*df_i['FI01'],df_i.index.astype(np.int64)/10**9))
        Q07=(integrate.trapz(df_i['h07']*df_i['FI01'],df_i.index.astype(np.int64)/10**9))    
        Q_combos={'STATE1' : [[Q_H2,Q16,Q12],[Q_C1,Q11,Q07]],
                  'STATE3' : [[Q_H1,Q11,Q07],[Q_C2,Q16,Q12]],
                  'STATE2a' : [[Q_H2_ST2a,Q16,Q12],[Q_C1_ST2a,Q11,Q07]],
                  'STATE2b' : [[Q_H2_ST2b,Q16,Q12],[Q_C1_ST2b,Q11,Q07]],
                  'STATE4a' : [[Q_H2_ST4a,Q16,Q12],[Q_C1_ST4a,Q11,Q07]],
                  'STATE4b' : [[Q_H2_ST4b,Q16,Q12],[Q_C1_ST4b,Q11,Q07]]}
        for k in range(0,len(possible_states)):
            if df_full_cols.loc[df_full_cols['unique_state_num']==i]['state'][0] == possible_states[k]:
                Q_combos[possible_states[k]][0][0].append(Q_combos[possible_states[k]][0][1]-Q_combos[possible_states[k]][0][2])
                Q_combos[possible_states[k]][1][0].append(Q_combos[possible_states[k]][1][1]-Q_combos[possible_states[k]][1][2])
            if df_full_cols.loc[df_full_cols['unique_state_num']==i]['state'][0] == 'STATE1':
                half_cycle_time.append((df_i.index[-1]-df_i.index[0]).total_seconds())        
    Q_C1_ave=0
    Q_H1_ave=0
    Q_H2_ave=0
    Q_C2_ave=0
    Q_C1_ST2a=[0]
    Q_C1_ST2b=[0]

    if len(Q_C1) != 0: Q_C1_ave=np.mean(Q_C1)
    if len(Q_H1) != 0: Q_H1_ave=np.mean(Q_H1)
    if len(Q_C2) !=0: Q_C2_ave=np.mean(Q_C2)
    if len(Q_H2) !=0: Q_H2_ave=np.mean(Q_H2)
    half_cycle_time_ave=np.mean(half_cycle_time)
    if len(Q_C1) !=0: Q_C1_std=np.std(Q_C1)
    if len(Q_H1) !=0: Q_H1_std=np.std(Q_H1)
    if len(Q_C2) !=0: Q_C2_std=np.std(Q_C2)
    if len(Q_H2) !=0: Q_H2_std=np.std(Q_H2)
    
    if len(Q_C1_ST2a) !=0: Q_C1_st2a_ave=np.mean(Q_C1_ST2a)
    if len(Q_C1_ST2b) !=0: Q_C1_st2b_ave=np.mean(Q_C1_ST2b)

    return [Q_C1_ave,Q_H1_ave,
            Q_C2_ave,Q_H2_ave,
            Q_C1_std,Q_H1_std,
            Q_C2_std,Q_H2_std,
            half_cycle_time_ave,
            Q_C1_st2a_ave,Q_C1_st2b_ave]
