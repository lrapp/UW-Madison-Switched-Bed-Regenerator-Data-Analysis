#ss density 0 < T < 1200 K

def ss_props(T):
    import numpy as np
    
    if (T<0 or T > 1200):
        print("T is below 0 or greater than 1200 K")
        return [np.nan, np.nan]
        
    
    a_all=[-1.4087,12.2486,-1879.464]
    b_all=[1.3982,-80.6422,3643.198]
    c_all=[0.2543,218.743,76.70125]
    d_all=[-0.6260,-308.854,-6176.028]
    e_all=[0.2334,239.5296,7437.6247]
    f_all=[0.4256,-89.9982,-4305.7217]
    g_all=[-0.4658,3.15315,1382.4627]
    h_all=[0.1650,8.44996,-237.22704]
    i_all=[-0.0199,-1.91368,17.05262]

#    x=2
#    c=10**(a_all[x])+b_all[x]*np.log10(T)+c_all[x]*(np.log10(T))**2+d_all[x]*(np.log10(T))**3+e_all[x]*(np.log10(T))**4+f_all[x]*(np.log10(T))**5+g_all[x]*(np.log10(T))**6+h_all[x]*(np.log10(T))**7+i_all[x]*(np.log10(T))**8
#    c=0.433294+0.000176199*T+7.68038e-8*T**2-8.83519e-11*T**3
    
    c=0.354425228 + 0.000672263976*T - 0.00000101732303*T**2 + 9.22897080E-10*T**3 - 3.33296910E-13*T**4
    
    rho=8044.05-0.422031*(T-273.15)
    
    return [c,rho]