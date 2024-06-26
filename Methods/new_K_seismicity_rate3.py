"""
Created on 2022.
@author: Louise Xiang

"""

###Librairies : 
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd


def seismicity_rate2(dn, df, dff) :
    
    """
    Compute seismicity rate of, can be :
        - predicted seismicity with EQ detection problem. --> df1.lam
        - predicted seismicity without EQ detection problem. --> df1.lam_tot
        - true seismicity with EQ detection problem. --> df.index
    
    Input : 
        dn : slide window in order to smooth our results 
        df : observed catalog (contains time data)
        dff : synthetic catalog (contain cumulative seismic events)
    """
    
    EQ_rate = np.zeros(len(df))

    for n in range(0+int(dn/2), len(df)-dn) :
        if (n >= (dn/2)) & (n <= (len(df) - (dn/2))) :
            n1 = int( n - (dn/2) )
            n2 = int( n + (dn/2) )
#             print("n =", n, "n1 =", n1, "n2 =", n2)

            df_t =  df.t.values[n1:n2]
            df_lam= dff.values[n1:n2] 

            duration = df_t[-1] - df_t[0]
            tot_seismes = df_lam[-1] - df_lam[0]

            taux_x = tot_seismes/duration
            EQ_rate[n] = taux_x

    return EQ_rate




























