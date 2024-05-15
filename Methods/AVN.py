###Librairies : 
import pandas as pd

def AVN(df) :
    """
    Finding AVN i.e Amatrice, Visso & Norcia data inside the table which contains all metadata.
    
    Input :
        df : table which contains all metadata.
    
    """
    ##Sorted df to see the highest magnitude values :
    df_m = df.sort_values('m', ascending=False)

    ##The 3 main EQ : 
    A = df_m.iloc[1]
    V = df_m.iloc[2]
    N = df_m.iloc[0]

    ##Contruct a new table :
    df_AVN = pd.concat([A, V, N], axis=1)
    df_AVN = df_AVN.transpose()
    
    return A, V, N, df_AVN