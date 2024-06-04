"""
Created on 2022.
@author: Louise Xiang

"""

###Librairies : 
import numpy as np 
import pandas as pd


def np_to_pd_3 (f) :
    
    """
    In our case : It is used for all files resulted from the inversion.
    Display data (inside a table) from matlab file to pandas version.
    /!\ 3 entry parameters inside the table : 
                - t : occurenec time of EQ.
                - cumulative EQ : predited EQ at the associated occurence time t.
                - nindex : only two possibilities : 
                               0 if EQ don't belong to the studied area, 
                               otherwise 1.
    
    Input :
        #path : path to the file.
        f : file name.
    """
    
    
    t, lam, ind = np.loadtxt(f, unpack = True)
    data = {"t":t, "lam":lam, "ind":ind}
    df = pd.DataFrame(data)
    
    return df
