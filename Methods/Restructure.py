###Librairies : 
import numpy as np 
import pandas as pd

###Restucture of data files : 
def np_to_pd_6 (path, f) :
    """
    In our case : It is used for data.dat which contains all metadata.
    Display data (inside a table) from matlab file to pandas version.
    /!\ 6 entry parameters inside the table : 
                - t (days) : time
                - m (ml preference) : magnitude
                - x (° preference) : latitude
                - y (° preference) : longitude
                - z (km) : depth
                - pi : EQ detection probability
    
    Input :
        path : path to the file.
        f : file name.
    """
    
    t, m, x, y, z, pi = np.loadtxt(path+f, unpack = True)
    data = {"t":t, "m":m, "x":x, "y":y, "z":z, "pi":pi}
    df = pd.DataFrame(data)
    
    return df

def np_to_pd_5 (path, f) :
    """
    In our case : It is used for data.dat which contains all metadata.
    Display data (inside a table) from matlab file to pandas version.
    /!\ 5 entry parameters inside the table : 
                - t (days) : time
                - m (ml preference) : magnitude
                - x (° preference) : latitude
                - y (° preference) : longitude
                - z (km) : depth
    
    Input :
        path : path to the file.
        f : file name.
    """
    
    t, m, x, y, z = np.loadtxt(path+f, unpack = True)
    data = {"t":t, "m":m, "x":x, "y":y, "z":z}
    df = pd.DataFrame(data)
    
    return df

def np_to_pd_4 (path, f) :
    
    """
    In our case : It is used for all files resulted from the inversion.
    Display data (inside a table) from matlab file to pandas version.
    /!\ 4 entry parameters inside the table : 
                - lambda : predicted seismicity without EQ detection probabilities correction
                - lambda_tot : prediected seismicity with EQ detection probabilities correction
                - K : productivity prefactors include K and Kj.
                - nrep : number of aftershocks.
    
    Input :
        path : path to the file.
        f : file name.
    """
    
    
    lam, lam_tot, K, nrep = np.loadtxt(path+f, unpack = True)
    data = {"lam":lam, "lam_tot":lam_tot, "K":K, "rep":nrep}
    df = pd.DataFrame(data)
    
    return df


##########Modified, no path anymore... 
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