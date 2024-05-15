import glob
import os
import subprocess
import numpy as np
import scipy.io as scp
import pandas as pd

from joblib import Parallel, delayed

#path_to_df0 = "/home/xianglo/Bureau/DATA_AMATRICE/DBSCAN_David/new_all_run_files_for_etas/catalog_all/"
#
path_to_df0 = "/data/local/xianglo/new_all_run_files_for_etas/catalog_all/"

#path_to_dbscan_data = "/home/xianglo/Bureau/DATA_AMATRICE/DBSCAN_David/new_all_run_files_for_etas/dbscan_data/"
#
path_to_dbscan_data = "/data/local/xianglo/new_all_run_files_for_etas/dbscan_data/"

#path_to_parameter_files = "/home/xianglo/Bureau/DATA_AMATRICE/DBSCAN_David/new_all_run_files_for_etas/dbscan_parameter_files_v0_d03_n50_every5/"
#path_to_parameter_files = "/data/local/xianglo/new_all_run_files_for_etas/dbscan_parameter_files_v0_d0_3_n50_every5/"
#
path_to_parameter_files = "/data/local/xianglo/new_all_run_files_for_etas/dbscan_parameter_files_v0_d03_n50_every1/"
#path_to_parameter_files = "/data/local/xianglo/new_all_run_files_for_etas/dbscan_parameter_files_v0_d03_n50_every5/"
if not os.path.exists(path_to_parameter_files):
    os.mkdir(path_to_parameter_files)
    
##Creating directory for output files:
#savedir = "/home/xianglo/Bureau/DATA_AMATRICE/DBSCAN_David/new_all_run_files_for_etas/output_dbscan_v0_d03_n50_every5/"
#savedir = "/data/local/xianglo/new_all_run_files_for_etas/output_dbscan_v0_d0_3_n50_every5/"
#
savedir = "/data/local/xianglo/new_all_run_files_for_etas/output_dbscan_v0_d03_n50_every1/"
#savedir = "/data/local/xianglo/new_all_run_files_for_etas/output_dbscan_v0_d03_n50_every5/"
if not os.path.exists(savedir):
    os.mkdir(savedir)

###Clusters from dbscan:
#v0_n50_d0_3 = scp.loadmat(path_to_dbscan_data+"v0_d0_3km_n50_every1.mat")["C300"].reshape(-1) ##au moins 50 séismes à < 300 m, v = 0km/d : 1 séisme sur 5 comme vertex ##from David ##total 312
#
v0_n50_d0_3 = scp.loadmat(path_to_dbscan_data+"v0_d03km_n50_every1.mat")["clusters"].reshape(-1) ##au moins 50 séismes à < 300 m, v = 0km/d : 1 séisme sur 1 comme vertex ##from me ##total 277
#v0_n50_d0_3 = scp.loadmat(path_to_dbscan_data+"v0_d03km_n50_every5.mat")["clusters"].reshape(-1) ##au moins 50 séismes à < 300 m, v = 0km/d : 1 séisme sur 5 comme vertex ##from me ##total 320

####Raw data:
df0 = pd.read_csv(path_to_df0+'raw_data_modified_version.csv')
###Structure : t, x, y, z, m, x_km, y_km

###Add new column:
df0["v0_n50_d0_3"] = v0_n50_d0_3

###Convert the column in integar:
df0 = df0.astype({"v0_n50_d0_3":"int"})

###Define parameters for loop:
pas = 1
list_i = np.arange(200, df0["v0_n50_d0_3"].max()+pas, pas)

def f(i):
    
    select_EQ_cluster_i = df0[df0["v0_n50_d0_3"] == i]
    t1,t2 = round(select_EQ_cluster_i.t.min(),3), round(select_EQ_cluster_i.t.max(),3)
    x1,x2 = round(select_EQ_cluster_i.x.min(),3), round(select_EQ_cluster_i.x.max(),3)
    y1,y2 = round(select_EQ_cluster_i.y.min(),3), round(select_EQ_cluster_i.y.max(),3)
    z1,z2 = round(select_EQ_cluster_i.z.min(),3), round(select_EQ_cluster_i.z.max(),3)

    ##Create parameters files:
    f = open(path_to_parameter_files+"dbscan_cluster"+str(i), "w")
    f.write('0.  366.\n' + "-5. 2. 5. \n" +str(x1) + ' ' + str(x2) + ' ' + str(y1) + ' ' + str(y2) + ' ' + str(z1) + ' ' + str(z2) + ' ' +'\n 2. 1.001 1e-5 1.80e-5 3. \n')
    f.close()
    
    ##Run inversion:
    process = subprocess.run(["./etas-amatrice-v5", str(path_to_df0)+"data.dat", str(path_to_parameter_files)+"dbscan_cluster"+str(i), str(savedir)+"dbscan_clust"+str(i)])
    
# Parallelize
Parallel(n_jobs=len(list_i), verbose=100)(delayed(f)(i) for i in list_i) 

# avec n_jobs = le nombre de valeurs prises par ma boucle


