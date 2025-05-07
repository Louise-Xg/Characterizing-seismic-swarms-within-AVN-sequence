# Open reasearch for "Characterizing seismic swarm activity in the aftershock zone at the 2016 Amatrice-Norcia seismic sequence"
@author:Louise Xiang

Date: 13/05/2024

To: Journal of Geophysical Research - Solid Earth

This work is part of the thesis of Louise Xiang, financed by a ministerial grant from the French Ministry of Higher Education, Research, and Innovation (MERSI).

## Open Research:
AGU requires that the underlying data and/or software or code needed to understand, evaluate, and build upon the reported research be available at the time of peer review and publication. Additionally, the code (e.g. Python, Jupyter Notebooks, R, MATLAB) should be made available in a free and open platform (e.g., Github) and preserved in a repository (e.g. Zenodo)

Therefore here I cite datasets, part of my scripts used in my study:

Datasets (tables): 
- Results obtained after applying Density-based clustering (3D): "v0_d03km_n50_every1.mat" 
- Results of final seismic swarms: "df_final_swarms_of_v0_d03_n50_every1.csv" (correspond to dataset s01 in the Supplementary Information)

Methods:
- Density-based clustering (3D) with matlab (@David Marsan): "cout_local.m"
                                                             "compute_pi_M.m"
- ETAS model (4D) including earthquake detection probability with C (@David Marsan): "etas-amatrice-v5.c"
- Determination of candidate seismic swarms with jupyter notebook: "Generating_seismicity_rate.ipynb"


