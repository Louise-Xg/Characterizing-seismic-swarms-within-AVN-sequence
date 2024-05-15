# Open reasearch for "Characterizing seismic swarm activity in the aftershock zone at the 2016 Amatrice-Norcia seismic sequence"
@Louise Xiang
Date: 13/05/2024
To: Journal of Geophysical Research - Solid Earth

## Open Research:
AGU requires that the underlying data and/or software or code needed to understand, evaluate, and build upon the reported research be available at the time of peer review and publication. Additionally, the code (e.g. Python, Jupyter Notebooks, R, MATLAB) used to perform any data analysis and to produce the manuscript’s figures should be made available in a free and open platform (e.g. Github) and preserved in a repository (e.g. Zenodo)

Therefore here I cite datasets, softwares and codes for analysis used in my study:

Datasets (tables): 
- Results obtained after applying Density-based clustering (3D): "v0_d03km_n50_every1.mat"
- Results of final seismic swarms: "df_final_swarms_of_v0_d03_n50_every1.csv"

Methods:
- Density-based clustering (3D) with matlab (@David Marsan): "cout_local.m"
                                                             "compute_pi_M.m"
- ETAS model (4D) including earthquake detection probability with C (@David Marsan): "etas-amatrice-v5.c"
- Determination of candidate seismic swarms with jupyter notebook: "automatisation_looking_for_swarms.ipynb"
- Statistical assessment of candidate seismic swarms with jupyter notebook: "take_a_look_at_swarms_of_v0_d03_n50_every1_review.ipynb"

Analyses:
- Migration analysis of seismic swarms with jupyter notebook: "Migration_behaviours.ipynb"
- Global statistics of final seismic swarms with jupyter notebook: "check_up_global_stats_of_v0_d03_n50_every1.ipynb"
- Interpretation of our results (figures mainly) with jupyter notebook: "plots_for_paper_about_swarms_v0_d03_n50_every1.ipynb"
                                                                        "map_results_for_v0_d03_n50_every1.ipynb"
                                                                        "cross_sections_of_v0_d03_n50_every1.ipynb"
                                                                        "Interpretation_swarms_N_S_parts_of_v0_d03_n50_every1.ipynb"
                                                                        "compute_pressure_variation.ipynb"

