# Garrett_Glitz_etal

Code for the preprint/submission: <i> Model sharing
in the human medial temporal lobe </i>; Glitz, Juechems, Summerfield & Garrett; 2021

This repository contains analysis scripts/functions for all analyses detailed in the main text of the manuscript and most supplementary analyses. 

<b> Requirements: </b> to run all scripts in this repository, you will need to have installed Matlab, Julia and R; 

<b> Data: </b> behavioural data is available in this folder under data/behavioural_data/; if you want to look at the second level SPM maps or re-run the non-crossvalidated RSA analyses or the encoding model analysis, you will need to download the neural data from this repository https://osf.io/zvkj3/ 

<b> Navigating this repository: </b> this repository is split up into several subfolders, the names of which detail the type of analyses within each. Learning_model refers to the computational model used to compute subjective transition probabilities. Each analysis script should contain a header which will tell you which paths to change and which inputs are required if you want to run the code on your computer.

Other: the plots in Figures 3d and 4 are recreated when you run the function /multivariate_analyses/RSA_analyses/without crossvalidation/RSA_onset.m with the data for the different brain regions from OSF. The boxplot in Figure 5b can be recreated when running multivariate_analyses/RSA_analyses/without crossvalidation/RSA_outcome.m. Figure 6a and b can be recreated by running valence_plots.R



