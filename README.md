# PredictiveComputationalPlatformEpigeneticTherapy2026
Codes associated with the paper "A  predictive computational platform for epigenetic therapy optimization in triple-negative breast cancer" by Bruno et al.

More precisely, this repository contains the code used to develop and simulate a computational platform for epigenetic therapy optimization in triple-negative breast cancer , as presented in the associated manuscript.
The platform integrates mechanistic modeling of chromatin dynamics and tumor progression with pharmacological inputs to study treatment effects, resistance mechanisms, and alternative therapeutic strategies.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Repository Structure

•	SB_Comp_Platform_Main_Model.R
Main computational model for combined epigenetic and signaling therapy (tazemetostat + ipatasertib) 

•	SB_Comp_Platform_Resist_Model.R
Extension of the model to study resistance mechanisms 

•	SB_Comp_Platform_Capiv.R
Alternative treatment scenario using capivasertib

•	posterior_samples.rds
Pre-estimated parameters for the tumor progression model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup

Before running any script, you must update the output directories in each file:
output_location <- "/Users/..."           # to be modified
output_location_figures <- "/Users/..."  # to be modified
Set these paths to your desired local directories where results and figures will be saved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Reproducibility

•	All simulations use fixed parameter values stored in posterior_samples.rds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DOI

This repository is archived on Zenodo: https://zenodo.org/badge/1200784060.svg)](https://doi.org/10.5281/zenodo.19410819
