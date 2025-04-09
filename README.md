This GitHub repository contains the code used for the DPMM algorithm and the Composite Individual-Level Model (C-ILM) framework, including both simulation studies and the real-world data application described in our ArXiv paper.

All code was developed and run in R version 4.3.1.

Code was run in R v4.3.1

File description:

- Simulation_data/: contains simulated epidemic data generated using Generate_epidemics.R and generate.cpp.
- log_likelihood_*.cpp: C++ code for computing the likelihood functions for DPMM and ILMs under basic SILM, C-ILM, and SEIR models.
- DPMM_cleaned.R: R implementation of the DPMM-C-ILM MCMC algorithm (Algorithm 1 in the paper).
- BasicILM.stan: Stan code for the basic SILM model.
- CILM.stan*: Stan code for C-ILM variants (Models M2, M3, and M4).
- CILM_cleaned.R: R code to run Stan-based C-ILM models (M2–M4).
- PPD_*.R: R scripts for generating posterior predictive distributions.
- FMD_Cumbria_Subset.csv: Subset of the 2001 Foot-and-Mouth Disease (FMD) outbreak data used in the case study.
- FMD_Application.R: Data processing and application of the model to the 2001 FMD dataset.
- Plots/: Code used to generate the figures included in the paper.
