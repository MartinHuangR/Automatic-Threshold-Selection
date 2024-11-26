# Automatic Threshold Selection
 
This repository contains code to reproduce results as given in the paper "Data-Adaptive and Automatic Stable Threshold Calibration for Stability Selection in Penalised Regression". 

## Simulations
We provide the three reproducible simulation results in "Artificial-Simulations.R", "Riboflavin-Simulations.R" and "Proteomics-Simulation.R". Note that running all three scripts will take at least 48 hours. Therefore we also provide the .Rdata files, which can be found in the "RData" folder.

To reproduce the figures, we provide an example given in "Example-Figure-Code.R".

All of these scripts require the use of user-defined functions given in "Functions.R".

## Data Availability

We provide the plasma proteomics dataset as "Proteomics.csv". The riboflavin dataset can be accessed through the R-package `hdi` and `data(riboflavin)`. 