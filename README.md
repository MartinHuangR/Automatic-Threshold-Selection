# Data-Adaptive Automatic Threshold Calibration for Stability Selection
 
This repository contains code to reproduce results as given in the paper "Data-Adaptive Automatic Threshold Calibration for Stability Selection". 

## Simulations
We provide the four reproducible simulation results in "Artificial-Simulations.R", "Riboflavin-Simulations.R", "Proteomics-Simulation.R", and FDR-Simulations.R. Note that running all four scripts will take at least 48 hours. Therefore we also provide the .Rdata files, which can be found in the "RData" folder.

To reproduce the figures, we provide an example given in "Example-Figure-Code.R".

All of these scripts require the use of user-defined functions given in "Functions.R".

Please note that we are aware of an issue whereby the `stabsel` function does not work in Positron. To avoid this issue, we recommend executing the code in RStudio.


## Data Availability

We provide the plasma proteomics dataset as "Proteomics.csv", and was sourced from Rumer et al. (2022). The riboflavin dataset can be accessed through the R-package `hdi` and `data(riboflavin)`. The artificial data generation process can be located in "Functions.R"
