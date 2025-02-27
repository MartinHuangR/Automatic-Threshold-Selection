# Data-Adaptive Automatic Threshold Calibration for Stability Selection
 
This repository contains code to reproduce results as given in the paper "Data-Adaptive Automatic Threshold Calibration for Stability Selection". 

## Simulations
We provide the four reproducible simulation results in "Artificial-Simulations.R", "Diabetes-Simulations.R", "Proteomics-Simulation.R", and FDR-Simulations.R. Note that running all four scripts will take at least 48 hours without parallelisation. Therefore we also provide the .Rdata files, which can be found in the "Data" folder.

To reproduce the figures, we provide RData files in the "Data" folder, and also provide code in their respective scripts. I.e. to reproduce figures in the Proteomics simulations, the code are found at the end of "Proteomics-Simulation.R".

All of these scripts require the use of user-defined functions given in "Functions.R".

Please note that we are aware of an issue whereby the `stabsel` function does not work in Positron. To avoid this issue, we recommend executing the code in RStudio.


## Data Availability

We provide the plasma proteomics dataset as "Proteomics.csv" in the folder "Data", and was sourced from Rumer et al. (2022). The diabetes dataset can be accessed through the R-package `lars` and `data(diabetes)`. The artificial data generation process can be located in "Functions.R"
