# Data-Adaptive Automatic Threshold Calibration for Stability Selection
 
This repository contains code to reproduce results as given in the paper "Data-Adaptive Automatic Threshold Calibration for Stability Selection". 

## Simulations
We provide the four reproducible simulation results in "Artificial-Simulations.R", "Diabetes-Simulations.R", "Proteomics-Simulation.R", and FDR-Simulations.R. Note that running all four scripts will take at least 48 hours without parallelisation. Therefore we also provide the .Rdata files, which can be found in the "Data" folder.

To reproduce the figures, we provide RData files in the "Data" folder, and also provide code in their respective scripts. I.e. to reproduce figures in the Proteomics simulations, the code are found at the end of "Proteomics-Simulation.R".

All of these scripts require the use of user-defined functions given in "Functions.R".

Please note that we are aware of an issue whereby the `stabsel` function does not work in Positron. To avoid this issue, we recommend executing the code in RStudio.


## Data Availability

We provide the plasma proteomics dataset as "Proteomics.csv" in the folder "Data", and was sourced from Rumer et al. (2022). The diabetes dataset can be accessed through the R-package `lars` and `data(diabetes)`. The artificial data generation process can be located in "Functions.R"

## Package References



Benjamin Hofner and Torsten Hothorn (2021). stabs: Stability Selection with Error Control, R package version 0.6-4 https://CRAN.R-project.org/package=stabs.

Breheny P, Huang J (2011). “Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection.” _Annals of Applied Statistics_, *5*(1), 232-253. doi:10.1214/10-AOAS388 <https://doi.org/10.1214/10-AOAS388>, <https://doi.org/10.1214/10-AOAS388>.

Eddelbuettel D, Francois R, Allaire J, Ushey K, Kou Q, Russell N, Ucar I, Bates D, Chambers J (2025). _Rcpp: Seamless R and C++ Integration_. R package version 1.0.14, <https://CRAN.R-project.org/package=Rcpp>.

Genz A, Bretz F (2009). _Computation of Multivariate Normal and t Probabilities_, series Lecture Notes in Statistics. Springer-Verlag, Heidelberg. ISBN 978-3-642-01688-2.

Hastie T, Efron B (2022). _lars: Least Angle Regression, Lasso and Forward Stagewise_. R package version 1.3, <https://CRAN.R-project.org/package=lars>.

Meschiari S (2022). _latex2exp: Use LaTeX Expressions in Plots_. R package version 0.9.6, <https://CRAN.R-project.org/package=latex2exp>.

Patterson E, Sesia M (2022). _knockoff: The Knockoff Filter for Controlled Variable Selection_. R package version 0.3.6, <https://CRAN.R-project.org/package=knockoff>.

Pedersen T (2024). _patchwork: The Composer of Plots_. R package version 1.3.0, <https://CRAN.R-project.org/package=patchwork>.

Ruben Dezeure, Peter Buehlmann, Lukas Meier and Nicolai Meinshausen (2015). High-Dimensional Inference: Confidence Intervals, p-values and R-Software hdi. Statistical Science 30 (4), 533-558.

Rudis B (2024). _hrbrthemes: Additional Themes, Theme Components and Utilities for 'ggplot2'_. R package version 0.8.7, <https://CRAN.R-project.org/package=hrbrthemes>.

Solymos P, Zawadzki Z (2023). _pbapply: Adding Progress Bar to '*apply' Functions_. R package version 1.7-2, <https://CRAN.R-project.org/package=pbapply>.

Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” _Journal of Open Source Software_, *4*(43), 1686. doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.