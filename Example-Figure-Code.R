# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code to reproduce Section 4 Figures                             #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2025)                      #                                    
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("Functions.R")
load("RData/r3sub.Rdata")
load("RData/r7sub.Rdata")
load("RData/S1hard.Rdata")
load("RData/S2hard.Rdata")
load("RData/S3hard.Rdata")
load("RData/S4hard.Rdata")
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# MCC Artificial Data
filtered = c("ATS", "Exclusion ATS",
             "Static 0.60","Static 0.75", "Static 0.90", "LASSO 1SE", "Knockoff", "SCAD")
repeats = 1000
c11hard = combine(S1.10hard, S1.5hard, S1.3hard, S1.1hard, 1, filtered = filtered) |> makeCluster()
c22hard = combine(S2.10hard, S2.5hard, S2.3hard, S2.1hard, 2, filtered = filtered) |> makeCluster()
c33hard = combine(S3.10hard, S3.5hard, S3.3hard, S3.1hard, 3, filtered = filtered) |> makeCluster()
c44hard = combine(S4.10hard, S4.5hard, S4.3hard, S4.1hard, 4, filtered = filtered) |> makeCluster()

c1hard = rbind(c11hard,c22hard) |> totplotnoaxis()
c2hard = rbind(c33hard,c44hard) |> totplot()
library(patchwork)
c1hard/c2hard
