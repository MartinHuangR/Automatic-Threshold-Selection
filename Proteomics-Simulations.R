# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code to reproduce Section 4 Proteomics                          #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2025)                      #                                  
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
#                                                                   #
# Note: This will take approximately 4 hours to run. I have provided#
# an Rdata file for my simulations. To access MCC use cleanMCC()    #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Each Rdata contains the two SNR values for each number of active variable settings.
load("Data/pro9.Rdata")
load("Data/pro4.Rdata")
source("Functions.R")

pro = read.csv("Data/Proteomics.csv")
pro = pro |> select(-sampleID)
X = pro |> as.matrix() |> scale()

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
set.seed(1)
active = 9; repeats = 1000; snr = 5; p = ncol(X)
beta = c(rep(1,active), rep(0, ncol(X) - active)) 
true = c(rep(1, active),rep(0,ncol(X) - active))
pro9.5 = replicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr, gaussian.knockoffs = T))

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
active = 9; repeats = 1000; snr = 10; p = ncol(X)
beta = c(rep(1,9), rep(0, ncol(X) - active)) 
true = c(rep(1, active),rep(0,ncol(X) - active))
pro9.10 = replicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr, gaussian.knockoffs = T))

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
active = 4; repeats = 1000; snr = 5; p = ncol(X)
beta = c(rep(1,4), rep(0, ncol(X) - active)) 
true = c(rep(1, active),rep(0,ncol(X) - active))
pro4.5 = replicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr, gaussian.knockoffs = T))

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
active = 4; repeats = 1000; snr = 10; p = ncol(X)
beta = c(rep(1,4), rep(0, ncol(X) - active)) 
true = c(rep(1, active),rep(0,ncol(X) - active))
pro4.10 = replicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr, gaussian.knockoffs = T))

