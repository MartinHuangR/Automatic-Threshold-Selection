# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code to reproduce Section 4 Riboflavin                          #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2025)                      #                                    
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
#                                                                   #
# Note: This will take approximately 1 hour to run. I have provided #
# an Rdata file for my simulations. To access MCC use cleanMCC()    #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Each Rdata contains the two SNR values for each number of active variable settings.
load("RData/r3sub.Rdata")
load("RData/r7sub.Rdata")
source("Functions.R")
data(riboflavin)
X = data.matrix(riboflavin$x)
y = riboflavin$y
yxdf <- as.data.frame(cbind(y, X))
names(yxdf)[1] <- "y"
X = yxdf[,-1]
X = scale(X)

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
set.seed(10)
X = X[,sample(1:ncol(X), size = 20, replace = F)]


#---#---#---#---#---#---#---#---#---#---#---#---#---#---
active = 3; repeats = 1000;  snr = 5; p = 20
beta = c(rep(1,active), rep(0, ncol(X) - active))
true = c(rep(1, active),rep(0,ncol(X) - active ))
r3.5.sub = replicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr))

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
active = 3; repeats = 1000;  snr = 10; p = 20
beta = c(rep(1,active), rep(0, ncol(X) - active))
true = c(rep(1, active),rep(0,ncol(X) - active ))
r3.10.sub = replicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr))

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
active = 7; repeats = 1000;  snr = 5; p = 20
beta = c(rep(1,active), rep(0, ncol(X) - active))
true = c(rep(1, active),rep(0,ncol(X) - active ))
r7.5.sub = replicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr))

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
active = 7; repeats = 1000;  snr = 10; p = 20
beta = c(rep(1,active), rep(0, ncol(X) - active))
true = c(rep(1, active),rep(0,ncol(X) - active ))
r7.10.sub = replicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr))

