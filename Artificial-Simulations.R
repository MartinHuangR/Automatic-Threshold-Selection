# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4 Artificial Simulations      #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2025)                      #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("Functions.R")
set.seed(1)


# Each Rdata file contains 4 SNRs for their respective simulation.
load("RData/S1hard.Rdata") 
load("RData/S2hard.Rdata")
load("RData/S3hard.Rdata")
load("RData/S4hard.Rdata")

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Simulation 1: SNR = 10
n = 20; p = 1000; active = 2; repeats = 1000; snr = 10
d = gendata(n = n, p = p, active = active)
true = c(rep(1, active), rep(0,p-active))
S1.10hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 1: SNR = 5
n = 20; p = 1000; active = 2; repeats = 1000; snr = 5
true = c(rep(1, active), rep(0,p - active))
S1.5hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 1: SNR = 3
n = 20; p = 1000; active = 2; repeats = 1000;  snr = 3
true = c(rep(1, active), rep(0,p - active ))
S1.3hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 1: SNR = 1
n = 20; p = 1000; active = 2; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S1.1hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))


#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Simulation 2: SNR = 10
n = 200; p = 200; active = 20; repeats = 1000;   snr = 10
true = c(rep(1, active), rep(0,p - active ))
d =  gendata(n = n, p = p, active = active)
S2.10hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 2: SNR = 5
n = 200; p = 200; active = 20; repeats = 1000;  snr = 5
true = c(rep(1, active), rep(0,p - active))
S2.5hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 2: SNR = 3
n = 200; p = 200; active = 20; repeats = 1000;  snr = 3
true = c(rep(1, active), rep(0,p - active ))
S2.3hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 2: SNR = 1
n = 200; p = 200; active = 20; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S2.1hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))


#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Simulation 3: SNR = 10
n = 500; p = 100; active = 20; repeats = 1000; snr = 10
true = c(rep(1, active), rep(0,p - active ))
d =  gendata(n = n, p = p, active = active)
S3.10hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 3: SNR = 5
n = 500; p = 100; active = 20; repeats = 1000; snr = 5
true = c(rep(1, active), rep(0,p - active ))
S3.5hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 3: SNR = 3
n = 500; p = 100; active = 20; repeats = 1000; snr = 3
true = c(rep(1, active), rep(0,p - active ))
S3.3hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 3: SNR = 1
n = 500; p = 100; active = 20; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S3.1hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Simulation 4: SNR = 10
n = 100; p = 500; active = 10; repeats = 1000; snr = 10
true = c(rep(1, active), rep(0,p - active ))
d =  gendata(n = n, p = p, active = active)
S4.10hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 4: SNR = 5
n = 100; p = 500; active = 10; repeats = 1000; snr = 5
true = c(rep(1, active), rep(0,p - active ))
S4.5hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 4: SNR = 3
n = 100; p = 500; active = 10; repeats = 1000; snr = 3
true = c(rep(1, active), rep(0,p - active ))
S4.3hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 4: SNR = 1
n = 100; p = 500; active = 10; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S4.1hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))


