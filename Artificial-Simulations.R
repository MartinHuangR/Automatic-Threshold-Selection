# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4 Artificial Simulations      #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2025)                      #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #  
# Note: This will take more than 24 hours without parallelisation.  #
#       I have provided the Rdata file for my simulations.          #
#       Please see Example-Figure-Code.R for code to visualise      #
#       results                                                     #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("Functions.R")


# Each Rdata file contains 4 SNRs for their respective simulation.
load("Data/S1hardR.Rdata") 
load("Data/S2hardR.Rdata")
load("Data/S3hardR.Rdata")
load("Data/S4hardR.Rdata")


source("Functions.R")

set.seed(1)

# Simulation 1: SNR = 10
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000; snr = 10
d = gendata(n = n, p = p, active = active)
true = c(rep(1, active), rep(0,p-active))
S1.10hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 1: SNR = 5
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000; snr = 5
true = c(rep(1, active), rep(0,p - active))
S1.5hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 1: SNR = 3
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000;  snr = 3
true = c(rep(1, active), rep(0,p - active ))
S1.3hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 1: SNR = 1
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S1.1hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))


save(S1.10hard, S1.5hard, S1.3hard, S1.1hard, file = "S1hardR.Rdata")
#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# Simulation 2: SNR = 10
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 10
true = c(rep(1, active), rep(0,p - active ))
d = gendata(n = n, p = p, active = active)
S2.10hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 4: SNR = 5
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 5
true = c(rep(1, active), rep(0,p - active ))
S2.5hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 4: SNR = 3
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 3
true = c(rep(1, active), rep(0,p - active ))
S2.3hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 4: SNR = 1
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S2.1hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

save( S2.10hard,  S2.5hard,  S2.3hard,  S2.1hard, file = "S2hardR.Rdata")
#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# Simulation 3: SNR = 10
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000;   snr = 10
true = c(rep(1, active), rep(0,p - active ))
d = gendata(n = n, p = p, active = active)
S3.10hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 2: SNR = 5
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000;  snr = 5
true = c(rep(1, active), rep(0,p - active))
S3.5hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 2: SNR = 3
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000;  snr = 3
true = c(rep(1, active), rep(0,p - active ))
S3.3hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 2: SNR = 1
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S3.1hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

save(S3.10hard, S3.5hard, S3.3hard ,S3.1hard, file = "S3hardR.Rdata")
#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# Simulation 4: SNR = 10
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 10
true = c(rep(1, active), rep(0,p - active ))
d = gendata(n = n, p = p, active = active)
S4.10hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 3: SNR = 5
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 5
true = c(rep(1, active), rep(0,p - active ))
S4.5hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 3: SNR = 3
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 3
true = c(rep(1, active), rep(0,p - active ))
S4.3hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Simulation 3: SNR = 1
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S4.1hard = replicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

save(S4.10hard, S4.5hard, S4.3hard, S4.1hard, file = "S4hardR.Rdata")

