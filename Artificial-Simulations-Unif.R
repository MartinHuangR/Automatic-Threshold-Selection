# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4.1                           #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2026)                      #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #  
#                                                                   #
# The following script provides artificial simulations with         #
# Uniform design                                                    #
#                                                                   #
# Note: I have provided the Rdata file for my simulations.          #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("Functions.R")


# Each Rdata file contains 4 SNRs for their respective simulation.
# load("Data/S2U.Rdata")

# Setting 2: snr = 0.5
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 0.5
true = c(rep(1, active), rep(0,p - active ))
d = gendata(n = n, p = p, active = active, distribution = "uniform")
S2.05hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 2: SNR = 1
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S2.1hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 2: snr = 2
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 2
true = c(rep(1, active), rep(0,p - active ))
S2.2hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 2: SNR = 3
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 3
true = c(rep(1, active), rep(0,p - active ))
S2.3hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))
save(S2.05hard.U, S2.2hard.U, S2.3hard.U, S2.1hard.U, file = paste0(Sys.Date(), "_S2U.RData"))

#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
# Figures
#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
methods = c("ATS", "EATS", "Static 0.75", "LASSO", "Knockoff", "SCAD")
# Prepare MCC
c22hard.U = combine(S2.05hard.U, S2.1hard.U, S2.2hard.U, S2.3hard.U,  2, filtered = filtered) |> makeCluster() |> filter(Method %in% methods)
c22NhardU = combineN(S2.05hard.U, S2.1hard.U, S2.2hard.U, S2.3hard.U,  2, filtered = filtered) |> makeCluster() |> filter(Method %in% methods)

c22hard.U$Dimension = "(Uniform):~n==100*`,`~p==500*`,`~`|`*beta[S]*`|`==10"
c22NhardU$Dimension = "(Uniform):~n==100*`,`~p==500*`,`~`|`*beta[S]*`|`==10"

# Plot MCC
a = c22hard.U |> totMEANplot()
b = c22NhardU |> NtotMEANplot(lim = 50)

#18x8
(a|b) + plot_layout(guides = "collect")
