# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4.1                           #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2025)                      #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #  
#                                                                   #
# Note: This will take more than 24 hours without parallelisation.  #
#       I have provided the Rdata file for my simulations.          #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("Functions.R")


# Each Rdata file contains 4 SNRs for their respective simulation.
# load("Data/S1U.RData")
# load("Data/S2U.Rdata")
# load("Data/S3U.Rdata")
# load("Data/S4U.Rdata")

# Setting 1: snr = 0.5
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000; snr = 0.5
d = gendata(n = n, p = p, active = active, distribution = "uniform")
true = c(rep(1, active), rep(0,p-active))
S1.05hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 1: SNR = 1
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S1.1hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 1: snr = 2
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000; snr = 2
true = c(rep(1, active), rep(0,p - active))
S1.2hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 1: SNR = 3
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000;  snr = 3
true = c(rep(1, active), rep(0,p - active ))
S1.3hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

save(S1.05hard.U, S1.2hard.U, S1.3hard.U, S1.1hard.U, file = paste0(Sys.Date(), "_S1U.RData"))
#---#---#---#---#---#---#---#---#---#---#---#---#---#---

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
#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# Setting 3: snr = 0.5
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000;   snr = 0.5
true = c(rep(1, active), rep(0,p - active ))
d = gendata(n = n, p = p, active = active, distribution = "uniform")
S3.05hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 3: SNR = 1
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S3.1hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 3: snr = 2
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000;  snr = 2
true = c(rep(1, active), rep(0,p - active))
S3.2hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 3: SNR = 3
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000;  snr = 3
true = c(rep(1, active), rep(0,p - active ))
S3.3hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))
save(S3.05hard.U, S3.2hard.U, S3.3hard.U ,S3.1hard.U, file = paste0(Sys.Date(), "_S3U.RData"))
#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# Setting 4: snr = 0.5
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 0.5
true = c(rep(1, active), rep(0,p - active ))
d = gendata(n = n, p = p, active = active, distribution = "uniform")
S4.05hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 4: SNR = 1
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S4.1hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 4: snr = 2
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 2
true = c(rep(1, active), rep(0,p - active ))
S4.2hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 4: SNR = 3
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 3
true = c(rep(1, active), rep(0,p - active ))
S4.3hard.U = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))
save(S4.05hard.U, S4.2hard.U, S4.3hard.U, S4.1hard.U, file = paste0(Sys.Date(), "_S4U.RData"))





#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
# Figures
#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
# MCC
c11hard.U = combine(S1.05hard.U, S1.1hard.U, S1.2hard.U, S1.3hard.U,  1, filtered = filtered) |> makeCluster()
c22hard.U = combine(S2.05hard.U, S2.1hard.U, S2.2hard.U, S2.3hard.U,  2, filtered = filtered) |> makeCluster()
c33hard.U = combine(S3.05hard.U, S3.1hard.U, S3.2hard.U, S3.3hard.U,  3, filtered = filtered) |> makeCluster()
c44hard.U = combine(S4.05hard.U, S4.1hard.U, S4.2hard.U, S4.3hard.U,  4, filtered = filtered) |> makeCluster()

c1hard.U = rbind(c11hard.U,c22hard.U) |> totplotnoaxis()
c2hard.U = rbind(c33hard.U,c44hard.U) |> totplot()

c1hard.U/c2hard.U
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Variables Selected
c11NhardU = combineN(S1.05hard.U, S1.1hard.U, S1.2hard.U, S1.3hard.U,  1, filtered = filtered) |> makeCluster()
c22NhardU = combineN(S2.05hard.U, S2.1hard.U, S2.2hard.U, S2.3hard.U,  2, filtered = filtered) |> makeCluster()
c33NhardU = combineN(S3.05hard.U, S3.1hard.U, S3.2hard.U, S3.3hard.U,  3, filtered = filtered) |> makeCluster()
c44NhardU = combineN(S4.05hard.U, S4.1hard.U, S4.2hard.U, S4.3hard.U,  4, filtered = filtered) |> makeCluster()

c1NhardU = rbind(c11NhardU,c22NhardU) |> totplotNnoaxis(lim = 50)
c2NhardU = rbind(c33NhardU,c44NhardU) |> totplotN(lim = 60)

c1NhardU/c2NhardU
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Recall
c11RecallhardU = combineRecall(S1.05hard.U,S1.1hard.U, S1.2hard.U, S1.3hard.U,  1, filtered = filtered) |> makeCluster()
c22RecallhardU = combineRecall(S2.05hard.U,S2.1hard.U, S2.2hard.U, S2.3hard.U,  2, filtered = filtered) |> makeCluster()
c33RecallhardU = combineRecall(S3.05hard.U,S3.1hard.U, S3.2hard.U, S3.3hard.U,  3, filtered = filtered) |> makeCluster()
c44RecallhardU = combineRecall(S4.05hard.U,S4.1hard.U, S4.2hard.U, S4.3hard.U,  4, filtered = filtered) |> makeCluster()

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
#Precision 
c11PrecisionhardU = combinePrecision(S1.05hard.U, S1.1hard.U,S1.2hard.U, S1.3hard.U, 1, filtered = filtered) |> makeCluster()
c22PrecisionhardU = combinePrecision(S2.05hard.U, S2.1hard.U,S2.2hard.U, S2.3hard.U, 2, filtered = filtered) |> makeCluster()
c33PrecisionhardU = combinePrecision(S3.05hard.U, S3.1hard.U,S3.2hard.U, S3.3hard.U, 3, filtered = filtered) |> makeCluster()
c44PrecisionhardU = combinePrecision(S4.05hard.U, S4.1hard.U,S4.2hard.U, S4.3hard.U, 4, filtered = filtered) |> makeCluster()

phard1 = prplot1(c11RecallhardU, c11PrecisionhardU, c22RecallhardU, c22PrecisionhardU)
phard2 = prplot1(c33RecallhardU, c33PrecisionhardU, c44RecallhardU, c44PrecisionhardU)

phard1
phard2
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Exclusion probability threshold
s1.05e.U = exclusion(S1.05hard.U)
s1.2e.U =  exclusion(S1.2hard.U)
s1.3e.U =  exclusion(S1.3hard.U)
s1.1e.U =  exclusion(S1.1hard.U)

s2.05e.U = exclusion(S2.05hard.U)
s2.2e.U =  exclusion(S2.2hard.U)
s2.3e.U =  exclusion(S2.3hard.U)
s2.1e.U =  exclusion(S2.1hard.U)

s3.05e.U =  exclusion(S3.05hard.U)
s3.2e.U =   exclusion(S3.2hard.U)
s3.3e.U =   exclusion(S3.3hard.U)
s3.1e.U =   exclusion(S3.1hard.U)

s4.05e.U = exclusion(S4.05hard.U)
s4.2e.U =  exclusion(S4.2hard.U)
s4.3e.U =  exclusion(S4.3hard.U)
s4.1e.U =  exclusion(S4.1hard.U)

excl = data.frame("Eta" = 
                    c(s1.05e.U, s1.1e.U, s1.2e.U, s1.3e.U ,
                      s2.05e.U, s2.1e.U, s2.2e.U, s2.3e.U ,
                      s3.05e.U, s3.1e.U, s3.2e.U, s3.3e.U ,
                      s4.05e.U, s4.1e.U, s4.2e.U, s4.3e.U ),
                  "Setting" = c(rep("(I):~n==20*`,`~p==1000*`,`~`|`*beta[S]*`|`==2",4000),
                                rep("(II):~n==100*`,`~p==500*`,`~`|`*beta[S]*`|`==10",4000),
                                rep("(III):~n==200*`,`~p==200*`,`~`|`*beta[S]*`|`==20",4000),
                                rep("(IV):~n==500*`,`~p==100*`,`~`|`*beta[S]*`|`==20",4000)),
                  "SNR" = rep(c(rep(0.5, 1000), rep(1,1000), rep(2, 1000), rep(3,1000)), 4)
) |> mutate(SNR = factor(SNR))

ggplot(excl, aes(x = SNR, y = Eta)) + geom_boxplot(alpha = 0.65) +
  facet_wrap(~Setting, labeller = label_parsed, ncol = 4) +
  theme_few_grid(base_size = 20) +
  ylab(TeX("Exclusion Probability Threshold $\\eta$")) + 
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  theme(legend.position = "none") 

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# ATS/EATS Estimated Pi

S1pi = rbind(extractPi(S1.05hard.U, 0.5, 1),
             extractPi(S1.1hard.U, 1, 1),
             extractPi(S1.2hard.U, 2, 1),
             extractPi(S1.3hard.U, 3, 1))

S2pi = rbind(extractPi(S2.05hard.U, 0.5, 2),
             extractPi(S2.1hard.U, 1, 2),
             extractPi(S2.2hard.U, 2, 2),
             extractPi(S2.3hard.U, 3, 2))

S3pi = rbind(extractPi(S3.05hard.U, 0.5, 3),
             extractPi(S3.1hard.U, 1, 3),
             extractPi(S3.2hard.U, 2, 3),
             extractPi(S3.3hard.U, 3, 3))

S4pi = rbind(extractPi(S4.05hard.U, 0.5, 4),
             extractPi(S4.1hard.U, 1, 4),
             extractPi(S4.2hard.U, 2, 4),
             extractPi(S4.3hard.U, 3, 4))

Spi = rbind(S1pi, S2pi, S3pi, S4pi)


ggplot(Spi, aes(x = value, fill = variable)) + geom_density(alpha = 0.7) +
  facet_grid(dimension~SNR, labeller = label_parsed) + xlab(TeX("Estimated $\\pi$")) + ylab("Frequency") +
  theme_few_grid(base_size = 20) + scale_fill_manual(values = c("#FC8D62", "#FFD92F")) + labs(fill = "Method")
