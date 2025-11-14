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
# load("Data/S1t.RData")
# load("Data/S2t.Rdata")
# load("Data/S3t.Rdata")
# load("Data/S4t.Rdata")

# Setting 1: snr = 0.5
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000; snr = 0.5
d = gendata(n = n, p = p, active = active, distribution = "t")
true = c(rep(1, active), rep(0,p-active))
S1.05hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 1: SNR = 1
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S1.1hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 1: snr = 2
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000; snr = 2
true = c(rep(1, active), rep(0,p - active))
S1.2hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 1: SNR = 3
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000;  snr = 3
true = c(rep(1, active), rep(0,p - active ))
S1.3hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

save(S1.05hard.t, S1.2hard.t, S1.3hard.t, S1.1hard.t, file = paste0(Sys.Date(), "_S1t.RData"))
#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# Setting 2: snr = 0.5
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 0.5
true = c(rep(1, active), rep(0,p - active ))
d = gendata(n = n, p = p, active = active, distribution = "t")
S2.05hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 2: SNR = 1
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S2.1hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 2: snr = 2
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 2
true = c(rep(1, active), rep(0,p - active ))
S2.2hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 2: SNR = 3
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 3
true = c(rep(1, active), rep(0,p - active ))
S2.3hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))
save(S2.05hard.t, S2.2hard.t, S2.3hard.t, S2.1hard.t, file = paste0(Sys.Date(), "_S2t.RData"))
#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# Setting 3: snr = 0.5
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000;   snr = 0.5
true = c(rep(1, active), rep(0,p - active ))
d = gendata(n = n, p = p, active = active, distribution = "t")
S3.05hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 3: SNR = 1
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S3.1hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 3: snr = 2
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000;  snr = 2
true = c(rep(1, active), rep(0,p - active))
S3.2hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 3: SNR = 3
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000;  snr = 3
true = c(rep(1, active), rep(0,p - active ))
S3.3hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))
save(S3.05hard.t, S3.2hard.t, S3.3hard.t ,S3.1hard.t, file = paste0(Sys.Date(), "_S3t.RData"))
#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# Setting 4: snr = 0.5
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 0.5
true = c(rep(1, active), rep(0,p - active ))
d = gendata(n = n, p = p, active = active, distribution = "t")
S4.05hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 4: SNR = 1
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S4.1hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 4: snr = 2
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 2
true = c(rep(1, active), rep(0,p - active ))
S4.2hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 4: SNR = 3
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 3
true = c(rep(1, active), rep(0,p - active ))
S4.3hard.t = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))
save(S4.05hard.t, S4.2hard.t, S4.3hard.t, S4.1hard.t, file = paste0(Sys.Date(), "_S4t.RData"))





#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
# Figures
#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
# MCC
c11hard.t = combine(S1.05hard.t, S1.1hard.t, S1.2hard.t, S1.3hard.t,  1, filtered = filtered) |> makeCluster()
c22hard.t = combine(S2.05hard.t, S2.1hard.t, S2.2hard.t, S2.3hard.t,  2, filtered = filtered) |> makeCluster()
c33hard.t = combine(S3.05hard.t, S3.1hard.t, S3.2hard.t, S3.3hard.t,  3, filtered = filtered) |> makeCluster()
c44hard.t = combine(S4.05hard.t, S4.1hard.t, S4.2hard.t, S4.3hard.t,  4, filtered = filtered) |> makeCluster()

c1hard.t = rbind(c11hard.t,c22hard.t) |> totplotnoaxis()
c2hard.t = rbind(c33hard.t,c44hard.t) |> totplot()

c1hard.t/c2hard.t
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Variables Selected
c11Nhardt = combineN(S1.05hard.t, S1.1hard.t, S1.2hard.t, S1.3hard.t,  1, filtered = filtered) |> makeCluster()
c22Nhardt = combineN(S2.05hard.t, S2.1hard.t, S2.2hard.t, S2.3hard.t,  2, filtered = filtered) |> makeCluster()
c33Nhardt = combineN(S3.05hard.t, S3.1hard.t, S3.2hard.t, S3.3hard.t,  3, filtered = filtered) |> makeCluster()
c44Nhardt = combineN(S4.05hard.t, S4.1hard.t, S4.2hard.t, S4.3hard.t,  4, filtered = filtered) |> makeCluster()

c1Nhardt = rbind(c11Nhardt,c22Nhardt) |> totplotNnoaxis(lim = 50)
c2Nhardt = rbind(c33Nhardt,c44Nhardt) |> totplotN(lim = 60)

c1Nhardt/c2Nhardt
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Recall
c11Recallhardt = combineRecall(S1.05hard.t,S1.1hard.t, S1.2hard.t, S1.3hard.t,  1, filtered = filtered) |> makeCluster()
c22Recallhardt = combineRecall(S2.05hard.t,S2.1hard.t, S2.2hard.t, S2.3hard.t,  2, filtered = filtered) |> makeCluster()
c33Recallhardt = combineRecall(S3.05hard.t,S3.1hard.t, S3.2hard.t, S3.3hard.t,  3, filtered = filtered) |> makeCluster()
c44Recallhardt = combineRecall(S4.05hard.t,S4.1hard.t, S4.2hard.t, S4.3hard.t,  4, filtered = filtered) |> makeCluster()

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
#Precision 
c11Precisionhardt = combinePrecision(S1.05hard.t, S1.1hard.t,S1.2hard.t, S1.3hard.t, 1, filtered = filtered) |> makeCluster()
c22Precisionhardt = combinePrecision(S2.05hard.t, S2.1hard.t,S2.2hard.t, S2.3hard.t, 2, filtered = filtered) |> makeCluster()
c33Precisionhardt = combinePrecision(S3.05hard.t, S3.1hard.t,S3.2hard.t, S3.3hard.t, 3, filtered = filtered) |> makeCluster()
c44Precisionhardt = combinePrecision(S4.05hard.t, S4.1hard.t,S4.2hard.t, S4.3hard.t, 4, filtered = filtered) |> makeCluster()

phard1 = prplot1(c11Recallhardt, c11Precisionhardt, c22Recallhardt, c22Precisionhardt)
phard2 = prplot1(c33Recallhardt, c33Precisionhardt, c44Recallhardt, c44Precisionhardt)

phard1
phard2
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Exclusion probability threshold
s1.05e.t = exclusion(S1.05hard.t)
s1.2e.t =  exclusion(S1.2hard.t)
s1.3e.t =  exclusion(S1.3hard.t)
s1.1e.t =  exclusion(S1.1hard.t)

s2.05e.t = exclusion(S2.05hard.t)
s2.2e.t =  exclusion(S2.2hard.t)
s2.3e.t =  exclusion(S2.3hard.t)
s2.1e.t =  exclusion(S2.1hard.t)

s3.05e.t =  exclusion(S3.05hard.t)
s3.2e.t =   exclusion(S3.2hard.t)
s3.3e.t =   exclusion(S3.3hard.t)
s3.1e.t =   exclusion(S3.1hard.t)

s4.05e.t = exclusion(S4.05hard.t)
s4.2e.t =  exclusion(S4.2hard.t)
s4.3e.t =  exclusion(S4.3hard.t)
s4.1e.t =  exclusion(S4.1hard.t)

excl = data.frame("Eta" = 
                    c(s1.05e.t, s1.1e.t, s1.2e.t, s1.3e.t ,
                      s2.05e.t, s2.1e.t, s2.2e.t, s2.3e.t ,
                      s3.05e.t, s3.1e.t, s3.2e.t, s3.3e.t ,
                      s4.05e.t, s4.1e.t, s4.2e.t, s4.3e.t ),
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

S1pi = rbind(extractPi(S1.05hard.t, 0.5, 1),
             extractPi(S1.1hard.t, 1, 1),
             extractPi(S1.2hard.t, 2, 1),
             extractPi(S1.3hard.t, 3, 1))

S2pi = rbind(extractPi(S2.05hard.t, 0.5, 2),
             extractPi(S2.1hard.t, 1, 2),
             extractPi(S2.2hard.t, 2, 2),
             extractPi(S2.3hard.t, 3, 2))

S3pi = rbind(extractPi(S3.05hard.t, 0.5, 3),
             extractPi(S3.1hard.t, 1, 3),
             extractPi(S3.2hard.t, 2, 3),
             extractPi(S3.3hard.t, 3, 3))

S4pi = rbind(extractPi(S4.05hard.t, 0.5, 4),
             extractPi(S4.1hard.t, 1, 4),
             extractPi(S4.2hard.t, 2, 4),
             extractPi(S4.3hard.t, 3, 4))

Spi = rbind(S1pi, S2pi, S3pi, S4pi)


ggplot(Spi, aes(x = value, fill = variable)) + geom_density(alpha = 0.7) +
  facet_grid(dimension~SNR, labeller = label_parsed) + xlab(TeX("Estimated $\\pi$")) + ylab("Frequency") +
  theme_few_grid(base_size = 20) + scale_fill_manual(values = c("#FC8D62", "#FFD92F")) + labs(fill = "Method")
