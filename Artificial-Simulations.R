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
# load("Data/S1.Rdata")
# load("Data/S2.Rdata")
# load("Data/S3.Rdata")
# load("Data/S4.Rdata")

# Setting 1: snr = 0.5
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000; snr = 0.5
d = gendata(n = n, p = p, active = active)
true = c(rep(1, active), rep(0,p-active))
S1.05hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 1: SNR = 1
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S1.1hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 1: snr = 2
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000; snr = 2
true = c(rep(1, active), rep(0,p - active))
S1.2hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 1: SNR = 3
set.seed(1)
n = 20; p = 1000; active = 2; repeats = 1000;  snr = 3
true = c(rep(1, active), rep(0,p - active ))
S1.3hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))
#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# Setting 2: snr = 0.5
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 0.5
true = c(rep(1, active), rep(0,p - active ))
d = gendata(n = n, p = p, active = active)
S2.05hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 2: SNR = 1
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S2.1hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 2: snr = 2
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 2
true = c(rep(1, active), rep(0,p - active ))
S2.2hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 2: SNR = 3
set.seed(1)
n = 100; p = 500; active = 10; repeats = 1000; snr = 3
true = c(rep(1, active), rep(0,p - active ))
S2.3hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))
#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# Setting 3: snr = 0.5
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000;   snr = 0.5
true = c(rep(1, active), rep(0,p - active ))
d = gendata(n = n, p = p, active = active)
S3.05hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 3: SNR = 1
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S3.1hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 3: snr = 2
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000;  snr = 2
true = c(rep(1, active), rep(0,p - active))
S3.2hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 3: SNR = 3
set.seed(1)
n = 200; p = 200; active = 20; repeats = 1000;  snr = 3
true = c(rep(1, active), rep(0,p - active ))
S3.3hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))
#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# Setting 4: snr = 0.5
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 0.5
true = c(rep(1, active), rep(0,p - active ))
d = gendata(n = n, p = p, active = active)
S4.05hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 4: SNR = 1
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 1
true = c(rep(1, active), rep(0,p - active ))
S4.1hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 4: snr = 2
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 2
true = c(rep(1, active), rep(0,p - active ))
S4.2hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

# Setting 4: SNR = 3
set.seed(1)
n = 500; p = 100; active = 20; repeats = 1000; snr = 3
true = c(rep(1, active), rep(0,p - active ))
S4.3hard = pbreplicate(repeats, simulationATS(X = d$X, beta = d$beta, true = true, p = p, snr = snr))

save(S1.05hard, S1.2hard, S1.3hard, S1.1hard, file = paste0(Sys.Date(), "_S1.RData"))
save(S2.05hard, S2.2hard, S2.3hard, S2.1hard, file = paste0(Sys.Date(), "_S2.RData"))
save(S3.05hard, S3.2hard, S3.3hard ,S3.1hard, file = paste0(Sys.Date(), "_S3.RData"))
save(S4.05hard, S4.2hard, S4.3hard, S4.1hard, file = paste0(Sys.Date(), "_S4.RData"))

#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
# Figures
#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
# MCC
c11hard = combine(S1.05hard, S1.1hard, S1.2hard, S1.3hard,  1, filtered = filtered) |> makeCluster()
c22hard = combine(S2.05hard, S2.1hard, S2.2hard, S2.3hard,  2, filtered = filtered) |> makeCluster()
c33hard = combine(S3.05hard, S3.1hard, S3.2hard, S3.3hard,  3, filtered = filtered) |> makeCluster()
c44hard = combine(S4.05hard, S4.1hard, S4.2hard, S4.3hard,  4, filtered = filtered) |> makeCluster()

# c1hard = rbind(c11hard,c22hard) |> totplotnoaxis()
# c2hard = rbind(c33hard,c44hard) |> totplot()

c1hardM = rbind(c11hard,c22hard) |> totMEANplotnoaxis()
c2hardM = rbind(c33hard,c44hard) |> totMEANplot()

# c1hard/c2hard
# 14/10
c1hardM/c2hardM + plot_layout(guides = "collect", axis_titles = "collect")
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Variables Selected
c11Nhard = combineN(S1.05hard, S1.1hard, S1.2hard, S1.3hard,  1, filtered = filtered) |> makeCluster()
c22Nhard = combineN(S2.05hard, S2.1hard, S2.2hard, S2.3hard,  2, filtered = filtered) |> makeCluster()
c33Nhard = combineN(S3.05hard, S3.1hard, S3.2hard, S3.3hard,  3, filtered = filtered) |> makeCluster()
c44Nhard = combineN(S4.05hard, S4.1hard, S4.2hard, S4.3hard,  4, filtered = filtered) |> makeCluster()

# c1nhard = rbind(c11Nhard,c22Nhard) |> totplotNnoaxis(lim = 50)
# c2nhard = rbind(c33Nhard,c44Nhard) |> totplotN(lim = 60)

c1nhardM = rbind(c11Nhard,c22Nhard) |> NtotMEANplotnoaxis(lim = 50)
c2nhardM = rbind(c33Nhard,c44Nhard) |> NtotMEANplot(lim = 60)

# c1nhard/c2nhard
c1nhardM/c2nhardM + plot_layout(guides = "collect", axis_titles = "collect")
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Recall
c11Recallhard = combineRecall(S1.05hard,S1.1hard, S1.2hard, S1.3hard,  1, filtered = filtered) |> makeCluster()
c22Recallhard = combineRecall(S2.05hard,S2.1hard, S2.2hard, S2.3hard,  2, filtered = filtered) |> makeCluster()
c33Recallhard = combineRecall(S3.05hard,S3.1hard, S3.2hard, S3.3hard,  3, filtered = filtered) |> makeCluster()
c44Recallhard = combineRecall(S4.05hard,S4.1hard, S4.2hard, S4.3hard,  4, filtered = filtered) |> makeCluster()

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
#Precision 
c11Precisionhard = combinePrecision(S1.05hard, S1.1hard,S1.2hard, S1.3hard, 1, filtered = filtered) |> makeCluster()
c22Precisionhard = combinePrecision(S2.05hard, S2.1hard,S2.2hard, S2.3hard, 2, filtered = filtered) |> makeCluster()
c33Precisionhard = combinePrecision(S3.05hard, S3.1hard,S3.2hard, S3.3hard, 3, filtered = filtered) |> makeCluster()
c44Precisionhard = combinePrecision(S4.05hard, S4.1hard,S4.2hard, S4.3hard, 4, filtered = filtered) |> makeCluster()

phard1 = prplot1(c11Recallhard, c11Precisionhard, c22Recallhard, c22Precisionhard)
phard2 = prplot1(c33Recallhard, c33Precisionhard, c44Recallhard, c44Precisionhard)

phard1
phard2
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Exclusion probability threshold
s1.05e = exclusion(S1.05hard)
s1.2e =  exclusion(S1.2hard)
s1.3e =  exclusion(S1.3hard)
s1.1e =  exclusion(S1.1hard)

s2.05e = exclusion(S2.05hard)
s2.2e =  exclusion(S2.2hard)
s2.3e =  exclusion(S2.3hard)
s2.1e =  exclusion(S2.1hard)

s3.05e =  exclusion(S3.05hard)
s3.2e =   exclusion(S3.2hard)
s3.3e =   exclusion(S3.3hard)
s3.1e =   exclusion(S3.1hard)

s4.05e = exclusion(S4.05hard)
s4.2e =  exclusion(S4.2hard)
s4.3e =  exclusion(S4.3hard)
s4.1e =  exclusion(S4.1hard)

excl = data.frame("Eta" = 
                    c(s1.05e, s1.1e, s1.2e, s1.3e ,
                      s2.05e, s2.1e, s2.2e, s2.3e ,
                      s3.05e, s3.1e, s3.2e, s3.3e ,
                      s4.05e, s4.1e, s4.2e, s4.3e ),
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

S1pi = rbind(extractPi(S1.05hard, 0.5, 1),
             extractPi(S1.1hard, 1, 1),
             extractPi(S1.2hard, 2, 1),
             extractPi(S1.3hard, 3, 1))

S2pi = rbind(extractPi(S2.05hard, 0.5, 2),
             extractPi(S2.1hard, 1, 2),
             extractPi(S2.2hard, 2, 2),
             extractPi(S2.3hard, 3, 2))

S3pi = rbind(extractPi(S3.05hard, 0.5, 3),
             extractPi(S3.1hard, 1, 3),
             extractPi(S3.2hard, 2, 3),
             extractPi(S3.3hard, 3, 3))

S4pi = rbind(extractPi(S4.05hard, 0.5, 4),
             extractPi(S4.1hard, 1, 4),
             extractPi(S4.2hard, 2, 4),
             extractPi(S4.3hard, 3, 4))

Spi = rbind(S1pi, S2pi, S3pi, S4pi)


ggplot(Spi, aes(x = value, fill = variable)) + geom_density(alpha = 0.7) +
  facet_grid(dimension~SNR, labeller = label_parsed) + xlab(TeX("Estimated $\\pi$")) + ylab("Frequency") +
  theme_few_grid(base_size = 20) + scale_fill_manual(values = c("#FC8D62", "#FFD92F")) + labs(fill = "Method")
