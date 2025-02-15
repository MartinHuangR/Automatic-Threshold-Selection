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
load("Data/S1.Rdata") 
load("Data/S2.Rdata")
load("Data/S3.Rdata")
load("Data/S4.Rdata")


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

#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
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




c11Nhard = combineN(S1.10hard, S1.5hard, S1.3hard, S1.1hard, 1, filtered = filtered) |> makeCluster()
c22Nhard = combineN(S2.10hard, S2.5hard, S2.3hard, S2.1hard, 2, filtered = filtered) |> makeCluster()
c33Nhard = combineN(S3.10hard, S3.5hard, S3.3hard, S3.1hard, 3, filtered = filtered) |> makeCluster()
c44Nhard = combineN(S4.10hard, S4.5hard, S4.3hard, S4.1hard, 4, filtered = filtered) |> makeCluster()

c1nhard = rbind(c11Nhard,c22Nhard) |> totplotNnoaxis(lim = 40)
c2nhard = rbind(c33Nhard,c44Nhard) |> totplotN(lim = 60)

library(patchwork)
c1nhard/c2nhard


# Recall
c11Recallhard = combineRecall(S1.10hard, S1.5hard, S1.3hard, S1.1hard, 1, filtered = filtered) |> makeCluster()
c22Recallhard = combineRecall(S4.10hard, S4.5hard, S4.3hard, S4.1hard, 2, filtered = filtered) |> makeCluster()
c33Recallhard = combineRecall(S2.10hard, S2.5hard, S2.3hard, S2.1hard, 3, filtered = filtered) |> makeCluster()
c44Recallhard = combineRecall(S3.10hard, S3.5hard, S3.3hard, S3.1hard, 4, filtered = filtered) |> makeCluster()


#Precision 
c11Precisionhard = combinePrecision(S1.10hard, S1.5hard, S1.3hard, S1.1hard, 1, filtered = filtered) |> makeCluster()
c22Precisionhard = combinePrecision(S4.10hard, S4.5hard, S4.3hard, S4.1hard, 2, filtered = filtered) |> makeCluster()
c33Precisionhard = combinePrecision(S2.10hard, S2.5hard, S2.3hard, S2.1hard, 3, filtered = filtered) |> makeCluster()
c44Precisionhard = combinePrecision(S3.10hard, S3.5hard, S3.3hard, S3.1hard, 4, filtered = filtered) |> makeCluster()


phard1 = prplot1(c11Recallhard, c11Precisionhard, c22Recallhard, c22Precisionhard)
phard2 = prplot1(c33Recallhard, c33Precisionhard, c44Recallhard, c44Precisionhard)

phard1
phard2


# Exclusion probability threshold
s1.10e = exclusion(S1.10hard)
s1.5e =  exclusion(S1.5hard)
s1.3e =  exclusion(S1.3hard)
s1.1e =  exclusion(S1.1hard)

s2.10e = exclusion(S2.10hard)
s2.5e =  exclusion(S2.5hard)
s2.3e =  exclusion(S2.3hard)
s2.1e =  exclusion(S2.1hard)

s3.10e =  exclusion(S3.10hard)
s3.5e =   exclusion(S3.5hard)
s3.3e =   exclusion(S3.3hard)
s3.1e =   exclusion(S3.1hard)

s4.10e = exclusion(S4.10hard)
s4.5e =  exclusion(S4.5hard)
s4.3e =  exclusion(S4.3hard)
s4.1e =  exclusion(S4.1hard)

excl = data.frame("Eta" = 
                    c(s1.10e, s1.5e, s1.3e, s1.1e,
                      s2.10e, s2.5e, s2.3e, s2.1e,
                      s3.10e, s3.5e, s3.3e, s3.1e,
                      s4.10e, s4.5e, s4.3e, s4.1e),
                  "Setting" = c(rep("(I):~n==20*`,`~p==1000*`,`~`|`*beta[S]*`|`==2",4000),
                                rep("(II):~n==100*`,`~p==500*`,`~`|`*beta[S]*`|`==10",4000),
                                rep("(III):~n==200*`,`~p==200*`,`~`|`*beta[S]*`|`==20",4000),
                                rep("(IV):~n==500*`,`~p==100*`,`~`|`*beta[S]*`|`==20",4000)),
                  "SNR" = rep(c(rep(10, 1000), rep(5,1000), rep(3, 1000), rep(1,1000)), 4)
) |> mutate(SNR = factor(SNR))

ggplot(excl, aes(x = SNR, y = Eta)) + geom_boxplot(alpha = 0.65) +
  facet_wrap(~Setting, labeller = label_parsed, ncol = 4) +
  theme_few_grid(base_size = 20) +
  ylab(TeX("Exclusion Probability Threshold $\\eta$")) + 
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  theme(legend.position = "none") 
