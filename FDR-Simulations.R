# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code to reproduce Section 4.4                                   #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2025)                      #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("Functions.R")

load("Data/FDR.I.RData")
load("Data/FDR.II.RData")
load("Data/FDR.III.RData")
load("Data/FDR.IV.RData")

# Setting 1
set.seed(1)
n = 20; p = 1000; active = 2; 
d = gendata(n = n, p = p, active = active)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 1
fIs1e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 3
fIs3e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 1
d = gendata(n = n, p = p, active = active)
fIs1e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 3
fIs3e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 1
fIs1e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 3
fIs3e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)

# Setting 2 

set.seed(1)
n = 100; p = 500; active = 10; 
d = gendata(n = n, p = p, active = active)

set.seed(1)
n = 100; p = 500; active = 10; snr = 1
fIIs1e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 100; p = 500; active = 10; snr = 3
fIIs3e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 100; p = 500; active = 10; snr = 1
d = gendata(n = n, p = p, active = active)
fIIs1e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 100; p = 500; active = 10; snr = 3
fIIs3e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 100; p = 500; active = 10; snr = 1
fIIs1e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)

set.seed(1)
n = 100; p = 500; active = 10; snr = 3
fIIs3e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)


# Setting 3
set.seed(1)
n = 200; p = 200; active = 20; 
d = gendata(n = n, p = p, active = active)

set.seed(1)
n = 200; p = 200; active = 20; snr = 1
fIIIs1e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 200; p = 200; active = 20; snr = 3
fIIIs3e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 200; p = 200; active = 20; snr = 1
d = gendata(n = n, p = p, active = active)
fIIIs1e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 200; p = 200; active = 20; snr = 3
fIIIs3e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 200; p = 200; active = 20; snr = 1
fIIIs1e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)

set.seed(1)
n = 200; p = 200; active = 20; snr = 3
fIIIs3e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)

# Setting 4
set.seed(1)
n = 500; p = 100; active = 20; snr = 1
d = gendata(n = n, p = p, active = active)

set.seed(1)
n = 500; p = 100; active = 20; snr = 1
fIVs1e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 500; p = 100; active = 20; snr = 3
fIVs3e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 500; p = 100; active = 20; snr = 1
d = gendata(n = n, p = p, active = active)
fIVs1e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 500; p = 100; active = 20; snr = 3
fIVs3e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 500; p = 100; active = 20; snr = 1
d = gendata(n = n, p = p, active = active)
fIVs1e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)

set.seed(1)
n = 500; p = 100; active = 20; snr = 3
fIVs3e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)

#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
# Table
#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
fI = rbind(fIs1e2, fIs3e2, fIs1e5, fIs3e5, fIs1e10, fIs3e10)  |> mutate(Setting = "I")
fII = rbind(fIIs1e2, fIIs3e2, fIIs1e5, fIIs3e5, fIIs1e10, fIIs3e10) |> mutate(Setting = "II")
fIII = rbind(fIIIs1e2, fIIIs3e2, fIIIs1e5, fIIIs3e5, fIIIs1e10, fIIIs3e10) |> mutate(Setting = "III")
fIV = rbind(fIVs1e2, fIVs3e2, fIVs1e5, fIVs3e5, fIVs1e10, fIVs3e10) |> mutate(Setting = "IV")

f = rbind(fI, fII, fIII, fIV) |> 
  group_by(Setting, snr, EV) |> 
  summarise(ebs = 1 - mean(false.selections > EV),
            correct.prop = mean(correct.selections.prop),
            avg.selected = mean(n.selected)) |> 
  mutate(correct.prop = round(correct.prop, 2),
         avg.selected = round(avg.selected,2))
