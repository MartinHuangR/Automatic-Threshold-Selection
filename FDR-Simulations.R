# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code to reproduce Section 4 FDR Study with Artificial Data      #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2025)                      #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
#                                                                   #
# Note: This will take around 24 hours to run. I have provided      #
# an Rdata file for my simulations.                                 #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("Functions.R")
load("RData/fdrI.RData")
load("RData/fdrII.RData")
load("RData/fdrIV.RData")
load("RData/fdrIII.RData")


# Setting I
set.seed(1)
n = 20; p = 1000; active = 2; snr = 10
d = gendata(n = n, p = p, active = active)
fIs10e1 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 1, LOOPS = 1000)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 5
d = gendata(n = n, p = p, active = active)
fIs5e1 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 1, LOOPS = 1000)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 10
d = gendata(n = n, p = p, active = active)
fIs10e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 5
d = gendata(n = n, p = p, active = active)
fIs5e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

# Setting II
set.seed(1)
n = 100; p = 500; active = 10; snr = 10
d = gendata(n = n, p = p, active = active)
fIIs10e1 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 1, LOOPS = 1000)

set.seed(1)
n = 100; p = 500; active = 10; snr = 5
d = gendata(n = n, p = p, active = active)
fIIs5e1 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 1, LOOPS = 1000)

set.seed(1)
n = 100; p = 500; active = 10; snr = 10
d = gendata(n = n, p = p, active = active)
fIIs10e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 100; p = 500; active = 10; snr = 5
d = gendata(n = n, p = p, active = active)
fIIs5e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

# Setting III
set.seed(1)
n = 200; p = 200; active = 20; snr = 10
d = gendata(n = n, p = p, active = active)
fIIIs10e1 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 1, LOOPS = 1000)

set.seed(1)
n = 200; p = 200; active = 20; snr = 5
d = gendata(n = n, p = p, active = active)
fIIIs5e1 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 1, LOOPS = 1000)

set.seed(1)
n = 200; p = 200; active = 20; snr = 10
d = gendata(n = n, p = p, active = active)
fIIIs10e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 200; p = 200; active = 20; snr = 5
d = gendata(n = n, p = p, active = active)
fIIIs5e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

# Setting IV
set.seed(1)
n = 500; p = 100; active = 10; snr = 10
d = gendata(n = n, p = p, active = active)
fIVs10e1 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 1, LOOPS = 1000)

set.seed(1)
n = 500; p = 100; active = 10; snr = 5
d = gendata(n = n, p = p, active = active)
fIVs5e1 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 1, LOOPS = 1000)

set.seed(1)
n = 500; p = 100; active = 10; snr = 10
d = gendata(n = n, p = p, active = active)
fIVs10e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 500; p = 100; active = 10; snr = 5
d = gendata(n = n, p = p, active = active)
fIVs5e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)



# Table results

fIII = rbind(fIIIs10e1, fIIIs5e1, fIIIs10e2, fIIIs5e2)
fIV = rbind(fIVs10e1,fIVs5e1,fIVs10e2, fIVs5e2)
fII = rbind(fIIs10e1,fIIs5e1,fIIs10e2, fIIs5e2)
fI = rbind(fIs10e1,fIs5e1,fIs10e2, fIs5e2)


fI = fI |> mutate(Setting = "I")
fII = fII |> mutate(Setting = "II")
fIII = fIII |> mutate(Setting = "III")
fIV = fIV |> mutate(Setting = "IV")

fdf = rbind(fI, fII, fIII, fIV)
fdf |> group_by(Setting, EV, snr) |> summarise(error = sum(false.selections > EV))