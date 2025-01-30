# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code to reproduce Section 4 FDR Study with Artificial Data      #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2025)                      #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
#                                                                   #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("Functions.R")
load("RData/FDR.I.RData")
load("RData/FDR.II.RData")
load("RData/FDR.III.RData")
load("RData/FDR.IV.RData")

# Setting 1
set.seed(1)
n = 20; p = 1000; active = 2; snr = 10
d = gendata(n = n, p = p, active = active)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 10
fIs10e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 5
fIs5e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 10
d = gendata(n = n, p = p, active = active)
fIs10e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 5
fIs5e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 10
fIs10e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)

set.seed(1)
n = 20; p = 1000; active = 2; snr = 5
fIs5e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)

# Setting 2 

set.seed(1)
n = 100; p = 500; active = 10; snr = 10
d = gendata(n = n, p = p, active = active)

set.seed(1)
n = 100; p = 500; active = 10; snr = 10
fIIs10e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 100; p = 500; active = 10; snr = 5
fIIs5e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 100; p = 500; active = 10; snr = 10
d = gendata(n = n, p = p, active = active)
fIIs10e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 100; p = 500; active = 10; snr = 5
fIIs5e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 100; p = 500; active = 10; snr = 10
fIIs10e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)

set.seed(1)
n = 100; p = 500; active = 10; snr = 5
fIIs5e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)


# Setting 3
set.seed(1)
n = 200; p = 200; active = 20; snr = 10
d = gendata(n = n, p = p, active = active)

set.seed(1)
n = 200; p = 200; active = 20; snr = 10
fIIIs10e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 200; p = 200; active = 20; snr = 5
fIIIs5e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 200; p = 200; active = 20; snr = 10
d = gendata(n = n, p = p, active = active)
fIIIs10e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 200; p = 200; active = 20; snr = 5
fIIIs5e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 200; p = 200; active = 20; snr = 10
fIIIs10e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)

set.seed(1)
n = 200; p = 200; active = 20; snr = 5
fIIIs5e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)

# Setting 4
set.seed(1)
n = 500; p = 100; active = 10; snr = 10
d = gendata(n = n, p = p, active = active)

set.seed(1)
n = 500; p = 100; active = 10; snr = 10
fIVs10e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 500; p = 100; active = 10; snr = 5
fIVs5e2 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 2, LOOPS = 1000)

set.seed(1)
n = 500; p = 100; active = 10; snr = 10
d = gendata(n = n, p = p, active = active)
fIVs10e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 500; p = 100; active = 10; snr = 5
fIVs5e10 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 10, LOOPS = 1000)

set.seed(1)
n = 500; p = 100; active = 10; snr = 10
d = gendata(n = n, p = p, active = active)
fIVs10e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)

set.seed(1)
n = 500; p = 100; active = 10; snr = 5
fIVs5e5 = fdr(X = d$X, beta = d$beta, p = p, snr = snr, EV = 5, LOOPS = 1000)

# Table Output
fI = rbind(fIs10e2, fIs5e2, fIs5e5, fIs10e5, fIs5e10, fIs10e10)  |> mutate(Setting = "I")
fII = rbind(fIIs10e2, fIIs5e2, fIIs5e5, fIIs10e5, fIIs5e10, fIIs10e10) |> mutate(Setting = "II")
fIII = rbind(fIIIs10e2, fIIIs5e2, fIIIs5e5, fIIIs10e5, fIIIs5e10, fIIIs10e10) |> mutate(Setting = "III")
fIV = rbind(fIVs10e2, fIVs5e2, fIVs5e5, fIVs10e5, fIVs5e10, fIVs10e10) |> mutate(Setting = "IV")

f = rbind(fI, fII, fIII, fIV) |> 
  group_by(Setting, snr, EV) |> 
  summarise(error = 1 - mean(false.selections > EV),
            correct.prop = mean(correct.selections.prop),
            avg.selected = mean(n.selected)) |> 
  mutate(correct.prop = round(correct.prop, 2))

