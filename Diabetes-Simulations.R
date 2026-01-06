# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code to reproduce Section 4.3                                   #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2026)                      #                                    
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #
#                                                                   #
# Note: This will take 4 hours without parallelisation.             #
#       I have provided the Rdata file for my simulations.          #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("Functions.R")
# load("Data/d10.Rdata")
# load("Data/d5.Rdata")

data("diabetes",package = "lars")
X = cbind(diabetes$x2) |> data.frame() |> scale()
y = diabetes$y

set.seed(1)
active = 5; repeats = 1000;  snr = 1; p = 64
beta = c(rep(1,active), rep(0, ncol(X) - active))
true = c(rep(1, active),rep(0,ncol(X) - active ))
d5.1 = pbreplicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr))

set.seed(1)
active = 5; repeats = 1000;  snr = 3; p = 64
beta = c(rep(1,active), rep(0, ncol(X) - active))
true = c(rep(1, active),rep(0,ncol(X) - active ))
d5.3 = pbreplicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr))


set.seed(1)
active = 10; repeats = 1000;  snr = 1; p = 64
beta = c(rep(1,active), rep(0, ncol(X) - active))
true = c(rep(1, active),rep(0,ncol(X) - active ))
d10.1 = pbreplicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr))

set.seed(1)
active = 10; repeats = 1000;  snr = 3; p = 64
beta = c(rep(1,active), rep(0, ncol(X) - active))
true = c(rep(1, active),rep(0,ncol(X) - active ))
d10.3 = pbreplicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr))

save(d5.1, d5.3, file = paste0(Sys.Date(), "_d5.RData"))
save(d10.1, d10.3, file = paste0(Sys.Date(), "_d10.RData"))

#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
# Figures
#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
# MCC
d5.1.df = d5.1 |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
d5.3.df = d5.3 |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
d5 = rbind(d5.1.df,d5.3.df) |> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
d5$Method = as.character(d5$Method)
d5$Method[d5$Method == "LASSO 1SE"] = "LASSO"
d5$Method[d5$Method == "Exclusion ATS"] = "EATS"
d5$Method = factor(d5$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
d5 = data.frame(d5) |> mutate(Dimension = "n==442*`,`~p==64*`,`~`|`*beta[S]*`|`==5") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()|> filter(!Method %in% c("Static 0.60", "Static 0.90"))

d10.1.df = d10.1 |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
d10.3.df = d10.3 |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
d10 = rbind(d10.1.df,d10.3.df)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
d10$Method = as.character(d10$Method)
d10$Method[d10$Method == "LASSO 1SE"] = "LASSO"
d10$Method[d10$Method == "Exclusion ATS"] = "EATS"
d10$Method = factor(d10$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
d10 = data.frame(d10) |> mutate(Dimension = "n==442*`,`~p==64*`,`~`|`*beta[S]*`|`==10") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()|> filter(!Method %in% c("Static 0.60", "Static 0.90"))
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Variables Selected
d5.1.df.N = d5.1 |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
d5.3.df.N = d5.3 |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
d5N = rbind(d5.1.df.N,d5.3.df.N)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
d5N$Method = as.character(d5N$Method)
d5N$Method[d5N$Method == "LASSO 1SE"] = "LASSO"
d5N$Method[d5N$Method == "Exclusion ATS"] = "EATS"
d5N$Method[d5N$Method == "All"] = "ATS"
d5N$Method = factor(d5N$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
d5N = data.frame(d5N) |> mutate(Dimension = "n==442*`,`~p==64*`,`~`|`*beta[S]*`|`==5") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster() |> filter(!Method %in% c("Static 0.60", "Static 0.90"))

d10.1.df.N = d10.1 |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
d10.3.df.N = d10.3 |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
d10N = rbind(d10.1.df.N,d10.3.df.N)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
d10N$Method = as.character(d10N$Method)
d10N$Method[d10N$Method == "LASSO 1SE"] = "LASSO"
d10N$Method[d10N$Method == "Exclusion ATS"] = "EATS"
d10N$Method[d10N$Method == "All"] = "ATS"
d10N$Method = factor(d10N$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
d10N = data.frame(d10N) |> mutate(Dimension = "n==442*`,`~p==64*`,`~`|`*beta[S]*`|`==10") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()|> filter(!Method %in% c("Static 0.60", "Static 0.90"))


# 14/10
((totMEANplotnoaxis(d10, proDiab = T) | totMEANplotnoaxis(d5,proDiab = T)) +  plot_layout(axis_titles = "collect", guides = "collect")) / ((NtotMEANplot(d10N,proDiab = T) + theme(legend.position = "none") | NtotMEANplot(d5N,proDiab = T)) + theme(legend.position = "none") +  plot_layout(guides = "collect", axis_titles = "collect")) 
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# ATS/EATS Estimated Pi
diaPi = rbind(extractPi(d5.1, snr = 1, data = "diabetes", active = 5),
              extractPi(d5.3, snr = 3, data = "diabetes", active = 5),
              extractPi(d10.1, snr = 1, data = "diabetes", active = 10),
              extractPi(d10.3, snr = 3, data = "diabetes", active = 10))

ggplot(diaPi, aes(x = value, fill = variable)) + geom_density(alpha = 0.7) +
  facet_grid(dimension~SNR, labeller = label_parsed) + xlab(TeX("Estimated $\\pi$")) + ylab("Frequency") +
  theme_few_grid(base_size = 20) + scale_fill_manual(values = c("#FC8D62", "#FFD92F")) + labs(fill = "Method")
