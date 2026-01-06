# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code to reproduce Section 4.2                                   #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2026)                      #                                    
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
#                                                                   #
# Note: This will take 8 hours without parallelisation.             #
#       I have provided the Rdata file for my simulations.          #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("Functions.R")

# load("Data/pro4.RData")
# load("Data/pro9.RData")

pro = read.csv("Data/Proteomics.csv")
pro = pro |> select(-sampleID)
X = pro |> as.matrix() |> scale()

set.seed(1)
active = 9; repeats = 1000; snr = 1; p = ncol(X)
beta = c(rep(1,active), rep(0, ncol(X) - active)) 
true = c(rep(1, active),rep(0,ncol(X) - active))
pro9.1 = pbreplicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr, gaussian.knockoffs = T))

set.seed(1)
active = 9; repeats = 1000; snr = 3; p = ncol(X)
beta = c(rep(1,9), rep(0, ncol(X) - active)) 
true = c(rep(1, active),rep(0,ncol(X) - active))
pro9.3 = pbreplicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr, gaussian.knockoffs = T))

set.seed(1)
active = 4; repeats = 1000; snr = 1; p = ncol(X)
beta = c(rep(1,4), rep(0, ncol(X) - active)) 
true = c(rep(1, active),rep(0,ncol(X) - active))
pro4.1 = pbreplicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr, gaussian.knockoffs = T))

set.seed(1)
active = 4; repeats = 1000; snr = 3; p = ncol(X)
beta = c(rep(1,4), rep(0, ncol(X) - active)) 
true = c(rep(1, active),rep(0,ncol(X) - active))
pro4.3 = pbreplicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr, gaussian.knockoffs = T))

save(pro9.1, pro9.3, file = paste0(Sys.Date(), "_pro9.RData"))
save(pro4.1, pro4.3, file = paste0(Sys.Date(), "_pro4.RData"))

#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
# Figures
#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
# MCC
pro9.1.df = pro9.1 |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
pro9.3.df = pro9.3 |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
pro9 = rbind(pro9.1.df,pro9.3.df)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
pro9$Method = as.character(pro9$Method)
pro9$Method[pro9$Method == "LASSO 1SE"] = "LASSO"
pro9$Method[pro9$Method == "Exclusion ATS"] = "EATS"
pro9$Method[pro9$Method == "All"] = "ATS"
pro9$Method = factor(pro9$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
pro9 = data.frame(pro9) |> mutate(Dimension = "n==91*`,`~p==721*`,`~`|`*beta[S]*`|`==9") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster() |> filter(!Method %in% c("Static 0.60", "Static 0.90"))

pro4.1.df = pro4.1 |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
pro4.3.df = pro4.3 |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
pro4 = rbind(pro4.1.df,pro4.3.df)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
pro4$Method = as.character(pro4$Method)
pro4$Method[pro4$Method == "LASSO 1SE"] = "LASSO"
pro4$Method[pro4$Method == "Exclusion ATS"] = "EATS"
pro4$Method[pro4$Method == "All"] = "ATS"
pro4$Method = factor(pro4$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
pro4 = data.frame(pro4) |> mutate(Dimension = "n==91*`,`~p==721*`,`~`|`*beta[S]*`|`==4") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster() |> filter(!Method %in% c("Static 0.60", "Static 0.90"))

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Variables Selected
pro9.1.df.N = pro9.1 |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
pro9.3.df.N = pro9.3 |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
pro9N = rbind(pro9.1.df.N,pro9.3.df.N)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
pro9N$Method = as.character(pro9N$Method)
pro9N$Method[pro9N$Method == "LASSO 1SE"] = "LASSO"
pro9N$Method[pro9N$Method == "Exclusion ATS"] = "EATS"
pro9N$Method[pro9N$Method == "All"] = "ATS"
pro9N$Method = factor(pro9N$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
pro9N = data.frame(pro9N) |> mutate(Dimension = "n==91*`,`~p==721*`,`~`|`*beta[S]*`|`==9") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()|> filter(!Method %in% c("Static 0.60", "Static 0.90"))

pro4.1.df.N = pro4.1 |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
pro4.3.df.N = pro4.3 |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
pro4N = rbind(pro4.1.df.N,pro4.3.df.N)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
pro4N$Method = as.character(pro4N$Method)
pro4N$Method[pro4N$Method == "LASSO 1SE"] = "LASSO"
pro4N$Method[pro4N$Method == "Exclusion ATS"] = "EATS"
pro4N$Method[pro4N$Method == "All"] = "ATS"
pro4N$Method = factor(pro4N$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
pro4N = data.frame(pro4N) |> mutate(Dimension = "n==91*`,`~p==721*`,`~`|`*beta[S]*`|`==4") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()|> filter(!Method %in% c("Static 0.60", "Static 0.90"))

((totMEANplotnoaxis(pro9, proDiab = T) | totMEANplotnoaxis(pro4, proDiab = T)) +  plot_layout(axis_titles = "collect", guides = "collect")) / ((NtotMEANplot(pro9N, proDiab = T) + theme(legend.position = "none") | NtotMEANplot(pro4N, proDiab = T)) + theme(legend.position = "none") +  plot_layout(guides = "collect", axis_titles = "collect")) 
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# ATS/EATS Estimated Pi
proPi = rbind(extractPi(pro4.1, snr = 1, data = "proteomics", active = 4),
              extractPi(pro4.3, snr = 3, data = "proteomics", active = 4),
              extractPi(pro9.1, snr = 1, data = "proteomics", active = 9),
              extractPi(pro9.3, snr = 3, data = "proteomics", active = 9))

ggplot(proPi, aes(x = value, fill = variable)) + geom_density(alpha = 0.7)  +
  facet_grid(dimension~SNR, labeller = label_parsed) + xlab(TeX("Estimated $\\pi$")) + ylab("Frequency") +
  theme_few_grid(base_size = 20) + scale_fill_manual(values = c("#FC8D62", "#FFD92F")) + labs(fill = "Method")
