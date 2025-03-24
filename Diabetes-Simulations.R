# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code to reproduce Section 4.3                                   #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2025)                      #                                    
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
  mutate(Dimension = factor(Dimension)) |> makeCluster()

d10.1.df = d10.1 |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
d10.3.df = d10.3 |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
d10 = rbind(d10.1.df,d10.3.df)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
d10$Method = as.character(d5$Method)
d10$Method[d10$Method == "LASSO 1SE"] = "LASSO"
d10$Method[d10$Method == "Exclusion ATS"] = "EATS"
d10$Method = factor(d10$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
d10 = data.frame(d10) |> mutate(Dimension = "n==442*`,`~p==64*`,`~`|`*beta[S]*`|`==10") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

d5plot = ggplot(d5,aes(x = Method, y = MCC, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_grid(Dimension ~SNR, labeller = label_parsed) +
  theme_few_grid(base_size = 20) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  xlab(element_blank()) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x= element_blank(),
        strip.text.x = element_blank()) + 
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

d10plot = ggplot(d10,aes(x = Method, y = MCC, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_grid(Dimension ~SNR, labeller = label_parsed) +
  theme_few_grid(base_size = 20) +
  xlab(element_blank()) + 
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x= element_blank()) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

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
  mutate(Dimension = factor(Dimension)) |> makeCluster()

d10.1.df.N = d10.1 |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
d10.3.df.N = d10.3 |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
d10N = rbind(d10.1.df.N,d10.3.df.N)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
d10N$Method = as.character(d5N$Method)
d10N$Method[d10N$Method == "LASSO 1SE"] = "LASSO"
d10N$Method[d10N$Method == "Exclusion ATS"] = "EATS"
d10N$Method[d10N$Method == "All"] = "ATS"
d10N$Method = factor(d10N$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
d10N = data.frame(d10N) |> mutate(Dimension = "n==442*`,`~p==64*`,`~`|`*beta[S]*`|`==10") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

d5nplot = ggplot(d5N,aes(x = Method, y = NN, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_grid(Dimension ~SNR, labeller = label_parsed) +
  ylab("Variables Selected") + 
  theme_few_grid(base_size = 20) +
  geom_hline( aes(yintercept = 5), linetype = "dashed", linewidth = 0.5) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  coord_cartesian(ylim = c(0, 25)) +
  xlab(element_blank()) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text.x = element_blank()) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

d10nplot = ggplot(d10N,aes(x = Method, y = NN, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_grid(Dimension ~SNR, labeller = label_parsed) +
  ylab("Variables Selected") + 
  theme_few_grid(base_size = 20) +
  coord_cartesian(ylim = c(0, 35)) +
  xlab(element_blank()) + 
  geom_hline( aes(yintercept = 10), linetype = "dashed", linewidth = 0.5) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x= element_blank(),
        strip.text.x = element_blank())   +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

d10plot/d5plot/d10nplot/d5nplot
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Precision
d5.1.df.Precision = d5.1 |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
d5.3.df.Precision = d5.3 |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
d5Precision = rbind(d5.1.df.Precision,d5.3.df.Precision)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
d5Precision$Method = as.character(d5Precision$Method)
d5Precision$Method[d5Precision$Method == "LASSO 1SE"] = "LASSO"
d5Precision$Method[d5Precision$Method == "Exclusion ATS"] = "EATS"
d5Precision$Method[d5Precision$Method == "All"] = "ATS"
d5Precision$Method = factor(d5Precision$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
d5Precision = data.frame(d5Precision) |> mutate(Dimension = "n==442*`,`~p==64*`,`~`|`*beta[S]*`|`==5") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

d10.1.df.Precision = d10.1 |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
d10.3.df.Precision = d10.3 |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
d10Precision = rbind(d10.1.df.Precision,d10.3.df.Precision)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
d10Precision$Method = as.character(d10Precision$Method)
d10Precision$Method[d10Precision$Method == "LASSO 1SE"] = "LASSO"
d10Precision$Method[d10Precision$Method == "Exclusion ATS"] = "EATS"
d10Precision$Method[d10Precision$Method == "All"] = "ATS"
d10Precision$Method = factor(d10Precision$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
d10Precision = data.frame(d10Precision) |> mutate(Dimension = "n==442*`,`~p==64*`,`~`|`*beta[S]*`|`==10") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

d5nplotprecision = ggplot(d5Precision,aes(x = Method, y = Precision, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_wrap(Dimension ~SNR, labeller = label_parsed, ncol = 4) +
  ylab("Precision") + 
  theme_few_grid(base_size = 20) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x= element_blank()) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

d10nplotprecision = ggplot(d10Precision,aes(x = Method, y = Precision, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_wrap(Dimension ~SNR, labeller = label_parsed, ncol = 4) +
  ylab("Precision") + 
  theme_few_grid(base_size = 20) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x= element_blank()) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Recall
d5.1.df.Recall = d5.1 |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
d5.3.df.Recall = d5.3 |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
d5Recall = rbind(d5.1.df.Recall,d5.3.df.Recall)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
d5Recall$Method = as.character(d5Recall$Method)
d5Recall$Method[d5Recall$Method == "LASSO 1SE"] = "LASSO"
d5Recall$Method[d5Recall$Method == "Exclusion ATS"] = "EATS"
d5Recall$Method[d5Recall$Method == "All"] = "ATS"
d5Recall$Method = factor(d5Recall$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
d5Recall = data.frame(d5Recall) |> mutate(Dimension = "n==442*`,`~p==64*`,`~`|`*beta[S]*`|`==5") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

d10.1.df.Recall = d10.1 |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
d10.3.df.Recall = d10.3 |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
d10Recall = rbind(d10.1.df.Recall,d10.3.df.Recall)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
d10Recall$Method = as.character(d5Recall$Method)
d10Recall$Method[d10Recall$Method == "LASSO 1SE"] = "LASSO"
d10Recall$Method[d10Recall$Method == "Exclusion ATS"] = "EATS"
d10Recall$Method[d10Recall$Method == "All"] = "ATS"
d10Recall$Method = factor(d10Recall$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
d10Recall = data.frame(d10Recall) |> mutate(Dimension = "n==442*`,`~p==64*`,`~`|`*beta[S]*`|`==10") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

pr = rbind(d10Recall, d5Recall)
pp = rbind(d10Precision, d5Precision)
prplot1 = ggplot(pr,aes(x = Method, y = Recall, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_grid(Dimension~SNR, labeller = label_parsed) +
  ylab("Recall") + 
  theme_few_grid(base_size = 20) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text.y = element_blank()) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

prplot2 = ggplot(pp,aes(x = Method, y = Precision, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_grid(Dimension~SNR, labeller = label_parsed) +
  ylab("Precision") + 
  theme_few_grid(base_size = 20) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

prplot1|prplot2
#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# ATS/EATS Estimated Pi
diaPi = rbind(extractPi(d5.1, snr = 1, data = "diabetes", active = 5),
              extractPi(d5.3, snr = 3, data = "diabetes", active = 5),
              extractPi(d10.1, snr = 1, data = "diabetes", active = 10),
              extractPi(d10.3, snr = 3, data = "diabetes", active = 10))

ggplot(diaPi, aes(x = value, fill = variable)) + geom_density(alpha = 0.7) +
  facet_grid(dimension~SNR, labeller = label_parsed) + xlab(TeX("Estimated $\\pi$")) + ylab("Frequency") +
  theme_few_grid(base_size = 20) + scale_fill_manual(values = c("#FC8D62", "#FFD92F")) + labs(fill = "Method")
