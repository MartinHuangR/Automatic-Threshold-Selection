# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code to reproduce Section 4.2                                   #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2025)                      #                                    
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
#                                                                   #
# Note: This will take 3 hours without parallelisation.             #
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
pro9 = data.frame(pro9) |> mutate(Dimension = "`|`*beta[S]*`|`==9") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

pro4.1.df = pro4.1 |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
pro4.3.df = pro4.3 |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
pro4 = rbind(pro4.1.df,pro4.3.df)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
pro4$Method = as.character(pro9$Method)
pro4$Method[pro4$Method == "LASSO 1SE"] = "LASSO"
pro4$Method[pro4$Method == "Exclusion ATS"] = "EATS"
pro4$Method[pro4$Method == "All"] = "ATS"
pro4$Method = factor(pro4$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
pro4 = data.frame(pro4) |> mutate(Dimension = "`|`*beta[S]*`|`==4") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

pro9plot = ggplot(pro9,aes(x = Method, y = MCC, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_grid(Dimension ~SNR, labeller = label_parsed) +
  theme_few_grid(base_size = 20) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  ylim(0, 1) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x= element_blank()) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

pro4plot = ggplot(pro4,aes(x = Method, y = MCC, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_grid(Dimension ~SNR, labeller = label_parsed) +
  theme_few_grid(base_size = 20) +
  ylim(0, 1) + 
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x= element_blank(),
        strip.text.x = element_blank())+
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

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
pro9N = data.frame(pro9N) |> mutate(Dimension = "`|`*beta[S]*`|`==9") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

pro4.1.df.N = pro4.1 |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
pro4.3.df.N = pro4.3 |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
pro4N = rbind(pro4.1.df.N,pro4.3.df.N)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
pro4N$Method = as.character(pro9N$Method)
pro4N$Method[pro4N$Method == "LASSO 1SE"] = "LASSO"
pro4N$Method[pro4N$Method == "Exclusion ATS"] = "EATS"
pro4N$Method[pro4N$Method == "All"] = "ATS"
pro4N$Method = factor(pro4N$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
pro4N = data.frame(pro4N) |> mutate(Dimension = "`|`*beta[S]*`|`==4") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

pro9nplot = ggplot(pro9N,aes(x = Method, y = NN, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_grid(Dimension ~SNR, labeller = label_parsed) +
  ylab("Variables Selected") + 
  theme_few_grid(base_size = 20) +
  coord_cartesian(ylim = c(0, 50)) + 
  geom_hline( aes(yintercept = 9), linetype = "dashed", linewidth = 0.5) +
  xlab(element_blank()) + 
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x= element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

pro4nplot = ggplot(pro4N,aes(x = Method, y = NN, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_grid(Dimension ~SNR, labeller = label_parsed) +
  ylab("Variables Selected") + 
  theme_few_grid(base_size = 20) +
  coord_cartesian(ylim = c(0, 60)) + 
  geom_hline( aes(yintercept = 4), linetype = "dashed", linewidth = 0.5) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  xlab(element_blank()) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text.x = element_blank()) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

pro9plot/pro4plot/pro9nplot/pro4nplot

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
# Precision
pro4.1.df.Precision = pro4.1 |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
pro4.3.df.Precision = pro4.3 |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
pro4Precision = rbind(pro4.1.df.Precision,pro4.3.df.Precision)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
pro4Precision$Method = as.character(pro4Precision$Method)
pro4Precision$Method[pro4Precision$Method == "LASSO 1SE"] = "LASSO"
pro4Precision$Method[pro4Precision$Method == "Exclusion ATS"] = "EATS"
pro4Precision$Method[pro4Precision$Method == "All"] = "ATS"
pro4Precision$Method = factor(pro4Precision$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
pro4Precision = data.frame(pro4Precision) |> mutate(Dimension = "`|`*beta[S]*`|`==4") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

pro9.1.df.Precision = pro9.1 |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
pro9.3.df.Precision = pro9.3 |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
pro9Precision = rbind(pro9.1.df.Precision,pro9.3.df.Precision)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
pro9Precision$Method = as.character(pro4Precision$Method)
pro9Precision$Method[pro9Precision$Method == "LASSO 1SE"] = "LASSO"
pro9Precision$Method[pro9Precision$Method == "Exclusion ATS"] = "EATS"
pro9Precision$Method[pro9Precision$Method == "All"] = "ATS"
pro9Precision$Method = factor(pro9Precision$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
pro9Precision = data.frame(pro9Precision) |> mutate(Dimension = "`|`*beta[S]*`|`==9") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

pro4nplotprecision = ggplot(pro4Precision,aes(x = Method, y = Precision, fill = Method)) + geom_boxplot(alpha = 0.65) +
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

pro9nplotprecision = ggplot(pro9Precision,aes(x = Method, y = Precision, fill = Method)) + geom_boxplot(alpha = 0.65) +
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
pro4.1.df.Recall = pro4.1 |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
pro4.3.df.Recall = pro4.3 |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
pro4Recall = rbind(pro4.1.df.Recall,pro4.3.df.Recall)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
pro4Recall$Method = as.character(pro4Recall$Method)
pro4Recall$Method[pro4Recall$Method == "LASSO 1SE"] = "LASSO"
pro4Recall$Method[pro4Recall$Method == "Exclusion ATS"] = "EATS"
pro4Recall$Method[pro4Recall$Method == "All"] = "ATS"
pro4Recall$Method = factor(pro4Recall$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
pro4Recall = data.frame(pro4Recall) |> mutate(Dimension = "`|`*beta[S]*`|`==4") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

pro9.1.df.Recall = pro9.1 |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==1")
pro9.3.df.Recall = pro9.3 |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==3")
pro9Recall = rbind(pro9.1.df.Recall,pro9.3.df.Recall)|> mutate(SNR = factor(SNR, levels = c("~SNR==1", "~SNR==3")))
pro9Recall$Method = as.character(pro9Recall$Method)
pro9Recall$Method[pro9Recall$Method == "LASSO 1SE"] = "LASSO"
pro9Recall$Method[pro9Recall$Method == "Exclusion ATS"] = "EATS"
pro9Recall$Method[pro9Recall$Method == "All"] = "ATS"
pro9Recall$Method = factor(pro9Recall$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
pro9Recall = data.frame(pro9Recall) |> mutate(Dimension = "`|`*beta[S]*`|`==9") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

pro4nplotRecall = ggplot(pro4Recall,aes(x = Method, y = Recall, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_wrap(Dimension ~SNR, labeller = label_parsed, ncol = 4) +
  ylab("Recall") + 
  theme_few_grid(base_size = 20) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

pro9nplotRecall = ggplot(pro9Recall,aes(x = Method, y = Recall, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_wrap(Dimension ~SNR, labeller = label_parsed, ncol = 4) +
  ylab("Recall") + 
  theme_few_grid(base_size = 20) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

pr = rbind(pro9Recall, pro4Recall)
pp = rbind(pro9Precision, pro4Precision)

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