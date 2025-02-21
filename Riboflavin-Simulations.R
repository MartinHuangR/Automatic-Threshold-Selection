# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code to reproduce Section 4                                     #
# in the paper: Data-Adaptive and Automatic Stable Threshold        #
# Calibration for Stability Selection in Penalised Regression       #
#                       (Huang et al. 2024)                         #                                     
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
#                                                                   #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


source("Functions.R")
data(riboflavin)
X = data.matrix(riboflavin$x)
y = riboflavin$y
yxdf <- as.data.frame(cbind(y, X))
names(yxdf)[1] <- "y"
X = yxdf[,-1]
X = scale(X)

set.seed(10)

X = X[,sample(1:ncol(X), size = 20, replace = F)]

active = 3; repeats = 1000;  snr = 5; p = 20
beta = c(rep(1,active), rep(0, ncol(X) - active))
true = c(rep(1, active),rep(0,ncol(X) - active ))
r3.5.sub = pbreplicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr))


active = 3; repeats = 1000;  snr = 10; p = 20
beta = c(rep(1,active), rep(0, ncol(X) - active))
true = c(rep(1, active),rep(0,ncol(X) - active ))
r3.10.sub = pbreplicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr))



active = 7; repeats = 1000;  snr = 5; p = 20
beta = c(rep(1,active), rep(0, ncol(X) - active))
true = c(rep(1, active),rep(0,ncol(X) - active ))
r7.5.sub = pbreplicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr))

active = 7; repeats = 1000;  snr = 10; p = 20
beta = c(rep(1,active), rep(0, ncol(X) - active))
true = c(rep(1, active),rep(0,ncol(X) - active ))
r7.10.sub = pbreplicate(repeats, simulationATS(X = X, beta = beta, true = true, p = p, snr = snr))

save(r7.5.sub, r7.10.sub, file = paste0(Sys.Date(), "_r7sub.RData"))
save(r3.5.sub, r3.10.sub, file = paste0(Sys.Date(), "_r3sub.RData"))

#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--
# load("Data/r7sub.Rdata")
# load("Data/r3sub.Rdata")
filtered = c("ATS", "Exclusion ATS",
             "Static 0.60","Static 0.75", "Static 0.90", "LASSO 1SE", "Knockoff", "SCAD")
repeats = 1000





r3.5.sub.df = r3.5.sub |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==5")
r3.10.sub.df = r3.10.sub |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==10")
r3sub = rbind(r3.5.sub.df,r3.10.sub.df)
r3sub$Method = as.character(r3sub$Method)
r3sub$Method[r3sub$Method == "LASSO 1SE"] = "LASSO"
r3sub$Method[r3sub$Method == "Exclusion ATS"] = "EATS"
r3sub$Method = factor(r3sub$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
r3sub = data.frame(r3sub) |> mutate(Dimension = "n==71*`,`~p==20*`,`~`|`*beta[S]*`|`==3") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()


r7.5.sub.df = r7.5.sub |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==5")
r7.10.sub.df = r7.10.sub |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==10")
r7sub = rbind(r7.5.sub.df,r7.10.sub.df)
r7sub$Method = as.character(r3sub$Method)
r7sub$Method[r7sub$Method == "LASSO 1SE"] = "LASSO"
r7sub$Method[r7sub$Method == "Exclusion ATS"] = "EATS"
r7sub$Method = factor(r7sub$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
r7sub = data.frame(r7sub) |> mutate(Dimension = "n==71*`,`~p==20*`,`~`|`*beta[S]*`|`==7") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()



r3subplot = ggplot(r3sub,aes(x = Method, y = MCC, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_grid(Dimension ~SNR, labeller = label_parsed) +
  theme_few_grid(base_size = 20) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  ylim(-0.32, 1) + 
  xlab(element_blank()) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x= element_blank()) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))



r7subplot = ggplot(r7sub,aes(x = Method, y = MCC, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_grid(Dimension ~SNR, labeller = label_parsed) +
  theme_few_grid(base_size = 20) +
  ylim(-0.32, 1) + 
  xlab(element_blank()) + 
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x= element_blank(),
        strip.text.x = element_blank())+
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))


r3.5.sub.df.N = r3.5.sub |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==5")
r3.10.sub.df.N = r3.10.sub |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==10")
r3subN = rbind(r3.5.sub.df.N,r3.10.sub.df.N)
r3subN$Method = as.character(r3subN$Method)
r3subN$Method[r3subN$Method == "LASSO 1SE"] = "LASSO"
r3subN$Method[r3subN$Method == "Exclusion ATS"] = "EATS"
r3subN$Method[r3subN$Method == "All"] = "ATS"
r3subN$Method = factor(r3subN$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
r3subN = data.frame(r3subN) |> mutate(Dimension = "n==71*`,`~p==20*`,`~`|`*beta[S]*`|`==3") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()



r7.5.sub.df.N = r7.5.sub |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==5")
r7.10.sub.df.N = r7.10.sub |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==10")
r7subN = rbind(r7.5.sub.df.N,r7.10.sub.df.N)
r7subN$Method = as.character(r3subN$Method)
r7subN$Method[r7subN$Method == "LASSO 1SE"] = "LASSO"
r7subN$Method[r7subN$Method == "Exclusion ATS"] = "EATS"
r7subN$Method[r7subN$Method == "All"] = "ATS"
r7subN$Method = factor(r7subN$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
r7subN = data.frame(r7subN) |> mutate(Dimension = "n==71*`,`~p==20*`,`~`|`*beta[S]*`|`==7") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()


r3nplot = ggplot(r3subN,aes(x = Method, y = NN, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_grid(Dimension ~SNR, labeller = label_parsed) +
  ylab("Variables Selected") + 
  theme_few_grid(base_size = 20) +
  geom_hline( aes(yintercept = 3), linetype = "dashed", linewidth = 0.5) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  coord_cartesian(ylim = c(0, 15)) + 
  xlab(element_blank()) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x= element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))



r7nplot = ggplot(r7subN,aes(x = Method, y = NN, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_grid(Dimension ~SNR, labeller = label_parsed) +
  ylab("Variables Selected") + 
  theme_few_grid(base_size = 20) +
  coord_cartesian(ylim = c(0, 15)) + 
  xlab(element_blank()) + 
  geom_hline( aes(yintercept = 7), linetype = "dashed", linewidth = 0.5) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text.x = element_blank()) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

library(patchwork)

r3subplot/r7subplot/r3nplot/r7nplot


# PRECISION AND RECALL

r3.5.sub.df.Precision = r3.5.sub |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==5")
r3.10.sub.df.Precision = r3.10.sub |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==10")
r3Precision = rbind(r3.5.sub.df.Precision,r3.10.sub.df.Precision)
r3Precision$Method = as.character(r3Precision$Method)
r3Precision$Method[r3Precision$Method == "LASSO 1SE"] = "LASSO"
r3Precision$Method[r3Precision$Method == "Shuffle 95%"] = "Exclusion ATS"
r3Precision$Method[r3Precision$Method == "All"] = "ATS"
r3Precision$Method = factor(r3Precision$Method, levels = c("ATS", "Exclusion ATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
r3Precision = data.frame(r3Precision) |> mutate(Dimension = "n==71*`,`~p==20*`,`~`|`*beta[S]*`|`==3") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()


r7.5.sub.df.Precision = r7.5.sub |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==5")
r7.10.sub.df.Precision = r7.10.sub |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==10")
r7Precision = rbind(r7.5.sub.df.Precision,r7.10.sub.df.Precision)
r7Precision$Method = as.character(r3Precision$Method)
r7Precision$Method[r7Precision$Method == "LASSO 1SE"] = "LASSO"
r7Precision$Method[r7Precision$Method == "Shuffle 95%"] = "Exclusion ATS"
r7Precision$Method[r7Precision$Method == "All"] = "ATS"
r7Precision$Method = factor(r7Precision$Method, levels = c("ATS", "Exclusion ATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
r7Precision = data.frame(r7Precision) |> mutate(Dimension = "n==71*`,`~p==20*`,`~`|`*beta[S]*`|`==7") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()


r3nplotprecision = ggplot(r3Precision,aes(x = Method, y = Precision, fill = Method)) + geom_boxplot(alpha = 0.65) +
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


r7nplotprecision = ggplot(r7Precision,aes(x = Method, y = Precision, fill = Method)) + geom_boxplot(alpha = 0.65) +
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



r3.5.sub.df.Recall = r3.5.sub |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==5")
r3.10.sub.df.Recall = r3.10.sub |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==10")
r3Recall = rbind(r3.5.sub.df.Recall,r3.10.sub.df.Recall)
r3Recall$Method = as.character(r3Recall$Method)
r3Recall$Method[r3Recall$Method == "LASSO 1SE"] = "LASSO"
r3Recall$Method[r3Recall$Method == "Shuffle 95%"] = "Exclusion ATS"
r3Recall$Method[r3Recall$Method == "All"] = "ATS"
r3Recall$Method = factor(r3Recall$Method, levels = c("ATS", "Exclusion ATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
r3Recall = data.frame(r3Recall) |> mutate(Dimension = "n==71*`,`~p==20*`,`~`|`*beta[S]*`|`==3") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()



r7.5.sub.df.Recall = r7.5.sub |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==5")
r7.10.sub.df.Recall = r7.10.sub |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==10")
r7Recall = rbind(r7.5.sub.df.Recall,r7.10.sub.df.Recall)
r7Recall$Method = as.character(r3Recall$Method)
r7Recall$Method[r7Recall$Method == "LASSO 1SE"] = "LASSO"
r7Recall$Method[r7Recall$Method == "Shuffle 95%"] = "Exclusion ATS"
r7Recall$Method[r7Recall$Method == "All"] = "ATS"
r7Recall$Method = factor(r7Recall$Method, levels = c("ATS", "Exclusion ATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
r7Recall = data.frame(r7Recall) |> mutate(Dimension = "n==71*`,`~p==20*`,`~`|`*beta[S]*`|`==7") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()


pr = rbind(r7Recall, r3Recall)
pp = rbind(r7Precision, r3Precision)
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

##############################################################################################################