# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code to reproduce Section 5 Figures 4 & 5                       #
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
load("RData/r3sub.Rdata")
load("RData/r7sub.Rdata")
#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# 3 Active Variables MCC
filtered = c("ATS", "Exclusion ATS",
             "Static 0.60","Static 0.75", "Static 0.90", "LASSO 1SE", "Knockoff", "SCAD")
r3.5.sub.df = r3.5.sub |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==5")
r3.10.sub.df = r3.10.sub |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==10")
r3sub = rbind(r3.5.sub.df,r3.10.sub.df)
r3sub$Method = as.character(r3sub$Method)
r3sub$Method[r3sub$Method == "LASSO 1SE"] = "LASSO"
r3sub$Method = factor(r3sub$Method, levels = c("ATS", "Exclusion ATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
r3sub = data.frame(r3sub) |> mutate(Dimension = "n==71*`,`~p==20*`,`~`|`*beta[S]*`|`==3") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

# 7 Active Variables MCC
r7.5.sub.df = r7.5.sub |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==5")
r7.10.sub.df = r7.10.sub |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==10")
r7sub = rbind(r7.5.sub.df,r7.10.sub.df)
r7sub$Method = as.character(r3sub$Method)
r7sub$Method[r7sub$Method == "LASSO 1SE"] = "LASSO"
r7sub$Method = factor(r7sub$Method, levels = c("ATS", "Exclusion ATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
r7sub = data.frame(r7sub) |> mutate(Dimension = "n==71*`,`~p==20*`,`~`|`*beta[S]*`|`==7") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

# 3 Active Variables Number of Variables Selected
r3.5.sub.df.N = r3.5.sub |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==5")
r3.10.sub.df.N = r3.10.sub |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==10")
r3subN = rbind(r3.5.sub.df.N,r3.10.sub.df.N)
r3subN$Method = as.character(r3subN$Method)
r3subN$Method[r3subN$Method == "LASSO 1SE"] = "LASSO"
r3subN$Method = factor(r3subN$Method, levels = c("ATS", "Exclusion ATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
r3subN = data.frame(r3subN) |> mutate(Dimension = "n==71*`,`~p==20*`,`~`|`*beta[S]*`|`==3") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

# 7 Active Variables Number of Variables Selected
r7.5.sub.df.N = r7.5.sub |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==5")
r7.10.sub.df.N = r7.10.sub |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = "~SNR==10")
r7subN = rbind(r7.5.sub.df.N,r7.10.sub.df.N)
r7subN$Method = as.character(r3subN$Method)
r7subN$Method[r7subN$Method == "LASSO 1SE"] = "LASSO"
r7subN$Method = factor(r7subN$Method, levels = c("ATS", "Exclusion ATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
r7subN = data.frame(r7subN) |> mutate(Dimension = "n==71*`,`~p==20*`,`~`|`*beta[S]*`|`==7") |> 
  mutate(Dimension = factor(Dimension)) |> makeCluster()

r3subplot = ggplot(r3sub,aes(x = Method, y = MCC, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_wrap(Dimension ~SNR, labeller = label_parsed, ncol = 4) +
  theme_few_grid(base_size = 20) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  ylim(-0.32, 1) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group = cluster, alpha = 2)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x= element_blank()) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

r7subplot = ggplot(r7sub,aes(x = Method, y = MCC, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_wrap(Dimension ~SNR, labeller = label_parsed, ncol = 4) +
  theme_few_grid(base_size = 20) +
  ylim(-0.32, 1) + 
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x= element_blank())+
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))


r3nplot = ggplot(r3subN,aes(x = Method, y = NN, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_wrap(Dimension ~SNR, labeller = label_parsed, ncol = 4) +
  ylab("Variables Selected") + 
  theme_few_grid(base_size = 20) +
  geom_hline( aes(yintercept = 3), linetype = "dashed", linewidth = 0.5) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  coord_cartesian(ylim = c(0, 15)) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))



r7nplot = ggplot(r7subN,aes(x = Method, y = NN, fill = Method)) + geom_boxplot(alpha = 0.65) +
  facet_wrap(Dimension ~SNR, labeller = label_parsed, ncol = 4) +
  ylab("Variables Selected") + 
  theme_few_grid(base_size = 20) +
  coord_cartesian(ylim = c(0, 15)) + 
  geom_hline( aes(yintercept = 7), linetype = "dashed", linewidth = 0.5) +
  stat_summary(fun ="mean", shape = 5, size = 0.5) + 
  stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))

# Plot All
patch = (r3subplot | r7subplot)/(r3nplot | r7nplot)
suppressWarnings(print(patch + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 40))))
