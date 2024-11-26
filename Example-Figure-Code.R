# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code to reproduce Section 5 Figures                             #
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

#---#---#---#---#---#---#---#---#---#---#---#---#---#---
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



ggplot(r3sub,aes(x = Method, y = MCC, fill = Method)) + geom_boxplot(alpha = 0.65) +
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
