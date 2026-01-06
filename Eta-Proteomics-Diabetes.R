# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code to reproduce Proteomics and Diabetes                       #
# Exclusion Probability Threshold figure                            #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2025)                      #                                    
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# load("Data/pro4.RData")
# load("Data/pro9.RData")
# load("Data/d10.Rdata")
# load("Data/d5.Rdata")

pro9.1e = exclusion(pro9.1)
pro9.3e =  exclusion(pro9.3)
pro4.1e =  exclusion(pro4.1)
pro4.3e =  exclusion(pro4.3)


exclpro = data.frame("Eta" =  c(pro4.1e, pro4.3e, pro9.1e, pro9.3e),
                  "Setting" = c(rep("n==91*`,`~p==721*`,`~`|`*beta[S]*`|`==4",2000),
                                rep("n==91*`,`~p==721*`,`~`|`*beta[S]*`|`==9",2000)),
                  "SNR" = c(rep(1, 1000), rep(3,1000), rep(1, 1000), rep(3,1000))) |>
  mutate(SNR = factor(SNR, levels = c(1,3)),
         Setting = factor(Setting, levels = c("n==91*`,`~p==721*`,`~`|`*beta[S]*`|`==9",
                                              "n==91*`,`~p==721*`,`~`|`*beta[S]*`|`==4")))

exclproplot = ggplot(exclpro, aes(x = SNR, y = Eta)) + geom_boxplot(alpha = 0.65) +
  facet_wrap(~Setting, labeller = label_parsed, ncol = 4) +
  theme_few_grid(base_size = 20) +
  ylab(TeX("Proteomics $\\eta$")) + 
  stat_summary(fun ="mean", shape = 5, size = 0.1) + 
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x= element_blank()) 

d5.1e = exclusion(d5.1)
d5.3e =  exclusion(d5.3)
d10.1e =  exclusion(d10.1)
d10.3e =  exclusion(d10.3)


excld = data.frame("Eta" =  c(d5.1e, d5.3e, d10.1e, d10.3e),
                     "Setting" = c(rep("n==442*`,`~p==64*`,`~`|`*beta[S]*`|`==5",2000),
                                   rep("n==442*`,`~p==64*`,`~`|`*beta[S]*`|`==10",2000)),
                     "SNR" = c(rep(1, 1000), rep(3,1000), rep(1, 1000), rep(3,1000))) |> mutate(SNR = factor(SNR, levels = c(1,3)))

excldplot = ggplot(excld, aes(x = SNR, y = Eta)) + geom_boxplot(alpha = 0.65) +
  facet_wrap(~Setting, labeller = label_parsed, ncol = 4) +
  theme_few_grid(base_size = 20) +
  ylab(TeX("Diabetes $\\eta$")) + 
  stat_summary(fun ="mean", shape = 5, size = 0.1) + 
  theme(legend.position = "none") 

exclprod = exclproplot/excldplot
exclprod
