# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4 Results                     #
# in the paper: Data-Adaptive Automatic Threshold Calibration       #
#  for Stability Selection (Huang et al. 2025)                      #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
library(tidyverse)
library(ncvreg)
library(Rcpp)
library(lars)
library(stabs)
library(mvtnorm)
library(knockoff)
library(latex2exp)
library(hdi)
library(patchwork)
library(hrbrthemes)
library(pbapply)

filtered = c("ATS", "Exclusion ATS",
             "Static 0.60","Static 0.75", "Static 0.90", "LASSO 1SE", "Knockoff", "SCAD")
repeats = 1000

simulationATS = function(X, true, p, beta,  snr = 10, gaussian.knockoffs = F){
  
  # Ensure SNR
  signal = sqrt(mean((as.matrix(X) %*% as.matrix(beta))^2))
  sigma = as.numeric(signal/sqrt(snr))
  
  # Compute Y with SNR
  Y = as.matrix(X)%*%as.matrix(beta) + rnorm(nrow(X), 0, sd = sigma) 
  
  # Stability Selection
  s = stabs::stabsel(x = X, y = Y, B = 100,
                     fitfun = stabs::lars.lasso, PFER = 5, cutoff = 0.75,
                     sampling.type = "MB")
  
  # Shuffle data
  idx = sample(1:nrow(X), replace = F)
  rX = X[idx,]
  idxPushed = c(tail(idx, 1), head(idx, -1))
  rY = Y[idxPushed] |> as.matrix(ncol = 1)
  
  # Exclusion Probability Threshold
  sMix = stabs::stabsel(x = rX, y = rY, B = 100,
                        fitfun = stabs::lars.lasso, PFER = 5, cutoff = 0.75,
                        sampling.type = "MB")
  
  sMix_prob = sort(sMix$max, decreasing = T)
  mix_exclusion = quantile(sMix_prob, 0.95)
  
  
  # Stable set for static values
  
  # 0.6 
  idx_chosen = sort(s$max, decreasing = T) >= 0.6
  q_static_selected = names(sort(s$max, decreasing = T))[idx_chosen]
  pred = makePred(q_static_selected, X,p)
  q_60_acc = acc(pred, true, p)
  
  # 0.75
  idx_chosen = sort(s$max, decreasing = T) >= 0.75
  q_static_selected = names(sort(s$max, decreasing = T))[idx_chosen]
  pred = makePred(q_static_selected, X,p)
  q_static_acc = acc(pred, true, p)
  
  
  # 0.9
  idx_chosen = sort(s$max, decreasing = T) >= 0.9
  q_static_selected = names(sort(s$max, decreasing = T))[idx_chosen]
  pred = makePred(q_static_selected, X,p)
  q_90_acc = acc(pred, true, p)
  
  
  
  
  # ATS 
  ATS = convert(s)
  if (length(ATS) == 2){ATS = c(ATS, ATS[1])}
  ATS = ATS |> getR()
  ATS_pi =  sort(s$max, decreasing = T)[ATS]
  ATS_selected = sort(s$max, decreasing = T)[1:ATS] |> names()
  pred = makePred(ATS_selected, X,p)
  ATS_acc = acc(pred, true, p)
  
  
  # Exclusion ATS  
  EATS =  convert(s)[1:(length(convert(s)))][convert(s)[1:(length(convert(s)))] >= 100*mix_exclusion] 
  if (length(EATS) == 2){EATS = c(EATS, EATS[1])}
  EATS = EATS |> getR()
  EATS_pi = sort(s$max, decreasing = T)[EATS]
  EATS_selected = sort(s$max, decreasing = T)[1:EATS] |> names()
  pred = makePred(EATS_selected, X,p)
  EATS_acc = acc(pred, true, p)
  
  # Knockoff
  if (gaussian.knockoffs == T){
    gaussian_knockoffs = function(X) knockoff::create.second_order(X, method='equi', shrink=T)
    kf = knockoff::knockoff.filter(X, Y, knockoffs = gaussian_knockoffs)
  }else{
    kf = knockoff::knockoff.filter(X, Y)
  }
  ko_selected = names(kf$selected)
  pred = makePred(ko_selected, X, p)
  ko_acc = acc(pred, true, p)
  
  
  # LASSO 
  l = cv.glmnet(X, Y)
  c1se = which(coef(l, lambda = l$lambda.1se)[-1] != 0)
  c1se_selected = colnames(X)[c1se]
  pred = makePred(c1se_selected, X, p)
  l1se_acc = acc(pred, true, p)
  
  cmin = which(coef(l, lambda = l$lambda.min)[-1] != 0)
  cmin_selected = colnames(X)[cmin]
  pred = makePred(cmin_selected, X, p)
  lmin_acc = acc(pred, true, p)
  
  
  # SCAD
  nr = cv.ncvreg(X,Y, penalty = "SCAD")
  nrcoef = which(coef(nr, lambda = nr$lambda.min)[-1] != 0)
  nr_selected = colnames(X)[nrcoef]
  pred = makePred(nr_selected, X, p)
  scad_acc = acc(pred, true, p)
  
  
  # Output
  list("0.75" = q_static_acc,
       "0.60" = q_60_acc,
       "0.90" = q_90_acc,
       "ATS" = ATS_acc,
       "Exclusion ATS" = EATS_acc,
       "LASSO 1SE" = l1se_acc,
       "LASSO MIN" = lmin_acc,
       "Knockoff" = ko_acc,
       "SCAD" = scad_acc,
       "Adaptive Exclusion Probability" = mix_exclusion,
       "EATS Pi" = EATS_pi,
       "ATS Pi" = ATS_pi)
}

exclusion = function(x){
  x["Adaptive Exclusion Probability",] |> unlist() |> as.vector()
}

extractPi = function(x, snr, setting = NULL, data = F, active = NULL){
  ep = x["EATS Pi",] |> unlist() |> as.vector()  
  ap = x["ATS Pi",] |> unlist() |> as.vector()  
  if (data == F){
    dimension = c("n = 20, p = 1000, active = 2",
                  "n = 100, p = 500, active = 10",
                  "n = 200, p = 200, active = 20",
                  "n = 500, p = 100, active = 20")
    labels  = c("(I):~n==20*`,`~p==1000*`,`~`|`*beta[S]*`|`==2",
                "(II):~n==100*`,`~p==500*`,`~`|`*beta[S]*`|`==10",
                "(III):~n==200*`,`~p==200*`,`~`|`*beta[S]*`|`==20",
                "(IV):~n==500*`,`~p==100*`,`~`|`*beta[S]*`|`==20")
    SNR = c(0.5,1,2, 3)
    SNRlabels = c("~SNR==0.5", "~SNR==1", "~SNR==2", "~SNR==3")
    
    d = data.frame("EATS" = ep, "ATS" = ap) |> reshape2::melt() |> 
      mutate(variable = factor(variable, levels = c("ATS", "EATS")))
    d$dimension = dimension[setting] 
    d = d |> mutate(dimension = factor(dimension, labels = labels[setting]),
                    SNR = SNR[SNR == snr]) 
    d$SNR = factor(d$SNR, levels = snr, labels = SNRlabels[SNR == snr])
  }else if (data == "proteomics"){
    labels = c("n==91*`,`~p==721*`,`~`|`*beta[S]*`|`==9",
                    "n==91*`,`~p==721*`,`~`|`*beta[S]*`|`==4")
    dimension = c("n = 91, p = 721, active = 9",
                  "n = 91, p = 721, active = 4")
    SNR = c(1,3)
    SNRlabels = c("~SNR==1", "~SNR==3")
    d = data.frame("EATS" = ep, "ATS" = ap) |> reshape2::melt() |> 
      mutate(variable = factor(variable, levels = c("ATS", "EATS")))
    d$dimension = ifelse(active == 9, dimension[1], dimension[2])
    d = d |> mutate(dimension = factor(dimension, labels = ifelse(active == 9, labels[1], labels[2])),
                    SNR = snr) 
    d$SNR = factor(d$SNR, levels = snr, labels = ifelse(snr == 1, SNRlabels[1], SNRlabels[2]))
  }else if(data == "diabetes"){
    labels = c("n==442*`,`~p==64*`,`~`|`*beta[S]*`|`==10",
               "n==442*`,`~p==64*`,`~`|`*beta[S]*`|`==5")
    dimension = c("n = 442, p = 64, active = 10",
                  "n = 442, p = 64, active = 5")
    SNR = c(1,3)
    SNRlabels = c("~SNR==1", "~SNR==3")
    d = data.frame("EATS" = ep, "ATS" = ap) |> reshape2::melt() |> 
      mutate(variable = factor(variable, levels = c("ATS", "EATS")))
    d$dimension = ifelse(active == 10, dimension[1], dimension[2])
    d = d |> mutate(dimension = factor(dimension, labels = ifelse(active == 10, labels[1], labels[2])),
                    SNR = snr) 
    d$SNR = factor(d$SNR, levels = snr, labels = ifelse(snr == 1, SNRlabels[1], SNRlabels[2]))
  }
  d
}



convert = function(s){
  return((as.vector(s$max) |> sort(decreasing = T))*100)
}

#### Data Cleaning Functions ###
makePred = function(q_selected, X, p){
  if (length(q_selected) == 0){
    pred = rep(0, p)
  }else {
    pred = rep(0,p)
    pred[which(colnames(X) %in% q_selected)] = 1
  }
  pred
}

acc = function(pred, true, p){
  tp = sum(pred == 1 & true == 1)  
  fp = sum(pred == 1 & true == 0)  
  tn = sum(pred == 0 & true == 0)  
  fn = sum(pred == 0 & true == 1)  
  list("TP" = tp, "FP" = fp, "TN" = tn, "FN" = fn)
}

extractMCC = function(ql){
  tp = sapply(ql, "[[", 1) 
  fp = sapply(ql, "[[", 2)
  tn = sapply(ql, "[[", 3)
  fn = sapply(ql, "[[", 4)
  MCC = numeric(length(tp))
  for (i in 1:length(tp)){
    MCC[i] = mltools::mcc(TP = tp[i],
                          FP = fp[i],
                          FN = fn[i],
                          TN = tn[i])
  }
  MCC
}

prep = function(l){
  tp = l[["TP"]]
  fp = l[["FP"]]
  tn = l[["TN"]]
  fn = l[["FN"]]
  c(tp,fp,tn,fn)
}

cleanMCC = function(x, repeats = 1000){
  z0.75 = lapply(x[1,], prep)
  z0.60 =  lapply(x[2,], prep)
  z0.90 =  lapply(x[3,], prep)
  ATS = lapply(x[4,],prep)
  EATS = lapply(x[5,], prep)
  lasso_1se = lapply(x[6,], prep)
  lasso_min = lapply(x[7,], prep)
  ko = lapply(x[8,],prep)
  scad = lapply(x[9,],prep)
  
  MCCdf = data.frame("MCC" = c(extractMCC(z0.75),
                               extractMCC(z0.60),
                               extractMCC(z0.90),
                               extractMCC(ATS),
                               extractMCC(EATS),
                               extractMCC(lasso_1se),
                               extractMCC(lasso_min),
                               extractMCC(ko),
                               extractMCC(scad)),
                     "Method" = c(rep("Static 0.75", repeats),
                                  rep("Static 0.60", repeats),
                                  rep("Static 0.90", repeats),
                                  rep("ATS", repeats),
                                  rep("Exclusion ATS", repeats),
                                  rep("LASSO 1SE", repeats),
                                  rep("LASSO MIN", repeats),
                                  rep("Knockoff", repeats),
                                  rep("SCAD", repeats)))
  
  MCCdf$Method = factor(MCCdf$Method, levels = c("ATS", "Exclusion ATS",
                                                 "Static 0.60", "Static 0.75", "Static 0.90",
                                                 "LASSO 1SE", "LASSO MIN", "Knockoff", "SCAD"))
  MCCdf
}

cleanN = function(x, repeats = 1000){
  z0.75 = lapply(x[1,], prep)
  z0.60 =  lapply(x[2,], prep)
  z0.90 =  lapply(x[3,], prep)
  ATS = lapply(x[4,],prep)
  EATS = lapply(x[5,], prep)
  lasso_1se = lapply(x[6,], prep)
  lasso_min = lapply(x[7,], prep)
  ko = lapply(x[8,],prep)
  scad = lapply(x[9,],prep)
  
  Ndf = data.frame("NN" = c(extractN(z0.75),
                            extractN(z0.60),
                            extractN(z0.90),
                            extractN(ATS),
                            extractN(EATS),
                            extractN(lasso_1se),
                            extractN(lasso_min),
                            extractN(ko),
                            extractN(scad)),
                   "Method" = c(rep("Static 0.75", repeats),
                                rep("Static 0.60", repeats),
                                rep("Static 0.90", repeats),
                                rep("ATS", repeats),
                                rep("Exclusion ATS", repeats),
                                rep("LASSO 1SE", repeats),
                                rep("LASSO MIN", repeats),
                                rep("Knockoff", repeats),
                                rep("SCAD", repeats)))
  
  Ndf$Method = factor(Ndf$Method, levels = c("ATS", "Exclusion ATS",
                                             "Static 0.60", "Static 0.75", "Static 0.90",
                                             "LASSO 1SE", "LASSO MIN", "Knockoff", "SCAD"))
  Ndf
}

cleanPrecision = function(x){
  z0.75 = lapply(x[1,], prep)
  z0.60 =  lapply(x[2,], prep)
  z0.90=  lapply(x[3,], prep)
  ats = lapply(x[4,],prep)
  eats = lapply(x[5,], prep)
  lasso_1se = lapply(x[6,], prep)
  lasso_min = lapply(x[7,], prep)
  ko = lapply(x[8,],prep)
  scad = lapply(x[9,],prep)
  
  P = data.frame("Precision" = c(extractPrecision(z0.75),
                                   extractPrecision(z0.60),
                                   extractPrecision(z0.90),
                                   extractPrecision(ats),
                                   extractPrecision(eats),
                                   extractPrecision(lasso_1se),
                                   extractPrecision(lasso_min),
                                   extractPrecision(ko),
                                   extractPrecision(scad)),
                   "Method" = c(rep("Static 0.75", repeats),
                                rep("Static 0.60", repeats),
                                rep("Static 0.90", repeats),
                                rep("ATS", repeats),
                                rep("Exclusion ATS", repeats),
                                rep("LASSO 1SE", repeats),
                                rep("LASSO MIN", repeats),
                                rep("Knockoff", repeats),
                                rep("SCAD", repeats)))
  
  P$Method = factor(P$Method, levels = c("ATS", "Exclusion ATS",
                                             "Static 0.60", "Static 0.75", "Static 0.90",
                                             "LASSO 1SE", "LASSO MIN", "Knockoff", "SCAD"))
  P
}

cleanRecall = function(x){
  z0.75 = lapply(x[1,], prep)
  z0.60 =  lapply(x[2,], prep)
  z0.90=  lapply(x[3,], prep)
  ats = lapply(x[4,],prep)
  eats = lapply(x[5,], prep)
  lasso_1se = lapply(x[6,], prep)
  lasso_min = lapply(x[7,], prep)
  ko = lapply(x[8,],prep)
  scad = lapply(x[9,],prep)
  
  R = data.frame("Recall" = c(extractRecall(z0.75),
                                extractRecall(z0.60),
                                extractRecall(z0.90),
                                extractRecall(ats),
                                extractRecall(eats),
                                extractRecall(lasso_1se),
                                extractRecall(lasso_min),
                                extractRecall(ko),
                                extractRecall(scad)),
                   "Method" = c(rep("Static 0.75", repeats),
                                rep("Static 0.60", repeats),
                                rep("Static 0.90", repeats),
                                rep("ATS", repeats),
                                rep("Exclusion ATS", repeats),
                                rep("LASSO 1SE", repeats),
                                rep("LASSO MIN", repeats),
                                rep("Knockoff", repeats),
                                rep("SCAD", repeats)))
  
  R$Method = factor(R$Method, levels = c("ATS", "Exclusion ATS",
                                             "Static 0.60", "Static 0.75", "Static 0.90",
                                             "LASSO 1SE", "LASSO MIN", "Knockoff", "SCAD"))
  R
}

extractPrecision = function(ql){
  fp = sapply(ql, "[[", 2)
  tp = sapply(ql, "[[", 1) 
  fn = sapply(ql, "[[", 4)
  tn = sapply(ql, "[[", 3)
  ans = c()
  for (i in 1:length(tn)){
    if ((tp[i] == 0 & fp[i] > 0) | (tp[i] == 0 & fn[i] > 0)){
      ans[i] = 0
    }else{
      ans[i] = tp[i]/(tp[i] + fp[i])
    }
  }
  ans
}

extractRecall = function(ql){
  fp = sapply(ql, "[[", 2)
  tp = sapply(ql, "[[", 1) 
  fn = sapply(ql, "[[", 4)
  tn = sapply(ql, "[[", 3)
  tp/(tp + fn)
}

extractN = function(ql){
  tp = sapply(ql, "[[", 1) 
  fp = sapply(ql, "[[", 2)
  tp + fp
}

### ATS Function ###
Rcpp::cppFunction('
int getR(const NumericVector& d) {
  int p = d.size();
  NumericVector lq(p, 0.0);
  NumericVector sigma2(p);
  for (int q = 0; q < p; q++) {
    NumericVector d1 = head(d, q + 1);
    NumericVector d2 = tail(d, p - (q + 1));
    double mu1 = mean(d1);
    double mu2 = mean(d2);
    sigma2[q] = (sum(pow(d1 - mu1, 2)) + sum(pow(d2 - mu2, 2))) / (p - 2);
    lq[q] = sum(dnorm(d1, mu1, sqrt(sigma2[q]), true)) +
      sum(dnorm(d2, mu2, sqrt(sigma2[q]), true));
  }
  return which_max(lq) + 1;
}
')

### Generating Data ###
extend <- function(alphabet) function(i) {
  base10toA <- function(n, A) {
    stopifnot(n >= 0L)
    N <- length(A)
    j <- n %/% N 
    if (j == 0L) A[n + 1L] else paste0(Recall(j - 1L, A), A[n %% N + 1L])
  }   
  vapply(i-1L, base10toA, character(1L), alphabet)
}

gendata = function(n, p, active){
  # Add letters as variable names
  MORELETTERS <- extend(LETTERS)     
  moreletters = MORELETTERS(c(1:p))
  
  # Generate random data from multivariate gaussian with Toeplitz covariance matrix
  toeshd = 0.5^abs(row(matrix(1:p, p, p)) - col(matrix(1:p, p, p)))
  x = mvtnorm::rmvnorm(n  = n, sigma = toeshd)
  
  # Randomly generate beta coefficient
  beta = rep(0,p)
  signs = sample(c(-1,1), active, replace = T)
  coefs = sample(1:3, active, replace = T)
  coefsigns = signs * coefs
  beta[1:active] = sample(coefsigns, active, replace = T)
  colnames(x) = moreletters
  names(beta) = moreletters[1:length(beta)]
  y <- x %*% beta
  
  list("Y" = y, "X" = x, "beta" = beta)
}

### Plotting ### 
makeCluster = function(c){
  c$cluster = case_when(c$Method == "All" ~ "All",
                        c$Method %in% c("Quarter", "Quarter Excl 0.2", "Quarter Excl 0.1",
                                        "Quarter Adaptive Excl 95%", "Quarter Adaptive Excl 99%",
                                        "Quarter Shuffle 95%") ~ "Quarter",
                        c$Method %in% c("Static 0.60", "Static 0.75", "Static 0.90") ~ "Static",
                        c$Method %in% c("LASSO 1SE", "LASSO MIN") ~ "LASSO",
                        .default = as.character(c$Method))
  c
}

theme_few_grid = function (base_size = 12, base_family = "") 
{
  gray <- "#4D4D4D"
  black <- "#000000"
  theme_bw(base_size = base_size, base_family = base_family) + 
    theme(line = element_line(colour = gray), rect = element_rect(fill = "white", 
                                                                  colour = NA), text = element_text(colour = black), 
          axis.ticks = element_line(colour = gray), legend.key = element_rect(colour = NA), 
          panel.border = element_rect(colour = gray), panel.grid = element_line(color = alpha("black", 0.05)), 
          strip.background = element_rect(fill = "#f7f9fc", colour = gray))
}

totplotnoaxis = function(TOT){
  ggplot(TOT,aes(x = Method, y = MCC, fill = Method)) + geom_boxplot(alpha = 0.65) +
    facet_grid(Dimension ~SNR, labeller = label_parsed) +
    theme_few_grid(base_size = 20) +
    stat_summary(fun ="mean", shape = 5, size = 0.5) + 
    stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x= element_blank()) +
    xlab(element_blank()) +
    scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))
}


totplot = function(TOT){
  ggplot(TOT,aes(x = Method, y = MCC, fill = Method)) + geom_boxplot(alpha = 0.65) +
    facet_grid(Dimension ~SNR, labeller = label_parsed) +
    theme_few_grid(base_size = 20) +
    stat_summary(fun ="mean", shape = 5, size = 0.5) + 
    stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.text.x = element_blank()) +
    xlab(element_blank()) + 
    scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))
}

combine = function(a,b,c,d = NULL, ref, filtered = NULL){
  if (is.null(filtered) == T){
    filtered = c("All", "Quarter Adaptive 95%",
                 "Quarter Adaptive 99%", "Quarter Shuffle 95%", "Shuffle 95%",
                 "Shuffle Weighted", "Static 0.60","Static 0.75", "Static 0.90", "LASSO 1SE", "LASSO MIN")
  }
  dimension = c("n = 20, p = 1000, active = 2",
                "n = 100, p = 500, active = 10",
                "n = 200, p = 200, active = 20",
                "n = 500, p = 100, active = 20")
  labels  = c("(I):~n==20*`,`~p==1000*`,`~`|`*beta[S]*`|`==2",
              "(II):~n==100*`,`~p==500*`,`~`|`*beta[S]*`|`==10",
              "(III):~n==200*`,`~p==200*`,`~`|`*beta[S]*`|`==20",
              "(IV):~n==500*`,`~p==100*`,`~`|`*beta[S]*`|`==20")
  SNR = c(0.5,1,2, 3)
  SNRlabels = c("~SNR==0.5", "~SNR==1", "~SNR==2", "~SNR==3")
  r = rbind(a |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[1], N = "A"),
            b |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[2], N = "B"),
            c |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[3], N = "C"),
            d |> cleanMCC() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[4], N = "C"))
  r$Method = as.character(r$Method)
  r$Method[r$Method == "LASSO 1SE"] = "LASSO"
  r$Method[r$Method == "Exclusion ATS"] = "EATS"
  r$Method[r$Method == "All"] = "ATS"
  r$Method = factor(r$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
  r = data.frame(r) |> mutate(Dimension = dimension[ref]) |>
    mutate(Dimension = factor(Dimension, labels = labels[ref]),
           SNR = factor(SNR, levels = c("0.5", "1", "2", "3"),
                        labels = c(SNRlabels[1],
                                   SNRlabels[2],
                                   SNRlabels[3],
                                   SNRlabels[4])))
  return(r)
}

combineN = function(a,b,c = NULL,d, ref,filtered = NULL){
  if (is.null(filtered) == T){
    filtered = c("ATS", "Exclusion ATS", "Static 0.60","Static 0.75", "Static 0.90", "LASSO 1SE", "LASSO MIN")
  }
  
  dimension = c("n = 20, p = 1000, active = 2",
                "n = 100, p = 500, active = 10",
                "n = 200, p = 200, active = 20",
                "n = 500, p = 100, active = 20")
  labels  = c("(I):~n==20*`,`~p==1000*`,`~`|`*beta[S]*`|`==2",
              "(II):~n==100*`,`~p==500*`,`~`|`*beta[S]*`|`==10",
              "(III):~n==200*`,`~p==200*`,`~`|`*beta[S]*`|`==20",
              "(IV):~n==500*`,`~p==100*`,`~`|`*beta[S]*`|`==20")
  SNR = c(0.5,1,2, 3)
  SNRlabels = c("~SNR==0.5", "~SNR==1", "~SNR==2", "~SNR==3")

  r = rbind(a |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[1], N = "A"),
            b |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[2], N = "B"),
            c |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[3], N = "C"),
            d |> cleanN() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[4], N = "C"))
  r$Method = as.character(r$Method)
  r$Method[r$Method == "LASSO 1SE"] = "LASSO"
  r$Method[r$Method == "Exclusion ATS"] = "EATS"
  r$Method[r$Method == "All"] = "ATS"
  r$Method = factor(r$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
  r = data.frame(r) |> mutate(Dimension = dimension[ref]) |>
    mutate(Dimension = factor(Dimension, labels = labels[ref]),
           SNR = factor(SNR, levels = c("0.5", "1", "2", "3"),
                        labels = c(SNRlabels[1],
                                   SNRlabels[2],
                                   SNRlabels[3],
                                   SNRlabels[4])))
  return(r)
}

totplotN = function(TOT, lim = max(TOT$NN)){
  str = unique(TOT$Dimension) |> as.character() |> strsplit(split = "|")
  active = numeric(length(str))
  for (i in 1:length(str)){
    if (str[[i]] |> last() == 0){
      active[i] = paste(str[[i]] |> nth(-2), str[[i]] |> last(), sep = "") |> as.numeric()
    }else{
      active[i] = as.numeric(str[[i]] |> last())
    }
  }
  TOTdummy = TOT |> group_by(Dimension) |> summarise(active = mean(NN))
  TOTdummy$active = active
  ggplot(TOT,aes(x = Method, y = NN, fill = Method)) + geom_boxplot(alpha = 0.65) +
    facet_grid(Dimension ~SNR, labeller = label_parsed) +
    ylab("Variables Selected") + 
    geom_hline(data = TOTdummy, aes(yintercept = active), linetype = "dashed", linewidth = 0.5) + 
    coord_cartesian(ylim = c(0, lim)) + 
    theme_few_grid(base_size = 20) +
    xlab(element_blank()) + 
    stat_summary(fun ="mean", shape = 5, size = 0.5) + 
    stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.text.x = element_blank()) + 
    scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))
}

totplotNnoaxis = function(TOT, lim = max(TOT$NN)){
  str = unique(TOT$Dimension) |> as.character() |> strsplit(split = "|")
  active = numeric(length(str))
  for (i in 1:length(str)){
    if (str[[i]] |> last() == 0){
      active[i] = paste(str[[i]] |> nth(-2), str[[i]] |> last(), sep = "") |> as.numeric()
    }else{
      active[i] = as.numeric(str[[i]] |> last())
    }
  }
  TOTdummy = TOT |> group_by(Dimension) |> summarise(active = mean(NN))
  TOTdummy$active = active
  ggplot(TOT,aes(x = Method, y = NN, fill = Method)) + geom_boxplot(alpha = 0.65) +
    facet_grid(Dimension ~SNR, labeller = label_parsed) +
    ylab("Variables Selected") +
    geom_hline(data = TOTdummy, aes(yintercept = active), linetype = "dashed", linewidth = 0.5) + 
    coord_cartesian(ylim = c(0, lim)) + 
    theme_few_grid(base_size = 20) +
    stat_summary(fun ="mean", shape = 5, size = 0.5) + 
    stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x= element_blank()) + 
    scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))
}

combineRecall = function(a,b,c = NULL,d, ref,filtered = NULL, ribo = F){
  if (is.null(filtered) == T){
    filtered = c("All", "Quarter Adaptive 95%",
                 "Quarter Adaptive 99%", "Quarter Shuffle 95%", "Shuffle 95%",
                 "Shuffle Weighted", "Static 0.60","Static 0.75", "Static 0.90", "LASSO 1SE", "LASSO MIN")
  }
  
  dimension = c("n = 20, p = 1000, active = 2",
                "n = 100, p = 500, active = 10",
                "n = 200, p = 200, active = 20",
                "n = 500, p = 100, active = 20")
  labels  = c("(I):~n==20*`,`~p==1000*`,`~`|`*beta[S]*`|`==2",
              "(II):~n==100*`,`~p==500*`,`~`|`*beta[S]*`|`==10",
              "(III):~n==200*`,`~p==200*`,`~`|`*beta[S]*`|`==20",
              "(IV):~n==500*`,`~p==100*`,`~`|`*beta[S]*`|`==20")
  SNR = c(0.5,1,2,3)
  SNRlabels = c("~SNR==0.5", "~SNR==1", "~SNR==2", "~SNR==3")
  
  r = rbind(a |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[1], N = "A"),
            b |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[2], N = "B"),
            c |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[3], N = "C"),
            d |> cleanRecall() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[4], N = "C"))
  r$Method = as.character(r$Method)
  r$Method[r$Method == "LASSO 1SE"] = "LASSO"
  r$Method[r$Method == "Exclusion ATS"] = "EATS"
  r$Method[r$Method == "All"] = "ATS"
  r$Method = factor(r$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
  r = data.frame(r) |> mutate(Dimension = dimension[ref]) |>
    mutate(Dimension = factor(Dimension, labels = labels[ref]),
           SNR = factor(SNR, levels = c("0.5", "1", "2", "3"),
                        labels = c(SNRlabels[1],
                                   SNRlabels[2],
                                   SNRlabels[3],
                                   SNRlabels[4])))
  return(r)
}

totplotRecall = function(TOT){
  ggplot(TOT,aes(x = Method, y = Recall, fill = Method)) + geom_boxplot() +
    facet_grid(Dimension ~SNR, labeller = label_parsed, scales = "free_y") +
    ylab("Recall") + 
    theme_few_grid(base_size = 20) +
    stat_summary(fun ="mean", shape = 5, size = 0.5) + 
    stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.text.x = element_blank()) + 
    scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))
}

combinePrecision = function(a,b,c = NULL,d, ref,filtered = NULL, ribo = F){
  if (is.null(filtered) == T){
    filtered = c("All", "Quarter Adaptive 95%",
                 "Quarter Adaptive 99%", "Quarter Shuffle 95%", "Shuffle 95%",
                 "Shuffle Weighted", "Static 0.60","Static 0.75", "Static 0.90", "LASSO 1SE", "LASSO MIN")
  }
  
  dimension = c("n = 20, p = 1000, active = 2",
                "n = 100, p = 500, active = 10",
                "n = 200, p = 200, active = 20",
                "n = 500, p = 100, active = 20")
  labels  = c("(I):~n==20*`,`~p==1000*`,`~`|`*beta[S]*`|`==2",
              "(II):~n==100*`,`~p==500*`,`~`|`*beta[S]*`|`==10",
              "(III):~n==200*`,`~p==200*`,`~`|`*beta[S]*`|`==20",
              "(IV):~n==500*`,`~p==100*`,`~`|`*beta[S]*`|`==20")
  SNR = c(0.5,1,2, 3)
  SNRlabels = c("~SNR==0.5", "~SNR==1", "~SNR==2", "~SNR==3")
  
  r = rbind(a |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[1], N = "A"),
            b |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[2], N = "B"),
            c |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[3], N = "C"),
            d |> cleanPrecision() |> dplyr::filter(Method %in% filtered) |> mutate(SNR = SNR[4], N = "C"))
  r$Method = as.character(r$Method)
  r$Method[r$Method == "LASSO 1SE"] = "LASSO"
  r$Method[r$Method == "Exclusion ATS"] = "EATS"
  r$Method[r$Method == "All"] = "ATS"
  r$Method = factor(r$Method, levels = c("ATS", "EATS", "Static 0.60", "Static 0.75", "Static 0.90", "LASSO","Knockoff","SCAD"))
  r = data.frame(r) |> mutate(Dimension = dimension[ref]) |>
    mutate(Dimension = factor(Dimension, labels = labels[ref]),
           SNR = factor(SNR, levels = c("0.5", "1", "2", "3"),
                        labels = c(SNRlabels[1],
                                   SNRlabels[2],
                                   SNRlabels[3],
                                   SNRlabels[4])))
  return(r)
}

totplotPrecision = function(TOT){
  ggplot(TOT,aes(x = Method, y = Precision, fill = Method)) + geom_boxplot() +
    facet_wrap(Dimension ~SNR, labeller = label_parsed, ncol = 4) +
    ylab("Recall") + 
    theme_few_grid(base_size = 20) +
    stat_summary(fun ="mean", shape = 5, size = 0.5) + 
    stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))
}

prplot1 = function(recall, precision, recall1, precision1){
  r = ggplot(recall, aes(x = Method, y = Recall, fill = Method)) +
    geom_boxplot() + facet_grid(Dimension~SNR, labeller = label_parsed) +
    ylab("Recall") + theme_few_grid(base_size = 20) +
    stat_summary(fun ="mean", shape = 5, size = 0.5) + 
    stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x= element_blank()) + 
    scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))
  
  p = ggplot(precision, aes(x = Method, y = Precision, fill = Method)) +
    geom_boxplot() + facet_grid(Dimension~SNR, labeller = label_parsed) +
    ylab("Precision") + theme_few_grid(base_size = 20) +
    xlab(element_blank()) + 
    stat_summary(fun ="mean", shape = 5, size = 0.5) + 
    stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          strip.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))                    
  
  rr = ggplot(recall1, aes(x = Method, y = Recall, fill = Method)) +
    geom_boxplot() + facet_grid(Dimension~SNR, labeller = label_parsed) +
    ylab("Recall") + theme_few_grid(base_size = 20) +
    stat_summary(fun ="mean", shape = 5, size = 0.5) + 
    stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x= element_blank(),
          strip.text.x = element_blank()) + 
    scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))
  
  pp = ggplot(precision1, aes(x = Method, y = Precision, fill = Method)) +
    geom_boxplot() + facet_grid(Dimension~SNR, labeller = label_parsed) +
    ylab("Precision") + theme_few_grid(base_size = 20) +
    stat_summary(fun ="mean", shape = 5, size = 0.5) + 
    stat_summary(fun= "mean", geom="line", linetype ="solid", linewidth = 0.5,  aes(group= cluster, alpha = 2)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.text.x = element_blank()) + 
    scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"))                    
  library(patchwork)
  
  r/p/rr/pp
}


### FDR Function ###
fdr = function(X, p, beta, snr, EV, LOOPS = 1000){
  n.selected = false.selections = correct.selections.prop = method = pi = c()
  for (i in 1:LOOPS){
    signal = sqrt(mean((as.matrix(X) %*% as.matrix(beta))^2))
    sigma = as.numeric(signal/sqrt(snr))
    
    # Compute Y with SNR
    Y = as.matrix(X)%*%as.matrix(beta) + rnorm(nrow(X), 0, sd = sigma) 
    
    # Stability Selection
    s = stabs::stabsel(x = X, y = Y, B = 100,
                       fitfun = stabs::lars.lasso, PFER = 5, cutoff = 0.75,
                       sampling.type = "MB")
    
    # Shuffle data
    idx = sample(1:nrow(X), replace = F)
    rX = X[idx,]
    idxPushed = c(tail(idx, 1), head(idx, -1))
    rY = Y[idxPushed] |> as.matrix(ncol = 1)
    
    # Exclusion Probability Threshold
    sMix = stabs::stabsel(x = rX, y = rY, B = 100,
                          fitfun = stabs::lars.lasso, PFER = 5, cutoff = 0.75,
                          sampling.type = "MB")
    
    sMix_prob = sort(sMix$max, decreasing = T)
    mix_exclusion = quantile(sMix_prob, 0.95)
    
    EATS =  convert(s)[1:(length(convert(s)))][convert(s)[1:(length(convert(s)))] >= 100*mix_exclusion] 
    if (length(EATS) == 2){EATS = c(EATS, EATS[1])}
    EATSix = EATS |> getR()
    pi.hat = EATS[EATSix]/100
    if (pi.hat <= 0.5){
      s = stabs::stabsel(x = X, y = Y, B = 100,
                         fitfun = stabs::lars.lasso, PFER = EV, cutoff = 0.501,
                         sampling.type = "MB")
    }else{
      s = stabs::stabsel(x = X, y = Y, B = 100,
                         fitfun = stabs::lars.lasso, PFER = EV, cutoff = pi.hat,
                         sampling.type = "MB")
    }
    
    n.selected[i] = length(s$selected)
    false.selections[i] = sum(!(s$selected %in% 1:active))
    correct.selections.prop[i] = sum((s$selected %in% 1:active))/active
    pi[i] = pi.hat
    method[i] = ifelse(pi.hat < 0.5, "CPSS", "SS")
    if (i %% 10 == 0){print(i)}
  }
  out = data.frame(n.selected,false.selections, correct.selections.prop, method, pi, EV, snr)
}