plot_elbow = function(n, p, active, snr = 3, LQ = F, upto = min(p, 30), true){
  d = gendata(n = n, p = p, active = active)
  signal = sqrt(mean((as.matrix(d$X) %*% as.matrix(d$beta))^2))
  sigma = as.numeric(signal/sqrt(snr))
  
  # Compute Y with SNR
  Y = as.matrix(d$X)%*%as.matrix(d$beta) + rnorm(nrow(d$X), 0, sd = sigma) 
  s = stabs::stabsel(x = d$X, y = Y, B = 100,
                     fitfun = stabs::lars.lasso, PFER = 5, cutoff = 0.75,
                     sampling.type = "MB")
  idx = sample(1:nrow(d$X), replace = F)
  rX = d$X[idx,]
  idxPushed = c(tail(idx, 1), head(idx, -1))
  rY = Y[idxPushed] |> as.matrix(ncol = 1)
  
  # Exclusion Probability Threshold
  sMix = stabs::stabsel(x = rX, y = rY, B = 100,
                        fitfun = stabs::lars.lasso, PFER = 5, cutoff = 0.75,
                        sampling.type = "MB")
  sMix_prob = sort(sMix$max, decreasing = T)
  mix_exclusion = quantile(sMix_prob, 0.95)
  EATS =  convert(s)[1:(length(convert(s)))][convert(s)[1:(length(convert(s)))] >= 100*mix_exclusion] 
  EATS = EATS |> getR()
  
  df1 = data.frame(sel.prob = s$max, vnames = factor(names(s$max), levels = names(s$max)))
  
  lq1 = data.frame(LQ = getLQ(convert(s)), vnames = factor(names(s$max), levels = names(s$max)))
  p1 = ggplot(lq1[1:upto,], aes(x = vnames, y = LQ)) + geom_line(group = 1) + theme_few_grid(base_size = 20) +
    labs(x = "Variable Name", y = "Likelihood") + geom_vline(xintercept = which.max(lq1$LQ), linetype = "dashed") +
    theme(axis.text.x = element_blank())
  
  
  p2 = ggplot(df1[1:upto, ], aes(x = reorder(vnames, sel.prob, decreasing = TRUE), y = sel.prob)) +
    geom_line(group = 1, color = "black") +
    geom_point(aes(shape = sel.prob >= sort(sel.prob, decreasing = TRUE)[EATS],
                   colour = vnames %in% names(d$beta)[1:active],
                   fill = vnames %in% names(d$beta)[1:active]), size = 3, stroke = 1) +
    geom_vline(xintercept = EATS, linetype = "dashed") +
    labs(x = "Variable Name", y = "Selection Probability") +
    scale_color_manual(name = "Truth",
                       labels = c("TRUE" = "Signal", "FALSE" = "Noise"),
                       values = c("TRUE" = "red", "FALSE" = "black")) +
    scale_fill_manual(name = "Truth",
                      labels = c("TRUE" = "Signal", "FALSE" = "Noise"),
                      values = c("TRUE" = "red", "FALSE" = "black")) +
    scale_shape_manual(name = "EATS-Selected",
                       labels = c("TRUE" = "Yes", "FALSE" = "No"),
                       values = c("TRUE" = 1, "FALSE" = 13)) + 
    theme_few_grid(base_size = 20) +
    theme(axis.text.x = element_blank(), legend.position = "none")
  
  list(p1,p2)
}

n = 20; p = 1000; active = 2
e1 = plot_elbow(n,p,active)

n = 100; p = 500; active = 10
e2 = plot_elbow(n,p,active)

n = 200; p = 200; active = 20
e3 = plot_elbow(n,p,active)

n = 500; p = 100; active = 20
e4 = plot_elbow(n,p,active)


e1[[2]] = e1[[2]] + labs(tag = "(I)")+ theme(plot.tag = element_text(size = 23))
e1[[1]] = e1[[1]] + labs(tag = "(I)")+ theme(plot.tag = element_text(size = 23))
e2[[2]] = e2[[2]] + labs(tag = "(II)")+ theme(plot.tag = element_text(size = 23))
e2[[1]] = e2[[1]] + labs(tag = "(II)")+ theme(plot.tag = element_text(size = 23))
e3[[2]] = e3[[2]] + labs(tag = "(III)")+ theme(plot.tag = element_text(size = 23))
e3[[1]] = e3[[1]] + labs(tag = "(III)")+ theme(plot.tag = element_text(size = 23))
e4[[2]] = e4[[2]] + labs(tag = "(IV)")+ theme(plot.tag = element_text(size = 23))
e4[[1]] = e4[[1]] + labs(tag = "(IV)")+ theme(plot.tag = element_text(size = 23))

pw1 = (e1[[2]] / e1[[1]])
pw2 = (e2[[2]] / e2[[1]])
pw3 = (e3[[2]] / e3[[1]])
pw4 = (e4[[2]] / e4[[1]])
pw = pw1 | pw2 | pw3 | pw4
pw 


