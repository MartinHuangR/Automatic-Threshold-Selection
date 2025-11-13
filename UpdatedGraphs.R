NtotMEANplot = function(TOT, lim = max(TOT$NN)){
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
  
  TOTm = TOT |> group_by(Method, SNR, Dimension) |> 
    summarise(mean = mean(NN), .groups = "drop")
  
  ggplot(TOTm, aes(x = Method, y = mean, group = SNR)) +
    geom_line(linetype = "dotted") + 
    ylab("Variables Selected") + 
    geom_hline(data = TOTdummy, aes(yintercept = active), linetype = "dashed", linewidth = 0.5) + 
    geom_point(aes(colour = Method, shape = SNR, fill = Method), size = 3, show.legend = TRUE) +  
    facet_grid(~Dimension, labeller = label_parsed) + 
    theme_few_grid(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
      legend.title = element_blank()) +
    scale_shape_manual(
      values = c(21,22,23,24), 
      labels = parse(text = c("SNR==0.5", "SNR==1", "SNR==2", "SNR==3"))  # parsed SNR labels
    ) +
    scale_colour_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"), guide = "none") +  # hide colour legend
    scale_fill_manual(values = c("#FC8D62", "#FFD92F","#A6D854","#A6D854","#A6D854","#8DA0CB","#8DA0CB","#8DA0CB"), guide = "none")  # hide fill legend
}

TOT = rbind(c11Nhard, c22Nhard)


load("Data/S1.Rdata")
load("Data/S2.Rdata")
load("Data/S3.Rdata")
load("Data/S4.Rdata")