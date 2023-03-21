#Plot the true growth over an interval against the modelled growth
delta_S_dot <- function(data1, data2, filename){
  plot <- ggplot_delta_S_dot(data1, data2)
  ggsave(filename, plot=plot, scale=1, width=100, height=100, units="mm")
  return(plot)
}

#Plot function for delta_S_dot()
ggplot_delta_S_dot <- function(data1, data2){
  plot1 <- ggplot(data=data1, aes(x,y, col=Model))+
    geom_point(aes(x=delta_S, y=G_hat, shape=Model), size=2) +
    scale_shape_manual(values = c(3, 2)) +
    scale_color_manual(values = c("blue", "red")) +
    xlab("True growth cm/yr") +
    ylab("Estimated growth cm/yr") +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) + 
    geom_abline(intercept=0, slope=1, linetype="dashed") +
    xlim(0, NA) +
    theme_classic() +
    theme(legend.position = c(.3, .8))
  
  plot2 <- plot1 + 
    geom_text(data = data2,
              aes(x=delta_S, y=G_hat, label = labs), 
              color = "blue", nudge_x = 0.1)
  
  return(plot2)
}
