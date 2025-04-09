#Figure 5 in Appendix - epicurves under seven scenarios
library(Rcpp)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(ggsci)
#set the file path
Visulization<-function(n){ #n is the # of population to plot
  datapath_prefix<-"./simulation_data/"
  spatials = c("homogeneous", "heterogeneous_small_K=3", "heterogeneous_small_K=5","heterogeneous_small_K=8","heterogeneous_large_K=3","heterogeneous_large_K=5", "heterogeneous_large_K=8")
  titles = c("CSR", "Low-variance Clustered, K=3", "Low-variance Clustered, K=5", "Low-variance Clustered, K=8", "High-variance Clustered, K=3", "High-variance Clustered, K=5", "High-variance Clustered, K=8")
  clusters = c("True Labels", "Spatial Clustering", "Spatio-temporal Clustering")
  plot_list <-vector('list', 21) #7*3
  # traceplot_list <-vector('list', 2)
  for(k in 1:length(spatials)){
    #import data
    datapath <- paste(datapath_prefix,spatials[k],"/",spatials[k],"_",toString(n),".csv",sep = "")
    data<-read.csv(datapath)
    
   #visualize clustering result
    plot_list[[(k-1)*3+1]] <- ggplot(data, aes(x = position_x, y = position_y, color=factor(cluster_id)))+
      geom_point(size = 5, alpha = 0.8)+
      scale_color_nejm()+
      theme(legend.position = "none", axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), panel.background = element_rect(fill = "white", color = "black"),
            panel.border = element_rect(color = "black", fill = NA))+
      ggtitle(paste(titles[k],clusters[1], sep = ", "))
    plot_list[[(k-1)*3+2]] <- ggplot(data, aes(x = position_x, y = position_y, color=factor(xymodes)))+
      geom_point(size = 5, alpha = 0.8)+
      scale_color_nejm()+
      theme(legend.position = "none", axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), panel.background = element_rect(fill = "white", color = "black"),
            panel.border = element_rect(color = "black", fill = NA))+
      ggtitle(paste(titles[k],clusters[2], sep = ", "))
    plot_list[[(k-1)*3+3]] <- ggplot(data, aes(x = position_x, y = position_y, color=factor(xytmodes)))+
      geom_point(size = 5, alpha = 0.8)+
      scale_color_nejm()+
      theme(legend.position = "none", axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), panel.background = element_rect(fill = "white", color = "black"),
            panel.border = element_rect(color = "black", fill = NA))+
      ggtitle(paste(titles[k],clusters[3], sep = ", "))
    #traceplots of variance parameters w_x, w_y, phi
    # tracedata_xy <- data.frame(
    #   Iteration = 1:chain_len,
    #   w_x = samples_xy[1:chain_len, 1],
    #   w_y = samples_xy[1:chain_len, 2]
    # )
    # tracedata_singlephi <- data.frame(
    #   Iteration = 1:chain_len,
    #   w_x = samples_singlephi[1:chain_len, 1],
    #   w_y = samples_singlephi[1:chain_len, 2],
    #   phi = samples_singlephi[1:chain_len, 3]
    # )
    # traceplot_list[[1]]<- ggplot() +
    #   geom_line(data = tracedata_xy, aes(x = Iteration, y = w_x), color = "red") +
    #   geom_line(data = tracedata_xy, aes(x = Iteration, y = w_y), color = "blue") +
    #   labs(title = "Traceplot", x = "Iteration", y = "Parameter Value") +
    #   theme_minimal(base_size = 15) +
    #   theme(
    #     plot.title = element_text(hjust = 0.5), 
    #     panel.grid = element_blank(), 
    #     panel.background = element_blank(), 
    #     axis.line = element_line(color = "black") 
    #   )
    # traceplot_list[[2]]<- ggplot() +
    #   geom_line(data = tracedata_singlephi, aes(x = Iteration, y = w_x), color = "red") +
    #   geom_line(data = tracedata_singlephi, aes(x = Iteration, y = w_y), color = "blue") +
    #   geom_line(data = tracedata_singlephi, aes(x = Iteration, y = phi), color = "green")
    # labs(title = "Traceplot", x = "Iteration", y = "Parameter Value") +
    #   theme_minimal(base_size = 15) +
    #   theme(
    #     plot.title = element_text(hjust = 0.5), 
    #     panel.grid = element_blank(), 
    #     panel.background = element_blank(), 
    #     axis.line = element_line(color = "black") )
    # combined_traceplots<-grid.arrange(grobs = traceplot_list, nrow = 1, ncol = 2, heights = c(3), weights = c(3, 3))    # ggsave(paste(outpath_prefix,spatials[k],"/",spatials[k],"_",toString(n),".png",sep = ""), combined_traceplots, units = "in", width = 12.8, height = 9.1)
  }
  combined_plots<-grid.arrange(grobs = plot_list, nrow = 7, ncol = 3)
  ggsave(paste("./Plots/Clustering_",toString(n),".pdf",sep = ""), combined_plots, units = "in", width = 19.2, height = 33.6)
}
for(n in 1:10){
  Visulization(n)
}