#Plot parameter estimation result
library(Rcpp)
library(ggplot2)
library(ggsci)
library(truncnorm)
library(gridExtra)
library(tidyr)
library(dplyr)
library(rstan)
library(latex2exp)
#Parameter Estimation plot
##data combination
#set the file path
data_combined_xyt<-function(spatialtype, modeltype, chain_len, N){
  if(modeltype == "StandardILM"){
    combined_data_a0 <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_a1 <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_beta <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    for (i in 1:N) {
      file_name <- paste0("./MCMC_CILM/StandardILM/", spatialtype,"/",spatialtype,"_",toString(i), ".csv", sep = "")
      data <- read.csv(file_name)
      combined_data_a0[[i]] <- data$a0
      combined_data_beta[[i]] <- data$beta
    }
    combined_data <- list(combined_data_a0, combined_data_beta)
  }
  if(modeltype == "M2"||modeltype == "M3"){
    combined_data_a0 <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_a1 <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_beta <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_tbeta <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    for (i in 1:N) {
      file_name <- paste0("./MCMC_CILM/", modeltype, "/", spatialtype,"/",spatialtype,"_singlephi_",toString(i), ".csv", sep = "")
      data <- read.csv(file_name)
      combined_data_a0[[i]] <- data$a0
      combined_data_beta[[i]] <- data$beta
      combined_data_tbeta[[i]] <- data$tbeta
    }
    combined_data <- list(combined_data_a0, combined_data_beta, combined_data_tbeta)
  }
  if(modeltype == "M4"){
    combined_data_a0 <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_beta <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_tbeta <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_delta <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    for (i in 1:N) {
      file_name <- paste0("./MCMC_CILM/M4/", spatialtype,"/",spatialtype,"_singlephi_",toString(i), ".csv", sep = "")
      data <- read.csv(file_name)
      combined_data_a0[[i]] <- data$a0
      combined_data_beta[[i]] <- data$beta
      combined_data_tbeta[[i]] <- data$tbeta
      combined_data_delta[[i]] <- data$delta
    }
    combined_data <- list(combined_data_a0, combined_data_beta, combined_data_tbeta, combined_data_delta)
  }
  return(combined_data)
}
##data summary
data_summary<-function(data, N){
  summary_stats <- apply(data, 2, function(x) {
    median <- median(x)
    lower_ci <- quantile(x, 0.025)
    upper_ci <- quantile(x, 0.975)
    return(c(median, lower_ci, upper_ci))
  })
  summary_stats <- t(summary_stats)
  summary_stats<-cbind(summary_stats, 1:N)
  colnames(summary_stats) <- c("median", "lower_ci", "upper_ci", "population")
  summary_stats<-as.data.frame(summary_stats)
  return(summary_stats)
}
##plot parameter estimation for basic ILM, M2, M3, M4
samples_plot_xyt<-function(chain_len, N){
  spatials = c("homogeneous", "heterogeneous_small_K=3", "heterogeneous_small_K=5","heterogeneous_small_K=8","heterogeneous_large_K=3","heterogeneous_large_K=5", "heterogeneous_large_K=8")
  titles = c("CSR", "Low-variance Clustered, K=3", "Low-variance Clustered, K=5", "Low-variance Clustered, K=8", "High-variance Clustered, K=3", "High-variance Clustered, K=5", "High-variance Clustered, K=8")
  models = c("StandardILM","M2","M3", "M4")
  truevalues = c(0.8, 2, NULL, NULL)
  parameters = c(", $\\alpha$", ", $\\beta$", ", $\\tilde{\\beta}$", ", $\\delta$")
  ranges = list(c(0,3.5), c(0, 4.0), c(0, 3.6), c(0, 10))
  for(m in 1:4){
    if(models[m] == "StandardILM"){
      plot_list <-vector('list', 14)#7*2
    }
    if(models[m] == "M2"||models[m] == "M3"){
      plot_list <-vector('list', 21)#7*3
    }
    if(models[m] == "M4"){
      plot_list <-vector('list', 28)#7*4
    }
    for(k in 1:7){
      combined_data_xyt<-data_combined_xyt(spatials[k], models[m], chain_len, N)
      for(c in 1:length(combined_data_xyt)){
        title = TeX(paste0(titles[k], parameters[c]))
        summary_stats_xyt<-data_summary(combined_data_xyt[[c]], N)
        plot<-ggplot(summary_stats_xyt, aes(x = population, y = median)) +
          geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.3, size =0.7) +
          geom_point(size = 4, shape = 16, color = "blue") +
          ggtitle(title)+
          theme(
            axis.title.x = element_text(size = 15), 
            axis.title.y = element_text(size = 15),
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
            panel.background = element_rect(fill = "white", color = "black"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
          ) +geom_hline(yintercept = truevalues[c], linetype = "dashed", color = "red", size = 0.7)+
          ylim(ranges[[c]])+
          labs(x = "Population", y = "Value")+scale_x_continuous(breaks = 1:15)
        plot_list[[(k-1)*length(combined_data_xyt)+c]]<-plot
      }
    }
    filepath_mcmcplot = paste("./Plots/", models[m],"_xyt_CIplot.pdf", sep = "")
    combined_plots<-grid.arrange(grobs = plot_list, nrow = 7, ncol = length(combined_data_xyt))
    ggsave(filepath_mcmcplot, combined_plots, units = "in", height = 33.6, width = 6.4*length(combined_data_xyt))
  }
}
samples_plot_xyt(1000,10)


data_combined_kt10<-function(spatialtype, modeltype, chain_len, N){
  if(modeltype == "basicILM"){
    combined_data_a0 <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_a1 <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_beta <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    for (i in 1:N) {
      file_name <- paste0("./MCMC/basicILM/", spatialtype,"/",spatialtype,"_",toString(i), ".csv", sep = "")
      data <- read.csv(file_name)
      combined_data_a0[[i]] <- data$a0
      combined_data_beta[[i]] <- data$beta
    }
    combined_data <- list(combined_data_a0, combined_data_beta)
  }
  if(modeltype == "M2"||modeltype == "M3"){
    combined_data_a0 <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_a1 <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_beta <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_tbeta <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    for (i in 1:N) {
      file_name <- paste0("./MCMC/", modeltype, "/", spatialtype,"/",spatialtype,"_kt10_",toString(i), ".csv", sep = "")
      data <- read.csv(file_name)
      combined_data_a0[[i]] <- data$a0
      combined_data_beta[[i]] <- data$beta
      combined_data_tbeta[[i]] <- data$tbeta
    }
    combined_data <- list(combined_data_a0, combined_data_beta, combined_data_tbeta)
  }
  if(modeltype == "M4"){
    combined_data_a0 <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_beta <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_tbeta <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    combined_data_delta <- data.frame(matrix(nrow=chain_len/2, ncol=N))
    for (i in 1:N) {
      file_name <- paste0("./MCMC/M4/", spatialtype,"/",spatialtype,"_xyt_",toString(i), ".csv", sep = "")
      data <- read.csv(file_name)
      combined_data_a0[[i]] <- data$a0
      combined_data_beta[[i]] <- data$beta
      combined_data_tbeta[[i]] <- data$tbeta
      combined_data_delta[[i]] <- data$delta
    }
    combined_data <- list(combined_data_a0, combined_data_beta, combined_data_tbeta, combined_data_delta)
  }
  return(combined_data)
}
samples_plot_kt10<-function(chain_len, N){
  spatials = c("homogeneous", "heterogeneous_small_K=3", "heterogeneous_small_K=5","heterogeneous_small_K=8","heterogeneous_large_K=3","heterogeneous_large_K=5", "heterogeneous_large_K=8")
  titles = c("CSR", "Low-variance Clustered, K=3", "Low-variance Clustered, K=5", "Low-variance Clustered, K=8", "High-variance Clustered, K=3", "High-variance Clustered, K=5", "High-variance Clustered, K=8")
  models = c("StandardILM","M2","M3", "M4")
  truevalues = c(0.8, 2, NULL, NULL)
  parameters = c(", $\\alpha$", ", $\\beta$", ", $\\tilde{\\beta}$", ", $\\delta$")
  for(m in 1:4){
    if(models[m] == "basicILM"){
      plot_list <-vector('list', 14)#7*2
    }
    if(models[m] == "M2"||models[m] == "M3"){
      plot_list <-vector('list', 21)#7*3
    }
    if(models[m] == "M4"){
      plot_list <-vector('list', 28)#7*4
    }
    for(k in 1:7){
      combined_data_kt10<-data_combined_kt10(spatials[k], models[m], chain_len, N)
      for(c in 1:length(combined_data_kt10)){
        title = TeX(paste0(titles[k], parameters[c]))
        summary_stats_kt10<-data_summary(combined_data_kt10[[c]], N)
        plot<-ggplot(summary_stats_kt10, aes(x = population, y = median)) +
          geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.1) +
          geom_point(size = 3, shape = 16, color = "blue") +
          ggtitle(title)+
          theme(
            axis.title.x = element_text(size = 15), 
            axis.title.y = element_text(size = 15),
            plot.title = element_text(hjust = 0.5, size = 18),
            panel.background = element_rect(fill = "white", color = "black"),
            panel.border = element_rect(color = "black", fill = NA)
          ) +geom_hline(yintercept = truevalues[c], linetype = "dashed", color = "red")+
          labs(x = "Population", y = "Value")+scale_x_continuous(breaks = 1:15)
        plot_list[[(k-1)*length(combined_data_kt10)+c]]<-plot
      }
    }
    filepath_mcmcplot = paste("./MCMC/", models[m], "/",models[m],"_kt10_CIplot.png", sep = "")
    combined_plots<-grid.arrange(grobs = plot_list, nrow = 7, ncol = length(combined_data_kt10))
    ggsave(filepath_mcmcplot, combined_plots, units = "in", height = 33.6, width = 6.4*length(combined_data_kt10))
  }
}

