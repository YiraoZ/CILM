#Plot Posterior Predictive Distribution
library(Rcpp)
library(ggplot2)
library(ggsci)
library(truncnorm)
library(gridExtra)
library(tidyr)
library(dplyr)
library(rstan)
library(latex2exp)
library(coda)
#PPD plot
#Plot PPD for all spatials(but only one dataset), and all models
PPD_plot_xyt<-function(t_start, t_end){
  spatials = c("homogeneous", "heterogeneous_small_K=3", "heterogeneous_small_K=5","heterogeneous_small_K=8","heterogeneous_large_K=3","heterogeneous_large_K=5", "heterogeneous_large_K=8")
  titles = c("CSR", "Low-variance Clustered, K=3", "Low-variance Clustered, K=5", "Low-variance Clustered, K=8", "High-variance Clustered, K=3", "High-variance Clustered, K=5", "High-variance Clustered, K=8")
  models = c("StandardILM","M2","M3","M4")
  models_rename = c("basic", "M2", "M3", "M4")
  plot_list <-vector('list', 28)#7*4
  start_from = paste0("start_from_",toString(t_start))
  tmax<-20
  for(k in 1:7){
    title = paste0(titles[k],", ",models_rename[1])
    filepath_epidata = paste0("./simulation_data/", spatials[k],"/",spatials[k],"_1.csv")
    filepath_PPD = paste0("./PPD/StandardILM/",spatials[k],"/",start_from, "/", spatials[k],"_1.csv")
    data<-read.csv(filepath_PPD)[-1,]
    data_true<-read.csv(filepath_epidata)
    df<-as.data.frame(0:tmax)
    colnames(df)<-c("time_of_infected")
    for(i in 1:100){
      inf_time <- as.data.frame(t(data[i,]),row.names = NULL)
      colnames(inf_time) <- c("time_of_infected")
      daily_cases <- df%>%left_join(inf_time %>%
                                      filter(time_of_infected != -1) %>%
                                      group_by(time_of_infected) %>%
                                      summarize(Cases = n()), by = c("time_of_infected"="time_of_infected"))%>%
        replace_na(list(Cases = 0))
      temp_df<-data.frame(x=daily_cases$time_of_infected, y=daily_cases$Cases)
      colname <-paste0("epicurve_",i)
      df[[colname]]<-temp_df$y
    }
    daily_cases <- df%>%left_join(data_true %>%
                                    filter(time_of_infected != -1) %>%
                                    group_by(time_of_infected) %>%
                                    summarize(Cases = n()), by = c("time_of_infected"="time_of_infected"))%>%
      replace_na(list(Cases = 0))
    temp_df<-data.frame(x=daily_cases$time_of_infected, y=daily_cases$Cases)
    df[["true"]]<-temp_df$y
    #calculate  medians, 95%CI
    medians <- apply(df[, -1], 1, median)
    HPDI <- t(apply(df[, -1], 1, function(row) {
      quantiles <- HPDinterval(as.mcmc(row), prob = 0.95)
      HPDI_low <- quantiles[1]
      HPDI_high <- quantiles[2]
      return(c(HPDI_low, HPDI_high))
    }))
    df[["median"]]<-medians
    df[["hpdi_low"]]<-HPDI[,1]
    df[["hpdi_high"]]<-HPDI[,2]
    plot_list[[(k-1)*4+1]]<-ggplot(df, aes(x = 1:(tmax+1))) + 
      geom_ribbon(aes(ymin = hpdi_low, ymax = hpdi_high), fill = "lightblue", alpha = 0.5) +  
      geom_line(aes(y = true), size = 1.2, color = "black") +  
      geom_line(aes(y = median), linetype = "dashed", size = 1, color = "navy") +  
      labs(title = title, x = "Time", y = "Incidence") +  
      theme(
        text = element_text(size = 15), 
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(color = "black", fill = NA)
      )    
    
    #plot_list[[(k-1)*4+1]]<-recordPlot()
    for(m in 2:4){
      title = paste0(titles[k],", ",models_rename[m])
      filepath_epidata = paste0("./simulation_data/", spatials[k],"/",spatials[k],"_1.csv")
      filepath_PPD = paste0("./PPD/", models[m],"/",spatials[k],"/",start_from, "/", spatials[k],"_singlephi_1.csv")
      data<-read.csv(filepath_PPD)[-1,]
      data_true<-read.csv(filepath_epidata)
      df<-as.data.frame(0:20)
      colnames(df)<-c("time_of_infected")
      for(i in 1:100){
        inf_time <- as.data.frame(t(data[i,]),row.names = NULL)
        colnames(inf_time) <- c("time_of_infected")
        daily_cases <- df%>%left_join(inf_time %>%
                                        filter(time_of_infected != -1) %>%
                                        group_by(time_of_infected) %>%
                                        summarize(Cases = n()), by = c("time_of_infected"="time_of_infected"))%>%
          replace_na(list(Cases = 0))
        temp_df<-data.frame(x=daily_cases$time_of_infected, y=daily_cases$Cases)
        colname <-paste0("epicurve_",i)
        df[[colname]]<-temp_df$y
      }
      daily_cases <- df%>%left_join(data_true %>%
                                      filter(time_of_infected != -1) %>%
                                      group_by(time_of_infected) %>%
                                      summarize(Cases = n()), by = c("time_of_infected"="time_of_infected"))%>%
        replace_na(list(Cases = 0))
      temp_df<-data.frame(x=daily_cases$time_of_infected, y=daily_cases$Cases)
      df[["true"]]<-temp_df$y
      #calculate  medians, 95%CI
      medians <- apply(df[, -1], 1, median)
      HPDI <- t(apply(df[, -1], 1, function(row) {
        quantiles <- HPDinterval(as.mcmc(row), prob = 0.95)
        HPDI_low <- quantiles[1]
        HPDI_high <- quantiles[2]
        return(c(HPDI_low, HPDI_high))
      }))
      df[["median"]]<-medians
      df[["hpdi_low"]]<-HPDI[,1]
      df[["hpdi_high"]]<-HPDI[,2]
      plot_list[[(k-1)*4+m]]<-ggplot(df, aes(x = 1:(tmax+1))) + 
        geom_ribbon(aes(ymin = hpdi_low, ymax = hpdi_high), fill = "lightblue", alpha = 0.5) +  
        geom_line(aes(y = true), size = 1.2, color = "black") +  
        geom_line(aes(y = median), linetype = "dashed", size = 1, color = "navy") +  
        labs(title = title, x = "Time", y = "Incidence") +  
        theme(
          text = element_text(size = 15), 
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.border = element_rect(color = "black", fill = NA)
        )    
    }
  }
  filepath_PPD = paste("./Plots/",start_from,"_xyt.pdf", sep = "")
  combined_plots<-grid.arrange(grobs = plot_list, nrow = 7, ncol = 4)
  ggsave(filepath_PPD, combined_plots, units = "in", height = 33.6, width = 25.6)
}

PPD_plot_xyt(1, 31)
PPD_plot_xyt(5, 31)
