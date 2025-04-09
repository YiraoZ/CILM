library(Rcpp)
library(ggplot2)
library(ggsci)
library(truncnorm)
library(gridExtra)
library(tidyr)
library(dplyr)
library(rstan)
library(latex2exp)
library(MASS)
library(coda)
sourceCpp('./PPD.cpp')
#Posterior Predictive Distribution(PPD)
PPD<-function(n_row, t_start, t_end, chain_len, N){#N is the number of repetition 
  spatials = c("homogeneous", "heterogeneous_small_K=3", "heterogeneous_small_K=5","heterogeneous_small_K=8","heterogeneous_large_K=3","heterogeneous_large_K=5", "heterogeneous_large_K=8")
  data_prefix<-"./simulation_data/"
  samples_basicILM_prefix<-"./MCMC/basic_ILM/"
  samples_M2_prefix<-"./MCMC/M2/"
  samples_M3_prefix<-"./MCMC/M3/"
  samples_M4_prefix<-"./MCMC/M4/"
  ppd_basicILM_prefix<-"./PPD/basic_ILM/"
  ppd_M2_prefix<-"./PPD/M2/"
  ppd_M3_prefix<-"./PPD/M3/"
  ppd_M4_prefix<-"./PPD/M4/"
  for(n in 1:N){
    for(k in 1:7){
      data<-read.csv(paste0(data_prefix,spatials[k],"/",spatials[k],"_",toString(n),".csv"))
      samples_basicILM<-read.csv(paste0(samples_basicILM_prefix,spatials[k],"/",spatials[k],"_",toString(n),".csv"))
      samples_M2<-read.csv(paste0(samples_M2_prefix,spatials[k],"/",spatials[k],"_",toString(n),".csv"))
      samples_M3<-read.csv(paste(samples_M3_prefix,spatials[k],"/",spatials[k],"_",toString(n),".csv"))
      samples_M4<-read.csv(paste0(samples_M4_prefix,spatials[k],"/",spatials[k],"_",toString(n),".csv"))
      
      ppd_basicILM <- as.data.frame(matrix(NA, nrow = n_row, ncol = nrow(data)))
      ppd_M2 <- as.data.frame(matrix(NA, nrow = n_row, ncol = nrow(data)))
      ppd_M3 <- as.data.frame(matrix(NA, nrow = n_row, ncol = nrow(data)))
      ppd_M4 <- as.data.frame(matrix(NA, nrow = n_row, ncol = nrow(data)))
      
      for(row in 1:n_row){
        s<-sample(1:(chain_len/2),1)
        ppd_basicILM[row,]<-PPD_StandardILM(data$time_of_infected,data$position_x, data$position_y, t_start, t_end, as.numeric(samples_basic[s,]))
        ppd_M2[row,]<-PPD_M2(data$time_of_infected,data$position_x, data$position_y, t_start, t_end, as.numeric(data$modes), as.numeric(samples_M2[s,]))
        ppd_M3[row,]<-PPD_M3(data$time_of_infected,data$position_x, data$position_y, t_start, t_end, as.numeric(data$modes), as.numeric(samples_M3[s,]))
        ppd_M4[row,]<-PPD_M4(data$time_of_infected,data$position_x, data$position_y, t_start, t_end, as.numeric(data$modes), as.numeric(samples_M4[s,]))
      }
      write.csv(ppd_basicILM,paste0(ppd_basicILM_prefix,spatials[k],"/start_from_",toString(t_start),"/",spatials[k],"_",toString(n),".csv"), row.names = FALSE)
      write.csv(ppd_M2,paste0(ppd_M2_prefix,spatials[k],"/start_from_",toString(t_start),"/",spatials[k],"_",toString(n),".csv"), row.names = FALSE)
      write.csv(ppd_M3,paste0(ppd_M3_prefix,spatials[k],"/start_from_",toString(t_start),"/",spatials[k],"_",toString(n),".csv"), row.names = FALSE)
      write.csv(ppd_M4,paste0(ppd_M4_prefix,spatials[k],"/start_from_",toString(t_start),"/",spatials[k],"_",toString(n),".csv"), row.names = FALSE)
    }
  }  
}