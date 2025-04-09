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
#set the file path
#Parameter Estimation using MCMC with Stan
##import the Stan files
basic_model<-stan_model(file = "./BasicILM.stan")
M2_model<-stan_model(file = "./CILM2.stan")
M3_model<-stan_model(file = "./CILM3.stan")
M4_model<-stan_model(file = "./CILM4.stan")
Basic_MCMC<-function(d, t_inf, filepath){
  stan_data <- list(N=100,
                    t_start=1,
                    t_end=31,
                    period=3,
                    d=d,
                    t_inf=t_inf
  )
  initf1 <- function() {
    list(a0 = 1, beta =2)}
  fit<-sampling(basic_model, stan_data, init = initf1, iter = 1000, chains = 4)
  print(summary(fit, pars=c('a0'))$summary)
  print(summary(fit, pars=c('beta'),)$summary)
  samples <-rstan::extract(fit)
  samples <-data.frame(a0 = samples$a0, beta = samples$beta)
  write.csv(samples, filepath, row.names = FALSE)
}
M2_MCMC<-function(d, K, t_inf, c, cluster_id, filepath){
  stan_data <- list(N=100,
                    K=K,
                    t_start=1,
                    t_end=31,
                    period=3,
                    d=d,
                    centroids = c,
                    t_inf=t_inf,
                    cluster_id = cluster_id
  )
  initf1 <- function() {
    list(a0 = 1, beta =2, tbeta = 2)
  }
  fit<-sampling(M2_model, stan_data, init = initf1, iter = 1000, chains = 4)
  print(summary(fit, pars=c('a0'))$summary)
  print(summary(fit, pars=c('beta'),)$summary)
  print(summary(fit, pars=c('tbeta'),)$summary)
  samples <-rstan::extract(fit)
  samples <-data.frame(a0 = samples$a0, beta = samples$beta, tbeta = samples$tbeta)
  write.csv(samples, filepath, row.names = FALSE)
}
M3_MCMC<-function(d, K, t_inf,cluster_id, filepath){
  stan_data <- list(N=100,
                    K=K,
                    t_start=1,
                    t_end=31,
                    period=3,
                    d=d,
                    t_inf=t_inf,
                    cluster_id = cluster_id
  )
  initf1 <- function() {
    list(a0 = 1, beta =2, tbeta = 2)
  }
  fit<-sampling(M3_model, stan_data, init = initf1, iter = 1000, chains = 4)
  print(summary(fit, pars=c('a0'))$summary)
  print(summary(fit, pars=c('beta'),)$summary)
  print(summary(fit, pars=c('tbeta'),)$summary)
  samples <-rstan::extract(fit)
  samples <-data.frame(a0 = samples$a0, beta = samples$beta, tbeta = samples$tbeta)
  write.csv(samples, filepath, row.names = FALSE)
}
M4_MCMC<-function(d, K, t_inf, c, cluster_id, filepath){
  stan_data <- list(N=100,
                    K=K,
                    t_start=1,
                    t_end=31,
                    period=3,
                    d=d,
                    centroids = c,
                    t_inf=t_inf,
                    cluster_id = cluster_id
  )
  initf1 <- function() {
    list(a0 = 1, beta =2, tbeta = 2, delta = 1)
  }
  fit<-sampling(M4_model, stan_data, init = initf1, iter = 1000, chains = 4)
  print(summary(fit, pars=c('a0'))$summary)
  print(summary(fit, pars=c('beta'),)$summary)
  print(summary(fit, pars=c('tbeta'),)$summary)
  print(summary(fit, pars=c('delta'),)$summary)
  samples <-rstan::extract(fit)
  samples <-data.frame(a0 = samples$a0, beta = samples$beta, tbeta = samples$tbeta, delta = samples$delta)
  write.csv(samples, filepath, row.names = FALSE)
}
##run basic ILM, M2, M3, M4 for all datasets
basic_ILM<-function(N){
  datasets <- c("homogeneous", "heterogeneous_small_K=3", "heterogeneous_small_K=5", "heterogeneous_small_K=8", "heterogeneous_large_K=3", "heterogeneous_large_K=5", "heterogeneous_large_K=8")
  for(k in 1:7){
    for(n in 1:N){
      data_path = paste("./simulation_data/",datasets[k],"/",datasets[k],"_",toString(n),".csv",sep = "")
      basic_path = paste("./MCMC/basicILM/", datasets[k],"/",datasets[k],"_", toString(n), ".csv", sep = "")     
      data<-read.csv(data_path)
      d<-cbind(data$position_x, data$position_y)
      t_inf<-data$time_of_infected
      Basic_MCMC(d, t_inf, basic_path)
    }
  }
}
M2_CILM<-function(N){
  spatials <- c("homogeneous", "heterogeneous_small_K=3", "heterogeneous_small_K=5", "heterogeneous_small_K=8", "heterogeneous_large_K=3", "heterogeneous_large_K=5", "heterogeneous_large_K=8")
  for(s in 1:7){
    for(n in 1:N){
      data_path = paste("./simulation_data/",spatials[s],"/",spatials[s],"_",toString(n),".csv",sep = "")
      #M2 paths
      # M2_path_xy <- paste("./MCMC_CILM/M2/", spatials[s],"/",spatials[s],"_xy_", n, ".csv", sep = "")
      M2_path_xyt <- paste("./MCMC/M2/", spatials[s],"/",spatials[s],"_xyt_", n, ".csv", sep = "")
      M2_path_k3<-paste("./MCMC/M2/", spatials[s],"/",spatials[s],"_k",3,"_", n, ".csv", sep = "")
      M2_path_kt3<-paste("./MCMC/M2/", spatials[s],"/",spatials[s],"_kt",3,"_", n, ".csv", sep = "")
      M2_path_k5<-paste("./MCMC/M2/", spatials[s],"/",spatials[s],"_k",5,"_", n, ".csv", sep = "")
      M2_path_kt5<-paste("./MCMC/M2/", spatials[s],"/",spatials[s],"_kt",5,"_", n, ".csv", sep = "")
      M2_path_k8<-paste("./MCMC/M2/", spatials[s],"/",spatials[s],"_k",8,"_", n, ".csv", sep = "")
      M2_path_kt8<-paste("./MCMC/M2/", spatials[s],"/",spatials[s],"_kt",8,"_", n, ".csv", sep = "")
      M2_path_k10<-paste("./MCMC/M2/", spatials[s],"/",spatials[s],"_k",10,"_", n, ".csv", sep = "")
      M2_path_kt10<-paste("./MCMC/M2/",spatials[s],"/",spatials[s],"_kt",10,"_", n, ".csv", sep = "")
      
      data<-read.csv(data_path)
      d<-cbind(data$position_x, data$position_y)
      t_inf<-data$time_of_infected
      # data$xymodes<-as.numeric(factor(data$xymodes))
      data$xytmodes<-as.numeric(factor(data$xytmodes))
      data$k3<-as.numeric(factor(data$k3))
      data$kt3<-as.numeric(factor(data$kt3))
      data$k5<-as.numeric(factor(data$k5))
      data$kt5<-as.numeric(factor(data$kt5))
      data$k8<-as.numeric(factor(data$k8))
      data$kt8<-as.numeric(factor(data$kt8))
      data$k10<-as.numeric(factor(data$k10))
      data$kt10<-as.numeric(factor(data$kt10))
      
      # centroids_xy <- aggregate(d ~ xymodes, data = data, FUN = mean)
      # centroids_xy <- cbind(centroids_xy$V1, centroids_xy$V2)
      centroids_xyt <- aggregate(d ~ xytmodes, data = data, FUN = base::mean)
      centroids_xyt <- cbind(centroids_xyt$V1, centroids_xyt$V2)
      centroids_3 <- aggregate(d ~ k3, data = data, FUN = mean)
      centroids_3 <- cbind(centroids_3$V1, centroids_3$V2)
      centroids_t3 <- aggregate(d ~ kt3, data = data, FUN = base::mean)
      centroids_t3 <- cbind(centroids_t3$V1, centroids_t3$V2)
      centroids_5 <- aggregate(d ~ k5, data = data, FUN = mean)
      centroids_5 <- cbind(centroids_5$V1, centroids_5$V2)
      centroids_t5 <- aggregate(d ~ kt5, data = data, FUN = base::mean)
      centroids_t5 <- cbind(centroids_t5$V1, centroids_t5$V2)
      centroids_8 <- aggregate(d ~ k8, data = data, FUN = mean)
      centroids_8 <- cbind(centroids_8$V1, centroids_8$V2)
      centroids_t8 <- aggregate(d ~ kt8, data = data, FUN = base::mean)
      centroids_t8 <- cbind(centroids_t8$V1, centroids_t8$V2)
      centroids_10 <- aggregate(d ~ k10, data = data, FUN = mean)
      centroids_10 <- cbind(centroids_10$V1, centroids_10$V2)
      centroids_t10 <- aggregate(d ~ kt10, data = data, FUN = base::mean)
      centroids_t10 <- cbind(centroids_t10$V1, centroids_t10$V2)
      
      # K_xy = length(unique(data$xymodes))
      K_xyt = length(unique(data$xytmodes))
      
      # Basic_MCMC(d, t_inf, basic_path)
      # M2_MCMC(d, K_xy, t_inf, centroids_xy, data$xymodes, M2_path_xy)
      M2_MCMC(d, K_xyt, t_inf, centroids_xyt, data$xytmodes, M2_path_xyt)
      M2_MCMC(d, 3, t_inf, centroids_3, data$k3, M2_path_k3)
      M2_MCMC(d, 3, t_inf, centroids_t3, data$kt3, M2_path_kt3)
      M2_MCMC(d, 5, t_inf, centroids_5, data$k5, M2_path_k5)
      M2_MCMC(d, 5, t_inf, centroids_t5, data$kt5, M2_path_kt5)
      M2_MCMC(d, 8, t_inf, centroids_8, data$k8, M2_path_k8)
      M2_MCMC(d, 8, t_inf, centroids_t8, data$kt8, M2_path_kt8)
      M2_MCMC(d, 10, t_inf, centroids_10, data$k10, M2_path_k10)
      M2_MCMC(d, 10, t_inf, centroids_t10, data$kt10, M2_path_kt10)
    }
  }
  
}
M2_CILM(5)
M3_CILM<-function(N){
  spatials <- c("homogeneous", "heterogeneous_small_K=3", "heterogeneous_small_K=5", "heterogeneous_small_K=8", "heterogeneous_large_K=3", "heterogeneous_large_K=5", "heterogeneous_large_K=8")
  for(s in 1:7){
    for(n in 1:N){
      data_path = paste("./simulation_data/",spatials[s],"/",spatials[s],"_",toString(n),".csv",sep = "")
      #M2 paths
      # M2_path_xy <- paste("./MCMC_CILM/M2/", spatials[s],"/",spatials[s],"_xy_", n, ".csv", sep = "")
      M3_path_xyt <- paste("./MCMC/M3/", spatials[s],"/",spatials[s],"_xyt_", n, ".csv", sep = "")
      M3_path_k3<-paste("./MCMC/M3/", spatials[s],"/",spatials[s],"_k",3,"_", n, ".csv", sep = "")
      M3_path_kt3<-paste("./MCMC/M3/", spatials[s],"/",spatials[s],"_kt",3,"_", n, ".csv", sep = "")
      M3_path_k5<-paste("./MCMC/M3/", spatials[s],"/",spatials[s],"_k",5,"_", n, ".csv", sep = "")
      M3_path_kt5<-paste("./MCMC/M3/", spatials[s],"/",spatials[s],"_kt",5,"_", n, ".csv", sep = "")
      M3_path_k8<-paste("./MCMC/M3/", spatials[s],"/",spatials[s],"_k",8,"_", n, ".csv", sep = "")
      M3_path_kt8<-paste("./MCMC/M3/", spatials[s],"/",spatials[s],"_kt",8,"_", n, ".csv", sep = "")
      M3_path_k10<-paste("./MCMC/M3/", spatials[s],"/",spatials[s],"_k",10,"_", n, ".csv", sep = "")
      M3_path_kt10<-paste("./MCMC/M3/",spatials[s],"/",spatials[s],"_kt",10,"_", n, ".csv", sep = "")
      
      data<-read.csv(data_path)
      d<-cbind(data$position_x, data$position_y)
      t_inf<-data$time_of_infected
      # data$xymodes<-as.numeric(factor(data$xymodes))
      data$xytmodes<-as.numeric(factor(data$xytmodes))
      data$k3<-as.numeric(factor(data$k3))
      data$kt3<-as.numeric(factor(data$kt3))
      data$k5<-as.numeric(factor(data$k5))
      data$kt5<-as.numeric(factor(data$kt5))
      data$k8<-as.numeric(factor(data$k8))
      data$kt8<-as.numeric(factor(data$kt8))
      data$k10<-as.numeric(factor(data$k10))
      data$kt10<-as.numeric(factor(data$kt10))
      
      # centroids_xy <- aggregate(d ~ xymodes, data = data, FUN = mean)
      # centroids_xy <- cbind(centroids_xy$V1, centroids_xy$V2)
      centroids_xyt <- aggregate(d ~ xytmodes, data = data, FUN = base::mean)
      centroids_xyt <- cbind(centroids_xyt$V1, centroids_xyt$V2)
      centroids_3 <- aggregate(d ~ k3, data = data, FUN = mean)
      centroids_3 <- cbind(centroids_3$V1, centroids_3$V2)
      centroids_t3 <- aggregate(d ~ kt3, data = data, FUN = base::mean)
      centroids_t3 <- cbind(centroids_t3$V1, centroids_t3$V2)
      centroids_5 <- aggregate(d ~ k5, data = data, FUN = mean)
      centroids_5 <- cbind(centroids_5$V1, centroids_5$V2)
      centroids_t5 <- aggregate(d ~ kt5, data = data, FUN = base::mean)
      centroids_t5 <- cbind(centroids_t5$V1, centroids_t5$V2)
      centroids_8 <- aggregate(d ~ k8, data = data, FUN = mean)
      centroids_8 <- cbind(centroids_8$V1, centroids_8$V2)
      centroids_t8 <- aggregate(d ~ kt8, data = data, FUN = base::mean)
      centroids_t8 <- cbind(centroids_t8$V1, centroids_t8$V2)
      centroids_10 <- aggregate(d ~ k10, data = data, FUN = mean)
      centroids_10 <- cbind(centroids_10$V1, centroids_10$V2)
      centroids_t10 <- aggregate(d ~ kt10, data = data, FUN = base::mean)
      centroids_t10 <- cbind(centroids_t10$V1, centroids_t10$V2)
      
      # K_xy = length(unique(data$xymodes))
      K_xyt = length(unique(data$xytmodes))
      
      # Basic_MCMC(d, t_inf, basic_path)
      # M2_MCMC(d, K_xy, t_inf, centroids_xy, data$xymodes, M2_path_xy)
      M3_MCMC(d, K_xyt, t_inf, data$xytmodes, M2_path_xyt)
      M3_MCMC(d, 3, t_inf, data$k3, M3_path_k3)
      M3_MCMC(d, 3, t_inf, data$kt3, M3_path_kt3)
      M3_MCMC(d, 5, t_inf, data$k5, M3_path_k5)
      M3_MCMC(d, 5, t_inf, data$kt5, M3_path_kt5)
      M3_MCMC(d, 8, t_inf, data$k8, M3_path_k8)
      M3_MCMC(d, 8, t_inf, data$kt8, M3_path_kt8)
      M3_MCMC(d, 10, t_inf, data$k10, M3_path_k10)
      M3_MCMC(d, 10, t_inf, data$kt10, M3_path_kt10)
    }
  }
  
}
M3_CILM(5)
M4_CILM<-function(N){
  spatials <- c("homogeneous", "heterogeneous_small_K=3", "heterogeneous_small_K=5", "heterogeneous_small_K=8", "heterogeneous_large_K=3", "heterogeneous_large_K=5", "heterogeneous_large_K=8")
  for(s in 1:7){
    for(n in 1:N){
      data_path = paste("./simulation_data/",spatials[s],"/",spatials[s],"_",toString(n),".csv",sep = "")
      # M4 paths
      M4_path_xyt <- paste("./MCMC/M4/", spatials[s],"/",spatials[s],"_xyt_", n, ".csv", sep = "")
      M4_path_k3<-paste("./MCMC/M4/", spatials[s],"/",spatials[s],"_k",3,"_", n, ".csv", sep = "")
      M4_path_kt3<-paste("./MCMC/M4/", spatials[s],"/",spatials[s],"_kt",3,"_", n, ".csv", sep = "")
      M4_path_k5<-paste("./MCMC/M4/", spatials[s],"/",spatials[s],"_k",5,"_", n, ".csv", sep = "")
      M4_path_kt5<-paste("./MCMC/M4/", spatials[s],"/",spatials[s],"_kt",5,"_", n, ".csv", sep = "")
      M4_path_k8<-paste("./MCMC/M4/", spatials[s],"/",spatials[s],"_k",8,"_", n, ".csv", sep = "")
      M4_path_kt8<-paste("./MCMC/M4/", spatials[s],"/",spatials[s],"_kt",8,"_", n, ".csv", sep = "")
      M4_path_k10<-paste("./MCMC/M4/", spatials[s],"/",spatials[s],"_k",10,"_", n, ".csv", sep = "")
      M4_path_kt10<-paste("./MCMC/M4/",spatials[s],"/",spatials[s],"_kt",10,"_", n, ".csv", sep = "")
      
      data<-read.csv(data_path)
      d<-cbind(data$position_x, data$position_y)
      t_inf<-data$time_of_infected
      
      data$xytmodes<-as.numeric(factor(data$xytmodes))
      data$k3<-as.numeric(factor(data$k3))
      data$kt3<-as.numeric(factor(data$kt3))
      data$k5<-as.numeric(factor(data$k5))
      data$kt5<-as.numeric(factor(data$kt5))
      data$k8<-as.numeric(factor(data$k8))
      data$kt8<-as.numeric(factor(data$kt8))
      data$k10<-as.numeric(factor(data$k10))
      data$kt10<-as.numeric(factor(data$kt10))
      
      centroids_xyt <- aggregate(d ~ xytmodes, data = data, FUN = base::mean)
      centroids_xyt <- cbind(centroids_xyt$V1, centroids_xyt$V2)
      centroids_3 <- aggregate(d ~ k3, data = data, FUN = mean)
      centroids_3 <- cbind(centroids_3$V1, centroids_3$V2)
      centroids_t3 <- aggregate(d ~ kt3, data = data, FUN = base::mean)
      centroids_t3 <- cbind(centroids_t3$V1, centroids_t3$V2)
      centroids_5 <- aggregate(d ~ k5, data = data, FUN = mean)
      centroids_5 <- cbind(centroids_5$V1, centroids_5$V2)
      centroids_t5 <- aggregate(d ~ kt5, data = data, FUN = base::mean)
      centroids_t5 <- cbind(centroids_t5$V1, centroids_t5$V2)
      centroids_8 <- aggregate(d ~ k8, data = data, FUN = mean)
      centroids_8 <- cbind(centroids_8$V1, centroids_8$V2)
      centroids_t8 <- aggregate(d ~ kt8, data = data, FUN = base::mean)
      centroids_t8 <- cbind(centroids_t8$V1, centroids_t8$V2)
      centroids_10 <- aggregate(d ~ k10, data = data, FUN = mean)
      centroids_10 <- cbind(centroids_10$V1, centroids_10$V2)
      centroids_t10 <- aggregate(d ~ kt10, data = data, FUN = base::mean)
      centroids_t10 <- cbind(centroids_t10$V1, centroids_t10$V2)
      
      K_xyt = length(unique(data$xytmodes))
      
      M4_MCMC(d, K_xyt, t_inf, centroids_xyt, data$xytmodes, M4_path_xyt)
      M4_MCMC(d, 3, t_inf, centroids_3, data$k3, M4_path_k3)
      M4_MCMC(d, 3, t_inf, centroids_t3, data$kt3, M4_path_kt3)
      M4_MCMC(d, 5, t_inf, centroids_5, data$k5, M4_path_k5)
      M4_MCMC(d, 5, t_inf, centroids_t5, data$kt5, M4_path_kt5)
      M4_MCMC(d, 8, t_inf, centroids_8, data$k8, M4_path_k8)
      M4_MCMC(d, 8, t_inf, centroids_t8, data$kt8, M4_path_kt8)
      M4_MCMC(d, 10, t_inf, centroids_10, data$k10, M4_path_k10)
      M4_MCMC(d, 10, t_inf, centroids_t10, data$kt10, M4_path_kt10)
    }
  }
}
M4_CILM(5)
