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
setwd("/Users/yiraozhang/Desktop/UC_PhD/DPMM/cleaned")
#cpp files
sourceCpp('./log_likelihood_DPMM_cleaned.cpp')
#calculate probability based on U
cal_prob <- function(U){
  U<-as.numeric(U)
  N<-length(U)
  prob <- U*c(1, cumprod(1-U)[1:(N-1)])
  return(prob)
}
#DPMM_MCMC: DPMM apply to (x, y, t)
DPMM_MCMC<-function(x, y, t, M, chain_len){#M is the max cluster number, chain_len is the max iteration number
  N<-length(x)
  x<-min(t)+(max(t) - min(t)) * (x - min(x)) /(max(x) - min(x))
  y<-min(t)+(max(t) - min(t)) * (y - min(y)) /(max(y) - min(y))
  #prior parameters
  a<-1
  b<-0.25
  #proposal parameters
  sd_x<-2
  sd_y<-2
  sd_phi<-2
  #create the output lists
  v<-data.frame(matrix(nrow = chain_len, ncol = 3))#omega_x, omega_y and phi
  c<-data.frame(matrix(nrow = chain_len, ncol = 4*M))#c^x, c^y, c^t and theta
  g<-data.frame(matrix(nrow = chain_len, ncol = N))#cluster membership of N individuals
  U<-data.frame(matrix(nrow = chain_len, ncol = M))#sticks
  r<-data.frame(matrix(nrow = chain_len, ncol = 1))
  #Initialization 
  v[1,]<-c(1,1,1)
  c[1,]<-0.5
  g[1,]<-sample(1:M,N,replace=TRUE)
  r[1,]<-rgamma(1,a,b)
  U[1,]<-rbeta(M,1,r[1,])
  prob<-cal_prob(U[1,])
  #MCMC
  pb <- txtProgressBar(min = 1, max = chain_len, style = 3)
  for(iter in 2:chain_len){
    #update variance-related parameters
    v[iter,]<-v[iter-1,]
    ##update w_x
    log_lkhd_crt<-log_likelihood_full(x, y, t, as.numeric(c[iter-1,]), as.numeric(v[iter,]), as.numeric(prob))
    v[iter,1]<-rtruncnorm(1, a=0, b=Inf, mean=v[iter-1,1], sd=sd_x)
    log_lkhd_nxt<-log_likelihood_full(x, y, t, as.numeric(c[iter-1,]), as.numeric(v[iter,]), as.numeric(prob))
    log_ratio<-(log_lkhd_nxt+log(dgamma(v[iter,1],1.5,1))+log(dtruncnorm(v[iter-1,1],a=0, b=Inf,mean=v[iter,1],sd=sd_x)))-(log_lkhd_crt+log(dgamma(v[iter-1,1],1.5,1))+log(dtruncnorm(v[iter,1],a=0, b=Inf,mean=v[iter-1,1],sd=sd_x)))
    if(is.na(log_ratio)||log(runif(1))>log_ratio){
      v[iter,1]<-v[iter-1,1]
    }
    else{
      log_lkhd_crt<-log_lkhd_nxt
    }
    ##update w_y
    v[iter,2]<-rtruncnorm(1, a=0, b=Inf, mean=v[iter-1,2], sd=sd_y)
    log_lkhd_nxt<-log_likelihood_full(x, y, t, as.numeric(c[iter-1,]), as.numeric(v[iter,]), as.numeric(prob))
    log_ratio<-(log_lkhd_nxt+log(dgamma(v[iter,2],1.5,1))+log(dtruncnorm(v[iter-1,2],a=0, b=Inf,mean=v[iter,2],sd=sd_y)))-(log_lkhd_crt+log(dgamma(v[iter-1,2],1.5,1))+log(dtruncnorm(v[iter,2],a=0, b=Inf,mean=v[iter-1,2],sd=sd_y)))
    if(is.na(log_ratio)||log(runif(1))>log_ratio){
      v[iter,2]<-v[iter-1,2]
    }
    else{
      log_lkhd_crt<-log_lkhd_nxt
    }
    ##update phi
    v[iter,3]<-rtruncnorm(1, a=0, b=Inf, mean=v[iter-1,3], sd=sd_phi)
    log_lkhd_nxt<-log_likelihood_full(x, y, t, as.numeric(c[iter-1,]), as.numeric(v[iter,]), as.numeric(prob))
    log_ratio<-(log_lkhd_nxt+log(dgamma(v[iter,3],1.5,1))+log(dtruncnorm(v[iter-1,3],a=0, b=Inf,mean=v[iter,3],sd=sd_phi)))-(log_lkhd_crt+log(dgamma(v[iter-1,3],1.5,1))+log(dtruncnorm(v[iter,3],a=0, b=Inf,mean=v[iter-1,3],sd=sd_phi)))
    if(is.na(log_ratio)||log(runif(1))>log_ratio){
      v[iter,3]<-v[iter-1,3]
    }
    #update assignments
    for(i in 1:N){
      prob_i<-log(prob)
      for(j in 1:M){
        prob_i[j]<-prob_i[j]+log(hNB2_pmf(t[i],c[iter-1, 4*(j-1)+3],v[iter,3],c[iter-1, 4*(j-1)+4]))-(((x[i]-c[iter-1,4*(j-1)+1])/v[iter,1])^2)/2-(((y[i]-c[iter-1,4*(j-1)+2])/v[iter,2])^2)/2
      }
      prob_i <- exp(prob_i-max(prob_i))# no need to normalize
      g[iter,i] <- sample(1:M, 1, prob =prob_i)
    }
    #update U
    U[iter,]<-rep(1, M)
    for(k in 1:(M-1)){
      U[iter,k]<-rbeta(1, sum(g[iter,]==k)+1, sum(g[iter,]>k)+r[iter-1,])
    }
    prob<-cal_prob(U[iter,])
    #update hyoerparameter gamma
    U[iter,][U[iter,]>0.9999]<-0.9999
    r[iter,]<-rgamma(1, sum(U[iter,]>0)+a-1, b-sum(log(1-U[iter,-M])))
    #update cluster-based c and theta
    c[iter,]<-c[iter-1,]
    for(j in 1:M){
      #update means
      #make sure the number in jth cluster is greater than 0
      n_m<-sum(g[iter,] == j)
      if(n_m>0){
        ##update c^x
        x_m_avg <- mean(x[g[iter,] == j])
        c[iter,(j-1)*4+1]<-rnorm(1, mean = x_m_avg, sd = v[iter,1]/sqrt(n_m))
        ##update c^y
        y_m_avg <- mean(y[g[iter, ] == j])
        c[iter, (j-1)*4+2]<-rnorm(1, mean = y_m_avg, sd = v[iter,2]/sqrt(n_m))
        ##update c^t
        log_lkhd_t_m_crt<-log_likelihood_conditional_t(t, as.numeric(c[iter, ]), as.numeric(v[iter,]), as.numeric(g[iter,]),j)
        c[iter,(j-1)*4+3]<-runif(1, min(t), max(t))
        log_lkhd_t_m_nxt<-log_likelihood_conditional_t(t, as.numeric(c[iter, ]), as.numeric(v[iter,]), as.numeric(g[iter,]),j)
        log_ratio<-log_lkhd_t_m_nxt-log_lkhd_t_m_crt
        if(is.na(log_ratio)||log(runif(1))>log_ratio){
          c[iter,(j-1)*4+3]<-c[iter-1,(j-1)*4+3]
        }
        else{
          log_lkhd_t_m_crt<-log_lkhd_t_m_nxt
        }
        #update theta
        c[iter,(j-1)*4+4]<-rbeta(1, 2, 2)
        log_lkhd_t_m_nxt<-log_likelihood_conditional_t(t, as.numeric(c[iter, ]), as.numeric(v[iter,]), as.numeric(g[iter,]),j)
        log_ratio<-log_lkhd_t_m_nxt-log_lkhd_t_m_crt
        if(is.na(log_ratio)||log(runif(1))>log_ratio){
          c[iter,(j-1)*4+4]<-c[iter-1,(j-1)*4+4]
        }
      }
    }
    setTxtProgressBar(pb, iter)
  }
  res <- list("g"=g,"v"=v, "c"=c,"U"=U, "r" = r)
  close(pb)
  return(res)
}
#DPMM_MCMC_spatial: DPMM apply to spatial data (x, y)
DPMM_MCMC_spatial<-function(x, y, M, chain_len){#M is the max cluster number, chain_len is the max iteration number
  N<-length(x)
  #prior parameters
  a<-1
  b<-0.25
  #proposal parameters
  sd_x<-1
  sd_y<-1
  #create the output samples
  v<-data.frame(matrix(nrow = chain_len, ncol = 2))#omega_x, omega_y
  c<-data.frame(matrix(nrow = chain_len, ncol = 2*M))#c^x, c^y
  g<-data.frame(matrix(nrow = chain_len, ncol = N))#cluster membership of N individuals
  U<-data.frame(matrix(nrow = chain_len, ncol = M))#sticks
  r<-data.frame(matrix(nrow = chain_len, ncol = 1))
  #Initialization 
  v[1,]<-c(1,1)
  c[1,]<-0.5
  g[1,]<-sample(1:M,N,replace=TRUE)
  r[1,]<-rgamma(1,a,b)
  U[1,]<-rbeta(M,1,r[1,])
  prob<-cal_prob(U[1,])
  #MCMC
  pb <- txtProgressBar(min = 1, max = chain_len, style = 3)
  for(iter in 2:chain_len){
    #update variance-related parameters 
    v[iter,]<-v[iter-1,]
    ##update w_x
    log_lkhd_crt<-log_likelihood_xy(x, y, as.numeric(c[iter-1,]), as.numeric(v[iter,]), as.numeric(prob))
    v[iter,1]<-rtruncnorm(1, a=0, b=Inf, mean=v[iter-1,1], sd=sd_x)
    log_lkhd_nxt<-log_likelihood_xy(x, y, as.numeric(c[iter-1,]), as.numeric(v[iter,]), as.numeric(prob))
    log_ratio<-(log_lkhd_nxt+log(dgamma(v[iter,1],1.5,1))+log(dtruncnorm(v[iter-1,1],a=0, b=Inf,mean=v[iter,1],sd=sd_x)))-(log_lkhd_crt+log(dgamma(v[iter-1,1],1.5,1))+log(dtruncnorm(v[iter,1],a=0, b=Inf,mean=v[iter-1,1],sd=sd_x)))
    if(is.na(log_ratio)||log(runif(1))>log_ratio){
      v[iter,1]<-v[iter-1,1]
    }
    else{
      log_lkhd_crt<-log_lkhd_nxt
    }
    ##update w_y
    v[iter,2]<-rtruncnorm(1, a=0, b=Inf, mean=v[iter-1,2], sd=sd_y)
    log_lkhd_nxt<-log_likelihood_xy(x, y, as.numeric(c[iter-1,]), as.numeric(v[iter,]), as.numeric(prob))
    log_ratio<-(log_lkhd_nxt+log(dgamma(v[iter,2],1.5,1))+log(dtruncnorm(v[iter-1,2],a=0, b=Inf,mean=v[iter,2],sd=sd_y)))-(log_lkhd_crt+log(dgamma(v[iter-1,2],1.5,1))+log(dtruncnorm(v[iter,2],a=0, b=Inf,mean=v[iter-1,2],sd=sd_y)))
    if(is.na(log_ratio)||log(runif(1))>log_ratio){
      v[iter,2]<-v[iter-1,2]
    }
    #update assignments
    for(i in 1:N){
      prob_i<-log(prob)
      for(j in 1:M){
        prob_i[j]<-prob_i[j]-(((x[i]-c[iter-1,2*(j-1)+1])/v[iter,1])^2)/2-(((y[i]-c[iter-1,2*(j-1)+2])/v[iter,2])^2)/2
      }
      prob_i <- exp(prob_i-max(prob_i))# no need to normalize
      g[iter,i] <- sample(1:M, 1, prob =prob_i)
    }
    #update U
    U[iter,]<-rep(1, M)
    for(k in 1:(M-1)){
      U[iter,k]<-rbeta(1, sum(g[iter,]==k)+1, sum(g[iter,]>k)+r[iter-1,])
    }
    prob<-cal_prob(U[iter,])
    #update r
    U[iter,][U[iter,]>0.9999]<-0.9999
    r[iter,]<-rgamma(1, sum(U[iter,]>0)+a-1, b-sum(log(1-U[iter,-M])))
    #update c
    c[iter,]<-c[iter-1,]
    for(j in 1:M){
      #update means
      #make sure the number in jth cluster is greater than 0
      n_m<-sum(g[iter,] == j)
      if(n_m>0){
        ##update c^x
        x_m_avg <- mean(x[g[iter,] == j])
        c[iter,(j-1)*2+1]<-rnorm(1, mean = x_m_avg, sd = v[iter,1]/sqrt(n_m))
        ##update c^y
        y_m_avg <- mean(y[g[iter, ] == j])
        c[iter, (j-1)*2+2]<-rnorm(1, mean = y_m_avg, sd = v[iter,2]/sqrt(n_m))
      }
    }
    setTxtProgressBar(pb, iter)
  }
  res <- list("g"=g,"v"=v, "c"=c,"U"=U, "r" = r)
  close(pb)
  return(res)
}
#Clustering using DPMM for (x,y,t), DPMM for (x,y)
#extract output statistics: mode of assignments
DPMM_mode<-function(g){
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  g_modes <- sapply(g, Mode)
  return(g_modes)
}
DPMM_xyt<-function(chain_len, N){#N is the # of datasets under each scenario
  datapath_prefix<-"./simulation_data/"
  spatials = c("homogeneous", "heterogeneous_small_K=3", "heterogeneous_small_K=5","heterogeneous_small_K=8","heterogeneous_large_K=3","heterogeneous_large_K=5", "heterogeneous_large_K=8")
  for(k in 1:length(spatials)){
    for(n in 1:N){
      datapath <- paste(datapath_prefix,spatials[k],"/",spatials[k],"_",toString(n),".csv",sep = "")
      data<-read.csv(datapath)
      t<-data$time_of_infected+1
      M<-30
      res<-DPMM_MCMC(data$position_x, data$position_y, t, M, chain_len)
      data$xytmode<-DPMM_mode(res$g[(chain/2):chain])
      write.csv(data, datapath, row.names = FALSE)
    }
  }
}
DPMM_xyt(2000,10)
DPMM_xy<-function(chain_len, N){#N is the # of datasets under each scenario
  datapath_prefix<-"./simulation_data/"
  spatials = c("homogeneous", "heterogeneous_small_K=3", "heterogeneous_small_K=5","heterogeneous_small_K=8","heterogeneous_large_K=3","heterogeneous_large_K=5", "heterogeneous_large_K=8")
  for(k in 1:length(spatials)){
    for(n in 1:N){
      datapath <- paste(datapath_prefix,spatials[k],"/",spatials[k],"_",toString(n),".csv",sep = "")
      data<-read.csv(datapath)
      M<-30
      res<-DPMM_MCMC_spatial(data$position_x, data$position_y, M, chain_len)
      data$xymode<-DPMM_mode(res$g[(chain/2):chain])
      write.csv(data, datapath, row.names = FALSE)
    }
  }
}
DPMM_xy(2000,10)
#Clustering using K-means when K = 3, 5, 8, 10
datapath_prefix<-"./simulation_data/"
spatials = c("homogeneous", "heterogeneous_small_K=3", "heterogeneous_small_K=5","heterogeneous_small_K=8","heterogeneous_large_K=3","heterogeneous_large_K=5", "heterogeneous_large_K=8")
for(s in 1:7){
  for(n in 1:10){
    datapath<-paste0(datapath_prefix, spatials[s],"/", spatials[s], "_", n, ".csv")
    data<-read.csv(datapath)
    xyt<-data.frame(data$position_x, data$position_x, data$time_of_infected)
    xy<-data.frame(data$position_x, data$position_y)#kmeans only consider xy
    #scaled_data <- scale(xyt)
    for(k in c(3, 5, 8, 10)){
      k_res <- kmeans(xy, centers = k, iter.max = 100)
      kt_res <- kmeans(xyt, centers = k, iter.max = 100)
      data[[paste0("k",k)]]<-k_res$cluster
      data[[paste0("kt",k)]]<-kt_res$cluster
      write.csv(data, datapath, row.names = FALSE)
    }
  }
}



