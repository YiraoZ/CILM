data {
  int<lower=1> N;//number of individuals
  int<lower=1> K;//number of clusters
  int<lower=0> t_start;//start time
  int<lower=0> t_end;//end time
  int<lower=0> period;//infection period
  matrix[N, 2] d;//position matrix
  matrix[K, 2] centroids;//centroids position matrix
  vector[N] t_inf;//table of infected time
  int<lower=0> cluster_id[N];//cluster assignments
}
transformed data{
  matrix[N,N] T_dis;
  for(i in 1:N){
    for(j in i:N){
      if(i==j){
        T_dis[i, j] = 0;
      }else{
      T_dis[i, j] = distance(d[i],d[j]);
      T_dis[j, i] = T_dis[i, j];
      }
    }
  }
  matrix[K, K] C_dis;//distances between centroids
  for(i in 1:K){
    for(j in i:K){
      if(i==j){
        C_dis[i, j] = 0;
      }else{
      C_dis[i, j] = distance(centroids[i],centroids[j]);
      C_dis[j, i] = C_dis[i, j];
      }
    }
  }
}
parameters{
  real<lower=0> a0;
  real<lower=0> beta;
  real<lower=0> tbeta;
}
model{
  a0 ~ gamma(1.5, 1);
  beta ~ gamma(2, 3);
  tbeta ~ gamma(2, 3);
  for(t in t_start:t_end){
    vector[N] p = rep_vector(0, N);
    vector[K] num_inf = rep_vector(0, K);//number of infected in each cluster
    for(n in 1:N){//count the number of infected in each cluster
      if(t_inf[n]!=-1 && t_inf[n]<t && t_inf[n]+period >=t){
        num_inf[cluster_id[n]] += 1;
      }
    }
    for(n in 1:N){
      if(t_inf[n]==-1 || t_inf[n]>=t){
        for(j in 1:N){
          if(j!=n && t_inf[j]<t && t_inf[j]+period >=t && cluster_id[j]==cluster_id[n]){
            p[n] += pow(T_dis[n, j], -beta);
          }
        }
        for(k in 1:K){
          if(k!=cluster_id[n]){
            p[n] += num_inf[k] * pow(C_dis[cluster_id[n], k], -tbeta);
          }
        }
        if(t_inf[n]==t && p[n]!=0){
          target += log1m_exp(-(a0) * p[n]);
        }
        if(t_inf[n]>t || t_inf[n]==-1){
          target += -(a0) * p[n];
      }
    }
  }
}
}
