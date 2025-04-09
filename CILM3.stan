data {
  int<lower=1> N;//number of individuals
  int<lower=1> K;//number of clusters
  int<lower=0> t_start;//start time
  int<lower=0> t_end;//end time
  int<lower=0> period;//infection period
  matrix[N, 2] d;//position matrix
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
    vector[N] kernel = rep_vector(0, N);
    vector[N] spark = rep_vector(0, N);
    vector[K] num_inf = rep_vector(0, K);//number of infected in each cluster
    matrix[K, 2] centroids = rep_matrix(0, K, 2);//centroids position matrix
    for(n in 1:N){//count the number of infected in each cluster
      if(t_inf[n]!=-1 && t_inf[n]<t && t_inf[n]+period >=t){
        num_inf[cluster_id[n]] += 1;
        centroids[cluster_id[n], 1]+=d[n, 1];
        centroids[cluster_id[n], 2]+=d[n, 2];
      }
    }
    for(k in 1:K){//divide by the number of each cluster
      if(num_inf[k]!=0){
        centroids[k, 1]/=num_inf[k];
        centroids[k, 2]/=num_inf[k];
      }
    }
    for(n in 1:N){
      if(t_inf[n]==-1 || t_inf[n]>=t){
        for(j in 1:N){
          if(j!=n && t_inf[j]<t && t_inf[j]+period >=t && cluster_id[j]==cluster_id[n]){
            kernel[n] += pow(T_dis[n, j], -beta);
          }
        }
        for(k in 1:K){
          if(k!=cluster_id[n] && num_inf[k]!=0){
            spark[n] += num_inf[k] * pow(distance(centroids[cluster_id[n]], centroids[k]), -tbeta);
          }
        }
        if(t_inf[n]==t && kernel[n]+spark[n]!=0){
          target += log1m_exp(-(a0) * (kernel[n]+spark[n]));
        }
        if(t_inf[n]>t || t_inf[n]==-1){
          target += -(a0) * (kernel[n]+spark[n]);
      }
    }
  }
}
}






