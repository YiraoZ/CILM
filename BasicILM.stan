data {
  int<lower=1> N;//number of individuals
  int<lower=0> t_start;//start time
  int<lower=0> t_end;//end time
  int<lower=0> period;//infection period
  matrix[N,2] d;//position matrix
  vector[N] t_inf;//table of infected time
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
}
model{
  a0 ~ gamma(1.5, 1);
  beta ~ gamma(2, 3);
  for(t in t_start:t_end){
    vector[N] p = rep_vector(0, N);
    for(n in 1:N){
      if(t_inf[n]==-1 || t_inf[n]>=t){
        for(j in 1:N){
          if(j!= n && t_inf[j]!=-1 && t_inf[j]<t && t_inf[j]+period >=t){
            p[n] += pow(T_dis[n, j], -beta);
          }
        }
        if(t_inf[n]==t && p[n]!=0){
          target += log1m_exp(-a0 * p[n]);
        }
        if(t_inf[n]>t || t_inf[n]==-1){
          target += -a0 * p[n];
        }
      }
    }
  }
}



