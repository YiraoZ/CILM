#include <Rcpp.h>
#include <random>
#include <unordered_set>
#include <map>
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
vector<int> PPD_StandardILM(vector<int> inf_time, vector<double> x, vector<double> y, int t_start, int t_end, vector<double> parameters){
  double a0 = parameters[0];
  double beta = parameters[1];
  int N = inf_time.size();//population size
  //copy inf time before t_start
  vector<int> ppd(N,-1);
  for(int i=0; i<N; i++){
    if(inf_time[i]<t_start){
      ppd[i] = inf_time[i];
    }
  }
  //simulate ppd
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<double> status_dis(0.0, 1.0);
  for(int t=t_start; t<t_end; t++){
    for(int i=0; i<N; i++){
      if(ppd[i]==-1){
        double kernel = 0;
        //calculate kernel
        for(int j = 0; j<N; j++){
          if(i!=j && ppd[j]!=-1 && ppd[j]<t && ppd[j]+3>=t){
            double d = sqrt(pow((x[i] - x[j]),2) + pow((y[i]-y[j]), 2));
            kernel += pow(d, -beta);
          }
        }
        double prob = 1-exp(-a0*kernel);
        double e = status_dis(gen);
        if(e < prob){
          ppd[i] = t;
        }
      }
    }
  }
  return ppd;
}
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
vector<int> PPD_M2(vector<int> inf_time, vector<double> x, vector<double> y, int t_start, int t_end, vector<int> assign, vector<double> parameters) {
    double a0 = parameters[0];
    double beta = parameters[1];
    double tbeta = parameters[2];
    int N = inf_time.size();//population size
    auto max_cluster_num = max_element(assign.begin(), assign.end());
    int K = *max_cluster_num + 1;//number of clusters
    //Clustering
    map<int, vector<int> > clusters;
    for (int m = 0; m < K; m++) {
      clusters[m] = vector<int>();
    }
    for(int i = 0; i<N; i++){
      clusters[assign[i]].push_back(i);
    }
    //calculate centroids
    vector<vector<double> > centroids(K, vector<double>(2, 0.0));
    for(int m=0; m<K;m++){
      for(int i: clusters[m]){
        centroids[m][0]+=x[i];
        centroids[m][1]+=y[i];
      }
      if(clusters[m].size()!=0){
        centroids[m][0]/=clusters[m].size();
        centroids[m][1]/=clusters[m].size();
      }
    }
    //copy inf time before t_start
     vector<int> ppd(N,-1);
    for(int i=0; i<N; i++){
      if(inf_time[i]<t_start){
        ppd[i] = inf_time[i];
      }
    }
    //simulate ppd
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> status_dis(0.0, 1.0);
    for(int t=t_start; t<t_end; t++){
      //calculate the inf number in each cluster
      vector<int> inf_num(K,0);
      for(int i=0; i<N; i++){
        if(ppd[i]!=-1 && ppd[i]<t && ppd[i]+3>=t){
          inf_num[assign[i]]++;
        }
      }
      for(int i=0; i<N; i++){
        if(ppd[i]==-1){
          double kernel = 0;
          double spark = 0;
          int k = assign[i];
          //calculate kernel
          for(int j: clusters[k]){
            if(i!=j && ppd[j]!=-1 && ppd[j]<t && ppd[j]+3>=t){
              double d = sqrt(pow((x[i] - x[j]),2) + pow((y[i]-y[j]), 2));
              kernel += pow(d, -beta);
            }
          }
          //calculate spark function
          for(int m=0; m<K; m++){
            if(k!=m){
              double d = sqrt(pow((centroids[k][0]-centroids[m][0]),2)+pow((centroids[k][1]-centroids[m][1]),2));
              spark += inf_num[m]*pow(d, -tbeta);
            }
          }
          double prob = 1-exp(-a0*(kernel+spark));
          double e = status_dis(gen);
          if(e < prob){
            ppd[i] = t;
          }
        }
      }
    }
     return ppd;
}
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
vector<int> PPD_M3(vector<int> inf_time, vector<double> x, vector<double> y, int t_start, int t_end, vector<int> assign, vector<double> parameters) {
  double a0 = parameters[0];
  double beta = parameters[1];
  double tbeta = parameters[2];
  int N = inf_time.size();//population size
  auto max_cluster_num = max_element(assign.begin(), assign.end());
  int K = *max_cluster_num + 1;//number of clusters
  //Clustering
  map<int, vector<int> > clusters;
  for (int m = 0; m < K; m++) {
    clusters[m] = vector<int>();
  }
  for(int i = 0; i<N; i++){
    clusters[assign[i]].push_back(i);
  }
  //copy inf time before t_start
  vector<int> ppd(N,-1);
  for(int i=0; i<N; i++){
    if(inf_time[i]<t_start){
      ppd[i] = inf_time[i];
    }
  }
  //simulate ppd
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<double> status_dis(0.0, 1.0);
  for(int t=t_start; t<t_end; t++){
    //calculate the inf number in each cluster
    vector<int> inf_num(K,0);
    // for(int i=0; i<N; i++){
    //   if(ppd[i]!=-1 && ppd[i]<t && ppd[i]+3>=t){
    //     inf_num[assign[i]]++;
    //   }
    // }
    //calculate centroids
    vector<vector<double> > centroids(K, vector<double>(2, 0.0));
    for(int m=0; m<K;m++){
      for(int i: clusters[m]){
        if(ppd[i]!= -1 && ppd[i]<t && ppd[t]+3>=t){
          inf_num[m]++;
          centroids[m][0]+=x[i];
          centroids[m][1]+=y[i];
        }
      }
      if(inf_num[m]!=0){
        centroids[m][0]/=inf_num[m];
        centroids[m][1]/=inf_num[m];
      }
    }
    for(int i=0; i<N; i++){
      if(ppd[i]==-1){
        double kernel = 0;
        double spark = 0;
        int k = assign[i];
        //calculate kernel
        for(int j: clusters[k]){
          if(i!=j && ppd[j]!=-1 && ppd[j]<t && ppd[j]+3>=t){
            double d = sqrt(pow((x[i] - x[j]),2) + pow((y[i]-y[j]), 2));
            kernel += pow(d, -beta);
          }
        }
        //calculate spark function
        for(int m=0; m<K; m++){
          if(k!=m && inf_num[m]!=0){
            double d = sqrt(pow((centroids[k][0]-centroids[m][0]),2)+pow((centroids[k][1]-centroids[m][1]),2));
            spark += inf_num[m]*pow(d, -tbeta);
          }
        }
        double prob = 1-exp(-a0*(kernel+spark));
        double e = status_dis(gen);
        if(e < prob){
          ppd[i] = t;
        }
      }
    }
  }
  return ppd;
}
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
vector<int> PPD_M4(vector<int> inf_time, vector<double> x, vector<double> y, int t_start, int t_end, vector<int> assign, vector<double> parameters) {
  double a0 = parameters[0];
  double beta = parameters[1];
  double tbeta = parameters[2];
  double delta = parameters[3];
  int N = inf_time.size();//population size
  auto max_cluster_num = max_element(assign.begin(), assign.end());
  int K = *max_cluster_num + 1;//number of clusters
  //Clustering
  map<int, vector<int> > clusters;
  for (int m = 0; m < K; m++) {
    clusters[m] = vector<int>();
  }
  for(int i = 0; i<N; i++){
    clusters[assign[i]].push_back(i);
  }
  //calculate centroids
  vector<vector<double> > centroids(K, vector<double>(2, 0.0));
  for(int m=0; m<K;m++){
    for(int i: clusters[m]){
      centroids[m][0]+=x[i];
      centroids[m][1]+=y[i];
    }
    if(clusters[m].size()!=0){
      centroids[m][0]/=clusters[m].size();
      centroids[m][1]/=clusters[m].size();
    }
  }
  //copy inf time before t_start
  vector<int> ppd(N,-1);
  for(int i=0; i<N; i++){
    if(inf_time[i]<t_start){
      ppd[i] = inf_time[i];
    }
  }
  //simulate ppd
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<double> status_dis(0.0, 1.0);
  for(int t=t_start; t<t_end; t++){
    //calculate the inf number in each cluster
    vector<int> inf_num(K,0);
    for(int i=0; i<N; i++){
      if(ppd[i]!=-1 && ppd[i]<t && ppd[i]+3>=t){
        inf_num[assign[i]]++;
      }
    }
    for(int i=0; i<N; i++){
      if(ppd[i]==-1){
        double kernel = 0;
        double spark = 0;
        int k = assign[i];
        //calculate kernel
        for(int j: clusters[k]){
          if(i!=j && ppd[j]!=-1 && ppd[j]<t && ppd[j]+3>=t){
            double d = sqrt(pow((x[i] - x[j]),2) + pow((y[i]-y[j]), 2));
            kernel += pow(d, -beta);
          }
        }
        //calculate spark function
        for(int m=0; m<K; m++){
          if(k!=m && inf_num[m]!=0){
            double d = sqrt(pow((centroids[k][0]-centroids[m][0]),2)+pow((centroids[k][1]-centroids[m][1]),2));
            spark += pow(inf_num[m], delta)*pow(d, -tbeta);
          }
        }
        double prob = 1-exp(-a0*(kernel+spark));
        double e = status_dis(gen);
        if(e < prob){
          ppd[i] = t;
        }
      }
    }
  }
  return ppd;
}
