#include <Rcpp.h>
#include <unordered_set>
#include <map>
using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
double log_likelihood_M2(vector<int> inf_time, vector<double> x, vector<double> y, int t_start, int t_end, vector<int> assign, vector<double> parameters){
    double a0 = parameters[0];
    double beta = parameters[1];
    double tbeta = parameters[2];
    int N = inf_time.size();
    auto max_cluster_num = max_element(assign.begin(), assign.end());
    int K = *max_cluster_num + 1;//number of clusters
    //clustering
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
    //calculate log likelihood day by day
    double log_lkhd = 0;
    for(int t=t_start; t<=t_end; t++){
        //count the infectious individuals in each cluster
        vector<int> inf_num(K,0);
        for(int i=0; i<N; i++){
            if(inf_time[i]!=-1 && inf_time[i]<t && inf_time[i]+3>=t){
                inf_num[assign[i]]++;
            }
        }
        //loop through all individuals
        for(int i=0; i<N; i++){
            if(inf_time[i]==-1 || inf_time[i]>=t) {//make sure the individual is susceptible at time t-1
                double kernel = 0;
                double spark = 0;
                int k = assign[i];
                //calculate kernel
                for(int j: clusters[k]){
                    if(i!=j && inf_time[j]!=-1 && inf_time[j]<t && inf_time[j]+3>=t){
                        double d = sqrt(pow((x[i] - x[j]),2) + pow((y[i]-y[j]), 2));
                        kernel += pow(d, -beta);
                    }
                }
                //calculate spark
                for(int m=0; m<K; m++){
                    if(k!=m && inf_num[m]!=0){
                        double d = sqrt(pow((centroids[k][0]-centroids[m][0]),2)+pow((centroids[k][1]-centroids[m][1]),2));
                        spark += inf_num[m]*pow(d, -tbeta);
                    }
                }
                double prob = 1-exp(-a0*(kernel+spark));
                if(inf_time[i]==t){
                    log_lkhd+=log(prob);
                }
                else{
                    log_lkhd+=log(1-prob);
                }
            }
        }
    }
    return log_lkhd;
}
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
double log_likelihood_M3(vector<int> inf_time, vector<double> x, vector<double> y, int t_start, int t_end, vector<int> assign, vector<double> parameters){
    double a0 = parameters[0];
    double beta = parameters[1];
    double tbeta = parameters[2];
    int N = inf_time.size();
    auto max_cluster_num = max_element(assign.begin(), assign.end());
    int K = *max_cluster_num + 1;//number of clusters
    //clustering
    map<int, vector<int> > clusters;
    for (int m = 0; m < K; m++) {
        clusters[m] = vector<int>();
    }
    for(int i = 0; i<N; i++){
        clusters[assign[i]].push_back(i);
    }
    //calculate log likelihood day by day
    double log_lkhd = 0;
    for(int t=t_start; t<=t_end; t++){
        //count the infectious individuals in each cluster
        vector<int> inf_num(K,0);
        // for(int i=0; i<N; i++){
        //     if(inf_time[i]!=-1 && inf_time[i]<t && inf_time[i]+3>=t){
        //         inf_num[assign[i]]++;
        //     }
        // }
        //calculate centroids
        vector<vector<double> > centroids(K, vector<double>(2, 0.0));
        for(int m=0; m<K;m++){
            for(int i: clusters[m]){
                if(inf_time[i]!= -1 && inf_time[i]<t && inf_time[t]+3>=t){
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
        //loop through all individuals
        for(int i=0; i<N; i++){
            if(inf_time[i]==-1 || inf_time[i]>=t) {//make sure the individual is susceptible at time t-1
                double kernel = 0;
                double spark = 0;
                int k = assign[i];
                //calculate kernel
                for(int j: clusters[k]){
                    if(i!=j && inf_time[j]!=-1 && inf_time[j]<t && inf_time[j]+3>=t){
                        double d = sqrt(pow((x[i] - x[j]),2) + pow((y[i]-y[j]), 2));
                        kernel += pow(d, -beta);
                    }
                }
                //calculate spark
                for(int m=0; m<K; m++){
                    if(k!=m && inf_num[m]!=0){
                        double d = sqrt(pow((centroids[k][0]-centroids[m][0]),2)+pow((centroids[k][1]-centroids[m][1]),2));
                        spark += inf_num[m]*pow(d, -tbeta);
                    }
                }
                double prob = 1-exp(-a0*(kernel+spark));
                if (prob <= 0) prob = 1e-10;
                if(inf_time[i]==t){
                    log_lkhd+=log(prob);
                }
                else{
                    log_lkhd+=log(1-prob);
                }
            }
        }
    }
    return log_lkhd;
}
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
double log_likelihood_M4(vector<int> inf_time, vector<double> x, vector<double> y, int t_start, int t_end, vector<int> assign, vector<double> parameters){
    double a0 = parameters[0];
    double beta = parameters[1];
    double tbeta = parameters[2];
    double delta = parameters[3];
    int N = inf_time.size();
    auto max_cluster_num = max_element(assign.begin(), assign.end());
    int K = *max_cluster_num + 1;//number of clusters
    //clustering
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
    //calculate log likelihood day by day
    double log_lkhd = 0;
    for(int t=t_start; t<=t_end; t++){
        //count the infectious individuals in each cluster
        vector<int> inf_num(K,0);
        for(int i=0; i<N; i++){
            if(inf_time[i]!=-1 && inf_time[i]<t && inf_time[i]+3>=t){
                inf_num[assign[i]]++;
            }
        }
        //loop through all individuals
        for(int i=0; i<N; i++){
            if(inf_time[i]==-1 || inf_time[i]>=t) {//make sure the individual is susceptible at time t-1
                double kernel = 0;
                double spark = 0;
                int k = assign[i];
                //calculate kernel
                for(int j: clusters[k]){
                    if(i!=j && inf_time[j]!=-1 && inf_time[j]<t && inf_time[j]+3>=t){
                        double d = sqrt(pow((x[i] - x[j]),2) + pow((y[i]-y[j]), 2));
                        kernel += pow(d, -beta);
                    }
                }
                //calculate spark
                for(int m=0; m<K; m++){
                    if(k!=m && inf_num[m]!=0){
                        double d = sqrt(pow((centroids[k][0]-centroids[m][0]),2)+pow((centroids[k][1]-centroids[m][1]),2));
                        spark += pow(inf_num[m], delta)*pow(d, -tbeta);
                    }
                }
                double prob = 1-exp(-a0*(kernel+spark));
                if(inf_time[i]==t){
                    log_lkhd+=log(prob);
                }
                else{
                    log_lkhd+=log(1-prob);
                }
            }
        }
    }
    return log_lkhd;
}
