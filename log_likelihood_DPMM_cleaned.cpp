#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
double NB2_pmf(int y, double mu, double phi) {
  // calculate combination C(y + phi - 1, y)
  double combination = tgamma(y+phi)/(tgamma(y+1)*tgamma(phi));
  // calculate pmf
  double pmf_value = combination * pow((mu/(mu+phi)), y) * pow((phi/(mu+phi)), phi);
  return pmf_value;
}
// [[Rcpp::export]]
double hNB2_pmf(int y, double mu, double phi, double theta) {
  if(y == 0){
    return theta;
  }
  else{
    // calculate pmf at y = 0
    double pmf_0 = NB2_pmf(0, mu, phi);
    double pmf_value = (1-theta)*NB2_pmf(y, mu, phi)/(1-pmf_0);
    return pmf_value;
  }
}
// [[Rcpp::export]]
double log_likelihood_full(vector<double> x, vector<double> y, vector<int> t, vector<double> c, vector<double> v, vector<double> prob){
  //single phi
  //theta length is 3: w_x, w_y, phi
  //centroids length is 4: mu_x, mu_y, mu_t, theta(hurdle)
  double log_lkhd = 0;
  double w_x = v[0];
  double w_y = v[1];
  double phi = v[2];
  int N = t.size();
  int M = prob.size();
  double sqrt_det_sigma = w_x*w_y;
  for(int i = 0; i<N; i++){
    double lkhd_cluster = 0;
    for(int j = 0; j<M; j++){
      vector<double> diff(2);
      diff[0] = x[i] - c[4*j];
      diff[1] = y[i] - c[4*j+1];
      double mu = c[4*j+2];
      double hurdle = c[4*j+3];
      lkhd_cluster += prob[j]*hNB2_pmf(t[i],mu,phi,hurdle)*exp(-0.5*(pow(diff[0]/w_x,2)+pow(diff[1]/w_y,2)))/(sqrt_det_sigma*sqrt(pow(2*M_PI,2)));
    }
    if(lkhd_cluster == 0){
      return -INFINITY;
    }
    log_lkhd += log(lkhd_cluster);
  }
  return log_lkhd;
}
// [[Rcpp::export]]
double log_likelihood_conditional(vector<double> x, vector<double> y, vector<int> t, vector<double> c, vector<double> v, vector<int> g){
  //with the cluster known
  //single phi
  //theta length is 3: w_x, w_y, phi
  //centroids length is 4: mu_x, mu_y, mu_t, theta(hurdle)
  double log_lkhd = 0;
  double w_x = v[0];
  double w_y = v[1];
  double phi = v[2];
  int N = t.size();
  double sqrt_det_sigma = w_x*w_y;
  for(int i=0; i<N; i++){
    int j = g[i];//j is the cluster that i individual belongs to 
    vector<double> diff(2);
    diff[0] = x[i] - c[4*(j-1)];
    diff[1] = y[i] - c[4*(j-1)+1];
    double mu = c[4*(j-1)+2];
    double hurdle = c[4*(j-1)+3];
    log_lkhd += log(hNB2_pmf(t[i],mu,phi,hurdle))-0.5*(pow(diff[0]/w_x,2)+pow(diff[1]/w_y,2))-log(sqrt_det_sigma*sqrt(pow(2*M_PI,2)));
  }
  return log_lkhd;
}
// [[Rcpp::export]]
double log_likelihood_conditional_t(vector<int> t, vector<double> c, vector<double> v, vector<int> g, int j){
  //with the cluster known
  //single phi
  //centroids length is 4: mu_x, mu_y, mu_t, theta(hurdle)
  double log_lkhd = 0;
  double phi = v[2];
  int N = t.size();
  for(int i=0; i<N; i++){
    if(j==g[i]){
      double mu = c[4*(j-1)+2];
      double hurdle = c[4*(j-1)+3];
      log_lkhd += log(hNB2_pmf(t[i],mu,phi,hurdle));
    }
  }
  return log_lkhd;
}
// [[Rcpp::export]]
double log_likelihood_xy(vector<double> x, vector<double> y, vector<double> c, vector<double> v, vector<double> prob){
  double log_lkhd = 0;
  double w_x = v[0];
  double w_y = v[1];
  int N = x.size();
  int M = prob.size();
  double sqrt_det_sigma = w_x*w_y;
  //double sqrt_det_sigma = w_x*w_y;
  for(int i = 0; i<N; i++){
    double lkhd_cluster = 0;
    for(int j = 0; j<M; j++){
      vector<double> diff(2);
      diff[0] = x[i] - c[2*j];
      diff[1] = y[i] - c[2*j+1];
      lkhd_cluster += prob[j]*exp(-0.5*(pow(diff[0]/w_x,2)+pow(diff[1]/w_y,2)))/(sqrt_det_sigma*sqrt(pow(2*M_PI,2)));
    }
    log_lkhd += log(lkhd_cluster);
  }
  return log_lkhd;
}


