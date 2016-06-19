#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
double cal_auc1(NumericVector ps_hat, NumericVector yr, NumericVector x, 
                NumericVector R){
  int n = ps_hat.size();
  double auc_nu = 0;
  double auc_de = 0;
  double tt[n] , I12[n] ;
    for (int i=0; i<n; i++) {
      for(int j=0; j<n; j++){
        tt[j] = R[i]/ps_hat[i] * yr[i] * (1 - yr[j]) * R[j]/ps_hat[j];
        if(x[i] == x[j]) I12[j] = 0.5;
        if(x[i] > x[j]) I12[j] = 1;
        if(x[i] < x[j]) I12[j] = 0;
        auc_nu += (tt[j] * I12[j]);
        auc_de += (tt[j]);
      }
    }
    double auc_A1 = auc_nu/auc_de;
    return auc_A1;
}


// [[Rcpp::export]]
double cal_auc2(NumericVector ps_hat, NumericVector yr, NumericVector x, 
                NumericVector R){
  int n = ps_hat.size();
  double auc_nu = 0;
  double auc_de = 0;
  double tt[n] , I12[n] ;
    for (int i=0; i<n; i++) {
      for(int j=0; j<n; j++){
        tt[j] = R[i]/ps_hat[i] * yr[i] * (1 - yr[j] * R[j]/ps_hat[j])  ;
        if(x[i] == x[j]) I12[j] = 0.5;
        if(x[i] > x[j]) I12[j] = 1;
        if(x[i] < x[j]) I12[j] = 0;
        auc_nu += (tt[j] * I12[j]);
        auc_de += (tt[j]);
      }
    }
    double auc_A2 = auc_nu/auc_de;
    return auc_A2;
}


// [[Rcpp::export]]
NumericVector cal_F1(NumericVector ps_hat, double auc, NumericVector yr, NumericVector x, 
                NumericVector R){
  int n = ps_hat.size();
  NumericVector F1(n);
  double tt[n] , I12[n];
    for (int i=0; i<n; i++) {
      F1[i] = 0.0;
      for(int j=0; j<n; j++){
        tt[j] = R[i]/ps_hat[i] * yr[i] * (1 - yr[j]) * R[j]/ps_hat[j] ;
        if(x[i] == x[j]) I12[j] = 0.5;
        if(x[i] > x[j]) I12[j] = 1;
        if(x[i] < x[j]) I12[j] = 0;
        F1[i] += tt[j] * (I12[j] - auc);
        
        tt[j] = R[j]/ps_hat[j] * yr[j] * (1 - yr[i]) * R[i]/ps_hat[i] ;
        if(x[i] == x[j]) I12[j] = 0.5;
        if(x[j] > x[i]) I12[j] = 1;
        if(x[j] < x[i]) I12[j] = 0;
        F1[i] += tt[j] * (I12[j] - auc);
      }
      F1[i] = F1[i]/n;
    }
    return (F1);
}


// [[Rcpp::export]]
NumericVector cal_F2(NumericVector ps_hat, double auc, NumericVector yr, NumericVector x, 
                NumericVector R){
  int n = ps_hat.size();
  double tt[n] , I12[n];
  NumericVector F2(n);
    for (int i=0; i<n; i++) {
      F2[i] = 0.0;
      for(int j=0; j<n; j++){
        tt[j] = R[i]/ps_hat[i] * yr[i] * (1 - yr[j] * R[j]/ps_hat[j]) ;
        if(x[i] == x[j]) I12[j] = 0.5;
        if(x[i] > x[j]) I12[j] = 1;
        if(x[i] < x[j]) I12[j] = 0;
        if(j == i) tt[j] = 0;
        F2[i] += tt[j] * (I12[j] - auc);
        tt[j] = R[j]/ps_hat[j] * yr[j] * (1 - yr[i] * R[i]/ps_hat[i])  ;
        if(x[i] == x[j]) I12[j] = 0.5;
        if(x[j] > x[i]) I12[j] = 1;
        if(x[j] < x[i]) I12[j] = 0;
        if(j == i) tt[j] = 0;
        F2[i] += tt[j] * (I12[j] - auc);
      }
       F2[i] = F2[i]/n;
    }
    return F2;
}

// [[Rcpp::export]]
NumericVector cal_gamma1(NumericVector ps_hat, double auc, NumericVector yr, NumericVector x, 
                NumericVector R, NumericMatrix vsub){
  int n = ps_hat.size();
  int k = vsub.ncol();
  NumericVector gamma(k+3);
  double tt1[n], tt2[n], I12[n];
  
  for (int s=0; s<k+3; s++) gamma[s] = 0;
  
    for (int i=0; i<n; i++) {
      for(int j=0; j<n; j++){
        tt1[j] = R[i]*(1/ps_hat[i] - 1) * yr[i] * (1 - yr[j]) * R[j]/ps_hat[j] ;
        tt2[j] = R[i]/ps_hat[i] * yr[i] * (1 - yr[j]) * R[j] * (1/ps_hat[j] - 1);
        if(x[i] == x[j]) I12[j] = 0.5;
        if(x[i] > x[j]) I12[j] = 1;
        if(x[i] < x[j]) I12[j] = 0;
        I12[j] = I12[j] - auc;
        gamma[0] += I12[j] * (tt1[j] + tt2[j]);
        gamma[1] += I12[j] * (tt1[j] * x[i] + tt2[j] * x[j]);
        for(int s =0; s<k; s++) gamma[2+s] += I12[j] * (tt1[j] * vsub(i, s) + tt2[j] * vsub(j, s));
        gamma[k+2] += I12[j] * (tt1[j] * yr[i] + tt2[j] * yr[j]);
      }
    }
    
    //for (int s=0; s<k+3; s++) gamma[s] = gamma[s]/(n*n);
    return (gamma);
}


// [[Rcpp::export]]
NumericVector cal_gamma2(NumericVector ps_hat, double auc, NumericVector yr, NumericVector x, 
                NumericVector R, NumericMatrix vsub){
  int n = ps_hat.size();
  int k = vsub.ncol();
  NumericVector gamma(k+3);
  double tt1[n], tt2[n], I12[n];
  
  for (int s=0; s<k+3; s++) gamma[s] = 0;
  
    for (int i=0; i<n; i++) {
      for(int j=0; j<n; j++){
        tt1[j] = R[i]*(1/ps_hat[i] - 1) * yr[i] * (1 - yr[j] * R[j]/ps_hat[j]) ;
        tt2[j] = R[i]/ps_hat[i] * yr[i] * (1 - yr[j] * R[j] * (1/ps_hat[j]) - 1);
        if(x[i] == x[j]) I12[j] = 0.5;
        if(x[i] > x[j]) I12[j] = 1;
        if(x[i] < x[j]) I12[j] = 0;
        I12[j] = I12[j] - auc;
        if(j == i) I12[j] = 0;
        gamma[0] += I12[j] * (tt1[j] + tt2[j]);
        gamma[1] += I12[j] * (tt1[j] * x[i] + tt2[j] * x[j]);
        for(int s =0; s<k; s++) gamma[2+s] += I12[j] * (tt1[j] * vsub(i, s) + tt2[j] * vsub(j, s));
        gamma[k+2] += I12[j] * (tt1[j] * yr[i] + tt2[j] * yr[j]);
      }
    }
    
    //for (int s=0; s<k+3; s++) gamma[s] = gamma[s]/(n*n);
    return (gamma);
}


// [[Rcpp::export]]
NumericVector cal_gamma1_v2(NumericVector ps_hat, double auc, NumericVector yr, 
                NumericVector R, NumericMatrix cova){
  int n = ps_hat.size();
  int k = cova.ncol();
  NumericVector gamma(k+2);
  double tt1[n], tt2[n], I12[n], x[n];
  
  for (int s=0; s<k+2; s++) gamma[s] = 0;
  for (int i=0; i<n; i++) x[i] = cova(i, 0);
  
    for (int i=0; i<n; i++) {
      for(int j=0; j<n; j++){
        tt1[j] = R[i]*(1/ps_hat[i] - 1) * yr[i] * (1 - yr[j]) * R[j]/ps_hat[j] ;
        tt2[j] = R[i]/ps_hat[i] * yr[i] * (1 - yr[j]) * R[j] * (1/ps_hat[j] - 1);
        if(x[i] == x[j]) I12[j] = 0.5;
        if(x[i] > x[j]) I12[j] = 1;
        if(x[i] < x[j]) I12[j] = 0;
        I12[j] = I12[j] - auc;
        gamma[0] += I12[j] * (tt1[j] + tt2[j]);
        for(int s = 0; s<k; s++) gamma[1+s] += I12[j] * (tt1[j] * cova(i, s) + tt2[j] * cova(j, s));
        gamma[k+1] += I12[j] * (tt1[j] * yr[i] + tt2[j] * yr[j]);
      }
    }
    
    //for (int s=0; s<k+3; s++) gamma[s] = gamma[s]/(n*n);
    return (gamma);
}


// [[Rcpp::export]]
NumericVector cal_gamma2_v2(NumericVector ps_hat, double auc, NumericVector yr, 
                NumericVector R, NumericMatrix cova){
  int n = ps_hat.size();
  int k = cova.ncol();
  NumericVector gamma(k+2);
  double tt1[n], tt2[n], I12[n], x[n];
  
  for (int s=0; s<k+2; s++) gamma[s] = 0;
  for (int i=0; i<n; i++) x[i] = cova(i, 0);
  
    for (int i=0; i<n; i++) {
      for(int j=0; j<n; j++){
        tt1[j] = R[i]*(1/ps_hat[i] - 1) * yr[i] * (1 - yr[j] * R[j]/ps_hat[j]) ;
        tt2[j] = R[i]/ps_hat[i] * yr[i] * (1 - yr[j] * R[j] * (1/ps_hat[j]) - 1);
        if(x[i] == x[j]) I12[j] = 0.5;
        if(x[i] > x[j]) I12[j] = 1;
        if(x[i] < x[j]) I12[j] = 0;
        I12[j] = I12[j] - auc;
        if(j == i) I12[j] = 0;
        gamma[0] += I12[j] * (tt1[j] + tt2[j]);
        for(int s =0; s<k; s++) gamma[1+s] += I12[j] * (tt1[j] * cova(i, s) + tt2[j] * cova(j, s));
        gamma[k+1] += I12[j] * (tt1[j] * yr[i] + tt2[j] * yr[j]);
      }
    }
    
    //for (int s=0; s<k+3; s++) gamma[s] = gamma[s]/(n*n);
    return (gamma);
}
