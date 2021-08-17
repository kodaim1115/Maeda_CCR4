// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(mvtnorm)]]
#include <RcppArmadillo.h>
//#include <RcppArmadilloExtensions/sample.h>
//#include <progress.hpp>
#include <mvtnormAPI.h>
#include<Rmath.h>
using namespace Rcpp;

// [[Rcpp::export]]
double calc_SSE(arma::mat Y, arma::mat W, arma::mat mu, int L){

  double SSE=0;
  int nl;
  
  for(int l=0;l<L;l++){
    nl = arma::sum(W.col(l));
    arma::uvec clustid(nl); clustid.fill(0);
    clustid = arma::find(W.col(l)==1);
    for(int i=0;i<nl;i++){
      SSE += arma::sum(square(Y.row(clustid(i))-mu.col(l).t()));
    }
  } 
  return SSE;
}

// [[Rcpp::export]]
double calc_SSE_v2(arma::mat Y, arma::mat W,int L){
  
  double SSE=0;
  int nl;
  int K = Y.n_cols;
  arma::rowvec mu(K);
  
  for(int l=0;l<L;l++){
    nl = arma::sum(W.col(l));
    arma::uvec clustid(nl); clustid.fill(0);
    clustid = arma::find(W.col(l)==1);
    mu.fill(0);
    for(int i=0;i<nl;i++){
      mu += Y.row(clustid(i));
    }
    mu = mu/nl;
    for(int i=0;i<nl;i++){
      SSE += arma::sum(square(Y.row(clustid(i))-mu));
    }
  } 
  return SSE;
}