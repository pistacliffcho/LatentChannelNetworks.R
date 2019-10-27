#ifndef UTILS
#define UTILS

#include <Rcpp.h>
#include <Rmath.h>
#include <vector>
#include <iostream>
#define vec std::vector
using namespace Rcpp;


/**
 * Threadsafe Matrix
 **/
class Mat{
public:
  vec<double> vals;
  int nRows;
  int nCols;
  double& operator()(int i, int j){ return(vals[ i + j * nRows ]); }
  Mat(){
    nRows = 0;
    nCols = 0;
  }
  Mat(NumericMatrix mat){
    nRows = mat.rows();
    nCols = mat.cols();
    vals.resize(nRows * nCols);
    for(int i = 0; i < (nRows * nCols); i++){
      vals[i] = mat[i];
    }
  }
  NumericMatrix getNumMat(){
    NumericMatrix ans(nRows, nCols);
    for(int i = 0; i < (nRows * nCols); i++){
      ans[i] = vals[i];
    }
    return(ans);
  }
  
  vec<double> colSums(){
    vec<double> ans(nCols);
    double curVal;
    int this_col_flat;
    for(int j = 0; j < nCols; j++){
      curVal = 0.0;
      this_col_flat = j * nRows;
      for(int i = 0; i < nRows; i++){
        curVal += vals[ i + this_col_flat ];
      }
      ans[j] = curVal;
    }
    return(ans);
  }
  
  vec<double> colMeans(){
    vec<double> ans = colSums();
    for(int j = 0; j < nCols; j++){
      ans[j] = ans[j] / (double(nRows));
    }
    return(ans);
  }
  
  Mat copy(){
    Mat ans;
    ans.nRows = nRows;
    ans.nCols = nCols;
    ans.vals.resize(nCols * nRows);
    for(int i = 0; i < nRows * nCols; i++){
      ans.vals[i] = vals[i];
    }
    return(ans);
  }
};


/**
 * Basic functions
 */

double my_abs(double x){
  if( x > 0.0){ return(x); }
  else{ return(-x); }
}

double max(double a, double b){
  if(a > b){
    return(a);
  }
  return(b);
}

/*
 * Computation functions
 */

double compute_err(Mat &m1, Mat &m2){
  double err = 0; 
  double this_err;
  for(int i = 0; i < m1.nRows; i++){
    for(int j = 0; j < m1.nCols; j++){
      this_err = my_abs(m1(i,j) - m2(i,j));
      err = max(err, this_err);
//      err += this_err;
    }
  }
  return(err);
}

void checkInd(int i, int max){
  if(i < 0) stop("negative index");
  if(i >= max) stop("i >= max");
}
  

  
void fillPosVector(int i, vec<int> &pos_vec, Mat &theta){
  int posCount = 0;
  int J = theta.nCols;
  for(int j = 0; j < J; j++){
    if(theta(i,j) > 0){ posCount++; }
  }
  pos_vec.resize(posCount);
  int pos_ind = -1;
  for(int j = 0; j < J; j++){
    if(theta(i,j) > 0){
      pos_ind++;
      pos_vec[pos_ind] = j;
    }
  }
}  

void setPosInds(vec<vec<int> > &pos_vecs, Mat &theta){
  int nRows = theta.nRows;
  pos_vecs.resize(nRows);
  for(int i = 0; i < nRows; i++){
    fillPosVector(i, pos_vecs[i], theta);
  }
}

double min(double a, double b){
  if(a < b) return(a);
  return(b);
}


/***
 * 
 * Tools for Over Relataxtion of EM algorithm,
 * not currently used!!
 * 
 ***/

// Monotonic Overrelaxation EM augmenteation step of Yu 2012
// Specifically for probabilities
void emRelaxedProb(Mat &theta_new, Mat &theta_old, double w){
  int nRows = theta_old.nRows;
  int nCols = theta_old.nCols;
  double new_val;
  for(int i = 0; i < nRows; i++){
    for(int j = 0; j < nCols; j++){
      new_val = (1.0 + w) * theta_new(i,j) - w * theta_old(i,j);
      new_val = max(0.0, new_val);
      new_val = min(1.0, new_val);
      theta_new(i,j) = new_val;
    }
  }
}

// Relaxation augmentation for non-negative parameters 
void emRelaxedPos(Mat &theta_new, Mat &theta_old, double w){
  int nRows = theta_old.nRows;
  int nCols = theta_old.nCols;
  double new_val;
  for(int i = 0; i < nRows; i++){
    for(int j = 0; j < nCols; j++){
      new_val = (1.0 + w) * theta_new(i,j) - w * theta_old(i,j);
      new_val = max(0.0, new_val);
      theta_new(i,j) = new_val;
    }
  }
}

NumericVector predict_lat_edges(IntegerVector r_i, IntegerVector r_j, 
                      NumericMatrix pmat, std::string model){
  int mod;
  if(model == "LCN"){ mod = 1;}
  else if(model == "BKN"){ mod = 2;}
  else{ stop("model not recognized"); }
  
  int K = pmat.ncol();
  int nRows = pmat.nrow();
  
  int n = r_i.size();
  int nj = r_j.size();
  int ci, cj;
  if(n != nj) stop("length(i) != length(j)");
  NumericVector ans(n);
  double this_par = 0.0;
  for(int i = 0; i < n; i++){
    if(mod == 1){ this_par = 1.0; }
    if(mod == 2){ this_par = 0.0; }
    ci = r_i[i] - 1; cj = r_j[i] - 1;
    if(ci < 0 | ci >= nRows | cj < 0 | cj >= nRows){ stop("invalid node id"); }
    for(int k = 0; k < K; k++){
      if(mod == 1){ this_par *= (1.0 - pmat(ci, k) * pmat(cj, k)); }
      if(mod == 2){ this_par += pmat(ci, k) * pmat(cj, k); }
    }
    if(mod == 1){ ans[i] = 1.0 - this_par; }
    if(mod == 2){ ans[i] = this_par; }
  }
  return(ans);
}

NumericMatrix predict_crossedge(IntegerVector r_i, IntegerVector r_j, 
                                NumericMatrix pmat, std::string model){
  int n_i = r_i.size();
  int n_j = r_j.size();
  
  NumericMatrix ans(n_i, n_j);
  IntegerVector rep_i(n_j);
  for(int i = 0; i < n_i; i++){
    rep_i = rep_i * 0;
    rep_i = rep_i + r_i[i];
    NumericVector this_row = predict_lat_edges(rep_i, r_j, pmat, model);
    ans.row(i) = this_row;
  }
  return(ans);
}


NumericMatrix chanConnect(IntegerVector r_i, IntegerVector r_j, 
                          NumericMatrix pmat, std::string model){
  if(model == "BKN"){stop("chanConnect not implemented for BKN yet");}
  NumericVector exp_connects = predict_lat_edges(r_i, r_j, pmat, model);
  int nChans = pmat.ncol();
  int nRows = exp_connects.size();
  NumericMatrix ans(nRows, nChans);
  double this_prob;
  int ci,cj;
  for(int i = 0; i < nRows; i++){
    this_prob = exp_connects[i];
    ci = r_i[i] - 1;
    cj = r_j[i] - 1;
    if(this_prob == 0){
      for(int j = 0; j < nChans; j++){ ans(i,j) = 0.0; }
    }
    else{
      for(int j = 0; j < nChans; j++){
        ans(i,j) = pmat(ci, j) * pmat(cj, j) / this_prob;
      }
    }
  }
  return(ans);
}

#endif
