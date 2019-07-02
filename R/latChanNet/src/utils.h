#ifndef UTILS
#define UTILS

#include <Rcpp.h>
#include <Rmath.h>
#include <vector>
#include <iostream>
#define vec std::vector
using namespace Rcpp;


/*
 * Threadsafe Matrix
 */

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
    }
  }
  return(err);
}


#endif
