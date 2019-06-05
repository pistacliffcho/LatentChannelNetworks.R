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
  
  vec<double> colMeans(){
    vec<double> ans(nCols);
    double curVal;
    for(int j = 0; j < nCols; j++){
      curVal = 0.0;
      for(int i = 0; i < nRows; i++){
        curVal += vals[ i + j * nRows ];
      }
      ans[j] = curVal / (double(nRows));
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



#endif
