#ifndef LATCHAN
#define LATCHAN

#include <Rcpp.h>
#include <Rmath.h>
#include <vector>
#include <iostream>
#include <RcppParallel.h>
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
  double get(int i, int j){ return(vals[ i + j * nRows ]); }
  void set(double x, int i, int j){ vals[i + j * nRows] = x; }
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
        curVal += get(i,j);
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


class LCN{
public:
  int nNodes;
  int dim;
  vec<vec<int> >edgeList;
  Mat pmat;
  Mat new_pmat;
  vec<double> pbar;
  
  vec<vec<double> > cache_probs;

  double err;
  
  void ingestEdges(List lst);
  void initializeCache();
  void par_init_cache();
  double edgeProb(int i, int j);
  void one_update(int i, int k, double pTol);
  void one_iter(double pTol);
  void par_one_iter();
  double llk();
  List cache_em(int max_its, bool use_par, double tol, double pTol);
  LCN(List edgeList, NumericMatrix input_pmat);
  NumericMatrix get_pmat();
  void set_pmat(NumericMatrix m);
  NumericVector computeTheta(int i, int j);
  NumericVector expectedConnections(int i);
};

double my_abs(double);
double max(double, double);

double LCN::llk(){
  double ans = 0.0;
  vec<bool> hasEdge(nNodes, false);
  int nEdges;
  double this_edgeProb;
  for(int i = 0; i < nNodes; i++){
    nEdges = edgeList[i].size();
    for(int ii = 0; ii < nEdges; ii++){
      hasEdge[ edgeList[i][ii] ] = true;
    }
    
    for(int j = 0; j < nNodes; j++){
      if(i == j){ continue; }
      this_edgeProb = edgeProb(i,j);
      if(hasEdge[j] == true){ ans += log(this_edgeProb); }
      else{ ans += log(1.0 - this_edgeProb); }
    }
    
    for(int ii = 0; ii < nEdges; ii++){
      hasEdge[ edgeList[i][ii] ] = FALSE;
    }
    
  }
  return(ans);
}

NumericMatrix deepcopy(NumericMatrix);

LCN::LCN(List input_edgeList, NumericMatrix input_pmat){
  pmat = Mat(input_pmat);
  new_pmat = pmat.copy();
  pbar = pmat.colMeans();
  dim = input_pmat.cols();
  nNodes = input_pmat.rows();
  ingestEdges(input_edgeList);
  initializeCache();
}

void LCN::initializeCache(){
  int j, n_these_edges;
  for(int i = 0; i < nNodes; i++){
    n_these_edges = edgeList[i].size();
    for(int ii = 0; ii < n_these_edges; ii++){
      j = edgeList[i][ii];
      cache_probs[i][ii] = edgeProb(i,j);
    }
  }
}

double LCN::edgeProb(int i, int j){
  double pNoEdge = 1.0;
  double pik, pjk;
  for(int k = 0; k < dim; k++){
    pik = pmat.get(i,k);
    pjk = pmat.get(j,k);
    pNoEdge = pNoEdge * (1.0 - pik * pjk );
  }
  double ans = 1.0 - pNoEdge;
  return(ans);
}

void LCN::ingestEdges(List lst){
  nNodes = lst.length();
  edgeList.resize(nNodes);
  cache_probs.resize(nNodes);
  int this_length, j;
  for(int i = 0; i < nNodes; i++){
    IntegerVector theseEdges = lst[i];
    this_length = theseEdges.length();
    edgeList[i].resize(this_length);
    cache_probs[i].resize(this_length);
    for(int ii = 0; ii < this_length; ii++){
      // Switching from 1 to zero index
      j = theseEdges[ii] - 1;
      edgeList[i][ii] = j;
    }
  }
}


void LCN::one_iter(double pTol){
  for(int k = 0; k < dim; k++){
    for(int i = 0; i < nNodes; i++){
      one_update(i,k,pTol);
    }
    pbar = pmat.colMeans();
    double this_err;
    for(int i = 0; i < pmat.vals.size(); i++){
      this_err = my_abs(pmat.vals[i] - new_pmat.vals[i]);
      err = max(err, this_err);
    }
    pmat = new_pmat.copy();
    initializeCache();
  }
}

void LCN::one_update(int i, int k, double pTol){
  double pik = pmat.get(i,k);
  if( (pik < pTol) || (1.0 - pik < pTol) ){
    return;
  }
  int this_J_tot = edgeList[i].size(); 
  if(this_J_tot == 0){
    pmat.set(0.0, i,k);
    return;
  }
  double edgeContribution = 0.0;
  double noEdgeContribution = nNodes * pik * 
    (1.0 - pbar[k]) - pik * (1.0 - pik);
  int j;
  double pjk, this_edgeP, pikpjk;
  int* jPtr = &edgeList[i][0];
  double* epPtr = &cache_probs[i][0];
  for(int j_cnt = 0; j_cnt < this_J_tot; j_cnt++){
    j = jPtr[j_cnt];
    pjk = pmat.get(j,k);
    pikpjk = pik * pjk;
    noEdgeContribution += pikpjk;
    this_edgeP = epPtr[j_cnt];
    edgeContribution +=  (pikpjk + (pik - pikpjk) * 
      ( 1.0 - (1.0 - this_edgeP) / (1.0 - pikpjk))  )/this_edgeP;
  }
  noEdgeContribution -= double(this_J_tot) * pik;
  
  double pikNew = (edgeContribution + noEdgeContribution) / (double(nNodes) - 1.0);
  new_pmat.set(pikNew, i,k);
}


double my_abs(double x){
  if( x > 0.0){ return(x); }
  else{ return(-x); }
}

double max(double a, double b){
  if(a > b){ return(a);}
  return(b);
}

List LCN::cache_em(int max_its, bool use_par, double tol, double pTol){
  int iter = 0;
  err = tol + 1.0;
  while( (iter < max_its) & (err > tol) ){
    iter++;
    // Error is updated *inside* one_iter
    err = 0.0;
    if(use_par){ par_one_iter(); }
    else{ one_iter(pTol); }
  }
  if(iter == max_its){
    Rprintf("Warning: maximum iterations reached\n");
  }
  
  List ans = List::create(Named("err") = err, 
                          Named("its") = iter);
  return(ans);
}

NumericMatrix LCN::get_pmat(){
  NumericMatrix ans = pmat.getNumMat();
  return(ans);
}

void LCN::set_pmat(NumericMatrix m){
  pmat = Mat(m);
  initializeCache();
}


NumericVector LCN::computeTheta(int i, int j){
  if(i < 0 | i >= nNodes){ stop("i out of bounds");}
  if(j < 0 | j >= nNodes){ stop("j out of bounds");}
  NumericVector ans(dim);
  double this_edgeProb = edgeProb(i,j);
  for(int k = 0; k < dim; k++){
    ans[k] = pmat.get(i,k) * pmat.get(j,k) / this_edgeProb;
  }
  return(ans);
}

NumericVector LCN::expectedConnections(int i){
  NumericVector ans(dim);
  int this_j;
  int nNodes = edgeList[i].size();
  for(int ii = 0; ii < nNodes; ii++){
    this_j = edgeList[i][ii];
    NumericVector this_theta = computeTheta(i, this_j);
    for(int k = 0; k < dim; k++){
      ans[k] = ans[k] + this_theta[k];
    }
  }
  return(ans);
}

struct ParInitCache : public RcppParallel::Worker{
  LCN* this_lcn;
  ParInitCache(LCN* mod){ this_lcn = mod; }
  void operator()(std::size_t begin, std::size_t end){
    int start = begin; int stop = end;
    int n_these_edges, j;
    for(int i = start; i < stop; i++){
      n_these_edges = this_lcn->edgeList[i].size();
      for(int ii = 0; ii < n_these_edges; ii++){
        j = this_lcn->edgeList[i][ii];
        this_lcn->cache_probs[i][ii] = this_lcn->edgeProb(i,j);
      }
    }
  }
};

void LCN::par_init_cache(){
  ParInitCache parCache(this);
  parallelFor(0, nNodes, parCache);
}

struct ParIter : public RcppParallel::Worker{
  LCN* this_lcn;
  ParIter(LCN* mod){ this_lcn = mod; }
  void operator()(std::size_t begin, std::size_t end){
    int nRow = this_lcn->pmat.nRows;
    int flat_index = begin;
    int i, k;
    while(flat_index < end){
      i = flat_index % nRow;
      k = floor(flat_index / nRow);
      this_lcn->one_update(i,k, 0.0000001);
      flat_index++;
    }
  }
};

void LCN::par_one_iter(){
  ParIter parIt(this);
  parallelFor(0, pmat.vals.size(), parIt);
  pbar = pmat.colMeans();
  double this_err;
  for(int i = 0; i < pmat.vals.size(); i++){
    this_err = my_abs(pmat.vals[i] - new_pmat.vals[i]);
    err = max(err, this_err);
  }
  pmat = new_pmat.copy();
  par_init_cache();
}


RCPP_MODULE(LCN){
  class_<LCN>("LCN")
  .constructor<List, NumericMatrix>("Args: EdgeList, Initial p-mat")
  .method("cache_em", &LCN::cache_em)
  .method("llk", &LCN::llk)
  .method("get_pmat", &LCN::get_pmat)
  .method("set_pmat", &LCN::set_pmat)
  .method("computeTheta", &LCN::computeTheta)
  .method("expectedConnections", &LCN::expectedConnections);
}

#endif
