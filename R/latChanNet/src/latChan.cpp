#include <Rcpp.h>
#include <Rmath.h>
#include <vector>
#include <iostream>
#include "utils.h"
#include <RcppParallel.h>
#define vec std::vector
using namespace Rcpp;

/***
 * Latent Channel Class
 ***/

class LCN{
public:
  int nNodes;
  int dim;
  double tol, pTol;
  vec<vec<int>>edgeList;
  Mat pmat;
  vec<double> pbar;
  
  vec<vec<double>> cache_probs;
  vec<vec<double*>> cache_map;
  
  double err;
  
  void ingestEdges(List lst);
  void initializeNode(int);
  void initializeCache();
  void parInitCache();
  LCN(List edgeList, NumericMatrix input_pmat);
  
  double edgeProb(int i, int j);
  double node_llk(int i);
  double llk();

  double update_pik(int i, int k);
  void one_ecm_update(int i, int k);
  void one_ecm_iter();
  void one_em_iter();
  void one_par_em_iter();
  List em(int max_its, int type, double tol, double pTol);
  
  NumericMatrix get_pmat();
  void set_pmat(NumericMatrix m);
  NumericVector computeTheta(int i, int j);
  NumericVector expectedConnections(int i);
};


/***
 * Initialization methods
 ***/

LCN::LCN(List input_edgeList, NumericMatrix input_pmat){
  pmat = Mat(input_pmat);
  pbar = pmat.colMeans();
  dim = input_pmat.cols();
  nNodes = input_pmat.rows();
  ingestEdges(input_edgeList);
  initializeCache();
  pTol = 0.00000001;
  tol = 0.0001;
}

void LCN::initializeNode(int i){
  int j, n_these_edges;
  n_these_edges = edgeList[i].size();
  for(int ii = 0; ii < n_these_edges; ii++){
    j = edgeList[i][ii];
    cache_probs[i][ii] = edgeProb(i,j);
  }
}

void LCN::initializeCache(){
  for(int i = 0; i < nNodes; i++){ initializeNode(i); }
}


int findTransposeInd(int i, int j, List rEdgeList){
  int r_i = i + 1;
  IntegerVector these_edges = rEdgeList[j];
  int n_edges = these_edges.length();
  for(int ii = 0; ii < n_edges; ii++){
    if(these_edges[ii] == r_i){
      return(ii);
    }
  }
  stop("Lookup unsuccessful!!");
  return(0);
}

void LCN::ingestEdges(List lst){
  nNodes = lst.length();
  edgeList.resize(nNodes);
  cache_probs.resize(nNodes);
  cache_map.resize(nNodes);
  int this_length, j, t_ind;
  for(int i = 0; i < nNodes; i++){
    IntegerVector theseEdges = lst[i];
    this_length = theseEdges.length();
    edgeList[i].resize(this_length);
    cache_probs[i].resize(this_length);
    cache_map[i].resize(this_length);
    for(int ii = 0; ii < this_length; ii++){
      // Switching from 1 to zero index
      j = theseEdges[ii] - 1;
      edgeList[i][ii] = j;
    }
  }
  for(int i = 0; i < nNodes; i++){
    IntegerVector theseEdges = lst[i];
    this_length = theseEdges.length();
    for(int ii = 0; ii < this_length; ii++){
      j = theseEdges[ii] - 1;
      t_ind = findTransposeInd(i,j,lst);
      cache_map[i][ii] = &(cache_probs[j][t_ind]);
    }
  }
}


/***
 * Likelihood functions
 ***/

double LCN::edgeProb(int i, int j){
  double pNoEdge = 1.0;
  double pik, pjk;
  for(int k = 0; k < dim; k++){
    pik = pmat(i,k);
    pjk = pmat(j,k);
    pNoEdge = pNoEdge * (1.0 - pik * pjk );
  }
  double ans = 1.0 - pNoEdge;
  return(ans);
}

double LCN::node_llk(int i){
  int nEdges = edgeList[i].size();
  vec<int> hasEdge(nNodes, 0);
  for(int ii = 0; ii < nEdges; ii++){
    hasEdge[ edgeList[i][ii] ] = 1;
  }
  double this_edgeProb;
  double ans = 0.0;
  for(int j = 0; j < nNodes; j++){
    if(i == j){ continue; }
    this_edgeProb = edgeProb(i,j);
    if(hasEdge[j] == 1){ ans += log(this_edgeProb); }
    else{ ans += log(1.0 - this_edgeProb); }
  }
  return(ans);
}

double LCN::llk(){
  double ans = 0.0;
  for(int i = 0; i < nNodes; i++){
    ans += node_llk(i);
  }
  return(ans);
}

/***
 * EM Algorithm: shared between Serial and Parallel methods
 ***/

double LCN::update_pik(int i, int k){
  double pik = pmat(i,k);
  if( (pik < pTol) || (1.0 - pik < pTol) ){
    return(pik);
  }
  int this_J_tot = edgeList[i].size(); 
  if(this_J_tot == 0){
    pmat(i,k) = 0.0;
    return(0.0);
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
    pjk = pmat(j,k);
    pikpjk = pik * pjk;
    noEdgeContribution += pikpjk;
    this_edgeP = epPtr[j_cnt];
    edgeContribution +=  (pikpjk + (pik - pikpjk) * 
      ( 1.0 - (1.0 - this_edgeP) / (1.0 - pikpjk))  ) / this_edgeP;
  }
  noEdgeContribution -= double(this_J_tot) * pik;
  
  double ans = (edgeContribution + noEdgeContribution) 
    / (double(nNodes) - 1.0);
  
  return(ans);
}

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

/***
 * EM Algorithm: Parallel methods
 ***/

// Structure for updating edge probabilities in parallel 
struct ParInitCache : public RcppParallel::Worker{
  LCN* mod_ptr;
  ParInitCache(LCN* mod){ mod_ptr = mod; }
  void operator()(size_t begin, size_t end){
    for(int i = begin; i < end; i++){
      mod_ptr->initializeNode(i);
    }
  }
};

void LCN::parInitCache(){
  ParInitCache cache_worker(this);
  RcppParallel::parallelFor(0, nNodes, cache_worker);
}

struct ParEMIter : public RcppParallel::Worker{
  LCN* mod;
  Mat* pmat_new;
  ParEMIter(LCN* mod_ptr, Mat* pmat_new_ptr){
    mod = mod_ptr;
    pmat_new = pmat_new_ptr;
  }
  void operator()(size_t begin, size_t end){
    int i,k;
    for(int flat_i = begin; flat_i < end; flat_i++){
      i = flat_i % mod->nNodes;
      k = int(flat_i / mod->nNodes);
      pmat_new->operator()(i,k) = mod->update_pik(i,k);
    }
  }
};

void LCN::one_par_em_iter(){
  Mat pmat_new = pmat.copy();
  int tot_pars = nNodes * dim;
  ParEMIter par_em_worker(this, &pmat_new);
  RcppParallel::parallelFor(0, tot_pars, par_em_worker);
  err = compute_err(pmat, pmat_new);
  pmat = pmat_new;
  pbar = pmat.colMeans();
  parInitCache();
}


/***
 * EM Algorithm: Serial methods
 ***/

void LCN::one_em_iter(){
  Mat pmat_new = pmat.copy();
  for(int i = 0; i < nNodes; i++){
    for(int k = 0; k < dim; k++){
      pmat_new(i,k) = update_pik(i,k);
    }
  }
  err = compute_err(pmat, pmat_new);
  pmat = pmat_new;
  pbar = pmat.colMeans();
  initializeCache();
}

/***
 * ECM Algorithm: Serial (only option) methods
 ***/

void LCN::one_ecm_iter(){
  for(int k = 0; k < dim; k++){
    for(int i = 0; i < nNodes; i++){
      one_ecm_update(i,k);
    }
  }
}



void LCN::one_ecm_update(int i, int k){
  double pikOld = pmat(i,k);
  double pikNew = update_pik(i,k);
  if(pikNew == pikOld){return;}
  pmat(i,k) = pikNew;
  double pikDiff = pikNew - pikOld;
  err = max(err, my_abs(pikDiff));
  pbar[k] += pikDiff/(double(nNodes));
  
  double oldEdgeP,newEdgeP, oldPNoEdge, newPNoEdge, pjk;
  int j;
  double* epPtr = &cache_probs[i][0];
  int this_J_tot = edgeList[i].size();
  for(int j_cnt = 0; j_cnt < this_J_tot; j_cnt++){
    j = edgeList[i][j_cnt]; 
    pjk = pmat(j,k);
    oldEdgeP = epPtr[j_cnt];
    oldPNoEdge = 1.0 - oldEdgeP;
    newPNoEdge = oldPNoEdge * (1.0 - pikNew * pjk) / (1.0 - pikOld*pjk);
    newEdgeP = 1.0 - newPNoEdge;
    cache_probs[i][j_cnt] = newEdgeP;
    (*cache_map[i][j_cnt]) = newEdgeP;
  }
}

// EM algorithm method. type allows to dictate whether to use 
// ECM (type = 1) 
// EM  (type = 2)
// Parallel EM (type = 3)
List LCN::em(int max_its, int type, double rtol, double rpTol){
  tol = rtol;
  pTol = rpTol;
  int iter = 0;
  err = tol + 1.0;
  while( (iter < max_its) & (err > tol) ){
    iter++;
    // Error is updated *inside* one_iter
    err = 0.0;
    if(type == 1){ one_ecm_iter(); }
    if(type == 2){ one_em_iter(); }
    if(type == 3){ one_par_em_iter(); }
  }
  if(iter == max_its){
    Rprintf("Warning: maximum iterations reached\n");
  }
  
  List ans = List::create(Named("err") = err, 
                          Named("its") = iter);
  return(ans);
}


/***
 * Misc. R interface functions
 ***/

// extract fitted p-matrix
NumericMatrix LCN::get_pmat(){
  NumericMatrix ans = pmat.getNumMat();
  return(ans);
}

// set p-matrix
void LCN::set_pmat(NumericMatrix m){
  int nRows = m.rows();
  int nCols = m.cols();
  if(nRows != nNodes){ stop("nRows != nNodes");}
  if(nCols != dim){ stop("nCols != dim");}
  pmat = Mat(m);
  initializeCache();
}

// Compute channel connection probability, conditional on observed edge
NumericVector LCN::computeTheta(int i, int j){
  if(i < 0 | i >= nNodes){ stop("i out of bounds");}
  if(j < 0 | j >= nNodes){ stop("j out of bounds");}
  NumericVector ans(dim);
  double this_edgeProb = edgeProb(i,j);
  for(int k = 0; k < dim; k++){
    ans[k] = pmat(i,k) * pmat(j,k) / this_edgeProb;
  }
  return(ans);
}

// Compute expected number of connections through each channel 
// for a given node, conditional on observed graph
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


// Exposing class to R
RCPP_MODULE(LCN){
  class_<LCN>("LCN")
  .constructor<List, NumericMatrix>("Args: EdgeList, Initial p-mat")
  .method("em", &LCN::em)
  .method("llk", &LCN::llk)
  .method("get_pmat", &LCN::get_pmat)
  .method("set_pmat", &LCN::set_pmat)
  .method("computeTheta", &LCN::computeTheta)
  .method("expectedConnections", &LCN::expectedConnections);
}