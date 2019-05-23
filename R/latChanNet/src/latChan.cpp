#include <Rcpp.h>
#include <Rmath.h>
#include <vector>
#include <iostream>
#define vec std::vector
using namespace Rcpp;


class LCN{
public:
  int nNodes;
  int dim;
  vec<vec<int>>edgeList;
  NumericMatrix pmat;
  NumericVector pbar;
  
  vec<vec<double>> cache_probs;
  vec<vec<double*>> cache_map;
  
  double err;
  
  void ingestEdges(List lst);
  void initializeCache();
  double edgeProb(int i, int j);
  void one_update(int i, int k, double pTol);
  void one_iter(double pTol);
  double llk();
  List cache_em(int max_its, double tol, double pTol);
  LCN(List edgeList, NumericMatrix input_pmat);
  NumericMatrix get_pmat();
};

double my_abs(double);
double max(double, double);

double LCN::llk(){
  double ans = 0.0;
  LogicalVector hasEdge(nNodes, FALSE);
  int nEdges;
  double this_edgeProb;
  for(int i = 0; i < nNodes; i++){
    nEdges = edgeList[i].size();
    for(int ii = 0; ii < nEdges; ii++){
      hasEdge[ edgeList[i][ii] ] = TRUE;
    }
    
    for(int j = 0; j < nNodes; j++){
      if(i == j){ continue; }
      this_edgeProb = edgeProb(i,j);
      if(hasEdge[j] == TRUE){ ans += log(this_edgeProb); }
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
  pmat = deepcopy(input_pmat);
  pbar = colMeans(pmat);
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
    pik = pmat(i,k);
    pjk = pmat(j,k);
    pNoEdge = pNoEdge * (1.0 - pik * pjk );
  }
  double ans = 1.0 - pNoEdge;
  return(ans);
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

void LCN::one_iter(double pTol){
  for(int k = 0; k < dim; k++){
    for(int i = 0; i < nNodes; i++){
      one_update(i,k,pTol);
    }
  }
}

void LCN::one_update(int i, int k, double pTol){
  double pik = pmat(i,k);
  if( (pik < pTol) || (1.0 - pik < pTol) ){
    return;
  }
  int this_J_tot = edgeList[i].size(); 
  if(this_J_tot == 0){
    pmat(i,k) = 0.0;
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
    pjk = pmat(j,k);
    pikpjk = pik * pjk;
    noEdgeContribution += pikpjk;
    this_edgeP = epPtr[j_cnt];
    edgeContribution +=  (pikpjk + (pik - pikpjk) * 
      ( 1.0 - (1.0 - this_edgeP) / (1.0 - pikpjk))  )/this_edgeP;
  }
  noEdgeContribution -= double(this_J_tot) * pik;
  
  double pikNew = (edgeContribution + noEdgeContribution) 
    / (double(nNodes) - 1.0);
  pmat(i,k) = pikNew;
  double pikDiff = pikNew - pik;
  err = max(err, my_abs(pikDiff));
  pbar[k] += pikDiff/(double(nNodes));
  
  double oldEdgeP,newEdgeP, oldPNoEdge, newPNoEdge;
  for(int j_cnt = 0; j_cnt < this_J_tot; j_cnt++){
    j = edgeList[i][j_cnt]; 
    pjk = pmat(j,k);
    oldEdgeP = epPtr[j_cnt];
    oldPNoEdge = 1.0 - oldEdgeP;
    newPNoEdge = oldPNoEdge * (1.0 - pikNew * pjk) / (1.0 - pik*pjk);
    newEdgeP = 1.0 - newPNoEdge;
    cache_probs[i][j_cnt] = newEdgeP;
    (*cache_map[i][j_cnt]) = newEdgeP;
  }
}

NumericMatrix deepcopy(NumericMatrix m){
  int nRows = m.rows();
  int nCols = m.cols();
  NumericMatrix ans(nRows, nCols);
  for(int j = 0; j < nCols; j++){
    for(int i = 0; i < nRows; i++){
      ans(i,j) = m(i,j);
    }
  }
  return(ans);
}

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

List LCN::cache_em(int max_its, double tol, double pTol){
  int iter = 0;
  err = tol + 1.0;
  while( (iter < max_its) & (err > tol) ){
    iter++;
    // Error is updated *inside* one_iter
    err = 0.0;
    one_iter(pTol);
  }
  if(iter == max_its){
    Rprintf("Warning: maximum iterations reached\n");
  }
  
  List ans = List::create(Named("err") = err, 
                          Named("its") = iter);
  return(ans);
}

NumericMatrix LCN::get_pmat(){
  NumericMatrix ans = deepcopy(pmat);
  return(ans);
}


RCPP_MODULE(LCN){
  class_<LCN>("LCN")
  .constructor<List, NumericMatrix>("Args: EdgeList, Initial p-mat")
  .method("cache_em", &LCN::cache_em)
  .method("llk", &LCN::llk)
  .method("get_pmat", &LCN::get_pmat);
}