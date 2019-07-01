/***
 * 
 * Implementation of Ball, Karrer and Newman (2011)
 * model, a related model to the latent channel networks
 */

#include <Rcpp.h>
#include <Rmath.h>
#include <math.h>
#include <vector>
#include <iostream>
#include "utils.h"
#include <RcppParallel.h>
#define vec std::vector
using namespace Rcpp;


class BKN{
public:
  int nNodes;
  int dim;
  double tol, pTol;
  vec<vec<vec<int> > >edgeCounts;
  
  /**
   * Initialization tools
   **/
  BKN(List edgeCountList, NumericMatrix input_pmat);
  void ingestEdges(List lst);
  
  /**
   * Likelihood tools
   **/
  
  double node_llk(int i);
  double llk();
  double meanEdges(int i, int j);
  
  
  /** 
   * Fields + methods for BKN's algorithm
   **/
  vec<vec<double> >qList;
  vec<double> qSqrtSums;
  Mat theta_mat;
  Mat QMat;
  
  void update_qSqrtSums();
  void update_qSqrtSums(int i);
  
  void updateQ(int);
  void updateQ();
  double qijz(int i, int j, int k, double meanEdges);

  void update_ti(int i);
  void one_em();
  void par_one_em();
  List em(int max_it, double tol, bool par, int type);

  /**
   * LCN-style EM
   **/
  void update_ti2(int i);
  void one_em2();
  vec<double> thetaColSums;
  
  /**
   * Querying tools
   **/  
  double expectedDegree(int i);
  NumericMatrix get_theta();
};

void BKN::one_em2(){
  thetaColSums = theta_mat.colSums();
  for(int i = 0; i < nNodes; i++){ update_ti2(i); }
}

void BKN::update_ti2(int i){
  int n_edges = edgeCounts[i].size();
  vec<double> denom(n_edges, 0.0);
  int this_j;
  double tik, tjk;
  for(int k = 0; k < dim; k++){
    tik = theta_mat(i,k);
    for(int ii = 0; ii < n_edges; ii++){
      this_j = edgeCounts[i][ii][0];
      tjk = theta_mat(this_j,k);
      denom[ii] += tik * tjk;
    }
  }
  
  double tik_new;
  for(int k = 0; k < dim; k++){
    tik_new = 0.0;
    tik = theta_mat(i,k);
    double this_A;
    for(int ii = 0; ii < n_edges; ii++){
      this_j = edgeCounts[i][ii][0];
      tjk = theta_mat(this_j,k);
      this_A = edgeCounts[i][ii][1];
      tik_new += this_A * tjk / denom[ii];
    }
    theta_mat(i,k) = tik_new * tik / thetaColSums[k];
  }
}


List BKN::em(int max_it, double tol, bool par, int type){
  double err = tol + 1.0;
  int iter = 0;
  Mat theta_old = theta_mat.copy();
  while( (iter < max_it) & (err > tol) ){
    if(type == 1){
      if(!par){ one_em(); }
      else{ par_one_em(); }
    }
    if(type == 2){
      one_em2();
    }
    err = compute_err(theta_old, theta_mat);
    theta_old = theta_mat.copy();
    iter++;
  }
  if(iter == max_it){ Rcout << "Warning: max iterations reached\n"; }
  List ans = List::create(Named("err") = err, 
                          Named("its") = iter);
  return(ans);
}

double BKN::expectedDegree(int i){
  i--;
  if(i < 0 || i >= nNodes){ stop("invalid i"); }
  double ans = 0.0;
  for(int j = 0; j < nNodes; j++){
    ans += meanEdges(i,j);
  }
  return(ans);
}

NumericMatrix BKN::get_theta(){
  NumericMatrix ans = theta_mat.getNumMat();
  return(ans);
}

BKN::BKN(List edgeCountList, NumericMatrix input_pmat){
  theta_mat = Mat(input_pmat);
  dim = input_pmat.cols();
  nNodes = input_pmat.rows();
  ingestEdges(edgeCountList);
  pTol = 0.00000001;
  tol = 0.0001;
  QMat = theta_mat.copy();
}

void BKN::ingestEdges(List lst){
  nNodes = lst.length();
  edgeCounts.resize(nNodes);
  qList.resize(nNodes);
  int this_length, j, A;
  for(int i = 0; i < nNodes; i++){
    IntegerMatrix theseEdges = lst[i];
    this_length = theseEdges.nrow();
    edgeCounts[i].resize(this_length);
    qList[i].resize(this_length);
    for(int ii = 0; ii < this_length; ii++){
      edgeCounts[i][ii].resize(2);
      // Switching from 1 to zero index
      j = theseEdges(ii,0) - 1;
      A = theseEdges(ii,1);
      edgeCounts[i][ii][0] = j;
      edgeCounts[i][ii][1] = A;
    }
  }
}

void BKN::update_qSqrtSums(int i){
  int this_n, this_j, this_A;
  double meanSum;
  this_n = edgeCounts[i].size();
  for(int k = 0; k < dim; k++){ QMat(i,k) = 0.0; }
  for(int ii = 0; ii < this_n; ii++){
    this_j = edgeCounts[i][ii][0];
    this_A = edgeCounts[i][ii][1];
    meanSum = meanEdges(i, this_j);
    for(int k = 0; k < dim; k++){
      QMat(i, k) = QMat(i,k) + this_A * qijz(i, this_j, k, meanSum);
    }
  }
}


void BKN::update_qSqrtSums(){
  qSqrtSums.resize(dim);
  for(int k = 0; k < dim; k++){ qSqrtSums[k] = 0.0; }
  for(int i = 0; i < nNodes; i++){ update_qSqrtSums(i); }
  for(int i = 0; i < nNodes; i++){
    for(int k = 0; k < dim; k++){
      qSqrtSums[k] += QMat(i,k);
    }
  }
  for(int k = 0; k < dim; k++){
    qSqrtSums[k] = sqrt(qSqrtSums[k]);
  }
}

void BKN::update_ti(int i){
  int this_n = edgeCounts[i].size();
  int this_j, this_A;
  double these_meanEdges;
  vec<double> temp_meanEdges(dim);
  for(int k = 0; k < dim; k++){ temp_meanEdges[k] = 0.0; }
  for(int ii = 0; ii < this_n; ii++){
    this_j = edgeCounts[i][ii][0];
    this_A = edgeCounts[i][ii][1];
    these_meanEdges = meanEdges(i, this_j);
    for(int k = 0; k < dim; k++){
      temp_meanEdges[k] += qijz(i, this_j, k, these_meanEdges) * this_A;
    }
  }
  for(int k = 0; k < dim; k++){
    theta_mat(i,k) = temp_meanEdges[k] / qSqrtSums[k];
  }
}


double BKN::meanEdges(int i, int j){
  double ans = 0.0;
  for(int k = 0; k < dim; k++){
    ans += theta_mat(i,k) * theta_mat(j,k);
  }
  return(ans);
}

double BKN::qijz(int i, int j, int k, double meanEdges){
  double ans = theta_mat(i,k) * theta_mat(j,k) / meanEdges;
  return(ans);
}

void BKN::updateQ(int i){
  int j, n_these_edges;
  n_these_edges = edgeCounts[i].size();
  for(int ii = 0; ii < n_these_edges; ii++){
    j = edgeCounts[i][ii][0];
    qList[i][ii] = meanEdges(i,j);
  }
}

void BKN::updateQ(){
  for(int i = 0; i < nNodes; i++){ updateQ(i); }
}



double BKN::node_llk(int i){
  int nEdges = edgeCounts[i].size();
  vec<double> these_edgeCounts(nNodes, 0);
  for(int ii = 0; ii < nEdges; ii++){
    these_edgeCounts[ edgeCounts[i][ii][0] ] = edgeCounts[i][ii][1];
  }
  double this_meanEdges, this_x;
  double ans = 0.0;
  for(int j = 0; j < nNodes; j++){
    if(i == j){ continue; }
    this_meanEdges = meanEdges(i,j);
    this_x = these_edgeCounts[j];
    ans += R::dpois(this_x, this_meanEdges, true);
  }
  return(ans);
}

double BKN::llk(){
  double ans = 0.0;
  for(int i = 0; i < nNodes; i++){
    ans += node_llk(i);
  }
  return(ans);
}



/*
 * Parallel Tools
 * 
 */

struct ThetaUpdater : public RcppParallel::Worker{
  BKN* bkn_mod;
  void operator()(size_t begin, size_t end){
    for(int i = begin; i < end; i++){ bkn_mod->update_ti(i); }
  }
  ThetaUpdater(BKN* ptr){
    bkn_mod = ptr;
  }
};

struct QMat_Updater : public RcppParallel::Worker{
  BKN* bkn_mod;
  QMat_Updater(BKN* ptr){ bkn_mod = ptr; }
  
  void operator()(size_t begin, size_t end){
    for(int i = begin; i < end; i++){ 
      bkn_mod->update_qSqrtSums(i); 
    }
  }
};

struct QUpdater : public RcppParallel::Worker{
  BKN* bkn_mod;
  void operator()(size_t begin, size_t end){
    for(int i = begin; i < end; i++){ bkn_mod->updateQ(i); }
  }
  QUpdater(BKN* ptr){ bkn_mod = ptr; }
};

void BKN::one_em(){
  updateQ();
  update_qSqrtSums();
  for(int i = 0; i < nNodes; i++){ update_ti(i); }
}

void BKN::par_one_em(){
  QUpdater qup(this);
  QMat_Updater qsqup(this);
  ThetaUpdater tup(this);
  
  RcppParallel::parallelFor(0, nNodes, qup);
  
  qSqrtSums.resize(dim);
  for(int k = 0; k < dim; k++){ qSqrtSums[k] = 0.0; }
  RcppParallel::parallelFor(0, nNodes, qsqup);
  for(int i = 0; i < nNodes; i++){
    for(int k = 0; k < dim; k++){
      qSqrtSums[k] += QMat(i,k);
    }
  }
  
  for(int k = 0; k < dim; k++){ qSqrtSums[k] = sqrt(qSqrtSums[k]); }
  
  RcppParallel::parallelFor(0, nNodes, tup);
}

