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
  vec<vec<double> >qList;
  vec<double> qSqrtSums;
  vec<double> temp_meanEdges;
  Mat theta_mat;
  
  void update_qSqrtSums();
  
  void ingestEdges(List lst);
  void updateQ(int);
  void updateQ();
//  void parInitCache();
  BKN(List edgeCountList, NumericMatrix input_pmat);
  
  double meanEdges(int i, int j);
  double qijz(int i, int j, int k, double meanEdges);
  
  double node_llk(int i);
  double llk();
  
  void update_ti(int i);
  void one_em();
  List em(int max_it, double tol);
  
  double expectedDegree(int i);
  NumericMatrix get_theta();
};

List BKN::em(int max_it, double tol){
  double err = tol + 1.0;
  int iter = 0;
  Mat theta_old = theta_mat.copy();
  while( (iter < max_it) & (err > tol) ){
    one_em();
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


void BKN::one_em(){
  updateQ();
  update_qSqrtSums();
  for(int i = 0; i < nNodes; i++){ update_ti(i); }
}

void BKN::update_qSqrtSums(){
  qSqrtSums.resize(dim);
  int this_n, this_j, this_A;
  double meanSum;
  for(int k = 0; k < dim; k++){ qSqrtSums[k] = 0.0; }
  for(int i = 0; i <nNodes; i++){
    this_n = edgeCounts[i].size();
    for(int ii = 0; ii < this_n; ii++){
      this_j = edgeCounts[i][ii][0];
      this_A = edgeCounts[i][ii][1];
      meanSum = meanEdges(i, this_j);
      for(int k = 0; k < dim; k++){
        qSqrtSums[k] += this_A * qijz(i, this_j, k, meanSum);
      }
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
  temp_meanEdges.resize(dim);
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
