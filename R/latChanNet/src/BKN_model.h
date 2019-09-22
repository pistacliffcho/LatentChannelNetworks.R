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


/**
 * Class for missing edges
 * j_mat: index for looking up parameters
 * j_lookup:index for looking up on edgeCountList
 **/
class MissingEdge{
public:
  int j_mat;
  int j_lookup;
};


/**
 * Class for information on edges
 * Because model is Poisson, want both connected node ID and count
 **/
class Edge{
public:
  int j;
  double cnt;
};

/**
 * Class for BKN model
 **/
class BKN{
public:
  int nNodes;
  int dim;
  double tol, pTol;
  vec<vec<Edge> >edgeCounts;
  vec<vec<MissingEdge> >unknownEdges;

  vec<vec<int> >posInds;
  
  /**
   * Initialization tools
   **/
  
  BKN(List edgeCountList, 
      NumericMatrix input_pmat, 
      List missingList);
  void ingestEdges(List lst);
  void ingestUnknownEdges(List unknownEdges);
  void set_theta(NumericMatrix theta);
  void resize(NumericMatrix p);
  
    
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
  void imputeMissingEdges(int i);
  void imputeMissingEdges(){
    for(int i = 0; i < nNodes; i++){ imputeMissingEdges(i); }
  }
  List em(int max_it, double tol, double pTol, bool par);

  /**
   * Querying tools
   **/  
  double expectedDegree(int i);
  NumericMatrix get_theta();
};

/**
 * EM call: controls which EM method used
 **/

List BKN::em(int max_it, double tol, 
             double ptol, bool par){
  pTol = ptol;
  double err = tol + 1.0;
  int iter = 0;
  Mat theta_old = theta_mat.copy();
  while( (iter < max_it) & (err > tol) ){
    R_CheckUserInterrupt();
    setPosInds(posInds, theta_mat);
    imputeMissingEdges();
    if(!par){ one_em(); }
    else{ par_one_em(); }
    err = compute_err(theta_old, theta_mat);
    theta_old = theta_mat.copy();
    iter++;
  }
//  if(iter == max_it){ Rcout << "Warning: max iterations reached\n"; }
  List ans = List::create(Named("err") = err, 
                          Named("its") = iter);
  return(ans);
}



/**
 * EM as presented in BKN 2011
 **/

void BKN::one_em(){
  updateQ();
  update_qSqrtSums();
  for(int i = 0; i < nNodes; i++){ update_ti(i); }
}


void BKN::update_qSqrtSums(int i){
  int this_n, this_j, this_A, k;
  double meanSum;
  vec<int>* pos_vec = &posInds[i];
  int nPos = (*pos_vec).size();
  this_n = edgeCounts[i].size();
  for(int k = 0; k < dim; k++){ QMat(i,k) = 0.0; }
  for(int ii = 0; ii < this_n; ii++){
    this_j = edgeCounts[i][ii].j;
    this_A = edgeCounts[i][ii].cnt;
    meanSum = meanEdges(i, this_j);
    for(int k_ind = 0; k_ind < nPos; k_ind++){
      k = (*pos_vec)[k_ind];
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
  int this_j, this_A, k;
  vec<int>* pos_vec = &posInds[i];
  int nPos = pos_vec->size();
  double these_meanEdges, new_t;
  vec<double> temp_meanEdges(dim);
  for(int k = 0; k < dim; k++){ temp_meanEdges[k] = 0.0; }
  for(int ii = 0; ii < this_n; ii++){
    this_j = edgeCounts[i][ii].j;
    this_A = edgeCounts[i][ii].cnt;
    these_meanEdges = meanEdges(i, this_j);
    for(int k_ind = 0; k_ind < nPos; k_ind++){
      k = (*pos_vec)[k_ind];
      temp_meanEdges[k] += qijz(i, this_j, k, these_meanEdges) * this_A;
    }
  }
  for(int k_ind = 0; k_ind < nPos; k_ind++){
    k = (*pos_vec)[k_ind];
    if(qSqrtSums[k] < 0.0001){
      theta_mat(i,k) = 0.0;
      continue;
    } 
    new_t = temp_meanEdges[k] / qSqrtSums[k];
    if(new_t < pTol){ new_t = 0.0; }
    theta_mat(i,k) = new_t;
  }
}


double BKN::meanEdges(int i, int j){
  checkInd(i, nNodes);
  checkInd(j, nNodes);
  double ans = 0.0;
  vec<int>* pos_vec;
  if(posInds[i].size() > posInds[j].size()){ pos_vec = &posInds[j]; }
  else{ pos_vec = &posInds[i]; }
  int nPosInds = pos_vec->size();
  int k;
  for(int k_ind = 0; k_ind < nPosInds; k_ind++){
    k = (*pos_vec)[k_ind];
    ans += theta_mat(i,k) * theta_mat(j,k);
  }
  return(ans);
}

double BKN::qijz(int i, int j, int k, double meanEdges){
  if(meanEdges < 0.001){
    double dk = dim;
    double ans = 1.0 / dk;
    return(ans);
  }
  double ans = theta_mat(i,k) * theta_mat(j,k) / meanEdges;
  return(ans);
}

void BKN::updateQ(int i){
  int j, n_these_edges;
  n_these_edges = edgeCounts[i].size();
  for(int ii = 0; ii < n_these_edges; ii++){
    j = edgeCounts[i][ii].j;
    qList[i][ii] = meanEdges(i,j);
  }
}

void BKN::updateQ(){
  for(int i = 0; i < nNodes; i++){ updateQ(i); }
}

void BKN::imputeMissingEdges(int i){
  int nMissing = unknownEdges[i].size();
  int j_mat, j_lookup;
  double this_mean;
  for(int ii = 0; ii < nMissing; ii++){
    j_mat = unknownEdges[i][ii].j_mat;
    this_mean = meanEdges(i, j_mat);
    j_lookup = unknownEdges[i][ii].j_lookup;
    edgeCounts[i][j_lookup].cnt = this_mean;
  }
}


/**
 * Model querying
 **/

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

/**
 * Model initialization
 **/

void BKN::resize(NumericMatrix p){
  int new_nrow = p.rows();
  if(new_nrow != nNodes){ stop("wrong number of rows"); }
  dim = p.cols();
  theta_mat = Mat(p);
  setPosInds(posInds, theta_mat);
  QMat = theta_mat.copy();
}


BKN::BKN(List edgeCountList, 
         NumericMatrix input_pmat, 
         List missingList){
  theta_mat = Mat(input_pmat);
  setPosInds(posInds, theta_mat);
  dim = input_pmat.cols();
  nNodes = input_pmat.rows();
  ingestEdges(edgeCountList);
  ingestUnknownEdges(missingList);
  pTol = 0.00000001;
  tol = 0.0001;
  QMat = theta_mat.copy();
}

/**
 * Ingest edge list
 * Note: unknown edges must be on this list
 **/
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
      // Switching from 1 to zero index
      j = theseEdges(ii,0) - 1;
      A = theseEdges(ii,1);
      edgeCounts[i][ii].j = j;
      edgeCounts[i][ii].cnt = A;
    }
  }
}

/**
 * Ingest unknown edges
 * All edges *must* be included on full edgeList
 **/
void BKN::ingestUnknownEdges(List RunknownEdges){
  if(RunknownEdges.size() != nNodes){ 
    stop("Input unknownEdge list wrong length"); 
  }
  unknownEdges.resize(nNodes);
  int this_size, n_obs, j_ind, this_j_lookup;
  for(int i = 0; i < nNodes; i++){
    IntegerVector these_edges = RunknownEdges[i];
    this_size = these_edges.size();
    unknownEdges[i].resize(this_size);
    for(int ii = 0; ii < this_size; ii++){
      j_ind = these_edges[ii] - 1;
      if(j_ind < 0){ 
        Rcout << "i = " << i << " ii = " << ii <<"\n";
        stop("j_ind < 0"); 
      }
      unknownEdges[i][ii].j_mat = j_ind;
      n_obs = edgeCounts[i].size();
      this_j_lookup = -1;
      for(int iii = 0; iii < n_obs; iii++){
        if(edgeCounts[i][iii].j == j_ind){
          this_j_lookup = iii;
          break;
        }
      }
      if(this_j_lookup == -1){stop("Did not find matching edge in edgeCounts from edgeList!");}
      unknownEdges[i][ii].j_lookup = this_j_lookup;
    }
  }
}

void BKN::set_theta(NumericMatrix new_theta){
  if(new_theta.rows() != theta_mat.nRows || 
     new_theta.cols() != theta_mat.nCols){
    stop("Incorrect dimensions for new theta");
  }
  
  for(int i = 0; i < theta_mat.nRows; i++){
    for(int j = 0; j < theta_mat.nCols; j++){
      theta_mat(i,j) = new_theta(i,j);
    }
  }
  setPosInds(posInds, theta_mat);
}


/**
 * Log-likelihood tools
 **/

double BKN::node_llk(int i){
  // Filling in non-zero edges
  // Note: this includes unknown edges, which need to be skipped
  int nEdges = edgeCounts[i].size();
  int this_j;
  vec<double> these_edgeCounts(nNodes, 0);
  for(int ii = 0; ii < nEdges; ii++){
    this_j = edgeCounts[i][ii].j;
    these_edgeCounts[ this_j ] = edgeCounts[i][ii].cnt;
  }

  // Making vector indicating missingness
  int nMissingEdges = unknownEdges[i].size();
  vec<int> isMissing(nNodes, 0);
  for(int ii = 0; ii < nMissingEdges; ii++){
    this_j = unknownEdges[i][ii].j_mat;
    isMissing[ this_j ] = 1;
  }
  double this_meanEdges, this_x;
  double ans = 0.0;
  for(int j = 0; j < nNodes; j++){
    // Skipping unknown edges
    if(isMissing[j] == 1){ continue; }
    this_meanEdges = meanEdges(i,j);
    this_x = these_edgeCounts[j];
//    if(i == j){
//      ans += 2.0 * R::dpois(this_x, this_meanEdges, true);
//    }
//    else{
      ans += R::dpois(this_x, this_meanEdges, true);
//    }
  }
  ans = ans/2.0;
  return(ans);
}

double BKN::llk(){
  double ans = 0.0;
  for(int i = 0; i < nNodes; i++){
    ans += node_llk(i);
  }
  return(ans);
}



/**
 * Parallel Tools
 **/

struct ThetaUpdater : public RcppParallel::Worker{
  BKN* bkn_mod;
  void operator()(size_t begin, size_t end){
    for(size_t i = begin; i < end; i++){ bkn_mod->update_ti(i); }
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

