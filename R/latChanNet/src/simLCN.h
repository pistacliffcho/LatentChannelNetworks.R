#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

double edgeProb(int i, int j, NumericMatrix p_mat){
  int nRows = p_mat.rows();
  if( i < 0 | j < 0 | i >= nRows | j >= nRows){
    stop("invalid indices");
  }
  double p_noEdge = 1.0;
  int nCols = p_mat.cols();
  for(int k  = 0; k < nCols; k++){
    p_noEdge *= (1.0 - p_mat(i,k) * p_mat(j,k));
  }
  double ans = 1.0 - p_noEdge;
  return(ans);
}

NumericMatrix simLCN(NumericMatrix p_mat){
  int nRows = p_mat.rows();
  int nCols = p_mat.cols();
  double this_p;
  for(int i = 0; i < nRows; i++){
    for(int j = 0; j < nCols; j++){
      this_p = p_mat(i,j);
      if(this_p < 0 | this_p > 1){
        stop("invalid probability supplied");
      }
    }
  }
  std::vector<int> i_inds(0);
  std::vector<int> j_inds(0);
  double this_prob;
  for(int i = 0; i < nRows; i++){
    for(int j = 0; j < i; j++){
      this_prob = edgeProb(i,j, p_mat);
      if(this_prob > runif(1, 0, 1)[0]){
        i_inds.push_back(i + 1);
        j_inds.push_back(j + 1);
      }
    }
  }
  int nEdges = i_inds.size();
  NumericMatrix ans(nEdges, 2);
  for(int i = 0; i < nEdges; i++){
    ans(i,0) = i_inds[i];
    ans(i,1) = j_inds[i];
  }
  return(ans);
}
