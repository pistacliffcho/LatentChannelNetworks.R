#include "latChan.h"
#include "BKN_model.h"

// LCN MODEL
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

// BKN MODEL
RCPP_MODULE(BKN){
  class_<BKN>("BKN")
  .constructor<List, NumericMatrix>("Args: EdgeCountList, Initial theta-mat")
  .method("llk", &BKN::llk)
  .method("one_em", &BKN::one_em)
  .method("get_theta", &BKN::get_theta);
}