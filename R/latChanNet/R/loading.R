#' @useDynLib latChanNet
#' @import ComplexHeatmap grid circlize Rcpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom mltools auc_roc

loadModule("LCN", TRUE)
loadModule("BKN", TRUE)
