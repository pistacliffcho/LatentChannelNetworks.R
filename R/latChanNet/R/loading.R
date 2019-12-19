#' @useDynLib latChanNet
#' @import grid 
#' @import methods stats utils
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom Rcpp loadModule cpp_object_initializer

loadModule("LCN", TRUE)
loadModule("BKN", TRUE)
