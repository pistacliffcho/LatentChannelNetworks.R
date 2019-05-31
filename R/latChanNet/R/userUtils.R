#' @title Expected Channel Usage
#' @description Compute expected channel usage between two connected nodes.
#' @param i Index of first node
#' @param j Index of second node
#' @param lcn_mod LCN model
#' @details Computes the expected channel usage between two nodes 
#' \strong{conditional on the two nodes sharing an edge.}
#' @export
computeTheta = function(i, j, lcn_mod){
  # Converting R indices to C indices
  c_i = i - 1; c_j = j - 1
  ans = lcn_mod$computeTheta(c_i,c_j)
  return(ans)
}


#' @title Compute expected connections for node
#' @description Computes the estimated expected connections through each latent channel
#' @param i Index of node
#' @param lcn_mod LCN model
#' @export
computeExpConnects = function(i, lcn_mod){
  # Converting R index to C indices
  c_i = i - 1
  ans = lcn_mod$expectedConnections(c_i)
  names(ans) = paste("Channel", seq_along(ans))
  return(ans)
}