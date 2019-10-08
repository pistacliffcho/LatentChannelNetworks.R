#' Reorder channels by variance explained by group
#' @param grp Vector of group for each node
#' @param pmat matrix of parameters for nodes
#' @noRd
reorder_nodes = function(grp, pmat){
  ord_grp_char = as.character(grp)
  use_grps = !(ord_grp_char %in% c("Missing", "Other"))
  grp2use = ord_grp_char[use_grps]
  pmat2use = pmat[use_grps,]
  ans = col_var_order(pmat2use, grp2use)
  return(ans)
}


weighted_coefs = function(x, y){
  fit = lm(y ~ x)
  coefs = coef(fit)
  cnts = table(x)
  return(list(coefs, cnts) )
}

wsd = function(x, w){
  w_mean = sum(x * w) / sum(w)
  w_mean2 = sum(x * x * w) / sum(w)
  ans = sqrt(w_mean2 - w_mean^2)
  return(ans)
}

col_var = function(pmat, grp){
  grp_fact = factor(grp)
  ans = NULL
  for(i in 1:ncol(pmat) ){
    wc = weighted_coefs(grp_fact, pmat[,i])
    ans[i] = wsd(wc[[1]], wc[[2]])
  }
  return(ans)
}

col_var_order = function(pmat, grp){
  colVars = col_var(pmat, grp)
  colOrder = order(colVars, decreasing = T)
  return(colOrder)
}

#' @param prob_cols Colors for probabilities
#' @param greater_col Color for values greater than 1
#' @param pars Parameters to be plotted
#' @noRd
make_colFxn = function(prob_cols, greater_col, pars){
  # Make color gradient function
  prob_points = seq(from = 0, to = 1, 
                    length.out = length(prob_cols))
  if(max(pars) > 1){
    cols2use = c(prob_cols, greater_col)
    use_points = c(prob_points, max(pars))
  }else{
    cols2use = prob_cols
    use_points = prob_points
  }
  colFxn = circlize::colorRamp2(use_points, 
                                colors = cols2use)
  return(colFxn)
}
