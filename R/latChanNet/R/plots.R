#' @title Build heatmap from model
#' @param lcn_mod LCN model 
#' @param grp Vector of group categories for each node
#' @param minGrpSize Minimum size of group in both. Smaller groups put in "other"
#' @param cols Colors for color gradient
#' @param reorderRows Should Channels be reorder by dependency on grp?
#' @param xlab X-axis label
#' @param ylab Y-asix label
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap
#' @export
heatmapLCN = function(lcn_mod, 
                      grp, 
                      minGrpSize = NULL,
                      cols = c("darkslateblue", "lightblue", 
                               "orange", "red"),
                      reorderRows = T,
                      xlab = "", ylab = "",
                      ...){
  pmat = lcn_mod$get_pmat()
  if(reorderRows){
    cnts = table(grp)
    biggestGrp_ind = which.max(cnts)
    biggestGrp = names(cnts)[biggestGrp_ind]
    
    biggestGrp_colSums = colSums(pmat[grp == biggestGrp,])
    new_col_order = order(biggestGrp_colSums, decreasing = T)
    pmat = pmat[,new_col_order]
  }
  
  
  if(!is.null(minGrpSize)){
    cnts = table(grp)
    is_too_small = cnts < minGrpSize
    too_small = which(is_too_small)
    for(i in seq_along(too_small)){
      this_ind = too_small[i]
      this_grp = names(cnts)[this_ind]
      grp[grp == this_grp] = "Other"
    }
  }
  
  ord = order(grp)
  ord_grp = grp[ord]
  pmat_ord = pmat[ord,]
  
  is_break = ord_grp[-1] != head(ord_grp, -1)
  is_break[1] = TRUE
  breaks = which(is_break)
  
  plot_grp_names = rep("", nrow(pmat_ord) ) 
  for(i in seq_len(length(breaks) - 1) ){
    brk_loc = round( ( breaks[i+1] + breaks[i] ) / 2)
    brk_name = ord_grp[brk_loc]
    plot_grp_names[brk_loc] = brk_name
  }
  last_loc = round( (tail(breaks,1) + length(ord_grp) )/2 )
  last_name = ord_grp[last_loc]
  plot_grp_names[last_loc] = last_name
  
  
  labels_row = seq_len(ncol(pmat_ord))
  colnames(pmat_ord) = labels_row
  rownames(pmat_ord) = plot_grp_names
  
  nColors = length(cols)
  colFxn = circlize::colorRamp2(seq(from = 0, to = max(pmat_ord), 
                                    length.out = nColors), 
                                colors = cols)
  
  p = ComplexHeatmap::Heatmap(t(pmat_ord), 
                          name = "",
                          cluster_rows = F, 
                          cluster_columns = F,
                          col = colFxn,
                          column_split = ord_grp,
                          column_title_side = "bottom",
                          column_title = ylab, 
                          row_title_side = "right",
                          row_title = xlab,
                          ...)  
 return(p) 
}
