LatPlot = setRefClass("LatPlot", 
                      fields = c("node_plot", 
                                 "node_info", 
                                 "meta_plot", 
                                 "meta_info", 
                                 "has_meta"), 
                      methods = c("initialize", "show"))

LatPlot$methods(
  initialize = function(node_list, meta_list = NULL){
    node_plot <<- node_list$plot
    node_info <<- node_list$info
    
    if(is.null(meta_list)){
      has_meta <<- FALSE
      meta_plot <<- NULL
      meta_info <<- NULL
    }
    else{
      has_meta <<- TRUE
      meta_plot <<- meta_list$plot
      meta_info <<- meta_list$info
    }
  }
)

LatPlot$methods(
  show = function(){
    if(!has_meta){
      print(node_plot)
    }
    else{
      print(node_plot + meta_plot)
    }
  }
)

#' @title Build heatmap from model
#' @param mod LatClass object 
#' @param grp Vector of group categories for each node
#' @param minGrpSize Minimum size of group in both. Smaller groups put in "other"
#' @param prob_cols Colors for color gradient of probability range
#' @param greater_col Color for color gradient beyond 1
#' @param reorderChannels Should Channels be reorder by dependency on grp?
#' @param xlab X-axis label
#' @param ylab Y-asix label
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap
#' @export
plot_net = function(mod, 
                      grp = NULL, 
                      metanames = NULL,
                      minGrpSize = NULL,
                      subset = NULL,
                      name = "",
                      ...){
  
  # Node plot
  node_info = hm_node_info(mod, grp = grp, 
                           minGrpSize = minGrpSize,
                           subset = subset)
  node_p = nodeinfo2ComplexHeatmap(node_info)
  nodeList = list(plot = node_p, info = node_info)
  if(is.null(metanames)){ 
    ans = LatPlot$new(nodeList)
    return(ans) 
  }
  
  # Meta plot
  meta_info = hm_meta_info(mod, metanames = metanames, 
                           colFxn = node_info$colFxn, 
                           chanOrder = node_info$chan_order)
  
  meta_p = meta2ComplexHeatmap(meta_info)
  metaList = list(plot = meta_p, info = meta_info)
  ans = LatPlot$new(nodeList, metaList)
  return(ans)
}


# Make heatmap info about nodes
hm_node_info = function(mod, 
                        grp = NULL, 
                        minGrpSize = NULL,
                        prob_cols = c("black", "white", "orange", "red"),
                        greater_col = "purple",
                        plotChannelNumber = T,
                        xlab = NULL, ylab = NULL,
                        sortColumns = T,
                        subset = NULL,
                        name = "",
                        ...){
  pmat = mod$get_pars()$nodes
  if(is.null(ylab)) ylab = "Channels"
  
  # Keeping track of original channel names
  if(plotChannelNumber){ colnames(pmat) = 1:ncol(pmat) }
  
  if(!is.null(subset)){
    if(!is.null(grp)){ grp = grp[subset] }
    pmat = pmat[subset,]
  }
  
  if(!is.null(minGrpSize)){
    grp = as.character(grp)
    cnts = table(grp)
    is_too_small = cnts < minGrpSize
    too_small = which(is_too_small)
    for(i in seq_along(too_small)){
      this_ind = too_small[i]
      this_grp = names(cnts)[this_ind]
      grp[grp == this_grp] = "Other"
    }
  }
  
  # Function for making colors
  colFxn = latChanNet:::make_colFxn(prob_cols, 
                                    greater_col, 
                                    pmat)
  
  # If grp is provided, need to sort and break up rows by grp
  if(!is.null(grp)){
    if(is.null(xlab) ){ xlab = "Nodes by Group" }
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
    rownames(pmat_ord) = plot_grp_names
  
    # Sorting columns by variance of category means
    if(sortColumns & !is.null(grp)){
      if(length(unique(ord_grp)) != 1){
        chan_order = reorder_nodes(ord_grp, pmat_ord)
        pmat_ord = pmat_ord[,chan_order]
      }
      else{
        chan_order = 1:ncol(pmat_ord)
      }
    }
  
    is_grp_name = plot_grp_names != ""
    grp_names = plot_grp_names[is_grp_name]
    grp_locs = which(is_grp_name)
    pmat_use = pmat_ord
  }
  else{
    if(is.null(xlab) ){ xlab = "Nodes" }
    pmat_use = pmat
    rownames(pmat_use) = rep("", nrow(pmat_use) )
    breaks = NULL
    grp_names = NULL
    grp_locs = NULL
    ord_grp = NULL
    chan_order = 1:ncol(pmat_use)
  }
  
  ans = list(mat = pmat_use, 
             colFxn = colFxn, 
             breaks = breaks, 
             grp_names = grp_names, 
             grp_locs = grp_locs, 
             xlab = xlab, ylab = ylab, 
             name = name, 
             ord_grp = ord_grp, 
             chan_order = chan_order) 
  return(ans)
}

# Take heatmap info and make Complex Heatmap
nodeinfo2ComplexHeatmap = function(hmi, 
                                 width = grid::unit(6, "cm")){
  p = ComplexHeatmap::Heatmap(t(hmi$mat), 
                              cluster_rows = F, 
                              cluster_columns = F, 
                              col = hmi$colFxn, 
                              column_title = hmi$xlab, 
                              column_title_side = "bottom",
                              row_title_side = "right", 
                              row_title = hmi$ylab,
                              name = hmi$name, 
                              column_split = hmi$ord_grp)
  return(p)
}

# Make meta heatmap meta info
hm_meta_info = function(mod,
                        metanames = NULL,
                        colFxn, 
                        chanOrder){
  if(is.null(metanames)){ return() }
  if(is.null(mod$metanames)){ stop("Model does not contain metadata") }
  if(any(!metanames %in% metanames) ){
    stop("invalid metanames. See mod$metanames")
  }
  
  meta_pmat = mod$get_pars()$meta
  meta_pmat = meta_pmat[metanames, ]
  meta_pmat = meta_pmat[,chanOrder]
  
  ans = list(mat = meta_pmat, 
             colFxn = colFxn, 
             xlab = "Metadata Nodes")
  return(ans)
}

# Make complex heatmap from meta info
meta2ComplexHeatmap = function(mi, 
                               width = grid::unit(6, "cm")){
  p = ComplexHeatmap::Heatmap(t(mi$mat), 
                              cluster_rows = F, 
                              cluster_columns = F, 
                              col = mi$colFxn,
                              width = width,
                              column_title = mi$xlab, 
                              column_title_side = "bottom",
                              show_heatmap_legend = F)
  return(p)
}
