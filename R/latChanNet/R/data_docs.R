#' Email data for EU Univeristy
#' 
#' An email network with metadata for an EU university. 
#' Nodes are professors, edges indicate an email was sent between the two nodes. 
#' Metadata is the department that each node belong to. 
#' 
#' @format List with two objects:
#' \describe{
#'     \item{edgeList}{A 16706 x 2 matrix of edges}
#'     \item{meta}{A 1005 x 1 data frame indicating department of each node}
#' }
"email_data"