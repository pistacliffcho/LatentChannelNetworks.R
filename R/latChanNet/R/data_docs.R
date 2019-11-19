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
#' @details 
#' Hao Yin, Austin R. Benson, Jure Leskovec, and David F. Gleich. "Local Higher-order Graph Clustering." In Proceedings of the 23rd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining. 2017.
#' 
#' J. Leskovec, J. Kleinberg and C. Faloutsos. Graph Evolution: Densification and Shrinking Diameters. ACM Transactions on Knowledge Discovery from Data (ACM TKDD), 1(1), 2007.
"email_data"