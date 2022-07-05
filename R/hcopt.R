#####
#This function performs hclust with optimal leaf ordering from S. Monti / D. Gusenleitner
#####
#' @export
hcopt <- function(d, HC=NULL, method = "ward.D", members = NULL){
  require("cba")
  if ( is.null(HC) ) {
    HC <- hclust(d,method=method,members=members)
  }
  #optimal leaf ordering
  ORD <- cba::order.optimal(d,merge=HC$merge)
  HC$merge <- ORD$merge
  HC$order <- ORD$order
  HC
}
