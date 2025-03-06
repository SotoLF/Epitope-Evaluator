# example_helper.R - Helper functions for example modules
#
# This module contains helper functions used by the example modules

#' Get intersect members for Up-Set plots
#'
#' @param data Data frame with peptide information
#' @param columns Columns to check for intersection
#' @return Matrix with intersection data
get_intersect_members <- function(data, columns) {
  if(length(columns) == 1) {
    return(as.matrix(data[,columns, drop = FALSE]))
  } else {
    result <- matrix(nrow = nrow(data), ncol = length(columns))
    colnames(result) <- columns
    
    for(i in 1:length(columns)) {
      result[,i] <- as.numeric(data[,columns[i]] == 1)
    }
    
    return(result)
  }
}