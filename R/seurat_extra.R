## ------------------------------------------------------------------------
#' Get available variable names (genes, identity classes, PCA embeddings, etc)
#'
#' @param object Seurat object
#' @return Returns a character vector of all eligible inputs to the `vars.all` argument of `FetchData`.
#' @export
AvailableData = function( object ){
  available_categorized = list( metadata = names( object@data.info ),
                                PCs = names(object@pca.x),
                                tsne = names(object@tsne.rot),
                                ICs = names(object@ica.rot),
                                genes = rownames( object@data ),
                                ident = "ident" )
  return( Reduce( f = union, x = available_categorized ) )
}


#' FetchData but with zeroes for unavailable genes
#'
#' @export
#' @param dge Seurat object
#' @param vars.all List of all variables to fetch. Missing entries are ignored.
#' @param ... Other arguments to pass to FetchData
#'
#' @details This function is stupid: if you ask for "PC1" and it's not available,
#' it will think you're asking for a non-expressed gene, so it will return zeroes.
FetchDataZeroPad = function( dge, vars.all, ... ){
  vars.all = vars.all[complete.cases(vars.all)]
  avail = intersect( vars.all, AvailableData( dge ) )
  unavail = setdiff( vars.all, AvailableData( dge ) )
  to_return  = FetchData( dge,  avail, ... ) 
  pad = as.data.frame( matrix(0,           
                              nrow = nrow( to_return ), 
                              ncol = length( unavail ),
                              dimnames = list( rownames( to_return ),               
                                               unavail) ) )
  to_return = cbind( to_return, pad )
  assertthat::are_equal( sort( vars.all ),   sort( colnames( to_return ) ) )   
  to_return = to_return[, vars.all, drop = F]
  assertthat::assert_that( is.data.frame( to_return ) )
  return( to_return )
}


