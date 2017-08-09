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


#' Subset data flexibly from a Seurat object.
#'
#' @param dge Seurat object
#' @param vars.use Variables to fetch for use in `predicate`.
#' @param predicate String to be parsed into an R expression and evaluated as an arg to `base::subset`.
#' @details Calls FetchData, subsets it as a dataframe using base::subset, and 
#' subsets the Seurat object correspondingly (using the df rownames.)
#' @export
SubsetDataFlex = function( dge, vars.use, predicate, zeropad = TRUE ){
  if( typeof(predicate) != "character"){
    print("predicate should be a character vector. It will be parsed into `subset` as an R expression.")
  }
  if(zeropad){
    df = FetchDataZeroPad(dge, vars.use) 
  } else {
    df = FetchData(dge, vars.use) 
  }
  cu = df %>% subset(eval(parse(text=predicate))) %>% rownames
  return( SubsetData(dge, cells.use = cu) )
}


#' Test for markers flexibly from a Seurat object.
#'
#' @param object Seurat object.
#' @param ident.use Variables to fetch for use in `predicate`.
#' @param ... Passed to Seurat::FindMarkers.
#' @details Calls Seurat::FindMarkers after replacing object@ident with 
#' `as.character(FetchData( object, ident.use )[[1]])`.
#' @export
FindMarkersFlex = function( object, ident.use, order_by = "avg_diff", ... ){
  object %<>% Seurat::SetIdent( ident.use = as.character(FetchData( object, ident.use )[[1]]) )
  x =  Seurat::FindMarkers( object, ... )
  x = x[order(x[[order_by]], decreasing = T), ]
  return( x )
}


