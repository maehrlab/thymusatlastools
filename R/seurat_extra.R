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
FetchDataZeroPad = function( dge, vars.all, warn = T, ... ){
  vars.all = vars.all[complete.cases(vars.all)]
  avail = intersect( vars.all, AvailableData( dge ) )
  unavail = setdiff( vars.all, AvailableData( dge ) )
  if(warn && length(unavail) > 0){
    warning("Some variables are not available. Returning zeroes.\n")
  }
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
#'
#' @export
#'
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

#' Merge two Seurat objects.
#'
#' By default, preserves any column of @data.info shared between both objects.  
#' You can also specify what variables to keep. They will be added to data.info in 
#' the output, warning and padding with zeroes if either object lacks any var in vars.keep.
#'  
#' @export
#' 
SeuratMerge = function( dge1, dge2, vars.keep = intersect(names(dge1@data.info), 
                                                          names(dge2@data.info) ) ){
  dge_all = list( dge1 = deseuratify_raw_data( dge1 ), 
                  dge2 = deseuratify_raw_data( dge2 ) ) %>%
    dge_merge_list %>% seuratify_thy_data 
  characterize_factor = function(x){ if(is.factor(x)) as.character(x) else x }
  all_factors_to_char = function(X) data.frame(lapply(X, characterize_factor), stringsAsFactors=FALSE)
  if( length(vars.keep) > 0 ){
    preserved_metadata = rbind( FetchDataZeroPad( dge1, vars.keep ) %>% all_factors_to_char, 
                                FetchDataZeroPad( dge2, vars.keep ) %>% all_factors_to_char )
    preserved_metadata %<>% as.data.frame
    rownames(preserved_metadata) = c( rownames(dge1@data.info),rownames(dge2@data.info)) 
    dge_all %<>% AddMetaData( preserved_metadata )     
  }
  return(dge_all)
}

#' Test for markers flexibly from a Seurat object.
#'
#' Calls FindMarkers with extra features.
#'
#' @param ident.use Fetched via FetchData to define the groups being tested. Should obey 
#' @param test.use Passed into FindMarkers unless it is "binomial_batch", in which case 
#'   it uses approximate p-values based on a binomial glmm with a random effect for batch (1|orig.ident). 
#'
#' All other parameters are passed into FindMarkers unless test.use=="binomial_batch", in 
#' which case I attempt to match the behavior of FindMarkers.
#' 
#' Output contains an extra column for q-values from p.adjust(..., method="fdr").
#'
#' @export
#'
FindMarkersFlex = function( object,
                            ident.use, ident.1, 
                            ident.2 = object %>% FetchData(ident.use) %>% extract2(1) %>% unique %>% setdiff(ident.1),
                            order_by_var = "avg_diff",
                            thresh.use = 0.25, 
                            test.use = "binomial_batch",
                            genes.use = object@data %>% rownames,
                            min.pct = 0.1, ... ){
  
  genes.use %<>% intersect(AvailableData(object))
  object %<>% Seurat::SetIdent( ident.use = as.character(FetchData( object, ident.use )[[1]]) )
  predicate = paste0( ident.use, " %in% c( '", ident.1, "', '", paste0(ident.2, collapse = "', '"), "' )" )
  object %<>% SubsetDataFlex( vars.use = ident.use, predicate )
  if( test.use == "binomial_batch" ){
    cat(" \n Computing summaries... \n")
    x = data.frame( gene = genes.use, stringsAsFactors = F )
    rownames( x ) = x$gene
    group_means = aggregate.nice( x  = FetchData(object, genes.use), 
                                  by = FetchData(object, ident.use), 
                                  FUN = mean ) %>% t
    group_pcts  = aggregate.nice( x  = FetchData(object, genes.use), 
                                  by = FetchData(object, ident.use), 
                                  FUN = prop_nz ) %>% t
    x$avg_diff  = group_means[, ident.1] - group_means[, ident.2]
    x$pct.1 = group_pcts[, ident.1]
    x$pct.2 = group_pcts[, ident.2]
    x = subset(x, abs(avg_diff) > thresh.use & ( pct.1 > min.pct | pct.2 > min.pct ) )
    cat(" Computing p-values... \n")
    get_p = function( gene ) {
      data = FetchData(object, c(gene, ident.use, "orig.ident"))
      data[[gene]]  %<>% is_greater_than(0)
      data[[ident.1]] = data[[ident.use]] %in% ident.1
      colnames(data) = make.names(colnames(data))
      mod = lme4::glmer(formula = paste0( colnames(data)[1], " ~ (1|orig.ident) + ", colnames(data)[4] ) , 
                        family = "binomial", data = data )
      mod_p = car::linearHypothesis( mod, hypothesis.matrix = paste0( make.names(ident.1), "TRUE = 0" ) )
      cat(".")
      return( mod_p$`Pr(>Chisq)`[[2]] )
    }
    x$p.value = parallel::mclapply( x$gene, 
                                    function(s) {
                                      tryCatch(get_p(s), error = function(e) NA) 
                                    }) %>% simplify2array()
    failures = is.na(x$p.value)
    cat("    ", sum( failures ), " failed tests out of ", nrow(x), 
        ". Setting failures to 1 for conservative FDR control. \n" )
    x$p.value[failures] = 1
  } else {
    x = Seurat::FindMarkers( object, ident.1 = ident.1, ident.2 = ident.2,                              
                             test.use = test.use,
                             genes.use = genes.use,   
                             thresh.use = thresh.use, 
                             min.pct = min.pct, ... )
  }
  x %<>% (plyr::rename)(c("p_val" = "p.value"))
  if( !is.null( x$p.value ) ){
    x$q.value = p.adjust( x$p.value, method = "fdr" )
  }
  x = x[order(x[[order_by_var]], decreasing = T), ]
  x$gene = rownames(x)
  return( x )
}


#' Sanitize gene names via `make.names`
#'
#' @export 
#'
SanitizeGenes = function( dge ){
  rownames( dge@raw.data )   %<>% make.names
  rownames( dge@data )       %<>% make.names
  rownames( dge@scale.data ) %<>% make.names
  names(    dge@var.genes )  %<>% make.names
  return( dge )
}


## ------------------------------------------------------------------------

#' Make a FACS-like plot from a single-cell rna-seq dataset.
#'
#' @param dge Seurat object
#' @param gene1 Horizontal axis on plot mimics this gene. Character, usually length 1 but possibly longer.
#' @param gene2 Vertical axis on plot mimics this gene. Character, usually length 1 but possibly longer. 
#' @param genesets_predetermined If FALSE, plot the sum of many genes similar to gene1 instead of gene1 alone (same 
#' for gene2). See ?get_similar_genes. If TRUE, plot the sum of only the genes given.
#' @param return_val If "all", returns a list with several internal calculations revealed.
#' If "plot", returns just a ggplot object. If "seurat", returns a Seurat object with gene scores added. 
#' @param cutoffs If given, divide plot into four quadrants and annotate with percentages. Numeric vector of length 2.
#' @param num_genes_add Each axis shows a simple sum of similar genes. This is how many (before removing overlap). Integer.
#'
#' @export 
#'
TACS = function( dge, gene1, gene2, genesets_predetermined = F, 
                 return_val = "plot", 
                 num_genes_add = 100, facet_by = NULL, cutoffs = NULL ){
  
  # Get gene sets to average
  if(genesets_predetermined){
    g1_similar = gene1
    g2_similar = gene2
  } else {
    g1_similar = get_similar_genes(dge, gene1, num_genes_add) %>% c( gene1, . )
    g2_similar = get_similar_genes(dge, gene2, num_genes_add) %>% c( gene2, . ) 
    shared = intersect(g1_similar, g2_similar)
    g1_similar %<>% setdiff(shared)
    g2_similar %<>% setdiff(shared)
  }
  
  # Average gene sets to get scores
  g1_score = rowMeans(FetchData(dge, g1_similar))
  g2_score = rowMeans(FetchData(dge, g2_similar))
  g1_score_name = paste0(gene1[1], "_score")
  g2_score_name = paste0(gene2[1], "_score")
  
  #Add scores as metadata. Extract with faceting var into plotting data.
  dge %<>% AddMetaData(g1_score, col.name = g1_score_name)
  dge %<>% AddMetaData(g2_score, col.name = g2_score_name)
  plot_df = FetchData(dge, c(g1_score_name, g2_score_name, facet_by))
  
  # Form plot
  p = ggplot(plot_df) + 
    geom_point( aes_string( x=g1_score_name, y=g2_score_name ) ) 
  p = p + expand_limits(y=0, x=0)
  # Facet to taste
  if(!is.null(facet_by)) {
    p = p + facet_wrap(as.formula(paste0("~", facet_by)))
  }
  
  # Add quadrants and percentages
  if( !is.null(cutoffs)){
    p = p + geom_vline(data = data.frame(xint=cutoffs[1]), aes(xintercept=xint))
    p = p + geom_hline(data = data.frame(yint=cutoffs[2]), aes(yintercept=yint))
    plot_df[[facet_by]] %<>% droplevels
    percentages = plot_df
    percentages[[g1_score_name]] %<>% is_greater_than(cutoffs[1])
    percentages[[g2_score_name]] %<>% is_greater_than(cutoffs[2])
    percentages %<>% table %>% melt 
    if(!is.null(facet_by)) {
      percentages = percentages[order(percentages[[facet_by]]), ]
      for( facet_level in seq_along(unique(plot_df[[facet_by]]))){
        percentages$value[1:4 + 4*(facet_level-1)] %<>% percentify()
      }
    } else {
      percentages$value %<>% percentify()
    }
    
    # Form annotation DF with correct facet and attempting sensible placement of percentages
    for( i in seq_along(percentages$value)){
      annot_df = data.frame(
        x = ifelse( percentages[i, g1_score_name], cutoffs[1]*2, cutoffs[1]*0.35) ,
        y = ifelse( percentages[i, g2_score_name], cutoffs[2]*2, cutoffs[2]*0.25) ,
        label = paste0( round(percentages$value[i], 1), "%") )
      if(!is.null(facet_by)) {
        annot_df[[facet_by]] = percentages[i, facet_by]
      }
      p = p + geom_text( data = annot_df, aes(x=x,y=y,label=label) )                
    }
  } else {
    percentages = NULL
  }
  if( return_val == "all" ){
    return( list( plot = p, 
                dge = dge, 
                score_names = c( g1_score_name, g2_score_name ), 
                genesets = list( g1_similar, g2_similar ),
                plot_data = plot_df,
                percentages = percentages ) )
  } else if( return_val == "seurat" ){
    return(dge)
  } else if( return_val == "plot" ){
    return( p )
  } else {
    stop(" return_val should be 'all', 'seurat', or 'plot'. ")
  }
}


