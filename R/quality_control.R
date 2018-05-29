## ------------------------------------------------------------------------
  
#' Check for female-specific and male-specific transcripts
#' 
#' @param raw_dge numeric matrix of raw molecule counts with genes in the rownames.
#'
#' @export
#'
check_xist_pure = function( raw_dge, rep_name, results_path ){

  dir.create.nice( results_path )
  print( paste0( "Checking for Xist in ", rep_name ) )
  # # try the gene as given, plus Capitalized, plus in UPPER CASE. Else return zeroes.
  friendly_gene_get = function( gene ) {
    genes_available = rownames( raw_dge )
    gene_present       = gene             %in% genes_available 
    gene_present_cap   = Capitalize(gene) %in% genes_available
    gene_present_upper = toupper(gene)    %in% genes_available
    
    if( gene_present       ){ return( raw_dge[gene,             ] ) }
    if( gene_present_cap   ){ return( raw_dge[Capitalize(gene), ] ) }
    if( gene_present_upper ){ return( raw_dge[toupper(gene),    ] ) }
    return( rep( 0, ncol( raw_dge ) ) )
  }
  
  data(Y_genes, package = "thymusatlastools") 
  Y_genes = intersect( c(Y_genes, toupper( Y_genes ) ), rownames( raw_dge ) )
  sink( file = file.path( results_path, paste0(rep_name, "_Xist_vs_Y_genes.txt" ) ) )
  {
    if( length( Y_genes ) == 0 ){
      print( "No y genes detected." ) 
      x = table( any_xist = friendly_gene_get( "XIST" ) > 0 ) 
      x %>% print
    } else {
      print( "Y_genes:" )
      print( Y_genes )
      any_y_genes_expressed = Reduce( f = "+", x = lapply( X = Y_genes, FUN = friendly_gene_get ) ) > 0
      x = table( any_y = any_y_genes_expressed, any_xist = friendly_gene_get( "XIST" ) > 0 ) 
      x %>% print
      x %>% chisq.test %>% print    
    }
  }
  sink()
  
  prop_xist = mean( friendly_gene_get( "XIST" ) > 0 )
  if( prop_xist > 0 ){
    cat( paste( "Detected Xist in ", 
             percentify( prop_xist ), 
             "percent of cells in replicate", 
             rep_name, "\n" ) )
  }
}


#' Assemble very basic summary stats: total UMIs, genes, and cells.
#'
#' @export
#'
save_depth_stats = function(results_path, samples = NULL, metadata = get_metadata()) {
  metadata %<>% subset( files_available == "yes" )
  if(!is.null( samples )){
    metadata %<>% subset( Sample_ID %in% samples )
  }
  all_runs = metadata$Sample_ID %>% as.list
  names( all_runs ) = all_runs
  get1 = function(sample) load_thymus_profiling_data(sample)[[1]]
  dimsum = function(X) c(dim(X), sum(X))
  basic_stats = data.frame( Reduce( rbind, lapply( all_runs, function(s) {dimsum( get1( s ) )} ) ) )
  colnames(basic_stats) = c( "transcripts", "cells", "UMIs" )
  rownames(basic_stats) = basic_stats$Sample_ID = unlist(all_runs)
  write.table( basic_stats, file = file.path( results_path, "basics_stats.txt" ),
               quote = F, row.names = F, col.names = T, sep = "\t")
  return(basic_stats)
}

#' Screen for female embryos using simple stats.
#'
#' @export
#'
check_xist_all = function( results_path, samples = NULL, metadata = get_metadata() ){
  if(!is.null(samples)){
    metadata %<>% subset( Sample_ID %in% samples )
  }
  # # Load the data; check for xist versus y chromosome genes
  all_runs = subset(metadata, files_available == "yes", select = "Sample_ID", drop = T) 
  for( rep_name in all_runs ){
    check_xist_pure( raw_dge = load_thymus_profiling_data( sample_ids = rep_name, test_mode = F )[[1]], 
                                       rep_name = rep_name, 
                                       results_path = file.path( results_path, "Xist_check" ) ) 
  }
}

## ------------------------------------------------------------------------
#' Consistently label replicates
#'
#' @details We use the convention that the sample ID, which goes into the orig.ident variable,
#' looks like "e13_5_<irregular>". This helps form a consistent way of labeling replicates.
custom.make.unique = function(s, ...){
  label_dupes = function( s2, my_seq, seqgen = seq_along, sep = "_", do_singletons = F ) {
    if ( ( length(s2) == 1 ) && !do_singletons ) {
      s2
    }  else{
      paste(s2, seqgen(s2), sep = sep)
    }
  }
  f = function(x) label_dupes(x, ...)
  ave(as.character(s), s, FUN = f)
} 



#' Plot a dataset's tSNE embedding faceted by experimental metadata (default replicate and eday).
#'
#' @export
faceted_tsne = function( dge, results_path, inner_factor = "orig.ident", outer_factor = "eday", colour = outer_factor, 
                         width = 20, height = 12 ){
  X = Y = FetchData(dge, unique( c( inner_factor, outer_factor, colour, "tSNE_1", "tSNE_2" ) ) )
  Y[[inner_factor]] = Y[[outer_factor]] = "all"
  replace_with_int_rank = function (x) as.numeric(as.factor(x))
  collapse_nested_factor = function( inner, outer ){
    ave(as.character(inner), outer, FUN = replace_with_int_rank )
  }
  plot_df = rbind(X, Y)
  plot_df[[inner_factor]] = collapse_nested_factor(inner = plot_df[[inner_factor]], 
                                                   outer = plot_df[[outer_factor]])
    p = ggplot(plot_df) + ggtitle("Cells stratified by replicate") + 
      facet_grid(paste0(inner_factor, "~",  outer_factor)) + 
      geom_point( aes_string(x = "tSNE_1", y = "tSNE_2", colour = colour ), size = 0.2) 
  fp = file.path(results_path, "tsne_replicates.pdf")
  cat( paste( "Saving faceted tSNE to", fp ) )
  ggsave( filename = fp, 
          plot = p, 
          width = width, height = height)
  return(p)
}

