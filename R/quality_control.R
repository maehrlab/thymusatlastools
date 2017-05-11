## ------------------------------------------------------------------------
check_all_scRNA = function( results_path ){
  metadata = read.csv( PATH_TO_METADATA, header = T, stringsAsFactors = F )
  
  # # Check to see the data is actually all there
  metadata$file_exists_auto = metadata$dge_path %>% file.exists 
  existence_status_correct = metadata$file_exists_auto == ( metadata$files_available == "yes" )
  if( !all( existence_status_correct ) ){
    print( "For these samples, the sheet says the files are available but they actually aren't, or vice versa." )
    print(" By row:" )
    cat( which( !existence_status_correct ) )
    print(" By Sample_ID:" )
    cat( metadata$Sample_ID[ !existence_status_correct ] )
  }
  
  # # Check for duplicate sample ID's and file names
  if( any( duplicated( metadata$Sample_ID ) ) ){
    print( "You have some duplicate sample IDs.")
    print( metadata$Sample_ID[ duplicated( metadata$Sample_ID ) ] )
  }
  file_dupes = metadata$raw_data_path[ duplicated( metadata$raw_data_path ) ]
  file_dupes = setdiff( file_dupes, c("") )
  if( !all( file_dupes %in% c("") ) ){
    print( "You have some duplicate files listed.")
    print( file_dupes )
  }
  
  # # Load the data; check for xist versus y chromosome genes
  all_runs = as.list(metadata$Sample_ID[metadata$file_exists_auto] ); names( all_runs ) = all_runs
  for( rep_name in names( all_runs ) ){
    all_runs[[rep_name]] = load_thymus_profiling_data( sample_ids = rep_name, test_mode = F, convert_all_to_mouse = F )[[1]]
    check_xist( raw_dge = all_runs[[rep_name]], 
                rep_name = rep_name, 
                results_path = file.path( results_path, "Xist_check" ) ) 
  }
  
  # # count transcripts, cells, umis, and reads
  basic_stats = data.frame( Reduce( rbind, lapply( all_runs, dim ) ) )
  rownames(basic_stats) = names( all_runs )
  colnames(basic_stats) = c( "transcripts", "cells" )
  basic_stats$umis_in_dge = as.numeric( lapply( all_runs, sum ) )
  basic_stats$Sample_ID = rownames( basic_stats )
  
  # # Format basic stats for seamless merging into metadata
  all_sampleids_df = read.csv( PATH_TO_METADATA, header = T, stringsAsFactors = F )[, 1, drop = F]
  rownames( all_sampleids_df ) = all_sampleids_df$Sample_ID
  basic_headers = c( "transcripts", "cells", "umis_in_dge" )
  all_sampleids_df[basic_headers] = NA
  all_sampleids_df[basic_stats$Sample_ID, basic_headers]  = basic_stats[, basic_headers]
  write.table( all_sampleids_df, file = file.path( results_path, "basics_stats.txt" ),
               quote = F, row.names = F, col.names = T, sep = "\t")
}

# # This function checks for female-specific and male-specific transcripts. 
# # It accepts a numeric matrix of raw molecule counts with genes in the rownames.
check_xist = function( raw_dge, rep_name, results_path ){
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
  
  Y_genes = read.table( file.path( PATH_TO_TABLES, "chrY_genes.txt"), stringsAsFactors = F )[[2]] %>% unique 
  Y_genes = intersect( c(Y_genes, toupper( Y_genes ) ), rownames( raw_dge ) )
  sink( file = file.path( results_path, paste0(rep_name, "_Xist_vs_Y_genes.txt" ) ) )
  {
    print( "Y_genes:" )
    print( Y_genes )
    any_y_genes_expressed = Reduce( f = "+", x = lapply( X = Y_genes, FUN = friendly_gene_get ) ) > 0
    x = table( any_y = any_y_genes_expressed, any_xist = friendly_gene_get( "XIST" ) > 0 ) 
    x %>% print
    x %>% chisq.test %>% print
  }
  sink()
  
  percentify = function( prop ) { return( 100*round( prop, 4 ) ) }
  prop_xist = mean( friendly_gene_get( "XIST" ) > 0 )
  if( prop_xist > 0 ){
    cat( paste( "Detected Xist in ", 
             percentify( prop_xist ), 
             "percent of cells in replicate", 
             rep_name, "\n" ) )
  }
}

# # This function checks for consistency across replicates as returned by `get_data_by_replicates`.
# # It plots avg expression for each gene and proportion expressing each gene.
scatterplot_replicates = function( results_path ){
  data_by_replicate = get_data_by_replicates()
  for( sample_type in names( data_by_replicate ) ){
    sample_ids = data_by_replicate[[sample_type]]
    reps = load_thymus_profiling_data( sample_ids = sample_ids, convert_all_to_mouse = F )
    all_genes = Reduce( x = lapply( reps, rownames ),f = union )
    # Initialize empty arrays
    mean_expr_by_gene = as.data.frame( matrix( 0, ncol = length( sample_ids ), nrow = length( all_genes ) ) )
    colnames( mean_expr_by_gene ) = sample_ids
    rownames( mean_expr_by_gene ) = all_genes
    prop_expr_by_gene = mean_expr_by_gene
    # Plot expression by gene both as log1p mean expression and proportion expressing.
    for( rep_id in sample_ids ){
      mean_expr_by_gene[ rownames(reps[[rep_id]]), rep_id ] = rowMeans( log1p( reps[[rep_id]]   ) )
      prop_expr_by_gene[ rownames(reps[[rep_id]]), rep_id ] = rowMeans(        reps[[rep_id]] > 0 ) 
    }
    dir.create.nice( file.path( results_path, "rep_check_total" ) )
    pdf( file.path( results_path, "rep_check_total", paste0( sample_type, ".pdf" ) ) ) 
    {
      plot_pairs( mean_expr_by_gene, main = "Total expression by gene" )
    } 
    dev.off()
    
    dir.create.nice( file.path( results_path, "rep_check_prop" ) )
    pdf( file.path( results_path, "rep_check_prop", paste0( sample_type, ".pdf" ) ) ) 
    {
      plot_pairs( prop_expr_by_gene, main = "Proportion expressing each gene" )
    } 
    dev.off()
    
  }
}

# # Another way of checking for consistency across replicates.
# # Plot a dataset's tSNE embedding, but faceted by replicate and embryonic day.
faceted_tsne = function( dge, results_path ){
  plot_df = rbind( FetchData(dge, c( "orig.ident", "eday", "tSNE_1", "tSNE_2" ) ),
                   mutate( FetchData(dge, c("orig.ident", "eday", "tSNE_1", "tSNE_2" ) ), 
                           orig.ident = "all", 
                           eday = "all") )
  old_levels = levels(plot_df$orig.ident)
  levels(plot_df$orig.ident) %<>% substring( 1, 3 ) %>% tolower %>% custom.make.unique( sep = ".5 rep ")
  write.table( cbind( old_levels, levels(plot_df$orig.ident) ) , quote = F, row.names = F )
  plot_df$rep = plot_df$orig.ident %>% substr( 11, 11 ) %>% as.character
  plot_df$rep [plot_df$rep ==""] = "1"
  p = ggplot(plot_df) + ggtitle("Cells stratified by replicate") + 
    facet_grid(rep~eday) + 
    geom_point(aes(x = tSNE_1, y = tSNE_2, colour = factor(eday)), size = 0.2) + 
    scale_color_manual(values = c(scales::hue_pal()(5), "gray"))
  fp = file.path(results_path, "tsne_replicates.pdf")
  cat( paste( "Saving faceted tSNE to", fp ) )
  ggsave( filename = fp, 
          plot = p, 
          width = 20, height = 12)
  return(p)
}
