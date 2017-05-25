## ------------------------------------------------------------------------

# Several different methods of variable gene selection, all based on outliers in mean-versus-sd plots.
# You can give either x and y cutoffs for these plots or the number of genes to select -- not both.
#' @export
var_gene_select = function( dge, results_path, test_mode = F, 
                            prop_genes_to_select = NULL,
                            num_genes_to_select = NULL,
                            excess_var_cutoff = NULL,
                            log_expr_cutoff = NULL,
                            method = "seurat" ){
  atat( method %in% c( "seurat", "spline", "knn" ) )
  dir.create.nice( results_path )
  # # Parse input
  if( !is.null( num_genes_to_select ) & is.null( prop_genes_to_select ) ){
    prop_genes_to_select = num_genes_to_select / dim(dge@data)[1]
  }
  cutoffs_given = !(is.null(excess_var_cutoff) & is.null(log_expr_cutoff))
  prop_was_given  = ! ( is.null(prop_genes_to_select) || is.na(prop_genes_to_select) )
  assertthat::assert_that( xor(prop_was_given, cutoffs_given) )
  
  if(test_mode){
    prop_was_given = T
    prop_genes_to_select = 0.015
  }
  
  # # Need to do this to initialize @mean.var even if not going to use those cutoffs
  if(is.null(log_expr_cutoff)){log_expr_cutoff = 0}
  if(is.null(excess_var_cutoff)){excess_var_cutoff = 0}
  dge = Seurat::MeanVarPlot( dge, x.low.cutoff = log_expr_cutoff, y.cutoff = excess_var_cutoff)
  if(method == "seurat"){
    # Variable gene selection method 1 -- Seurat built-in
    dge@mean.var$avg_log_exp = dge@mean.var$data.x
    dge@mean.var$variability = dge@mean.var$data.y
    dge@mean.var$excess_variability = dge@mean.var$data.norm.y
  }
  if(method == "spline"){
    # Variable gene selection method 2 -- residuals from spline curve
    sd_log_exp  = unlist( apply( FUN = sd  , dge@data, MARGIN = 1 ) )
    avg_log_exp = unlist( apply( FUN = mean, dge@data, MARGIN = 1 ) )
    gam_data = data.frame( y = sd_log_exp, x = avg_log_exp )
    s = mgcv:::s
    mean_var_reg = mgcv::gam( data = gam_data, formula = y~s(x), family = mgcv::scat() )
    excess_variability = standardize( sd_log_exp - mean_var_reg$fitted.values )
    dge@mean.var$avg_log_exp = avg_log_exp
    dge@mean.var$variability = sd_log_exp
    dge@mean.var$excess_variability = excess_variability
  }
  if(method == "knn"){
    # # Variable gene selection method 3 -- knn-based outlier detection
    avg_log_exp = unlist( apply( FUN = mean, dge@data, MARGIN = 1 ) )
    sd_log_exp  = unlist( apply( FUN = sd  , dge@data, MARGIN = 1 ) )
    p = outlier_labeled_scatterplot( data.frame( avg_log_exp = avg_log_exp,
                                                 sd_log_exp = sd_log_exp, 
                                                 gene = rownames( dge@data ) ) )
    dge@mean.var$avg_log_exp = avg_log_exp
    dge@mean.var$variability = sd_log_exp
    temp = p$plot_env$data$dist_to_nn
    dge@mean.var$excess_variability = standardize( temp )
  }
  if(prop_was_given){
    dge@var.genes = top_n_preserve_rownames( dge@mean.var, 
                                             n = round(prop_genes_to_select * nrow( dge@mean.var ) ),
                                             wt = excess_variability ) %>% rownames
  } else {
    dge@var.genes = subset( dge@mean.var, 
                            excess_variability > excess_var_cutoff & 
                              avg_log_exp > log_expr_cutoff ) %>% rownames
  }
  dge@mean.var$included = rownames( dge@mean.var ) %in% dge@var.genes
  
  
  # # Save mean_var_plots 
  p_reg = ggplot(dge@mean.var) + ggtitle("Mean versus coef. var for all genes") +
    geom_point(aes(avg_log_exp, variability, colour = included), alpha = 0.5) + 
    geom_vline(aes(xintercept = log_expr_cutoff))
    
  ggsave(file.path(results_path, "MeanVarPlot_reg.pdf"), p_reg)
  
  p_res = ggplot(dge@mean.var) + ggtitle("Mean versus excess var across genes") +
    geom_point(aes(avg_log_exp, excess_variability, colour = included), alpha = 0.5) + 
    geom_vline(aes(xintercept = log_expr_cutoff)) + 
    geom_hline(aes(yintercept = excess_var_cutoff))
    
  ggsave(file.path(results_path, "MeanVarPlot_res.pdf"), p_res)
  
  # # Save variable genes and parameters
  vgsrp = file.path( results_path, "var_gene_select" )
  dir.create.nice( vgsrp )
  cell_markers = get_rene_markers()$marker %>% harmonize_species(dge)
  variable_cell_markers = intersect( cell_markers$marker, as.character( dge@var.genes ) )
  variable_cell_markers = c(paste0(length(variable_cell_markers), "total"), variable_cell_markers)
  text2file(file.path(vgsrp, "markers_among_variable_genes.txt"), variable_cell_markers)
  totalstring = paste(length(as.character(dge@var.genes)), "total")
  var_genes = c(totalstring, as.character(dge@var.genes))
  text2file(file.path(vgsrp, "variable_genes.txt"), var_genes)
  if( is.null( num_genes_to_select  ) ) { num_genes_to_select  = "NULL" }
  if( is.null( prop_genes_to_select ) ) { prop_genes_to_select = "NULL" }
  if( is.null( excess_var_cutoff    ) ) { excess_var_cutoff    = "NULL" }
  if( is.null( log_expr_cutoff      ) ) { log_expr_cutoff      = "NULL" }
  vsp        = c( prop_genes_to_select,   num_genes_to_select,   excess_var_cutoff,   log_expr_cutoff)
  names(vsp) = c("prop_genes_to_select", "num_genes_to_select", "excess_var_cutoff", "log_expr_cutoff")
  text2file(file.path(vgsrp, "var_gene_selection_params.txt"), collapse_by_name(vsp))
  
  return(dge)
}

#' @export
save_complexity_plot = function(dge, result_path){
  dir.create.nice( file.path( result_path, "QC" ) )
  f = file.path(result_path, "QC/complexity.pdf")
  p = ggplot( data.frame( complexity = colSums( dge@raw.data > 0 ) ) ) + 
    ggtitle("Complexity") + geom_histogram(aes(x = complexity)) 
  ggsave(f, p)
  f = file.path(result_path, "QC/UMIs_by_cell.pdf")
  p = ggplot( data.frame( UMIs_by_cell = colSums( dge@raw.data ) ) ) + 
    ggtitle("UMIs per cell") + geom_histogram(aes(x = log10(UMIs_by_cell)))
  ggsave(f, p)
}


## ------------------------------------------------------------------------
#' @export
add_cc_score = function(dge, method = "average"){
  # Get cell cycle genes and expression data
  # ccDat   = read.table(file.path(PATH_TO_TABLES, "ms_cellcycleGO0007049.txt"),   sep="\t")
  macosko_cell_cycle_genes = get_macosko_cc_genes()
  phases = colnames( macosko_cell_cycle_genes )
  raw_dge = as.matrix( dge@data ) # calculate scores on log1p-scale data
  
  # # To reconcile the data with the predefined list:
  # # - get human orthologs when appropriate
  # # - get both Capital and UPPERCASE
  # # - intersect gene list with available genes
  geneset = as.list(phases); names(geneset) = phases
  for( phase in phases ){
    if ( !"species" %in% AvailableData( dge ) ){
      warning(paste("No species metadata present. Assuming mouse, but please add species metadata.", 
                    "Try `object = add_maerhlab_metadata(object, 'species')` with `add_maerhlab_metadata`",
                    "from either of the packages `thymusatlasdataprivate` or `thymusatlasdatapublic`." ) )
    } else if( any( Seurat::FetchData(dge, "species")[[1]] == "human" ) ){
      human_ortho = unlist( lapply( X = macosko_cell_cycle_genes[, phase], 
                                    FUN = get_ortholog, from = "human", to = "mouse" ) )
      human_ortho = human_ortho[!is.na( human_ortho )]
    } else {
      human_ortho = c()
    }
    geneset[[phase]] = union( Capitalize( macosko_cell_cycle_genes[, phase] ), 
                     toupper(    macosko_cell_cycle_genes[, phase] ) )
    geneset[[phase]]  = union( geneset[[phase]], human_ortho )
    geneset[[phase]]  = intersect( geneset[[phase]] , rownames( raw_dge ) )
  }
  if( method == "correlation" ){
    # produce a score for each phase
    pseudocells = as.list(phases); names( pseudocells ) = phases
    for( phase in phases ){
      pseudocells[[phase]] = rownames( raw_dge ) %in% geneset[[phase]] 
    }
    scores_mat = cor( as.matrix( as.data.frame( pseudocells ) ), as.matrix( raw_dge ) )
    rownames( scores_mat ) = colnames( macosko_cell_cycle_genes )
    colnames( scores_mat ) = colnames( raw_dge )
    scores_df = data.frame(t(scores_mat))
  } else if( method == "average" ){
    # produce a score for each phase
    scores_mat = matrix( 0, nrow = length( phases ), ncol = ncol( raw_dge ) )
    rownames( scores_mat ) = colnames( macosko_cell_cycle_genes )
    colnames( scores_mat ) = colnames( raw_dge )
    for( phase in phases ){
      scores_mat[phase, ] = apply( X = raw_dge[geneset[[phase]],], MARGIN = 2, FUN = mean )
    }
    scores_df = data.frame(t(scores_mat))
  } else {
    # take the union of the cc genes and do a PCA
    method_split = down_idx(strsplit(method, split = "_"))
    method_prefix = method_split[[1]]
    method_suffix = method_split[[2]]
    atat(method_prefix == c("pca"))
    num_pc = method_suffix %>% as.numeric %>% round
    macosko_cell_cycle_genes = Reduce( union, geneset )
    macosko_cell_cycle_genes = intersect(macosko_cell_cycle_genes, rownames(raw_dge))
    dge@var.genes = macosko_cell_cycle_genes
    scores_df = Seurat::PCA(dge, pc.genes = macosko_cell_cycle_genes, pcs.store = num_pc, do.print = F)@pca.rot
    dge@var.genes = ""
    names(scores_df) = paste0("cc_pc", seq_along(scores_df))
  }  
  if( any( names( scores_df ) %in% names( dge@data.info ) ) ){
    warning("Overwriting some metadata fields in Seurat object.")
  }
  dge = Seurat::AddMetaData(dge, scores_df)
  return(list(dge = dge, cc_score_names = names(scores_df)))
}

## ------------------------------------------------------------------------
#' @export
do_dim_red = function( dge, pc.use, ... ){
  dge = Seurat::PCA( dge, do.print = F, pcs.store = max( pc.use ) )
  if( test_mode ){
    pc.use = 1:4
    # # Solve problem downstream where t-SNE can't handle exact duplicates
    if( test_mode ) {
      dge@pca.rot = dge@pca.rot + matrix( 0.1*rnorm( prod( dim( dge@pca.rot ) ) ), nrow = nrow( dge@pca.rot ) )
    }
  } 
  dge = Seurat::RunTSNE( dge, dims.use = pc.use, ... )
  return( dge )
}
 

## ------------------------------------------------------------------------
#' @export
top_genes_by_pc = function(dge, results_path, test_mode, num_pc = 30){
  if(test_mode){return()} #Doesn't work with small data. Don't care; bypassing.
  dir.create.nice(file.path(results_path, "top_genes_for_each_pc") )
  for(pc in 1:num_pc){
    pdf(file.path(results_path, "top_genes_for_each_pc", paste0("pc", pc, ".pdf")))
    VizPCA(dge, pcs.use = pc)
    dev.off()
  }
}

## ------------------------------------------------------------------------
#' A wrapper for Seurat's clustering functions.
#' 
#' @param method is either "DBSCAN" or "SNN".
#' @param granularities_as_string should be numbers separated by commas and whitespace.
#'     The number of clusterings performed is the length of (the split and cleaned version of) `granularities_as_string`.
#'     The granularity number is used as the neighborhood radius in DBSCAN or the "resolution" in SNN.
#' @details 
#' For more information on density-based clustering (DBSCAN), look at the KDD-96 paper by Ester, Kriegel, Sander, and Xu.
#' "A density-based algorithm for discovering clusters in large spatial databases with noise."
#' For more information on the SNN-based clustering, look at (the parameter gamma in)
#' "A unified approach to mapping and clustering of bibliometric networks", 
#' Ludo Waltman, Nees Jan van Eck, and Ed C.M. Noyons

# # A postcondition of this wrapper is that cluster 1 is not a legit cluster.
# # It is either empty or the set of "rejects".
#' @export
cluster_wrapper = function(dge, results_path, test_mode, 
                           pc.use = NULL, 
                           method = c("DBSCAN"),
                           granularities_as_string = "6",
                           merge_small = F){
  
  granularities = as.numeric( trimws( strsplit( granularities_as_string, split = "," )[[1]] ) )
  dir.create.nice( file.path( results_path, method ) )
  for(my_resolution in granularities){
    if(method == "DBSCAN"){
      dge = Seurat::DBClustDimension(dge, reduction.use="tsne", G.use=my_resolution, set.ident = T)
    } else {
      atat( method == "SNN" )
      dge = Seurat::FindClusters(dge, 
                                 pc.use = pc.use, 
                                 resolution = my_resolution, 
                                 print.output = F)
      ident_no_1 = dge@ident %>% as.character %>% as.numeric
      ident_no_1[ ident_no_1==1 ] = max( ident_no_1 ) + 1
      ident_no_1 = as.character( ident_no_1 )
      dge = Seurat::SetIdent(dge, ident.use = ident_no_1)
    }      
    tsne_colored( dge, file.path(results_path, method), colour = "ident",
                  fig_name = paste0("res=", my_resolution, ".pdf"))
  }
  
  

  if(test_mode & 1==length(levels(dge@ident))){
    warning("In test mode, got 1 cluster. 
            Randomly splitting into 2 clusters to facilitate debugging of differential expression code.")
    dge = Seurat::SetIdent(dge, ident.use = sample(1:2, size = length(dge@ident), replace = T) %>% as.character)
  }
  cluster_summary(dge, results_path)
  return(dge)
}

# # This records summary information for each of the clusters (currently just size). 
# # It also makes an extra copy in `named_clusters.txt` that can be hand-edited to manually name the clusters.
#' @export
cluster_summary = function(dge, results_path){
  clus_summ = data.frame(table(dge@ident))
  names(clus_summ) = c("Cluster", "Num Cells")
  write.table(x = clus_summ, 
              file = file.path(results_path, "cluster_summary.txt"),
              sep = "\t", quote = F, row.names = F, col.names = T)
  clus_summ$cluster_name = "insert_name_here"
  clus_summ$color = rainbow(length(clus_summ[[1]]))
  write.table(x = clus_summ, 
              file = file.path(results_path, "named_clusters.txt"), 
              sep = "\t", quote = F, row.names = F, col.names = T)
  return( clus_summ )
}

# # This helps compare two different analyses of the same cells.
# # You put in the usual Seurat object and results path
# # plus another Seurat object that you want to compare against,
# # and names for each of them.
#' @export
compare_views = function(dge, results_path, comparator_dge, dge_name, comparator_name){
  figname = paste0("embedding=", dge_name, "|colors=", comparator_name, ".pdf")
  dge = Seurat::AddMetaData( dge, metadata = comparator_dge@ident[dge@cell.names] , col.name = "other_id" )
  ggsave(file.path(results_path, figname),
         custom_feature_plot(dge, colour = "other_id"),
         width = 5.5, height = 5)
}

## ------------------------------------------------------------------------
#' Return spline-smoothed expression plots over pseudotime.
#'
#' @export
time_series = function( dge, gene, colour = "eday", main = NULL, x = "pseudotime", col = Thanksgiving_colors ){
  atae( length( gene ), 1 )
  if( is.null(main)){ main = paste0( "Expression by ", x)}
    
  # Sanitize input -- `aes_string` chokes on a genes with hyphens (Nkx2-1)
  rownames( dge@data ) = make.names( rownames( dge@data ) )
  rownames( dge@raw.data ) = make.names( rownames( dge@raw.data ) )
  rownames( dge@scale.data ) = make.names( rownames( dge@scale.data ) )
  gene = make.names( gene )
  
  
  my_df = FetchDataZeroPad( dge, vars.all = c( x, gene, colour ) ) 
  atat( all( sapply( my_df, FUN = is.numeric)))
  s = mgcv:::s
  p = ggplot( my_df ) + ggtitle( main ) + 
    geom_smooth( aes_string( x=x, y=gene ), colour = "black",
                 method = mgcv::gam, formula = y ~ s(x),
                 method.args = list( family = mgcv::nb() )) +
    geom_point(  aes_string( x=x, y=gene, colour = colour ) ) 
  p = p + scale_y_continuous(labels = function(x) sprintf("%4.1f", x) )
  p = p + ggtitle( gene )
  if( !is.null( col ) ){ p = p + scale_color_gradientn( colours = col ) }
  return( p )
}

#' Save plots from `times_series`.
#'
#' @export
time_series_save = function( dge, 
                             results_path, 
                             gene,
                             x = "pseudotime",
                             types = c("pdf", "pdf_no_leg", "png_pdf_split", "pdf_no_cells"), 
                             width = 8,
                             height = 6,
                             colour = "eday",
                             ... ){
  types = tolower(types)
  # Sanitize input -- `aes_string` chokes on a genes with hyphens (Nkx2-1)
  rownames( dge@data ) = make.names( rownames( dge@data ) )
  rownames( dge@raw.data ) = make.names( rownames( dge@raw.data ) )
  rownames( dge@scale.data ) = make.names( rownames( dge@scale.data ) )
  gene = make.names( gene )
  
  
  p = time_series( dge, gene, colour = colour, x = x, ... )
  results_path = file.path( results_path, "time_series" )
  dir.create.nice( results_path )
  if( "pdf" %in% types ){
    ggsave( filename = file.path( results_path, paste0( gene, ".pdf") ),
          plot = p,
          width = width, height = height)
  } 
  if( any( c("pdf_noleg", "pdf_no_leg") %in% types ) ){
    ggsave( filename = file.path( results_path, paste0( gene, "_no_leg.pdf") ),
            plot = p + theme(legend.position="none"),
            width = width, height = height)
  }
  if( any( c( "png_pdf_split", "pdf_png_split" ) %in% types ) ){
    # PNG no axis tick labels, no axis labels, and no legend
    ggsave( filename = file.path( results_path, paste0( gene, ".png") ),
            plot = p + 
              theme(legend.position="none") +
              theme(axis.text.x  = element_blank(), 
                    axis.text.y  = element_blank()) + 
              xlab("") + ylab("") + ggtitle(""),
            width = width, height = height)
    
    # ==== PDF with no points ====
    # Copy plot and remove points
    p_no_pts = p
    p_no_pts$layers = p_no_pts$layers[1]
    # Add four points to get the right y axis and color legend
    p1 = which.max( FetchDataZeroPad( dge, gene )[[1]] )[1]
    p2 = which.min( FetchDataZeroPad( dge, gene )[[1]] )[1]
    p3 = which.max( FetchDataZeroPad( dge, colour )[[1]] )[1]
    p4 = which.min( FetchDataZeroPad( dge, colour )[[1]] )[1]
    p_no_pts = p_no_pts + geom_point( data = FetchDataZeroPad( dge, c( x, colour, gene ) )[c( p1, p2, p3, p4 ) , , drop = F],
                                      aes_string( x = x, y = gene, colour = colour ) )
    ggsave( filename = file.path( results_path, paste0( gene, "_few_pts.pdf") ),
            plot = p_no_pts ,
            width = width, height = height)
  } 
  if( "pdf_no_cells" %in% types ){
    p_no_pts = p
    p_no_pts$layers = p_no_pts$layers[1]
    ggsave( filename = file.path( results_path, paste0( gene, "_no_pts.pdf") ),
            plot = p_no_pts ,
            width = width, height = height)
  }
}

#' Make tSNE plots (or PCA, or Monocle; it's customizable)
#' 
#' I needed a little bit more flexibility than Seurat was giving me with the feature plots.
#' Note: this used to plot the ranks of the data, but now it doesn't unless you specify `use_rank = T`.
#' If you want cols.use to work for categorical variables, then it should be named with the variable's levels.
#' For a blank plot (default), set `colour = NULL`.
#' @export
custom_feature_plot = function(dge, colour = NULL, subset_id = NULL, axes = c("tSNE_1", "tSNE_2"),
                               alpha = 1, cols.use = blue_gray_red, use_rank = F, overplot_adjust = F, ...){
  
  # Sanitize input -- `aes_string` was choking on a gene with a hyphen (Nkx2-1)
  rownames( dge@data ) = make.names( rownames( dge@data ) )
  rownames( dge@raw.data ) = make.names( rownames( dge@raw.data ) )
  rownames( dge@scale.data ) = make.names( rownames( dge@scale.data ) )
  colour = make.names( colour )
  axes = make.names( axes )
  my_df = FetchDataZeroPad(dge, vars.all = c(axes, colour, "ident" ), use.raw = F)
  
  # # Omit some cells if user wants to omit them
  # # but keep the plotting window the same.
  if( !is.null( subset_id ) ){
    cells.use = as.character(my_df$ident) %in% as.character(subset_id)
  } else {
    cells.use = rownames(my_df)
  }
  p = ggplot() + geom_blank( aes_string( x = axes[1], y = axes[2] ), data = my_df )
  my_df = my_df[cells.use, ]
  
  # # Treat categorical variables one way and continuous one another way.
  # # For categorical, assign randomly-ordered diverging colors if none given or not enough given
  # # Convert to hexadecimal if any given as e.g. "red"
  is_categorical = (length(colour) > 0) && ( is.factor(my_df[[colour]]) | is.character(my_df[[colour]]) )
  if( overplot_adjust & is_categorical ){
    warning("Cannot adjust for overplotting with categorical variables due to color aggregation issues. 
            Continuing with `overplot_adjust=F`." )
    overplot_adjust = F
  }
  if( is_categorical ){
    is_default = ( length( cols.use ) == length( blue_gray_red ) ) && all( cols.use == blue_gray_red )
    if( is_default || length( cols.use ) < length( unique( my_df[[colour]] ) ) ){
      better_rainbow = scales::hue_pal()
      cols.use = ( my_df[[colour]] %>% unique %>% length %>% better_rainbow )
    } else if ( any( cols.use %in% colors() ) ){
      preserved_names = names(cols.use)
      cols.use = gplots::col2hex( cols.use )
      names(cols.use) = preserved_names
    }
    p = p + scale_color_manual(values = cols.use)    + scale_fill_manual(values = cols.use)
  } else { 
    if( (length( colour ) > 0) ){
      # Optional rank transformation
      if( use_rank ){
        my_df[[colour]] = rank(my_df[[colour]]) 
        p = p + labs(colour="Cell rank")
      } else {
        p = p + labs(colour="Log normalized expression")
      }
      # Set color scale by individual points, even if aggregating as in overplot_adjust
      my_limits = c(min(my_df[[colour]]), max(my_df[[colour]]))
      p = p + 
        scale_color_gradientn(colors = cols.use, limits=my_limits ) +
        scale_fill_gradientn( colors = cols.use, limits=my_limits ) 
    }
    p = p + xlab( axes[1] ) + ylab( axes[2] ) 
  }
  
 
  if( !overplot_adjust ){
    if( length( colour ) == 0 ){
      p = p + geom_point( aes_string(x = axes[1], y = axes[2]), colour = "grey25",
                          alpha = alpha, data = my_df,
                          size = 4 / log10( length( cells.use ) ) )  
    } else {
      p = p + geom_point(aes_string(x = axes[1], y = axes[2], colour = colour), 
                         alpha = alpha, data = my_df,
                         size = 4 / log10(length(cells.use))) 
      p = p + ggtitle( colour )
    }
  } else {
    if( length( colour ) == 0 ){
      p = p + geom_hex( aes_string( x = axes[1], y = axes[2], alpha = "..count.." ), fill = "grey25",
                    data = my_df )  
      p = p + ggtitle( axes_description )
    } else {
      hex_data = hexbin::hexbin(my_df)
      hex_data = data.frame( x = hex_data@xcm, 
                             y = hex_data@ycm, 
                             count = hex_data@count )
      names(hex_data)[1:2] = axes
      nearest_bin = FNN::get.knnx( query = my_df[axes], 
                                   data = hex_data[axes], 
                                   k = 1, algorithm = "cover_tree" )$nn.index %>% c
      bin_averages = aggregate.nice( my_df[[colour]], by = nearest_bin, FUN = mean )[, 1]
      hex_data[names(bin_averages), colour] = bin_averages
      p = p + geom_point( aes_string(x = axes[1], y = axes[2], size = "count", colour = colour ),
                          data = hex_data )
      hex_data = subset( hex_data, count > 20 )
      p = p + ggtitle( colour )
    }
  }
  return(p)
}

## ------------------------------------------------------------------------
#' Save plots from `custom_feature_plot`.
#'
#' @param dge Seurat object
#' @param results_path Where to save files.
#' @param colour A gene or type of metadata. Numeric zeroes plotted if `!is.element( colour, AvailableData( dge ) )`.
#' @param fig_name Figure gets named <fig_name>.pdf or <fig_name>.png or similar. If you put a name ending in ".png" or ".pdf", the extension is stripped off.
#' @param axes Character vector of length 2. Name of numeric variables available from `FetchData`.
#' @param axes_description Character. Used in file paths, so no spaces please.
#' @param alpha Numeric of length 1 between 0 and 1. Point transparency.
#' @param height Passed to ggsave.
#' @param width Passed to ggsave, but when you ask for a legend, it gets stretched a bit to make up for lost horizontal space.
#' @param types Atomic character vector; can be longer than 1 element. If contains "PDF", you get a PDF back. If "PDF_no_leg", you get a PDF with no legend. If "PNG_PDF_split", you get back the points and bare axes in a PNG, plus text-containing elements in a PDF with no points. By default, does all three. Matching is not case sensitive.
#' @param ... Additional arguments passed to `custom_feature_plot`.
#'
#' @export
tsne_colored = function(dge, results_path, colour = NULL, fig_name = NULL,
                        axes = c("tSNE_1", "tSNE_2"), axes_description = "TSNE", 
                        alpha = 1, height = 7, width = 8, 
                        types = c("PDF", "PDF_no_leg", "PNG_PDF_split"), ... ){
  
  # Sanitize input -- `aes_string` was choking on a gene with a hyphen (Nkx2-1)
  rownames( dge@data ) = make.names( rownames( dge@data ) )
  rownames( dge@raw.data ) = make.names( rownames( dge@raw.data ) )
  rownames( dge@scale.data ) = make.names( rownames( dge@scale.data ) )
  colour = make.names( colour )
  axes == make.names( axes )
  
  # More input cleaning
  types = tolower(types)
  if( is.null( fig_name ) ){ fig_name = colour }
  fig_name %<>% strip_suffix( ".pdf" )
  fig_name %<>% strip_suffix( ".png" )
  
  # Get plot
  p = custom_feature_plot(dge = dge, colour = colour, axes = axes, alpha = alpha, ...)
  
  # Save plots
  stretch = (1 + 0.025*max( nchar(colour), 0))
  dir.create.nice( file.path( results_path, axes_description ) )
  if( "pdf" %in% types ){
    ggsave( filename = file.path( results_path, axes_description, paste0(fig_name, ".pdf") ),
          plot = p,
          width = width*stretch, height = height)
  } 
  if( any( c("pdf_noleg", "pdf_no_leg") %in% types ) ){
    ggsave( filename = file.path( results_path, axes_description, paste0(fig_name, "_no_leg.pdf") ),
            plot = p + theme(legend.position="none"),
            width = width, height = height)
  }
  if( any( c( "png_pdf_split", "pdf_png_split" ) %in% types ) ){
    # PNG no axis tick labels, no axis labels, and no legend
    ggsave( filename = file.path( results_path, axes_description, paste0(fig_name, ".png") ),
            plot = p + 
              theme(legend.position="none") +
              theme(axis.text.x  = element_blank(), 
                    axis.text.y  = element_blank()) + 
              xlab("") + ylab("") + ggtitle(""),
            width = width*stretch, height = height)
    
    # ==== PDF with no points ====
    # Copy plot and remove points
    p_no_pts = p
    p_no_pts$layers = p_no_pts$layers[1]
    # Add two points to get the right color legend 
    if(length(colour)!=0){
      max_idx = which.max( FetchDataZeroPad(dge, colour)[[1]] )[1]
      min_idx = which.min( FetchDataZeroPad(dge, colour)[[1]] )[1]
      p_no_pts = p_no_pts + geom_point( data = FetchDataZeroPad( dge, c( axes, colour ) )[c( max_idx, min_idx ) , ],
                                        aes_string( x = axes[[1]], y = axes[[2]], colour = colour ) )
    }
    ggsave( filename = file.path( results_path, axes_description, paste0(fig_name, "_no_pts.pdf") ),
            plot = p_no_pts ,
            width = width, height = height)
  }

}

## ------------------------------------------------------------------------
#' Save summary plots: plain gray, eday, clusters, replicates, nUMI, nGene, and pseudotime if available.
#'
#' @export
misc_summary_info = function(dge, results_path, clusters_with_names = NULL,
                             axes = c("tSNE_1", "tSNE_2"), axes_description = "TSNE", alpha = 1,
                             ident.use = "eday" ){
  results_path = file.path( results_path, "summaries" )
  
  # # Plot summary, checking automatically whether the colour variable is available
  fplot = function( fig_name, colour, ... ){
    if( length( colour ) == 0 || colour %in% c( names( dge@data.info ), "ident" ) ){
      tsne_colored( dge = dge, results_path,
                    fig_name = fig_name, colour = colour, 
                    axes = axes, axes_description = axes_description, alpha = alpha, ...)
    } else {
      print( paste0( "Skipping summary of ", colour, " because it's not available." ) )
    }
  }
  
  fplot( fig_name = "plain_gray.pdf", colour = NULL )
  fplot( "replicates.pdf", "rep" )
  fplot( "cell_type.pdf" , "cell_type" )
  fplot( "clusters.pdf"  , "ident" )
  fplot( "samples.pdf"   , "orig.ident" )
  fplot( "nGenes.pdf"    , "nGenes" )
  fplot( "branch.pdf"    , "branch" )
  fplot( "day.pdf"       , "eday" )
  if( all( c("pseudotime", "eday") %in% AvailableData( dge ) ) ) {
    fplot( fig_name = "pseudotime.pdf" , colour = "pseudotime" )
    ggsave( filename = file.path( results_path, "pseudotime_by_eday_box.pdf"),
            plot = ggplot( FetchData( dge, c( "pseudotime", "eday" ) ), 
                           aes( y = pseudotime, x = factor( eday ) ) ) + geom_boxplot() )
  }
}


#' Save plots en masse.
#'
#' @param dge Seurat object with available t-SNE coords (or whatever's in `axes`) 
#' @param results_path: where to save the resulting plots
#' @param top_genes: deprecated; do not use
#' @param by_cluster: deprecated; do not use
#' @param gene_list: character vector consisting of gene names
#' @param gene_list_name: used in file paths so that you can call this function again with different `gene_list_name`
#    but the same results_path and it won't overwrite.
#' @param axes: any pair of numeric variables retrievable via FetchData. Defaults to `c("tSNE_1", "tSNE_2")`.
#' @param axes_description: used in file paths so that you can call this function again with different `axes_description` but the same `results_path` and it won't overwrite.
#' @param time_series: Uses `time_series` internally instead of `custom_feature_plot`. Changes defaults for
#' `axes` and `axes_description`.
#' @param alpha Transparency of points
#' @param ... Additional parameters are passed to `custom_feature_plot` or `time_series`
#' @export
save_feature_plots = function( dge, results_path, 
                               top_genes = NULL, 
                               by_cluster = NULL,
                               gene_list = NULL, 
                               gene_list_name = NULL, 
                               axes = NULL,
                               axes_description = NULL,
                               do_time_series = F,
                               alpha = 1, ... ){
  # # Adjust defaults sensibly
  if( do_time_series ){
    if( is.null( axes            ) )  { axes             = "pseudotime" }
    if( is.null( axes_description ) ) { axes_description = "pseudotime" }
  } else {
    if( is.null( axes             ) ) { axes = c( "tSNE_1", "tSNE_2" ) }
    if( is.null( axes_description ) ) { axes_description = "TSNE" }
  }
  
  # # Defaults to rene's markers if gene_list not given
  # # If gene_list is not given, gene_list_name is replaced with "rene_picks"
  # # gene_list_name defaults to "unknown" if only gene_list_name not given
  if( is.null( gene_list ) ){
    gene_list = get_rene_markers()$marker %>% harmonize_species(dge)
    if( !is.null( gene_list_name ) ){
      warning("Overwriting gene_list_name argument with 'rene_picks' since gene_list was not given.")
    }
    gene_list_name = "rene_picks"
  } else if(is.null(gene_list_name)){
    warning("Please fill in the gene_list_name argument. Defaulting to 'unknown'.")
    gene_list_name = "unknown"
  }
  
  if(!is.null(top_genes) || !is.null(by_cluster)){
    warning( paste ( "`top_genes` and `by_cluster` arguments have been deprecated.",
                     "If you want plots of cluster markers, use the new arg `gene_list_name`." ) )
  }
  
  # # Put all feature plots in one PDF
  no_data = c()
  feature_plots_path = file.path(results_path, "feature_plots", gene_list_name)
  dir.create.nice( feature_plots_path )
  dir.create.nice( file.path( feature_plots_path ) )
   
  gene_list = as.character( gene_list )
  for( gene_name in gene_list ){
    if( !do_time_series ){
      tsne_colored( dge, results_path = feature_plots_path, colour = gene_name, 
                    axes = axes, axes_description = axes_description, alpha = alpha, ... )
    } else {
      time_series_save( dge, results_path = feature_plots_path, gene = gene_name, ... )
    }
  } 
  cat( "Plots saved to", file.path( feature_plots_path ), "\n" )
}


# # Find genes with expression patterns similar to the genes you've specified.
# #
# # `dge` : a Seurat object with field `@scale.data` filled in.
# # `markers`: a character vector; giving gene names.
# # `n`: integer; number of results to return.
# # `anticorr` : allow negatively correlated genes; defaults to `FALSE`.
# # Given a Seurat object and a list of gene names, this function returns genes 
# # that are strongly correlated with those markers. 
# # Return value: character vector.
#' @export
get_similar_genes = function( dge, markers, n, anticorr = F ){
  data.use = dge@scale.data
  if(!all(markers %in% rownames(data.use))){ 
    warning("Some of your markers have no data available. Trying Various CASE Changes.")
    markers = unique( c( markers, toupper(markers), Capitalize( markers ) ) )
  }
  markers = intersect(markers,
                      rownames( data.use) ) 
  correlation = rowSums( data.use %*% t( data.use[markers, , drop = F]) ) 
  correlation = correlation[ setdiff( names( correlation ), markers ) ]
  if( anticorr ){
    similar_genes = names( sort( abs( correlation ), decreasing = T )[ 1:n ] )
  } else {
    similar_genes = names( sort( correlation, decreasing = T )[ 1:n ] )
  }
  return( similar_genes )  
}

## ------------------------------------------------------------------------
#' @export
do_enrichr = function( results_path, geneset, geneset_name, 
                       desired_db = c( "KEGG_2016", 
                                       "WikiPathways_2016",
                                       "Reactome_2016",
                                       "BioCarta_2016",
                                       "Panther_2016",
                                       "NCI-Nature_2016", 
                                       "GO_Biological_Process_2015" ),
                       N_ANNOT_PER_DB = 2 ){
  dir.create.nice( results_path )
  
  # # Get enrichr results and parse them
  output_table = enrichR::enrichGeneList( gene.list = geneset, databases = desired_db ) 
  output_table %<>% group_by( database ) %>% top_n( wt = -pval, n = N_ANNOT_PER_DB) %>% as.data.frame
  output_table$pval = NULL
  output_table %<>% mutate( log10_qval = round( log10( qval ), 1 ) )
  output_table$qval = NULL
  output_table = output_table[c(1,2,4,3)]
  # # Save raw table
  write.table( x = output_table, 
               file = file.path( results_path, paste0("annot_", geneset_name, "_raw.txt" ) ),
               sep = "\t", quote = F, row.names = T, col.names = T )
 
  write.table( x = geneset, 
               file = file.path( results_path, paste0("annot_", geneset_name, "_annotated_genes.txt" ) ),
               sep = "\t", quote = F, row.names = F, col.names = F )
  # # Don't include genes in pretty version
  output_table$genes  = NULL

  
  # # Set up color palette and print to file
  n_colors = length(unique(output_table$database)); 
  color_idx = output_table$database %>% factor(levels = desired_db, ordered = T) %>% as.integer
  my_cols = scales::hue_pal()(n_colors)[ color_idx ]
  theme_color_db = ttheme_minimal( core=list( bg_params = list( fill = my_cols ) ) )
  ggsave( tableGrob( d = output_table, theme = theme_color_db ), 
          file = file.path( results_path, paste0("annot_", geneset_name, "_color=database.pdf")  ), 
          width = 15, height = N_ANNOT_PER_DB*length(desired_db) / 3, limitsize = F )
}

## ------------------------------------------------------------------------

#' Helper function. Check if a list of markers is compatible with a given Seurat object 
#' so that genes and cluster assignments are present in both `marker_info` and
#' the Seurat object `dge`.
#'
#' If `desired_cluster_order` is given, `are_compatible` checks that it is 
#' free of duplicates and it is a superset of the identity values occurring in other inputs.
are_compatible = function( dge, marker_info, ident.use ){
  atat( all( c("gene", "cluster") %in% names( marker_info ) ) )
  factor_flag = FALSE
  for( i in seq_along( marker_info ) ) {
    factor_flag = factor_flag || is.factor( marker_info[[i]] )
    marker_info[[i]] = as.character( marker_info[[i]] )
  }
  if( factor_flag ){ warning( "Factor columns detected in marker_info." ) }
  
  dge_ident =  as.character( FetchData( dge, ident.use )[[1]] )
  gene_compat  = all( marker_info$gene    %in% rownames( dge@data ) )
  ident_compat = all( marker_info$cluster %in% dge_ident )
  if( !gene_compat ){
    warning("marker_info$cluster has ID's not available in Seurat object")
  }
  if( !ident_compat ){
    warning("marker_info$gene has genes not available in Seurat object")
  }
  return( gene_compat && ident_compat )
}

#' Check if a list of cell types is good to be used to order genes and cells in a heatmap.
#' @details This means:
#' - no duplicates
#' - up to order, should be equal to union of table cluster labels and dge cluster labels
#' - no missing values
fix_cluster_order = function(  dge, marker_info, ident.use, desired_cluster_order = NULL ){
  if( is.null( desired_cluster_order ) ){
    warning( "No cluster order specified. Ordering clusters stupidly." )
    desired_cluster_order = union( Seurat::FetchData( dge, ident.use )[[1]], 
                                   marker_info$cluster )
  }
  atae( typeof( desired_cluster_order ), "character" )
  all_celltypes = union( Seurat::FetchData(dge, ident.use)[[1]], marker_info$cluster )
  missing = setdiff( all_celltypes, desired_cluster_order)
  extra   = setdiff( desired_cluster_order, all_celltypes)
  if( length( missing ) > 0 ){
    warning("Adding missing cell types to `desired_cluster_order` from Seurat object or marker table.")
    desired_cluster_order = c(desired_cluster_order, missing)
  }
  if( length( extra ) > 0 ){
    warning("Removing extraneous cell types from `desired_cluster_order`.")
    desired_cluster_order = intersect( desired_cluster_order, all_celltypes )
  }
  if( anyDuplicated( desired_cluster_order ) ){
    warning("Omitting duplicates in `desired_cluster_order`.")
    desired_cluster_order = unique( desired_cluster_order )
  }
  atat(  !any( is.na( desired_cluster_order ) ) )
  return( desired_cluster_order )
}


#' Returns a new ordering of the cells (a list of cell names).
#' Cells will be ordered first by cluster and then (within clusters) by the sum of 
#' the expression levels of that cluster's markers.
#'
#' @param marker_info dataframe containing variables `gene` and `cluster`.
#' @param dge should be a Seurat object.
#' @param ident.use tells you what variable to pull from the Seurat object.
#' The first three args should be compatible according to `are_compatible()`.
#' Also, `desired_cluster_order` should be a superset of `Seurat::FetchData(dge, vars.all = ident.use)[[1]]`.
optimize_cell_order = function(dge, marker_info, ident.use, desired_cluster_order ){
  # # Check inputs
  atat( are_compatible( dge, marker_info, ident.use ) )
  desired_cluster_order = fix_cluster_order( dge, marker_info, ident.use, desired_cluster_order )
  # # Set up vector: gets sorted on values, but the item of interest is the resulting ordering of the names.
  cluster_labels_orig = factor( as.character( Seurat::FetchData(dge, ident.use)[[1]]), 
                                levels = desired_cluster_order, 
                                ordered = T ) 
  names( cluster_labels_orig ) = dge@cell.names
  
  # Put clusters together
  # clcibcn means cluster_labels_contiguous_indexed_by_cell_names
  clcibcn = sort(cluster_labels_orig)

  # Rearrange each cluster by its top markers
  for(cluster_id in unique( marker_info$cluster )){
    this_cluster_idx =       ( clcibcn == cluster_id )
    this_cluster_cells = names(clcibcn)[ this_cluster_idx ]
    markers.use = marker_info$gene[marker_info$cluster == cluster_id]
    if( sum( marker_info$cluster == cluster_id ) == 0){ next } # skip reordering if no markers available
    new_cell_order = Seurat::FetchData(dge, 
                                       vars.all = markers.use, 
                                       cells.use = this_cluster_cells) %>% 
      rowSums %>% sort %>% names
    clcibcn[this_cluster_idx] = clcibcn[new_cell_order]
  }
  
  # Make sure clusters are still together
  assertthat::are_equal(sort(clcibcn), clcibcn)
  return(names(clcibcn))
}

## ------------------------------------------------------------------------

#' Set up a sparse axis with few labels for many x/y values
#' Not an ideal interface (sorry!): for eight observations in groups of 5 and 3, 
#' the `labels` input must be like c("", "", "lab1", "", "", "", "lab2", "").
#' 
sparse_axis = function(labels, side, ...){
  for( i in seq_along( labels ) ){ 
    if( "" != labels[i] )
      axis(side = side, at = i / length( labels ), labels = labels[i], ... ) 
    par(las = 2)
  }
}

#' Save a big PDF file to `<results_path>/<main>.pdf` containing a heatmap of gene expression levels.
#' 
#' @param marker_info a dataframe containing variables `gene` and `cluster`.
#' @param dge a Seurat object
#' @param results_path should be a character such that `dir.exists( results_path )`
#' @param desired_cluster_order: a superset of `levels(dge@ident)`
#' @param do_key If `TRUE` a color legend gets added
#' @param cs a vector of color names
#' @param dendrogram If TRUE, a cell dendrogram is added. It is constrained to be compatible with `dge@ident`. The 
#' implementation changes considerably. This feature is underdeveloped and other inputs -- especially `desired_cluster_order` -- may not work properly.
#' @param test_mode If `TRUE`, use only 100 genes and 100 cells.
#' @details Each column is a cell and each row is a gene. Each gene is rescaled so that its peak expression is 1.
#' This facilitates comparison within genes and across cells, though it's bad for comparison across genes.
#' @export
save_heatmap = function( dge, results_path, marker_info, 
                         desired_cluster_order = NULL,
                         main = "heatmap",
                         ident.use = "ident",
                         do_key = F,  
                         cs = blue_yellow, 
                         dendrogram = F, 
                         test_mode = F ){
  if( do_key ){warning("Sorry, do_key is not implemented right now.")}
  if( dendrogram ) atat(is.null(desired_cluster_order))
  
  # # Check inputs
  atat( are_compatible( dge, marker_info, ident.use ) )
  desired_cluster_order = fix_cluster_order( dge, marker_info, ident.use, desired_cluster_order )
  
  # reorder genes and set sparse column labels
  marker_info$cluster = factor( as.character(marker_info$cluster), 
                                   levels = desired_cluster_order, 
                                   ordered = T)
  marker_info = marker_info[order(marker_info$cluster), ]
  gene_labels = rep("", length( marker_info$cluster ))
  tmc = table( marker_info$cluster )
  tmc = tmc[tmc > 0] # Don't want to include labels such as "doublets" if there are no corresponding markers.
  non_blanks = round(cumsum(tmc) - tmc/2)+1
  gene_labels[non_blanks] = names(tmc) 

  # reorder cells and set sparse row labels
  cell_order = optimize_cell_order( dge, marker_info = marker_info, ident.use = ident.use, 
                                    desired_cluster_order = desired_cluster_order )
  cell_clusters = Seurat::FetchData(dge, ident.use)[cell_order, 1] %>% as.character
  cell_labels = rep("",length(cell_clusters))
  tcc = table(cell_clusters)[desired_cluster_order]
  tcc = tcc[tcc > 0] # Don't want to include labels such as "doublets" if there are no corresponding cells.
  non_blanks = round(cumsum(tcc)-tcc/2)+1
  cell_labels[non_blanks] = names(tcc) 
  
  # sweep out max expression level and set colorscale
  norm_expr = t( apply(X = dge@data[marker_info$gene, cell_order], FUN = div_by_max, MARGIN = 1) )

  if( test_mode ){
    rand_idx = sample(1:min(dim(norm_expr)), size = 100, replace = F) %>% sort
    cell_labels = cell_labels[rand_idx]
    gene_labels = gene_labels[rand_idx]
    norm_expr = norm_expr[rand_idx, rand_idx]
    marker_info = marker_info[rand_idx, ]
  }
  
  num_breaks = length(cs) + 1
  max_num_clust = max( length(unique(marker_info$cluster)), length(unique(dge@ident)) )
  # # diverging
  categ_colors = colorspace::rainbow_hcl(n = max_num_clust)
  # # alternating gray
  # categ_colors = colorspace::sequential_hcl(n = 2, h = 0, c = 0, l = c(20, 80))[1 + mod(0:max_num_clust, 2)] 
  # time to rock and roll!
  names(categ_colors) = desired_cluster_order

  print("Making heatmap...")
  fname = paste0( main, ".pdf" )
  pdf( file.path( results_path, fname ) )
  {
    if(!dendrogram){
      image( z = t(norm_expr), col = cs, xaxt = "n", yaxt = "n", main = main, xlab = "Cells", ylab = "Genes" ) 
      sparse_axis( labels = gene_labels, side = 2, tick = F )
      sparse_axis( labels = cell_labels, side = 1, tick = F )
      # # fields::image.plot is a good basic setup with color key, if you can ever get the damn thing working.
    } else {
      gplots::heatmap.2( norm_expr, 
                         Rowv = F, 
                         Colv = T, 
                         dendrogram = "column",
                         symm = F, 
                         scale = "none", 
                         col = blue_yellow,
                         trace = "none",
                         labCol = cell_labels, xlab = "Cells", 
                         labRow = gene_labels, ylab = "Genes",
                         ColSideColors = categ_colors[ dge@ident[ colnames( norm_expr ) ] %>% as.character ] ,
                         RowSideColors = categ_colors[ marker_info$cluster %>% as.character ] )
    }
  }
  dev.off()
  print( paste0( "Heatmap saved as ", fname) ) 
}

## ------------------------------------------------------------------------
#' Make a heatmap with one column for each cluster in `unique( Seurat::FetchData(dge, ident.use)[[1]])` and 
#' one row for every gene in `genes_in_order`. 
#' 
#' If the cluster's expression values are stored in `x`, then `aggregator(x)` gets (normalized and) plotted.
#' Optional parameter `desired_cluster_order` gets coerced to character. Should be a permutation of 
#' `unique(Seurat::FetchData(dge, ident.use))`, though elements may be omitted.
#' @export
make_heatmap_for_table = function( dge, genes_in_order, 
                                   desired_cluster_order = NULL, 
                                   ident.use = "ident",
                                   labels = NULL, 
                                   aggregator = mean, 
                                   normalize = "row", 
                                   norm_fun = div_by_max,
                                   main = "Genes aggregated by cluster" ){

  # # Set up simple ident variable
  if(is.null(desired_cluster_order)){
    warning("No cell-type ordering given. Using arbitrary ordering.")
    desired_cluster_order = list(FetchData(dge, ident.use)[1, 1])
  }
  marker_info = data.frame( gene = genes_in_order, cluster = desired_cluster_order[[1]] )
  desired_cluster_order = fix_cluster_order( dge, marker_info, ident.use, desired_cluster_order )
  ident = FetchData(dge, ident.use) %>% 
    vectorize_preserving_rownames %>% 
    factor(levels = desired_cluster_order, ordered = T)
  ident = sort(ident)
  
  # # Sanitize input -- characters for genes, and no duplicate genes.
  genes_in_order = as.character( genes_in_order )
  if( anyDuplicated( genes_in_order ) ){
    warning( "Sorry, can't handle duplicate genes. Removing them." )
    genes_in_order = genes_in_order[ !duplicated( genes_in_order )]
  }
  if( !all( genes_in_order %in% AvailableData(dge) ) ){
    warning( "Some of those markers are not available." )
    genes_in_order = intersect( genes_in_order, AvailableData(dge) )
  }
  
  # # Get cluster mean expression for each gene and row normalize
  logscale_expression = Seurat::FetchData(dge, vars.all = genes_in_order)[names( ident ), ]
  expression_by_cluster = aggregate.nice( x = logscale_expression, by = list( ident ), FUN = aggregator )
  if( normalize == "row" ){
    expression_by_cluster = apply(X = expression_by_cluster, FUN = norm_fun, MARGIN = 2 ) %>% t
  } else if( normalize == "column" ){
    expression_by_cluster = apply(X = expression_by_cluster, FUN = norm_fun, MARGIN = 1 ) 
  } else if( normalize != "none"){
    warning('normalize should be one of "row", "column", or "none". Performing row normalization.')
    normalize = "row"
  } else {
    expression_by_cluster = t(expression_by_cluster)
  }
  
  # # Form matrix in shape of heatmap and then melt into ggplot
  plot_df_wide = cbind( as.data.frame( expression_by_cluster ) , gene = rownames(expression_by_cluster))
  plot_df_wide$y = 1:nrow(plot_df_wide)
  plot_df_long = reshape2::melt( plot_df_wide, 
                                 id.vars = c("gene", "y"), 
                                 value.name = "RelLogExpr")

  plot_df_long$RelLogExpr = plot_df_long$value
  plot_df_long$Cluster = factor( as.character( plot_df_long$variable ),
                                 levels = desired_cluster_order )
  plot_df_long$gene = factor( as.character( plot_df_long$gene ),
                                 levels = genes_in_order )

  p = ggplot( plot_df_long ) + ggtitle( main ) +
    geom_tile( aes(x = Cluster, y = gene, fill = RelLogExpr ) )
  p = p + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if( is.null( labels ) ){ 
    if( length( genes_in_order ) < 30 ){
      labels = "regular"
    } else if ( length( genes_in_order ) < 80 ){
      labels = "stagger"
    } else {
      labels = "none"
    }
  }
  
  if( labels=="none"){
    p = p + theme(axis.ticks.y=element_blank(), axis.text.y = element_blank()) 
  }
  # # Add staggered labels with invisible one farther out to make room
  if( labels == "stagger" ){
    p = p + theme(axis.ticks.y=element_blank(), axis.text.y = element_blank())
    plot_df_long$stagger_pos = rep( c(-0.75, 0), length.out = nrow( plot_df_long ) )
    invisible_label_row = plot_df_long[1, ]
    invisible_label_row$gene = ""
    invisible_label_row$stagger_pos = -1.5
    invisible_label_row$y = 1
    plot_df_long = rbind(plot_df_long, invisible_label_row)
    p = p + geom_text(data = plot_df_long, aes(x = stagger_pos, y = y, label = gene ))
  }
  print(p)
  
  return( p )
}


## ------------------------------------------------------------------------
#' @export
screen_receptor_ligand = function( is_expressed, results_path ){
  
  # # Get receptor-ligand pairs; annotate with tissues expressed; save
  ramilowski = get_ramilowski()
  ramilowski$Ligand.ApprovedSymbol = NULL
  ramilowski$Receptor.ApprovedSymbol = NULL
  ramilowski = subset( ramilowski, ligand_mouse %in% rownames(is_expressed) )
  ramilowski = subset( ramilowski, receptor_mouse %in% rownames(is_expressed) )
  ramilowski$ligand_cell_types = 
    is_expressed[ramilowski$ligand_mouse, ] %>% 
    apply( 1, which ) %>% 
    sapply( names ) %>% 
    sapply( paste, collapse = "_")
  ramilowski$receptor_cell_types = 
    is_expressed[ramilowski$receptor_mouse, ] %>% 
    apply( 1, which ) %>% 
    sapply( names ) %>% 
    sapply( paste, collapse = "_")
  write.table( ramilowski, file =  file.path( results_path, "Receptor_ligand_all.txt" ), 
               sep = "\t", row.names = F, col.names = T, quote = F )

  absent = union( ramilowski$ligand_mouse, ramilowski$receptor_mouse ) %>% setdiff( rownames( is_expressed ) )
  if( length( absent ) > 0 ){
    zeropad = matrix(F, ncol = is_expressed, nrow = length( absent ), 
                     dimnames = list( gene = absent, 
                                      celltype = colnames(is_expressed)) )
    is_expressed %<>% rbind( zeropad )
  }
  
  # # Get lists of receptors and ligands for each tissue pairing
  num_unique_ligands = matrix( 0, nrow = 3, ncol = 3, 
                               dimnames = list( lig_expr_tissue = c("BLD", "MES", "TEC"),
                                                rec_expr_tissue = c("BLD", "MES", "TEC") ) )
  dir.create.nice( file.path( results_path, "ligand_lists" ) )
  dir.create.nice( file.path( results_path, "receptor_lists" ) )
  for( lig_tissue in rownames(num_unique_ligands)){
    for( rec_tissue in colnames(num_unique_ligands)){
      eligible_subset = subset( ramilowski, 
                                is_expressed[ligand_mouse,   lig_tissue] & 
                                  is_expressed[receptor_mouse, rec_tissue] )
      num_unique_ligands[lig_tissue, rec_tissue] = length( unique( eligible_subset$ligand_mouse ) )
      write.table( unique(eligible_subset$ligand_mouse  ), row.names = F, col.names = F, quote = F,
                   paste0( results_path, "/ligand_lists/"  , lig_tissue, "_to_", rec_tissue, ".txt") )
      write.table( unique(eligible_subset$receptor_mouse), row.names = F, col.names = F, quote = F,
                   paste0( results_path, "/receptor_lists/", lig_tissue, "_to_", rec_tissue, ".txt") )
      
    }  
  }
  
  # Background lists composed of everything that got successfully converted to mouse ortholog
  # For pathway analysis with a background list
  ramilowski_orig = read.table( file.path( PATH_TO_TABLES, "LigandReceptor_Ramilowski2015_mouse.txt" ), 
                                header = T, sep="\t", stringsAsFactors = F )
  write.table( ramilowski_orig$ligand_mouse   %>% unique, 
               paste0( results_path, "/ligand_lists/background.txt"),
               row.names = F, col.names = F, quote = F)
  write.table( ramilowski_orig$receptor_mouse %>% unique, 
               paste0( results_path, "/receptor_lists/background.txt"),
               row.names = F, col.names = F, quote = F)
  
  
  # # Save summary to file
  # # Sink helps get the full dimnames
  sink( file.path( results_path, "num_unique_ligands.txt" ) )
  {
    print( num_unique_ligands )
  }
  sink()
  return()
}

## ------------------------------------------------------------------------

#' Quickly explore many parameter settings
#' @export
explore_embeddings = function(dge, results_path, all_params, test_mode = F){
  atat(is.data.frame( all_params ) )
  required_params = c( "cc_method",
                       "num_pc", 
                       "clust_granularities_as_string",
                       "plot_all_var_genes" )
  atat( all( required_params %in% names( all_params ) ) )
  atat( any( c( "excess_var_cutoff","log_expr_cutoff", "prop_genes_to_select" ) %in% names( all_params ) ) )
  
  # # Record the things you're gonna try
  dir.create.nice( results_path )
  write.table( all_params, file = file.path(results_path, "params_to_try.txt"), 
               quote = F, sep = "\t", row.names = F)
  
  # # Try all the things
  prev_param_row = NULL
  for(i in rownames( all_params ) ){
    param_row = all_params[i, ]; names(param_row) = names(all_params)
    print("Trying these settings:")
    print(param_row)
    rp_mini = file.path(results_path, collapse_by_name( all_params[i,] ))
    dir.create.nice(rp_mini)

    # # remove cc variation
    if( "extra_regressout" %in% names( param_row ) ){
      extra_vars = trimws( strsplit( param_row[["extra_regressout"]], "," )[[1]] )
    } else {
      extra_vars = c()
    }
    if(is.null(prev_param_row) || param_row[["cc_method"]] != prev_param_row[["cc_method"]]){
      cc_scores_out = add_cc_score(dge, method = param_row[["cc_method"]])
      dge = Seurat::RegressOut(object = cc_scores_out$dge, 
                               latent.vars = c(cc_scores_out$cc_score_names, extra_vars))
    }
    # # select genes; do dim red; cluster cells
    dge = var_gene_select( dge, results_path = rp_mini, test_mode,
                           excess_var_cutoff   = param_row[["excess_var_cutoff"]],
                           log_expr_cutoff     = param_row[["log_expr_cutoff"]],
                           prop_genes_to_select = param_row[["prop_genes_to_select"]],
                           method = param_row[["var_gene_method"]])
    if( !is.null( param_row[[ "TF_only" ]] ) && param_row[[ "TF_only" ]] ){
      dge@var.genes %<>% intersect(get_mouse_tfs())
    }

    dge = Seurat::PCA(dge, pc.genes = dge@var.genes, do.print = F) 
    pc.use = 1:param_row[["num_pc"]]
    dge = Seurat::RunTSNE(dge, dims.use = pc.use, do.fast = T) 
    dge = cluster_wrapper(dge, results_path = rp_mini, test_mode = test_mode, 
                          method = param_row[["clust_method"]],
                          granularities_as_string = param_row[["clust_granularities_as_string"]],
                          pc.use = 1:param_row[["num_pc"]])
                             
    # # save plots and summaries 
    misc_summary_info( dge, results_path = rp_mini)
    saveRDS( dge, file.path( rp_mini, "dge.data") ) 
    
    top_genes_by_pc(   dge, results_path = rp_mini, test_mode)
    fm = param_row[["find_markers_thresh"]]
    if( !is.null( fm ) ){
      de_genes = find_de(dge, results_path = rp_mini, ident.use = "ident", thresh.use = fm )
      top_markers = get_top_de_genes(de_genes, results_path = rp_mini, 
                                     test_mode = test_mode, 
                                     top_n = 10, thresh_df = NULL)
      save_feature_plots(dge, results_path = rp_mini, gene_list = top_markers$gene, gene_list_name = "top_markers")
      # save_heatmap(dge, results_path = rp_mini, marker_info = de_genes,    main = "heatmap_all_de_genes")
      save_heatmap(dge, results_path = rp_mini, marker_info = top_markers, main = "heatmap_top_markers")
    }
    save_feature_plots(dge, results_path = rp_mini)
    pavg = param_row[["plot_all_var_genes"]]
    if( !is.null( pavg ) && pavg  ){
      save_feature_plots(dge, results_path = rp_mini, gene_list = dge@var.genes, gene_list_name = "var_genes")
    }
    prev_param_row = param_row
  }
  return(dge)
}


