## ------------------------------------------------------------------------

#' Characterize genes by behavior over pseudotime, returning cluster assignments and p values. 
#' 
#' @param dge should be a seurat object with a field "pseudotime". 
#' The field `dge@data` is accessed for expression levels -- for Eric's objects, the units will be log2(1+CP10K).
#' @param results_path is a character vector showing where to dump the output.
#' @param num_periods_initial_screen Cells are partitioned into this many pseudotime periods (equal number of cells in each). Initial screening is based on a piecewise linear model where expression is constant within these bins.
#' @param pval_cutoff Genes are screened by p-value to avoid too much computationally expensive smoothing.
#' @param do_adjust Logical -- if TRUE, then apply BH correction. 
#' @param abcissae_kmeans Gene expression is fed into k-means as a series of predictions at successive time points.
#' The arguments says how many time points to predict and feed in (if length one) or what time points (if longer).
#' @param num_clusters Genes are partitioned into this many modules. If NULL (default) the value is selected via the gap statistics and their SEs using the method in the original gap statistic paper:
#' 
#' Tibshirani, R., Walther, G. and Hastie, T. (2001). Estimating the number 
#' of data clusters via the Gap statistic. Journal of the Royal Statistical Society B, 63, 411–423.
#'                                 
#' There's one adjustment: this function will never use just one cluster. It will issue a warning and use 2.
#'
#' @return A list with elements:
#' \itemize{
#' \item dge: the Seurat object
#' \item gene_stats: genes with effect sizes and cluster labels.
#' \item smoothers: fitted regression models, one for each gene.
#' \item cluster_mod: output from stats::kmeans 
#' \item gap_stats: output from cluster::clusGap
#' }
#'
#' 
#' This function helps explore gene dynamics over pseudotime. It goes through three main steps:
#' \itemize{
#' \item find genes that respond strongly to pseudotime.
#' \item smooth those genes' expression to form an overall pseudotime trend.
#' \item cluster genes based on smoothed expression patterns that have been shifted/scaled to the unit interval.
#' }
#'
#' @export
#'
smooth_and_cluster_genes = function( dge, results_path, genes.use = get_mouse_tfs(),
                                     num_periods_initial_screen = 20, 
                                     pval_cutoff = 0.05, do_adjust = T, 
                                     abcissae_kmeans = 20,
                                     num_clusters = NULL ){
  
  dir.create.nice( results_path )
  
  # # Filter genes based on qvals from piecewise linear model 
  cell_data = Seurat::FetchData( dge , "pseudotime" )
  nbreaks = num_periods_initial_screen+1
  breaks = quantile(cell_data$pseudotime, probs = (1:nbreaks - 0.5)/nbreaks)
  cell_data$period = cut( x = cell_data$pseudotime, 
                          breaks = breaks, 
                          labels = as.character( 1:num_periods_initial_screen ),
                          include.lowest = T, ordered_result = T)
  dge = Seurat::AddMetaData( dge, cell_data[, "period"] %>% matrixify_preserving_rownames, "period" )
  ssq = function(x) sum((x - mean(x))^2)
  genes.use = FetchDataZeroPad( dge, genes.use ) %>% apply(2, prop_nz) %>% is_greater_than(0.05) %>% which %>% names
  expr_mat  = FetchDataZeroPad( dge, genes.use ) %>% as.matrix
  n = length(dge@cell.names)
  df1 = n - 1
  df2 = n - num_periods_initial_screen
  mss_constant  = apply( expr_mat, 2, ssq ) %>% divide_by(df1)
  mss_piecewise = aggregate_nice( expr_mat, by = cell_data$period, FUN = ssq  ) %>% t %>% rowSums %>% divide_by(df2)
  get_p = function(x) pf( x, df1 = df1, df2 = df2, lower.tail = F )
  pvals = (mss_constant  / mss_piecewise) %>% sapply( get_p )
  names(pvals) = genes.use
  if( do_adjust ){
    pvals %<>% p.adjust()
  }
  genes_included = names(which(pvals < pval_cutoff ))
  
  # Smooth gene expression over pseudotime
  pt_data = FetchData( dge, c( "pseudotime", genes_included ) )
  s = mgcv:::s
  fast_smooth = function( i ) { 
    gene = genes_included[i]
    expr = pt_data[[gene]] 
    pseudotime = pt_data[["pseudotime"]]
    if( i %% 50 == 0 ){
      cat("Smoothing gene ", i, " of ", length( genes_included ), "\n" )
    }
    return(  mgcv::gam( expr ~ s( pseudotime ), family = mgcv::nb() ) )
  }
  smoothers = lapply( seq_along( genes_included ), fast_smooth )
  names( smoothers ) = genes_included

  # # Generate kmeans features
  if( length( abcissae_kmeans ) == 0 ) { abcissae_kmeans = 20 }
  if( length( abcissae_kmeans ) == 1 ) {
    abcissae_kmeans = quantile(cell_data$pseudotime, (1:abcissae_kmeans - 0.5)/abcissae_kmeans )
  }
  kmeans_features = sapply( smoothers, predict, type = "response", 
                            newdata = data.frame( pseudotime = abcissae_kmeans ) )
  kmeans_features = t(kmeans_features)
  
  # # Rescale each gene to the unit interval
  kmeans_features = apply( kmeans_features, MARGIN = 1, FUN = function( x ) (x - min(x)) ) %>% t
  kmeans_features = apply( kmeans_features, MARGIN = 1, FUN = div_by_max ) %>% t
  atae(ncol(kmeans_features), length(abcissae_kmeans))
  atae(nrow(kmeans_features), length(genes_included))
  rownames(kmeans_features) = genes_included
  
  cat("  Calculating gap statistics... \n")
  gap_stats = cluster::clusGap( x = kmeans_features, FUNcluster = kmeans, K.max = 25, spaceH0 = "scaledPCA", 
                                iter.max = length(abcissae_kmeans) * 2 )
  {
    pdf( file.path( results_path, "gap_stats.pdf" ) )
    plot( gap_stats )
    dev.off()
  }
  if( is.null( num_clusters ) ) {
    num_clusters = cluster::maxSE(gap_stats$Tab[, 3], gap_stats$Tab[, 4], "Tibs2001SEmax") 
    if( num_clusters == 1 ){
      warning("Based on gap statistics, Tibshirani et alii suggest one cluster, but that's boring so I'll use two.")
      num_clusters = 2
    }
  }
  cat("  Continuing with ", num_clusters,  " clusters.\n")
  cluster_mod = kmeans( kmeans_features, centers = num_clusters, iter.max = 500 )
  
  # Reorder clusters by hierarchical clustering 
  peak_position = apply( cluster_mod$centers, 1, function( x ) mean( which( x > 0.75 ) ) )
  preferred_ordering = order( peak_position )

  old_given_new = preferred_ordering
  new_given_old = function( k ){ which(k==old_given_new)}
  cluster_mod$cluster  = sapply( cluster_mod$cluster, new_given_old )
  cluster_mod$centers  = cluster_mod$centers [old_given_new, ]
  cluster_mod$withinss = cluster_mod$withinss[old_given_new]
  cluster_mod$size     = cluster_mod$size    [old_given_new]
  rownames(cluster_mod$centers) = NULL
  
  # # Assemble data for export
  cat("Exporting cluster assignments and pvals for genes under study... \n")
  to_return = data.frame( gene = genes_included,
                          pval = pvals[genes_included], 
                          cluster = cluster_mod$cluster[genes_included] )
  to_return = to_return[order( to_return$cluster, to_return$pval ), ]
  write.table( to_return, file.path(results_path, "gene_pt_dependence_stats.txt"), 
               sep = "\t", quote = F, row.names = F, col.names = T)
  
  cat("Done.\n")
  return( list( dge = dge,
                # smoothing/filtering
                gene_stats = to_return,
                smoothers = smoothers, 
                # clustering
                abcissae_kmeans = abcissae_kmeans,
                kmeans_features = kmeans_features,
                cluster_mod = cluster_mod,
                gap_stats = gap_stats ) )
}

#' Draw overlaid line plots of gene clusters from output of `smooth_and_cluster_genes`.
#'
#' @param facet_ncol For resulting facet plots, number of columns.
#' @export
facet_plot_gene_clusters = function( dge, results_path, 
                                     cluster_mod, 
                                     kmeans_features,
                                     abcissae_kmeans, 
                                     facet_ncol = NULL ){
  cluster_mod$cluster = factor( cluster_mod$cluster )
  num_clusters = length(levels(cluster_mod$cluster))
  centers = aggregate_nice( kmeans_features, by = cluster_mod$cluster, FUN = mean )
  data_wide = data.frame( rbind( centers, kmeans_features ) )
  data_wide$is_center = F; data_wide$is_center[1:num_clusters] = T
  data_wide$cluster = c( levels(cluster_mod$cluster), as.character( cluster_mod$cluster ) )
  data_wide$gene = rownames( data_wide )
  
  data_long = reshape2::melt( data_wide, id.vars = c( "is_center", "cluster", "gene" ))
  data_long$unit_scaled_expression = data_long$value
  data_long$pseudotime = abcissae_kmeans[ data_long$variable ]
  data_long$cluster %<>% factor(levels = rev(sort(unique(data_long$cluster))), ordered = T)
  p_faceted_clusters = ggplot() + ggtitle( "Major gene clusters" ) +
    geom_line( data = subset( data_long, !is_center), alpha = 0.6,
               aes( x = pseudotime, y = unit_scaled_expression, group = gene, colour = cluster ) ) +
    geom_line( data = subset( data_long, is_center), alpha = 1, colour = "black",
               aes( x = pseudotime, y = unit_scaled_expression, group = gene  ) ) +
    facet_wrap( ~cluster, ncol = facet_ncol ) + 
    theme(strip.background = element_blank(), legend.position = "none") +
    ylab("Unit-scaled expression") + 
    scale_colour_grey()

  ggsave( filename = file.path(results_path, "faceted_gene_clusters.pdf"), 
          p_faceted_clusters, 
          width = 5, height = 6 )
  return( p_faceted_clusters )
}

#' Draw heatmaps of gene clusters from output of `smooth_and_cluster_genes`.
#'
#' @param dge Seurat object with raw data and pseudotime metadata used to train smoothers. 
#' If metadata `simple_branch` is present, function is hardwired to look for "mTEC", "cTEC", "branchpoint", and "progenitor"
#' and make a heatmap similar to figure 2 in http://dx.doi.org/10.1101/122531. 
#' @param results_path Where to save plots and files.
#' @param cluster_mod K-means output with cluster labels for the genes in `smoothers` and also with cluster centers.
#' @param smoothers List of regression models for the genes in `cluster_mod` and `gene_stats`.
#' @param gap_size White bars separating clusters are formed by adding fake genes. gap_size is how many fake genes per bar.
#' @param genes_use Gene names or anything in `AvailableData(dge)`.
#' @param genes_to_label Subset of `genes_use` to write tick labels for.
#' @export
heatmap_gene_clusters = function( dge, results_path, 
                                  cluster_mod, 
                                  smoothers,
                                  gap_size = NULL,
                                  genes_use = NULL,
                                  genes_to_label = NULL ){  
  cell_data = dge %>% FetchData("pseudotime")
  num_clusters = length( unique( cluster_mod$cluster ) )
  atae( names( smoothers ), names(cluster_mod$cluster) )
  
  if( is.null( genes_use ) ){
    genes_use = names(cluster_mod$cluster)
  } else {
    cluster_mod$cluster = cluster_mod$cluster[genes_use]
    smoothers = smoothers[genes_use]
  }
  
  # # set up wide-format data
  if( "simple_branch" %in% AvailableData( dge ) ){
    abcissae_heatmap = FetchData(dge, "pseudotime")[[1]] 
    #avoids exact duplicates; need to use these as dimnames
    abcissae_heatmap = abcissae_heatmap + rnorm( n = length(abcissae_heatmap), mean = 0, sd = 1e-8 )
    abcissae_heatmap %<>% sort
  } else {
    abcissae_heatmap = seq( min(cell_data$pseudotime), max(cell_data$pseudotime), length.out = 100 )
  }
  data_wide_heat = sapply( smoothers, predict, type = "response", 
                           newdata = data.frame( pseudotime = abcissae_heatmap ) ) 
  rownames(data_wide_heat) = abcissae_heatmap
  get_expression_peak = function(x){
    weighted.mean( x = seq_along( x ), w = x * ( standardize(x)>0 ) )
  }
  peak_expression = apply( data_wide_heat, 2, get_expression_peak )
  names(peak_expression) = names(smoothers)
  
  data_wide_heat = apply(data_wide_heat, 2, standardize) %>% t %>% as.data.frame
  data_wide_heat$gene = names(smoothers)
  data_wide_heat$cluster = cluster_mod$cluster

  # # Add bars between clusters by adding dummy genes ranked last in each cluster
  if( is.null(gap_size) ){
    gap_size = ifelse( length(genes_use) < 50, 0, ceiling( length(smoothers) / 100 ) )
  }
  scaffold = rep( 1:nrow(cluster_mod$centers), each = gap_size )
  if( gap_size > 0 ){
    dummy_genes = matrix( NA, 
                          ncol = length( abcissae_heatmap ), 
                          nrow = length( scaffold ) ) %>% as.data.frame
    dummy_genes$gene = paste0("DUMMY_GENE_", scaffold ) %>% make.unique
    dummy_genes$cluster =                    scaffold
    colnames(dummy_genes) = colnames(data_wide_heat)
    data_wide_heat = rbind( data_wide_heat, dummy_genes )
    peak_expression = c( peak_expression, rep( Inf, length(scaffold) ) ) #used for ranking
  }
  
  # # Add vertical bar for branching heatmap
  if( "simple_branch" %in% AvailableData( dge ) ){
    pt_by_branch = aggregate_nice( x =  FetchData(dge, "pseudotime"), 
                                   by = FetchData(dge, "simple_branch"), 
                                   FUN = mean )
    atae( pt_by_branch["progenitor", ], 0, tol = 1e-8 )
    atae( pt_by_branch["branchpoint", ], 0, tol = 1e-8 )
    atat( pt_by_branch["mTEC", ] < 0 )
    atat( pt_by_branch["cTEC", ] > 0 )
    branch_counts = FetchData(dge, "simple_branch") %>% table
    boundary = branch_counts["mTEC"] + branch_counts["progenitor"] / 2 + branch_counts["branchpoint"] / 2
    dummy_cells = matrix(NA, nrow = nrow( data_wide_heat ), ncol = ncol( data_wide_heat )/40 )
    middle_pt = abcissae_heatmap[boundary + 0:1 ]
    colnames(dummy_cells) = seq( max(middle_pt) + 1e-8, 
                                 min(middle_pt) - 1e-8, length.out = ncol( dummy_cells ) )
    is_left = 1:ncol(data_wide_heat) < boundary
    data_wide_heat = cbind( data_wide_heat[,  is_left],
                            dummy_cells,
                            data_wide_heat[, !is_left] )
  }
  
  # # Melt data and order rows
  data_wide_heat$gene = factor( data_wide_heat$gene, 
                                levels = data_wide_heat$gene[order(data_wide_heat$cluster, peak_expression)] , 
                                ordered = T)
  data_wide_heat = data_wide_heat[order(data_wide_heat$gene), ]
  data_long_heat = reshape2::melt(data_wide_heat, id.vars = c("gene", "cluster"))
  data_long_heat %<>% plyr::rename( replace = c("variable" = "pseudotime", "value" = "rel_log_expr") )

  # # Smoothed heatmap
  p_heat_clusters = ggplot(data = data_long_heat ) + 
    ggtitle( "Genes clustered by temporal expression pattern" ) +
    geom_raster( aes( x = pseudotime, y = gene, fill = rel_log_expr ), 
                 interpolate = T ) +
    scale_fill_gradientn( colors = blue_gray_red, na.value="white" ) + 
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x=element_blank()) 
  
  # # Label selected genes
  rest = setdiff( data_wide_heat$gene, genes_to_label )
  ylabels = c( genes_to_label, rep("", length(rest))  )
  names(ylabels) = c( genes_to_label, rest )
  p_heat_clusters = p_heat_clusters + 
    scale_y_discrete( labels = ylabels ) + labs( y = "Genes" ) 

  # # X axis is either by cells (branching) or by pseudotime (not branching).
  if( "simple_branch" %in% AvailableData( dge ) ){
    p_heat_clusters = p_heat_clusters + xlab("Cells")  
  } else {
    p_heat_clusters = p_heat_clusters + xlab("Pseudotime")
  }
  
  # # Add cluster colorbar
  boundary_genes = which(1==diff(sort(data_wide_heat$cluster)))
  boundary_genes = c(0, boundary_genes, length( data_wide_heat$gene ) ) 
  right_edge = max(as.numeric(data_long_heat$pseudotime))
  left_edge  = min(as.numeric(data_long_heat$pseudotime))
  for( i in 1:num_clusters ){
    p_heat_clusters = p_heat_clusters +
      annotate( "rect", 
                xmin = left_edge - 4 * (right_edge - left_edge) / 100,
                xmax = left_edge -     (right_edge - left_edge) / 100,
                ymin = boundary_genes[i  ] + gap_size,
                ymax = boundary_genes[i+1],
                fill = scales::hue_pal()(num_clusters)[i]) + 
      annotate( "text", 
                x = left_edge - 8 * (right_edge - left_edge) / 100, 
                y = (boundary_genes[i  ] + gap_size + boundary_genes[i+1]) / 2, 
                label = i) 
    
  }
  
  # # Add cell types
  branch_colors = c(scales::hue_pal()(2), "purple", "gray" )
  names(branch_colors) = c("mTEC", "cTEC", "progenitor", "branchpoint" )
  if( "simple_branch" %in% AvailableData( dge ) ){
    X = FetchData(dge, c("simple_branch", "pseudotime")) 
    X = X[order(X$pseudotime), ]
    X = X[["simple_branch"]] %>% as.character %>% rle #run-length encoding keeps from overwhelming ggplot
    X$positions = c(0, cumsum(X$lengths) )
    X$positions %<>% div_by_max
    X$positions = X$positions*right_edge
    for( i in 1:length(X$values) ){
      cell_label_data = 
        data.frame( 
          xmin = X$positions[i], 
          xmax = X$positions[i+1], 
          ymin = 2*ceiling( length(smoothers) / 100 ),
          ymax = 0 )
      p_heat_clusters = p_heat_clusters + 
        geom_rect( data = cell_label_data,
                   aes( xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax ), 
                   fill = branch_colors[X$values[i]] ) 
      
      dummy_cells_pseudotime = as.numeric(colnames(dummy_cells)) * right_edge / max(X$positions)
      p_heat_clusters = p_heat_clusters + 
        annotate( geom = "rect", fill = "white", 
                  xmin = min( dummy_cells_pseudotime ), 
                  xmax = max( dummy_cells_pseudotime ), 
                  ymin = 2*ceiling( length(smoothers) / 100 ), 
                  ymax = 0 ) 
    }                
  }

  ggsave( filename = file.path(results_path, "heatmapped_gene_clusters.png"), 
          p_heat_clusters, width = 6, height = 7 )
  ggsave( filename = file.path(results_path, "heatmapped_gene_clusters.pdf"), 
          p_heat_clusters, width = 6, height = 7 )
  return( p_heat_clusters )
}

#' Plot a dendrogram, labeling merges given by `converter` with colored rectangles. 
#'
#' @details If converter is c("a", "a", "b", "c"), you get edges labeled as a 
#' four-colour muted rainbow next to a grayscale with light, light, medium, and dark.
#' This is to convey that the first two tips have been somehow grouped.
#' @export
 plot_dendro_with_rect = function( hc, converter, main = "Dendrogram" ){
    p = ggdendro::ggdendrogram(hc)
    p = p + annotate( geom = "tile", x = 1:length(hc$order), 
                      y = -2, 
                      fill = scales::hue_pal()(length(hc$order)) ) 
    p = p + annotate( geom = "tile", x = 1:length(hc$order), 
                      y = -1, 
                      fill = scales::grey_pal()(3)[factor(converter[hc$order])] ) 
    p = p + ggtitle(main)
    print(p)
    return(p)
  }
#' Plot a dendrogram and also cut it to merge input into `num_desired` groups.
#'
#' @export
dendrogram_merge_points = function( X, num_desired, results_path, 
                                    FUN = function(x) hclust(dist(x), method = "ward.D2"), 
                                    REORDER_FUN = function(hc) as.hclust(stats::reorder(as.dendrogram(hc), 1:nrow(X))), 
                                    PLOT_FUN = plot_dendro_with_rect,
                                    CUT_FUN = stats::cutree,
                                    main = "Dendrogram",
                                    return_hc = F,
                                    ... ){
  hc = FUN( X, ... )
  hc = REORDER_FUN( hc )
  converter = CUT_FUN(hc, num_desired)
  converter = setNames( letters[converter], names(converter) )
  atae(as.character(hc$labels), names(converter))
  hc$labels = paste0( hc$labels, " (", converter[hc$labels], ")" ) 
  {
    pdf( file.path( results_path, paste0( main, ".pdf" ) ) )
    PLOT_FUN( hc, converter, main = main )
    dev.off()
  }
  PLOT_FUN( hc, converter, main = main )
  
  if(return_hc){
    return( hc )
  }
  return( converter )
}


