## ------------------------------------------------------------------------

#' Run Monocle, DPT, or PCA-based pseudotime analysis on a Seurat object.
#'
#' @param dge Seurat object
#' @param results_path Where to save tables and plots
#' @param earliest_day Numeric. Used to root the Monocle trajectory; gets compared to `FetchData(dge, "eday")[[1]]` via `==`.
#' @param mp Passed to call_monocle_on_seurat as monocle_params arg.
#' @param ... Additional parameters passed to call_monocle_on_seurat.
#'
#' The statistical inference is outsourced to Monocle via call_monocle_on_seurat (or Seurat via pc_as_pt).
#' @export
master_pt = function( dge, results_path, method = "monocle", earliest_day = NULL,
                      reset_var_genes = T,
                      mp = list( reset_var_genes       = reset_var_genes, 
                                 log_scale_expr_thresh = 0.1,
                                 excess_disp           = 1,
                                 num_mature_types      = NULL,
                                 reduction_method      = "DDRTree" ), ...  ){
  # # Run monocle & transfer results
  if( method=="PCA" ){
    dge = pc_as_pt( dge )$dge
  } else if (method=="DPT") {
    if( !requireNamespace("destiny", quietly = TRUE) | packageVersion("destiny") < "2" ){
      stop("Option DPT requires the package 'destiny' version >=2.5.0.")
    }
    if( reset_var_genes ){
      dge %<>% MeanVarPlot
    }
    #dm = destiny::DiffusionMap( data = t( as.matrix( dge@scale.data ) ) )
    dm = destiny::DiffusionMap( data = t( as.matrix( dge@scale.data[dge@var.genes, ] ) ) )
    dpt_out = destiny::DPT( dm, ... )
    dge = add_pseudotime_to_seurat( dge, pt_obj = dpt_out, pt_method = "dpt" )
  } else {
    mobj = call_monocle_on_seurat( dge, results_path, monocle_params = mp, earliest_day = earliest_day, ... )
    dge = add_pseudotime_to_seurat( dge, pt_obj = mobj )
  }
  return( dge )
}


#' Take a principal component (first by default) and return a list of results gained by interpreting that PC as a pseudotime axis.
#'
#' @param dge: Seurat object to be used. Should have variable genes and PCA fields already filled in. 
#' Jackstraw and t-SNE are optional but preferred.
#' @param orient_var: pseudotime is flipped to correlate positively with this. Anything available 
#' from `Seurat::FetchData(dge, orient_var)` is fair game... as long as it's numeric.
#' @param pc.use: which principal component to take
#' 
#' @return List with named elements:
#' `$dge`: Seurat object with a `pseudotime` metadatum filled in
#' `$gene_corrs`: dataframe with genes ordered by correlation with pseudotime.
#'    Columns are `gene`, `corr`, and (if available from Jackstraw) `p_value`.
pc_as_pt = function( dge, pc.use = 1, orient_var = "eday" ){
  
  # # Construct pseudotime
  pt = Seurat::FetchData(dge, paste0("PC", pc.use)) %>% vectorize_preserving_rownames
  
  # # Get correlations with genes and sort through them
  dge = ProjectPCA( dge, do.print = F )
  gene_corrs = dge@pca.x.full[, "PC1", drop = F] %>% vectorize_preserving_rownames %>% sort
  gene_corrs = data.frame( gene = names(gene_corrs), 
                           corr = gene_corrs ) 
  if( prod( dim( dge@jackStraw.empP.full ) ) > 0 ){
    # Need to index by gene because of earlier sorting
    gene_corrs$p_value = dge@jackStraw.empP.full[gene_corrs$gene, "PC1"]
  } 
  rownames(gene_corrs) = gene_corrs$gene
  
  # # Flip everything if necessary
  ov = Seurat::FetchData( dge, orient_var )[[1]]
  atat(is.numeric(ov))
  wrong_way = cor( pt, ov ) < 0
  if( wrong_way ){ orient = function(x) -x } else { orient = function(x) x }
  pt                  %<>% orient
  dge@pca.rot   [, 1] %<>% orient
  dge@pca.x     [, 1] %<>% orient
  dge@pca.x.full[, 1] %<>% orient
  gene_corrs$corr %<>% orient
  
  # # Add pt as metadata and return
  dge = Seurat::AddMetaData( dge, pt, "pseudotime" )
  return( list( dge = dge, gene_corrs = gene_corrs ) )
}

## ------------------------------------------------------------------------
#' Extract data from a Seurat object and run Monocle, returning a Monocle object.
#'
#'
#' @export
call_monocle_on_seurat = function( dge, results_path, monocle_params, earliest_day = NULL ){
  
  attach( monocle_params )
  
  # # Make (m)onocle (obj)ect
  raw_dge = deseuratify_raw_data( dge )
  rownames(dge@mean.var) %<>% make.names
  geneInfo = data.frame( dge@mean.var[rownames(raw_dge), ] )
  geneInfo$gene = rownames(geneInfo)
  adf_class = getClassDef("AnnotatedDataFrame", where, package = "Biobase")
  pd = new(adf_class, data = dge@data.info[colnames(raw_dge), ])
  fd = new(adf_class, data = geneInfo)
  mobj = monocle::newCellDataSet( as.matrix( raw_dge ), 
                                  phenoData = pd, 
                                  featureData=fd,
                                  lowerDetectionLimit=1,
                                  expressionFamily=VGAM::negbinomial.size( ) )
  
  # # Select variable genes or get from Seurat object
  mobj = monocle::detectGenes( mobj, min_expr = 0 )
  mobj = BiocGenerics::estimateSizeFactors( mobj )
  mobj = BiocGenerics::estimateDispersions( mobj )
  if(reset_var_genes){
    disp_table = monocle::dispersionTable( mobj )
    ordering_genes = subset(disp_table,
                            ( mean_expression >= log_scale_expr_thresh ) & 
                              ( dispersion_empirical >= excess_disp * dispersion_fit ) )[["gene_id"]]
    cc_genes_macosko = get_macosko_cc_genes()
    cc_genes_all = Reduce(union, cc_genes_macosko) %>% sapply( get_ortholog ) %>% unique %>% na.omit
    ordering_genes %<>% as.character %>% setdiff( cc_genes_all )
    expressed_genes = row.names(subset(Biobase::fData(mobj), num_cells_expressed >= 10))
    ordering_genes = intersect(expressed_genes,ordering_genes)
    mobj = monocle::setOrderingFilter(mobj,ordering_genes)
  } else {
    mobj = monocle::setOrderingFilter(mobj, dge@var.genes )
  }
  
  # # Save variable genes and parameters
  vgsrp = file.path( results_path, "var_gene_select" )
  dir.create.nice( vgsrp )
  gd = mobj@featureData@data
  var_genes = gd$gene[gd$use_for_ordering]
  cell_markers = get_rene_markers()
  variable_cell_markers = intersect( Capitalize(cell_markers$marker), Capitalize(as.character(var_genes)) )
  variable_cell_markers = c(paste0(length(variable_cell_markers), "total"), variable_cell_markers)
  text2file( file.path( vgsrp, "markers_among_variable_genes_monocle.txt" ), variable_cell_markers )
  totalstring = paste(length(as.character(var_genes)), "total")
  var_genes = c(totalstring, as.character(var_genes))
  text2file(file.path(vgsrp, "variable_genes_monocle.txt"), var_genes)
  vsp        = c( excess_disp,   log_scale_expr_thresh)
  names(vsp) = c("excess_disp", "log_scale_expr_thresh")
  text2file(file.path(vgsrp, "var_gene_selection_params_monocle.txt"), collapse_by_name(vsp))
  
  
  # # Do Monocle dimension reduction, correcting for cell cycle but not batch effects 
  # # (For us, batch effects are often completely nested inside of embryonic day, 
  # # so we can't adjust for them without removing the temporal signal.)
  form_str = paste0( "~", paste0( CC_PHASES, collapse = " + " ) )
  library(DDRTree) # required to prevent  "error in DDRTree_reduce_dim: cannot convert to environment"
  mobj = monocle::reduceDimension(mobj,residualModelFormulaStr=form_str)
  mobj = monocle::orderCells(mobj, num_paths = num_mature_types)
  
  if ( is.null( earliest_day ) ) {
    earliest_day = Biobase::phenoData(mobj)[["eday"]] %>% as.character %>% as.numeric %>% min
  } 
  
  # # This code roots the lineage tree using a horrible horrible hack!  
  # # Monocle crashes when you give orderCells a root_state that is an internal node of the ddrtree graph.
  # # I can't tell which f***ing nodes are internal. So I pick the one I want based on embryonic day,
  # # then rank the states by
  # # their proximity to it in one of the 2d embeddings. I try rooting the tree at the best
  # # ones first, stopping when something works.
  pct_eq_10_5 = function(x) mean( x == earliest_day )
  eday  = Biobase::phenoData(mobj)[["eday"]] %>% as.character %>% as.numeric
  state = Biobase::phenoData(mobj)[["State"]] %>% as.character %>% as.numeric
  coords = t( mobj@reducedDimS )
  eday_by_state     = aggregate_nice( x = eday,   by = state, FUN = pct_eq_10_5 )
  centroid_by_state = aggregate_nice( x = coords, by = state, FUN = mean )
  root_state = which.max( eday_by_state ) 
  distance_sq = function(x, y) sum((x - y)*(x - y))
  dist_sq_to_root = apply( centroid_by_state, 1, FUN = function(x) distance_sq(x, centroid_by_state[root_state,]) )
  names( dist_sq_to_root ) = rownames( centroid_by_state )
  dist_sq_to_root = sort( dist_sq_to_root, decreasing = F )
  for( i in seq_along( dist_sq_to_root ) ){
    candidate_state = names( dist_sq_to_root )[[i]]
    e = tryCatch( expr = { mobj = orderCells( mobj, root_state = candidate_state ); e = "success" },
                  error=function( err ){ "failure" } )
    if( e != "failure" ){
      break
    }
  }
  
  
  detach( monocle_params )
  
  return(mobj)  
}



#' Transfer data from a pseudotime modeling object to a Seurat object.
#' 
#' @details This function takes Seurat object and a finished object from a pseudotime analysis package
#' It transfers info from the latter to the former, guaranteeing that the Seurat object
#' will have complete metadata fields of the following names:
#' `pseudotime`, `branch`, `branch_viz_1`, `branch_viz_2`
#' @export
add_pseudotime_to_seurat = function(dge, pt_obj, pt_method = "monocle" ){
  if( pt_method == "monocle" ){
    to_add = data.frame(
      branch_viz_1 = pt_obj@reducedDimS[1, dge@cell.names], 
      branch_viz_2 = pt_obj@reducedDimS[2, dge@cell.names],
      branch = as.character( pt_obj@phenoData[["State"]] ),
      pseudotime = pt_obj@phenoData[["Pseudotime"]]       ) 
    atae( dge@cell.names,  rownames( pt_obj@phenoData@data ) )
    rownames( to_add ) = dge@cell.names
  } else if (pt_method == "dpt") {
    to_add = data.frame(
      branch_viz_1 = pt_obj@dm@eigenvectors[, 1],
      branch_viz_2 = pt_obj@dm@eigenvectors[, 2],
      branch_viz_3 = pt_obj@dm@eigenvectors[, 3],
      branch = as.character( pt_obj@branch [, 1] ),
      pseudotime = pt_obj[, 1]        )
    to_add$pseudotime = to_add$pseudotime*sign(cor(to_add$pseudotime, FetchData(dge, "eday")))
    rownames( to_add ) = dge@cell.names
    warning("Removing root cell from DPT output due to weird gap from rest of dataset.")
    bad_idx = which( to_add$pseudotime == 0 )
    to_add = to_add[-bad_idx, ]
    dge = Seurat::SubsetData( dge, cells.use = rownames( to_add ) )
  } else {
    warning("Only transfer from monocle and dpt has been implemented. Returning seurat object untouched.")
  }
  dge = Seurat::AddMetaData(dge, to_add)

  return( dge )
}


#' Project cells from one Seurat object onto the pseudotime axis (or any other continuous variable) from another Seurat object.
#'
#' @param dge_train Seurat object with 'pseudotime' field included.
#' @param dge_test Seurat object. Cells will be assigned a pseudotime value based on `dge_train`.
#' @param genes.use Defaults to `dge_train@var.genes`. If you use method="PCA", use the same genes here.
#' If you aren't sure, use method="KNN".
#' @param regressout Please put whatever you used for 'latent.vars' in RegressOut when you processed 'dge_train'.
#' @param to_project Character vector such that FetchData( dge_train, to_project) contains numeric columns.
#' Try "pseudotime" or c("tSNE_1", "tSNE_2").
#' @param k Number of nearest neighbors to use.
#' @param pc.use which principal components to use for nearest neighbor search.
#'
#' @export
#'
ProjectCells = function( dge_train, dge_test, 
                         genes.use = dge_train@var.genes,
                         regressout = c( "IG1.S", "S", "G2.M", "M", "M.G1" ),  
                         to_project = "pseudotime",
                         k = 25,
                         pc.use = 1:25, 
                         test_mode = F ){
  
  if(any(pc.use > ncol(dge_train@pca.x))){
    stop("Not enough principal components available in trainset. Compute moar. \n")
  }
  
  # Assemble data
  train_genes = FetchDataZeroPad( dge_train, genes.use ) %>% as.matrix
  test_genes  = FetchDataZeroPad( dge_test,  genes.use ) %>% as.matrix
  train_cc = FetchDataZeroPad( dge_train, regressout ) %>% as.matrix %>% cbind(1)
  test_cc  = FetchDataZeroPad( dge_test,  regressout ) %>% as.matrix %>% cbind(1)
  
  # # Regress out cc from testset using trainset params
  cc_coeffs = solve( t(train_cc) %*% train_cc, ( t(train_cc) %*% train_genes ) )
  train_genes = train_genes - train_cc %*% cc_coeffs
  test_genes  = test_genes  - test_cc  %*% cc_coeffs
  
  # # Standardize all data based on trainset statistics
  gene_means = apply( train_genes, 2, mean )
  gene_sds   = apply( train_genes, 2, sd )
  atae( length(gene_means), length(genes.use))
  train_genes = train_genes %>% as.matrix %>%
    sweep( 2, STATS = gene_means, "-") %>%
    sweep( 2, STATS = gene_sds,   "/")
  test_genes = test_genes %>% as.matrix %>%
    sweep( 2, STATS = gene_means, "-") %>%
    sweep( 2, STATS = gene_sds,   "/")
  
  # # Project into principal subspace of training data
  projection_mat = as.matrix(dge_train@pca.x[paste0("PC", pc.use)])
  train_pca_embeddings = train_genes %*% projection_mat
  test_pca_embeddings  = test_genes  %*% projection_mat
  
  # # Check
  if( test_mode ){
    plot( FetchData(dge_train, "PC1")[[1]], train_pca_embeddings[, 1], 
          main = "Existing versus recomputed PC1 projections of cells", 
          xlab = "existing", ylab = "recomputed")
  }
  
  # # Interpolate via average over nearest neighbors
  train_response = FetchDataZeroPad( dge_train, to_project )
  test_response = data.frame(matrix(NA, nrow = nrow(test_genes), ncol = length(to_project)))
  colnames(test_response) = to_project
  rownames(test_response) = rownames(test_genes)
  neighbors = FNN::get.knnx( data = train_pca_embeddings, query = test_pca_embeddings, k = k )
  for(i in 1:nrow(test_response)){
    test_response[i, ] = colMeans(train_response[neighbors$nn.index[i, ], ])
  }

  dge_test %<>% AddMetaData( metadata = test_response )
  return( dge_test )
}

## ------------------------------------------------------------------------

#' Align simple (y-shaped) lineage trees represented as points in 2D space.
#' 
#' Output from Monocle is arbitrarily rotated and reflected.
#' For bootstrap analysis where I call Monocle on subsets of the data, it's useful to 
#' align all of the results so that cells are in roughly the same positions across bootstrap samples. 
#' This program does that for any two matrices mat_ref and mat_samp.
#' Rows are cells and columns are dimensions of some embedding.
#'
#' The reference and the query can have different numbers of columns.
#' They should have the same number of rows, or else 
#' the property all( rownames( mat_samp ) %in% rownames( mat_ref ) ).

#' The solution is based on linear regression. 
#' The return value is a projection of each column of mat_ref onto the column space of the 
#' augmented matrix [mat_samp | 1]. In other terms, mat_samp is translated and linearly combined to 
#' minimize the l2 distance to each column of mat_ref, then returned.
#' This function always returns a numeric matrix with ncol = ncol(mat_ref) and nrow = nrow(mat_samp).
#' @export
align_embedding_to_reference = function(mat_samp, mat_ref, do.plot = F){
  
  # # Handle vectors
  if( is.null( dim ( mat_samp ) ) ) { mat_samp = matrixify_preserving_rownames( mat_samp ) }
  if( is.null( dim ( mat_ref  ) ) ) { mat_ref  = matrixify_preserving_rownames( mat_ref  ) }
  
  # # Handle dataframes
  mat_samp = as.matrix( mat_samp )
  mat_ref  = as.matrix( mat_ref  )

  # # Check dimensions
  if( ncol( mat_samp ) > nrow( mat_samp ) ||
      ncol( mat_ref  ) > nrow( mat_ref  ) ){
    warning("Did you forget to transpose?")
  }
  
  # # Ensure rownames present and compatible between sample and reference
  if( is.null( rownames( mat_samp ) ) |  is.null( rownames( mat_ref ) ) ){
    if( nrow( mat_samp ) != nrow( mat_ref ) ){ 
      stop( "Either number of rows should be equal or rows should be named with `all( rownames( mat_samp ) %in% rownames( mat_ref ) )`." ) 
    } else {
      rownames( mat_samp ) = as.character( 1:nrow( mat_samp ) )
      rownames( mat_ref )  = as.character( 1:nrow( mat_samp ) )
    }
  }
  atat( !is.null( rownames( mat_samp ) ) )
  atat( !is.null( rownames( mat_ref  ) ) )
  atat( all( rownames( mat_samp ) %in% rownames( mat_ref ) ) )
  mat_ref = mat_ref[ rownames( mat_samp ) , , drop = F] # drop = F GODDAMMIT

  # # Do math
  mat_return = matrix( lm( mat_ref ~ mat_samp )$fitted.values, ncol = ncol( mat_ref ) )
  
  # # guarantee nice named-matrix output
  rownames( mat_return ) = rownames( mat_samp ) 
  atat( nrow( mat_return ) == nrow( mat_samp ) )
  atat( ncol( mat_return ) == ncol( mat_ref ) )
  atat( is.numeric( mat_return ) )
  
  # # For testing
  if(do.plot){
    mat_ref     = data.frame( mat_ref );    mat_ref$type    = "input_ref"
    mat_samp    = data.frame( mat_samp );   mat_samp$type   = "input_sample"
    mat_return  = data.frame( mat_return ); mat_return$type = "return_val"
    df_all = rbind( mat_ref, mat_samp, mat_return )
    p = ggplot(df_all) + ggtitle("return_val and input_ref should be aligned") +
      geom_point(aes(X1, X2, colour = type), alpha = 0.2)
    print(p)
  } 
  return( mat_return )
}

# # Brief test script
# # One dimensional test
# refscale = 1:5; names(refscale) = letters[1:5]
# atae( matrixify_preserving_rownames( refscale ), 
#       align_embedding_to_reference(mat_ref = refscale, mat_samp = -10*refscale + 4 ),
#       tolerance = 0.000001,
#       check.attributes = F)
# # # Two dimensional test
# noisy_l  = cbind(-99:100, c( rep(1, 100), 1:100 ) - 50 ) + matrix( 2*rnorm(400), ncol = 2 )
# rownames( noisy_l ) = make.unique( letters[(1:200 %% 25) + 1] )
# noisy_l_rot = noisy_l %*% matrix( c(1, 1, -1, 1), nrow = 2 ) + matrix( 2*rnorm(400), ncol = 2 )
# tmp = align_embedding_to_reference(mat_ref  = noisy_l,
#                                    mat_samp = noisy_l, 
#                                    do.plot = T)
# tmp = align_embedding_to_reference(mat_ref  = noisy_l,
#                                    mat_samp = noisy_l_rot, 
#                                    do.plot = T)


