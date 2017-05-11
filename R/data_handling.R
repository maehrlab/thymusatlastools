## ------------------------------------------------------------------------
#' Rescale every cell to a certain amount of UMIs, where
#' that amount is selected by rounding up the median UMI count up to the next power of 10.
#'
#'
normalize_cpx_amt = function(dge, results_path, do.plot = T){

  umis_by_cell = apply( dge, 2, sum )
  assertthat::are_equal(length(umis_by_cell), length(colnames(dge)))
  magnitude = umis_by_cell %>% median %>% log10 %>% ceiling
  to.return = 10^magnitude 
  
  # # Plot results
  if( do.plot ){
      dir.create.nice( file.path( results_path, "QC" ) )
    pdf( file.path( results_path, "QC", "total_umis_by_cell.pdf"))
    {
      hist(log10(umis_by_cell), breaks = 40,
           xlab = "log10 UMI count", main = "Number of UMIs by cell")
      abline(v = to.return)
    }
    dev.off()
  }
  
  return(to.return)
}
demo1 = matrix(2000, nrow = 5, ncol = 3)
assertthat::are_equal(normalize_cpx_amt(demo1, results_path = "~/Desktop/scRNA_junk", do.plot = F), 10000)


#' Make arrays into Seurat objects.
#' 
#' Keeps all genes expressed in at least
#' (by default) 3 cells and keeps all cells with at least 1000 genes. It reports some 
#' summary figures, plotting number of genes by cell, num UMIs by cell, and number of cells by
#' gene. 
seuratify_thy_data = function(raw_dge, results_path, test_mode = F, 
                              min.genes = 1000, min.cells = 3, do.plot = T){
  atat(1 < raw_dge %>% dim %>% min %>% min)

  total_desired = normalize_cpx_amt( raw_dge, results_path )
  raw_dge_norm = apply( X = raw_dge, MARGIN = 2, FUN = div_by_sum ) * total_desired
  raw_dge_norm = as( raw_dge_norm, "sparseMatrix" )
  seurat_dge = Setup( new( "seurat", raw.data = log2( 1+raw_dge_norm ) ),
                      min.cells = min.cells * (1 - test_mode), # 0 if test_mode; else 3
                      min.genes = min.genes * (1 - test_mode), # 0 if test_mode; else min.genes
                      is.expr=0, 
                      do.logNormalize = F, 
                      project = "thymus_scRNAseq", 
                      names.delim = "\\|",
                      names.field = 2 )
  
  # # Make sure the raw data is just UMI counts and update the nUMI field, which will be
  # # filled incorrectly given the above
  seurat_dge@raw.data = raw_dge[ seurat_dge@data %>% rownames, 
                                 seurat_dge@data %>% colnames ]
  seurat_dge %<>% AddMetaData( col.name = "nUMI", 
                               metadata = setNames( colSums( seurat_dge@raw.data ),
                                                    colnames( seurat_dge@raw.data ) ) )
  
  # # Plot results
  if(do.plot){
    dir.create.nice( file.path(results_path, "QC" ) )
    genes_by_cell = apply( raw_dge, 2, nnz )
    atat( length( genes_by_cell ) == ncol( raw_dge ) ) # num cells = num cols
    pdf(file.path(results_path, "QC", "total_genes_by_cell.pdf"))
    {
      hist(log10(genes_by_cell), breaks = 40,
           xlab = "log10 gene count", main = "Number of genes by cell")
      abline(v = log10( min.genes ) )
    }
    dev.off()
    
    cells_by_gene = apply( raw_dge, 1, nnz )
    atat( length( cells_by_gene ) == nrow( raw_dge )[1] ) # num genes = num rows
    pdf(file.path(results_path, "QC", "total_cells_by_gene.pdf"))
    {
      hist(log10(cells_by_gene), breaks = 40,
           xlab = "log10 cell count", main = "Number of cells by gene")
      abline(v = log10( min.cells ) )
    }
    dev.off()
    
    print( paste0("There are ", 
                  length(cells_by_gene), " genes and ", 
                  length(genes_by_cell), " cells in the raw data."))
    print( paste0( length( cells_by_gene ) - nrow( seurat_dge@data ), " genes and ", 
                   length( genes_by_cell ) - ncol( seurat_dge@data ), " cells were excluded from the Seurat object."))
    
  }

  return(seurat_dge)
}


#' Extract mostly-raw (medium rare?) data from a Seurat object.
#' 
#' Seurat preserves the raw data exactly. Sometimes that's not ideal.
#' This function helps you get rawish data that have undergone the same QC filters as the Seurat scale.data,
#' so some cells and genes are filtered out. 
#' But, the numbers are UMI counts, integers, not logged or with any of that normalization BS.
#' This function guarantees output with `colnames(output) == dge@cell.names`.
#' 
#' Because the `@raw.data` slot was filled in wrong in some of my Seurat objects,
#' this can use `load_thymus_profiling_data` to get the raw data.
#' To toggle this behavior, set `retrieve_anew = {T,F}`.
deseuratify_raw_data = function( seurat_dge, retrieve_anew = F ){
  if( !retrieve_anew ){
    raw_dge = seurat_dge@raw.data
  } else {
    raw_dge = load_thymus_profiling_data( sample_ids = unique( FetchData(seurat_dge, "orig.ident" )[[1]]) ) %>% dge_merge_list
  }
  desired_genes = rownames( seurat_dge@scale.data )
  acceptable_genes = intersect( desired_genes, rownames( raw_dge ) )
  missing_genes    = setdiff(   desired_genes, rownames( raw_dge ) )
  raw_dge = as.matrix( raw_dge )[acceptable_genes, seurat_dge@cell.names]
  
  ##Zero-pad to ensure all genes from scale.data are present
  my_zeroes = matrix( 0, ncol = ncol( raw_dge ), nrow = length( missing_genes ) )
  rownames( my_zeroes ) = missing_genes
  colnames( my_zeroes ) = colnames( raw_dge )
  raw_dge = rbind( raw_dge, my_zeroes )
  raw_dge = raw_dge[desired_genes, ]
  
  atae(raw_dge, round( raw_dge ) )
  atae(desired_genes, rownames( raw_dge ) )
  return( raw_dge )
}

## ------------------------------------------------------------------------
#' Merge a list of digital gene expression matrices.
#'
#' @param dge_list list of matrices with genes as rows and cell barcodes as columns. Duplicate barcodes across datasets cause errors; I recommend you append the sample ID.
#' @details `dge_merge_list` converts the digital gene expression matrices to dataframes with genes as columns, merges them, then converts the result back, all without disturbing the gene labels. A gene will be included if it appears in any of the datasets. If a gene appears in one dataset but not another, zeroes will be filled in for missing expression levels.
dge_merge_list = function(dge_list){
  # Allow duplicate cells but not duplicate genes
  all_genes = Reduce( f = union, x = lapply( dge_list, rownames ) )
  all_cells = Reduce( f = c,     x = lapply( dge_list, colnames ) )
  if(anyDuplicated(all_cells)){ warning( "Duplicate cell barcodes present." )}
  new_dge = as.data.frame( matrix( 0, nrow = length( all_genes ), 
                                      ncol = length( all_cells ) ) )
  rownames( new_dge ) = all_genes
  colnames( new_dge ) = all_cells
  for( ii in seq_along( dge_list ) ){
    dge = dge_list[[ii]]
    new_dge[rownames(dge), colnames(dge)] = dge
    print( paste( "Finished merging dge matrix", names(dge_list)[[ii]] ) )
  }
  #Results should have at least as many genes as the input with the most genes
  #the number of barcodes should equal the sum of the total barcodes.
  input_dimensions = Reduce(rbind, lapply(dge_list, dim))
  assertthat::assert_that(dim(new_dge)[1] >= max(input_dimensions[,1]))
  assertthat::assert_that(dim(new_dge)[2] == sum(input_dimensions[,2]))
  assertthat::assert_that(!is.null(colnames(new_dge)))
  assertthat::assert_that(!is.null(rownames(new_dge)))
  return( new_dge )
}


