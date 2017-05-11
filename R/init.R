package_code_r = c( "setup.R", 
                    "utilities.R",
                    "seurat_extra.R",
                    "data_handling.R",
                    "quality_control.R",
                    "exploratory_pipeline.R",
                    "classifier.R",
                    "pseudotime.R" )
package_code_rmd = paste0( package_code_r, "md" )
for(i in seq_along(package_code_r)){
  knitr::purl(file.path( "R", package_code_rmd[[i]] ), output = file.path( "R", package_code_r  [[i]] ) )
  source(                                                       file.path( "R", package_code_r  [[i]] ) )
}

