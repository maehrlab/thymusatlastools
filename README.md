##README

####Paper Reproduction 

This folder contains an R package with code to reproduce the analyses in (the as-yet-unpublished Maehrlab scRNA thymus atlas manuscript). It complements the R package `thymusatlasdatapublic` and the free-standing scripts in `ekernf01/thymusatlasanalysis`. To get started, open up `main.Rmd` from `ekernf01/thymusatlasanalysis`.


####New work

The R package contained in this folder is built around [Seurat](http://satijalab.org/seurat/), and it adds considerable functionality. Notable features: supervised classification of cell types and an interface with Monocle. Consult the package index and documentation (`?knn_classifier`, `?custom_feature_plot`, `?FetchDataZeroPad`, `?get_ortholog`) for more information.