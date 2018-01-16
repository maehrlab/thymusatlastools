library(magrittr)
# How I made this test dataset: 30 cells of each cluster from the thymus atlas.
# test_data_path  = "~/Desktop/software_projects/thymusatlastools/tests/testdata/thymus_test.Rdata"
# thymus_test_with_junk =
#   freezr::inventory_get(inv_location = "~/Desktop/scRNA_data_analysis/atlas_paper_2017/scRNA_redo_12_p0/results/",
#                         tag = "overview_clean_labeled"  ) %>%
#   readRDS %>%
#   SubsetData(max.cells.per.ident = 30)
# thymus_test = thymus_test_with_junk
#   deseuratify_raw_data %>%
#   seuratify_thy_data
# thymus_test@data.info = thymus_test_with_junk@data.info
# saveRDS(thymus_test, test_data_path)
#
# load a minimal example data set (subset of thymus atlas)
thymus_test = readRDS("../testdata/thymus_test.Rdata")

testthat::test_that( "AvailableData runs and gives back genes and essentials in a bare-bones object", {
  should_have = c( "nGene"    ,     "nUMI"      ,    "orig.ident" , "ident", "eday",
                   rownames(thymus_test@data))
  atat(all(should_have %in% AvailableData(thymus_test)))
})


testthat::test_that( "FetchDataZeroPad works for things we have and things we lack", {
  have = c("nGene"    ,     "Foxn1" )
  lack = "iq34gobafrovurbryeubv78b3"
  expect_warning(       FetchDataZeroPad( thymus_test,          lack ) )
  X = suppressWarnings( FetchDataZeroPad( thymus_test, c( have, lack ) ) )
  expect_equal( X[1:2], FetchData( thymus_test, have ) )
  expect_equal( X[[3]]  , rep( 0, nrow( X ) ) )
  expect_equal( 3, ncol( X ) )
})


testthat::test_that( "FetchDataZeroPad works for things we have and things we lack", {
  queries = c( "nGene"    ,     "Foxn1"      ,    "iq34gobafrovurbryeubv78b3" )
  expect_warning(      FetchDataZeroPad(thymus_test, queries))
  X = suppressWarnings(FetchDataZeroPad(thymus_test, queries))
  atat(all(0==X$iq34gobafrovurbryeubv78b3))
  atae(ncol(X), 3)
})

FetchDataZeroPad

testthat::test_that( "SubsetDataFlex works", {
  thymus_test %<>% SubsetDataFlex("eday", "eday==12.5")
  thymus_test@data.info$eday %>% table %>% equals(12.5)
})

testthat::test_that( "FindMarkersFlex matches FindMarkers", {
  gu = thymus_test@data %>% rownames %>% sample(100)
  X = FindMarkersFlex(thymus_test,
                      ident.use = "ident", test.use = "bimod",
                      ident.1 = "C1", ident.2 = "C2", genes.use = gu )
  Y = FindMarkers    (thymus_test,
                      ident.1 = "C1", ident.2 = "C2", genes.use = gu )
  Y = Y[order(-Y$avg_diff), ]
  expect_equal(X[1:4], Y)
})


testthat::test_that( "SeuratPie runs", {
  SeuratPie(thymus_test, ident.use = "ident" )
})


testthat::test_that( "TACS runs with and without faceting and cutoffs", {
  TACS(thymus_test, gene1 = "Epcam", gene2 = "Ptprc" )
  TACS(thymus_test, gene1 = "Epcam", gene2 = "Ptprc", facet_by = "eday" )
  TACS(thymus_test, gene1 = "Epcam", gene2 = "Ptprc", facet_by = "eday", cutoffs = c(0.7, 0.7) )
  TACS(thymus_test, gene1 = "Epcam", gene2 = "Ptprc", facet_by = "eday", cutoffs = c(0.7, 0.7), density = T )
})

testthat::test_that( "", {
  TACS(thymus_test, gene1 = "Epcam", gene2 = "Ptprc" )
})


