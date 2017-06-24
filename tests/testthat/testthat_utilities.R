testthat::test_that( "aggregate_nice works for means and sums", {
  test1 = data.frame(1:6, 2, 3)
  atae( aggregate_nice( test1, letters[1:6], mean ), test1 )
  atat( all( aggregate_nice( test1, rep("A", 6), sum ) == matrix( c(21, 12, 18), nrow = 1) ) )
})


testthat::test_that( "replace_with_int_rank works for integers and characters", {
  atat(all(replace_with_int_rank(1:5) == 1:5))
  atat(all(replace_with_int_rank(c("a", "b", "a", "c", "c")) == c(1, 2, 1, 3, 3)))
  atat(all(replace_with_int_rank(c("0", "2", "0", "4", "4")) == c(1, 2, 1, 3, 3)))
})

testthat::test_that( "replace_with_int_rank works for integers and characters", {
  assertthat::are_equal(Capitalize(c("FOXN1", "PSMB11")), c("Foxn1", "Psmb11") )
})


testthat::test_that( "strip_suffix works for removing '.pdf'", {
  atae( strip_suffix("blah.pdf", ".pdf"), "blah")
  atae( strip_suffix("blah",     ".pdf"), "blah")
})


testthat::test_that( "get_preimage works for simple examples", {
  atae( get_preimage( map = setNames( LETTERS, letters) ), 
        as.list( setNames( letters, LETTERS) )  )
  to_invert = setNames(      c("A", "B", "DUPE", "DUPE", "DUPE2", "DUPE2"), 
                             nm = c("a", "b", "c",    "d", "e",    "f") )
  inverse = list(A = "a", 
                 B = "b",
                 DUPE  = c( "c", "d" ) ,
                 DUPE2 = c( "e", "f" ) )
  atae( get_preimage( map = to_invert ), inverse ) 
})

