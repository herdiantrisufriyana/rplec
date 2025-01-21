library(testthat)

test_that("ipla validates inputs correctly", {
  # Mock data for valid input
  aging <- get("aging", envir = environment())
  ga <- get("ga", envir = environment())
  phenotype <- get("phenotype", envir = environment())
  
  case <- "Case"
  control <- "Control"
  method <- "Mann-Whitney U"
  from <- 10
  to <- 20
  
  # Test valid input
  expect_error(ipla(aging, ga, phenotype, case, control, method, from, to), NA)
  
  # Test 'aging' is not a data frame
  expect_error(
    ipla(as.matrix(aging), ga, phenotype, case, control, method, from, to),
    "'aging' must be a data frame."
  )
  
  # Test 'aging' missing 'output' column
  aging_missing_output <- aging[, -1, drop = FALSE]
  expect_error(
    ipla(aging_missing_output, ga, phenotype, case, control, method, from, to),
    "'aging' must have a column named 'output'."
  )
  
  # Test 'ga' is not a data frame
  expect_error(
    ipla(aging, as.matrix(ga), phenotype, case, control, method, from, to),
    "'ga' must be a data frame."
  )
  
  # Test 'ga' missing 'GA' column
  ga_missing_ga_column <- ga[, -1, drop = FALSE]
  expect_error(
    ipla(
      aging, ga_missing_ga_column, phenotype, case, control, method, from, to
    ),
    "'ga' must have a column named 'GA'."
  )
  
  # Test 'phenotype' is not a data frame
  expect_error(
    ipla(aging, ga, as.matrix(phenotype), case, control, method, from, to),
    "'phenotype' must be a data frame if provided."
  )
  
  # Test 'phenotype' missing 'phenotype' column
  phenotype_missing_column <- phenotype[, -1, drop = FALSE]
  expect_error(
    ipla(aging, ga, phenotype_missing_column, case, control, method, from, to),
    "'phenotype' must have a column named 'phenotype'."
  )
  
  # Test invalid 'case'
  expect_error(
    ipla(aging, ga, phenotype, case = 123, control, method, from, to),
    "'case' must be a single character string."
  )
  
  # Test invalid 'control'
  expect_error(
    ipla(aging, ga, phenotype, case, control = 123, method, from, to),
    "'control' must be a single character string."
  )
  
  # Test invalid 'method'
  expect_error(
    ipla(
      aging, ga, phenotype, case, control, method = "InvalidMethod", from, to
    ),
    "'method' must be either 'Mann-Whitney U' or 'Permutation' if provided."
  )
  
  # Test invalid 'from'
  expect_error(
    ipla(aging, ga, phenotype, case, control, method, from = 50, to),
    "'from' must be a single integer between 5 and 44."
  )
  
  # Test invalid 'to'
  expect_error(
    ipla(aging, ga, phenotype, case, control, method, from, to = 50),
    "'to' must be a single integer between 5 and 44."
  )
  
  # Test 'from' greater than 'to'
  expect_error(
    ipla(aging, ga, phenotype, case, control, method, from = 20, to = 10),
    "'from' cannot be greater than 'to'."
  )
  
  # Test valid input with null 'method'
  expect_error(
    ipla(aging, ga, phenotype, case, control, method = NULL, from, to), NA
  )
  
  # Test valid input with null 'from' and 'to'
  expect_error(
    ipla(aging, ga, phenotype, case, control, method, from = NULL, to = NULL)
    , NA
  )
})