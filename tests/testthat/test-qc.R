library(testthat)

test_that("qc validates inputs correctly", {
  # Create mock data for valid input
  beta <- download_beta_values_case()
  beta <- beta[, 1:2, drop = FALSE]
  norm_beta <- bmiq_norm_450k(beta)
  dnam_ga <- plec(norm_beta)
  
  ga <- get("ga", envir = environment())
  ga <- ga[colnames(beta), , drop = FALSE]
  
  phenotype <- get("phenotype", envir = environment())
  phenotype <- phenotype[colnames(beta), , drop = FALSE]
  
  # Test valid input
  expect_error(qc(dnam_ga, ga, phenotype), NA)
  
  # Test dnam_ga is not a data frame
  expect_error(
    qc(as.matrix(dnam_ga), ga, phenotype),
    "'dnam_ga' must be a data frame."
  )
  
  # Test dnam_ga missing 'output' column
  dnam_ga_missing_output <- dnam_ga[, -1, drop = FALSE]
  expect_error(
    qc(dnam_ga_missing_output, ga, phenotype),
    "'dnam_ga' must have a column named 'output'."
  )
  
  # Test ga is not a data frame
  expect_error(
    qc(dnam_ga, as.matrix(ga), phenotype),
    "'ga' must be a data frame."
  )
  
  # Test ga missing 'GA' column
  ga_missing_ga_column <- ga[, -1, drop = FALSE]
  expect_error(
    qc(dnam_ga, ga_missing_ga_column, phenotype),
    "'ga' must have a column named 'GA'."
  )
  
  # Test phenotype is not a data frame
  expect_error(
    qc(dnam_ga, ga, as.matrix(phenotype)),
    "'phenotype' must be a data frame."
  )
  
  # Test phenotype missing 'phenotype' column
  phenotype_missing_column <- phenotype[, -1, drop = FALSE]
  expect_error(
    qc(dnam_ga, ga, phenotype_missing_column),
    "'phenotype' must have a column named 'phenotype'."
  )
  
  # Test valid input with null phenotype
  expect_error(qc(dnam_ga, ga, phenotype = NULL), NA)
})