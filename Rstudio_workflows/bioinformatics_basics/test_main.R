#!/usr/bin/Rscript
# if you change the name of your script, this line must be changed as well
source("main.R")
library(testthat)

test_that("loading csv works using load_expression()", {
  result_tib <- load_expression("data/example_intensity_data_subset.csv")
  expect_equal(dim(result_tib), c(1000, 36))
  expect_true(is_tibble(result_tib))
})

test_that("the correct rows are filtered out in filter_15()", {
  test_tib <- tibble(probe=c('1_s_at', '2_s_at', '3_s_at', '4_s_at'),
                     GSM1=c(1.0, 3.95, 4.05, 0.5),
                     GSM2=rep(1.6, 4),
                     GSM3=rep(2.5, 4),
                     GSM4=rep(3.99, 4),
                     GSM5=rep(3.0, 4),
                     GSM6=rep(1.0, 4),
                     GSM7=rep(0.5, 4))
  expect_equal(c(filter_15(test_tib)$probe), c("2_s_at", "3_s_at"))
})

test_that("affy ids can be converted to HGNC names properly using affy_to_hgnc()", {
  # biomaRt super buggy so we can try to account for not connecting well
  success <- FALSE
  attempt <- 0
  while (!success ) {
    try({
      attempt <- attempt + 1
      tryCatch({
        response <- affy_to_hgnc(tibble('1553551_s_at'))},
        warning = function(w) stop("Could not connect to ENSEMBL. Reattempting connection")
      )
      Sys.sleep(1)
      success <- TRUE
    })
    if (attempt == 3) stop("Process halted. Function faced non-connection errors. Try again later.")
  }
  expect_equal(response$hgnc_symbol, c("MT-ND1", "MT-TI", "MT-TM", "MT-ND2"))
})

test_that("reduce_data() is correctly changing the size and shape of the tibble", {
  t_tibble <- tibble(probe = c("1_s_at", "2_s_at", "3_s_at"),
                     GSM1 = c(9.5, 7.6, 5.5),
                     GSM2 = c(9.7, 7.2, 2.9),
                     GSM3 = c(6.9, 4.3, 6.8))
  names <- tibble(affy_hg_u133_plus_2 = c("1_s_at", "2_s_at", "3_s_at"),
                  hgnc_symbol = c("A-REAL-GENE", "SONIC", "UTWPU"))
  good <- c("A-REAL-GENE")
  bad <- c("SONIC")
  reduce_test <- reduce_data(t_tibble, names, good, bad)
  result <- tibble(probe = c("1_s_at", "2_s_at"),
                   hgnc_symbol = c("A-REAL-GENE", "SONIC"),
                   gene_set = c("good", "bad"),
                   GSM1 = c(9.5, 7.6),
                   GSM2 = c(9.7, 7.2),
                   GSM3 = c(6.9, 4.3))
  expect_equal(reduce_test, result)
})

test_that("convert_to_long() correctly converts from long to wide format", {
  
  # wide format tibble
  wide_tib <- tibble(probe = c("202274_at", "202541_at", "202542_s_at", "203919_at"),
                         hgnc = c("ACTG2", "AIMP1", "AIMP1", "TCEA2"),
                         gene_set = rep("good", 4),
                         GSM1 = c(8.05, 8.40, 9.55, 4.44),
                         GSM2 = c(7.74, 7.11, 8.48, 5.39))
  
  # expected long format tibble
  expected_tib <- tibble(probe = rep(c("202274_at", "202541_at", "202542_s_at", "203919_at"), each = 2),
                     hgnc = rep(c("ACTG2", "AIMP1", "AIMP1", "TCEA2"), each = 2),
                     gene_set = rep(rep("good", 4), each = 2),
                     sample = rep(c("GSM1", "GSM2"), 4),
                     value = c(8.05, 7.74, 8.40, 7.11, 9.55, 8.48, 4.44, 5.39))
  
  # Convert to long format
  long_format <- convert_to_long(wide_tib)
  
  # Test if the converted long format matches the expected tibble
  expect_equal(long_format, expected_tib)
})