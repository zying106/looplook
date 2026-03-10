# tests/testthat/test-utils.R

test_that("read_robust_general pushes all I/O boundary conditions", {
  tmp_tsv <- tempfile(fileext = ".txt")
  write.table(data.frame(Gene = c("A", "B"), Val1 = 1:2, Val2 = 3:4), tmp_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  res_tsv <- read_robust_general(tmp_tsv, header = TRUE, row_name = 1, min_cols = 2)
  expect_s3_class(res_tsv, "data.frame")

  tmp_csv <- tempfile(fileext = ".csv")
  write.csv(data.frame(Col1 = 1:3, Col2 = 4:6, Col3 = 7:9), tmp_csv, row.names = FALSE)
  res_csv <- read_robust_general(tmp_csv, header = TRUE, row_name = NULL, min_cols = 2)
  expect_equal(ncol(res_csv), 3)

  tmp_dirty <- tempfile(fileext = ".txt")
  writeLines(c("# This is a comment", "Gene\tExp1\tExp2", "BRD4\t10\t20", "MYC\tNA\t30"), tmp_dirty)
  res_dirty <- read_robust_general(tmp_dirty, header = TRUE, row_name = 1, min_cols = 2)
  expect_equal(nrow(res_dirty), 2)

  expect_error(read_robust_general("fake_file.txt"), "not found")
  expect_error(read_robust_general(tmp_tsv, min_cols = 100), "insufficient columns")

  unlink(c(tmp_tsv, tmp_csv, tmp_dirty))
})
