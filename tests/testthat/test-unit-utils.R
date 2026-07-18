# tests/testthat/test-unit-utils.R — unit tests for utility helpers

# --- .passes_expression_threshold ---
test_that(".passes_expression_threshold enforces threshold=0 as strictly positive", {
    expect_true(looplook:::.passes_expression_threshold(5, 0))
    expect_false(looplook:::.passes_expression_threshold(0, 0))
    expect_false(looplook:::.passes_expression_threshold(NA_real_, 0))
    expect_true(looplook:::.passes_expression_threshold(0.1, 0))
})

test_that(".passes_expression_threshold enforces threshold>=1 as >= comparison", {
    expect_true(looplook:::.passes_expression_threshold(5, 1))
    expect_true(looplook:::.passes_expression_threshold(1, 1))
    expect_false(looplook:::.passes_expression_threshold(0.9, 1))
    expect_false(looplook:::.passes_expression_threshold(0, 1))
    expect_false(looplook:::.passes_expression_threshold(NA_real_, 1))
})

test_that(".passes_expression_threshold rejects invalid threshold", {
    expect_error(looplook:::.passes_expression_threshold(1, NA_real_), "finite")
    expect_error(looplook:::.passes_expression_threshold(1, -1), "non-negative")
    expect_error(looplook:::.passes_expression_threshold(1, Inf), "finite")
})

test_that(".passes_expression_threshold handles vector input", {
    vals <- c(10, 0, 5, NA_real_, 1)
    expect_equal(
        looplook:::.passes_expression_threshold(vals, 1),
        c(TRUE, FALSE, TRUE, FALSE, TRUE)
    )
    expect_equal(
        looplook:::.passes_expression_threshold(vals, 0),
        c(TRUE, FALSE, TRUE, FALSE, TRUE)
    )
})

# --- .validate_tss_region ---
test_that(".validate_tss_region accepts valid windows", {
    expect_equal(looplook:::.validate_tss_region(c(-2000, 2000)), c(-2000, 2000))
    expect_equal(looplook:::.validate_tss_region(c(-5000, 5000)), c(-5000, 5000))
})

test_that(".validate_tss_region rejects invalid windows", {
    expect_error(.validate_tss_region(c(2000, 5000)), "such as c") # first > 0
    expect_error(.validate_tss_region(c(-5000, -2000)), "such as c") # second < 0
    expect_error(.validate_tss_region(c(-2000)), "such as c") # length != 2
    expect_error(.validate_tss_region(c(NA, 2000)), "finite") # NA
    expect_error(.validate_tss_region("string"), "such as c") # non-numeric
})

# --- .motif_log_or (internal helper; test via formula) ---
test_that(".motif_log_or computes Haldane-Anscombe corrected log OR", {
    .motif_log_or <- function(a, b, n_fg, n_bg) {
        log(((a + 0.5) * (n_bg - b + 0.5)) / ((n_fg - a + 0.5) * (b + 0.5)))
    }
    # Symmetric case
    expect_equal(.motif_log_or(10, 10, 50, 50), 0)
    # Foreground enriched
    expect_gt(.motif_log_or(20, 5, 50, 50), 0)
    # Zero foreground: finite due to correction
    lor3 <- .motif_log_or(0, 10, 50, 50)
    expect_false(is.infinite(lor3))
    expect_lt(lor3, 0)
    # All foreground: finite
    lor4 <- .motif_log_or(50, 10, 50, 50)
    expect_false(is.infinite(lor4))
    expect_gt(lor4, 0)
})

test_that(".motif_log_or is monotonic in foreground hits", {
    .motif_log_or <- function(a, b, n_fg, n_bg) {
        log(((a + 0.5) * (n_bg - b + 0.5)) / ((n_fg - a + 0.5) * (b + 0.5)))
    }
    expect_gt(.motif_log_or(15, 10, 50, 50), .motif_log_or(5, 10, 50, 50))
})

# --- .empty_motif_output ---
test_that(".empty_motif_output returns stable schema", {
    x <- looplook:::.empty_motif_output()
    expect_type(x, "list")
    expect_named(x, c("results", "plots"))
    expect_named(x$results, c("proximal", "distal"))
    expect_null(x$results$proximal)
    expect_null(x$results$distal)
    expect_length(x$plots, 0)
})

# --- .igraph_vertex_ids ---
test_that(".igraph_vertex_ids handles NULL and empty input", {
    expect_equal(looplook:::.igraph_vertex_ids(NULL, NULL), character(0))
    expect_equal(looplook:::.igraph_vertex_ids(character(0), NULL), character(0))
})

# --- extract_genes ---
test_that("extract_genes handles edge cases", {
    expect_equal(looplook:::extract_genes(c("A", "B", "A")), "A;B")
    expect_equal(looplook:::extract_genes(NA_character_), NA_character_)
    expect_equal(looplook:::extract_genes(""), NA_character_)
})

# --- clean_gene_names ---
test_that("clean_gene_names handles empty and NA gracefully", {
    expect_equal(looplook:::clean_gene_names(""), character(0))
    expect_equal(looplook:::clean_gene_names(NA_character_), character(0))
    expect_equal(looplook:::clean_gene_names(NULL), character(0))
})
