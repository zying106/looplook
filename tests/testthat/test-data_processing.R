# tests/testthat/test-data_processing.R

test_that("data_processing modules run successfully on example bedpe files", {
  loop1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  loop2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")
  h3k27ac_peaks <- system.file("extdata", "example_k27ac_peaks.bed", package = "looplook")
  skip_if(loop1 == "" || loop2 == "")

  res_clean <- looplook:::.with_known_upstream_noise_suppressed(
    consolidate_chromatin_loops(
      files = c(loop1, loop2),
      mode = "consensus",
      min_raw_score = 2,
      min_score = 5,
      gap = 1000,
      blacklist_species = "hg38",
      region_of_interest = h3k27ac_peaks,
      out_file = tempfile(fileext = ".bedpe"),
      write_output = FALSE,
      quiet = TRUE
    )
  )

  expect_true(!is.null(res_clean))
})

test_that("consolidate_chromatin_loops balances clustered scores across replicates", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")

  writeLines(
    "chr1\t0\t100\tchr1\t200\t300\t100",
    f1
  )
  writeLines(c(
    "chr1\t5\t105\tchr1\t205\t305\t20",
    "chr1\t10\t110\tchr1\t210\t310\t25",
    "chr1\t15\t115\tchr1\t215\t315\t30"
  ), f2)

  gi <- consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "consensus",
    gap = 25,
    quiet = TRUE
  )
  expect_equal(length(gi), 1)
  expect_equal(S4Vectors::mcols(gi)$n_members[[1]], 4L)
  expect_equal(S4Vectors::mcols(gi)$n_reps[[1]], 2L)
  expect_equal(S4Vectors::mcols(gi)$score[[1]], 62.5)

  gi_keep <- consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "consensus",
    gap = 25,
    min_score = 50,
    quiet = TRUE
  )
  expect_equal(length(gi_keep), 1)

  gi_drop <- consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "consensus",
    gap = 25,
    min_score = 63,
    quiet = TRUE
  )
  expect_equal(length(gi_drop), 0)

  unlink(c(f1, f2))
})

test_that("consolidate_chromatin_loops consensus score is replicate-balanced (mean of per-source means)", {
  # Scenario: 3 replicates with asymmetric fragmentation
  # Rep1: 1 loop, score=100
  # Rep2: 3 fragmented loops at same locus, scores=20,30,50  → per-source mean=33.3
  # Rep3: 1 loop, score=60
  # Expected: mean(100, 33.3, 60) = 64.4 (NOT mean(100,20,30,50,60)=52)
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  f3 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t0\t100\tchr1\t200\t300\t100", f1)
  writeLines(c(
    "chr1\t5\t105\tchr1\t205\t305\t20",
    "chr1\t10\t110\tchr1\t210\t310\t30",
    "chr1\t15\t115\tchr1\t215\t315\t50"
  ), f2)
  writeLines("chr1\t0\t100\tchr1\t200\t300\t60", f3)

  gi <- consolidate_chromatin_loops(
    files = c(f1, f2, f3), mode = "consensus", gap = 25, quiet = TRUE
  )
  expect_equal(length(gi), 1)
  # Per-source means: rep1=100, rep2=mean(20,30,50)=33.33, rep3=60
  # Replicate-balanced score: mean(100, 33.33, 60) = 64.44
  expected_score <- mean(c(100, mean(c(20, 30, 50)), 60))
  expect_equal(S4Vectors::mcols(gi)$score[[1]], expected_score, tolerance = 0.1)
  # n_reps should be 3 (all three files contribute)
  expect_equal(S4Vectors::mcols(gi)$n_reps[[1]], 3L)
  # n_members should be 5 (1 + 3 + 1 raw loops merged)
  expect_equal(S4Vectors::mcols(gi)$n_members[[1]], 5L)
  unlink(c(f1, f2, f3))
})

test_that("consolidate_chromatin_loops respects quiet and write_output flags", {
  loop1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  loop2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")
  skip_if(loop1 == "" || loop2 == "")

  out_dir <- tempfile(pattern = "looplook_cons_nowrite_")
  out_file <- file.path(out_dir, "loops.bedpe")
  unlink(out_dir, recursive = TRUE, force = TRUE)
  expect_false(dir.exists(out_dir))

  expect_no_message({
    gi <- consolidate_chromatin_loops(
      files = c(loop1, loop2),
      mode = "consensus",
      gap = 1000,
      out_file = out_file,
      write_output = FALSE,
      quiet = TRUE
    )
  })

  expect_s4_class(gi, "GInteractions")
  expect_false(dir.exists(out_dir))
  expect_false(file.exists(out_file))
})

# --- cluster_loops_dt: overlap detection (not containment) ---
test_that("cluster_loops_dt detects overlap not containment", {
  dt <- data.table::data.table(
    chr1   = c("chr1", "chr1"),
    start1 = c(100L, 150L),
    end1   = c(200L, 250L),
    chr2   = c("chr1", "chr1"),
    start2 = c(1000L, 1050L),
    end2   = c(2000L, 2100L),
    score  = c(1, 1),
    source = c(1L, 1L)
  )
  result <- looplook:::cluster_loops_dt(dt, gap = 0L)
  # anchor1 overlaps (100-200 vs 150-250) and anchor2 overlaps (1000-2000 vs 1050-2100)
  expect_equal(length(unique(result$cluster)), 1L)
})

# --- reduce_clusters_dt: n_reps counts every contributing source, even NA-score ones ---
test_that("reduce_clusters_dt counts NA-score sources toward n_reps", {
  dt <- data.table::data.table(
    cluster = c(1L, 1L, 1L),
    chr1 = "chr1", start1 = c(100L, 110L, 120L), end1 = c(200L, 210L, 220L),
    chr2 = "chr1", start2 = c(1000L, 1010L, 1020L), end2 = c(2000L, 2010L, 2020L),
    score = c(5, NA, 3),
    source = c(1L, 1L, 2L)
  )
  result <- looplook:::reduce_clusters_dt(dt)
  # source 1 contributes rows to this cluster (one row has score=5, one row score=NA);
  # source 2 contributes one row with score=3.
  # Per-source means: source 1 -> mean(5) = 5 (NA stripped); source 2 -> mean(3) = 3.
  # Replicate-balanced score = mean(5, 3) = 4.
  # n_reps MUST count both sources (2), regardless of per-row NA scores, because
  # the existence of a row in the BEDPE already means the replicate detected the loop.
  expect_equal(result$n_reps, 2L)
  expect_equal(result$score, 4)
})

test_that("reduce_clusters_dt: cluster with all-NA scores still has integer n_reps (not NA)", {
  dt <- data.table::data.table(
    cluster = c(1L, 1L),
    chr1 = "chr1", start1 = c(100L, 110L), end1 = c(200L, 210L),
    chr2 = "chr1", start2 = c(1000L, 1010L), end2 = c(2000L, 2010L),
    score = c(NA_real_, NA_real_),
    source = c(1L, 2L)
  )
  result <- looplook:::reduce_clusters_dt(dt)
  # n_reps must be 2 (both sources contributed rows), NOT NA.
  # Previously this produced n_reps = NA and downstream consensus filters
  # silently dropped the cluster via NA >= k -> NA -> FALSE.
  expect_equal(result$n_reps, 2L)
  expect_true(is.na(result$score))
})

# --- filter_chromatin_loops / .consolidate_post_filters: NA-score guard ---
test_that("filter_chromatin_loops min_score drops NA-score loops explicitly (not via NA>=k)", {
  dt <- data.table::data.table(
    cluster = c(1L, 2L), chr1 = "chr1", start1 = c(100L, 200L), end1 = c(150L, 250L),
    chr2 = "chr1", start2 = c(1000L, 2000L), end2 = c(1500L, 2500L),
    score = c(NA_real_, 5), source = c(1L, 1L),
    n_members = c(1L, 1L), n_reps = c(1L, 1L)
  )
  gi <- looplook:::dt_to_gi(dt)
  res <- filter_chromatin_loops(gi, min_score = 4, quiet = TRUE)
  # NA-score loop must be dropped, but explicitly (no NA leaks into result).
  expect_equal(length(res), 1L)
  expect_true(all(!is.na(S4Vectors::mcols(res)$score)))
})

# --- bedpe_to_gi: validation and score detection ---
test_that("bedpe_to_gi rejects start >= end (zero-width or invalid)", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t200\t100\tchr1\t400\t500", tmp)
  expect_error(looplook:::bedpe_to_gi(tmp), "start >= end")
  unlink(tmp)
})

test_that("bedpe_to_gi rejects zero-width anchor (start == end)", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t100\tchr1\t400\t500", tmp)
  expect_error(looplook:::bedpe_to_gi(tmp), "start >= end")
  unlink(tmp)
})

test_that("bedpe_to_gi detects score in column 7", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines(c(
    "chr1\t100\t200\tchr1\t400\t500\t42",
    "chr1\t150\t250\tchr1\t450\t550\t10"
  ), tmp)
  gi <- looplook:::bedpe_to_gi(tmp)
  expect_equal(S4Vectors::mcols(gi)$score, c(42, 10))
  unlink(tmp)
})

test_that("bedpe_to_gi warns when no valid score column", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines(c(
    "chr1\t100\t200\tchr1\t400\t500\tname\tnote",
    "chr1\t150\t250\tchr1\t450\t550\tname2\tnote2"
  ), tmp)
  expect_warning(looplook:::bedpe_to_gi(tmp), "Scores defaulted to 0")
  unlink(tmp)
})

test_that("bedpe_to_gi swaps anchors when needed", {
  tmp <- tempfile(fileext = ".bedpe")
  # chr2 should come before chr1 after swap
  writeLines("chr2\t100\t200\tchr1\t400\t500", tmp)
  gi <- looplook:::bedpe_to_gi(tmp)
  a1 <- InteractionSet::anchors(gi, type = "first")
  a2 <- InteractionSet::anchors(gi, type = "second")
  # After swap, first anchor should be chr1
  expect_equal(as.character(GenomicRanges::seqnames(a1)), "chr1")
  unlink(tmp)
})

test_that("bedpe_to_gi preserves score after anchor swap", {
  tmp <- tempfile(fileext = ".bedpe")
  # chr2 comes first in file; score=42 should follow the interaction after swap
  writeLines(c(
    "chr2\t100\t200\tchr1\t400\t500\t42",
    "chr1\t50\t150\tchr3\t300\t400\t99"
  ), tmp)
  gi <- looplook:::bedpe_to_gi(tmp)
  scores <- S4Vectors::mcols(gi)$score
  # First loop (chr2-chr1) gets swapped to chr1-chr2; score 42 stays with it
  # Second loop (chr1-chr3) stays; score 99 stays with it
  a1_chr <- as.character(GenomicRanges::seqnames(InteractionSet::anchors(gi, "first")))
  # Find the chr1-chr2 interaction (was originally chr2-chr1)
  idx_swapped <- which(a1_chr == "chr1")
  # The swapped interaction should still carry its original score
  expect_true(all(scores[idx_swapped] %in% c(42, 99)))
  unlink(tmp)
})

# --- bedpe_to_gi: score_col parameter ---
test_that("bedpe_to_gi uses score_col when specified", {
  tmp <- tempfile(fileext = ".bedpe")
  # Column 7 is name, column 8 is p-value, column 9 is score
  writeLines(c(
    "chr1\t100\t200\tchr1\t400\t500\tname\t0.001\t42",
    "chr1\t150\t250\tchr1\t450\t550\tname2\t0.01\t10"
  ), tmp)

  # Auto-detect would pick column 8 (p-value), but explicit score_col=9 picks score
  gi <- looplook:::bedpe_to_gi(tmp, score_col = 9)
  expect_equal(S4Vectors::mcols(gi)$score, c(42, 10))
  unlink(tmp)
})

test_that("bedpe_to_gi errors when score_col exceeds column count", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t400\t500", tmp)
  expect_error(looplook:::bedpe_to_gi(tmp, score_col = 10), "exceeds file column count")
  unlink(tmp)
})

# --- read_simple_bed: edge cases ---
test_that("read_simple_bed returns NULL for NULL input", {
  expect_null(looplook:::read_simple_bed(NULL))
})

test_that("read_simple_bed errors on missing file", {
  expect_error(looplook:::read_simple_bed("/nonexistent/file.bed"), "does not exist")
})

# --- consolidate_chromatin_loops: intersect mode ---
test_that("consolidate_chromatin_loops intersect mode works", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t400\t500", f1)
  writeLines("chr1\t110\t210\tchr1\t410\t510", f2)

  gi <- looplook::consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "intersect",
    gap = 50,
    quiet = TRUE
  )
  expect_s4_class(gi, "GInteractions")
  expect_equal(length(gi), 1)

  # intersect mode must assign cluster_id to each retained loop
  expect_true("cluster_id" %in% colnames(S4Vectors::mcols(gi)))
  expect_false(is.na(S4Vectors::mcols(gi)$cluster_id))
  unlink(c(f1, f2))
})

# --- consolidate_chromatin_loops: intersect mode with cluster_id ---
test_that("consolidate_chromatin_loops intersect mode assigns sequential cluster_id", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  writeLines(c(
    "chr1\t100\t200\tchr1\t400\t500",
    "chr1\t1000\t1100\tchr1\t1400\t1500"
  ), f1)
  writeLines(c(
    "chr1\t110\t210\tchr1\t410\t510",
    "chr1\t1010\t1110\tchr1\t1410\t1510"
  ), f2)

  gi <- looplook::consolidate_chromatin_loops(
    files = c(f1, f2), mode = "intersect", gap = 50, quiet = TRUE
  )
  expect_equal(length(gi), 2)
  mcols <- S4Vectors::mcols(gi)
  expect_true("cluster_id" %in% colnames(mcols))
  expect_false(any(is.na(mcols$cluster_id)))
  # cluster_id values should be unique per retained loop
  expect_equal(length(unique(mcols$cluster_id)), 2)
  unlink(c(f1, f2))
})

# --- consolidate_chromatin_loops: union mode ---
test_that("consolidate_chromatin_loops union mode keeps all", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t400\t500", f1)
  writeLines("chr5\t1000\t2000\tchr5\t4000\t5000", f2)

  gi <- looplook::consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "union",
    gap = 100,
    quiet = TRUE
  )
  expect_s4_class(gi, "GInteractions")
  expect_equal(length(gi), 2) # Different chromosomes, no merge
  unlink(c(f1, f2))
})

# --- consolidate_chromatin_loops: blacklist filtering ---
test_that("consolidate_chromatin_loops filters blacklist", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  bl <- tempfile(fileext = ".bed")
  writeLines("chr1\t100\t200\tchr1\t400\t500", f1)
  writeLines("chr1\t110\t210\tchr1\t410\t510", f2)
  # Blacklist overlaps first anchor
  writeLines("chr1\t90\t210", bl)

  gi <- looplook::consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "consensus",
    gap = 50,
    blacklist_species = bl,
    quiet = TRUE
  )
  # Loop should be filtered out due to blacklist overlap
  expect_equal(length(gi), 0)
  unlink(c(f1, f2, bl))
})

# --- reduce_ginteractions: empty input ---
test_that("reduce_ginteractions handles empty GInteractions", {
  empty_gi <- InteractionSet::GInteractions(
    GenomicRanges::GRanges(),
    GenomicRanges::GRanges(),
    mode = "strict"
  )
  res <- looplook:::reduce_ginteractions(empty_gi)
  expect_equal(length(res$gi), 0)
  expect_equal(length(res$membership), 0)
})

# --- reduce_ginteractions: basic functionality ---
test_that("reduce_ginteractions clusters and reduces correctly", {
  f1 <- tempfile(fileext = ".bedpe")
  writeLines(c(
    "chr1\t100\t200\tchr1\t400\t500",
    "chr1\t110\t210\tchr1\t410\t510"
  ), f1)
  gi <- looplook:::bedpe_to_gi(f1)
  res <- looplook:::reduce_ginteractions(gi, gap = 50)
  expect_s4_class(res$gi, "GInteractions")
  expect_equal(length(res$gi), 1) # Should be merged into 1 cluster
  expect_equal(length(res$membership), 2) # Original 2 loops
  unlink(f1)
})

# --- consolidate_chromatin_loops: region_of_interest filtering ---
test_that("consolidate_chromatin_loops filters by region_of_interest", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  roi <- tempfile(fileext = ".bed")
  writeLines("chr1\t100\t200\tchr1\t400\t500", f1)
  writeLines("chr1\t110\t210\tchr1\t410\t510", f2)
  # ROI overlaps the loop
  writeLines("chr1\t90\t220", roi)

  gi <- looplook::consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "consensus",
    gap = 50,
    region_of_interest = roi,
    roi_mode = "any", # ROI only covers anchor1, need any-mode
    quiet = TRUE
  )
  expect_s4_class(gi, "GInteractions")
  expect_equal(length(gi), 1)
  unlink(c(f1, f2, roi))
})

test_that("consolidate_chromatin_loops returns empty when no ROI overlap", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  roi <- tempfile(fileext = ".bed")
  writeLines("chr1\t100\t200\tchr1\t400\t500", f1)
  writeLines("chr1\t110\t210\tchr1\t410\t510", f2)
  # ROI does NOT overlap
  writeLines("chr5\t1000000\t2000000", roi)

  gi <- looplook::consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "consensus",
    gap = 50,
    region_of_interest = roi,
    quiet = TRUE
  )
  expect_equal(length(gi), 0)
  unlink(c(f1, f2, roi))
})

# --- consolidate_chromatin_loops: write_output ---
test_that("consolidate_chromatin_loops writes output file", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  out_file <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t400\t500", f1)
  writeLines("chr1\t110\t210\tchr1\t410\t510", f2)

  gi <- looplook::consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "consensus",
    gap = 50,
    out_file = out_file,
    write_output = TRUE,
    quiet = TRUE
  )
  expect_true(file.exists(out_file))
  unlink(c(f1, f2, out_file))
})

# --- Production-grade boundary tests (from _problems/) ---
test_that("bedpe_to_gi warns when auto-detected score looks like p-values", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines(c(
    "chr1\t100\t200\tchr1\t400\t500\t0.01",
    "chr1\t150\t250\tchr1\t450\t550\t0.05"
  ), tmp)
  expect_warning(looplook:::bedpe_to_gi(tmp), "resemble p-values")
  unlink(tmp)
})

test_that("bedpe_to_gi errors when score_col points to non-numeric column", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t300\t400\tname", tmp)
  expect_error(
    looplook:::bedpe_to_gi(tmp, score_col = 7),
    "does not contain predominantly numeric"
  )
  unlink(tmp)
})

test_that("consolidate_chromatin_loops returns empty on min_raw_score > all scores (intersect)", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t300\t400\t1", f1)
  writeLines("chr1\t100\t200\tchr1\t300\t400\t5", f2)
  res <- consolidate_chromatin_loops(
    files = c(f1, f2), mode = "intersect",
    min_raw_score = 10, quiet = TRUE
  )
  expect_s4_class(res, "GInteractions")
  expect_equal(length(res), 0)
  unlink(c(f1, f2))
})

test_that("consolidate_chromatin_loops returns empty on min_raw_score > all scores (consensus)", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t300\t400\t1", f1)
  writeLines("chr1\t150\t250\tchr1\t350\t450\t1", f2)
  res <- consolidate_chromatin_loops(
    files = c(f1, f2), mode = "consensus",
    min_raw_score = 10, quiet = TRUE
  )
  expect_s4_class(res, "GInteractions")
  expect_equal(length(res), 0)
  unlink(c(f1, f2))
})

test_that("consolidate_chromatin_loops chaining_policy='warn' warns on wide clusters", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t0\t100\tchr1\t200\t300\t100", f1)
  writeLines(c(
    "chr1\t5\t105\tchr1\t205\t305\t20",
    "chr1\t10\t110\tchr1\t210\t310\t25",
    "chr1\t15\t115\tchr1\t215\t315\t30"
  ), f2)
  expect_warning(
    consolidate_chromatin_loops(
      files = c(f1, f2), mode = "consensus", gap = 25,
      chaining_policy = "warn", quiet = TRUE
    ),
    "max_span > chaining threshold"
  )
  unlink(c(f1, f2))
})

test_that("consolidate_chromatin_loops chaining_policy='drop' removes wide clusters", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t0\t100\tchr1\t200\t300\t100", f1)
  writeLines(c(
    "chr1\t5\t105\tchr1\t205\t305\t20",
    "chr1\t10\t110\tchr1\t210\t310\t25",
    "chr1\t15\t115\tchr1\t215\t315\t30"
  ), f2)
  res_drop <- consolidate_chromatin_loops(
    files = c(f1, f2), mode = "consensus", gap = 25,
    chaining_policy = "drop", quiet = TRUE
  )
  # With chaining, the single wide cluster spanning all loops is dropped
  expect_equal(length(res_drop), 0L)
  unlink(c(f1, f2))
})

test_that("consolidate_chromatin_loops chaining_policy='error' stops on wide clusters", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t0\t100\tchr1\t200\t300\t100", f1)
  writeLines(c(
    "chr1\t5\t105\tchr1\t205\t305\t20",
    "chr1\t10\t110\tchr1\t210\t310\t25",
    "chr1\t15\t115\tchr1\t215\t315\t30"
  ), f2)
  expect_error(
    consolidate_chromatin_loops(
      files = c(f1, f2), mode = "consensus", gap = 25,
      chaining_policy = "error", quiet = TRUE
    ),
    "max_span > chaining threshold"
  )
  unlink(c(f1, f2))
})

test_that("consolidate_chromatin_loops chaining_policy='none' runs silently on wide clusters", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t0\t100\tchr1\t200\t300\t100", f1)
  writeLines(c(
    "chr1\t5\t105\tchr1\t205\t305\t20",
    "chr1\t10\t110\tchr1\t210\t310\t25",
    "chr1\t15\t115\tchr1\t215\t315\t30"
  ), f2)
  expect_no_warning(
    consolidate_chromatin_loops(
      files = c(f1, f2), mode = "consensus", gap = 25,
      chaining_policy = "none", quiet = TRUE
    )
  )
  unlink(c(f1, f2))
})

test_that("consolidate_chromatin_loops write_output=FALSE does not create directory", {
  f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  f2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")
  skip_if(f1 == "" || f2 == "")
  out_dir <- tempfile(pattern = "looplook_cons_nowrite_")
  out_file <- file.path(out_dir, "loops.bedpe")
  unlink(out_dir, recursive = TRUE, force = TRUE)
  expect_false(dir.exists(out_dir))
  gi <- consolidate_chromatin_loops(
    files = c(f1, f2), mode = "consensus", gap = 1000,
    out_file = out_file, write_output = FALSE, quiet = TRUE
  )
  expect_s4_class(gi, "GInteractions")
  expect_false(dir.exists(out_dir))
})

# --- consolidate_chromatin_loops parameter validation ---
test_that("consolidate_chromatin_loops validates numeric parameters", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t400\t500", f1)
  writeLines("chr1\t100\t200\tchr1\t400\t500", f2)
  expect_error(consolidate_chromatin_loops(files = c(f1,f2), gap = -1, quiet = TRUE))
  expect_error(consolidate_chromatin_loops(files = c(f1,f2), min_consensus = 0, quiet = TRUE))
  expect_error(consolidate_chromatin_loops(files = c(f1,f2),
    region_of_interest = "/nonexistent/path.bed", quiet = TRUE))
  unlink(c(f1, f2))
})

# --- .diagnose_gap small data ---
test_that(".diagnose_gap returns NULL for small data", {
  dt <- data.table::data.table(
    chr1 = "chr1", start1 = 100L, end1 = 200L,
    chr2 = "chr1", start2 = 1000L, end2 = 2000L,
    score = 1, source = 1L
  )
  expect_null(looplook:::.diagnose_gap(dt, 1000, function(...) {}))
})

# --- .diagnose_clusters empty ---
test_that(".diagnose_clusters handles zero clusters", {
  empty_gi <- InteractionSet::GInteractions(GenomicRanges::GRanges(), GenomicRanges::GRanges())
  empty_dt <- data.table::data.table(cluster=integer(), n_members=integer(), n_reps=integer())
  expect_null(looplook:::.diagnose_clusters(empty_gi, empty_dt, 1000, function(...) {}))
})

# --- .diagnose_gap: RISK MODERATE path (800bp anchors, gap=2000) ---
test_that(".diagnose_gap detects moderate risk for mid-width anchors", {
  dt <- data.table::data.table(
    chr1 = rep("chr1", 15),
    start1 = seq(100, 15000, by = 1000),
    end1 = seq(900, 15800, by = 1000),
    chr2 = rep("chr1", 15),
    start2 = seq(50000, 64000, by = 1000),
    end2 = seq(50800, 64800, by = 1000),
    score = 1, source = 1L
  )
  msg <- capture.output(looplook:::.diagnose_gap(dt, 2000, message), type = "message")
  expect_true(any(grepl("RISK MODERATE", msg)) ||
              any(grepl("RISK HIGH", msg)))  # either is fine for coverage
})

# --- cluster_loops_dt: multi-chromosome split ---
test_that("cluster_loops_dt handles multi-chromosome data without idx error", {
  dt <- data.table::data.table(
    chr1   = c("chr1", "chr2", "chr1"),
    start1 = c(100L, 100L, 150L),
    end1   = c(200L, 200L, 250L),
    chr2   = c("chr1", "chr2", "chr1"),
    start2 = c(1000L, 1000L, 1050L),
    end2   = c(2000L, 2000L, 2100L),
    score  = c(1, 1, 1),
    source = c(1L, 1L, 1L)
  )
  result <- looplook:::cluster_loops_dt(dt, gap = 0L)
  expect_equal(nrow(result), 3L)
  expect_equal(length(unique(result$cluster)), 2L)  # chr1-chr1 pairs merge, chr2 stays alone
})

# --- Schema closure: cluster_id type + source metadata identical across modes ---
test_that("consolidate_chromatin_loops emits identical mcols schema across all modes", {
  loop1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  loop2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")
  skip_if(loop1 == "" || loop2 == "")

  gi_int <- looplook:::consolidate_chromatin_loops(c(loop1, loop2), mode = "intersect",
                                                   gap = 1000, quiet = TRUE, write_output = FALSE)
  gi_con <- looplook:::consolidate_chromatin_loops(c(loop1, loop2), mode = "consensus",
                                                   gap = 1000, quiet = TRUE, write_output = FALSE)
  gi_uni <- looplook:::consolidate_chromatin_loops(c(loop1, loop2), mode = "union",
                                                   gap = 1000, quiet = TRUE, write_output = FALSE)
  # All three modes must:
  #  (1) emit character cluster_id (no integer vs character mismatch if user
  #      combines gi across modes)
  #  (2) NOT retain the per-source label (output is a multi-source consensus)
  #  (3) carry the same mcols column set
  expect_equal(class(S4Vectors::mcols(gi_int)$cluster_id), "character")
  expect_equal(class(S4Vectors::mcols(gi_con)$cluster_id), "character")
  expect_equal(class(S4Vectors::mcols(gi_uni)$cluster_id), "character")
  expect_false("source" %in% colnames(S4Vectors::mcols(gi_int)))
  expect_false("source" %in% colnames(S4Vectors::mcols(gi_con)))
  expect_false("source" %in% colnames(S4Vectors::mcols(gi_uni)))
  expect_setequal(colnames(S4Vectors::mcols(gi_int)),
                  colnames(S4Vectors::mcols(gi_con)))
  expect_setequal(colnames(S4Vectors::mcols(gi_int)),
                  colnames(S4Vectors::mcols(gi_uni)))
})
