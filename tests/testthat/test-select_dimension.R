################################################################################
##                     Unit Tests: select_dimension                            ##
################################################################################

test_that("select_dimension returns a vesalius_assay with expected slots", {
    skip_if_not_installed("vesalius")
    skip_if_not_installed("spdep")
    skip_if_not_installed("oneiric")

    library(vesalius)
    library(oneiric)

    set.seed(42)
    layout <- simulate_spatial(n_cells = 500, n_territories = 3,
                               pattern = "circle")[[1]]
    counts <- matrix(rpois(500 * 30, lambda = 2), nrow = 30, ncol = 500)
    colnames(counts) <- layout$barcodes
    rownames(counts) <- paste0("Gene", 1:30)
    # Inject signal
    for (t in 0:2) {
        g_start <- t * 10 + 1
        g_end   <- (t + 1) * 10
        counts[g_start:g_end, layout$Territory == t] <-
            counts[g_start:g_end, layout$Territory == t] + 5
    }
    coords <- layout[, c("barcodes", "x", "y")]

    ves <- build_vesalius_assay(coords, counts, verbose = FALSE)
    ves <- generate_embeddings(ves, dim_reduction = "PCA",
                               nfeatures = 30, verbose = FALSE)
    ves <- smooth_image(ves, sigma = 2, iter = 3, verbose = FALSE)

    result <- select_dimension(ves, dimensions = 1:5,
                                metric = "morans_I", verbose = FALSE)

    # Returns a vesalius_assay
    expect_s4_class(result, "vesalius_assay")

    # Has the expected meta slots
    expect_true("dim_scores" %in% names(result@meta))
    expect_true("selected_dims" %in% names(result@meta))

    # dim_scores is a data.frame with correct columns
    scores <- result@meta$dim_scores
    expect_s3_class(scores, "data.frame")
    expect_true(all(c("dimension", "score") %in% colnames(scores)))
    expect_equal(nrow(scores), 5L)

    # Scores are sorted descending
    expect_true(all(diff(scores$score) <= 0))

    # selected_dims is an integer vector, subset of tested dims
    sel <- result@meta$selected_dims
    expect_type(sel, "integer")
    expect_true(all(sel %in% 1:5))
    expect_true(length(sel) >= 1)
})


test_that("all three metrics run without error", {
    skip_if_not_installed("vesalius")
    skip_if_not_installed("spdep")
    skip_if_not_installed("entropy")
    skip_if_not_installed("oneiric")

    library(vesalius)
    library(oneiric)

    set.seed(99)
    layout <- simulate_spatial(n_cells = 300, n_territories = 2,
                               pattern = "circle")[[1]]
    counts <- matrix(rpois(300 * 20, lambda = 2), nrow = 20, ncol = 300)
    colnames(counts) <- layout$barcodes
    rownames(counts) <- paste0("G", 1:20)
    counts[1:10, layout$Territory == 0] <-
        counts[1:10, layout$Territory == 0] + 6
    coords <- layout[, c("barcodes", "x", "y")]

    ves <- build_vesalius_assay(coords, counts, verbose = FALSE)
    ves <- generate_embeddings(ves, dim_reduction = "PCA",
                               nfeatures = 20, verbose = FALSE)

    for (m in c("morans_I", "gearys_C", "spatial_entropy")) {
        res <- select_dimension(ves, dimensions = 1:3, metric = m,
                                verbose = FALSE)
        expect_s4_class(res, "vesalius_assay")
        expect_equal(nrow(res@meta$dim_scores), 3L)
    }
})


test_that("invalid inputs are caught gracefully", {
    skip_if_not_installed("vesalius")

    expect_error(select_dimension("not_an_assay"),
                 "vesalius_assay")
})


test_that("dimensions exceeding embedding width produce a warning", {
    skip_if_not_installed("vesalius")
    skip_if_not_installed("spdep")
    skip_if_not_installed("oneiric")

    library(vesalius)
    library(oneiric)

    set.seed(7)
    layout <- simulate_spatial(n_cells = 200, n_territories = 2,
                               pattern = "circle")[[1]]
    counts <- matrix(rpois(200 * 10, lambda = 2), nrow = 10, ncol = 200)
    colnames(counts) <- layout$barcodes
    rownames(counts) <- paste0("G", 1:10)
    coords <- layout[, c("barcodes", "x", "y")]

    ves <- build_vesalius_assay(coords, counts, verbose = FALSE)
    ves <- generate_embeddings(ves, dim_reduction = "PCA",
                               nfeatures = 10, verbose = FALSE)

    # PCA on 10 genes will have <= 10 dims; asking for 1:50 should warn
    expect_warning(
        select_dimension(ves, dimensions = 1:50, verbose = FALSE),
        "dropping"
    )
})


test_that("redundancy filter reduces selected dimensions", {
    skip_if_not_installed("vesalius")
    skip_if_not_installed("spdep")
    skip_if_not_installed("oneiric")

    library(vesalius)
    library(oneiric)

    set.seed(42)
    layout <- simulate_spatial(n_cells = 500, n_territories = 3,
                               pattern = "circle")[[1]]
    counts <- matrix(rpois(500 * 40, lambda = 2), nrow = 40, ncol = 500)
    colnames(counts) <- layout$barcodes
    rownames(counts) <- paste0("Gene", 1:40)
    for (t in 0:2) {
        g_start <- t * 10 + 1
        g_end   <- (t + 1) * 10
        counts[g_start:g_end, layout$Territory == t] <-
            counts[g_start:g_end, layout$Territory == t] + 5
    }
    coords <- layout[, c("barcodes", "x", "y")]

    ves <- build_vesalius_assay(coords, counts, verbose = FALSE)
    ves <- generate_embeddings(ves, dim_reduction = "PCA",
                               nfeatures = 40, verbose = FALSE)
    ves <- smooth_image(ves, sigma = 2, iter = 3, verbose = FALSE)

    res_no_filt <- select_dimension(ves, dimensions = 1:10,
                                     filter_redundant = FALSE,
                                     verbose = FALSE)
    res_filt    <- select_dimension(ves, dimensions = 1:10,
                                     filter_redundant = TRUE,
                                     cor_threshold = 0.8,
                                     verbose = FALSE)

    # Redundancy filter should keep <= as many dims
    expect_true(length(res_filt@meta$selected_dims) <=
                length(res_no_filt@meta$selected_dims))
})