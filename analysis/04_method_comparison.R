# ==============================================================================
# 04_method_comparison.R
# Purpose: Compare Vesalius + select_dimension against other spatial
#          territory identification tools on synthetic data.
# Tools:   Vesalius (ours), BayesSpace, Giotto, SpaGCN
#          (Install each as needed; script will skip unavailable tools.)
# Output:  analysis/results/method_comparison.csv
# ==============================================================================

devtools::load_all(".")

library(oneiric)

if (!requireNamespace("mcclust", quietly = TRUE)) {
    install.packages("mcclust", repos = "https://cloud.r-project.org")
}
library(mcclust)
library(mcclust)
library(mclust)

set.seed(42)

dir.create("analysis/results", recursive = TRUE, showWarnings = FALSE)

# --- Configuration -----------------------------------------------------------

n_cells        <- 1000
n_replicates   <- 10
architectures  <- list(
    list(name = "simple_3",   n_terr = 3, pattern = "circle"),
    list(name = "complex_6",  n_terr = 6, pattern = "circle"),
    list(name = "complex_8",  n_terr = 8, pattern = "circle")
)

# --- Helper: generate synthetic data -----------------------------------------

make_data <- function(arch, seed) {
    set.seed(seed)
    layout <- simulate_spatial(n_cells = n_cells,
                               n_territories = arch$n_terr,
                               pattern = arch$pattern)[[1]]
    n_actual    <- nrow(layout)
    n_sig_genes <- arch$n_terr * 10
    n_noise     <- 50
    total       <- n_sig_genes + n_noise

    counts <- matrix(rpois(n_actual * total, lambda = 2),
                     nrow = total, ncol = n_actual)
    colnames(counts) <- layout$barcodes
    rownames(counts) <- paste0("Gene", seq_len(total))

    for (t in seq(0, arch$n_terr - 1)) {
        g_start <- t * 10 + 1
        g_end   <- (t + 1) * 10
        idx     <- which(layout$Territory == t)
        counts[g_start:g_end, idx] <- counts[g_start:g_end, idx] + 5
    }

    coords <- layout[, c("barcodes", "x", "y")]
    list(coords = coords, counts = counts, layout = layout,
         total_genes = total)
}


# --- Method runners (each returns a cluster vector) --------------------------

run_vesalius_baseline <- function(dat, n_terr) {
    # Vesalius with ALL dims (no selection)
    ves <- build_vesalius_assay(dat$coords, dat$counts, verbose = FALSE)
    ves <- generate_embeddings(ves, dim_reduction = "PCA",
                               nfeatures = dat$total_genes, verbose = FALSE)
    ves <- smooth_image(ves, sigma = 2, iter = 5, verbose = FALSE)
    embed <- vesalius:::check_embedding_selection(ves, "last", 10)
    km <- kmeans(embed, centers = n_terr, nstart = 10, iter.max = 100)
    km$cluster
}


run_vesalius_selected <- function(dat, n_terr) {
    # Vesalius + our dimension selection
    ves <- build_vesalius_assay(dat$coords, dat$counts, verbose = FALSE)
    ves <- generate_embeddings(ves, dim_reduction = "PCA",
                               nfeatures = dat$total_genes, verbose = FALSE)
    ves <- smooth_image(ves, sigma = 2, iter = 5, verbose = FALSE)
    ves <- select_dimension(ves, dimensions = 1:min(n_terr + 5, 15),
                            metric = "morans_I", verbose = FALSE)
    sel <- ves@meta$selected_dims
    embed <- vesalius:::check_embedding_selection(ves, "last",
                                                  max(sel))
    km <- kmeans(embed[, sel, drop = FALSE], centers = n_terr,
                 nstart = 10, iter.max = 100)
    km$cluster
}


run_vesalius_selected_redundancy <- function(dat, n_terr) {
    # Vesalius + dimension selection + redundancy filter
    ves <- build_vesalius_assay(dat$coords, dat$counts, verbose = FALSE)
    ves <- generate_embeddings(ves, dim_reduction = "PCA",
                               nfeatures = dat$total_genes, verbose = FALSE)
    ves <- smooth_image(ves, sigma = 2, iter = 5, verbose = FALSE)
    ves <- select_dimension(ves, dimensions = 1:min(n_terr + 5, 15),
                            metric = "morans_I",
                            filter_redundant = TRUE,
                            cor_threshold = 0.85,
                            verbose = FALSE)
    sel <- ves@meta$selected_dims
    embed <- vesalius:::check_embedding_selection(ves, "last",
                                                  max(sel))
    km <- kmeans(embed[, sel, drop = FALSE], centers = n_terr,
                 nstart = 10, iter.max = 100)
    km$cluster
}


# Add wrappers for BayesSpace, SpaGCN, etc. here as they become available.
# Each should take (dat, n_terr) and return a cluster vector of length
# nrow(dat$layout).
#
# run_bayesspace <- function(dat, n_terr) { ... }
# run_spagcn     <- function(dat, n_terr) { ... }


# --- Registry of methods to test ---------------------------------------------

methods <- list(
    list(name = "Vesalius_all_dims",           fn = run_vesalius_baseline),
    list(name = "Vesalius_selected",           fn = run_vesalius_selected),
    list(name = "Vesalius_selected_redundancy", fn = run_vesalius_selected_redundancy)
    # list(name = "BayesSpace",                fn = run_bayesspace),
    # list(name = "SpaGCN",                    fn = run_spagcn)
)


# --- Main comparison loop ----------------------------------------------------

results <- data.frame()

for (arch in architectures) {
    message(sprintf("\n=== Architecture: %s (%d territories) ===",
                    arch$name, arch$n_terr))

    for (rep_i in seq_len(n_replicates)) {
        dat <- make_data(arch, seed = rep_i * 100 + arch$n_terr)
        truth <- dat$layout$Territory

        for (m in methods) {
            clusters <- tryCatch(
                m$fn(dat, arch$n_terr),
                error = function(e) {
                    message(sprintf("  SKIP %s rep %d: %s",
                                    m$name, rep_i, e$message))
                    NULL
                }
            )
            if (is.null(clusters)) next

            ari <- adjustedRandIndex(clusters, truth)
            vi  <- vi.dist(clusters, truth)

            row <- data.frame(
                architecture = arch$name,
                n_territories = arch$n_terr,
                replicate    = rep_i,
                method       = m$name,
                ARI          = round(ari, 4),
                VI           = round(vi, 4),
                stringsAsFactors = FALSE
            )
            results <- rbind(results, row)
        }
    }
}

# --- Save results -------------------------------------------------------------
write.csv(results, "analysis/results/method_comparison.csv",
          row.names = FALSE)

message(sprintf("\nMethod comparison complete: %d rows saved.", nrow(results)))

# --- Summary ------------------------------------------------------------------
message("\n=== MEAN ARI BY METHOD & ARCHITECTURE ===")
agg <- aggregate(ARI ~ method + architecture, data = results,
                 FUN = function(x) sprintf("%.3f (sd=%.3f)",
                                            mean(x), sd(x)))
print(agg)

message("\n=== OVERALL MEAN ARI BY METHOD ===")
agg2 <- aggregate(ARI ~ method, data = results,
                  FUN = function(x) sprintf("%.3f (sd=%.3f)",
                                             mean(x), sd(x)))
print(agg2)