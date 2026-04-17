# ==============================================================================
# 01_synthetic_benchmark.R
# Purpose: Benchmark select_dimension across architectures, metrics,
#          dim reduction methods, and complexity levels on synthetic data.
# Output:  analysis/results/synthetic_benchmark.csv
# ==============================================================================

devtools::load_all(".")

library(oneiric)
library(mclust)

set.seed(42)

# --- Configuration -----------------------------------------------------------

architectures <- list(
    list(name = "3_territories", n_terr = 3, pattern = "circle"),
    list(name = "5_territories", n_terr = 5, pattern = "circle"),
    list(name = "8_territories", n_terr = 8, pattern = "circle")
    # Add "stripes" or "layered" once you verify oneiric supports them
)

metrics        <- c("morans_I", "gearys_C", "spatial_entropy")
dim_methods    <- c("PCA")
# Add "UMAP" once confirmed it runs stably on your Mac
# VAE requires additional setup — add later on HPC

n_cells        <- 1000
n_replicates   <- 5      # multiple replicates for variance estimates
signal_strengths <- c(3, 5, 8)   # fold-change over background

dir.create("analysis/results", recursive = TRUE, showWarnings = FALSE)

# --- Helper: generate one synthetic dataset -----------------------------------

make_synthetic <- function(arch, signal, seed) {
    set.seed(seed)
    layout <- simulate_spatial(n_cells = n_cells,
                               n_territories = arch$n_terr,
                               pattern = arch$pattern)[[1]]
    n_actual      <- nrow(layout)
    n_signal_genes <- arch$n_terr * 10
    n_noise_genes  <- 50
    total_genes    <- n_signal_genes + n_noise_genes

    counts <- matrix(rpois(n_actual * total_genes, lambda = 2),
                     nrow = total_genes, ncol = n_actual)
    colnames(counts) <- layout$barcodes
    rownames(counts) <- paste0("Gene", seq_len(total_genes))

    for (t in seq(0, arch$n_terr - 1)) {
        g_start <- t * 10 + 1
        g_end   <- (t + 1) * 10
        idx     <- which(layout$Territory == t)
        counts[g_start:g_end, idx] <- counts[g_start:g_end, idx] + signal
    }

    coords <- layout[, c("barcodes", "x", "y")]
    list(coords = coords, counts = counts, layout = layout,
         total_genes = total_genes)
}


# --- Main benchmark loop -----------------------------------------------------

results <- data.frame()

for (arch in architectures) {
    for (sig in signal_strengths) {
        for (rep_i in seq_len(n_replicates)) {
            seed <- rep_i * 100 + sig
            dat  <- make_synthetic(arch, signal = sig, seed = seed)

            for (dr in dim_methods) {
                # Build + embed
                ves <- tryCatch({
                    v <- build_vesalius_assay(dat$coords, dat$counts,
                                             verbose = FALSE)
                    v <- generate_embeddings(v, dim_reduction = dr,
                                            nfeatures = dat$total_genes,
                                            verbose = FALSE)
                    v <- smooth_image(v, sigma = 2, iter = 5,
                                     verbose = FALSE)
                    v
                }, error = function(e) {
                    message(sprintf("SKIP %s/%s/sig=%d/rep=%d: %s",
                                    arch$name, dr, sig, rep_i, e$message))
                    NULL
                })
                if (is.null(ves)) next

                dims_to_test <- seq_len(min(arch$n_terr + 5, 15))

                for (met in metrics) {
                    ves_scored <- tryCatch({
                        select_dimension(ves,
                                         dimensions = dims_to_test,
                                         metric     = met,
                                         verbose    = FALSE)
                    }, error = function(e) {
                        message(sprintf("METRIC ERROR %s: %s", met,
                                        e$message))
                        NULL
                    })
                    if (is.null(ves_scored)) next

                    sel_dims <- ves_scored@meta$selected_dims

                    # Clustering accuracy: selected vs all dims
                    embed <- vesalius:::check_embedding_selection(
                        ves, "last", max(dims_to_test))

                    # Baseline: all tested dims
                    km_all <- kmeans(embed[, dims_to_test, drop = FALSE],
                                    centers = arch$n_terr,
                                    nstart = 10, iter.max = 100)
                    ari_all <- adjustedRandIndex(km_all$cluster,
                                                 dat$layout$Territory)

                    # Selected dims only
                    if (length(sel_dims) >= 1) {
                        km_sel <- kmeans(embed[, sel_dims, drop = FALSE],
                                         centers = arch$n_terr,
                                         nstart = 10, iter.max = 100)
                        ari_sel <- adjustedRandIndex(km_sel$cluster,
                                                      dat$layout$Territory)
                    } else {
                        ari_sel <- NA
                    }

                    row <- data.frame(
                        architecture    = arch$name,
                        n_territories   = arch$n_terr,
                        signal_strength = sig,
                        replicate       = rep_i,
                        dim_reduction   = dr,
                        metric          = met,
                        n_dims_tested   = length(dims_to_test),
                        n_dims_selected = length(sel_dims),
                        selected_dims   = paste(sel_dims, collapse = ";"),
                        ari_all_dims    = round(ari_all, 4),
                        ari_selected    = round(ari_sel, 4),
                        improvement     = round(ari_sel - ari_all, 4),
                        stringsAsFactors = FALSE
                    )
                    results <- rbind(results, row)
                }
            }
        }
    }
}

# --- Save results -------------------------------------------------------------
write.csv(results, "analysis/results/synthetic_benchmark.csv",
          row.names = FALSE)

message(sprintf("\nBenchmark complete: %d rows saved to analysis/results/synthetic_benchmark.csv",
                nrow(results)))

# --- Quick summary to console -------------------------------------------------
message("\n=== SUMMARY BY METRIC & DIM REDUCTION ===")
agg <- aggregate(cbind(ari_all_dims, ari_selected, improvement) ~
                     metric + dim_reduction,
                 data = results, FUN = mean, na.rm = TRUE)
print(agg)

message("\n=== SUMMARY BY ARCHITECTURE ===")
agg2 <- aggregate(cbind(ari_all_dims, ari_selected, improvement) ~
                      architecture + signal_strength,
                  data = results, FUN = mean, na.rm = TRUE)
print(agg2)