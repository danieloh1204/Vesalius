# ==============================================================================
# 06_consensus_metrics.R
# Purpose: Test whether combining multiple spatial metrics via consensus
#          produces more robust dimension selection than any single metric.
# Strategies: intersection (all agree), union (any says keep),
#             majority vote (2 of 3 agree)
# Output:  analysis/results/consensus_metrics.csv
# ==============================================================================

devtools::load_all(".")

library(oneiric)
library(mclust)

set.seed(42)

dir.create("analysis/results", recursive = TRUE, showWarnings = FALSE)

# --- Configuration -----------------------------------------------------------

architectures <- list(
    list(name = "3_territories", n_terr = 3, pattern = "circle"),
    list(name = "5_territories", n_terr = 5, pattern = "circle"),
    list(name = "8_territories", n_terr = 8, pattern = "circle")
)

metrics      <- c("morans_I", "gearys_C", "spatial_entropy")
n_cells      <- 1000
n_replicates <- 5
signal       <- 5

# --- Helper ------------------------------------------------------------------

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
        counts[g_start:g_end, idx] <- counts[g_start:g_end, idx] + signal
    }

    coords <- layout[, c("barcodes", "x", "y")]
    list(coords = coords, counts = counts, layout = layout,
         total_genes = total)
}


compute_ari <- function(embed, dims, n_terr, truth) {
    if (length(dims) == 0) return(NA)
    km <- kmeans(embed[, dims, drop = FALSE], centers = n_terr,
                 nstart = 10, iter.max = 100)
    adjustedRandIndex(km$cluster, truth)
}

# --- Main loop ----------------------------------------------------------------

results <- data.frame()

for (arch in architectures) {
    message(sprintf("\n=== %s ===", arch$name))

    for (rep_i in seq_len(n_replicates)) {
        dat <- make_data(arch, seed = rep_i * 100 + arch$n_terr)

        ves <- tryCatch({
            v <- build_vesalius_assay(dat$coords, dat$counts, verbose = FALSE)
            v <- generate_embeddings(v, dim_reduction = "PCA",
                                     nfeatures = dat$total_genes,
                                     verbose = FALSE)
            smooth_image(v, sigma = 2, iter = 5, verbose = FALSE)
        }, error = function(e) NULL)
        if (is.null(ves)) next

        dims_to_test <- seq_len(min(arch$n_terr + 5, 15))
        embed <- vesalius:::check_embedding_selection(ves, "last",
                                                       max(dims_to_test))
        truth <- dat$layout$Territory

        # Baseline
        ari_baseline <- compute_ari(embed, dims_to_test, arch$n_terr, truth)

        # Get selected dims from each metric
        selected_by_metric <- list()
        for (met in metrics) {
            ves_scored <- tryCatch(
                select_dimension(ves, dimensions = dims_to_test,
                                 metric = met, verbose = FALSE),
                error = function(e) NULL
            )
            if (!is.null(ves_scored)) {
                selected_by_metric[[met]] <- ves_scored@meta$selected_dims
            }
        }

        if (length(selected_by_metric) < 3) next

        # --- Build consensus strategies ---

        all_dims_selected <- unlist(selected_by_metric)
        dim_counts <- table(all_dims_selected)

        # Intersection: dimension must be selected by ALL metrics
        intersection_dims <- as.integer(
            names(dim_counts[dim_counts == length(metrics)]))

        # Union: dimension selected by ANY metric
        union_dims <- as.integer(names(dim_counts))

        # Majority: dimension selected by >= 2 metrics
        majority_dims <- as.integer(
            names(dim_counts[dim_counts >= 2]))

        # --- Compute ARI for each strategy ---

        strategies <- list(
            list(name = "baseline_all",     dims = dims_to_test),
            list(name = "morans_I_only",    dims = selected_by_metric[["morans_I"]]),
            list(name = "gearys_C_only",    dims = selected_by_metric[["gearys_C"]]),
            list(name = "entropy_only",     dims = selected_by_metric[["spatial_entropy"]]),
            list(name = "intersection",     dims = intersection_dims),
            list(name = "union",            dims = union_dims),
            list(name = "majority_vote",    dims = majority_dims)
        )

        for (s in strategies) {
            valid_dims <- s$dims[s$dims <= ncol(embed)]
            ari <- compute_ari(embed, valid_dims, arch$n_terr, truth)

            row <- data.frame(
                architecture  = arch$name,
                n_territories = arch$n_terr,
                replicate     = rep_i,
                strategy      = s$name,
                n_dims        = length(valid_dims),
                dims_used     = paste(sort(valid_dims), collapse = ";"),
                ARI           = round(ari, 4),
                stringsAsFactors = FALSE
            )
            results <- rbind(results, row)
        }
    }
}

# --- Save results -------------------------------------------------------------
write.csv(results, "analysis/results/consensus_metrics.csv",
          row.names = FALSE)

message(sprintf("\nConsensus analysis complete: %d rows saved.", nrow(results)))

# --- Summary ------------------------------------------------------------------
message("\n=== MEAN ARI BY STRATEGY & ARCHITECTURE ===")
agg <- aggregate(cbind(ARI, n_dims) ~ strategy + architecture,
                 data = results, FUN = mean, na.rm = TRUE)
agg <- agg[order(agg$architecture, -agg$ARI), ]
print(agg, row.names = FALSE)

message("\n=== OVERALL MEAN ARI BY STRATEGY ===")
overall <- aggregate(ARI ~ strategy, data = results,
                     FUN = function(x) sprintf("%.3f (sd=%.3f)",
                                                mean(x, na.rm = TRUE),
                                                sd(x, na.rm = TRUE)))
overall <- overall[order(-as.numeric(gsub(" .*", "",  overall$ARI))), ]
print(overall, row.names = FALSE)