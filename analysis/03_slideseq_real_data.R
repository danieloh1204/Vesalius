# ==============================================================================
# 03_slideseq_real_data.R
# Purpose: Run select_dimension on real Slide-seq hippocampus data.
#          Dev mode uses 3k-bead local subset. For paper, run on HPC
#          with full dataset.
# Output:  analysis/results/slideseq_scores.csv
#          analysis/results/slideseq_selected_dims.rds
# ==============================================================================

devtools::load_all(".")

library(ggplot2)

set.seed(42)

dir.create("analysis/results", recursive = TRUE, showWarnings = FALSE)
dir.create("analysis/figures", recursive = TRUE, showWarnings = FALSE)

# --- Configuration -----------------------------------------------------------
# Toggle when moving to HPC:
DATA_PATH <- "data/slide_seq_3k_macbook.rds"     # local dev
# DATA_PATH <- "data/slide_seq_full.rds"          # HPC full run

DEV_MODE  <- TRUE   # Set FALSE on HPC for full analysis

# --- 1. Load data ------------------------------------------------------------
message("Loading Slide-seq data from: ", DATA_PATH)
slide_data <- readRDS(DATA_PATH)

counts <- slide_data$counts
coords <- slide_data$coords

coords_df <- as.data.frame(coords)
colnames(coords_df)[1:2] <- c("x", "y")
coords_df$barcodes <- rownames(coords)
coords_df <- coords_df[, c("barcodes", "x", "y")]

message(sprintf("Dataset: %d beads, %d genes",
                nrow(coords_df), nrow(counts)))

# --- 2. Process with Vesalius ------------------------------------------------
message("Building Vesalius assay...")
ves <- build_vesalius_assay(coords_df, counts, verbose = FALSE)

# Test multiple dim reduction methods
dim_methods <- c("PCA")
n_features  <- ifelse(DEV_MODE, 300, 2000)
dims_to_test <- if (DEV_MODE) 1:10 else 1:30
metrics     <- c("morans_I", "gearys_C", "spatial_entropy")

all_results <- data.frame()

for (dr in dim_methods) {
    message(sprintf("\n--- Processing: %s ---", dr))

    ves_dr <- tryCatch({
        v <- generate_embeddings(ves, dim_reduction = dr,
                                 nfeatures = n_features, verbose = FALSE)
        smooth_image(v, sigma = 2, iter = 5, verbose = FALSE)
    }, error = function(e) {
        message(sprintf("SKIP %s: %s", dr, e$message))
        NULL
    })
    if (is.null(ves_dr)) next

    for (met in metrics) {
        message(sprintf("  Scoring with %s...", met))

        ves_scored <- tryCatch(
            select_dimension(ves_dr, dimensions = dims_to_test,
                             metric = met, verbose = FALSE),
            error = function(e) {
                message(sprintf("  ERROR %s: %s", met, e$message))
                NULL
            }
        )
        if (is.null(ves_scored)) next

        scores <- ves_scored@meta$dim_scores
        scores$dim_reduction <- dr
        scores$metric        <- met
        scores$selected      <- scores$dimension %in%
                                    ves_scored@meta$selected_dims

        all_results <- rbind(all_results, scores)
    }
}

# --- 3. Save results ---------------------------------------------------------
write.csv(all_results, "analysis/results/slideseq_scores.csv",
          row.names = FALSE)

message(sprintf("\nScoring complete: %d rows saved.", nrow(all_results)))

# --- 4. Summary --------------------------------------------------------------
message("\n=== SELECTED DIMENSIONS BY METHOD & METRIC ===")
sel_summary <- all_results[all_results$selected == TRUE, ]
agg <- aggregate(dimension ~ dim_reduction + metric, data = sel_summary,
                 FUN = function(x) paste(sort(x), collapse = ", "))
colnames(agg)[3] <- "selected_dims"
print(agg)

# --- 5. Consensus across metrics ---------------------------------------------
message("\n=== CONSENSUS: Dimensions selected by ALL metrics ===")
for (dr in unique(all_results$dim_reduction)) {
    sub <- all_results[all_results$dim_reduction == dr &
                       all_results$selected == TRUE, ]
    if (nrow(sub) == 0) next
    freq <- table(sub$dimension)
    consensus <- as.integer(names(freq[freq == length(metrics)]))
    message(sprintf("  %s consensus dims: [%s]", dr,
                    paste(consensus, collapse = ", ")))
}