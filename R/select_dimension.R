################################################################################
##                        Latent Space Dimension Selection                     ##
################################################################################

#' Select Best Latent Space Dimensions
#'
#' Scores each latent-space dimension by a spatial information metric and
#' returns the subset of dimensions that carry significant spatial structure.
#' Optionally clusters redundant dimensions so that only one representative
#' per spatial pattern is retained.
#'
#' @param vesalius_assay A vesalius_assay object that already contains
#'   embeddings (and ideally has been smoothed).
#' @param dimensions Integer vector indicating which dimensions to evaluate
#'   (e.g. \code{1:30}). Defaults to \code{1:10}.
#' @param metric Character string. One of \code{"morans_I"},
#'   \code{"gearys_C"}, or \code{"spatial_entropy"}.
#' @param k Integer. Number of nearest neighbours used to build the spatial
#'   weights matrix (default 6).
#' @param embedding Character string describing which embedding to use.
#'   Default is \code{"last"}.
#' @param threshold Numeric between 0 and 1. Relative gap size (as a fraction
#'   of the total score range) used to draw the signal/noise cutoff. Larger
#'   values are more conservative. Default 0.15.
#' @param filter_redundant Logical. If \code{TRUE} (default \code{FALSE}),
#'   performs hierarchical clustering on the selected dimensions and keeps one
#'   representative per cluster.
#' @param cor_threshold Numeric. Correlation threshold for the redundancy
#'   filter (default 0.9). Dimensions with absolute pairwise correlation above
#'   this value are grouped together.
#' @param entropy_k Integer. Number of k-means centres used to discretize
#'   continuous values when \code{metric = "spatial_entropy"}. Default 5.
#' @param verbose Logical. Print progress messages? Default \code{TRUE}.
#' @return The input \code{vesalius_assay} with the following additions
#'   accessible via \code{vesalius_assay@meta}:
#'   \describe{
#'     \item{\code{dim_scores}}{Data frame with columns \code{dimension} and
#'       \code{score}, sorted descending.}
#'     \item{\code{selected_dims}}{Integer vector of recommended dimensions.}
#'   }
#' @details
#' For Moran's I and Geary's C the score is computed directly on the
#' continuous embedding values (no discretisation). For spatial entropy the
#' values are discretised via k-means before building a spatial co-occurrence
#' matrix whose Shannon entropy is computed.
#'
#' All scores are oriented so that \strong{higher = more spatial structure}.
#' The automatic cutoff sorts scores in descending order and looks for the
#' largest relative gap; dimensions above that gap are kept.
#'
#' The optional redundancy filter computes pairwise absolute Pearson
#' correlations among the selected dimensions, cuts a hierarchical
#' clustering tree at \code{cor_threshold}, and retains the highest-scoring
#' member from each cluster.
#'
#' @examples
#' \dontrun{
#' ves <- build_vesalius_assay(coords, counts)
#' ves <- generate_embeddings(ves, dim_reduction = "PCA")
#' ves <- smooth_image(ves, sigma = 2, iter = 5)
#' ves <- select_dimension(ves, dimensions = 1:15)
#' }
#' @export
select_dimension <- function(vesalius_assay,
                             dimensions     = 1:10,
                             metric         = "morans_I",
                             k              = 6L,
                             embedding      = "last",
                             threshold      = 0.50,
                             filter_redundant = FALSE,
                             cor_threshold  = 0.9,
                             entropy_k      = 5L,
                             verbose        = TRUE) {

    # ------------------------------------------------------------------
    # 0. Input validation
    # ------------------------------------------------------------------
    if (!inherits(vesalius_assay, "vesalius_assay")) {
        stop("Input must be a vesalius_assay object.")
    }
    metric <- match.arg(metric, c("morans_I", "gearys_C", "spatial_entropy"))
    if (metric == "spatial_entropy" &&
        !requireNamespace("entropy", quietly = TRUE)) {
        stop("Please install the 'entropy' package for spatial_entropy metric.")
    }
    if (!requireNamespace("spdep", quietly = TRUE)) {
        stop("Please install the 'spdep' package.")
    }

    dimensions <- as.integer(dimensions)
    if (any(dimensions < 1)) stop("Dimensions must be positive integers.")

    # ------------------------------------------------------------------
    # 1. Extract coordinates & build spatial weights
    # ------------------------------------------------------------------
    if (verbose) message("Extracting coordinates and building spatial network...")
    coord_df <- vesalius:::get_tiles(vesalius_assay)
    coord_df <- coord_df[coord_df$origin == 1, ]
    coords   <- as.matrix(coord_df[, c("x", "y")])

    # Tiny jitter to avoid spdep duplicate-point errors
    set.seed(42)
    coords <- coords + matrix(runif(length(coords), 0, 1e-5),
                              ncol = 2)

    knn_nb   <- spdep::knn2nb(spdep::knearneigh(coords, k = k))
    listw    <- spdep::nb2listw(knn_nb, style = "W", zero.policy = TRUE)

    # ------------------------------------------------------------------
    # 2. Extract embedding matrix & validate requested dimensions
    # ------------------------------------------------------------------
    max_dim  <- max(dimensions)
    embed    <- tryCatch(
        vesalius:::check_embedding_selection(vesalius_assay,
                                             embedding, max_dim),
        error = function(e) {
            # Patrick's function errors if max_dim exceeds available dims.
            # Fall back to fetching whatever is available.
            vesalius:::check_embedding_selection(vesalius_assay,
                                                 embedding, 1)
        }
    )
    if (max_dim > ncol(embed)) {
        bad <- dimensions[dimensions > ncol(embed)]
        warning(sprintf(
            "Embedding has %d dimensions; dropping requested dims: %s",
            ncol(embed), paste(bad, collapse = ", ")))
        dimensions <- dimensions[dimensions <= ncol(embed)]
        if (length(dimensions) == 0) stop("No valid dimensions remaining.")
    }

    # ------------------------------------------------------------------
    # 3. Score each dimension
    # ------------------------------------------------------------------
    if (verbose) message(sprintf("Scoring %d dimensions with %s...",
                                  length(dimensions), metric))

    scores <- vapply(dimensions, function(d) {
        vals <- embed[, d]
        score <- switch(metric,
            morans_I = {
                res <- spdep::moran.test(vals, listw = listw,
                                         randomisation = FALSE,
                                         zero.policy = TRUE)
                res$estimate["Moran I statistic"]
            },
            gearys_C = {
                res <- spdep::geary.test(vals, listw = listw,
                                         randomisation = FALSE,
                                         zero.policy = TRUE)
                # Geary's C: 0 = perfect positive SA, 1 = random, 2 = negative
                # Invert so higher = more structure
                2 - res$estimate["Geary C statistic"]
            },
            spatial_entropy = {
                # Discretise for entropy calculation
                km  <- suppressWarnings(
                    stats::kmeans(vals, centers = entropy_k,
                                  iter.max = 100, nstart = 10))
                seg <- km$cluster
                # Spatial co-occurrence matrix
                center_cl   <- seg[rep(seq_along(seg),
                                        lengths(knn_nb))]
                neighbor_cl <- seg[unlist(knn_nb)]
                co_occ      <- table(center_cl, neighbor_cl)
                ent         <- entropy::entropy(co_occ, method = "ML")
                # Lower entropy = more structure; negate for consistency
                -ent
            }
        )
        unname(score)
    }, numeric(1))

    result_df <- data.frame(dimension = dimensions, score = scores)
    result_df <- result_df[order(-result_df$score), ]
    rownames(result_df) <- NULL

    # ------------------------------------------------------------------
    # 4. Automatic cutoff via largest gap
    # ------------------------------------------------------------------
    if (verbose) message("Determining cutoff...")
    selected <- .find_cutoff(result_df$score, result_df$dimension, threshold)

    # ------------------------------------------------------------------
    # 5. Optional redundancy filter
    # ------------------------------------------------------------------
    if (filter_redundant && length(selected) > 1) {
        if (verbose) message("Filtering redundant dimensions...")
        selected <- .filter_redundant(embed, selected, result_df,
                                       cor_threshold)
    }

    # ------------------------------------------------------------------
    # 6. Store results back into the vesalius_assay
    # ------------------------------------------------------------------
    if (verbose) {
        message(sprintf("Selected %d dimensions: [%s]",
                        length(selected),
                        paste(sort(selected), collapse = ", ")))
    }

    vesalius_assay@meta$dim_scores    <- result_df
    vesalius_assay@meta$selected_dims <- sort(selected)

    return(vesalius_assay)
}


# ======================================================================
# Internal helpers (not exported)
# ======================================================================

#' Largest-gap cutoff
#'
#' Given scores sorted descending and their matching dimensions, finds the
#' largest gap (relative to the total range) and keeps everything above.
#' @noRd
.find_cutoff <- function(sorted_scores, sorted_dims, threshold) {
    n <- length(sorted_scores)
    if (n <= 2) return(sorted_dims)

    total_range <- sorted_scores[1] - sorted_scores[n]
    if (total_range == 0) return(sorted_dims)       # all identical

    gaps     <- -diff(sorted_scores)                 # positive decreases
    rel_gaps <- gaps / total_range

    # Find the first gap exceeding the threshold
    cut_idx <- which(rel_gaps >= threshold)
    if (length(cut_idx) == 0) {
        # No clear gap — keep all
        return(sorted_dims)
    }

    # Keep everything at or above the gap
    keep_up_to <- cut_idx[1]
    sorted_dims[seq_len(keep_up_to)]
}


#' Redundancy filter via hierarchical clustering on correlation
#' @noRd
.filter_redundant <- function(embed_matrix, selected, score_df,
                               cor_threshold) {
    sub <- embed_matrix[, selected, drop = FALSE]
    cor_mat  <- abs(stats::cor(sub))
    dist_mat <- stats::as.dist(1 - cor_mat)
    hc       <- stats::hclust(dist_mat, method = "complete")
    clusters <- stats::cutree(hc, h = 1 - cor_threshold)

    # From each cluster, keep the dimension with the highest score
    reps <- vapply(unique(clusters), function(cl) {
        members <- selected[clusters == cl]
        sub_df  <- score_df[score_df$dimension %in% members, ]
        sub_df$dimension[which.max(sub_df$score)]
    }, integer(1))

    sort(reps)
}