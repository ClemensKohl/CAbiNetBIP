# Test for biMAP optimization
# This test verifies that the optimized sparse matrix implementation
# produces the same kNN indices and distances as the original implementation

library(testthat)
library(Matrix)

# Original biMAP implementation (for comparison)
run_biMAP_original_dist <- function(
  obj,
  caobj,
  k_knn = min(
    c(Matrix::rowSums(obj@inc), Matrix::colSums(obj@inc)),
    na.rm = TRUE
  )
) {
  n_cells <- nrow(obj@inc)
  n_genes <- ncol(obj@inc)
  n_total <- n_cells + n_genes

  # Create square adjacency matrix: [cells, genes] x [cells, genes]
  square_adj <- Matrix::sparseMatrix(
    i = integer(0),
    j = integer(0),
    x = numeric(0),
    dims = c(n_total, n_total)
  )

  # Set row/column names
  rownames(square_adj) <- c(rownames(obj@inc), colnames(obj@inc))
  colnames(square_adj) <- c(rownames(obj@inc), colnames(obj@inc))

  # Fill in the bipartite connections
  # Upper right: cells to genes
  square_adj[1:n_cells, (n_cells + 1):n_total] <- obj@inc
  # Lower left: genes to cells (transpose)
  square_adj[(n_cells + 1):n_total, 1:n_cells] <- Matrix::t(obj@inc)

  dp <- caobj@std_coords_cols %*% t(caobj@prin_coords_rows)
  dp <- max(dp) - dp
  rownames(dp) <- rownames(caobj@std_coords_cols)
  colnames(dp) <- rownames(caobj@prin_coords_rows)

  # Filter to match incidence matrix dimensions
  dp <- dp[rownames(obj@inc), colnames(obj@inc), drop = FALSE]

  # Create square distance matrix
  square_dist <- rbind(
    cbind(
      matrix(
        Inf, # Distance between cells (not connected)
        ncol = n_cells,
        nrow = n_cells
      ),
      dp # Cell-gene distances
    ),
    cbind(
      t(dp), # Gene-cell distances (transpose)
      matrix(
        Inf, # Distance between genes (not connected)
        ncol = n_genes,
        nrow = n_genes
      )
    )
  )

  rownames(square_dist) <- c(rownames(obj@inc), colnames(obj@inc))
  colnames(square_dist) <- c(rownames(obj@inc), colnames(obj@inc))
  diag(square_dist) <- 0

  # Check k_knn against minimum connections for both cells and genes
  cell_connections <- Matrix::rowSums(obj@inc)
  gene_connections <- Matrix::colSums(obj@inc)
  min_connections <- min(c(cell_connections, gene_connections), na.rm = TRUE)

  if (is.finite(min_connections) && k_knn > min_connections) {
    k_knn <- min_connections
  }

  # Create kNN index and distance matrices for UMAP
  knn_indices <- matrix(0, nrow = n_total, ncol = k_knn)
  knn_distances <- matrix(Inf, nrow = n_total, ncol = k_knn)
  rownames(knn_indices) <- rownames(square_adj)
  rownames(knn_distances) <- rownames(square_adj)

  # For each node, find its k nearest connected neighbors
  for (i in seq_len(n_total)) {
    # Get connected neighbors (non-zero adjacency)
    connected_neighbors <- which(square_adj[i, ] > 0)

    if (length(connected_neighbors) > 0) {
      # Get distances to connected neighbors
      neighbor_distances <- square_dist[i, connected_neighbors]

      # Sort by distance and take top k
      k_to_take <- min(k_knn, length(connected_neighbors))
      sorted_idx <- order(neighbor_distances)[1:k_to_take]

      # Store indices and distances
      knn_indices[i, 1:k_to_take] <- connected_neighbors[sorted_idx]
      knn_distances[i, 1:k_to_take] <- neighbor_distances[sorted_idx]
    }
  }

  return(list(knn_indices = knn_indices, knn_distances = knn_distances))
}

# Optimized implementation (extract the kNN computation part)
run_biMAP_optimized_dist <- function(
  obj,
  caobj,
  k_knn = min(
    c(Matrix::rowSums(obj@inc), Matrix::colSums(obj@inc)),
    na.rm = TRUE
  )
) {
  n_cells <- nrow(obj@inc)
  n_genes <- ncol(obj@inc)
  n_total <- n_cells + n_genes

  # Compute distance matrix between cells and genes
  dp <- caobj@std_coords_cols %*% t(caobj@prin_coords_rows)
  dp <- max(dp) - dp
  rownames(dp) <- rownames(caobj@std_coords_cols)
  colnames(dp) <- rownames(caobj@prin_coords_rows)

  # Filter to match incidence matrix dimensions
  dp <- dp[rownames(obj@inc), colnames(obj@inc), drop = FALSE]

  # Check k_knn against minimum connections for both cells and genes
  cell_connections <- Matrix::rowSums(obj@inc)
  gene_connections <- Matrix::colSums(obj@inc)
  min_connections <- min(c(cell_connections, gene_connections), na.rm = TRUE)

  if (is.finite(min_connections) && k_knn > min_connections) {
    k_knn <- min_connections
  }

  # Create kNN index and distance matrices for UMAP efficiently using sparse ops
  knn_indices <- matrix(0, nrow = n_total, ncol = k_knn)
  knn_distances <- matrix(Inf, nrow = n_total, ncol = k_knn)

  # Set row names
  node_names <- c(rownames(obj@inc), colnames(obj@inc))
  rownames(knn_indices) <- node_names
  rownames(knn_distances) <- node_names

  # Use sparse matrix summary to get all non-zero connections efficiently
  sparse_summary <- Matrix::summary(obj@inc)

  # Only proceed if we have connections and k_knn > 0
  if (nrow(sparse_summary) > 0 && k_knn > 0) {
    # Create a list of connections for each cell
    cell_connections <- split(sparse_summary$j, sparse_summary$i)

    # Process cells using sparse operations
    for (i in seq_len(n_cells)) {
      connected_genes <- cell_connections[[as.character(i)]]

      if (!is.null(connected_genes) && length(connected_genes) > 0) {
        # Get distances to connected genes
        neighbor_distances <- dp[i, connected_genes]

        # Sort by distance and take top k
        k_to_take <- min(k_knn, length(connected_genes))
        sorted_idx <- order(neighbor_distances)[seq_len(k_to_take)]

        # Store indices (offset by n_cells for genes) and distances
        knn_indices[i, seq_len(k_to_take)] <- n_cells +
          connected_genes[sorted_idx]
        knn_distances[i, seq_len(k_to_take)] <- neighbor_distances[sorted_idx]
      }
    }

    # Create a list of connections for each gene (transpose operation)
    gene_connections <- split(sparse_summary$i, sparse_summary$j)

    # Process genes using sparse operations
    for (j in seq_len(n_genes)) {
      connected_cells <- gene_connections[[as.character(j)]]

      if (!is.null(connected_cells) && length(connected_cells) > 0) {
        # Get distances to connected cells (transpose)
        neighbor_distances <- dp[connected_cells, j]

        # Sort by distance and take top k
        k_to_take <- min(k_knn, length(connected_cells))
        sorted_idx <- order(neighbor_distances)[seq_len(k_to_take)]

        # Store indices and distances
        gene_idx <- n_cells + j
        knn_indices[gene_idx, seq_len(k_to_take)] <- connected_cells[sorted_idx]
        knn_distances[gene_idx, seq_len(k_to_take)] <- neighbor_distances[
          sorted_idx
        ]
      }
    }
  }

  return(list(knn_indices = knn_indices, knn_distances = knn_distances))
}

# Helper function to create real test data with caching
create_real_test_data <- function() {
  # Define cache file path
  cache_file <- file.path("data", "test_bimap_data.rds")

  # Load from cache if it exists
  if (file.exists(cache_file)) {
    cat("Loading cached test data...\n")
    return(readRDS(cache_file))
  }

  # Skip test if required packages are not available
  skip_if_not_installed("scRNAseq")
  skip_if_not_installed("scater")
  skip_if_not_installed("APL")

  cat("Creating new test data (this may take a moment)...\n")

  # Load required packages
  suppressPackageStartupMessages({
    library(scRNAseq)
    library(scater)
    library(APL)
  })

  # Get a small dataset - Zeisel brain data
  sce <- ZeiselBrainData()

  # Subset to reasonable size for testing: 300 genes, 200 cells
  set.seed(42)
  selected_genes <- sample(nrow(sce), 300)
  selected_cells <- sample(ncol(sce), 200)
  sce_small <- sce[selected_genes, selected_cells]

  # Basic preprocessing to get logcounts
  sce_small <- logNormCounts(sce_small)

  # Extract logcount matrix
  logcount_matrix <- SingleCellExperiment::logcounts(sce_small)

  # Run correspondence analysis on the logcount matrix
  cacomp_result <- APL::cacomp(
    logcount_matrix,
    princ_coords = 3,
    top = nrow(logcount_matrix),
    dims = 20
  )

  # Run CAbiNetBIP clustering using the same cacomp result
  caclust_result <- run_caclust_bip(
    cacomp_result,
    k = 50,
    min_edges = 1,
    MNN = FALSE,
    algorithm = "leiden",
    leiden_pack = "igraph",
    save_dists = FALSE,
    handle_isolated = "remove"
  )

  # Create result object
  result <- list(caclust = caclust_result, cacomp = cacomp_result)

  # Save to cache
  saveRDS(result, cache_file)
  cat("Test data cached for future use.\n")

  return(result)
}

test_that("biMAP optimization produces same kNN results", {
  # Create toy data
  real_data <- create_real_test_data()
  caclust_obj <- real_data$caclust
  cacomp_obj <- real_data$cacomp

  k_knn <- 10 # Reasonable k for real data

  # Run both implementations
  original_result <- run_biMAP_original_dist(
    caclust_obj,
    cacomp_obj,
    k_knn = k_knn
  )
  optimized_result <- run_biMAP_optimized_dist(
    caclust_obj,
    cacomp_obj,
    k_knn = k_knn
  )

  # Test that dimensions are identical
  expect_equal(
    dim(original_result$knn_indices),
    dim(optimized_result$knn_indices)
  )
  expect_equal(
    dim(original_result$knn_distances),
    dim(optimized_result$knn_distances)
  )

  # Test that row names are identical
  expect_equal(
    rownames(original_result$knn_indices),
    rownames(optimized_result$knn_indices)
  )
  expect_equal(
    rownames(original_result$knn_distances),
    rownames(optimized_result$knn_distances)
  )

  # Test that the kNN indices are identical
  expect_equal(original_result$knn_indices, optimized_result$knn_indices)

  # Test that the kNN distances are identical (with tolerance for floating point)
  expect_equal(
    original_result$knn_distances,
    optimized_result$knn_distances,
    tolerance = 1e-10
  )
})

test_that("optimized biMAP runs without errors on larger example", {
  # Create slightly larger toy data
  set.seed(123)
  n_cells <- 10
  n_genes <- 15

  # Create more complex sparse matrix
  n_connections <- 40
  inc_data <- Matrix::sparseMatrix(
    i = sample(1:n_cells, n_connections, replace = TRUE),
    j = sample(1:n_genes, n_connections, replace = TRUE),
    x = rep(1, n_connections),
    dims = c(n_cells, n_genes)
  )

  rownames(inc_data) <- paste0("cell_", 1:n_cells)
  colnames(inc_data) <- paste0("gene_", 1:n_genes)

  # Create caclust object
  caclust_obj <- new("caclust")
  caclust_obj@inc <- inc_data
  caclust_obj@cell_clusters <- factor(sample(1:3, n_cells, replace = TRUE))
  names(caclust_obj@cell_clusters) <- rownames(inc_data)
  caclust_obj@gene_clusters <- factor(sample(1:4, n_genes, replace = TRUE))
  names(caclust_obj@gene_clusters) <- colnames(inc_data)

  # Create cacomp object
  cacomp_obj <- new("cacomp")
  cacomp_obj@std_coords_cols <- matrix(
    runif(n_cells * 5, -2, 2),
    nrow = n_cells,
    ncol = 5
  )
  cacomp_obj@prin_coords_rows <- matrix(
    runif(n_genes * 5, -2, 2),
    nrow = n_genes,
    ncol = 5
  )
  rownames(cacomp_obj@std_coords_cols) <- rownames(inc_data)
  rownames(cacomp_obj@prin_coords_rows) <- colnames(inc_data)

  k_knn <- 10

  # Test that both implementations run without errors
  expect_no_error(
    original_result <- run_biMAP_original_dist(
      caclust_obj,
      cacomp_obj,
      k_knn = k_knn
    )
  )
  expect_no_error(
    optimized_result <- run_biMAP_optimized_dist(
      caclust_obj,
      cacomp_obj,
      k_knn = k_knn
    )
  )

  # Verify results are identical
  expect_equal(original_result$knn_indices, optimized_result$knn_indices)
  expect_equal(
    original_result$knn_distances,
    optimized_result$knn_distances,
    tolerance = 1e-10
  )
})

test_that("full biMAP pipeline works with optimized implementation", {
  # Test that the actual biMAP function works with the optimization
  real_data <- create_real_test_data()
  caclust_obj <- real_data$caclust
  cacomp_obj <- real_data$cacomp

  k_knn <- 10

  # This should run without any errors now
  expect_no_error({
    result <- run_biMAP(
      caclust_obj,
      cacomp_obj,
      k_knn = k_knn,
      rand_seed = 42,
      method = "dist"
    )
  })

  # Check that biMAP coordinates were generated
  expect_true(nrow(result@bimap) > 0)
  expect_true(all(
    c("x", "y", "name", "type", "cluster") %in% colnames(result@bimap)
  ))

  # Check that we have both cells and genes
  expect_true("cell" %in% result@bimap$type)
  expect_true("gene" %in% result@bimap$type)
})

test_that("biMAP throws error for disconnected nodes", {
  # Create toy data with disconnected nodes
  n_cells <- 3
  n_genes <- 3

  # Create matrix where some nodes have no connections
  inc_data <- Matrix::sparseMatrix(
    i = c(1, 2), # Only cells 1 and 2 have connections, cell 3 is isolated
    j = c(1, 2), # Only genes 1 and 2 have connections, gene 3 is isolated
    x = c(1, 1),
    dims = c(n_cells, n_genes)
  )

  rownames(inc_data) <- paste0("cell_", 1:n_cells)
  colnames(inc_data) <- paste0("gene_", 1:n_genes)

  # Create caclust object
  caclust_obj <- new("caclust")
  caclust_obj@inc <- inc_data
  caclust_obj@cell_clusters <- factor(c(1, 1, 2), levels = c(1, 2))
  names(caclust_obj@cell_clusters) <- rownames(inc_data)
  caclust_obj@gene_clusters <- factor(c(1, 2, 2), levels = c(1, 2))
  names(caclust_obj@gene_clusters) <- colnames(inc_data)

  # Create cacomp object
  cacomp_obj <- new("cacomp")
  cacomp_obj@std_coords_cols <- matrix(
    runif(n_cells * 3, -1, 1),
    nrow = n_cells,
    ncol = 3
  )
  cacomp_obj@prin_coords_rows <- matrix(
    runif(n_genes * 3, -1, 1),
    nrow = n_genes,
    ncol = 3
  )
  rownames(cacomp_obj@std_coords_cols) <- rownames(inc_data)
  rownames(cacomp_obj@prin_coords_rows) <- colnames(inc_data)

  # Should throw an error about disconnected nodes
  expect_error(
    run_biMAP(
      caclust_obj,
      cacomp_obj,
      k_knn = 2,
      rand_seed = 42,
      method = "dist"
    ),
    "Cannot compute biMAP.*no connections"
  )
})

