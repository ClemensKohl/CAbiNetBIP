#' @include classes.R
NULL

#' Run UMAP embedding for cell-gene graph built up by caclust.
#'
#' @description
#' This function takes cacomp and caclust object as input to calculate
#' UMAP embedding of cell-gene graph in several different ways:
#' * 'dist'(Default): run UMAP on the distance matrix of cell-gene
#' bipartite graph built up by caclust.
#' * 'spectral': run UMAP on the selected eigenvectors of cell-gene graph
#' laplacian (only eligiable when algorithm is set as 'spectral' in 'caclust'
#' function)
#' * 'ca': run UMAP on the singular vectors from Correspondence Analysis.
#' @rdname run_biMAP
#' @param obj results from biclustering of class "caclust"
#' @param caobj A cacomp object with principal and standard coordinates
#' calculated. Only needs to be supplied when using method "ca".
#' @param k_umap integer. Number of nearest neighbours to use to compute UMAP.
#' @param k_knn integer. Number of nearest neighbours to use use from the
#' bipartite graph for construciton of the pre-computed kNN graph for UMAP.
#' @param rand_seed integer. Random seed for UMAP.
#' @param method Can be either "dist", "spectral" or "ca". When using "dist"
#' or ca", a "cacomp" object has to be provided for `caobj`.
#' @param use_inc TRUE/FALSE. This parameter only works when method == 'ca'.
#' If TRUE, it will calculate bimap embedding of genes and cells contained in
#' the bipartite graph together with the given 'features' if any.
#' If FALSE, only calculate bimap embedding of cells contained in the bipartite
#' graph together with the given 'fearures' if any.
#' @param features character/vector of gene names.
#' This parameter only works when method == 'ca'.
#' 'ca' method allows visualization of genes which are absent from biparitte
#' graph, so the 'features' could be genes detected by our biclustering
#' algorithm, and genes which goes through CA analysis.
#' Setting use_inc = FALSE allows users to have a visualization of the
#' features genes user defined/interested.
#'
#' @return
#' caclust object with biMAP coordinates stored in the `bimap` slot.
#'
#' @md
run_biMAP <- function(
  obj,
  caobj = NULL,
  k_umap = 30,
  k_knn = min(Matrix::rowSums(obj@inc)),
  rand_seed = 2358,
  method = c("dist", "spectral", "ca"),
  use_inc = TRUE,
  features = NULL
) {
  stopifnot(is(obj, "caclust"))
  method <- match.arg(method)

  cellc <- names(cell_clusters(obj))
  genec <- names(gene_clusters(obj))

  if (method == "dist") {
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

    if (k_knn > min(Matrix::rowSums(obj@inc))) {
      k_knn <- min(Matrix::rowSums(obj@inc))
      cat(
        "Reducing knn_k to ",
        k_knn,
        " as it is larger than the minimum number of neighbours in the incidence matrix.\n"
      )
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

    bip_knn_umap <- umap::umap.knn(
      indexes = knn_indices,
      distances = knn_distances
    )

    custom_config <- umap::umap.defaults
    custom_config$random_state <- rand_seed

    input_data <- rbind(
      caobj@std_coords_cols[rownames(obj@inc), , drop = FALSE],
      caobj@prin_coords_rows[colnames(obj@inc), , drop = FALSE]
    )

    sce_umap <- umap::umap(
      input_data,
      config = custom_config,
      n_neighbors = k_umap,
      knn = bip_knn_umap
    )

    umap_coords <- as.data.frame(sce_umap$layout)
    rownames(umap_coords) <- rownames(input_data)
  } else if (method == "spectral") {
    eigen <- get_eigen(obj)

    if (is.empty(eigen)) {
      stop("Spectral clustering not run.")
    }
    custom_config <- umap::umap.defaults
    custom_config$random_state <- rand_seed

    caclust_umap <- umap::umap(
      eigen,
      config = custom_config,
      n_neighbors = k_umap,
      metric = "cosine"
    )

    umap_coords <- as.data.frame(caclust_umap$layout)
  } else if (method == "ca") {
    stopifnot(!is.null(caobj))
    stopifnot(is(caobj, "cacomp"))

    if (isTRUE(use_inc)) {
      inc <- get_inc(obj)
      selected_items <- c(rownames(inc), colnames(inc))
    } else {
      selected_items <- c(rownames(caobj@V), rownames(caobj@U))
      cellc <- rownames(caobj@V)
      genec <- rownames(caobj@U)
    }

    if (!is.null(features)) {
      ix <- features %in% rownames(caobj@U)
      if (sum(ix) < length(features)) {
        warning(
          "only ",
          sum(ix),
          " out of ",
          length(features),
          " features are found from the caobj."
        )
      }

      selected_items <- c(selected_items, features[ix])
      selected_items <- unique(selected_items)

      genec <- c(genec, features[ix])
      genec <- unique(genec)
    }

    # eigen = rbind(caobj@V, caobj@U)
    # eigen <- rbind(caobj@std_coords_cols, caobj@prin_coords_rows)
    eigen <- rbind(caobj@prin_coords_cols, caobj@std_coords_rows)

    custom.config <- umap::umap.defaults
    custom.config$random_state <- rand_seed

    eigen <- eigen[rownames(eigen) %in% selected_items, ]

    caclust_umap <- umap::umap(
      eigen,
      config = custom.config,
      metric = "cosine",
      n_neighbors = k_umap
    )
    umap_coords <- as.data.frame(caclust_umap$layout)
  } else {
    stop()
  }

  colnames(umap_coords) <- c("x", "y")
  umap_coords$name <- rownames(umap_coords)

  umap_coords$type <- "none"
  umap_coords$type[umap_coords$name %in% cellc] <- "cell"
  umap_coords$type[umap_coords$name %in% genec] <- "gene"

  umap_coords$cluster <- "none"

  cell_idx <- stats::na.omit(match(names(cell_clusters(obj)), umap_coords$name))
  gene_idx <- stats::na.omit(match(names(gene_clusters(obj)), umap_coords$name))

  umap_coords$cluster[cell_idx] <- cell_clusters(obj)
  umap_coords$cluster[gene_idx] <- gene_clusters(obj)
  umap_coords$cluster <- factor(
    umap_coords$cluster,
    levels = sort(as.numeric(unique(umap_coords$cluster)), decreasing = FALSE)
  )

  umap_coords <- umap_coords %>% dplyr::arrange(desc(type))

  obj@bimap <- umap_coords
  return(obj)
}


# #' Add cacomp obj results to SingleCellExperiment object
# #' @param sce SingleCellExperiment object
# #' @param umap_coords data.frame with coordinates of genes and cells
# #' @param  biMAP_meta_name name not listed in colData(sce), rowData(sce), or metadata(sce)
# #' @export
# #'
# add_biMAP_sce <- function(sce, umap_coords, biMAP_meta_name = 'biMAP'){
#
#   S4Vectors::metadata(sce)[[biMAP_meta_name]] <- umap_coords
#
#   return(sce)
# }

#' Compute biMAP
#'
#' @description
#' The function takes either a `caclust` or `SingleCellExperiment`as input and
#' stores the biMAP in the "bimap" slot in the caclust object. If a
#' SingleCellExperiment was provided the caclust object is stored in its
#' metadata.
#' @name biMAP
#' @rdname biMAP
#' @param obj A caclust object or SingleCellExperiment object
#' @param method Can be either "dist" or "spectral".
#' @inheritParams run_biMAP
#' @param ... Further arguments
#' @details
#' The biMAP cell and gene embeddings can be calculated via different methods
#' as controlled by the parameter `method`:
#' * 'dist'(Default): run UMAP on the distance matrix of cell-gene SNN graph
#' built up by caclust, which is '1-adj(SNN)'.
#' * 'spectral': run UMAP on the selected eigenvectors of cell-gene graph
#' laplacian (only eligiable when algorithm is set as 'spectral' in 'caclust'
#' function)
#' @return
#' A caclust object or SingleCellExperiment object.
#'
#' @md
#' @export
setGeneric(
  "biMAP",
  function(
    obj,
    caobj,
    k_umap = 30,
    k_knn = min(Matrix::rowSums(obj@inc)),
    rand_seed = 2358,
    method = "dist",
    use_inc = TRUE,
    features = NULL,
    ...
  ) {
    standardGeneric("biMAP")
  }
)


#' @rdname biMAP
#' @export
setMethod(
  f = "biMAP",
  signature(obj = "caclust"),
  function(
    obj,
    caobj,
    k_umap = 30,
    k_knn = min(Matrix::rowSums(obj@inc)),
    rand_seed = 2358,
    method = "dist",
    use_inc = TRUE,
    features = NULL,
    ...
  ) {
    stopifnot(method %in% c("dist", "spectral"))

    obj <- run_biMAP(
      obj = obj,
      caobj = caobj,
      k_umap = k_umap,
      k_knn = k_knn,
      rand_seed = rand_seed,
      method = method,
      use_inc = use_inc,
      features = features
    )
    return(obj)
  }
)

#' @rdname biMAP
#' @param caclust_meta_name The name of caclust object stored in metadata
#' (SingleCellExperiment object).
#' @export
setMethod(
  f = "biMAP",
  signature(obj = "SingleCellExperiment"),
  function(
    obj,
    caobj,
    k_umap = 30,
    k_knn = min(Matrix::rowSums(obj@inc)),
    rand_seed = 2358,
    method = "dist",
    use_inc = TRUE,
    features = NULL,
    ...,
    caclust_meta_name = "caclust"
  ) {
    stopifnot(method %in% c("dist", "spectral"))

    if (isFALSE(caclust_meta_name %in% names(S4Vectors::metadata(obj)))) {
      stop(
        "The caclust_meta_name in not found in metadata(sce obj), change meta_name"
      )
    }

    caclust_obj <- S4Vectors::metadata(obj)[[caclust_meta_name]]

    caclust_obj <- run_biMAP(
      obj = caclust_obj,
      caobj = caobj,
      k_umap = k_umap,
      k_knn = k_knn,
      rand_seed = rand_seed,
      method = method,
      use_inc = use_inc,
      features = features
    )

    S4Vectors::metadata(obj)[[caclust_meta_name]] <- caclust_obj
    # TODO
    # allow adding multi-bimap coordinate slots to caclust with slot
    # names 'biMAP_'+algorithm, eg. 'biMAP_SNNdist'

    return(obj)
  }
)


#' Compute biMAP basedon CA results
#'
#' @description
#' The function takes either a `caclust` and a `cacomp` object as input and
#' computes the biMAP embedings for cells and genes on the basis of the singular
#' vectors of CA.
#' @name ca_biMAP
#' @rdname ca_biMAP
#'
#' @param obj A `caclust` object or `SingleCellExperiment` object with `caclust`
#' and `cacomp` objects stored in the metadata.
#' @param caobj A `cacomp` object.
#' @inheritParams run_biMAP
#' @param ... Further arguments
#' @details
#' The biMAP embeddings are computed on the basis of the singular vectors
#' from Correspondence Analysis.
#' If a `SingleCellExperiment` object with `caclust` and `cacomp`
#' objects stored is provided the argument `caobj` is not required.
#'
#' @return
#' an caclust object or SingleCellExperiment objects
#' @export
setGeneric(
  "ca_biMAP",
  function(
    obj,
    caobj,
    k_umap = 30,
    rand_seed = 2358,
    use_inc = TRUE,
    features = NULL,
    ...
  ) {
    standardGeneric("ca_biMAP")
  }
)


#' @rdname ca_biMAP
#' @export
setMethod(
  f = "ca_biMAP",
  signature(obj = "caclust", caobj = "cacomp"),
  function(
    obj,
    caobj,
    k_umap = 30,
    rand_seed = 2358,
    use_inc = TRUE,
    features = NULL,
    ...
  ) {
    obj <- run_biMAP(
      obj = obj,
      caobj = caobj,
      k_umap = k_umap,
      rand_seed = rand_seed,
      method = "ca",
      use_inc = use_inc,
      features = features
    )
    return(obj)
  }
)


#' @rdname ca_biMAP
#' @param cacomp_meta_name the name of cacomp object stored in
#' metadata(obj)
#' @param caclust_meta_name the name of caclust object stored in metadata(obj)
#' @export
setMethod(
  f = "ca_biMAP",
  signature(obj = "SingleCellExperiment"),
  function(
    obj,
    caobj = NULL,
    k_umap = 30,
    rand_seed = 2358,
    use_inc = TRUE,
    features = NULL,
    ...,
    caclust_meta_name = "caclust",
    cacomp_meta_name = "CA"
  ) {
    correct <- check_caobj_sce(obj, cacomp_meta_name = cacomp_meta_name)

    if (isFALSE(correct)) {
      stop(
        "No 'CA' dimension reduction object found. ",
        "Please run cacomp(sce_obj, top, coords = FALSE, ",
        "return_input=TRUE) first."
      )
    }

    caobj <- APL::as.cacomp(obj)

    caclust_obj <- S4Vectors::metadata(obj)[[caclust_meta_name]]

    caclust_obj <- run_biMAP(
      obj = caclust_obj,
      caobj = caobj,
      k_umap = k_umap,
      rand_seed = rand_seed,
      method = "ca",
      use_inc = use_inc,
      features = features
    )

    S4Vectors::metadata(obj)[[caclust_meta_name]] <- caclust_obj
    # TODO
    # allow adding multi-bimap coordinate slots to caclust with slot
    # names 'biMAP_'+algorithm, eg. 'biMAP_SNNdist'

    return(obj)
  }
)
