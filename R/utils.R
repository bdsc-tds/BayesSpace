#' Compute pairwise distances between all spots and return list of neighbors
#' for each spot.
#'
#' @param positions (n x 2) matrix of spot coordinates.
#' @param radius The maximum distance for two spots to be considered neighbors.
#' @param method Distance metric to use.
#'
#' @return List df_j, where \code{df_j[[i]]} is a vector of zero-indexed
#'   neighbors of i.
#'
#' @keywords internal
#' @importFrom stats dist
find_neighbors <- function(
    positions, radius,
    method = c("manhattan", "euclidean")) {
  method <- match.arg(method)

  pdist <- as.matrix(stats::dist(positions, method = method))
  neighbors <- (pdist <= radius & pdist > 0)
  df_j <- sapply(
    seq_len(nrow(positions)),
    function(x) as.vector(which(neighbors[x, ])) - 1
  )

  df_j
}

#' Estimate the distance between two neighboring spots
#'
#' Fit linear models between each image pixel coordinate and its corresponding
#' array coordinate to estimate the pixel distance between two spots along
#' each axis. Add these distances to estimate the L1 distance between two
#' spots, then add a small buffer.
#'
#' @param sce SingleCellExperiment (must include array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres
#'   in colData)
#' @param scale.factor Scale estimated L1 difference up by this amount.
#'
#' @return doubles xdist, ydist
#'
#' @keywords internal
#' @importFrom stats lm coef
.compute_interspot_distances <- function(sce) {
  cols <- c("array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  assert_that(all(cols %in% colnames(colData(sce))))

  dists <- list()
  dists$xdist <- coef(lm(sce$pxl_col_in_fullres ~ sce$array_col))[2]
  dists$ydist <- coef(lm(sce$pxl_row_in_fullres ~ sce$array_row))[2]

  dists
}

#' Find the mode
#'
#' Used for finding the most frequent cluster for each z
#'
#' @param x Numeric vector
#'
#' @return mode Numeric scalar, most frequent element in x
#'
#' @keywords internal
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Prepare cluster/deconvolve inputs from SingleCellExperiment object
#'
#' @return List of PCs, names of columns with x/y positions, and inter-spot
#'   distances
#'
#' @keywords internal
#'
#' @importFrom SingleCellExperiment reducedDimNames
#' @importFrom purrr imap
#' @importFrom assertthat assert_that
.prepare_inputs <- function(
    sce, subspots, spot_neighbors, use.dimred = list(PCA = seq_len(15)), use.subspot.dimred = NULL,
    jitter_prior = 0.3, calc.init = TRUE, q = NULL, init = NULL,
    init.method = c("spatialCluster", "mclust", "kmeans"), positions = NULL,
    position.cols = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    xdist = NULL, ydist = NULL, scalef = 1, platform = c("Visium", "VisiumHD", "ST"),
    verbose = FALSE) {
  inputs <- list()

  ## PCs on spot-level (to be enhanced)
  spotPCs <- .check_dimred(sce, use.dimred)[[1]]
  inputs$d2enhance <- sum(vapply(
    names(spotPCs$name),
    function(x) length(spotPCs$name[[x]]),
    FUN.VALUE = integer(1)
  ))

  ## PCs on subspot-level (fixed)
  subspotPCs <- .check_dimred(sce, use.subspot.dimred)
  if (!is.null(subspotPCs)) subspotPCs <- subspotPCs[[1]]

  n.spotPCs <- nrow(spotPCs$PCs)

  ## PCs to be enhanced on subspot-level
  .PCs2enhance <- spotPCs$PCs[rep(seq_len(n.spotPCs), subspots), ]
  rownames(.PCs2enhance) <- paste(
    sce$barcode,
    rep(seq_len(subspots), each = n.spotPCs),
    sep = ":"
  )

  ## PCs to be fixed on subspot-level
  .PCs2fix <- NULL
  if (!is.null(subspotPCs)) .PCs2fix <- subspotPCs$PCs

  ## The amount of jittering (the variance) for the prior distribution
  inputs$c <- jitter_prior * 1 / (2 * mean(diag(cov(spotPCs$PCs))))

  ## pxl coordinates of spots
  if (is.null(positions)) {
    positions <- as.matrix(colData(sce)[position.cols])
  }
  colnames(positions) <- c("x", "y")

  ## pxl coordinates of subspots
  .positions <- positions[rep(seq_len(n.spotPCs), subspots), ]

  ## Compute inter-spot distances (for neighbor finding)
  ## This should only be necessary for Visium enhancement since switching to
  ## array-coordinate-based neighbor finding
  if (is.null(xdist) && is.null(ydist)) {
    dists <- .compute_interspot_distances(sce)
    dists <- imap(dists, function(d, n) ifelse(is.null(get(n)), d, get(n)))

    inputs <- c(inputs, dists)
  } else {
    inputs$xdist <- xdist
    inputs$ydist <- ydist
  }

  ## Compute the pxl coordinates of subspots.
  shift <- .make_subspots(platform, xdist, ydist, scalef)
  shift_long <- shift$shift[rep(seq_len(subspots), each = n0), ]
  positions2[, "x"] <- positions2[, "x"] + shift_long[, "x"]
  positions2[, "y"] <- positions2[, "y"] + shift_long[, "y"]

  inputs$dist <- shift$dist

  # Prepare colData for subspots
  .cdata <- .make_subspot_coldata(as.data.frame(.positions), sce, NULL, platform, 3)

  ## Reorder rows of PCs to be fixed on subspot-level
  if (!is.null(subspotPCs)) .PCs2fix <- .PCs2fix[.cdata$barcode, , drop = FALSE]

  ## Combine PCs to be enhanced and to be fixed
  if (is.null(subspotPCs)) {
    inputs$PCs <- .PCs2enhance
  } else {
    inputs$PCs <- cbind(.PCs2enhance, .PCs2fix)
  }

  ## Initialize cluster assignments (use spatialCluster by default)
  if (calc.init && is.null(init)) {
    if (verbose) {
      message("Initializing clusters...")
    }

    init.method <- match.arg(init.method)
    if (init.method == "spatialCluster") {
      msg <- paste0(
        "Must run spatialCluster on sce before enhancement ",
        "if using spatialCluster to initialize."
      )
      assert_that("spatial.cluster" %in% colnames(colData(sce)), msg = msg)
      init <- sce$spatial.cluster
      inputs$init <- rep(init, subspots)
    } else {
      assert_that(!is.null(q) && q > 0)

      inputs$init <- .init_cluster(inputs$PCs, q, init, init.method)
    }
  }

  ## Create an SCE object for subspot
  colnames(.PCs2enhance) <- vapply(
    strsplit(colnames(.PCs2enhance), "_"),
    function(x) x[length(x)],
    FUN.VALUE = character(1)
  )
  if (!is.null(subspotPCs)) {
    colnames(.PCs2fix) <- vapply(
      strsplit(colnames(.PCs2fix), "_"),
      function(x) x[length(x)],
      FUN.VALUE = character(1)
    )
  }

  inputs$sce <- SingleCellExperiment(
    assays = list(),
    rowData = rowData(sce),
    colData = .cdata,
    reducedDims = compact(list(
      "PCA" = .PCs2enhance,
      "image" = .PCs2fix
    ))
  )

  # Add metadata to new SingleCellExperiment object
  metadata(inputs$sce)$BayesSpace.data <- list()
  metadata(inputs$sce)$BayesSpace.data$platform <- platform
  metadata(inputs$sce)$BayesSpace.data$is.enhanced <- TRUE

  inputs
}

#' Create minimal \code{SingleCellExperiment} for documentation examples.
#'
#' @param nrow Number of rows of spots
#' @param ncol Number of columns of spots
#' @param n_genes Number of genes to simulate
#' @param n_PCs Number of principal components to include
#'
#' @return A SingleCellExperiment object with simulated counts, corresponding
#'   logcounts and PCs, and positional data in \code{colData}. Spots are
#'   distributed over an (\code{nrow} x \code{ncol}) rectangle.
#'
#' @details
#' Inspired by scuttle's \code{mockSCE()}.
#'
#' @examples
#' set.seed(149)
#' sce <- exampleSCE()
#'
#' @importFrom stats rnorm
#' @importFrom SingleCellExperiment SingleCellExperiment logcounts
#' @importFrom scater logNormCounts
#' @importFrom stats prcomp rnbinom runif
#' @importFrom S4Vectors metadata<-
#'
#' @export
exampleSCE <- function(nrow = 8, ncol = 12, n_genes = 100, n_PCs = 10) {
  n_spots <- nrow * ncol

  ## Borrowed from scuttle::mockSCE()
  mu <- 2^runif(n_genes, 2, 10)
  size <- 1 / (100 / mu + 0.5)
  counts <- matrix(rnbinom(n_genes * n_spots, mu = mu, size = size), ncol = n_spots)
  rownames(counts) <- paste0("gene_", seq_len(n_genes))
  colnames(counts) <- paste0("spot_", seq_len(n_spots))

  ## Make array coordinates - filled rectangle
  cdata <- list()
  cdata$array_row <- rep(seq_len(nrow), each = ncol)
  cdata$array_col <- rep(seq_len(ncol), nrow)
  cdata <- as.data.frame(do.call(cbind, cdata))

  ## Scale and jitter image coordinates
  scale.factor <- rnorm(1, 8)
  cdata$pxl_row_in_fullres <- scale.factor * cdata$row + rnorm(n_spots)
  cdata$pxl_col_in_fullres <- scale.factor * cdata$col + rnorm(n_spots)

  ## Make SCE
  ## note: scater::runPCA throws warning on our small sim data, so use prcomp
  sce <- SingleCellExperiment(assays = list(counts = counts), colData = cdata)
  sce <- logNormCounts(sce)
  reducedDim(sce, "PCA") <- prcomp(t(logcounts(sce)))$x[, seq_len(n_PCs)]

  ## Add cluster labels for clusterPlot() examples
  sce$spatial.cluster <- floor(runif(ncol(sce), 1, 4))

  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- "ST"
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

  sce
}

#' Download a processed sample from our S3 bucket
#'
#' Datasets are cached locally using \code{BiocFileCache}. The first time using
#' this function, you may need to consent to creating a BiocFileCache directory
#' if one does not already exist.
#'
#' @param dataset Dataset identifier
#' @param sample Sample identifier
#' @param cache If true, cache the dataset locally with \code{BiocFileCache}.
#'   Otherwise, download directly from our S3 bucket. Caching saves time on
#'   subsequent loads, but consumes disk space.
#'
#' @return sce A SingleCellExperiment with positional information in colData and
#'   PCs based on the top 2000 HVGs
#'
#' @details
#' The following datasets are available via \code{getRDS}.
#' | Dataset  | Sample(s) |
#' | ------------- | ------------- |
#' | 2018_thrane_melanoma | ST_mel1_rep2 |
#' | 2020_maynard_prefrontal-cortex  | 151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673, 151674, 151675, 151676  |
#' | 2020_ji_squamous-cell-carcinoma | P4_rep1 |
#' | 2020_10X-IDC | IDC1 |
#' | 2020_10X-demo_ovarian-cancer | whole_transcriptome |
#'
#' @md
#'
#' @examples
#' sce <- getRDS("2018_thrane_melanoma", "ST_mel1_rep2", cache = FALSE)
#'
#' @export
#' @importFrom RCurl url.exists
#' @importFrom utils download.file
#' @importFrom assertthat assert_that
#' @importFrom BiocFileCache BiocFileCache bfcrpath
getRDS <- function(dataset, sample, cache = TRUE) {
  url <- "https://fh-pi-gottardo-r-eco-public.s3.amazonaws.com/SpatialTranscriptomes/%s/%s.rds"
  url <- sprintf(url, dataset, sample)
  assert_that(url.exists(url), msg = "Dataset/sample not available")

  if (cache) {
    bfc <- BiocFileCache()
    local.path <- bfcrpath(bfc, url)
  } else {
    local.path <- tempfile(fileext = ".rds")
    download.file(url, local.path, quiet = TRUE, mode = "wb")
  }

  ret <- readRDS(local.path)

  # Rename columns of colData of `ret` for compatibility reasons.
  if (any(c("row", "col") %in% colnames(colData(ret)))) {
    col.names <- colnames(colData(ret))
    col.names <- gsub("row", "array_row", col.names)
    col.names <- gsub("col", "array_col", col.names)
    colnames(colData(ret)) <- col.names
  }

  ret
}

#' Access BayesSpace metadata
#'
#' @param sce SingleCellExperiment
#' @param name Metadata name
#'
#' @return Requested metadata
#'
#' @keywords internal
.bsData <- function(sce, name, default = NULL, warn = FALSE) {
  if (!exists("BayesSpace.data", metadata(sce))) {
    return(default)
  }

  bsData <- metadata(sce)[["BayesSpace.data"]]
  if (is.null(name) || is.na(name)) {
    bsData
  } else if (exists(name, bsData)) {
    bsData[[name]]
  } else {
    if (warn) {
      default.name <- ifelse(is.null(default), "NULL", default)
      warning(name, " not found in BayesSpace metadata. Using default: ", default.name)
    }
    default
  }
}

#' Check for reducedDim features
#'
#' @param sce SingleCellExperiment
#' @param name Names of reduced dimensions
#'
#' @return A list of combined reducedDim features
#'
#' @keywords internal
#' @importFrom purrr keep discard compact
.check_dimred <- function(sce, name) {
  ret <- list()

  if (is.character(name)) {
    names(name) <- name
    name <- sapply(
      name,
      function(x) -1,
      simplify = FALSE
    )
  }

  reduced_dim <- compact(sapply(
    names(name),
    function(x) {
      if (x %in% reducedDimNames(sce)) {
        if (is.null(name[[x]]) || (length(name[[x]]) == 1 && name[[x]] == 0)) {
          return(0)
        }

        if (length(name[[x]]) == 1 && name[[x]] == -1) {
          return(seq_len(dim(reducedDim(sce, x))[2]))
        }

        return(keep(name[[x]], function(y) y <= dim(reducedDim(sce, x))[2]))
      }

      NULL
    },
    simplify = FALSE
  ))
  left_over <- discard(names(name), function(x) x %in% names(reduced_dim))
  names(left_over) <- left_over

  metadata <- compact(sapply(
    left_over,
    function(x) {
      if (x %in% names(metadata(sce)[["BayesSpace.data"]])) {
        if (is.null(name[[x]]) || (length(name[[x]]) == 1 && name[[x]] == 0)) {
          return(0)
        }

        if (length(name[[x]]) == 1 && name[[x]] == -1) {
          return(seq_len(dim(metadata(sce)[["BayesSpace.data"]][[x]])[2]))
        }

        return(keep(name[[x]], function(y) y <= dim(metadata(sce)[["BayesSpace.data"]][[x]])[2]))
      }

      NULL
    },
    simplify = FALSE
  ))
  left_over <- discard(left_over, function(x) x %in% names(metadata))

  if (!is.null(reduced_dim) && length(reduced_dim) > 0) {
    ret$reduced_dim <- list(
      name = reduced_dim,
      func = reducedDim
    )
  }

  if (!is.null(metadata) && length(metadata) > 0) {
    ret$metadata <- list(
      name = metadata,
      func = .bsData
    )
  }

  if (!is.null(left_over) && length(left_over) > 0) {
    warning(paste0(
      "Cannot find the following reduced dimensions: ",
      paste0(left_over, collapse = ", ")
    ))
  }

  if (length(ret) > 0) {
    ret <- discard(
      sapply(
        names(ret),
        function(x) {
          c(
            ret[[x]],
            list(
              PCs = do.call(
                cbind,
                lapply(
                  names(ret[[x]]$name),
                  function(y) {
                    if (is.null(ret[[x]]$name[[y]]) || (length(ret[[x]]$name[[y]]) == 1 && ret[[x]]$name[[y]] == 0)) {
                      return(NULL)
                    }

                    .Y <- ret[[x]]$func(sce, y)
                    colnames(.Y) <- paste(y, colnames(.Y), sep = "_")

                    .valid_spots <- rowSums(sapply(
                      colnames(sce),
                      function(z) grepl(paste0("^", z, ".*"), rownames(.Y))
                    ))

                    .Y[rownames(.Y)[.valid_spots > 0], ret[[x]]$name[[y]], drop = FALSE]
                  }
                )
              )
            )
          )
        },
        simplify = FALSE
      ),
      function(x) is.null(x$PCs)
    )
  }

  if (is.null(ret) || length(ret) == 0) {
    ret <- NULL
  }

  ret
}

#' Convert a list into vectors for easier output.
#'
#' @param X A list.
#'
#' @return A vector converted from the input list \code{X}.
#'
#' @keywords internal
#'
#' @importFrom assertthat assert_that
.list2vec <- function(X, sep = "=", collapse = ",", use_names = TRUE) {
  assert_that(!missing(X) && !is.null(X) && is.list(X))

  if (is.null(names(X)) && use_names) {
    names(X) <- X
  }

  if (use_names) {
    ret <- sapply(names(X), function(x) paste(x, X[x], sep = sep), USE.NAMES = FALSE)
  } else {
    ret <- sapply(X, function(x) X[x], USE.NAMES = FALSE)
  }

  if (!is.null(collapse) && length(collapse) > 0) {
    ret <- paste(ret, collapse = collapse)
  }

  ret
}
