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
.prepare_inputs <- function(
    sce, use.dimred = "PCA", d = 15,
    positions = NULL, position.cols = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    xdist = NULL, ydist = NULL) {
    inputs <- list()

    if (!(use.dimred %in% reducedDimNames(sce))) {
        stop("reducedDim \"", use.dimred, "\" not found in input SCE.")
    }

    PCs <- reducedDim(sce, use.dimred)
    d <- min(d, ncol(PCs))
    inputs$PCs <- PCs[, seq_len(d)]

    if (is.null(positions)) {
        positions <- as.matrix(colData(sce)[position.cols])
    }

    colnames(positions) <- c("x", "y")
    inputs$positions <- positions

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
