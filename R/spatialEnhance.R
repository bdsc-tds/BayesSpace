#' Enhance spot resolution
#'
#' Enhanced clustering of a spatial expression dataset to subspot resolution.
#'
#' @param sce A SingleCellExperiment object containing the spatial data.
#' @param q The number of clusters.
#' @param platform Spatial transcriptomic platform. Specify 'Visium' for hex
#'   lattice geometry or 'ST' and 'VisiumHD' for square lattice geometry.
#'   Specifying this parameter is optional when analyzing SingleCellExperiments
#'   processed using \code{\link{readVisium}}, \code{\link{spatialPreprocess}},
#'   or \code{\link{spatialCluster}}, as this information is included in their
#'   metadata.
#' @param use.dimred Name of a reduced dimensionality result in
#'   \code{reducedDims(sce)}. If provided, cluster on these features directly.
#' @param d Number of top principal components to use when clustering.
#' @param nsubspots.per.edge Number of subspots per edge of the square. Only
#'   valid when \code{platform} is 'ST' or 'VisiumHD'.
#' @param init Initial cluster assignments for spots.
#' @param init.method If \code{init} is not provided, cluster the top \code{d}
#'   PCs with this method to obtain initial cluster assignments.
#' @param model Error model. ('normal' or 't')
#' @param nrep The number of MCMC iterations.
#' @param gamma Smoothing parameter. (Values in range of 1-3 seem to work well.)
#' @param mu0 Prior mean hyperparameter for mu. If not provided, mu0 is set to
#'   the mean of PCs over all spots.
#' @param lambda0 Prior precision hyperparam for mu. If not provided, lambda0
#'   is set to a diagonal matrix \eqn{0.01 I}.
#' @param alpha Hyperparameter for Wishart distributed precision lambda.
#' @param beta Hyperparameter for Wishart distributed precision lambda.
#' @param save.chain If true, save the MCMC chain to an HDF5 file.
#' @param chain.fname File path for saved chain. Tempfile used if not provided.
#' @param burn.in Number of iterations to exclude as burn-in period. The MCMC
#'   iterations are currently thinned to every 100; accordingly \code{burn.in}
#'   is rounded down to the nearest multiple of 100. If a value no larger than 1
#'   is set, it is considered as a percentage. It is always considered as
#'   percentage for \code{adjustClusterLabels}.
#' @param jitter.scale Controls the amount of jittering. Small amounts of
#'   jittering are more likely to be accepted but result in exploring the space
#'   more slowly. We suggest tuning \code{jitter.scale} so that Ychange is on
#'   average around 25\%-40\%. Ychange can be accessed via \code{mcmcChain()}.
#'   Alternatively, set it to 0 to activate adaptive MCMC.
#' @param adapt.before Adapting the MCMC chain before the specified number
#'   or proportion of iterations (by default equal to \code{burn.in}; set to 0
#'   to always adapt). Only valid when \code{jitter.scale} is 0.
#' @param jitter.prior Scale factor for the prior variance, parameterized as the
#'   proportion (default = 0.3) of the mean variance of the PCs.
#'   We suggest making \code{jitter.prior} smaller if the jittered values are
#'   not expected to vary much from the overall mean of the spot.
#' @param cores The number of threads to use. The results are invariate to the
#'   value of \code{cores}.
#' @param verbose Log progress to stderr.
#' @param test.cores Either a list of, or a maximum number of cores to test. In
#'   the latter case, a list of values (power of 2) will be created
#' @param test.times Times to repeat the benchmarking with microbenchmark.
#' @param ... Arguments for \code{spatialEnhance} (except for cores).
#'
#' @return
#' \code{spatialEnhance} returns a new SingleCellExperiment object.
#'   By default, the \code{assays} of this object are empty, and the enhanced
#'   resolution PCs are stored as a reduced dimensionality result accessible
#'   with \code{reducedDim(sce, 'PCA')}.
#'
#' \code{coresTune} returns the output of \code{microbenchmark}.
#'
#' \code{adjustClusterLabels} adjusts the cluster labels from the MCMC samples
#'   via \code{burn.in}, the percentage of samples to drop. The MCMC chain
#'   must be retained.
#'
#' @details
#' The enhanced \code{SingleCellExperiment} has most of the properties of the
#'   input SCE - \code{rowData}, \code{colData}, \code{reducedDims} - but does
#'   not include expression data in \code{counts} or \code{logcounts}. To impute
#'   enhanced expression vectors, please use [enhanceFeatures()] after
#'   running \code{spatialEnhance}.
#'
#' The \code{colData} of the enhanced \code{SingleCellExperiment} includes the
#'   following columns to permit referencing the subspots in spatial context and
#'   linking back to the original spots:
#'   \itemize{
#'   \item \code{spot.idx}: Index of the spot this subspot belongs to (with
#'     respect to the input SCE).
#'   \item \code{subspot.idx}: Index of the subspot within its parent spot.
#'   \item \code{spot.row}: Array row of the subspot's parent spot.
#'   \item \code{spot.col}: Array col of the subspot's parent spot.
#'   \item \code{array_row}: Array row of the subspot. This is the parent spot's row
#'     plus an offset based on the subspot's position within the spot.
#'   \item \code{array_col}: Array col of the subspot. This is the parent spot's col
#'     plus an offset based on the subspot's position within the spot.
#'   \item \code{pxl_row_in_fullres}: Pixel row of the subspot. This is the parent spot's
#'     row plus an offset based on the subspot's position within the spot.
#'   \item \code{pxl_col_in_fullres}: Pixel col of the subspot. This is the parent spot's
#'     col plus an offset based on the subspot's position within the spot.
#'   }
#'
#' @examples
#' set.seed(149)
#' sce <- exampleSCE()
#' sce <- spatialCluster(sce, 7, nrep = 100, burn.in = 10)
#' enhanced <- spatialEnhance(sce, 7, nrep = 100, burn.in = 10)
#'
#' @seealso \code{\link{spatialCluster}} for clustering at the spot level
#'   before enhancing, \code{\link{clusterPlot}} for visualizing the cluster
#'   assignments, \code{\link{enhanceFeatures}} for imputing enhanced
#'   expression, and \code{\link{mcmcChain}} for examining the full MCMC chain
#'   associated with the enhanced clustering.
#'   .
#'
#' @name spatialEnhance
NULL

#' Wrapper around C++ \code{iterate_deconv()} function
#'
#' @return List of enhancement parameter values at each iteration
#'
#' @keywords internal
#' @importFrom stats cov
#' @importFrom purrr discard
deconvolve <- function(Y, positions, xdist, ydist, scalef, q, spot_neighbors, init, nrep = 1000, thin = 100,
                       model = "normal", platform = c("Visium", "VisiumHD", "ST"), nsubspots.per.edge = 3, verbose = TRUE,
                       jitter.scale = 5, jitter.prior = 0.01, adapt.before = 100, mu0 = colMeans(Y), gamma = 2,
                       lambda0 = diag(0.01, nrow = ncol(Y)), alpha = 1, beta = 0.01, cores = 1) {
  d <- ncol(Y)
  n0 <- nrow(Y)
  Y <- as.matrix(Y)
  c <- jitter.prior * 1 / (2 * mean(diag(cov(Y))))

  positions <- as.matrix(positions)
  colnames(positions) <- c("x", "y")

  platform <- match.arg(platform)
  subspots <- ifelse(platform == "Visium", 6, nsubspots.per.edge^2)

  init1 <- rep(init, subspots)
  Y2 <- Y[rep(seq_len(n0), subspots), ]
  positions2 <- positions[rep(seq_len(n0), subspots), ]
  
  shift <- .make_subspots(platform, xdist, ydist, nsubspots.per.edge = nsubspots.per.edge)
  shift_long <- shift$shift[rep(seq_len(subspots), each = n0), ]
  positions2[, "x"] <- positions2[, "x"] + shift_long[, "x"]
  positions2[, "y"] <- positions2[, "y"] + shift_long[, "y"]
  n <- nrow(Y2)

  if (verbose) {
    message("Fitting model...")
  }
  tdist <- (model == "t")
  out <- iterate_deconv(
    subspot_positions = positions2,
    dist = as.numeric(shift$dist),
    spot_neighbors = spot_neighbors,
    Y = Y2, tdist = tdist, nrep = nrep, thin = thin, n = n, n0 = n0,
    d = d, gamma = gamma, q = q, init = init1, subspots = subspots, verbose = verbose,
    jitter_scale = jitter.scale, adapt_before = adapt.before, c = c, mu0 = mu0,
    lambda0 = lambda0, alpha = alpha, beta = beta, thread_num = cores
  )

  # The indices of neighbors are 1-based.
  out$df_j <- apply(
    out$df_j,
    1,
    function(x) paste(sort(discard(x, function(y) y == 0)), collapse = ",")
  )

  out$positions <- positions2
  out
}

#' Define offsets and Manhattan distances for each subspot layout.
#'
#' Hex spots are divided into 6 triangular subspots, square spots are divided
#' into 9 squares. Offsets are relative to the spot center. A unit corresponds
#' to the diameter of a spot.
#'
#' Manhattan distance is used here instead of Euclidean to avoid numerical
#' issues.
#'
#' @param platform The platform from which the data comes.
#' @param scalef Scale factors of Visium data.
#' @return Matrix of x and y offsets, one row per subspot
#'
#' @keywords internal
#'
#' @importFrom assertthat assert_that
.make_subspots <- function(
    platform, xdist, ydist, force = FALSE, nsubspots.per.edge = 3, tolerance = 1.05
) {
  if (platform == "Visium") {
    if (abs(xdist) >= abs(ydist) && !force) {
      stop("Unable to find neighbors of subspots. Please raise an issue to maintainers.")
    }

    shift <- rbind(
      expand.grid(
        list(
          x = c(1 / 3, -1 / 3),
          y = c(1 / 3, -1 / 3)
        )
      ),
      expand.grid(
        list(
          x = c(2 / 3, -2 / 3),
          y = 0
        )
      )
    )
  } else if (platform %in% c("VisiumHD", "ST")) {
    vec <- .make_square_vec(nsubspots.per.edge, tolerance)

    shift <- expand.grid(
      list(
        x = vec$vec,
        y = vec$vec
      )
    )

    dist <- vec$dist
  } else {
    stop("Only data from Visium, VisiumHD and ST currently supported.")
  }

  shift[, c("x", "y")] <- as.data.frame(t(
    t(as.matrix(shift[, c("x", "y")])) * c(xdist, ydist)
  ))

  if (platform == "Visium") {
    dist <- max(rowSums(abs(shift))) * tolerance
  }

  list(
    shift = shift,
    dist = dist
  )
}

#' @keywords internal
.make_square_vec <- function(nsubspots.per.edge, tolerance = 1.05) {
  stopifnot(tolerance > 1)

  list(
    vec = seq(
      1 / (2 * nsubspots.per.edge),
      1,
      by = 1 / nsubspots.per.edge
    ) - 1 / 2,
    dist = tolerance / nsubspots.per.edge
  )
}

#' #' Define offsets for each subspot layout.
#' #'
#' #' Hex spots are divided into 6 triangular subspots, square spots are divided
#' #' into 9 squares. Offsets are relative to the spot center.
#' #'
#' #' @param n_subspots_per Number of subspots per spot
#' #' @return Matrix of x and y offsets, one row per subspot
#' #'
#' #' @keywords internal
#' .make_subspot_offsets <- function(n_subspots_per) {
#'   if (n_subspots_per == 6) {
#'     rbind(
#'       expand.grid(
#'         list(
#'           x = c(1 / 3, -1 / 3),
#'           y = c(1 / 3, -1 / 3)
#'         )
#'       ),
#'       expand.grid(
#'         list(
#'           x = c(2 / 3, -2 / 3),
#'           y = 0
#'         )
#'       )
#'     )
#'     # } else if (n_subspots_per == 7) {
#'     #     rbind(expand.grid(c(1/3, -1/3), c(1/3, -1/3)), expand.grid(c(2/3, -2/3, 0), 0))
#'   } else if (n_subspots_per == 9) {
#'     rbind(expand.grid(c(1 / 3, -1 / 3, 0), c(1 / 3, -1 / 3, 0)))
#'   } else {
#'     stop("Only 6 and 9 subspots currently supported.")
#'   }
#' }

#' Add subspot labels and offset row/col locations before making enhanced SCE.
#'
#' Subspots are stored as (1.1, 2.1, 3.1, ..., 1.2, 2.2, 3.2, ...)
#'
#' @param cdata Table of colData (imagerow and imagecol; from deconv$positions)
#' @param sce Original sce (to obtain number of spots and original row/col)
#' @param n_subspots_per Number of subspots per spot
#'
#' @return Data frame with added subspot names, parent spot indices, and offset
#'   row/column coordinates
#'
#' @keywords internal
#' @importFrom assertthat assert_that
.make_subspot_coldata <- function(
    cdata, sce, subspot_neighbors, platform, nsubspots.per.edge = 3
) {
  if (platform == "Visium") {
    n_subspots_per <- 6
    colnames(cdata) <- c("pxl_col_in_fullres", "pxl_row_in_fullres")
  } else if (platform %in% c("VisiumHD", "ST")) {
    n_subspots_per <- nsubspots.per.edge ^ 2
    colnames(cdata) <- c("array_col", "array_row")
  } else {
    stop("Only Visium, VisumHD, and ST are supported.")
  }

  n_spots <- ncol(sce)
  n_subspots <- nrow(cdata)
  assert_that(nrow(cdata) == n_spots * n_subspots_per)

  ## Index of parent spot is (subspot % n_spots)
  idxs <- seq_len(n_subspots)
  spot_idxs <- ((idxs - 1) %% n_spots) + 1
  subspot_idxs <- rep(seq_len(n_subspots_per), each = n_spots)
  cdata$spot.idx <- spot_idxs
  cdata$spot.neighbors <- sce$spot.neighbors[spot_idxs]
  cdata$subspot.idx <- subspot_idxs
  cdata$subspot.neighbors <- subspot_neighbors
  rownames(cdata) <- paste0("subspot_", spot_idxs, ".", subspot_idxs)

  # w.r.t. the coordinate system of array
  array_offsets <- .make_subspots(platform, 1, 1, TRUE, nsubspots.per.edge)$shift
  
  cdata$spot.row <- rep(sce$array_row, n_subspots_per)
  cdata$spot.col <- rep(sce$array_col, n_subspots_per)
  
  if (platform == "Visium") {
    cdata$array_col <- cdata$spot.col + rep(array_offsets[, "x"], each = n_spots)
    cdata$array_row <- cdata$spot.row + rep(array_offsets[, "y"], each = n_spots)
  }
  
  cols <- c(
    "spot.idx", "spot.neighbors",
    "subspot.idx", "subspot.neighbors",
    "spot.row", "spot.col", "array_row", "array_col"
  )
  
  if (platform %in% c("Visium", "VisiumHD")) {
    cdata$spot.pxl.row <- rep(sce$pxl_row_in_fullres, n_subspots_per)
    cdata$spot.pxl.col <- rep(sce$pxl_col_in_fullres, n_subspots_per)
    
    if (platform == "VisiumHD") {
      dist <- .compute_interspot_distances(sce)
      
      # w.r.t. the coordinate system of the full resolution image
      pxl_offsets <- .make_subspots(platform, dist$xdist, dist$ydist, nsubspots.per.edge = nsubspots.per.edge)$shift
      
      cdata$pxl_col_in_fullres <- cdata$spot.pxl.col + rep(pxl_offsets[, "x"], each = n_spots)
      cdata$pxl_row_in_fullres <- cdata$spot.pxl.row + rep(pxl_offsets[, "y"], each = n_spots)
    }
    
    cols <- c(
      cols,
      "spot.pxl.row", "spot.pxl.col",
      "pxl_row_in_fullres", "pxl_col_in_fullres"
    )
  }

  cdata[, cols]
}

#' @export
#' @rdname spatialEnhance
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim<-
#' @importFrom SummarizedExperiment rowData
#' @importFrom assertthat assert_that
spatialEnhance <- function(sce, q, platform = c("Visium", "VisiumHD", "ST"),
                           use.dimred = "PCA", d = 15, nsubspots.per.edge = 3,
                           init = NULL, init.method = c("spatialCluster", "mclust", "kmeans"),
                           model = c("t", "normal"), nrep = 100000, gamma = NULL,
                           mu0 = NULL, lambda0 = NULL, alpha = 1, beta = 0.01,
                           save.chain = FALSE, chain.fname = NULL, burn.in = 10000,
                           jitter.scale = 5, jitter.prior = 0.3, adapt.before = burn.in, cores = 1,
                           verbose = FALSE) {
  assert_that(nrep >= 100) # require at least one iteration after thinning
  assert_that(burn.in >= 0)
  if (burn.in >= nrep) {
    stop("Please specify a burn-in period shorter than the total number of iterations.")
  } else if (burn.in < 1) {
    burn.in <- as.integer(nrep * burn.in)
  }

  if (jitter.scale == 0) {
    assert_that(adapt.before >= 0)

    if (adapt.before >= nrep) {
      stop("Please specify a period for adaptive MCMC shorter than the total number of iterations.")
    } else if (adapt.before < 1) {
      adapt.before <- as.integer(nrep * adapt.before)
    }
  }

  ## Thinning interval; only every 100 iterations are kept to reduce memory
  ## This is temporarily hard-coded into the C++ code
  thin <- 100

  ## If user didn't specify a platform, attempt to parse from SCE metadata
  ## otherwise check against valid options
  if (length(platform) > 1) {
    platform <- .bsData(sce, "platform", match.arg(platform))
  } else {
    platform <- match.arg(platform)
  }

  if (platform == "Visium") {
    position.cols <- c("pxl_col_in_fullres", "pxl_row_in_fullres")
    xdist <- ydist <- NULL # Compute with .prepare_inputs
  } else if (platform %in% c("VisiumHD", "ST")) {
    position.cols <- c("array_col", "array_row")
    xdist <- ydist <- 1
  }

  inputs <- .prepare_inputs(sce,
    use.dimred = use.dimred, d = d,
    positions = NULL, position.cols = position.cols,
    xdist = xdist, ydist = ydist
  )

  ## Initialize cluster assignments (use spatialCluster by default)
  if (is.null(init)) {
    init.method <- match.arg(init.method)
    if (init.method == "spatialCluster") {
      msg <- paste0(
        "Must run spatialCluster on sce before enhancement ",
        "if using spatialCluster to initialize."
      )
      assert_that("spatial.cluster" %in% colnames(colData(sce)), msg = msg)
      init <- sce$spatial.cluster
    } else {
      init <- .init_cluster(inputs$PCs, q, init, init.method)
    }
  }

  ## Set model parameters
  model <- match.arg(model)
  if (is.null(mu0)) {
    mu0 <- colMeans(inputs$PCs)
  }
  if (is.null(lambda0)) {
    lambda0 <- diag(0.01, ncol(inputs$PCs))
  }
  if (is.null(gamma)) {
    if (platform == "Visium") {
      gamma <- 3
    } else if (platform %in% c("VisiumHD", "ST")) {
      gamma <- 2
    }
  }

  deconv <- deconvolve(inputs$PCs, inputs$positions,
    nrep = nrep, gamma = gamma,
    xdist = inputs$xdist, ydist = inputs$ydist, scalef = .bsData(sce, "scalef"),
    q = q, spot_neighbors = sce$spot.neighbors, init = init, model = model,
    platform = platform, nsubspots.per.edge = nsubspots.per.edge,
    verbose = verbose, jitter.scale = jitter.scale, jitter.prior = jitter.prior,
    adapt.before = adapt.before, mu0 = mu0, lambda0 = lambda0, alpha = alpha,
    beta = beta, cores = cores
  )

  ## Create enhanced SCE
  cdata <- .make_subspot_coldata(
    as.data.frame(deconv$positions), sce, deconv$df_j, platform, nsubspots.per.edge
  )
  enhanced <- SingleCellExperiment(
    assays = list(),
    rowData = rowData(sce), colData = cdata
  )

  ## Scale burn.in period to thinned intervals
  .burn.in <- burn.in %/% thin
  .nrep <- nrep %/% thin

  ## Average PCs, excluding burn-in
  deconv_PCs <- Reduce(`+`, deconv$Y[-seq_len(.burn.in + 1)]) / (.nrep - .burn.in)
  colnames(deconv_PCs) <- paste0("PC", seq_len(ncol(deconv_PCs)))
  reducedDim(enhanced, "PCA") <- deconv_PCs

  ## Choose modal cluster label, excluding burn-in
  message(
    "Calculating labels using iterations ", burn.in + 1,
    " through ", nrep, "."
  )
  zs <- deconv$z[seq(.burn.in + 2, .nrep + 1), ]
  if (.burn.in + 1 == .nrep) {
    labels <- matrix(zs, nrow = 1)
  } # if only one iteration kept, return it
  else {
    labels <- apply(zs, 2, Mode)
  } # else take modal assignment

  enhanced$spatial.cluster <- unname(labels)

  if (save.chain) {
    deconv <- .clean_chain(deconv, method = "enhance")
    params <- c("z", "mu", "lambda", "weights", "Y", "Ychange", "plogLik")
    metadata(enhanced)$chain.h5 <- .write_chain(deconv, chain.fname, params)
  }

  ## Add metadata to new SingleCellExperiment object
  metadata(enhanced)$BayesSpace.data <- list()
  metadata(enhanced)$BayesSpace.data$platform <- platform
  metadata(enhanced)$BayesSpace.data$is.enhanced <- TRUE

  enhanced
}

#' @export
#' @rdname spatialEnhance
#' @importFrom microbenchmark microbenchmark
#' @importFrom parallel detectCores
#' @importFrom purrr compact discard
coreTune <- function(sce, test.cores = detectCores(), test.times = 1, ...) {
  assert_that(
    length(test.cores) == length(unique(test.cores)),
    msg = paste0("Duplicate values found in 'test.cores'.")
  )

  args <- list(...)

  # Maximum 1000 iterations.
  if (("nrep" %in% names(args) && args["nrep"] > 1000) || !("nrep" %in% names(args))) {
    args["nrep"] <- 1000
    args["burn.in"] <- 100
  }
  
  eff.args <- discard(
    names(args),
    function(x) x %in% c("save.chain", "chain.fname", "cores")
  )
  
  if (length(test.cores) == 1) {
    cores <- as.integer(vapply(
      seq(
        from = 0,
        to = floor(log(test.cores, 2)),
        by = 1
      ),
      function(x) 2^x,
      FUN.VALUE = numeric(1),
      USE.NAMES = FALSE
    ))
  } else {
    cores <- test.cores
  }

  if (is.null(names(cores))) {
    names(cores) <- paste("c", cores, sep = "_")
  }

  message(paste0(
    "Number of cores to test: ",
    paste0(cores, collapse = ",")
  ))

  exprs <- sapply(
    cores,
    function(x) {
      bquote(do.call(
        spatialEnhance,
        c(
          sce = sce,
          args[eff.args],
          list(
            save.chain = FALSE,
            cores = .(x)
          )
        )
      ))
    },
    simplify = FALSE
  )

  microbenchmark(
    list = exprs,
    times = test.times
  )
}

#' @export
#' @rdname spatialEnhance
adjustClusterLabels <- function(sce, burn.in) {
  zsamples <- mcmcChain(sce, "z")

  assert_that(
    burn.in >= 0,
    burn.in < 1
  )
  burn.in <- as.integer((nrow(zsamples) - 1) * burn.in)

  zs <- zsamples[seq(burn.in + 2, nrow(zsamples)), ]
  if (burn.in + 2 == nrow(zsamples)) {
    labels <- matrix(zs, nrow = 1)
  } else {
    labels <- apply(zs, 2, Mode)
  }

  sce$spatial.cluster <- unname(labels)

  sce
}

#' @export
#' @rdname spatialEnhance
mapSubspot2Ref <- function(sce, sce.ref, cores = 1) {
  map_subspot2ref(
    as.matrix(colData(sce)[c("pxl_col_in_fullres", "pxl_row_in_fullres")]),
    as.matrix(colData(sce.ref)[c("pxl_col_in_fullres", "pxl_row_in_fullres")]),
    cores
  )
}

#' @export
#' @rdname spatialEnhance
computeCorr <- function(m1, m2, cores = 1) {
  compute_corr(
    m1,
    m2,
    cores
  )
}
