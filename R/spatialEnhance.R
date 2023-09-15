#' Enhance spot resolution
#'
#' Enhanced clustering of a spatial expression dataset to subspot resolution.
#'
#' @param sce A SingleCellExperiment object containing the spatial data.
#' @param q The number of clusters.
#' @param platform Spatial transcriptomic platform. Specify 'Visium' for hex
#'   lattice geometry or 'ST' for square lattice geometry. Specifying this
#'   parameter is optional when analyzing SingleCellExperiments processed using
#'   \code{\link{readVisium}}, \code{\link{spatialPreprocess}}, or
#'   \code{\link{spatialCluster}}, as this information is included in their
#'   metadata.
#' @param use.dimred A named list with vectors of numbers of top principal
#'   components to use from spot-level data when clustering, named after the
#'   names of several reduced dimensionality results in \code{reducedDims(sce)}
#'   or \code{metadata(sce)$BayesSpace.data}. They must share the same number
#'   of rows and row names. If provided, cluster on these features directly.
#' @param subspot.d A vector of principal components from the corresponding H&E
#'   image to use during the clustering on the subspot-level, if such data is
#'   available (by default \code{seq_len(5)}; set to \code{NULL} or \code{0}
#'   if no such data).
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
#'   is set, it is considered as a percentage.
#' @param jitter_scale Controls the amount of jittering. Small amounts of
#'   jittering are more likely to be accepted but result in exploring the space
#'   more slowly. We suggest tuning \code{jitter_scale} so that Ychange is on
#'   average around 25\%-40\%. Ychange can be accessed via \code{mcmcChain()}.
#'   Alternatively, set it to 0 to activate adaptive MCMC.
#' @param jitter_prior Scale factor for the prior variance, parameterized as the
#'   proportion (default = 0.3) of the mean variance of the PCs.
#'   We suggest making \code{jitter_prior} smaller if the jittered values are
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
#'   via \code{burn.in}. The MCMC chain must be retained.
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

#' Define offsets for each subspot layout.
#'
#' Hex spots are divided into 6 triangular subspots, square spots are divided
#' into 9 squares. Offsets are relative to the spot center.
#'
#' @param n_subspots_per Number of subspots per spot
#' @return Matrix of x and y offsets, one row per subspot
#'
#' @keywords internal
.make_subspot_offsets <- function(n_subspots_per) {
  if (n_subspots_per == 6) {
    rbind(expand.grid(c(1 / 3, -1 / 3), c(1 / 3, -1 / 3)), expand.grid(c(2 / 3, -2 / 3), 0))
    # } else if (n_subspots_per == 7) {
    #     rbind(expand.grid(c(1/3, -1/3), c(1/3, -1/3)), expand.grid(c(2/3, -2/3, 0), 0))
  } else if (n_subspots_per == 9) {
    rbind(expand.grid(c(1 / 3, -1 / 3, 0), c(1 / 3, -1 / 3, 0)))
  } else {
    stop("Only 6 and 9 subspots currently supported.")
  }
}

#' @export
#' @rdname spatialEnhance
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim<-
#' @importFrom SummarizedExperiment rowData
#' @importFrom assertthat assert_that
spatialEnhance <- function(sce, q, platform = c("Visium", "ST"),
                           use.dimred = list(PCA = seq_len(15)), subspot.d = seq_len(5),
                           init = NULL, init.method = c("spatialCluster", "mclust", "kmeans"),
                           model = c("t", "normal"), nrep = 200000, gamma = NULL,
                           mu0 = NULL, lambda0 = NULL, alpha = 1, beta = 0.01,
                           save.chain = FALSE, chain.fname = NULL, burn.in = 10000,
                           jitter_scale = 5, jitter_prior = 0.3, cores = 1, verbose = FALSE) {
  assert_that(nrep >= 100) # require at least one iteration after thinning
  assert_that(burn.in >= 0)
  if (burn.in >= nrep) {
    stop("Please specify a burn-in period shorter than the total number of iterations.")
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
  } else if (platform == "ST") {
    position.cols <- c("array_col", "array_row")
    xdist <- ydist <- 1
  }

  subspots <- ifelse(platform == "Visium", 6, 9)

  ## Prepare for inputs
  inputs <- .prepare_inputs(
    sce,
    subspots = subspots,
    use.dimred = use.dimred,
    use.subspot.dimred = list(subspot_image_feats_pcs = subspot.d),
    jitter_prior = jitter_prior,
    init = init, init.method = init.method,
    positions = NULL, position.cols = position.cols,
    xdist = xdist, ydist = ydist, platform = platform, verbose = verbose
  )

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
    } else if (platform == "ST") {
      gamma <- 2
    }
  }

  deconv <- deconvolve(inputs$PCs, inputs$positions,
    nrep = nrep, gamma = gamma,
    xdist = inputs$xdist, ydist = inputs$ydist, q = q, init = init, model = model,
    platform = platform, verbose = verbose, jitter_scale = jitter_scale,
    jitter_prior = jitter_prior, mu0 = mu0, lambda0 = lambda0, alpha = alpha,
    beta = beta, cores = cores
  )

  ## Create enhanced SCE
  n_subspots_per <- ifelse(platform == "Visium", 6, 9)
  cdata <- .make_subspot_coldata(deconv$positions, sce, n_subspots_per)
  enhanced <- SingleCellExperiment(
    assays = list(),
    rowData = rowData(sce), colData = cdata
  )

  ## Scale burn.in period to thinned intervals, and
  ## add one to skip initialization values stored before first iteration
  burn.in <- (burn.in %/% thin) + 1

  ## Average PCs, excluding burn-in
  deconv_PCs <- Reduce(`+`, deconv$Y[-seq_len(burn.in)]) / (length(deconv$Y) - burn.in)
  colnames(deconv_PCs) <- colnames(reducedDim(inputs$sce, "PCA"))
  reducedDim(inputs$sce, "PCA") <- deconv_PCs

  ## Choose modal cluster label, excluding burn-in
  message(
    "Calculating labels using iterations ", (burn.in - 1) * thin,
    " through ", nrep, "."
  )
  zs <- deconv$z[seq(burn.in, (nrep %/% thin) + 1), ]
  if (burn.in == (nrep %/% thin) + 1) {
    labels <- matrix(zs, nrow = 1)
  } else {
    labels <- apply(zs, 2, Mode)
  }

  inputs$sce$spatial.cluster <- unname(labels)

  if (save.chain) {
    deconv$d2enhance <- inputs$d2enhance
    deconv <- .clean_chain(deconv, method = "enhance")
    params <- c("z", "mu", "lambda", "weights", "Y", "Ychange", "plogLik")
    metadata(inputs$sce)$chain.h5 <- .write_chain(deconv, chain.fname, params)
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
  eff.args <- discard(
    names(args),
    function(x) x %in% c("save.chain", "chain.fname", "cores")
  )

  # Maximum 1000 iterations.
  if (("nrep" %in% names(args) && args["nrep"] > 1000) || !("nrep" %in% names(args))) {
    args["nrep"] <- 1000
    args["burn.in"] <- 100
  }

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
  n_iter <- nrow(zsamples) - 1 # this is technically n_iters / 100

  assert_that(burn.in >= 0)
  if (burn.in < 1) {
    burn.in <- as.integer(n_iter * burn.in)
  }

  zs <- zsamples[seq(burn.in, n_iter + 1), ]
  if (burn.in == n_iter + 1) {
    labels <- matrix(zs, nrow = 1)
  } else {
    labels <- apply(zs, 2, Mode)
  }

  sce$spatial.cluster <- unname(labels)

  sce
}
