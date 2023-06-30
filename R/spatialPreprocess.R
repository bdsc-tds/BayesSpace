#' Preprocess a spatial dataset for BayesSpace
#'
#' Adds metadata required for downstream analyses, and (optionally) performs PCA
#' on log-normalized expression of top HVGs.
#'
#' @param sce SingleCellExperiment to preprocess
#' @param platform Spatial sequencing platform. Used to determine spot layout
#'   and neighborhood structure (Visium = hex, ST = square).
#' @param n.PCs Number of principal components to compute. We suggest using the
#'   top 15 PCs in most cases.
#' @param n.HVGs Number of highly variable genes to run PCA upon.
#' @param skip.PCA Skip PCA (if dimensionality reduction was previously
#'   computed.)
#' @param log.normalize Whether to log-normalize the input data with scater. May
#'   be omitted if log-normalization previously computed.
#' @param assay.type Name of assay in \code{sce} containing normalized counts.
#'   Leave as "logcounts" unless you explicitly pre-computed a different
#'   normalization and added it to \code{sce} under another assay. Note that we
#'   do not recommend running BayesSpace on PCs computed from raw counts.
#' @param h2o.max.mem The maximum amount of memory that h2o uses.
#' @param n.PCs.image The number of PCs to be extracted from the image (if 
#'   available).
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying which
#'   algorithm should be used to perform the PCA. By default, an exact PCA is
#'   performed, as current spatial datasets are generally small (<10,000 spots).
#'   To perform a faster approximate PCA, please specify
#'   \code{FastAutoParam()} and set a random seed to ensure
#'   reproducibility.
#'
#' @return SingleCellExperiment with PCA and BayesSpace metadata
#'
#' @examples
#' sce <- exampleSCE()
#' sce <- spatialPreprocess(sce)
#'
#' @export
#' @importFrom scater logNormCounts runPCA
#' @importFrom scran modelGeneVar getTopHVGs
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom BiocSingular ExactParam
spatialPreprocess <- function(sce, platform=c("Visium", "ST"),
                              n.PCs=15, n.HVGs=2000, skip.PCA=FALSE,
                              log.normalize=TRUE, assay.type="logcounts",
                              h2o.max.mem="5g", n.PCs.image=5,
                              BSPARAM=ExactParam(), ...) {
    ## Set BayesSpace metadata
  if (is.null(metadata(sce)$BayesSpace.data)) {
    metadata(sce)$BayesSpace.data <- list()
  }
  metadata(sce)$BayesSpace.data$platform <- match.arg(platform)
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
  # metadata(sce)$BayesSpace.data$use_dimred <- use.dimred
  # metadata(sce)$BayesSpace.data$d <- n.PCs
  }
    metadata(sce)$BayesSpace.data$platform <- match.arg(platform)
    metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
    # metadata(sce)$BayesSpace.data$use_dimred <- use.dimred
    # metadata(sce)$BayesSpace.data$d <- n.PCs

    ## Run PCA on HVGs, log-normalizing if necessary
    if (!skip.PCA) {
        if (log.normalize)
            sce <- logNormCounts(sce)
   
        dec <- modelGeneVar(sce, assay.type=assay.type)
        top <- getTopHVGs(dec, n=n.HVGs)
        sce <- runPCA(sce, subset_row=top, ncomponents=n.PCs, 
                      exprs_values=assay.type, BSPARAM=BSPARAM)
        rowData(sce)[["is.HVG"]] <- (rownames(sce) %in% top)
    }
    
    ## If the H&E image is provided, run VAE and PCA.
    ## For spot
    if (!is.null(metadata(sce)$BayesSpace.data$spot_image)) {
      ## Get features extracted by VAE.
      metadata(sce)$BayesSpace.data$spot_image_feats <- extractImageFeatures(
        metadata(sce)$BayesSpace.data$spot_image,
        h2o.max.mem,
        ...
      )
      
      ## Get rid of the images to save memory.
      metadata(sce)$BayesSpace.data$spot_image <- NULL
      
      ## Get PCs from VAE features.
      reducedDim(sce, "image") <- scater::calculatePCA(
        metadata(sce)$BayesSpace.data$spot_image_feats,
        ncomponents = n.PCs.image,
        ntop = dim(metadata(sce)$BayesSpace.data$spot_image_feats)[1],
        BSPARAM = BSPARAM
      )
    }
    
    ## For subspot
    if (!is.null(metadata(sce)$BayesSpace.data$subspot_image)) {
      ## Get features extracted by VAE.
      metadata(sce)$BayesSpace.data$subspot_image_feats <- extractImageFeatures(
        metadata(sce)$BayesSpace.data$subspot_image,
        h2o.max.mem,
        ...
      )
      
      ## Get rid of the images to save memory.
      metadata(sce)$BayesSpace.data$subspot_image <- NULL
      
      ## Get PCs from VAE features.
      reducedDim(sce, "image") <- scater::calculatePCA(
        metadata(sce)$BayesSpace.data$subspot_image_feats,
        ncomponents = n.PCs.image,
        ntop = dim(metadata(sce)$BayesSpace.data$subspot_image_feats)[1],
        BSPARAM = BSPARAM
      )
    }

    sce
}

#' @importFrom h2o h2o.init as.h2o h2o.deeplearning h2o.deepfeatures h2o.shutdown
extractImageFeatures <- function(images, h2o.max.mem="5g", ...) {
  h2o.init(
    max_mem_size = h2o.max.mem,
    ...
  )
  
  features <- as.h2o(t(images))
  vae.model <- h2o.deeplearning(x = seq_along(features),
                               training_frame = features,
                               autoencoder = T,
                               hidden = 64,
                               activation = 'Tanh')
  img.feats <- t(as.matrix(h2o.deepfeatures(vae.model, features, layer = 1)))
  
  h2o.shutdown(prompt = FALSE)
  
  colnames(img.feats) <- colnames(images)
  img.feats
}
