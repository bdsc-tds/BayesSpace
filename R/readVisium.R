#' Load a Visium spatial dataset as a SingleCellExperiment.
#'
#' @param dirname Path to spaceranger output directory (e.g. "sampleID/outs/").
#'   This directory must contain the counts matrix and feature/barcode TSVs in
#'   \code{filtered_feature_bc_matrix/} for \code{readVisium}, or in
#'   \code{filtered_feature_bc_matrix.h5} for \code{read10Xh5}. Besides, it
#'   must also contain a file for spot positions at
#'   \code{spatial/tissue_positions_list.csv} (before Space Ranger V2.0) or
#'   \code{spatial/tissue_positions.csv} (since Space Ranger V2.0).
#'   (To understand the output directory, refer to the corresponding
#'   \href{https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/overview}{10X Genomics help page}.)
#' @param fname File name of the h5 file. It should be inside \code{dirname}.
#'   (By default "filtered_feature_bc_matrix.h5")
#' @param fullres.image Single H&E brightfield image in either TIFF or JPG
#'   format; used as input to "spaceranger count".
#' @param tile.image.dir Path to the directory of tile images.
#' @param scale.factor.fname File name of the scale factor.
#'
#' @return SingleCellExperiment containing the counts matrix in \code{counts}
#'   and spatial data in \code{colData}. Array coordinates for each spot are
#'   stored in columns \code{array_row} and \code{array_col}, while image
#'   coordinates are stored in columns \code{pxl_row_in_fullres} and
#'   \code{pxl_col_in_fullres}.
#'
#' @details We store two variables associated with downstream BayesSpace
#'   functions in a list called \code{BayesSpace.data} in the
#'   SingleCellExperiment's \code{metadata}.
#'   \itemize{
#'     \item \code{platform} is set to "Visium", and is used to determine spot
#'       layout and neighborhood structure.
#'     \item \code{is.enhanced} is set to \code{FALSE} to denote the object
#'       contains spot-level data.
#'   }
#'
#' @examples
#' \dontrun{
#' sce <- readVisium("path/to/outs/")
#' }
#'
#' @name readVisium
NULL

#' @export
#' @importFrom utils read.csv read.table
#' @importFrom scater uniquifyFeatureNames
#' @importFrom Matrix readMM colSums
#' @importFrom SingleCellExperiment SingleCellExperiment counts
#' @importFrom S4Vectors metadata metadata<-
#' @rdname readVisium
readVisium <- function(dirname, fullres.image = NULL, tile.image.dir = NULL,
                       init.backend = TRUE, shutdown.backend = TRUE,
                       cores = 1, num.spots = -1) {
  spatial_dir <- file.path(dirname, "spatial")
  matrix_dir <- file.path(dirname, "filtered_feature_bc_matrix")

  if (!dir.exists(matrix_dir)) {
    stop("Matrix directory does not exist:\n  ", matrix_dir)
  }
  if (!dir.exists(spatial_dir)) {
    stop("Spatial directory does not exist:\n  ", spatial_dir)
  }

  rowData <- read.table(file.path(matrix_dir, "features.tsv.gz"), header = FALSE, sep = "\t")
  colnames(rowData) <- c("gene_id", "gene_name", "feature_type")
  rowData <- rowData[, c("gene_id", "gene_name")]
  rownames(rowData) <- scater::uniquifyFeatureNames(rowData$gene_id, rowData$gene_name)

  .counts <- Matrix::readMM(file.path(matrix_dir, "matrix.mtx.gz"))
  barcodes <- read.table(file.path(matrix_dir, "barcodes.tsv.gz"), header = FALSE, sep = "\t")
  colData <- .read_spot_pos(spatial_dir, barcodes)
  colnames(.counts) <- barcodes$V1
  rownames(.counts) <- rownames(rowData)
  .counts <- .counts[, rownames(colData)]

  sce <- SingleCellExperiment(
    assays = list(counts = .counts),
    rowData = rowData,
    colData = colData
  )

  # Remove spots with no reads for all genes.
  sce <- sce[, Matrix::colSums(counts(sce)) > 0]

  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- "Visium"
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

  if (!is.null(fullres.image)) {
    if (!file.exists(fullres.image)) {
      stop("Full resolution image does not exist:\n  ", fullres.image)
    }

    scale_factor_fname <- file.path(spatial_dir, "scalefactors_json.json")
    if (!file.exists(scale_factor_fname)) {
      stop("Scale factor file does not exist:\n  ", scale_factor_fname)
    }

    if (num.spots <= 0) {
      num.spots <- dim(sce)[2]
    } else {
      sce <- sce[, seq_len(num.spots)]
    }

    flattened_tiles <- create_tiles(
      sce, fullres.image, scale_factor_fname,
      tile.image.dir,
      num.spots, init.backend, shutdown.backend,
      cores
    )

    metadata(sce)$BayesSpace.data$spot_image <- flattened_tiles$spot
    metadata(sce)$BayesSpace.data$subspot_image <- flattened_tiles$subspot
  }

  sce
}

#' @export
#' @importFrom Matrix sparseMatrix colSums readMM
#' @importFrom SingleCellExperiment SingleCellExperiment counts
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom rhdf5 h5read
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by mutate select n case_when
#' @importFrom tibble column_to_rownames
#' @rdname readVisium
read10Xh5 <- function(dirname, fname = "filtered_feature_bc_matrix.h5",
                      fullres.image = NULL, tile.image.dir = NULL,
                      init.backend = TRUE, shutdown.backend = TRUE,
                      cores = 1, num.spots = -1) {
  spatial_dir <- file.path(dirname, "spatial")
  h5_file <- file.path(dirname, fname)

  if (!dir.exists(spatial_dir)) {
    stop("Spatial directory does not exist:\n  ", spatial_dir)
  }

  if (!file.exists(h5_file)) {
    stop("H5 file does not exist:\n  ", h5_file)
  }

  colData <- .read_spot_pos(spatial_dir)

  non.zero.indices <- .extract_indices(
    h5read(h5_file, "matrix/indices"),
    h5read(h5_file, "matrix/indptr")
  )

  rowData <- data.frame(
    gene_id = h5read(h5_file, "matrix/features/id"),
    gene_name = h5read(h5_file, "matrix/features/name")
  ) %>%
    group_by(
      gene_name
    ) %>%
    mutate(
      idx = 1:n(),
      row_name = case_when(
        max(idx) > 1 ~ paste(gene_name, gene_id, sep = "_"),
        TRUE ~ gene_name
      )
    ) %>%
    column_to_rownames("row_name") %>%
    select(
      -idx
    )

  .counts <- sparseMatrix(
    i = non.zero.indices$i,
    j = non.zero.indices$j,
    x = h5read(h5_file, "matrix/data"),
    dims = h5read(h5_file, "matrix/shape"),
    dimnames = list(
      rownames(rowData),
      h5read(h5_file, "matrix/barcodes")
    ),
    index1 = FALSE
  )
  .counts <- .counts[, rownames(colData)]

  sce <- SingleCellExperiment(
    assays = list(
      counts = .counts
    ),
    rowData = rowData,
    colData = colData
  )

  # Remove spots with no reads for all genes.
  sce <- sce[, Matrix::colSums(counts(sce)) > 0]

  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- "Visium"
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

  if (!is.null(fullres.image)) {
    if (!file.exists(fullres.image)) {
      stop("Full resolution image does not exist:\n  ", fullres.image)
    }

    scale_factor_fname <- file.path(spatial_dir, "scalefactors_json.json")
    if (!file.exists(scale_factor_fname)) {
      stop("Scale factor file does not exist:\n  ", scale_factor_fname)
    }

    if (num.spots <= 0) {
      num.spots <- dim(sce)[2]
    } else {
      sce <- sce[, seq_len(num.spots)]
    }

    flattened_tiles <- create_tiles(
      sce, fullres.image, scale_factor_fname,
      tile.image.dir,
      num.spots, init.backend, shutdown.backend,
      cores
    )

    metadata(sce)$BayesSpace.data$spot_image <- flattened_tiles$spot
    metadata(sce)$BayesSpace.data$subspot_image <- flattened_tiles$subspot
  }

  sce
}

#' @export
#' @importFrom rjson fromJSON
#' @rdname readVisium
create_tiles <- function(sce, fullres.image, scale.factor.fname,
                         tile.image.dir = NULL, num.spots = dim(sce)[2],
                         init.backend = TRUE, shutdown.backend = TRUE,
                         cores = 1) {
  if (!file.exists(fullres.image)) stop("Full resolution image does not exist:\n", fullres.image)
  if (!file.exists(scale.factor.fname)) stop("Scale factor file does not exist:\n", scale.factor.fname)

  if (is.null(tile.image.dir)) tile.image.dir <- tempdir(check = TRUE)
  message(paste0("Tiles are stored in ", tile.image.dir))

  .barcodes <- as.vector(colData(sce)[seq_len(num.spots), "barcode"])
  flattened_tiles <- get_spot_subspot_tiles_from_image(
    .barcodes,
    as.matrix(colData(sce)[seq_len(num.spots), c("pxl_col_in_fullres", "pxl_row_in_fullres")]),
    fromJSON(file = scale.factor.fname)$spot_diameter_fullres / 2,
    fullres.image, tile.image.dir, init.backend, shutdown.backend, cores
  )

  colnames(flattened_tiles$spot) <- .barcodes
  colnames(flattened_tiles$subspot) <- flattened_tiles$subspot_barcodes

  flattened_tiles
}

#' @export
#' @importFrom rhdf5 h5read h5createFile h5createGroup h5write
#' @rdname readVisium
counts2h5 <- function(dirname) {
  h5_file <- file.path(dirname, "filtered_feature_bc_matrix.h5")
  spatial_dir <- file.path(dirname, "spatial")
  matrix_dir <- file.path(dirname, "filtered_feature_bc_matrix")

  if (file.exists(h5_file)) {
    stop("H5 file exists:\n ", h5_file)
  }

  if (!dir.exists(matrix_dir)) {
    stop("Matrix directory does not exist:\n  ", matrix_dir)
  }
  if (!dir.exists(spatial_dir)) {
    stop("Spatial directory does not exist:\n  ", spatial_dir)
  }

  rowData <- read.table(file.path(matrix_dir, "features.tsv.gz"), header = FALSE, sep = "\t")
  colnames(rowData) <- c("gene_id", "gene_name", "feature_type")

  .counts <- readMM(file.path(matrix_dir, "matrix.mtx.gz"))
  counts <- matrix(
    as.integer(as.matrix(.counts)),
    nrow = dim(.counts)[1]
  )
  .counts <- as(.counts, "CsparseMatrix")
  barcodes <- read.table(file.path(matrix_dir, "barcodes.tsv.gz"), header = FALSE, sep = "\t")
  colData <- .read_spot_pos(spatial_dir, barcodes)

  h5createFile(h5_file)

  h5createGroup(h5_file, "matrix")

  h5write(barcodes[[1]], h5_file, "matrix/barcodes")
  h5write(counts[counts > 0], h5_file, "matrix/data")

  h5createGroup(h5_file, "matrix/features")
  h5write("genome", h5_file, "/matrix/features/_all_tag_keys")
  h5write(rowData$feature_type, h5_file, "/matrix/features/feature_type")
  h5write(rowData$gene_id, h5_file, "/matrix/features/id")
  h5write(rowData$gene_name, h5_file, "/matrix/features/name")

  h5write(.counts@i, h5_file, "matrix/indices")
  h5write(.counts@p, h5_file, "matrix/indptr")
  h5write(dim(counts), h5_file, "matrix/shape")

  NULL
}

#' Load spot positions.
#'
#' @param dirname Path to spaceranger outputs of spatial pipeline, i.e., "outs/spatial".
#'     This directory must contain a file for the spot positions at
#'     \code{tissue_positions_list.csv} (before Space Ranger V2.0) or
#'     \code{tissue_positions.csv} (since Space Ranger V2.0).
#'
#' @return Data frame of spot positions.
#'
#' @keywords internal
#'
#' @importFrom utils read.csv
#' @importFrom tibble as_tibble
#' @importFrom dplyr inner_join
#' @importFrom tidyr uncount
.read_spot_pos <- function(dirname, barcodes = NULL) {
  if (file.exists(file.path(dirname, "tissue_positions_list.csv"))) {
    message("Inferred Space Ranger version < V2.0")
    colData <- read.csv(file.path(dirname, "tissue_positions_list.csv"), header = FALSE)
    colnames(colData) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  } else if (file.exists(file.path(dirname, "tissue_positions.csv"))) {
    message("Inferred Space Ranger version >= V2.0")
    colData <- read.csv(file.path(dirname, "tissue_positions.csv"))
  } else {
    stop("No file for spot positions found in ", dirname)
  }

  if (!is.null(barcodes)) {
    colData <- inner_join(
      colData,
      barcodes,
      by = c("barcode" = "V1")
    )
  }

  # Sanity check.
  if (
    abs(cor(colData$array_row, colData$pxl_row_in_fullres)) < abs(cor(colData$array_row, colData$pxl_col_in_fullres)) &&
      abs(cor(colData$array_col, colData$pxl_col_in_fullres)) < abs(cor(colData$array_col, colData$pxl_row_in_fullres))
  ) {
    message("Warning! The coordinates with indices and pixels do not match. Swapping the pixel values between the row and column...")

    colData <- transform(
      colData,
      pxl_row_in_fullres = pxl_col_in_fullres,
      pxl_col_in_fullres = pxl_row_in_fullres
    )
  }

  rownames(colData) <- colData$barcode
  colData <- colData[colData$in_tissue > 0, ]
  return(colData)
}

#' Extract row and column indices of the count matrix from h5 file.
#'
#' @param idx Row index of corresponding element in the non-zero count matrix.
#' @param new.start Index of the start of each column corresponding to
#'     \code{idx} and the non-zero count matrix.
#' @param zero.based Whether the \code{} and \code{} are zero-based or not.
#'     (By default is TRUE)
#'
#' @return List of row (i) and column (j) indices of the non-zero elements
#'     in the count matrix.
#'
#' @keywords internal
#'
#' @importFrom tibble as_tibble
#' @importFrom tidyr uncount
.extract_indices <- function(idx, new.start, zero.based = TRUE) {
  if (length(idx) < 1) {
    return(NULL)
  }

  idx.cnts <- do.call(
    rbind,
    lapply(
      seq_len(length(new.start))[-1],
      function(x) c(x - ifelse(zero.based, 2, 1), new.start[[x]] - new.start[[x - 1]])
    )
  )
  colnames(idx.cnts) <- c("id", "n")

  return(
    list(
      i = idx,
      j = as.integer(uncount(as_tibble(idx.cnts), n)[[1]]),
      new.start = new.start
    )
  )
}
