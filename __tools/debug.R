Rcpp::compileAttributes()
# devtools::clean_dll()
roxygen2::roxygenise()
devtools::clean_dll()
devtools::load_all()

# Sys.setenv(OMP_NUM_THREADS = 20)

# library(tictoc)
# library(BayesSpace, lib.loc = "/users/skang/R/rocker-rstudio/4.2")
# packageVersion("BayesSpace")

# sce <- readRDS("__tools/results/test.rds")
sce <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/skang/Programming/Git/Xenium_data_exploration/Visium_data/adaptive_mcmc/enhanced_yes.rds")
adjustClusterLabels(sce, 0.2)
# set.seed(149)
# sce <- spatialCluster(
#   sce, q=15, use.dimred = list(image = 1),
#   platform="Visium",
#   init.method="kmeans", model="t", gamma=3,
#   nrep=3000, burn.in=300, save.chain=FALSE
# )

# ret <- pairwiseComp(
#   list(
#     no_img_feat = readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/skang/Programming/Git/Xenium_data_exploration/Visium_data/img_feat/__results/clustered.rds"),
#     img_feat_1 = readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/skang/Programming/Git/Xenium_data_exploration/Visium_data/img_feat/__results/clustered_img_1.rds"),
#     img_feat_2 = readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/skang/Programming/Git/Xenium_data_exploration/Visium_data/img_feat/__results/clustered_img_2.rds"),
#     img_feat_5 = readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/skang/Programming/Git/Xenium_data_exploration/Visium_data/img_feat/__results/clustered_img_5.rds")
#   ),
#   func.comp = aricode::ARI,
#   func.data = function(x) colData(x)$spatial.cluster
# )

# plotHeatmap(ret, name = "ARI", plot.diag = FALSE)

# set.seed(149)
# sce_enhanced <- spatialEnhance(
#   sce,
#   q = 15, platform = "Visium",
#   d = 15,
#   model = "t", gamma = 3,
#   jitter_prior = 0.3, jitter_scale = 3.5,
#   nrep = 200, burn.in = 20,
#   save.chain = TRUE, verbose = TRUE, cores = 2,
#   chain.fname = "__tools/results/adaptive_mcmc/no.h5"
# )

# set.seed(149)
# sce_enhanced <- spatialEnhance(
#   sce,
#   q = 15, platform = "Visium",
#   d = 15,
#   model = "t", gamma = 3,
#   jitter_prior = 0.3, jitter_scale = 0,
#   nrep = 20000, burn.in = 2000,
#   save.chain = TRUE, verbose = TRUE, cores = 2,
#   chain.fname = "__tools/results/adaptive_mcmc/yes.h5"
# )

# img.plots <- imageFeaturePlot(sce, datatype = "both", res = "both",
#                               d = c(pca = 5, vae = 16),
#                               display = list(pca = c(3, 2), vae = c(4, 4)), color = NA)

# num.spots <- -1

# sce_img <- readVisium(
#   "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/skang/Programming/Git/Xenium_data_exploration/Visium_data/data",
#   "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Xenium_Preprint_Data/CytAssist_FFPE_Human_Breast_Cancer/tissue_image.tif",
#   init.backend = TRUE, shutdown.backend = TRUE, num.spots = num.spots, cores = 4)

# sce_2 <- readVisium("__tools/data_2", "__tools/data/tissue_image.tif", "__tools/results_2", FALSE, TRUE, 8, num.spots)


# sce_2 <- read10Xh5("__tools/data")

# n.PCs = 15
#
# set.seed(1)
# sce <- spatialPreprocess(
#   sce,
#   platform = "Visium",
#   n.PCs = n.PCs,
#   log.normalize = TRUE,
#   h2o.max.mem = "500g",
#   h2o.hidden.layer.size = 64,
#   port = 50505
# )
# saveRDS(sce, "__tools/results/preprocessed.rds")

# set.seed(149)
# sce <- spatialCluster(
#   sce, q=15, use.dimred = c(image = 5),
#   platform="Visium",
#   init.method="mclust", model="t", gamma=3,
#   nrep=1000, burn.in=100, save.chain=FALSE
# )
#
# saveRDS(sce, "__tools/results/clustered.rds")

# tmp <- coreTune(
#   sce, 20, 2,
#   q = 15, platform = "Visium",
#   use.dimred = c(PCA = 7),
#   subspot.d = 3,
#   model = "t", gamma = 3,
#   jitter_prior = 0.3, jitter_scale = 3.5,
#   nrep = 500, burn.in = 50,
#   save.chain = FALSE, verbose = TRUE
# )

# set.seed(149)
# sce_enhanced <- spatialEnhance(sce,
#     q = 15, platform = "Visium",
#     use.dimred = list(PCA = seq_len(7)),
#     subspot.d = c(3, 2),
#     model = "t", gamma = 3,
#     jitter_prior = 0.3, jitter_scale = 3.5,
#     nrep = 500, burn.in = 50,
#     save.chain = FALSE, verbose = TRUE, cores = 1
# )

# sce <- readVisium(
#   "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/skang/Programming/Git/Xenium_data_exploration/Visium_data/data"
# )

# sce <- spatialPreprocess(
#   sce,
#   platform = "Visium",
#   log.normalize = TRUE
# )

# set.seed(149)
# sce <- spatialCluster(
#   sce,
#   q = 7,
#   platform = "Visium",
#   init.method = "mclust", model = "t", gamma = 3,
#   nrep = 1000, burn.in = 100, save.chain = FALSE
# )

# saveRDS(sce, "__tools/results/test.rds")

# set.seed(149)
# tic()
# sce_enhanced <- spatialEnhance(
#   sce,
#   q = 7, platform = "Visium",
#   d = 7,
#   model = "t", gamma = 3,
#   jitter_prior = 0.3, jitter_scale = 15,
#   nrep = 10000, burn.in = 1000,
#   save.chain = FALSE, # TRUE,
#   # chain.fname = "__tools/results/chain.h5",
#   verbose = TRUE, cores = 4
# )
# toc()
#
#
#
# detach("package:BayesSpace", unload = TRUE)
# library(BayesSpace, lib.loc = "/Users/se3594/Rlibs")
# packageVersion("BayesSpace")
