Rcpp::compileAttributes()
roxygen2::roxygenise()
devtools::clean_dll()
devtools::load_all()

Sys.setenv(OMP_NUM_THREADS = 20)

# library(tictoc)
# library(BayesSpace, lib.loc = "/users/skang/R/rocker-rstudio/4.2")

# sce <- readRDS("__tools/results/preprocessed.rds")

# img.plots <- imageFeaturePlot(sce, datatype = "both", res = "both",
#                               d = c(pca = 5, vae = 16),
#                               display = list(pca = c(3, 2), vae = c(4, 4)), color = NA)

# num.spots <- -1

# tic()
# sce <- readVisium(
#   "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/skang/Programming/Git/Xenium_data_exploration/Visium_data/data",
#   "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Xenium_Preprint_Data/CytAssist_FFPE_Human_Breast_Cancer/tissue_image.tif",
#   init.backend = TRUE, shutdown.backend = TRUE, num.spots = num.spots)
# toc()

# tic()
# sce_2 <- readVisium("__tools/data_2", "__tools/data/tissue_image.tif", "__tools/results_2", FALSE, TRUE, 8, num.spots)
# toc()


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
#     use.dimred = c(PCA = 7),
#     subspot.d = 3,
#     model = "t", gamma = 3,
#     jitter_prior = 0.3, jitter_scale = 3.5,
#     nrep = 500, burn.in = 50,
#     save.chain = FALSE, verbose = TRUE, cores = 2
# )

# library(ggplot2)
# library(tictoc)
#
# nrep <- 10000
#
# detach("package:BayesSpace", unload = TRUE)
# library(BayesSpace, lib.loc = "/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library")
# packageVersion("BayesSpace")
# set.seed(100)
# sce_old = spatialPreprocess(
#   getRDS("2018_thrane_melanoma", "ST_mel1_rep2", cache=FALSE),
#   platform = "ST",
#   n.PCs = 7
# )
#
# set.seed(149)
# sce_old <- spatialCluster(
#   sce_old, q=4, platform="ST", d=7,
#   init.method="mclust", model="t", gamma=2,
#   nrep=1000, burn.in=100, save.chain=FALSE
# )
# clusterPlot(sce_old)
# ggsave("__tools/old.jpg", device = "jpg")
#
# set.seed(149)
# sce_old_enhanced <- spatialEnhance(sce_old, q=4, platform="ST", d=7,
#                                     model="t", gamma=2,
#                                     jitter_prior=0.3, jitter_scale=3.5,
#                                     nrep=nrep, burn.in=100,
#                                     save.chain=FALSE)
# clusterPlot(sce_old_enhanced)
# ggsave("__tools/old_enhanced.jpg", device = "jpg")
#
#
#
# detach("package:BayesSpace", unload = TRUE)
# library(BayesSpace, lib.loc = "/Users/se3594/Rlibs")
# packageVersion("BayesSpace")
set.seed(100)
sce_new <- spatialPreprocess(
  getRDS("2018_thrane_melanoma", "ST_mel1_rep2", cache = FALSE),
  platform = "ST",
  n.PCs = 7
)

set.seed(149)
sce_new <- spatialCluster(
  sce_new,
  q = 4, platform = "ST", d = 7,
  init.method = "mclust", model = "t", gamma = 2,
  nrep = 1000, burn.in = 100, save.chain = FALSE
)
# clusterPlot(sce_new)
# ggsave("__tools/new.jpg", device = "jpg")

# tmp <- coreTune(sce_new,
#   test.cores = c(1, 2, 4, 8, 16), test.times = 2,
#   q = 4, platform = "ST", d = 7,
#   model = "t", gamma = 2,
#   jitter_prior = 0.3, jitter_scale = 3.5,
#   nrep = 10000, burn.in = 100, verbose = TRUE
# )

set.seed(149)
# tic()
sce_new_enhanced <- spatialEnhance(sce_new,
  q = 4, platform = "ST", d = 7,
  model = "t", gamma = 2,
  jitter_prior = 0.3, jitter_scale = 3.5,
  nrep = 2000, burn.in = 100,
  save.chain = FALSE, cores = 2, verbose = TRUE
)
# toc()
clusterPlot(sce_new_enhanced)
# ggsave("__tools/new_enhanced_1.jpg", device = "jpg")

set.seed(149)
sce_new_enhanced_2 <- spatialEnhance(sce_new,
  q = 4, platform = "ST", d = 7,
  model = "t", gamma = 2,
  jitter_prior = 0.3, jitter_scale = 3.5,
  nrep = 2000, burn.in = 100,
  save.chain = FALSE, cores = 4, verbose = TRUE
)
clusterPlot(sce_new_enhanced_2)
#
# set.seed(149)
# tic()
# sce_new_enhanced_2 <- spatialEnhance(sce_new, q=4, platform="ST", d=7,
#                                    model="t", gamma=2,
#                                    jitter_prior=0.3, jitter_scale=3.5,
#                                    nrep=nrep, burn.in=100,
#                                    save.chain=FALSE, cores = 4, verbose = TRUE)
# toc()
# clusterPlot(sce_new_enhanced_2)
# ggsave("__tools/new_enhanced_4.jpg", device = "jpg")
