analysis.round <- "1st"
num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25
regressOut_mode <- NULL
features_to_regressOut <- NULL
use.sctransform <- TRUE
vars.to.regress <- c("percent.mt")
num.dim.integration <- 25 
num.dim.cluster <- 25
my_random_seed <- 42
cluster.resolution <- 0.5


filter.config.params <- list(
  default = list(nFeatureRNAfloor = NULL,
                 nFeatureRNAceiling = NULL,
                 nCountRNAfloor = NULL, 
                 nCountRNAceiling = NULL,
                 pct_mitofloor = NULL, 
                 pct_mitoceiling = 10,
                 pct_ribofloor = NULL, 
                 pct_riboceiling = NULL,
                 ambientRNA_thres = 0.5,
                 log10GenesPerUMI_thres = NULL,
                 MINCELLS  = 50,
                 MINGENES  = 5)
)