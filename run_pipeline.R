gc()
rm(list = ls())

my_random_seed <- 42
set.seed(my_random_seed)

# if ("DoubletFinder" %in% installed.packages() == FALSE){
#   # remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
#   remotes::install_github('chris-mcginnis-ucsf/DoubletFinder@3b420df')
# }

# remove.packages("glmGamPoi")
# BiocManager::install('glmGamPoi', update = FALSE)
# remove.packages("DoubletFinder")
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder@3b420df')

# if ("glmGamPoi" %in% installed.packages() == FALSE){
#   BiocManager::install('glmGamPoi', update = FALSE)
# }


sctoolbox.dir <- "/home/hieunguyen/src/sctoolboxs"
pipeline.version <- "v0.2"
source(file.path(sctoolbox.dir, "single_cell_pipelines", pipeline.version, "import_libraries.R"))
source(file.path(sctoolbox.dir, "single_cell_pipelines", pipeline.version, "helper_functions.R"))
source(file.path(sctoolbox.dir, "single_cell_pipelines", pipeline.version, "scRNA_GEX_pipeline.R"))


PROJECT <- "CancerSCEM2"

config.version <- "default"

path.to.anno.contigs <- NULL
path.to.count.clonaltype <- NULL
filtered.barcodes <- NULL
path.to.s3a <- NULL

PROJECT.with.version <- sprintf("%s_%s", PROJECT, config.version)

# outdir <- "/media/HNSD01/outdir/cancerscem2"
outdir <- "/media/HNSD02/outdir/cancerscem2"

path.to.main.output <- file.path(outdir, PROJECT.with.version)

path.to.main.input <- file.path("/media/HNSD01/raw_data", PROJECT)

all.input.files <- Sys.glob(file.path(path.to.main.input, "counts", "*.tsv"))
print(sprintf("Number of input files: %s", length(all.input.files)))

finished.samples <- Sys.glob(file.path("/media/HNWD02/outdir/cancerscem2/CancerSCEM2_default/1st_round", "*", "s8a_output", "*")) 
finished.samples <- to_vec(
  for (item in finished.samples) str_replace(str_split(item, "/")[[1]][[8]], "_1st_round", "")
)
print(sprintf("Number of finished pipeline files: %s", length(finished.samples)))

all.input.files <- to_vec(
  for (item in all.input.files) if (basename(item) %in% finished.samples == FALSE) item
)
print(sprintf("Number of to-run input files: %s", length(all.input.files)))

path.to.project.src <- file.path("/home/hieunguyen/src/cancerscem2")
source(file.path(path.to.project.src, "config.R"))

analysis.round <- "1st"

filter.thresholds <- filter.config.params[[config.version]]

MINCELLS <- filter.thresholds$MINCELLS
MINGENES <- filter.thresholds$MINGENES

input.stage_lst <- basename(all.input.files)
# names(input.stage_lst) <- to_vec (for (item in basename(all.input.files)) str_replace(item, ".counts.matrix.tsv", ""))
names(input.stage_lst) <- basename(all.input.files)

save.RDS <- list(s1 = TRUE,
                 s2 = TRUE,
                 s3 = TRUE,
                 s4 = TRUE,
                 s5 = TRUE,
                 s6 = TRUE,
                 s7 = TRUE,
                 s8 = TRUE,
                 s8a = TRUE,
                 s9 = FALSE)

sw <- list(s1 = "on",
           s2 = "on",
           s3 = "on",
           s4 = "on",
           s5 = "on",
           s6 = "on",
           s7 = "on",
           s8 = "off",
           s8a = "on",
           s9 = "off")

rerun <- list(s1 = FALSE, 
              s2 = FALSE,
              s3 = FALSE,
              s4 = FALSE,
              s5 = FALSE,
              s6 = FALSE,
              s7 = FALSE,
              s8 = FALSE,
              s8a = FALSE,
              s9 = FALSE)

filter.thresholds <- filter.config.params[[config.version]]
remove_doublet <- FALSE
path.to.10X.doublet.estimation <- file.path(sctoolbox.dir, "/resources/DoubletEstimation10X.csv")
with.VDJ <- FALSE
DE.test <- "wilcox"

path2src <- file.path(sctoolbox.dir, "single_cell_pipelines", pipeline.version, "processes_src")

for (sample.id in names(input.stage_lst)){
# sample.id <- "AML-029-27-1E.counts.matrix.tsv"
  print(sprintf("Working on sample.id %s", sample.id))
  path.to.output <- file.path(path.to.main.output, sprintf("%s_round", analysis.round))
  dir.create(path.to.output, showWarnings = FALSE)
    
  input.method <- "from_count_matrix"
  path.to.output <- file.path(path.to.output, sprintf("%s_%s_round", sample.id, analysis.round))
  dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(path.to.output, sprintf("%s_%s_round", sample.id, analysis.round), "logs"), showWarnings = FALSE, recursive = TRUE)
  
  path.to.s8a.output <- file.path(path.to.output, "s8a_output/CancerSCEM2_default.output.s8a.rds")
  if (file.exists(path.to.s8a.output) == FALSE){
    tryCatch(
      expr = {
        s.obj <- run_pipeline_GEX(path2src=path2src,
                                  path2input=file.path(path.to.main.input, "counts"),
                                  path.to.logfile.dir=file.path(path.to.output, sprintf("%s_%s_round", sample.id, analysis.round), "logs"),
                                  stage_lst=input.stage_lst[sample.id],
                                  path.to.10X.doublet.estimation=path.to.10X.doublet.estimation,
                                  MINCELLS=MINCELLS,
                                  MINGENES=MINGENES,
                                  PROJECT=PROJECT.with.version,
                                  remove_doublet=remove_doublet,
                                  save.RDS=save.RDS,
                                  path.to.output=path.to.output,
                                  rerun=rerun,
                                  DE.test="DE.test",
                                  num.PCA=num.PCA,
                                  num.PC.used.in.UMAP=num.PC.used.in.UMAP,
                                  num.PC.used.in.Clustering=num.PC.used.in.Clustering,
                                  use.sctransform=use.sctransform,
                                  filtered.barcodes=filtered.barcodes,
                                  filter.thresholds=filter.thresholds,
                                  path.to.anno.contigs=path.to.anno.contigs,
                                  path.to.count.clonaltype=path.to.count.clonaltype,
                                  input.method = input.method,
                                  my_random_seed = my_random_seed,
                                  with.VDJ = with.VDJ, 
                                  path.to.s3a.source = path.to.s3a, 
                                  regressOut_mode = regressOut_mode,
                                  features_to_regressOut = features_to_regressOut,
                                  sw = sw,
                                  vars.to.regress = vars.to.regress,
                                  cluster.resolution = cluster.resolution)     
        # ***** clean up, remove all files in S1 to S7, keep only s8a output. 
        for (i in seq(1, 7)){
          system(sprintf("rm -rf %s", file.path(path.to.output, sprintf("s%s_output", i))))
        }
      },
      error = function(e) {
        message(sprintf("Error with sample %s", sample.id))
        return(NA) 
      }
    )    
  } else {
    print(sprintf("File %s", path.to.s8a.output))
  }
}

#### ALWAYS REMEMBER TO SAVE SESSIONINFO !!!!!!
writeLines(capture.output(sessionInfo()), file.path(path.to.output, sprintf("%s_sessionInfo.txt", PROJECT)))
