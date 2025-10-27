gc()
rm(list = ls())

library(Seurat)
library(tidyverse)
library(dplyr)
library(stringr)
library(comprehenr)

path.to.pipeline.src <- "/home/hieunguyen/src/sctoolboxs/single_cell_pipelines/v0.2"
source(file.path(path.to.pipeline.src, "processes_src", "s8_integration_and_clustering_SeuratV5.selectedGenes.R"))

path.to.main.src <- "/home/hieunguyen/src/cancerscem2"
path.to.celltype.annot <- "/media/HNSD01/raw_data/CancerSCEM2/celltype_annotation"
outdir <- "/media/HNSD01/outdir/cancerscem2"

output.version <- "20251026"

# candidate.cohorts <- c("LUSC", "LUAD", "NSCLC", "STAD", "CRC", "HCC", "BRCA")
candidate.cohorts <- c("NSCLC")

for (input.cohort in candidate.cohorts){
  print("---------------------------------------------------------------")
  print(sprintf("***** WORKING ON COHORT %s *****", input.cohort))
  print("---------------------------------------------------------------")
  
  path.to.main.output <- file.path(outdir, output.version)
  path.to.01.output <- file.path(path.to.main.output, "01_output", input.cohort)
  dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)
  
  dataset.metadata <- read.csv(file.path(path.to.main.src, "processed_metadata.csv"))
  
  cohort.metadata <- subset(dataset.metadata, dataset.metadata$ORGAN == input.cohort & 
                              dataset.metadata$status == "finished")
  
  # save this version of metadata to output
  writexl::write_xlsx(cohort.metadata, file.path(path.to.01.output, "integrated_metadata.xlsx"))
  
  normal.metadata <- subset(cohort.metadata, cohort.metadata$Sample.Type == "Normal")
  tumor.metadata <- subset(cohort.metadata, cohort.metadata$Sample.Type == "Tumour")
  pbmc.metadata <- subset(cohort.metadata, cohort.metadata$Sample.Type == "PBMC")
  
  all.samples <- c(normal.metadata$Sample.ID, tumor.metadata$Sample.ID, pbmc.metadata$Sample.ID)
  
  print(sprintf("Number of samples in NORMAL: %s", nrow(normal.metadata)))
  print(sprintf("Number of samples in TUMOR: %s", nrow(tumor.metadata)))
  print(sprintf("Number of samples in PBMC: %s", nrow(pbmc.metadata)))
  
  if (file.exists(file.path(path.to.01.output, sprintf("raw_merge_dataset_%s.rds", input.cohort))) == FALSE){
    data.list <- list()
    
    for (sample.id in all.samples){
      print(sprintf("working on sample %s", sample.id))
      input.path <- subset(dataset.metadata, dataset.metadata$Sample.ID == sample.id)$path
      tmp.s.obj <- readRDS(str_replace(input.path, "Volumes", "media"))
      
      tmp.celltypedf <- read.csv(file.path(path.to.celltype.annot, sprintf("%s.cell.type.txt", sample.id)), sep = "\t")
      colnames(tmp.celltypedf) <- c("barcode", "celltype")
      tmp.celltypedf <- tmp.celltypedf %>% rowwise() %>%
        mutate(barcode = str_replace_all(barcode, "[.]", "-"))
      
      tmp.metadata <- tmp.s.obj@meta.data %>% rownames_to_column("barcode") %>%
        rowwise() %>%
        mutate(barcode2 = str_replace_all(barcode, "[.]", "-"))
      
      tmp.metadata <- merge(tmp.metadata, tmp.celltypedf, by.x = "barcode2", by.y = "barcode") %>%
        column_to_rownames("barcode")
      tmp.metadata <- tmp.metadata[row.names(tmp.s.obj@meta.data), ]
      tmp.metadata$SampleID <- sample.id
      
      tmp.s.obj <- AddMetaData(object = tmp.s.obj, col.name = "celltype", metadata = tmp.metadata$celltype)
      tmp.s.obj <- AddMetaData(object = tmp.s.obj, col.name = "SampleID", metadata = tmp.metadata$SampleID)
      
      # DimPlot(object = tmp.s.obj, reduction = "RNA_UMAP", label = TRUE, label.box = TRUE, group.by = "celltype") +
      #   ggtitle(sample.id)
      
      # save data to the summary data.list object. 
      data.list[[sample.id]] <- tmp.s.obj
    }
    s.obj <- merge(data.list[[1]], data.list[2: length(data.list)])
    saveRDS(s.obj, file.path(path.to.01.output, sprintf("raw_merge_dataset_%s.rds", input.cohort)))    
  } else {
    print("reading in saved merged object ...")
    s.obj <- readRDS(file.path(path.to.01.output, sprintf("raw_merge_dataset_%s.rds", input.cohort)))
    print("finished reading in saved merged object ...")
  }
  
  # num.PCA <- 25
  # num.PC.used.in.UMAP <- 25
  # num.PC.used.in.Clustering <- 25
  # regressOut_mode <- NULL
  # features_to_regressOut <- NULL
  # use.sctransform <- TRUE
  # vars.to.regress <- c("percent.mt")
  # cluster.resolution <- 0.5
  # 
  # PROJECT <- input.cohort
  # if (file.exists(file.path(path.to.01.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT))) == FALSE){
  #   DefaultAssay(s.obj) <- "RNA"
  #   s.obj <- JoinLayers(s.obj)
  #   s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = s.obj, 
  #                                                        save.RDS.s8 = TRUE,
  #                                                        path.to.output = path.to.01.output,
  #                                                        use.sctransform = TRUE,
  #                                                        num.PCA = num.PCA,
  #                                                        num.PC.used.in.UMAP = num.PC.used.in.UMAP,
  #                                                        num.PC.used.in.Clustering = num.PC.used.in.Clustering,
  #                                                        cluster.resolution = cluster.resolution,
  #                                                        vars.to.regress = vars.to.regress, 
  #                                                        remove.genes = NULL)
  # }
}

