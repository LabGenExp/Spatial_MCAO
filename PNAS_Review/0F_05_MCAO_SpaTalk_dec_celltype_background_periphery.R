# OF_05_MCAO_SpaTalk_dec_celltype_background_periphery

# Hello, this R script is made to run some of the SpaTalk time heavy computations on the background. 
# This script matches selected code chunks from 1DP_05_MCAO_spatial_cell_cell_interactions_SpaTalk.
# We will process just the lesions periphery

# day-to-day libraries
library(tidyverse)
library(dplyr)
library(Seurat)
library(ggplot2)
library(openxlsx)
library(magrittr)
library(stringr)
library(patchwork)

# deconvolution and CCI
library(spacexr)
library(SpaTalk)

# set for SpaTalk to operate properly
options(Seurat.object.assay.version = 'v3')

# Let's measure the runtime
start.time <- Sys.time()
cat(paste0("Initiating the computations. Local time is: ", start.time))

{
  # define sections to process
  sections <- c("1DPI", "3DPI", "7DPI")
  cat(crayon::bgBlack("Processing the following sections:", sections, "\n"))
}

# Deconvolution to single-cell resolution
if(file.exists(file.path("data", "SpaTalk_RefZeng2023_MCAO_spatial_periphery.Rds"))){
  cat(crayon::silver("Skipping computing single-cell data. The file already exists.", "\n"))
} else {
  #load the data
  load(file.path("data", "1DP_05_MCAO_dec_celltype_background_periphery.Rdata"))
  cat(crayon::silver(paste0("Objects were loaded. Initiating the deconvolution for single-cell resolution.", "\n")))
  
  # Assign deconvolution data to the SpaTalk object
  dec.spatalk.list <- sections %>% lapply(\(x){
    cat(crayon::yellow(paste0("Processing section: ", x, "\n")))
    # use tryCatch to capture if there are errors
    tryCatch({
      # define core usage
      cores_to_use <- (detectCores() * 0.2) %>% round(0)
      
      # create spatalk object with the single-cell spatial interaction data
      dec_spatalk <- dec_celltype(
        object = spatalk.list[[x]],
        sc_data = result.list[["sc_data"]],
        sc_celltype = result.list[["sc_celltype"]],
        if_use_normalize_data = TRUE,
        if_doParallel = TRUE,
        use_n_cores = cores_to_use,
        method = 2,
        dec_result = prep.list[[x]][["result"]] # or run the rctd deconvolution anew with dec_result = NULL
      )
      cat(crayon::yellow(paste0("Completed section: ", x, "\n")))
      return(dec_spatalk)
    }, error = function(e) {
      message(paste0("Error in section: ", x, " - ", e$message))
    })
  })
  
  # Check if all sections were processed successfully
  if (length(dec.spatalk.list) == length(sections)) {
    # save the output
    names(dec.spatalk.list) <- sections
    saveRDS(dec.spatalk.list, file = file.path("data", "SpaTalk_RefZeng2023_MCAO_spatial_periphery.Rds"))
    cat(crayon::silver("Deconvolution into single-cells is completed. The output was saved.", "\n"))
  } else {
    cat(crayon::red("Some sections failed to process. Review the logs for details.", "\n"))
  }
}

# Infer cell-cell communications. Find ligand-receptor pairs.
if(file.exists(file.path("data", "SpaTalk_RefZeng2023_MCAO_spatial_lr_periphery.Rds"))){
  cat(crayon::silver("Skipping computing L-R pairings. The file already exists.", "\n"))
} else {
  # load the data if not present in the environment
  if(!exists("dec.spatalk.list")){
    dec.spatalk.list <- readRDS(file = file.path("data", "SpaTalk_RefZeng2023_MCAO_spatial_periphery.Rds"))
    cat(crayon::silver(paste0("Objects loaded. Computing L-R pairs.", "\n")))
  }
  
  # compute the lr pairs
  dec.spatalk.list_lr <- sections %>% lapply(\(x){
    cat(crayon::yellow(paste0("Processing section: ", x, "\n")))
    # define cores usage
    cores_to_use <- (detectCores() * 0.2) %>% round(0)
  
    temp <- find_lr_path(
      object = dec.spatalk.list[[x]],
      lrpairs = lrpairs,
      pathways = pathways,
      if_doParallel = T,
      use_n_cores = cores_to_use)
    
    cat(crayon::yellow(paste0("Completed section: ", x, "\n")))
    return(temp)
  })
  
  # save the output
  saveRDS(dec.spatalk.list_lr, file = file.path("data", "SpaTalk_RefZeng2023_MCAO_spatial_lr_periphery.Rds"))
  cat(crayon::silver("Searching LR pairs is done. The output was saved.", "\n"))
}


# Get all cell-cell interaction pairs
# load the file if not present in environment
if(!exists("dec.spatalk.list_lr")){
  dec.spatalk.list_lr <- readRDS(file = file.path("data", "SpaTalk_RefZeng2023_MCAO_spatial_lr_periphery.Rds"))
}

# set up the for loop objects
cat(crayon::silver("Starting with celltype CCI.", "\n"))
for(section in seq_along(dec.spatalk.list_lr)){
  
  section_name <- sections[section]
  cat(crayon::yellow(paste0("Processing section: ", section_name, "\n")))
  
  full_spatalk <- dec_cci_all(
    object = dec.spatalk.list_lr[[section]], 
    if_doParallel = T, 
    use_n_cores = 12
  )
  cat(crayon::yellow(paste0("Completed section: ", section_name, ". Saving the file.", "\n")))
  
  # save the output, section by section
  file_name <- paste0("SpaTalk_RefZeng2023_MCAO_spatial_full_", section_name, "_periphery.Rdata")
  save(full_spatalk, file = file.path("data", file_name))
  
  # clean the memory and environment
  rm(full_spatalk, file_name, section_name)
  gc()
}


end.time <- Sys.time()
time.taken <- end.time - start.time
cat(crayon::magenta(paste0("Job done. Time taken was ", time.taken, " hours.", "\n")))
