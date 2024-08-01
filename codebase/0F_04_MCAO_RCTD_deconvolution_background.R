# RCTD deconvolution on the background

# Hello, in this script we will run RCTD deconvolution on the background.

## libraries
library(tidyverse)
library(tibble)
library(magrittr)
library(stringr)
library(Matrix)
library(spacexr)

# lists
plot.list <- list()
result.list <- list()

# load 
load(file = file.path("data", "1DP_04_MCAO_rctd_zeng2023_background.Rdata"))
cat(crayon::silver("Data were loaded. Creating RCTD objects.", "\n"))

# create RCTD objects
myRCTD.list <- list()
for(replicate in sections){
  
  cat(crayon::yellow(paste0("Processing section:", replicate, "\n")))
  
  myRCTD.list[[replicate]] <- create.RCTD(
    spatialRNA = puck[[replicate]], 
    reference = reference, 
    max_cores = 32)
  
  cat(crayon::yellow(paste0("Finished processing section:", replicate, "\n")))
}

rm(replicate)

# deconvolute
cat(crayon::silver("Objects are ready. Starting with RCTD.", "\n"))
celltypes_levels <- reference@cell_types %>% levels
gc()

for(i in sections){
  
  cat(crayon::yellow(paste0("Processing section:", i, "\n")))
  
  myRCTD.list[[i]] <- run.RCTD(
    myRCTD.list[[i]], 
    doublet_mode = "full") ## deconvo itself; mode = "full" as we expect multiple cell types present in each spot.
  
  cat(crayon::yellow(paste0("Finished processing section:", i, "\n")))
  
}

rm(i)

# save the output
saveRDS(myRCTD.list, file = "data/RCTDlist_scRefZeng2023_fullmode.Rds")
save(myRCTD.list, sections, celltypes_levels, file = "data/RCTDlist_scRefZeng2023_fullmode.Rdata")

cat(crayon::silver("RCTD finished. The outputs were saved.", "\n"))