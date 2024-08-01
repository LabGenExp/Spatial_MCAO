# SOP: Helper functions
# author: Daniel Zucha

# ---
# Hello! 
# This script contains helper functions that I encountered using repeatedly. Feel free to modify to your needs! ###


# tidyverse_priority
## resolving conflicts such that tidyverse has priority over other loaded packages
tidyverse_priority <- function() {
  library(conflicted)
  conflicts <- conflict_scout()
  conflicts_tidyverse <- purrr::map(conflicts, ~any(. %in% tidyverse_packages())) %>% unlist
  conflicts <- conflicts[conflicts_tidyverse]
  iwalk(conflicts, ~conflict_prefer(.y, dplyr::intersect(.x, tidyverse_packages())[[1]]))
}


# ggsave: Tiff with 300 DPI
## saving tiffs with constant output parameters. 
## The functions checks for your working space path 'ws' and save the plot there. 
## Adjust height and width in cm.
## not the preferred variant to save is the functino png_save_show(), as it is more flexible and pngs take up less space.
ggsave_tiff <- function(x, ws_location,  plotname, height, width){
  ggsave(plot = x, device = "tiff", dpi = 400, bg = "transparent", filename = paste0(ws_location, "/", plotname, ".tiff"), units = "cm", height = height, width = width)
}

# openxlsx read in function
## xlxs package is great for saving lists into separate tabs of an excel file. 
## Until now however, reconstructing the list back from the multi-tab excel file required more than just a single line of code. Not anymore...
read_all_sheets = function(file, ...) {
  require(openxlsx)
  sheet_list <- lapply(getSheetNames(file), function(sheet){
    read.xlsx(file, sheet = sheet, colNames = T, rowNames = F)
  })
  names(sheet_list) <- getSheetNames(file)
  return(sheet_list)
}

# Variance explained plot
## selecting the correct number of PCs to include in analysis is an important step. 
## Here, I wrote a small function that output a lollipop chart identifying PC contributions to variance explained.
PC_var_explained <- function(Seurat = seurat){
  pca <- Seurat[["pca"]]
  
  ## Get the total variance:
  total_variance <- Seurat@reductions$pca@misc$total.variance
  eigValues = (pca@stdev)^2  ## EigenValues
  varExplained = (eigValues / total_variance)*100
  
  print(paste0(round(varExplained %>% sum, digits = 2), " %"))
  
  df <- data.frame(
    "PC" = Seurat@reductions[["pca"]]@feature.loadings %>% colnames %>% factor(levels = unique(.)),
    "varExplained" = varExplained
  )
  
  p1 <- ggplot(data = df, mapping = aes(x = PC, y = varExplained, color = varExplained)) + 
    geom_point(size = 3) + 
    scale_color_gradient2(low = "#FFF7FB", mid = "#D0D1E6", high = "#023858", midpoint = mean(varExplained)) +
    geom_segment(aes(x=PC, xend=PC, y=0, yend=varExplained)) +
    labs(y= "Variance explained [%]", x="PC components") + 
    theme_light() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p1)
}


# Frequencies barplot
## Display total numbers of populations/conditions as barplot
total_barplot <- function(seurat, groupvar, ylim, x_label = "", y_label = "Cell counts \n"){
  df <- seurat[[groupvar]] %>% 
    table(dnn = ("Celltypes")) %>% 
    sort(decreasing = T) %>% 
    as.data.frame()
  
  return.plot <- df %>% ggplot(data = ., aes(x = Celltypes, y = Freq, label = Freq, fill = Celltypes)) +
    geom_bar(stat = "identity", width = 0.9, alpha = 0.75) + 
    geom_text(hjust = -0.05, vjust = 0.5, color = "black", size = 10, angle = 90) +
    ggtitle(label = NULL) +
    scale_fill_manual(values = col.list[[groupvar]]) +
    ylim(0, ylim) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 105, hjust = 1, vjust = 0, size = 35, color = "black"), 
          axis.text.y = element_blank(),
          axis.title = element_text(size = 25, color = 'black', face = "bold", angle = 90),
          plot.background = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "none") +  
    xlab(x_label) + ylab(y_label)
  
  return(return.plot)
}

# Side-by-side Barplot of proportions
## To create a Barplot that displays the proportions of a categorical value (e.g. Condition) across a grouping variable (var, e.g. celltype).
sidebyside_barplot_percent <- function(seurat, var, group.by){
  ids <- seurat@meta.data[[group.by]] %>% levels
  categories <- seurat@meta.data[[var]] %>% levels
  counts <- matrix(nrow=length(ids), ncol=length(categories)); rownames(counts) <- ids; colnames(counts) <- categories
  
  for (i in seq_along(ids)) {
    for (j in seq_along(categories)) {
      count <- seurat@meta.data[seurat@meta.data[[group.by]] == ids[i] & seurat@meta.data[[var]] == categories[j], ] %>%
        nrow()
      counts[i, j] <- count
    }
  }; rm(i,j)
  counts <- counts/rowSums(counts)
  counts_mm_1 <- reshape2::melt(counts, id = "rownames") ## percent count proportions
  
  
  plot <-  ggplot(data = counts_mm_1, aes(x = Var2, y = value, fill = Var1)) + ## percent count proportions
    geom_col(position = "dodge") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 105, hjust = 1, vjust = 0, size = 35, color = "black"), 
          axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 18, color = "black"),
          axis.title = element_text( size = 25),
          plot.background = element_blank(),
          panel.background = element_blank()) +
    ylab('Proportion of cells\n') +
    xlab(label = NULL) + ggtitle(label = NULL) + NoLegend() + ylim(0,1) +
    scale_fill_manual(values = col.list[[group.by]])
  
  return(plot)
}


## Density plot for Seurat metadata
density_seurat <- function(seu, feature, groups, groups_fill_colors){
  
  stopifnot(groups %in% colnames(seu@meta.data) | feature %in% colnames(seu@meta.data))
  
  plot <- ggplot(data = seu@meta.data, aes(.data[[feature]], fill = .data[[groups]])) + 
    geom_density(alpha = 0.8, color = "white", linewidth = 0.5) +
    scale_fill_manual(values = groups_fill_colors) + 
    theme_mk +
    remove_grid +
    theme(legend.position = c(0.8, 0.75)) +
    labs(x = NULL, y = "Density", fill = NULL)
  
  return(plot)
}


## Seurat's enhanced violin plot (VlnPlot) with boxplots and optional significance testing
enh_vlnplot <- function(seu, feature, grouping, colors, compare_means = TRUE, stat_test = "t.test", paired = F, ref.group = NULL, ...){
  plot <- seu %>% VlnPlot(
    features = feature, 
    group.by = grouping, 
    pt.size = 0, 
    cols = colors, 
    ...) &
    geom_boxplot(width=0.1, fill="white", outlier.size = 0.2) & # add boxplot inside the violins
    xlab(NULL) &
    ylab(NULL)
     
    # statistical test if wanted
  if(isTRUE(compare_means)){
    require(ggpubr)
    plot <- plot & ggpubr::stat_compare_means(
      method = stat_test, 
      ref.group = ref.group, 
      paired = paired, 
      label = "p.signif")
  }
    return(plot)
}


## Correlation of selected metadata (`celltypes`) against a metadata variable (`correlate_against`)
cor_celltype_metadata <- function(seu, celltypes, correlate_against, grouping) {
  require(glue)
  # initiate an empty dataframe to store results
  correlation_results <- data.frame(celltype = character(), 
                                    correlation_coefficient = numeric(),
                                    p_value = numeric(), 
                                    stringsAsFactors = FALSE)
  
  # correlate the values
  for(celltype in celltypes) {
    stopifnot(c(celltype, correlate_against) %in% colnames(seu@meta.data))
    
    cor_test_result <- cor.test(seu@meta.data[[correlate_against]], seu@meta.data[[celltype]], method = "pearson")
    
    ## append the results to the dataframe
    correlation_results <- rbind(correlation_results, 
                                 data.frame(Celltype = celltype, 
                                            correlation_coefficient = cor_test_result$estimate %>% round(digits = 2), 
                                            p_value = cor_test_result$p.value))
  }
  
  # correct for multiple comparisons
  correlation_results %<>% mutate(pvaladj = p.adjust(p_value, method = "fdr") %>% round(digits = 3)) 
  
  # plot the results
  plot_list <- list()
  
  
  for(celltype in celltypes) {
    coef_r <- correlation_results %>% filter(celltype == Celltype) %>% pull(correlation_coefficient)
    coef_p <- correlation_results %>% filter(celltype == Celltype) %>% pull(pvaladj)
    
    p <- ggplot(data = seu@meta.data, aes(x = .data[[celltype]], y = .data[[correlate_against]])) + 
      geom_point(aes(color = grouping, alpha = 0.8)) + 
      scale_color_manual(values = col.list$pals_tableau20) +
      geom_smooth(method = "lm", se = FALSE, color = "#3d3d3d", linewidth = 0.5) + 
      theme_minimal() +
      remove_grid + 
      theme_mk + 
      labs(title = glue("{celltype} \n 
    r={coef_r}   P={coef_p}"), 
    x = "Score", 
    y = paste0(correlate_against)) + NoLegend() + remove_grid
    
    plot_name <- paste0("ScatterPlot_corr_", correlate_against, "_vs_", celltype)
    ## append to the spatial list
    plot_list[[plot_name]] <- p
  }
  return(plot_list)
}

## custom SpatialPlot for metadata. Driven mostly as the crop = F in the Seurat's leaves too much white space, or crop = T distorts the image
stplot <- \(seu, md.feature, section, colors = NULL, image = NULL){
  # seu = seurat object with spatiald data
  # md.feature = metadata feature to plot
  # section = a feature that has identical values in names(seu@images) and seu@meta.data
  # colors = color palette
  # image = only specific image to be plotted?
  
  md <- seu@meta.data
  
  # how many images are to be plotted?
  if(is.null(image)){
    section.levels <- md[[section]] %>% levels
  } else {
    section.levels <- image
  }
  
  # prep the data for plotting
  data.list <- lapply(section.levels, \(level){
    
    # filter for the metadata in one section
    md %<>% filter(.data[[section]] %in% level)
    md %<>% rownames_to_column("barcodes")
    
    # add coordinates
    pos <- spatial.seurat@images[[level]]@coordinates[, c("imagerow", "imagecol")] %>% 
      rename(
        x = imagecol, 
        y = imagerow
      ) %>% 
      mutate(y = y %>% "*"(-1)) # turn the y coordinates upside-down
    pos %<>% rownames_to_column("barcodes")
    
    # add the coordinates to metadata
    md %<>% right_join(pos) %>% column_to_rownames("barcodes")
    
    
    # plot 
    st_plot <- ggplot(data = md, aes(x, y, fill = .data[[md.feature]])) + 
      geom_point(shape = 21, size = 0.5, stroke = 0) + 
      theme_minimal() + 
      remove_grid
    
    if(class(md[[md.feature]]) %in%  c("factor", "character")) {
      st_plot <- st_plot + scale_fill_manual(values = colors)
    } else if (class(md[[md.feature]]) %in%  c("numeric", "integer")) {
      limits <- c(
        min(md[[md.feature]]),
        max(md[[md.feature]])
      )
      st_plot <- st_plot + scale_fill_gradientn(colours = colors, limits = limits)
    } else {
      message("Sorry, I do not know this data class type.")
    }
    
    return(st_plot)
    
  })
}


# min-max normalization
## custom function to bound quantitative values to range(0,1). Good for normalizing i.e. Module Scores 
min_max_normalize <- \(x){(x - min(x)) / (max(x) - min(x))}

# ModuleScoring for bulk data
## Implemented from https://github.com/HerpelinckT/geneset-modulescoring
AddGeneSetScore <- function(
    dds,
    features,
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    name = 'Set',
    seed = 123
) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  object <- counts(dds)
  object.old <- object
  object <- object %||% object.old
  
  features <- list(features)
  features.old <- features
  
  if (is.null(x = features)) {
    stop("Missing input feature list")
  }
  features <- lapply(
    X = features,
    FUN = function(x) {
      missing.features <- setdiff(x = x, y = rownames(x = object))
      if (length(x = missing.features) > 0) {
        warning(
          "The following features are not present in the object: ",
          paste(missing.features, collapse = ", ")
        ) 
        warning(
          paste0("\n ",
                 paste(missing.features, collapse = ", "),
                 " dropped for calculating the geneset score."
          )
        )
      }
      return(intersect(x = x, y = rownames(x = object)))
    }
  )
  
  geneset.length <- length(x = features)
  
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = object[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = geneset.length)
  for (i in 1:geneset.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
          size = ctrl,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = object)
  )
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = object[features.use, ])
  }
  features.scores <- matrix(
    data = numeric(length = 1L),
    nrow = geneset.length,
    ncol = ncol(x = object)
  )
  for (i in 1:geneset.length) {
    features.use <- features[[i]]
    data.use <- object[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:geneset.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  
  range01 <- lapply(
    X = features.scores.use,
    FUN = function(x) {
      range01 <- (x-min(x))/(max(x)-min(x))
    }
  )
  
  range01 <- as.data.frame(x = range01)
  rownames(x = range01) <- colnames(object)
  
  colData(dds) <- cbind(colData(dds), range01)
  
  return(dds)
}

# Compute co-expression in spatial data
compute_coexpression <- \(dataset, ligand, receptor, pair_name){
  # compute expression value of LR pair (the minimum of the two)
  mtx <- dataset %>% 
    GetAssayData(assay = "Spatial", layer = "data") %>% 
    .[c(ligand, receptor), ] %>% 
    t %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Barcodes") %>% 
    dplyr::rename(
      exp_ligand = .data[[ligand]], 
      exp_receptor = .data[[receptor]]) %>% 
    rowwise() %>% 
    mutate(!!pair_name := case_when(
      exp_ligand == 0 | exp_receptor == 0 ~ 0, 
      TRUE ~ min(c(exp_ligand, exp_receptor))
    )) %>% 
    ungroup() %>% 
    select(-exp_ligand, -exp_receptor)
  
  # add to spatial.seurat data
  dataset@meta.data %<>% rownames_to_column(var = "Barcodes")
  dataset@meta.data %<>% right_join(y = mtx, by = "Barcodes")
  dataset@meta.data %<>% 
    mutate(barcode = Barcodes) %>% 
    column_to_rownames(var = "Barcodes")
  dataset@meta.data[[pair_name]] %>% summary %>% print
  
  return(dataset)
}

## predefined theme for gg plotting
theme_mk <- theme_bw() + theme(
  text = element_text(size = 7),
  axis.text.x = element_text(size = 7),
  axis.text.y = element_text(size = 7),
  legend.text = element_text(size = 7),
  strip.text = element_text(size = 7),
  plot.title = element_text(hjust = 0.5, size = 7),
  legend.key.size = unit(2.5, "mm"),
  line = element_line(size = 0.3),
  legend.background=element_blank()
)

## remove grid from ggplot backgrounds
remove_grid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## define plot layout size
plot_layout_mm <- function(..., width, height) {
  plot_layout(..., widths = unit(width, "mm"), heights = unit(height, "mm"))
}


## a helper function to show and save a plot if its size is defined.
png_save_show <- function(plot, file, show_plot = TRUE, ...) {
  bro::bro_ggsave_paged(gg = plot,
                        filename = file,
                        units = "mm",
                        ...)
  if(show_plot) {
    knitr::include_graphics(file)  
  }
}

# rotate images in spatial seurat
## a helper function to rotate /flip images in spatial seurat object.
# Helper function to rotate matrices
# Helper function to rotate matrices
rotimat <- function(foo, rotation) {
  # Check if input is a matrix
  if (!is.matrix(foo)) {
    cat("Input is not a matrix\n")
    return(foo)
  }
  
  # Check if the rotation argument is valid
  if (!(rotation %in% c("180", "Hf", "Vf", "R90", "L90"))) {
    cat("Rotation should be either L90, R90, 180, Hf, or Vf\n")
    return(foo)
  }
  
  # Apply the specified rotation
  if (rotation == "180") {
    foo <- foo %>%
      .[, dim(.)[2]:1] %>%  # Reverse columns
      .[dim(.)[1]:1, ]  # Reverse rows
  }
  
  if (rotation == "Hf") {  # Horizontal flip
    foo <- foo %>%
      .[, dim(.)[2]:1]  # Reverse columns
  }
  
  if (rotation == "Vf") {  # Vertical flip
    foo <- foo %>%
      .[dim(.)[1]:1, ]  # Reverse rows
  }
  
  if (rotation == "L90") {  # Left 90 degrees
    foo <- t(foo)  # Transpose
    foo <- foo %>%
      .[dim(.)[1]:1, ]  # Reverse rows after transpose
  }
  
  if (rotation == "R90") {  # Right 90 degrees
    foo <- t(foo)  # Transpose
    foo <- foo %>%
      .[, dim(.)[2]:1]  # Reverse columns after transpose
  }
  
  return(foo)
}

# Function to rotate images in a Seurat object
rotateSeuratImage <- function(seuratVisumObject, slide = "slice1", rotation = "Vf") {
  # Check if the rotation argument is valid
  if (!(rotation %in% c("180", "Hf", "Vf", "L90", "R90"))) {
    cat("Rotation should be either 180, L90, R90, Hf, or Vf\n")
    return(NULL)
  } else {
    # Retrieve the Seurat object
    seurat.visium <- seuratVisumObject
    
    # Retrieve the image array
    ori.array <- (seurat.visium@images)[[slide]]@image
    
    # Get the image dimensions
    img.dim <- dim(ori.array)[1:2] / (seurat.visium@images)[[slide]]@scale.factors$lowres
    new.mx <- list()
    
    # Transform the image array for each RGB channel
    for (rgb_idx in 1:3) {
      each.mx <- ori.array[,,rgb_idx]
      each.mx.trans <- rotimat(each.mx, rotation)
      new.mx <- c(new.mx, list(each.mx.trans))
    }
    
    # Construct new RGB image array
    new.X.dim <- dim(each.mx.trans)[1]
    new.Y.dim <- dim(each.mx.trans)[2]
    new.array <- array(c(new.mx[[1]], new.mx[[2]], new.mx[[3]]), 
                       dim = c(new.X.dim, new.Y.dim, 3))
    
    # Swap old image with new image
    seurat.visium@images[[slide]]@image <- new.array
    
    # Change the tissue pixel-spot index
    img.index <- (seurat.visium@images)[[slide]]@coordinates
    
    # Swap index based on rotation type
    if (rotation == "Hf") {
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[2] - img.index$imagecol
    }
    
    if (rotation == "Vf") {
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[1] - img.index$imagerow
    }
    
    if (rotation == "180") {
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[1] - img.index$imagerow
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[2] - img.index$imagecol
    }
    
    if (rotation == "L90") {
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[2] - img.index$imagecol
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.index$imagerow
    }
    
    if (rotation == "R90") {
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.index$imagecol
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[1] - img.index$imagerow
    }
    
    return(seurat.visium)
  }
}



# Give some color to your life! ... and analyses
## custom color list
library(RColorBrewer)

col.list <- list()
big_col_palette_f <- function(){
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  return(col_vector)
}; col.list[["big_col_palette"]] <- big_col_palette_f()
col.list[["large_col_palette"]] <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]


col.list[["ctrl_cols"]] <- c("Isocortex L1-4" = "#EBCC2A", "Isocortex L5" = "#E58601", "Isocortex L6" = "#70ACE0", "GLS" = "#AD5248", "Lateral ventricle"  = "#DD8D29", "Hippocampus" = "#E6A0C4", "Fiber Tracts" = "#6F8E78", "Caudoputamen" = "#F8AFA8", "Thalamus" = "#7294D4", "Hypothalamus" =  "#E1AF00", "Amygdalar nuclei" = "#46ACC8", "Cortical subplate" = "#899DA4", "Piriform area" = "#C27D38")

col.list[["cols_mono_short"]] <- c("AMY" = "#D0D0CA", "CP" = "#A1A19C", "CS" = "#B1B1AC", "FT" = "#F0EFE2", "GLS" = "#A1A19C", "HIP" = "#A1A19C", "HY" = "#B1B1AC", "ISD1c" = "#F8EFF1", "ISD1p" = "#E7A9B1", "ISD3c" = "#F1D2D4", "ISD3p" = "#AF3039",  "ISD7c" = "#C65E6A", "ISD7p" =   "#7B1113", "CTX1-4" = "#B1B1AC", "CTX5" = "#D0D0CA", "CTX6" = "#C1C0BB", "lCTX4-5" = "#C1C0BB", "lCTX6" = "#D0D0CA", "LV" = "#B1B1AC", "PAL" = "#D0D0CA", "PIR" = "#C1C0BB", "TH" = "#C1C0BB")
col.list[["cols_mono_short_v2"]] <- c("AMY" = "#D0D0CA", 
                                      "CC" = "#F0EFE2", 
                                      "CP" = "#A1A19C", 
                                      "CS" = "#B1B1AC", 
                                      "FT" = "#F0EFE2", 
                                      "GLS" = "#A1A19C", 
                                      "HIP" = "#A1A19C", 
                                      "HY" = "#B1B1AC", 
                                      "ISD3c" = "#F1D2D4", 
                                      "ISD3p" = "#AF3039",  
                                      "CTX1-4" = "#B1B1AC", 
                                      "CTX5-6" = "#D0D0CA",
                                      "V" = "#B1B1AC", 
                                      "PAL" = "#D0D0CA", 
                                      "PIR" = "#C1C0BB", 
                                      "TH" = "#C1C0BB")

col.list[["cols_mono"]] <- c("Amygdalar area" = "#D0D0CA", 
                             "Caudoputamen" = "#A1A19C", 
                             "Cortical subplate" = "#B1B1AC", 
                             "Fiber tracts" = "#F0EFE2", 
                             "GLS" = "#A1A19C", 
                             "Hippocampus" = "#A1A19C", 
                             "Hypothalamus" = "#B1B1AC", 
                             "Ischemic area D1 center" = "#F8EFF1", 
                             "Ischemic area D1 periphery" = "#E7A9B1", 
                             "Ischemic area D3 center" = "#F1D2D4", 
                             "Ischemic area D3 periphery" = "#AF3039",  
                             "Ischemic area D7 center" = "#C65E6A", 
                             "Ischemic area D7 periphery" =   "#7B1113", 
                             "Isocortex L1-4" = "#B1B1AC", 
                             "Isocortex L5" = "#D0D0CA", 
                             "Isocortex L6" = "#C1C0BB", 
                             "Lateral isocortex L4/5" = "#C1C0BB", 
                             "Lateral isocortex L6" = "#D0D0CA", 
                             "Lateral ventricle" = "#B1B1AC", 
                             "Pallidum" = "#D0D0CA", 
                             "Piriform area" = "#C1C0BB", 
                             "Thalamus" = "#C1C0BB");

col.list[["Timepoint"]] <- c("sham" = "#525475",
                             "D1" = "#F1D7AB",
                             "D3" = "#E6B188", 
                             "D7" = "#D18080")

col.list[["Phase"]] <- c("G1" = "#1B9E77", 
                         "G2M" = "#D95F02", 
                         "S" = "#7570B3")

col.list$alloursn <- c(
  'Neuro 1' = "#2F535C",
  'Glutamatergic 1' = '#99c7e3',
  'Glutamatergic 2' = '#8abfdf',
  'Glutamatergic 3' = '#7ab6da',
  'Glutamatergic 4' = '#6baed6',
  'Glutamatergic 5' = '#5ca6d2',
  'Glutamatergic 6' = '#4c9dcd',
  'Glutamatergic 7' = '#3d95c9',
  'Glutamatergic 8' = '#358bbe',
  'Glutamatergic 9' = '#307faf',
  'Glutamatergic 10' = '#2c749f',
  'GABAergic 1' = '#2fba5d',
  'GABAergic 2' = '#2baa55',
  'GABAergic 3' = '#279b4d',
  'GABAergic 4' = '#238b45',
  'GABAergic 5' = '#1f7b3d',
  'GABAergic 6' = '#1b6c36',
  'NB' = '#D0D1E6',
  'EPEN' = "#016C59", 
  'ASTRO' = "#02818A",
  'ASTRO/OLIGO' = "#014636",
  'OPC' = "#BFD3E6",
  'OLIGO' = "#8C96C6",
  'MG' = "#ba2f47",  # #eb5a5e
  'VLMC' = "#DE77AE",
  'PER/Endo' = "#FEE0B6"
)

col.list[["Celltypes"]] <- c(
  'Astrocytes' = "#02818A",
  'EndothelialCells' = "#ed8b03", 
  'EpendymalCells' = "#016C59",
  'Microglia' = "#CB181D",
  'Neuroblasts' = '#D0D1E6',
  'NeuronsGABA' = '#238b45',
  'NeuronsGLUT' = '#4c9dcd',
  'OPCs' = "#BFD3E6",
  'OLs' = "#8C96C6",
  'VLMCs' = "#c67503",
  'debris' = "grey60"
)

col.list[["CelltypesDetailed"]] <- c(
  'Astrocytes' = "#02818A",
  'Astrocytes2' = "#014636",
  'EndothelialCells' = "#ed8b03", 
  'EpendymalCells' = "#016C59",
  'Microglia1' = '#EF3B2C',
  'Microglia2' = "#CB181D",
  'Neuroblasts' = '#D0D1E6',
  'NeuronsGABA1' = '#2fba5d',
  'NeuronsGABA2' = '#2baa55',
  'NeuronsGABA3' = '#279b4d',
  'NeuronsGABA4' = '#238b45',
  'NeuronsGABA5' = '#1f7b3d',
  'NeuronsGABA6' = '#1b6c36',
  'NeuronsGLUT1' = '#99c7e3',
  'NeuronsGLUT2' = '#8abfdf',
  'NeuronsGLUT3' = '#7ab6da',
  'NeuronsGLUT4' = '#6baed6',
  'NeuronsGLUT5' = '#5ca6d2',
  'NeuronsGLUT6' = '#4c9dcd',
  'NeuronsGLUT7' = '#3d95c9',
  'NeuronsGLUT8' = '#358bbe',
  'NeuronsGLUT9' = '#307faf',
  'OPCs' = "#BFD3E6",
  'OLs1' = "#8C96C6",
  'OLs2' = "#4D3D70",
  'VLMCs' = "#c67503"
)

col.list[["astro_cols"]] <- c(
  "Neuroblasts" = "#D0D1E6", 
  "Ependymal" = "#016C59", 
  "ChoroidPlexus" = "#67A9CF", 
  "TE Astrocytes" = "#3690C0", 
  "DE Astrocytes" = "#02818A", 
  "Activated Astrocytes" = "#A6BDDB", 
  "Reactive Astrocytes" = "#014636"
)

col.list[["mg_cols"]] <- c(
  "Homeostatic" = "#FC9272",
  "Reactive" = "#A50F15"
)


col.list[["intOLSubtypes"]] <- c(
  "OPCs" = '#D0D1E6', 
  "NFOL" = "#BFD3E6",
  "MOL2" = "#9EBCDA",
  "MOL5/6" = "#8C96C6",
  "MOL_IFN" = "#A23C8A",
  "MOL_DA1" = "#4D3D70",
  "MOL_DA2" = "#211E61"
)

col.list[["Condition"]] <- c("Sham" = "#525475",
                             "Ctrl" = "#525475",
                             "1DPI" = "#F1D7AB",
                             "2DPI" = "#E6B188",
                             "3DPI" = "#E6B188", 
                             "7DPI" = "#D18080", 
                             "10DPI" = "#D18080", 
                             "21DPI" = "#AC4C53")

col.list[["BrainAreas"]] <- c("Cortex" = "grey60",
                             "Telencephalon" = "#48759E",
                             "Diencephalon" = "#00A087CC", 
                             "LesionCore" = "#F1D2D4", 
                             "LesionPeriphery" = "#AF3039",
                             "Lesion" = "#AF3039"
                             )

col.list[["Phase"]] <- c("G1" = "#1B9E77", 
                         "G2M" = "#D95F02", 
                         "S" = "#7570B3")

col.list[['milich.cols']] <- c(
  'Astrocytes' = "#02818A",
  'Ependymal' = "#016C59",
  'Oligodendrocytes' = "#8C96C6",
  'OPCs' = "#BFD3E6",
  'Endothelial' = "#ed8b03",
  'Mural cells' = "#c67503",
  'Fibroblasts' = "#BDBDBD",
  'Microglia Homeostatic' = '#EF3B2C',
  'Microglia Chemotactic' = "#CB181D",
  'Microglia Proliferating' = "#A50F15",
  'Microglia Reactive' = "#67000D",
  'Neurons' = '#5ca6d2',
  'T cells' = "#FEE6CE",
  'Neutrophils' = "#FDD0A2",
  'Dendritic' = "#FDAE6B",
  'Monocytes' = "#FD8D3C",
  'Macrophages' = "#F16913"
)

col.list[['zeng.cols']] <- c(
  'Astrocytes' = "#02818A",
  'EpendymalCells' = "#016C59",
  'EndothelialCells' = "#ed8b03",
  "ChoroidPlexusCells" = "#78AE99",
  'DendriticCells' = "#FDAE6B",
  'Fibroblasts' = "#BDBDBD",
  'Granulocytes' = "#FDD0A2",
  'Microglia' = "#A50F15",
  'Neuroblast' = '#5ca6d2', 
  'Neurons' = '#307faf',
  'OLs' = "#8C96C6",
  'OPCs' = "#BFD3E6",
  'Pericytes' = "#D66597", 
  'PeripheralMyeloidCells' = "#F16913", 
  "ProliferatingCells" = "#8C68A4",
  'Tcells' = "#FEE6CE",
  'VSMCs' = "#c67503"
)


col.list[["gradient_blue_to_red"]] <- c("#D01B1B", "#FF4242", "#FFFFFF", "#e7f9ff", "#95D2EC", "#47abd8") %>% rev
col.list[["gradient_grey_to_red"]] <- c("#C2CCD5", "#FCE283", "#F2B83A", "#ED872D", "#E7622F", "#D93934")
col.list[["gradient_blue_yellow_red"]] <- c("#48759E", "#709FC1", "#e7f9ff", "#EAD76E", "#E2AA60", "#CE5C5E", "#AC4C53")
col.list[["Retro_FiveDiscrete"]] <- c("#79A69B", "#CC945C", "#5F4C62", "#BAD4C6", "#B86453")
col.list[["Retro_SixDiscrete"]] <- c("#434C6D", "#6EA89E", "#E6D7C1", "#F4C05B", "#B65F4D", "#553B23")
col.list[["Retro_ThreeDiscrete"]] <- c("#6EA89E","#F4C05B", "#B65F4D")
col.list[["Retro_TenDiscreet"]] <- c("#E64B35CC", "#4DBBD5CC", "#00A087CC", "#3C5488CC", "#F39B7FCC", "#8491B4CC", "#91D1C2CC", "#DC0000CC", "#7E6148CC", "#B09C85CC")

## large discrete color palettes
col.list[["pals_tableau20"]] <- pals::tableau20()

