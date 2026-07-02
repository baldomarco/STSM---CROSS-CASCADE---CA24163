
#                                       Dr. Marco Baldo, MSc
# 
#                                11/06/2026  CZU and SAS - CROSS-CASCADE COST Action CA24163


#            This script will read and manipulate the IBL Polish Forest Institute forest sampling plots for the CROSS-CASCADE 
#            COST Action STSM for the study region of the Bialowieza National Park and prepare the input tables required to synthesize 
#            the forest and environmental structures into the iLand model of forest dynamics and disturbances simulations

#            This a clean and automatized script to synthesize the whole Bottoms-Up database into input tables readable in iLand model

library(geosphere)
library(sf)
library(ggpubr)
library(ggplot2)
library(writexl)
library(dplyr)
library(tidyverse)

# Load the data
ALS_clean_alive <- readRDS("C:/iLand/20230901_Bottoms_Up/plot_init/R/stsm_roma/ALS_clean_alive.rds")

# Check where aren't present trees names: deciduous as species.. 
plots_no_deciduous <- ALS_clean_alive %>%
  dplyr::group_by(plotid) %>%
  dplyr::filter(!any(species == "deciduous")) %>%
  dplyr::distinct(plotid) %>%
  dplyr::pull(plotid)
# Save the list they are 641 on 668
writeLines(plots_no_deciduous, "C:/iLand/20230901_Bottoms_Up/plot_init/R/stsm_roma/ALS_plots_without_deciduous.txt")

# Categorize the species to start the replace of these trees with the most common broadl species in the plot
broadleaf_species <- setdiff(unique(ALS_clean_alive$species), c("piab", "pisy", "deciduous"))

# Find the most common broadl species in the plots
broadleaf_mode <- ALS_clean_alive %>%
  dplyr::filter(species %in% broadleaf_species) %>%
  dplyr::group_by(plotid, species) %>%
  dplyr::tally() %>%
  dplyr::group_by(plotid) %>%
  dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
  dplyr::select(plotid, species)

# Create a directory to store the plots
plot_dir <- "C:/iLand/20230901_Bottoms_Up/plot_init/plots/New folder/"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# Use the one here correcting it in case you divide the plots in sites/regions

# Process each site and save into the new folder in .xlsx format
#unique_sites <- unique(ALS_clean_alive$siteID)
#for (site in unique_sites) {
#  site_data <- ALS_clean_alive %>%
#    dplyr::filter(siteID == site)
#  output_dir <- file.path(plot_dir, site)
#  if (!dir.exists(output_dir)) dir.create(output_dir)
#  process_site(site_data, output_dir)
#}

# Process each plot and save into the new folder in .xlsx format 
unique_plots <- unique(ALS_clean_alive$plotid)

for (plot in unique_plots) {
  subset_df <- ALS_clean_alive %>%
    dplyr::filter(plotid == plot) %>%
    tidyr::drop_na(x, y)
  
  if ("deciduous" %in% subset_df$species) {
    replacement <- broadleaf_mode$species[match(plot, broadleaf_mode$plotid)]
    subset_df$species[subset_df$species == "deciduous"] <- replacement
  }
  
  if (nrow(subset_df) == 0) next
  
  subset_df$x <- as.numeric(subset_df$x)
  subset_df$y <- as.numeric(subset_df$y)
  
  write_xlsx(subset_df,
             file.path(plot_dir, paste0("plot_", plot, ".xlsx")))
}


# You can add the section for creating and writing the RU and Stand grids within the loop after calculating the corner points of the plot. Here's the modified script with the added section
# To name the output files based on the plot ID, you can extract the plot ID from the Excel file's name and use it to generate the output file names. Here's how you can modify the script to achieve this

# Load the required libraries
# Load the required libraries
library(readxl)
library(dplyr)
library(raster)
library(fields)

# Define the directory path where your Excel files are located
directory_path <- "C:/iLand/20230901_Bottoms_Up/plot_init/plots/clean_plot/"

# List all Excel files in the directory
file_paths <- list.files(directory_path, pattern = "\\.xlsx$", full.names = TRUE)

# Create a function to process each Excel file
process_excel_file <- function(file_path) {
  # Read the data from the Excel file
  data <- read_excel(file_path)
  
  # Extract the plot ID from the "plotid" column
  plot_id <- unique(data$plotid)
  
  # Calculate the minimum and maximum x and y coordinates of the trees
  min_x <- min(as.numeric(data$x), na.rm = TRUE)
  max_x <- max(as.numeric(data$x), na.rm = TRUE)
  min_y <- min(as.numeric(data$y), na.rm = TRUE)
  max_y <- max(as.numeric(data$y), na.rm = TRUE)
  
  # Plot corners = raw tree bounding box, NO buffer
  corner1 <- c(min_x, min_y)
  corner2 <- c(max_x, min_y)
  corner3 <- c(max_x, max_y)
  corner4 <- c(min_x, max_y)
  
  cat("Corner 1: ", corner1[1], ",", corner1[2], "\n")
  cat("Corner 2: ", corner2[1], ",", corner2[2], "\n")
  cat("Corner 3: ", corner3[1], ",", corner3[2], "\n")
  cat("Corner 4: ", corner4[1], ",", corner4[2], "\n")
  
  #----------------------------------------------------------------
  out.dataroot <- "C:/iLand/20230901_Bottoms_Up/plot_init/gis/"
  
  # Single shared origin used for BOTH the grid files and the tree shift below
  x.coord.corner <- as.integer(corner1[1])
  y.coord.corner <- as.integer(corner1[2])
  
  # RU GRID
  xn <- 1
  yn <- 1
  RU.values <- rep(110, xn * yn)
  RU.grid <- matrix(RU.values, ncol = xn)
  
  # STAND GRID
  xn_ <- 10
  yn_ <- 10
  Stand.values <- rep(110, xn_ * yn_)
  Stand.grid <- matrix(Stand.values, ncol = xn_)
  
  set.panel(2, 2)
  par(mar = c(2, 4, 2, 4))
  
  #---------------------------------------ENVIRONMENT-------------------------------------
  RU.grid.file <- paste0(out.dataroot, "RU_grid_", plot_id, ".asc")
  
  write.table(paste("NCOLS", xn, sep="\t"), file = RU.grid.file, append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(paste("NROWS", yn, sep="\t"), file = RU.grid.file, append = TRUE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(paste("XLLCORNER", x.coord.corner, sep="\t"), file = RU.grid.file, append = TRUE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(paste("YLLCORNER", y.coord.corner, sep="\t"), file = RU.grid.file, append = TRUE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(paste("CELLSIZE", "100", sep="\t"), file = RU.grid.file, append = TRUE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(paste("NODATA_value", "-9999", sep="\t"), file = RU.grid.file, append = TRUE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(RU.grid, file = RU.grid.file, append = TRUE, quote = FALSE, sep = "\t",
              eol = "\n", na = "-9999", dec = ".", row.names = FALSE, col.names = FALSE)
  
  #---------------------------------------- STAND------------------
  S.grid.file <- paste0(out.dataroot, "Stand_grid_", plot_id, ".asc")
  
  write.table(paste("NCOLS", xn_, sep="\t"), file = S.grid.file, append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(paste("NROWS", yn_, sep="\t"), file = S.grid.file, append = TRUE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(paste("XLLCORNER", x.coord.corner, sep="\t"), file = S.grid.file, append = TRUE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(paste("YLLCORNER", y.coord.corner, sep="\t"), file = S.grid.file, append = TRUE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(paste("CELLSIZE", "10", sep="\t"), file = S.grid.file, append = TRUE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(paste("NODATA_value", "-9999", sep="\t"), file = S.grid.file, append = TRUE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(Stand.grid, file = S.grid.file, append = TRUE, quote = FALSE, sep = "\t",
              eol = "\n", na = "-9999", dec = ".", row.names = FALSE, col.names = FALSE)
  
  #------------------------------- TREE TABLE ---------------------------
  # select FROM `data` (the table just read in); real height column is "h"
  desired_columns <- dplyr::select(data, x, y, species, dbh, h)
  colnames(desired_columns) <- c("x", "y", "species", "dbh", "height")
  
  sp <- unique(desired_columns$species)
  
  desired_columns$height <- as.numeric(desired_columns$height)          # meters
  desired_columns$dbh <- round(as.numeric(desired_columns$dbh) / 10, 1) # mm -> cm
  
  #--------------------------------------------------------------------
  # Shift trees so the plot's own bottom-left corner becomes (0,0).
  # Same x.coord.corner/y.coord.corner used for the grid files above,
  # so trees and grids share exactly one consistent origin.
  desired_columns$x <- as.numeric(desired_columns$x) - x.coord.corner
  desired_columns$y <- as.numeric(desired_columns$y) - y.coord.corner
  
  # Clip tiny rounding overflow at the edges to stay within 0-100
  desired_columns$x <- pmin(pmax(desired_columns$x, 0), 100)
  desired_columns$y <- pmin(pmax(desired_columns$y, 0), 100)
  
  desired_columns$x <- round(desired_columns$x, 2)
  desired_columns$y <- round(desired_columns$y, 2)
  
  #-------------------------------------------------------------------
  out.dataroot <- "C:/iLand/20230901_Bottoms_Up/plot_init/gis/init/"
  
  Plot.grid.file <- paste0(out.dataroot, plot_id, "_init.txt")
  
  write.table(desired_columns, file = Plot.grid.file,
              append = FALSE, quote = FALSE, sep = ";", eol = "\n", na = "NA",
              dec = ".", row.names = FALSE, col.names = TRUE)
  
  return(desired_columns)
}

# Use lapply to process all Excel files
processed_data_list <- lapply(file_paths, process_excel_file)

# If you need to access specific processed data, you can do so like this:
# processed_data <- processed_data_list[[index]]  # Replace 'index' with the desired file's index

# Now you have a list of processed data, one for each Excel file in the directory

#-------------------------------------------------------------------------------
#                                       END
#
#-------------------------------------------------------------------------------
