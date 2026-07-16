#                           Dr. Marco Baldo, MSc
#           11/06/2026  CZU and SAS - CROSS-CASCADE COST Action CA22136
#
#   Replicates FID 2015 circular 500m2 plots to 1ha (100x100m) squares assuming
#       homogeneous forest structure, then writes iLand tree init tables.
#       Adapted from 03_tree_initialization.R + _TreeReplicationForData.R

rm(list=ls()) 

library(dplyr)
library(sf)
library(tidyr)
library(writexl)

#===============================================================================
# --- Load data ---
#===============================================================================
FID_2015_clean_alive <- readRDS("C:/iLand/20230901_Bottoms_Up/plot_init/R/stsm_roma/FID_2015_clean_alive.rds")

# Check NAs in single tree coordinates and species name
na_check <- FID_2015_clean_alive %>%
  summarise(
    na_x       = sum(is.na(x)),
    na_y       = sum(is.na(y)),
    na_species = sum(is.na(species))
  )
na_check

# Species missing - Just to know where because later we will fill the appropriate species name
FID_2015_clean_alive %>% filter(is.na(species)) %>% select(plotid, treeid, sp_ibl, sp_name, species)

# Test to see the NA in DBH & H and H within the tree init source data. 
# This is important to avoid errors in the simulation
n_before <- nrow(FID_2015_clean_alive)
na_dbh   <- sum(is.na(FID_2015_clean_alive$dbh))
na_h     <- sum(is.na(FID_2015_clean_alive$h))

na_dbh_trees <- FID_2015_clean_alive %>% filter(is.na(dbh)) %>% select(plotid, treeid)
na_h_trees   <- FID_2015_clean_alive %>% filter(is.na(h))   %>% select(plotid, treeid)

table(na_dbh_trees$plotid)
table(na_h_trees$plotid)


# Drop NA in DBH
FID_2015_clean_alive <- FID_2015_clean_alive %>% tidyr::drop_na(dbh, h)

message(sprintf("Removed %d NA dbh, %d NA h out of %d total rows (%.1f%% / %.1f%%)",
                na_dbh, na_h, n_before,
                100 * na_dbh / n_before,
                100 * na_h   / n_before))


#===============================================================================
# --- Convert iLand Species ----
#===============================================================================
sp <- unique(FID_2015_clean_alive$sp_name)
#write.csv(sp, file.path("C:/P/DMP_CROSS_CASCADE/03_rawdata/Species conversion.csv"), row.names = FALSE)
sp_con <- read.csv("C:/P/DMP_CROSS_CASCADE/03_rawdata/Species conversion.csv")
            
sp_NA<-subset(FID_2015_clean_alive, is.na(FID_2015_clean_alive$sp_name)) 
sp_NA

# First conversion based on the species not present in iLand sp parameters db                    
FID_2015_clean_alive <- FID_2015_clean_alive %>%
  mutate(species = recode(
    species,
    "tipe" = "tipl",
    "saal" = "saca",
    "safr" = "saca",
    "savi" = "saca",
    "masi" = "prse",
    "pyco" = "prse",
    "ulla" = "ulgl",
    "ulca" = "ulgl",
    .default = species
  ))

# Replace the NA
FID_2015_clean_alive <- FID_2015_clean_alive %>%
  mutate(
    species = case_when(
      is.na(species) & sp_ibl == "SL"  ~ "prse",
      is.na(species) & sp_ibl == "WIP" ~ "saca",
      TRUE ~ species
    )
  )


# --- Constants ---
nominal_plot_area_m2 <- 500
plot_radius_m        <- sqrt(nominal_plot_area_m2 / pi)   # 12.6157 m
square_side_m         <- 100
scaling_factor        <- (square_side_m^2) / nominal_plot_area_m2   # 20
center_xy              <- square_side_m / 2                          # 50

# Outer 1ha square, fixed for every plot
outer_polygon <- st_polygon(list(matrix(
  c(0,0, square_side_m,0, square_side_m,square_side_m, 0,square_side_m, 0,0),
  ncol = 2, byrow = TRUE)))

# Exclusion zone = the real 500m2 circle, same for every plot (fixed design, not tree-hull based)
exclusion_circle <- st_buffer(st_point(c(center_xy, center_xy)), dist = plot_radius_m)
valid_sampling_area <- st_difference(outer_polygon, exclusion_circle)

# --- Function: replicate one plot's trees to fill 1 ha ---
replicate_plot_to_1ha <- function(plot_df) {

  # Center original trees on (50,50) using their own centroid as reference
  cx <- mean(plot_df$x, na.rm = TRUE)
  cy <- mean(plot_df$y, na.rm = TRUE)

  original_trees <- plot_df %>%
    mutate(
      treeid = as.character(treeid),
      x = x - cx + center_xy,
      y = y - cy + center_xy,
      Source = "Original"
    )

  # Bootstrap-resample attributes (species, dbh, h, v, ...) with replacement
  n_orig <- nrow(plot_df)
  n_replicate <- round(n_orig * (scaling_factor - 1))   # e.g. 19x original count

  replicated_trees <- plot_df %>%
    slice_sample(n = n_replicate, replace = TRUE) %>%
    mutate(Source = "Replicated")

  # New coordinates for replicated trees: random points in square minus the real plot circle
  new_points  <- st_sample(valid_sampling_area, size = n_replicate, type = "random")
  new_coords  <- st_coordinates(new_points)
  replicated_trees$x <- new_coords[, "X"]
  replicated_trees$y <- new_coords[, "Y"]

  # Unique treeids for the replicates
  replicated_trees$treeid <- paste0(replicated_trees$treeid, "_r", seq_len(nrow(replicated_trees)))

  combined <- bind_rows(original_trees, replicated_trees)
  combined
}

# --- Apply per plot, write per-plot xlsx (same downstream format as 03_tree_initialization.R) ---
plot_dir <- "C:/iLand/20230901_Bottoms_Up/plot_init/plots/New folder/"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

unique_plots <- unique(FID_2015_clean_alive$plotid)

one_ha_list <- list()

for (plot in unique_plots) {
  subset_df <- FID_2015_clean_alive %>%
    dplyr::filter(plotid == plot) %>%
    tidyr::drop_na(x, y)

  if (nrow(subset_df) == 0) next

  one_ha_plot <- replicate_plot_to_1ha(subset_df)
  one_ha_list[[plot]] <- one_ha_plot

  write_xlsx(one_ha_plot, file.path(plot_dir, paste0("plot_", plot, "_1ha.xlsx")))
}


# --- Verification plot for one example plot ---
example_plot <- one_ha_list[[unique_plots[230]]]
ggplot2::ggplot(example_plot, ggplot2::aes(x = x, y = y, color = Source)) +
  ggplot2::geom_point(ggplot2::aes(size = dbh), alpha = 0.7) +
  ggplot2::coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
  ggplot2::labs(title = paste("Simulated 1ha plot -", unique_plots[230])) +
  ggplot2::theme_bw()

#-------------------------------------------------------------------------------
#            TREE INIT TABLE WRITER
#            No bounding-box/corner logic needed: coordinates already live
#            in local 0-100 space by construction. RU/Stand grids come from ALS.
#-------------------------------------------------------------------------------

init_dir <- "C:/iLand/20230901_Bottoms_Up/plot_init/gis/init/"
if (!dir.exists(init_dir)) dir.create(init_dir, recursive = TRUE)

write_tree_init <- function(one_ha_plot, plot_id) {
  desired_columns <- dplyr::select(one_ha_plot, x, y, species, dbh, h)
  colnames(desired_columns) <- c("x", "y", "species", "dbh", "height")

  desired_columns$height <- as.numeric(desired_columns$height)            # meters
  desired_columns$dbh    <- round(as.numeric(desired_columns$dbh) / 10, 1) # mm -> cm

  # Safety clip only (sampling is already bounded to [0,100])
  desired_columns$x <- round(pmin(pmax(desired_columns$x, 0), 100), 2)
  desired_columns$y <- round(pmin(pmax(desired_columns$y, 0), 100), 2)

  out_file <- file.path(init_dir, paste0(plot_id, "_init.txt"))
  write.table(desired_columns, file = out_file,
              append = FALSE, quote = FALSE, sep = ";", eol = "\n", na = "NA",
              dec = ".", row.names = FALSE, col.names = TRUE)
}

for (plot in names(one_ha_list)) {
  write_tree_init(one_ha_list[[plot]], plot)
}

#-------------------------------------------------------------------------------
#            VERIFICATION: did replication preserve forest structure?
#            Compare original 500m2 plot (scaled x20 to per-ha) vs the
#            actual built 1ha replicate, across all plots. BA in m2/ha,
#            N in stems/ha, using standard DBH classes (cm).
#-------------------------------------------------------------------------------

dbh_breaks <- c(0, 10, 20, 30, 40, 50, 60, 80, Inf)
dbh_labels <- c("0-10","10-20","20-30","30-40","40-50","50-60","60-80","80+")

structural_attributes <- function(df, area_m2) {
  df <- df %>% mutate(dbh_cm = as.numeric(dbh) / 10)  # mm -> cm
  scale_to_ha <- 10000 / area_m2

  tibble(
    N_ha    = nrow(df) * scale_to_ha,
    BA_ha   = sum(pi * (df$dbh_cm / 200)^2, na.rm = TRUE) * scale_to_ha,
    dbh_mean = mean(df$dbh_cm, na.rm = TRUE)
  )
}

dbh_class_counts <- function(df, area_m2) {
  scale_to_ha <- 10000 / area_m2
  df %>%
    mutate(dbh_cm = as.numeric(dbh) / 10,
           dbh_class = cut(dbh_cm, breaks = dbh_breaks, labels = dbh_labels, right = FALSE)) %>%
    count(dbh_class, .drop = FALSE) %>%
    mutate(n_ha = n * scale_to_ha)
}

comparison_list <- list()

for (plot in names(one_ha_list)) {
  orig_df <- FID_2015_clean_alive %>% dplyr::filter(plotid == plot) %>% tidyr::drop_na(x, y)
  rep_df  <- one_ha_list[[plot]]

  orig_attr <- structural_attributes(orig_df, nominal_plot_area_m2) %>%
    rename_with(~paste0(., "_orig"))
  rep_attr  <- structural_attributes(rep_df, square_side_m^2) %>%
    rename_with(~paste0(., "_rep"))

  comparison_list[[plot]] <- bind_cols(plotid = plot, orig_attr, rep_attr)
}

comparison_df <- bind_rows(comparison_list)

# --- Save post-processing data ---
out_rds <- "C:/iLand/20230901_Bottoms_Up/plot_init/postprocess/"
if (!dir.exists(out_rds)) dir.create(out_rds, recursive = TRUE)


saveRDS(one_ha_list,
        file = file.path(out_rds, "one_ha_list.rds"))

saveRDS(comparison_df,
        file = file.path(out_rds, "comparison_df.rds"))


# --- Error metrics (bias, MAE, RMSE, R2) treating scaled-original as the reference ---
error_metrics <- function(obs, pred) {
  tibble(
    bias = mean(pred - obs, na.rm = TRUE),
    mae  = mean(abs(pred - obs), na.rm = TRUE),
    rmse = sqrt(mean((pred - obs)^2, na.rm = TRUE)),
    r2   = cor(obs, pred, use = "complete.obs")^2
  )
}

ba_metrics <- error_metrics(comparison_df$BA_ha_orig, comparison_df$BA_ha_rep)
n_metrics  <- error_metrics(comparison_df$N_ha_orig,  comparison_df$N_ha_rep)

cat("Basal area (m2/ha) — original(scaled x20) vs replicate:\n"); print(ba_metrics)
cat("Stem density (N/ha) — original(scaled x20) vs replicate:\n"); print(n_metrics)

# --- Gauch et al. (2003) MSD decomposition: MSD = SB + NU + LC ---
msd_decomposition <- function(obs, pred) {
  n <- length(obs)
  msd <- mean((pred - obs)^2, na.rm = TRUE)
  sb  <- (mean(pred, na.rm = TRUE) - mean(obs, na.rm = TRUE))^2         # squared bias
  b   <- coef(lm(pred ~ obs))[2]                                        # regression slope
  nu  <- (b - 1)^2 * (sum((obs - mean(obs))^2) / n)                     # non-unity slope
  r2  <- cor(obs, pred, use = "complete.obs")^2
  lc  <- (1 - r2) * (sum((pred - mean(pred))^2) / n)                    # lack of correlation
  tibble(MSD = msd, SB = sb, NU = nu, LC = lc,
         SB_pct = sb/msd*100, NU_pct = nu/msd*100, LC_pct = lc/msd*100)
}

cat("MSD decomposition — Basal area:\n")
print(msd_decomposition(comparison_df$BA_ha_orig, comparison_df$BA_ha_rep))
cat("MSD decomposition — Stem density:\n")
print(msd_decomposition(comparison_df$N_ha_orig, comparison_df$N_ha_rep))

# --- DBH class distribution check for one example plot ---
example_id <- unique_plots[240]
orig_classes <- dbh_class_counts(FID_2015_clean_alive %>% dplyr::filter(plotid == example_id), nominal_plot_area_m2)
rep_classes  <- dbh_class_counts(one_ha_list[[example_id]], square_side_m^2)

dbh_compare <- bind_rows(
  orig_classes %>% mutate(Source = "Original x20"),
  rep_classes  %>% mutate(Source = "1ha Replicate")
)

ggplot2::ggplot(dbh_compare, ggplot2::aes(x = dbh_class, y = n_ha, fill = Source)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::labs(title = paste("DBH class distribution -", example_id),
                x = "DBH class (cm)", y = "Stems / ha") +
  ggplot2::theme_bw()

#-------------------------------------------------------------------------------
#                                       END
#-------------------------------------------------------------------------------
