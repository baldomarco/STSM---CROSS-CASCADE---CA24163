# Marco Baldo - contact: baldo@fld.czu.cz
# Verification: FID2015 (500m2 source) vs FID2015 1ha-replicated (one_ha_list) vs ALS2015 (already 1ha)
# Trees + tree volume: numeric and visual, on a plot selection (not all ~668 together)
# Also: does simulation FID2015 t0 match source data? And does simulated stand age match FID plot-description age?

#===============================================================================
# --- START ----
#===============================================================================

rm(list = ls())
library(tidyverse)

base_dir               <- "C:/P/DMP_CROSS_CASCADE"
raw_process_data_path  <- file.path(base_dir, "03_rawdata/02_process_storage")
sim_root               <- "C:/iLand/2026/Bialowieza/RDS"
figure_path             <- file.path(base_dir, "04_work/03_analysis/02_figure")
table_path              <- file.path(base_dir, "04_work/03_analysis/03_table")

# NOTE: update this to wherever 03_tree_initialization.R actually saves the object.
# It is NOT written by 01_data_organization.R -- that script only writes per-plot xlsx.
one_ha_list_path <- file.path("C:/Users/baldo/Documents/GitHub/STSM---CROSS-CASCADE---CA24163/02_iLand_experiment/postprocess/one_ha_list.rds")

#===============================================================================
# --- LOAD DATA ----
#===============================================================================

FID_2015_clean_alive <- readRDS(file.path(raw_process_data_path, "FID_order_rds/FID_2015_clean_alive.rds"))
ALS_clean_alive       <- readRDS(file.path(raw_process_data_path, "ALS_order_rds/ALS_clean_alive.rds"))
one_ha_list           <- readRDS(one_ha_list_path)

# plot description info (field-assessed stand age lives here, "plots_info" sheet)
FID_2015_2019_Info_clean <- readRDS(file.path(raw_process_data_path, "FID_2015_2019_Info_clean.rds"))

# simulation, FID_V1 run (adjust scenario if you want ALS_V1 / V2 instead)
plot_variables_all_FID_V1_FULL <- readRDS(file.path(sim_root, "FID_V1_Full/plot_variables_all_FID_V1_FULL.rds"))
lnd_scen_FID_V1_FULL            <- readRDS(file.path(sim_root, "FID_V1_Full/lnd_scen_FID_V1_FULL.rds"))

#===============================================================================
# --- SANITY CHECKS ----
#===============================================================================

stopifnot(is.list(one_ha_list), !is.data.frame(one_ha_list))

need_tree_cols <- c("plotid", "x", "y", "species", "dbh", "h", "dead")
chk_cols <- function(df, cols, label) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(paste0(label, " missing columns: ", paste(miss, collapse = ", ")))
}
chk_cols(FID_2015_clean_alive, c(need_tree_cols, "v"), "FID_2015_clean_alive")
chk_cols(ALS_clean_alive, need_tree_cols, "ALS_clean_alive")
chk_cols(one_ha_list[[1]], c(need_tree_cols, "v", "Source"), "one_ha_list[[1]]")

has_als_volume <- "v" %in% names(ALS_clean_alive)
if (!has_als_volume) {
  message("ALS_clean_alive has no 'v' (volume) column -> ALS is EXCLUDED from volume comparisons below, ",
          "and included only in tree-count / DBH-structure comparisons.")
}

#===============================================================================
# --- SPECIES GROUPS (kept consistent with 06a_Error_Analysis_general.R) ----
#===============================================================================

species_to_remove <- c("piab", "pisy", "abal", "lade", "psme", "pini", "pice")

#===============================================================================
# --- GENERIC PLOT-LEVEL METRIC FUNCTION (same shape as metrics_one_var in 06a) ----
# moved up here so PART A0 can use it too
#===============================================================================

metrics_one_var <- function(pred_df, obs_df, pred_col, obs_col, label, var_label) {
  
  pred_sel <- pred_df %>% dplyr::select(plotid, pred = dplyr::all_of(pred_col))
  obs_sel  <- obs_df  %>% dplyr::select(plotid, obs  = dplyr::all_of(obs_col))
  
  if (anyDuplicated(pred_sel$plotid) != 0) stop(paste0(label, ": pred_df has duplicate plotid rows."))
  if (anyDuplicated(obs_sel$plotid)  != 0) stop(paste0(label, ": obs_df has duplicate plotid rows."))
  
  joined <- dplyr::inner_join(pred_sel, obs_sel, by = "plotid")
  if (nrow(joined) == 0) {
    warning(paste0("No overlap for ", label))
    return(list(table = joined, summary = tibble::tibble()))
  }
  
  per_plot <- joined %>%
    dplyr::mutate(
      rel_diff_pct = (pred - obs) / (obs + 0.001) * 100,
      sq_diff      = (pred - obs)^2
    )
  
  r <- tryCatch(cor(per_plot$pred, per_plot$obs, use = "complete.obs"), error = function(e) NA_real_)
  if (is.na(r)) r <- 0
  
  summary_row <- tibble::tibble(
    comparison         = label,
    variable            = var_label,
    n_plots             = nrow(per_plot),
    rel_difference_pct  = median(per_plot$rel_diff_pct, na.rm = TRUE),
    #rel_difference_pct  = mean(per_plot$rel_diff_pct, na.rm = TRUE),
    rel_difference_sd   = sd(per_plot$rel_diff_pct, na.rm = TRUE),
    pooled_rel_diff_pct = sum(per_plot$pred - per_plot$obs, na.rm = TRUE) / (sum(per_plot$obs, na.rm = TRUE) + 0.001) * 100,
    MSD                 = mean(per_plot$sq_diff, na.rm = TRUE),
    SB                  = (mean(per_plot$pred, na.rm = TRUE) - mean(per_plot$obs, na.rm = TRUE))^2,
    SDSD                = (sd(per_plot$pred, na.rm = TRUE) - sd(per_plot$obs, na.rm = TRUE))^2,
    LC                  = 2 * sd(per_plot$pred, na.rm = TRUE) * sd(per_plot$obs, na.rm = TRUE) * (1 - r),
    r                   = r
  )
  
  list(table = per_plot, summary = summary_row)
}

#===============================================================================
# PART A0 -- CROP ALS2015 TO THE SAME REAL FOOTPRINT AS FID2015 (500 m2 circle)
# Real detected ALS trees vs real measured FID trees, no replication/scaling
# involved on either side -- this is the actual "same trees, same species?" check.
#===============================================================================

nominal_plot_area_m2 <- 500
plot_radius_m         <- sqrt(nominal_plot_area_m2 / pi)   # 12.6157 m, same constant as 03_tree_initialization.R

plot_centers <- FID_2015_2019_Info_clean %>%
  dplyr::select(plotid, center_x = x, center_y = y) %>%
  dplyr::filter(!is.na(center_x), !is.na(center_y))

# --- Sanity check: are ALS coords and the field plot-center coords in the same frame? ---
# If the nearest ALS tree to the recorded center is much farther than plot_radius_m,
# the crop below is not trustworthy and the coordinate systems need reconciling first.
check_plot  <- plot_centers$plotid[1]
chk_center  <- plot_centers %>% dplyr::filter(plotid == check_plot)
chk_als     <- ALS_clean_alive %>% dplyr::filter(plotid == check_plot, dead == 0)
if (nrow(chk_als) > 0 && nrow(chk_center) > 0) {
  chk_dist <- sqrt((chk_als$x - chk_center$center_x)^2 + (chk_als$y - chk_center$center_y)^2)
  message("Sanity check, plot ", check_plot, ": nearest ALS tree is ", round(min(chk_dist), 1),
          " m from the field-recorded plot center (plot radius = ", round(plot_radius_m, 1), " m).")
  message("If this is much larger than the plot radius across most plots, ALS x/y and ",
          "FID_2015_2019_Info_clean x/y are not directly comparable -- stop and check the CRS/frame before trusting PART A0.")
}

# --- Crop ALS to the real 500m2 circle, per plot ---
als_footprint_cropped <- ALS_clean_alive %>%
  dplyr::filter(dead == 0) %>%
  dplyr::inner_join(plot_centers, by = "plotid") %>%
  dplyr::mutate(dist_to_center = sqrt((x - center_x)^2 + (y - center_y)^2)) %>%
  dplyr::filter(dist_to_center <= plot_radius_m) %>%
  dplyr::select(-center_x, -center_y, -dist_to_center)

# how many ALS trees survive the crop, per plot, vs the full 1ha ALS count -- should be a small fraction (~500/10000 = 5% of area)
als_crop_check <- als_footprint_cropped %>%
  dplyr::group_by(plotid) %>%
  dplyr::summarise(n_trees_ALS_cropped_500m2 = dplyr::n(), .groups = "drop") %>%
  dplyr::full_join(
    ALS_clean_alive %>% dplyr::filter(dead == 0) %>% dplyr::group_by(plotid) %>%
      dplyr::summarise(n_trees_ALS_1ha = dplyr::n(), .groups = "drop"),
    by = "plotid"
  ) %>%
  dplyr::mutate(pct_of_1ha_in_circle = n_trees_ALS_cropped_500m2 / n_trees_ALS_1ha * 100)

summary(als_crop_check$pct_of_1ha_in_circle)

#-------------------------------------------------------------------------------
# Numeric: trees, basal area, species richness -- FID2015 (real) vs ALS2015 (cropped, real)
# no volume here (ALS has none) -- basal area used as the size-structure proxy instead
#-------------------------------------------------------------------------------

footprint_metrics <- function(df) {
  df %>%
    dplyr::mutate(dbh_cm = dbh / 10, ba_m2 = pi * (dbh_cm / 200)^2) %>%
    dplyr::group_by(plotid) %>%
    dplyr::summarise(
      n_trees   = dplyr::n(),
      ba_total  = sum(ba_m2, na.rm = TRUE),
      n_species = dplyr::n_distinct(species),
      .groups = "drop"
    )
}

fid_footprint_metrics <- footprint_metrics(FID_2015_clean_alive %>% dplyr::filter(dead == 0))
als_footprint_metrics <- footprint_metrics(als_footprint_cropped)

trees_footprint_vs <- metrics_one_var(
  pred_df = als_footprint_metrics, obs_df = fid_footprint_metrics,
  pred_col = "n_trees", obs_col = "n_trees",
  label = "ALS_cropped500m2_vs_FID500m2", var_label = "n_trees"
)
trees_footprint_vs$summary

ba_footprint_vs <- metrics_one_var(
  pred_df = als_footprint_metrics, obs_df = fid_footprint_metrics,
  pred_col = "ba_total", obs_col = "ba_total",
  label = "ALS_cropped500m2_vs_FID500m2", var_label = "basal_area_m2"
)
ba_footprint_vs$summary

#-------------------------------------------------------------------------------
# Numeric: species composition -- tree counts per species, same footprint
#-------------------------------------------------------------------------------

species_composition <- dplyr::bind_rows(
  FID_2015_clean_alive %>% dplyr::filter(dead == 0) %>%
    dplyr::count(plotid, species, name = "n") %>% dplyr::mutate(source = "FID2015_500m2"),
  als_footprint_cropped %>%
    dplyr::count(plotid, species, name = "n") %>% dplyr::mutate(source = "ALS2015_cropped500m2")
)

#-------------------------------------------------------------------------------
# Visual: same-footprint comparison, sample plots
#-------------------------------------------------------------------------------

# NEED TO OPEN A PDF WRITER AND GIVE IT THE ROOT, THE NAME, AND THE SIZE
dataroot <- "C:/P/DMP_CROSS_CASCADE/04_work/03_analysis/02_figure/PDF_archive/"
pdf(paste0(dataroot, "20260715_02d_FID_ALS_real_sim_plot_&_age_comparison.pdf"), height=8, width=12)

# Plot selection - if you want to use a specific plot use target_plotid <- "KS194" instead of footprint_sample_plots
set.seed(123)
# Run without the set.seed (123) if you want to change the selection randomly
footprint_sample_plots <- sample(
  intersect(unique(FID_2015_clean_alive$plotid), unique(als_footprint_cropped$plotid)), 6
)

ggplot2::ggplot(
  species_composition %>% dplyr::filter(plotid %in% footprint_sample_plots),
  ggplot2::aes(x = source, y = n, fill = species)
) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(~ plotid, scales = "free_y") +
  ggplot2::labs(title = "Species composition, same 500m2 footprint: FID2015 (real) vs ALS2015 (cropped, real)",
                x = NULL, y = "N trees", fill = "Species") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

footprint_spatial <- dplyr::bind_rows(
  FID_2015_clean_alive %>%
    dplyr::filter(dead == 0, plotid %in% footprint_sample_plots) %>%
    dplyr::mutate(source = "FID2015_500m2") %>%
    dplyr::select(plotid, x, y, dbh, species, source),
  als_footprint_cropped %>%
    dplyr::filter(plotid %in% footprint_sample_plots) %>%
    dplyr::mutate(source = "ALS2015_cropped500m2") %>%
    dplyr::select(plotid, x, y, dbh, species, source)
)

# Make the overlapping trees comparison plots based on the selection
ggplot2::ggplot(footprint_spatial, ggplot2::aes(x = x, y = y, color = source)) +
  ggplot2::geom_point(ggplot2::aes(size = dbh), alpha = 0.6) +
  ggplot2::facet_wrap(~ plotid, scales = "free") +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
  ggplot2::labs(title = "Same 500m2 footprint, sample plots: FID2015 (real) vs ALS2015 (cropped, real)") +
  ggplot2::theme_bw()


# Make single trees comparison in the 500m2 plot FID-ALS ground data per species and dbh size 
library(patchwork)

p_fid <- footprint_spatial %>% dplyr::filter(source == "FID2015_500m2") %>%
  ggplot2::ggplot(ggplot2::aes(x = x, y = y, color = species)) +
  ggplot2::geom_point(ggplot2::aes(size = dbh), alpha = 0.7) +
  ggplot2::facet_wrap(~ plotid, scales = "free") +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
  ggplot2::labs(title = "FID2015_500m2") + ggplot2::theme_bw()

p_als <- footprint_spatial %>% dplyr::filter(source == "ALS2015_cropped500m2") %>%
  ggplot2::ggplot(ggplot2::aes(x = x, y = y, color = species)) +
  ggplot2::geom_point(ggplot2::aes(size = dbh), alpha = 0.7) +
  ggplot2::facet_wrap(~ plotid, scales = "free") +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
  ggplot2::labs(title = "ALS2015_cropped500m2") + ggplot2::theme_bw()

p_fid / p_als


# NOTE: FID x/y are plot-local (centered near each plot's own centroid per 03_tree_initialization.R)
# while ALS x/y here are still real-world CRS 2180 coordinates (not re-centered) -- the facet above
# uses free scales per plot so shape/spread are comparable, but the two point clouds will NOT sit
# on the same absolute x/y unless ALS is recentered by (center_x, center_y) from plot_centers first.
# Uncomment below to recenter ALS on (0,0) for a true overlay:
# footprint_spatial_centered <- footprint_spatial %>%
#   dplyr::left_join(plot_centers, by = "plotid") %>%
#   dplyr::mutate(x = dplyr::if_else(source == "ALS2015_cropped500m2", x - center_x, x),
#                 y = dplyr::if_else(source == "ALS2015_cropped500m2", y - center_y, y))

#===============================================================================
# PART A -- TREE COUNTS: FID500(scaled) vs FID1ha-replicated vs ALS1ha
#===============================================================================

fid500_counts <- FID_2015_clean_alive %>%
  dplyr::filter(dead == 0) %>%
  dplyr::group_by(plotid) %>%
  dplyr::summarise(n_trees_500m2 = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(n_trees_FID500_scaled_1ha = n_trees_500m2 * 20)

fid1ha_counts <- purrr::imap_dfr(
  one_ha_list,
  ~ dplyr::tibble(plotid = .y, n_trees_FID1ha_replicated = nrow(dplyr::filter(.x, dead == 0)))
)

als_counts <- ALS_clean_alive %>%
  dplyr::filter(dead == 0) %>%
  dplyr::group_by(plotid) %>%
  dplyr::summarise(n_trees_ALS_1ha = dplyr::n(), .groups = "drop")

counts_compare <- fid500_counts %>%
  dplyr::full_join(fid1ha_counts, by = "plotid") %>%
  dplyr::full_join(als_counts, by = "plotid") %>%
  dplyr::mutate(
    # should be ~1 : replication is built to conserve tree density, this is the internal consistency check
    replication_ratio_1ha_vs_500scaled = n_trees_FID1ha_replicated / n_trees_FID500_scaled_1ha,
    diff_FID1ha_vs_ALS      = n_trees_FID1ha_replicated - n_trees_ALS_1ha,
    rel_diff_FID1ha_vs_ALS_pct = diff_FID1ha_vs_ALS / (n_trees_ALS_1ha + 0.001) * 100
  )

summary(counts_compare$replication_ratio_1ha_vs_500scaled)
summary(counts_compare$rel_diff_FID1ha_vs_ALS_pct)

#===============================================================================
# PART B -- TREE VOLUME: FID500(scaled x20) vs FID1ha-replicated (NO further scaling)
# ALS excluded here (no v column) -- see has_als_volume check above
#===============================================================================

fid500_volume <- FID_2015_clean_alive %>%
  dplyr::filter(dead == 0, !is.na(v)) %>%
  dplyr::group_by(plotid) %>%
  dplyr::summarise(volume_m3ha_FID500_scaled = sum(v, na.rm = TRUE) * 20, .groups = "drop")

fid1ha_volume <- purrr::imap_dfr(one_ha_list, function(df, id) {
  df <- dplyr::filter(df, dead == 0)
  dplyr::tibble(
    plotid = id,
    volume_m3ha_FID1ha_replicated = sum(df$v, na.rm = TRUE)
  )
})

volume_compare <- fid500_volume %>%
  dplyr::full_join(fid1ha_volume, by = "plotid") %>%
  dplyr::mutate(
    diff_1ha_vs_500scaled     = volume_m3ha_FID1ha_replicated - volume_m3ha_FID500_scaled,
    rel_diff_1ha_vs_500scaled_pct = diff_1ha_vs_500scaled / (volume_m3ha_FID500_scaled + 0.001) * 100
  )

summary(volume_compare$rel_diff_1ha_vs_500scaled_pct)

vol_1ha_vs_500 <- metrics_one_var(
  pred_df = fid1ha_volume, obs_df = fid500_volume,
  pred_col = "volume_m3ha_FID1ha_replicated", obs_col = "volume_m3ha_FID500_scaled",
  label = "FID1ha_replicated_vs_FID500scaled", var_label = "volume_m3ha"
)
vol_1ha_vs_500$summary

trees_1ha_vs_ALS <- metrics_one_var(
  pred_df = fid1ha_counts, obs_df = als_counts,
  pred_col = "n_trees_FID1ha_replicated", obs_col = "n_trees_ALS_1ha",
  label = "FID1ha_replicated_vs_ALS1ha", var_label = "n_trees"
)
trees_1ha_vs_ALS$summary

#===============================================================================
# --- VISUAL: numeric checks above, scatter + 1:1 line, all plots at once (points only, cheap) ----
#===============================================================================

ggplot2::ggplot(vol_1ha_vs_500$table, ggplot2::aes(x = obs, y = pred)) +
  ggplot2::geom_point(alpha = 0.4) +
  ggplot2::geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  ggplot2::labs(
    title = "Volume: FID2015 500m2 (x20 scaled) vs FID2015 1ha-replicated",
    x = "Volume m3/ha - FID500 scaled", y = "Volume m3/ha - FID 1ha replicated"
  ) +
  ggplot2::theme_bw()

ggplot2::ggplot(trees_1ha_vs_ALS$table, ggplot2::aes(x = obs, y = pred)) +
  ggplot2::geom_point(alpha = 0.4) +
  ggplot2::geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  ggplot2::labs(
    title = "Tree count: ALS2015 (1ha) vs FID2015 1ha-replicated",
    x = "N trees - ALS 1ha", y = "N trees - FID 1ha replicated"
  ) +
  ggplot2::theme_bw()

#===============================================================================
# PART C -- SELECTED-PLOT VISUALS (spatial maps), NOT all ~668 plots together
#===============================================================================

# Select a fix plot selection = in this case 6
set.seed(123)
sample_plots <- sample(intersect(names(one_ha_list), unique(ALS_clean_alive$plotid)), 6)

# --- FID2015: original vs replicated trees, per plot (user-provided base pattern, extended to 6 plots) ---
fid1ha_sample <- purrr::imap_dfr(one_ha_list[sample_plots], ~ dplyr::mutate(.x, plotid = .y))

# Plot 500m2 + replicates visualization
ggplot2::ggplot(fid1ha_sample, ggplot2::aes(x = x, y = y, color = Source)) +
  ggplot2::geom_point(ggplot2::aes(size = dbh), alpha = 0.7) +
  ggplot2::coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
  ggplot2::facet_wrap(~ plotid) +
  ggplot2::labs(title = "FID2015 1ha-replicated: original vs replicated trees, sample plots") +
  ggplot2::theme_bw()

# --- Same sample plots: FID2015 1ha-replicated vs ALS2015, side by side ---
als_sample <- ALS_clean_alive %>%
  dplyr::filter(plotid %in% sample_plots, dead == 0) %>%
  dplyr::left_join(plot_centers, by = "plotid") %>%
  dplyr::mutate(x = x - center_x + 50, y = y - center_y + 50, source = "ALS2015", species = as.character(species)) %>%
  dplyr::select(plotid, x, y, dbh, h, species, source)

fid1ha_sample_plain <- fid1ha_sample %>%
  dplyr::filter(dead == 0) %>%
  dplyr::mutate(source = "FID2015_1ha_replicated", species = as.character(species)) %>%
  dplyr::select(plotid, x, y, dbh, h, species, source)

spatial_compare <- dplyr::bind_rows(
  fid1ha_sample_plain,
  als_sample %>% dplyr::select(plotid, x, y, dbh, h, species, source)
)

# Plot comparison ALS - FID 1-ha scale overlapping trees dbh
ggplot2::ggplot(spatial_compare, ggplot2::aes(x = x, y = y, color = source)) +
  ggplot2::geom_point(ggplot2::aes(size = dbh), alpha = 0.5) +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::facet_wrap(~ plotid) +
  ggplot2::labs(title = "FID2015 1ha-replicated vs ALS2015, sample plots") +
  ggplot2::theme_bw()

# --- DBH distribution, sample plots, FID1ha replicated vs ALS ---
ggplot2::ggplot(spatial_compare, ggplot2::aes(x = dbh, fill = source)) +
  ggplot2::geom_histogram(alpha = 0.5, position = "identity", bins = 25) +
  ggplot2::facet_wrap(~ plotid, scales = "free_y") +
  ggplot2::labs(title = "DBH (mm) distribution, sample plots: FID2015 1ha-replicated vs ALS2015") +
  ggplot2::theme_bw()

# --- HEIGHT distribution, sample plots, FID1ha replicated vs ALS ---
ggplot2::ggplot(spatial_compare, ggplot2::aes(x = h, fill = source)) +
  ggplot2::geom_histogram(alpha = 0.5, position = "identity", bins = 25) +
  ggplot2::facet_wrap(~ plotid, scales = "free_y") +
  ggplot2::labs(title = "Height (m) distribution, sample plots: FID2015 1ha-replicated vs ALS2015") +
  ggplot2::theme_bw()

# Let's make the ALS-FID comparison tree by tree colored by species at 1ha scale
target_plotid <- sample_plots[1]

spatial_compare %>%
  dplyr::filter(plotid == target_plotid) %>%
  ggplot2::ggplot(ggplot2::aes(x = x, y = y, color = species)) +
  ggplot2::geom_point(ggplot2::aes(size = dbh), alpha = 0.7) +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::facet_wrap(~ source) +
  ggplot2::labs(title = paste("FID2015 1ha-replicated vs ALS2015 -", target_plotid)) +
  ggplot2::theme_bw()


# --- Explicit single-plot example (kept in original single-object style) ---
example_id  <- sample_plots[1]
example_plot <- one_ha_list[[example_id]]
ggplot2::ggplot(example_plot, ggplot2::aes(x = x, y = y, color = Source)) +
  ggplot2::geom_point(ggplot2::aes(size = h), alpha = 0.7) +
  ggplot2::coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
  ggplot2::labs(title = paste("Simulated 1ha plot -", example_id)) +
  ggplot2::theme_bw()

# Let's compare the total basal area by species
spatial_compare %>%
  dplyr::filter(plotid == target_plotid) %>%
  dplyr::mutate(ba_m2 = pi * (dbh/2000)^2) %>%
  dplyr::group_by(source, species) %>%
  dplyr::summarise(ba_m2 = sum(ba_m2, na.rm = TRUE), .groups = "drop") %>%
  ggplot2::ggplot(ggplot2::aes(x = source, y = ba_m2, fill = species)) +
  ggplot2::geom_col() +
  ggplot2::labs(title = paste("Basal area by species -", target_plotid), y = "Basal area (m2)", x = NULL) +
  ggplot2::theme_bw()


#===============================================================================
# PART D -- SIMULATION vs SOURCE DATA AT t0: does the initialization reproduce
# the FID2015 1ha-replicated trees/volume in the sqlite output at year 0?
#===============================================================================

sim_t0 <- plot_variables_all_FID_V1_FULL %>%
  dplyr::mutate(
    plotid = stringr::str_remove(run, "DB_") %>% stringr::str_remove("\\.sqlite")
  ) %>%
  dplyr::filter(year == 0)

# volume: use whichever column your plot_variables table actually carries; check names first
if (!"volume_m3" %in% names(sim_t0)) {
  message("plot_variables_all_FID_V1_FULL has no 'volume_m3' column at plot level -- ",
          "available columns are printed below, adjust pred_col accordingly:")
  print(names(sim_t0))
}

if ("volume_m3" %in% names(sim_t0)) {
  vol_sim_vs_fid1ha <- metrics_one_var(
    pred_df = sim_t0 %>% dplyr::select(plotid, volume_m3),
    obs_df  = fid1ha_volume,
    pred_col = "volume_m3", obs_col = "volume_m3ha_FID1ha_replicated",
    label = "sim_FID_V1_t0_vs_FID1ha_replicated", var_label = "volume_m3ha"
  )
  vol_sim_vs_fid1ha$summary
  
  ggplot2::ggplot(vol_sim_vs_fid1ha$table, ggplot2::aes(x = obs, y = pred)) +
    ggplot2::geom_point(alpha = 0.4) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    ggplot2::labs(
      title = "Volume at t0: simulation (FID_V1) vs FID2015 1ha-replicated init data",
      x = "Volume m3/ha - FID2015 1ha replicated (source)", y = "Volume m3/ha - simulation t0"
    ) +
    ggplot2::theme_bw()
}

#===============================================================================
# PART E -- AGE CHECK: FID2015 plot-description age vs simulation FID2015 age (t0)
#===============================================================================

# Confirmed field-age column: "domspeciesage2015" (age of the dominant species, plot level).
# Kept as a variable + fallback detection rather than hardcoded inline, in case naming changes.
age_col_2015 <- "domspeciesage2015"
if (!age_col_2015 %in% names(FID_2015_2019_Info_clean)) {
  age_candidates <- names(FID_2015_2019_Info_clean)[
    grepl("age", names(FID_2015_2019_Info_clean), ignore.case = TRUE) &
      grepl("2015", names(FID_2015_2019_Info_clean))
  ]
  if (length(age_candidates) == 0) {
    stop("'domspeciesage2015' not found and no fallback age-2015 column detected. Columns available:\n",
         paste(names(FID_2015_2019_Info_clean), collapse = ", "))
  }
  age_col_2015 <- age_candidates[1]
  message("'domspeciesage2015' not found -- falling back to '", age_col_2015, "'.")
}

message("Using '", age_col_2015, "' as the FID2015 plot-description age column.")

fid_age_2015 <- FID_2015_2019_Info_clean %>%
  dplyr::select(plotid, field_age = dplyr::all_of(age_col_2015)) %>%
  dplyr::filter(!is.na(field_age))

if (!"age" %in% names(sim_t0)) {
  stop("plot_variables_all_FID_V1_FULL has no 'age' column at year == 0 -- cannot run age check.")
}

sim_age_2015 <- sim_t0 %>%
  dplyr::select(plotid, sim_age = age) %>%
  dplyr::filter(!is.na(sim_age))

age_compare_result <- metrics_one_var(
  pred_df = sim_age_2015, obs_df = fid_age_2015,
  pred_col = "sim_age", obs_col = "field_age",
  label = "sim_FID_V1_t0_age_vs_FID2015_field_age", var_label = "stand_age"
)
age_compare_result$summary

# --- Numeric: full per-plot table for inspection ---
age_compare_result$table

# --- Visual: scatter with 1:1 line ---
ggplot2::ggplot(age_compare_result$table, ggplot2::aes(x = obs, y = pred)) +
  ggplot2::geom_point(alpha = 0.5) +
  ggplot2::geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  ggplot2::labs(
    title = "Stand age: FID2015 field description vs simulation t0 (FID_V1)",
    x = "Field-assessed age (years)", y = "Simulated age at t0 (years)"
  ) +
  ggplot2::theme_bw()

# --- Visual: distribution of the difference (sim - field) ---
ggplot2::ggplot(age_compare_result$table, ggplot2::aes(x = pred - obs)) +
  ggplot2::geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  ggplot2::labs(
    title = "Age difference (simulation t0 - FID2015 field age)",
    x = "Difference in years", y = "N plots"
  ) +
  ggplot2::theme_bw()

# --- Visual: paired boxplot, field vs sim, sample of plots only (readability) ---
age_sample_plots <- sample(age_compare_result$table$plotid, min(30, nrow(age_compare_result$table)))

age_compare_result$table %>%
  dplyr::filter(plotid %in% age_sample_plots) %>%
  tidyr::pivot_longer(cols = c(obs, pred), names_to = "source", values_to = "age") %>%
  dplyr::mutate(source = dplyr::recode(source, obs = "FID2015_field", pred = "sim_FID_V1_t0")) %>%
  ggplot2::ggplot(ggplot2::aes(x = source, y = age, group = plotid)) +
  ggplot2::geom_line(alpha = 0.3) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::labs(title = "Age pairing, 30-plot sample: FID2015 field age vs simulation t0",
                x = NULL, y = "Age (years)") +
  ggplot2::theme_bw()

# ---- Close the PDF if open
STOP!
print("CARFUL NOT RUN THE NEXT")  
#---
  
dev.off()

#===============================================================================
# --- SAVE SUMMARY TABLES ----
#===============================================================================

final_check_summary <- dplyr::bind_rows(
  trees_footprint_vs$summary,
  ba_footprint_vs$summary,
  vol_1ha_vs_500$summary,
  trees_1ha_vs_ALS$summary,
  if (exists("vol_sim_vs_fid1ha")) vol_sim_vs_fid1ha$summary,
  age_compare_result$summary
)

final_check_summary

write.csv(final_check_summary, file.path(table_path, "07_FID_ALS_1ha_volume_age_check_summary.csv"), row.names = FALSE)
saveRDS(counts_compare,  file.path(table_path, "07_counts_compare.rds"))
saveRDS(volume_compare,  file.path(table_path, "07_volume_compare.rds"))
saveRDS(age_compare_result$table, file.path(table_path, "07_age_compare_table.rds"))
