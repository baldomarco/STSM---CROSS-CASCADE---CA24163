#   #
#   #
#   #
#   #

#===============================================================================
# --- START ----
#===============================================================================

# Clean the environment from previous analysed data
rm(list=ls()) 

library(tidyverse)

#===============================================================================
# --- UPLOAD DATA ----
#===============================================================================

lnd_scen_ALS_V1_FULL <- readRDS("C:/iLand/2026/Bialowieza/RDS/ALS_V1_Full/lnd_scen_ALS_V1_FULL.rds")
lnd_scen_ALS_V2_FULL <- readRDS("C:/iLand/2026/Bialowieza/RDS/ALS_V2_Full/lnd_scen_ALS_V2_FULL.rds")
lnd_scen_FID_V1_FULL <- readRDS("C:/iLand/2026/Bialowieza/RDS/FID_V1_Full/lnd_scen_FID_V1_FULL.rds")
lnd_scen_FID_V2_FULL <- readRDS("C:/iLand/2026/Bialowieza/RDS/FID_V2_Full/lnd_scen_FID_V2_FULL.rds")
ALS_clean_alive <-      readRDS("C:/P/DMP_CROSS_CASCADE/03_rawdata/02_process_storage/ALS_order_rds/ALS_clean_alive.rds")
FID_2019_clean_alive <- readRDS("C:/P/DMP_CROSS_CASCADE/03_rawdata/02_process_storage/FID_order_rds/FID_2019_clean_alive.rds")
FID_2015_clean_alive <- readRDS("C:/P/DMP_CROSS_CASCADE/03_rawdata/02_process_storage/FID_order_rds/FID_2015_clean_alive.rds")
plot_variables_all_ALS_V1_FULL <- readRDS("C:/iLand/2026/Bialowieza/RDS/ALS_V1_Full/plot_variables_all_ALS_V1_FULL.rds")
plot_variables_all_ALS_V2_FULL <- readRDS("C:/iLand/2026/Bialowieza/RDS/ALS_V2_Full/plot_variables_all_ALS_V2_FULL.rds")
plot_variables_all_FID_V1_FULL <- readRDS("C:/iLand/2026/Bialowieza/RDS/FID_V1_Full/plot_variables_all_FID_V1_FULL.rds")
plot_variables_all_FID_V2_FULL <- readRDS("C:/iLand/2026/Bialowieza/RDS/FID_V2_Full/plot_variables_all_FID_V2_FULL.rds")

#===============================================================================
# --- PART A - SETTINGS ----
#===============================================================================

species_to_remove <- c(
  "piab", "pisy", "abal",
  "lade", "psme", "pini",
  "pice"
)

labs <- c(
  lai_sim     = "LAI",
  ba_broadl   = "Broadleaf BA (m²/ha)",
  trees_10_40 = "Trees (10–40 cm)",
  broadl_40   = "Broadleaf >40 cm",
  age         = "Stand age"
)

# Ground-truth data cannot supply lai_sim/age (no LAI or age field at tree
# level in FID/ALS), so MSD comparisons vs real data are restricted to these:
vars_ground_truth <- c("ba_broadl", "trees_10_40", "broadl_40")


#===============================================================================
# --- FUNCTION: TREE LEVEL INVENTORY -> PLOT LEVEL (REAL DATA) ----
# NOTE: dbh in FID_*/ALS_clean_alive is in millimetres -> convert to cm first
#===============================================================================

process_inventory <- function(df, year, source, scale = 1){
  
  df <- df %>%
    dplyr::mutate(
      dbh_cm        = dbh / 10,
      basal_area_m2 = pi * (dbh_cm / 200)^2
    )
  
  df %>%
    dplyr::group_by(plotid) %>%
    dplyr::summarise(
      
      ba_broadl =
        sum(basal_area_m2[!species %in% species_to_remove], na.rm = TRUE),
      
      trees_10_40 =
        sum(dbh_cm >= 10 & dbh_cm <= 40, na.rm = TRUE),
      
      broadl_40 =
        sum(dbh_cm > 40 & !species %in% species_to_remove, na.rm = TRUE),
      
      .groups = "drop"
      
    ) %>%
    dplyr::mutate(
      year   = year,
      source = source,
      dplyr::across(
        c(ba_broadl, trees_10_40, broadl_40),
        ~ .x * scale
      )
    )
}


#===============================================================================
# REAL DATA -- built dynamically from whatever FID_*/ALS objects exist
# FID plots are 500 m2 -> scale x20 to 1 ha; ALS is already 1 ha -> scale x1
#===============================================================================

real_sources <- list(
  list(df = FID_2015_clean_alive, year = 2015, source = "real_FID", scale = 20),
  list(df = FID_2019_clean_alive, year = 2019, source = "real_FID", scale = 20)
)

if (exists("FID_2022_clean_alive")) {
  real_sources <- append(
    real_sources,
    list(list(df = FID_2022_clean_alive, year = 2022, source = "real_FID", scale = 20))
  )
}

real_sources <- append(
  real_sources,
  list(list(df = ALS_clean_alive, year = 2015, source = "real_ALS", scale = 1))
)

real_df <- purrr::map_dfr(
  real_sources,
  ~ process_inventory(.x$df, .x$year, .x$source, .x$scale)
)

real_fid_years <- real_df %>%
  dplyr::filter(source == "real_FID") %>%
  dplyr::distinct(year) %>%
  dplyr::pull(year) %>%
  sort()


#===============================================================================
# --- SIMULATION DATA ----
# raw "year" is an offset from 2015 -> converted to calendar year for joins
#===============================================================================

sim_df <- dplyr::bind_rows(
  
  plot_variables_all_ALS_V1_FULL %>% dplyr::mutate(scenario = "sim_ALS_V1"),
  plot_variables_all_ALS_V2_FULL %>% dplyr::mutate(scenario = "sim_ALS_V2"),
  plot_variables_all_FID_V1_FULL %>% dplyr::mutate(scenario = "sim_FID_V1"),
  plot_variables_all_FID_V2_FULL %>% dplyr::mutate(scenario = "sim_FID_V2")
  
) %>%
  dplyr::mutate(
    plotid = stringr::str_remove(run, "DB_") %>%
      stringr::str_remove("\\.sqlite"),
    year = year + 2015
  )

sim_scenarios <- unique(sim_df$scenario)


#===============================================================================
# --- KEEP ONLY PLOTS PRESENT IN EVERY REAL SAMPLING YEAR + EVERY SIM SCENARIO ----
# (built dynamically so a missing FID_2022 etc. no longer breaks the pipeline)
#===============================================================================

fid_plot_sets <- purrr::map(
  real_fid_years,
  ~ real_df %>% dplyr::filter(source == "real_FID", year == .x) %>% dplyr::pull(plotid)
)

als_plots <- real_df %>% dplyr::filter(source == "real_ALS") %>% dplyr::pull(plotid)

sim_plot_sets <- purrr::map(
  sim_scenarios,
  ~ sim_df %>% dplyr::filter(scenario == .x) %>% dplyr::pull(plotid)
)

common_plots <- purrr::reduce(
  c(fid_plot_sets, list(als_plots), sim_plot_sets),
  intersect
)

real_df <- real_df %>% dplyr::filter(plotid %in% common_plots)
sim_df  <- sim_df  %>% dplyr::filter(plotid %in% common_plots)


#===============================================================================
# --- SINGLE COMBINED TABLE, ONE COLUMN PER RUN/SCENARIO ----
#===============================================================================

master_long <- dplyr::bind_rows(
  
  sim_df %>%
    dplyr::mutate(run = scenario) %>%
    dplyr::select(plotid, year, run, lai_sim, ba_broadl, trees_10_40, broadl_40, age),
  
  real_df %>%
    dplyr::mutate(run = source) %>%
    dplyr::select(plotid, year, run, ba_broadl, trees_10_40, broadl_40)
  
) %>%
  tidyr::pivot_longer(
    cols      = c(lai_sim, ba_broadl, trees_10_40, broadl_40, age),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  dplyr::filter(!is.na(value))

master_wide <- master_long %>%
  tidyr::pivot_wider(names_from = run, values_from = value)




#===============================================================================
# ---- MSD DECOMPOSITION & RELATIVE DIFFERENCES -----
#===============================================================================

metrics_by_plot <- function(pred_df, obs_df, vars_used){
  
  joined <- dplyr::inner_join(
    pred_df, obs_df,
    by = c("plotid", "year"),
    suffix = c("_pred", "_obs")
  )
  
  if (nrow(joined) == 0) return(tibble::tibble())
  
  purrr::map_dfr(vars_used, function(v){
    
    per_plot <- joined %>%
      dplyr::transmute(
        plotid,
        pred = .data[[paste0(v, "_pred")]],
        obs  = .data[[paste0(v, "_obs")]],
        rel_diff_pct = (pred - obs) / (obs + 0.001) * 100,
        sq_diff      = (pred - obs)^2
      )
    
    tibble::tibble(
      variable            = v,
      n_plots             = nrow(per_plot),
      rel_difference_pct  = mean(per_plot$rel_diff_pct, na.rm = TRUE),
      rel_difference_sd   = sd(per_plot$rel_diff_pct, na.rm = TRUE),
      MSD                 = mean(per_plot$sq_diff, na.rm = TRUE),
      SB                  = (mean(per_plot$pred, na.rm = TRUE) - mean(per_plot$obs, na.rm = TRUE))^2,
      SDSD                = (sd(per_plot$pred, na.rm = TRUE) - sd(per_plot$obs, na.rm = TRUE))^2,
      LC                  = {
        r <- tryCatch(cor(per_plot$pred, per_plot$obs, use = "complete.obs"),
                      error = function(e) NA_real_)
        if (is.na(r)) r <- 0
        2 * sd(per_plot$pred, na.rm = TRUE) * sd(per_plot$obs, na.rm = TRUE) * (1 - r)
      }
    )
  })
}


#===============================================================================
# COMPARISON FUNCTION (guards against empty joins)
#===============================================================================

#compare_sets <- function(pred, obs, label, vars_used = vars_ground_truth){
  
  joined <- dplyr::inner_join(
    pred, obs,
    by = c("plotid", "year"),
    suffix = c("_pred", "_obs")
  )
  
  if (nrow(joined) == 0) {
    warning(paste0("No overlapping plots for: ", label))
    return(tibble::tibble())
  }
  
  purrr::map_dfr(
    vars_used,
    function(v){
      metrics(joined[[paste0(v, "_pred")]], joined[[paste0(v, "_obs")]]) %>%
        dplyr::mutate(comparison = label, variable = v)
    }
  )
}


compare_sets <- function(pred, obs, label, vars_used = vars_ground_truth){
  metrics_by_plot(pred, obs, vars_used) %>%
    dplyr::mutate(comparison = label, .before = 1)
}

#===============================================================================
# FINAL COMPARISONS
# 1) ALS2015 vs FID2015 (real vs real)
# 2) each sim scenario vs real FID, at every real sampling year (temporal)
#===============================================================================

als_vs_fid <- compare_sets(
  pred  = real_df %>% dplyr::filter(source == "real_ALS", year == 2015),
  obs   = real_df %>% dplyr::filter(source == "real_FID", year == 2015),
  label = "ALS2015_vs_FID2015"
)

sim_vs_fid <- purrr::map_dfr(sim_scenarios, function(sc){
  purrr::map_dfr(real_fid_years, function(yr){
    compare_sets(
      pred  = sim_df  %>% dplyr::filter(scenario == sc, year == yr),
      obs   = real_df %>% dplyr::filter(source == "real_FID", year == yr),
      label = paste0(sc, "_vs_FID", yr)
    )
  })
})

results <- dplyr::bind_rows(als_vs_fid, sim_vs_fid)


#===============================================================================
# FINAL OUTPUT
#===============================================================================

final_summary_table <- results %>%
  dplyr::select(
    comparison, 
    variable, 
    n_plots,
    rel_difference_pct, 
    rel_difference_sd,
    MSD, 
    SB, 
    SDSD, 
    LC)

final_summary_table

#===============================================================================




#===============================================================================
# PART B - DIAGNOSIS TEST
#===============================================================================


#===============================================================================
# --- Single-variable metrics function — one variable, one comparison, at a time ----
#===============================================================================

metrics_one_var <- function(pred_df, obs_df, var, label){
  
  pred_sel <- pred_df %>% dplyr::select(plotid, pred = dplyr::all_of(var))
  obs_sel  <- obs_df  %>% dplyr::select(plotid, obs  = dplyr::all_of(var))
  
  # HARD GUARD: each source must have exactly one row per plot for this
  # comparison, otherwise inner_join fans out and starts pairing values
  # across different plots instead of matching the same plot.
  if (anyDuplicated(pred_sel$plotid) != 0) {
    stop(paste0(label, " / ", var, ": pred_df has duplicate plotid rows — ",
                "aggregate to one row per plot before comparing."))
  }
  if (anyDuplicated(obs_sel$plotid) != 0) {
    stop(paste0(label, " / ", var, ": obs_df has duplicate plotid rows — ",
                "aggregate to one row per plot before comparing."))
  }
  
  joined <- dplyr::inner_join(pred_sel, obs_sel, by = "plotid")
  
  if (nrow(joined) == 0) {
    warning(paste0("No overlap for ", label, " / ", var))
    return(list(table = joined, summary = tibble::tibble()))
  }
  
  per_plot <- joined %>%
    dplyr::mutate(
      rel_diff_pct = (pred - obs) / (obs + 0.001) * 100,
      sq_diff      = (pred - obs)^2
    )
  
  r <- tryCatch(cor(per_plot$pred, per_plot$obs, use = "complete.obs"),
                error = function(e) NA_real_)
  if (is.na(r)) r <- 0
  
  summary_row <- tibble::tibble(
    comparison         = label,
    variable           = var,
    n_plots            = nrow(per_plot),
    rel_difference_pct = mean(per_plot$rel_diff_pct, na.rm = TRUE),
    rel_difference_sd  = sd(per_plot$rel_diff_pct, na.rm = TRUE),
    MSD                = mean(per_plot$sq_diff, na.rm = TRUE),
    SB                 = (mean(per_plot$pred, na.rm = TRUE) - mean(per_plot$obs, na.rm = TRUE))^2,
    SDSD               = (sd(per_plot$pred, na.rm = TRUE) - sd(per_plot$obs, na.rm = TRUE))^2,
    LC                 = 2 * sd(per_plot$pred, na.rm = TRUE) * sd(per_plot$obs, na.rm = TRUE) * (1 - r)
  )
  
  list(table = per_plot, summary = summary_row)
}


#===============================================================================
# --- 1) ba_broadl — inspect plot-level values before trusting the summary -----
#===============================================================================

ba_ALS_vs_FID2015 <- metrics_one_var(
  pred_df = real_df %>% dplyr::filter(source == "real_ALS", year == 2015),
  obs_df  = real_df %>% dplyr::filter(source == "real_FID", year == 2015),
  var = "ba_broadl", label = "ALS2015_vs_FID2015"
)

A <- ba_ALS_vs_FID2015$table    # look here first: pred vs obs, plot by plot
ba_ALS_vs_FID2015$summary


#===============================================================================
# --- 2) trees_10_40 — inspect separately (this is the one that blew up to ~100000%) ----
#===============================================================================

trees_ALS_vs_FID2015 <- metrics_one_var(
  pred_df = real_df %>% dplyr::filter(source == "real_ALS", year == 2015),
  obs_df  = real_df %>% dplyr::filter(source == "real_FID", year == 2015),
  var = "trees_10_40", label = "ALS2015_vs_FID2015"
)
trees_ALS_vs_FID2015$table
trees_ALS_vs_FID2015$summary


#===============================================================================
# --- 3) broadl_40 — inspect separately (unstable when obs is 0/1 per plot) ----
#===============================================================================

broadl40_ALS_vs_FID2015 <- metrics_one_var(
  pred_df = real_df %>% dplyr::filter(source == "real_ALS", year == 2015),
  obs_df  = real_df %>% dplyr::filter(source == "real_FID", year == 2015),
  var = "broadl_40", label = "ALS2015_vs_FID2015"
)
broadl40_ALS_vs_FID2015$table
broadl40_ALS_vs_FID2015$summary


#===============================================================================
# --- 4) volume_m3 — real FID (from tree-level v, x20) vs each sim scenario, ----
#    at every available real sampling year. ALS excluded: no volume formula
#    confirmed yet for ALS tree list.
#===============================================================================

real_volume_sources <- list(
  list(df = FID_2015_clean_alive, year = 2015)
)
if (exists("FID_2019_clean_alive")) {
  real_volume_sources <- append(real_volume_sources, list(list(df = FID_2019_clean_alive, year = 2019)))
}
if (exists("FID_2022_clean_alive")) {
  real_volume_sources <- append(real_volume_sources, list(list(df = FID_2022_clean_alive, year = 2022)))
}

real_volume_all <- purrr::map_dfr(real_volume_sources, function(src){
  src$df %>%
    dplyr::filter(!is.na(v)) %>%
    dplyr::group_by(plotid) %>%
    dplyr::summarise(volume_m3 = sum(v, na.rm = TRUE) * 20, .groups = "drop") %>%
    dplyr::mutate(year = src$year)
})

sim_volume <- dplyr::bind_rows(
  lnd_scen_ALS_V1_FULL %>% dplyr::mutate(scenario = "sim_ALS_V1"),
  lnd_scen_ALS_V2_FULL %>% dplyr::mutate(scenario = "sim_ALS_V2"),
  lnd_scen_FID_V1_FULL %>% dplyr::mutate(scenario = "sim_FID_V1"),
  lnd_scen_FID_V2_FULL %>% dplyr::mutate(scenario = "sim_FID_V2")
) %>%
  dplyr::mutate(
    plotid = stringr::str_remove(run, "DB_") %>% stringr::str_remove("\\.sqlite"),
    year   = year + 2015
  ) %>%
  dplyr::group_by(plotid, year, scenario) %>%
  dplyr::summarise(volume_m3 = sum(volume_m3, na.rm = TRUE), .groups = "drop")

volume_results <- purrr::map_dfr(unique(sim_volume$scenario), function(sc){
  purrr::map_dfr(unique(real_volume_all$year), function(yr){
    metrics_one_var(
      pred_df = sim_volume      %>% dplyr::filter(scenario == sc, year == yr),
      obs_df  = real_volume_all %>% dplyr::filter(year == yr),
      var = "volume_m3",
      label = paste0(sc, "_vs_FID", yr)
    )$summary
  })
})

volume_results   # relative difference + full MSD decomposition, volume only


#===============================================================================
# Only after each block above looks sane: combine into one table
#===============================================================================

# results <- dplyr::bind_rows(
#   ba_ALS_vs_FID2015$summary,
#   trees_ALS_vs_FID2015$summary,
#   broadl40_ALS_vs_FID2015$summary,
#   volume_simALSV1_vs_FID2015$summary
# )

#===============================================================================




#===============================================================================
# --- PART C — SPECIES x VOLUME, PLOT BY PLOT: FID2015 (real) vs 4 sim scenarios ----
#===============================================================================

species_to_remove <- c("piab", "pisy", "abal", "lade", "psme", "pini", "pice")

#-------------------------------------------------------------------------------
# Real FID2015: tree-level volume (v, m3) by species, scaled 500m2 -> 1ha (x20)
#-------------------------------------------------------------------------------

real_vol_species <- FID_2015_clean_alive %>%
  dplyr::filter(!is.na(v)) %>%
  dplyr::group_by(plotid, species) %>%
  dplyr::summarise(volume_m3 = sum(v, na.rm = TRUE) * 20, .groups = "drop") %>%
  dplyr::mutate(run = "real_FID2015")

#-------------------------------------------------------------------------------
# Sim landscape tables: already per-ha, species x year x run
# year 0 = 2015 -> keep year == 0 to match FID2015
#-------------------------------------------------------------------------------

sim_vol_species <- dplyr::bind_rows(
  lnd_scen_ALS_V1_FULL %>% dplyr::mutate(run = "sim_ALS_V1"),
  lnd_scen_ALS_V2_FULL %>% dplyr::mutate(run = "sim_ALS_V2"),
  lnd_scen_FID_V1_FULL %>% dplyr::mutate(run = "sim_FID_V1"),
  lnd_scen_FID_V2_FULL %>% dplyr::mutate(run = "sim_FID_V2")
) %>%
  dplyr::filter(year == 0) %>%
  dplyr::mutate(
    plotid = stringr::str_remove(run, "DB_") # placeholder, corrected below
  )

# NOTE: in this table "run" already got overwritten by scenario name above,
# so plotid must come from the ORIGINAL sqlite filename column instead.
# Re-derive properly:

sim_vol_species <- dplyr::bind_rows(
  lnd_scen_ALS_V1_FULL %>% dplyr::mutate(scenario = "sim_ALS_V1"),
  lnd_scen_ALS_V2_FULL %>% dplyr::mutate(scenario = "sim_ALS_V2"),
  lnd_scen_FID_V1_FULL %>% dplyr::mutate(scenario = "sim_FID_V1"),
  lnd_scen_FID_V2_FULL %>% dplyr::mutate(scenario = "sim_FID_V2")
) %>%
  dplyr::filter(year == 0) %>%
  dplyr::mutate(
    plotid = stringr::str_remove(run, "DB_") %>% stringr::str_remove("\\.sqlite")
  ) %>%
  dplyr::select(plotid, species, volume_m3, run = scenario)

#-------------------------------------------------------------------------------
# Restrict to plots common to real + all 4 sims (reuse common_plots if present,
# else recompute from what's available here)
#-------------------------------------------------------------------------------

vol_common_plots <- purrr::reduce(
  list(
    unique(real_vol_species$plotid),
    sim_vol_species %>% dplyr::filter(run == "sim_ALS_V1") %>% dplyr::pull(plotid) %>% unique(),
    sim_vol_species %>% dplyr::filter(run == "sim_ALS_V2") %>% dplyr::pull(plotid) %>% unique(),
    sim_vol_species %>% dplyr::filter(run == "sim_FID_V1") %>% dplyr::pull(plotid) %>% unique(),
    sim_vol_species %>% dplyr::filter(run == "sim_FID_V2") %>% dplyr::pull(plotid) %>% unique()
  ),
  intersect
)

vol_species_all <- dplyr::bind_rows(real_vol_species, sim_vol_species) %>%
  dplyr::filter(plotid %in% vol_common_plots) %>%
  dplyr::mutate(
    species_group = ifelse(species %in% species_to_remove, species, "broadleaf_other")
  )

#-------------------------------------------------------------------------------
# --- Real FID volume, all years (2015, 2019, 2022 if present) ----
#-------------------------------------------------------------------------------

real_vol_species <- purrr::map_dfr(real_volume_sources, function(src){
  src$df %>%
    dplyr::filter(!is.na(v)) %>%
    dplyr::group_by(plotid, species) %>%
    dplyr::summarise(volume_m3 = sum(v, na.rm = TRUE) * 20, .groups = "drop") %>%
    dplyr::mutate(year = src$year, run = paste0("real_FID", src$year))
})

#-------------------------------------------------------------------------------
# Sim species volume, all years matching real FID years
#-------------------------------------------------------------------------------

sim_vol_species <- dplyr::bind_rows(
  lnd_scen_ALS_V1_FULL %>% dplyr::mutate(scenario = "sim_ALS_V1"),
  lnd_scen_ALS_V2_FULL %>% dplyr::mutate(scenario = "sim_ALS_V2"),
  lnd_scen_FID_V1_FULL %>% dplyr::mutate(scenario = "sim_FID_V1"),
  lnd_scen_FID_V2_FULL %>% dplyr::mutate(scenario = "sim_FID_V2")
) %>%
  dplyr::mutate(
    plotid = stringr::str_remove(run, "DB_") %>% stringr::str_remove("\\.sqlite"),
    year   = year + 2015
  ) %>%
  dplyr::filter(year %in% unique(real_vol_species$year)) %>%
  dplyr::select(plotid, species, volume_m3, year, run = scenario)

vol_species_all <- dplyr::bind_rows(
  real_vol_species %>% dplyr::select(plotid, species, volume_m3, year, run),
  sim_vol_species
)

vol_common_plots <- purrr::reduce(
  c(
    list(unique(real_vol_species$plotid)),
    purrr::map(unique(sim_vol_species$run), ~ sim_vol_species %>% dplyr::filter(run == .x) %>% dplyr::pull(plotid) %>% unique())
  ),
  intersect
)

vol_species_all <- vol_species_all %>%
  dplyr::filter(plotid %in% vol_common_plots)

#-------------------------------------------------------------------------------
# Species colors — restrict full palette/order to species actually present
#-------------------------------------------------------------------------------

cols.all <- c(
  "rops"="#e0e0e0", "acpl"="#A9A9A9",   "alin"="#696969", "alvi"="#2e2e2e",
  "bepe"="#fadfad", "bepu"= "#5F9EA0", "prse" = "#8B2323",
  "casa"="#7eeadf", "coav"="#20c6b6",
  "tipl"="#645394", "ulgl"="#311432",
  "saca"="#D8BFD8",  "soar"="#DDA0DD", "soau"="#BA55D3",
  "pice"="#D27D2D", "pini"="#a81c07",
  "algl"="#2ECBE9","tico"="#128FC8",  "potr"="#00468B","poni"="#5BAEB7",
  "frex"="#fe9cb5","cabe"="#fe6181","acps"="#fe223e",
  "lade"="#FFFE71","abal"="#FFD800", "pisy"="#A4DE02",
  "fasy"="#76BA1B", "piab"="#006600",
  "quro"="#FF7F00", "qupe"="#FF9900", "qupu"="#CC9900"
)

new_order_gg.all <- c(
  "alvi","alin", "acpl", "rops","bepe" ,"bepu", "prse", "coav", "casa", "ulgl", "tipl",  "soau", "soar", "saca",  "pini", "pice",
  "poni", "algl", "tico", "potr",  "frex","cabe", "acps",  "lade", "abal",  "qupu", "qupe","quro","pisy", "fasy", "piab"
)

species.we.have <- unique(vol_species_all$species)
cols            <- cols.all[names(cols.all) %in% species.we.have]
new_order_gg    <- new_order_gg.all[new_order_gg.all %in% species.we.have]

#-------------------------------------------------------------------------------
# Single plot only, explicit selection — now faceted by year (2015 + 2019)
#-------------------------------------------------------------------------------

example_id <- vol_common_plots[218]

vol_example <- vol_species_all %>%
  dplyr::filter(plotid == example_id) %>%
  dplyr::mutate(
    species = factor(species, levels = new_order_gg),
    year    = factor(year)
  )

ggplot2::ggplot(vol_example, ggplot2::aes(x = run, y = volume_m3, fill = species)) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(~ year, scales = "free_x") +
  ggplot2::scale_fill_manual(values = cols) +
  ggplot2::labs(
    title = paste("Species volume -", example_id),
    x = NULL, y = "Volume (m³/ha)", fill = "Species"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))








#===============================================================================
# ---- PART C-ALTERNATIVE — SPECIES x VOLUME, PLOT BY PLOT: FID2015 (real) vs 4 sim scenarios ----
#===============================================================================

species_to_remove <- c("piab", "pisy", "abal", "lade", "psme", "pini", "pice")

#-------------------------------------------------------------------------------
# Real FID2015: tree-level volume (v, m3) by species, scaled 500m2 -> 1ha (x20)
#-------------------------------------------------------------------------------

real_vol_species <- FID_2015_clean_alive %>%
  dplyr::filter(!is.na(v)) %>%
  dplyr::group_by(plotid, species) %>%
  dplyr::summarise(volume_m3 = sum(v, na.rm = TRUE) * 20, .groups = "drop") %>%
  dplyr::mutate(run = "real_FID2015")

#-------------------------------------------------------------------------------
# Sim landscape tables: already per-ha, species x year x run
# year 0 = 2015 -> keep year == 0 to match FID2015
#-------------------------------------------------------------------------------

sim_vol_species <- dplyr::bind_rows(
  lnd_scen_ALS_V1_FULL %>% dplyr::mutate(run = "sim_ALS_V1"),
  lnd_scen_ALS_V2_FULL %>% dplyr::mutate(run = "sim_ALS_V2"),
  lnd_scen_FID_V1_FULL %>% dplyr::mutate(run = "sim_FID_V1"),
  lnd_scen_FID_V2_FULL %>% dplyr::mutate(run = "sim_FID_V2")
) %>%
  dplyr::filter(year == 0) %>%
  dplyr::mutate(
    plotid = stringr::str_remove(run, "DB_") # placeholder, corrected below
  )

# NOTE: in this table "run" already got overwritten by scenario name above,
# so plotid must come from the ORIGINAL sqlite filename column instead.
# Re-derive properly:

sim_vol_species <- dplyr::bind_rows(
  lnd_scen_ALS_V1_FULL %>% dplyr::mutate(scenario = "sim_ALS_V1"),
  lnd_scen_ALS_V2_FULL %>% dplyr::mutate(scenario = "sim_ALS_V2"),
  lnd_scen_FID_V1_FULL %>% dplyr::mutate(scenario = "sim_FID_V1"),
  lnd_scen_FID_V2_FULL %>% dplyr::mutate(scenario = "sim_FID_V2")
) %>%
  dplyr::filter(year == 0) %>%
  dplyr::mutate(
    plotid = stringr::str_remove(run, "DB_") %>% stringr::str_remove("\\.sqlite")
  ) %>%
  dplyr::select(plotid, species, volume_m3, run = scenario)

#-------------------------------------------------------------------------------
# Restrict to plots common to real + all 4 sims (reuse common_plots if present,
# else recompute from what's available here)
#-------------------------------------------------------------------------------

vol_common_plots <- purrr::reduce(
  list(
    unique(real_vol_species$plotid),
    sim_vol_species %>% dplyr::filter(run == "sim_ALS_V1") %>% dplyr::pull(plotid) %>% unique(),
    sim_vol_species %>% dplyr::filter(run == "sim_ALS_V2") %>% dplyr::pull(plotid) %>% unique(),
    sim_vol_species %>% dplyr::filter(run == "sim_FID_V1") %>% dplyr::pull(plotid) %>% unique(),
    sim_vol_species %>% dplyr::filter(run == "sim_FID_V2") %>% dplyr::pull(plotid) %>% unique()
  ),
  intersect
)

vol_species_all <- dplyr::bind_rows(real_vol_species, sim_vol_species) %>%
  dplyr::filter(plotid %in% vol_common_plots) %>%
  dplyr::mutate(
    species_group = ifelse(species %in% species_to_remove, species, "broadleaf_other")
  )



#-------------------------------------------------------------------------------
# Single-plot stacked bar: one plotid at a time, bars = run, fill = species
#-------------------------------------------------------------------------------

plot_species_volume <- function(target_plotid){
  
  if (!target_plotid %in% vol_common_plots) {
    stop(paste0(target_plotid, " is not in vol_common_plots (missing from real or one of the sims)"))
  }
  
  vol_species_all %>%
    dplyr::filter(plotid == target_plotid) %>%
    ggplot2::ggplot(ggplot2::aes(x = run, y = volume_m3, fill = species)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(
      x = NULL, y = "Volume (m³/ha)", fill = "Species",
      title = paste0("Species volume — plot ", target_plotid, ": FID2015 (real) vs simulation scenarios")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))+
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      plot.background = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )
}

#---- Single Plot Volume example: sage: pick one plotid at a time, e.g. -----
plot_species_volume("KS004")
plot_species_volume(vol_common_plots[218])



#-------------------------------------------------------------------------------
# --- ALL PLOTS TOGETHER VISUALIZATION - Plot-by-plot bar chart ----
#-------------------------------------------------------------------------------

plot_ids_to_show <- sort(vol_common_plots)  # all common plots; subset if too many

p_species_volume <- vol_species_all %>%
  dplyr::filter(plotid %in% plot_ids_to_show) %>%
  ggplot2::ggplot(ggplot2::aes(x = run, y = volume_m3, fill = species)) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = cols) +
  ggplot2::facet_wrap(~ plotid, scales = "free_y") +
  ggplot2::labs(
    x = NULL, y = "Volume (m³/ha)", fill = "Species",
    title = "Species volume by plot: FID2015 (real) vs simulation scenarios"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

# Visualize all the plots in one graph
p_species_volume
