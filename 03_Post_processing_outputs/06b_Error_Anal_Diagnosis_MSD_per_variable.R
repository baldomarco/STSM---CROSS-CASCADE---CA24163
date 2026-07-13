library(tidyverse)

species_to_remove <- c("piab", "pisy", "abal", "lade", "psme", "pini", "pice")

#===============================================================================
# Single-variable metrics function — one variable, one comparison, at a time
#===============================================================================

metrics_one_var <- function(pred_df, obs_df, var, label){

  joined <- dplyr::inner_join(
    pred_df %>% dplyr::select(plotid, year, pred = dplyr::all_of(var)),
    obs_df  %>% dplyr::select(plotid, year, obs  = dplyr::all_of(var)),
    by = c("plotid", "year")
  )

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
# 1) ba_broadl — inspect plot-level values before trusting the summary
#===============================================================================

ba_ALS_vs_FID2015 <- metrics_one_var(
  pred_df = real_df %>% dplyr::filter(source == "real_ALS", year == 2015),
  obs_df  = real_df %>% dplyr::filter(source == "real_FID", year == 2015),
  var = "ba_broadl", label = "ALS2015_vs_FID2015"
)
ba_ALS_vs_FID2015$table    # look here first: pred vs obs, plot by plot
ba_ALS_vs_FID2015$summary


#===============================================================================
# 2) trees_10_40 — inspect separately (this is the one that blew up to ~100000%)
#===============================================================================

trees_ALS_vs_FID2015 <- metrics_one_var(
  pred_df = real_df %>% dplyr::filter(source == "real_ALS", year == 2015),
  obs_df  = real_df %>% dplyr::filter(source == "real_FID", year == 2015),
  var = "trees_10_40", label = "ALS2015_vs_FID2015"
)
trees_ALS_vs_FID2015$table
trees_ALS_vs_FID2015$summary


#===============================================================================
# 3) broadl_40 — inspect separately (unstable when obs is 0/1 per plot)
#===============================================================================

broadl40_ALS_vs_FID2015 <- metrics_one_var(
  pred_df = real_df %>% dplyr::filter(source == "real_ALS", year == 2015),
  obs_df  = real_df %>% dplyr::filter(source == "real_FID", year == 2015),
  var = "broadl_40", label = "ALS2015_vs_FID2015"
)
broadl40_ALS_vs_FID2015$table
broadl40_ALS_vs_FID2015$summary


#===============================================================================
# 4) volume_m3 — real FID (from tree-level v, x20) vs sim (aggregated landscape)
#    ALS excluded here: no volume formula confirmed yet for ALS tree list
#===============================================================================

real_volume <- FID_2015_clean_alive %>%
  dplyr::filter(!is.na(v)) %>%
  dplyr::group_by(plotid) %>%
  dplyr::summarise(volume_m3 = sum(v, na.rm = TRUE) * 20, .groups = "drop") %>%
  dplyr::mutate(year = 2015)

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

volume_simALSV1_vs_FID2015 <- metrics_one_var(
  pred_df = sim_volume %>% dplyr::filter(scenario == "sim_ALS_V1", year == 2015),
  obs_df  = real_volume,
  var = "volume_m3", label = "sim_ALS_V1_vs_FID2015"
)
volume_simALSV1_vs_FID2015$table
volume_simALSV1_vs_FID2015$summary


#===============================================================================
# Only after each block above looks sane: combine into one table
#===============================================================================

# results <- dplyr::bind_rows(
#   ba_ALS_vs_FID2015$summary,
#   trees_ALS_vs_FID2015$summary,
#   broadl40_ALS_vs_FID2015$summary,
#   volume_simALSV1_vs_FID2015$summary
# )
