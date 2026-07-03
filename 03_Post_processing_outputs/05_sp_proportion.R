#===============================================================================
# Species proportion (basal area based) per PLOT and per SITE
# - Sites are not encoded in plotid (all plots are "KS###"), so site groups
#   are detected automatically from plot coordinates (spatial clustering).
# - If no spatial cluster structure is found, falls back to one site = whole
#   landscape.
# - Includes the "deciduous" -> most-common-broadleaf-in-plot correction.
#===============================================================================

library(dplyr)
library(tidyr)
library(writexl)
library(ggplot2)
library(dbscan)     # install.packages("dbscan") if needed

#------------------------------------------------------------------------------
# 0. Paths
#------------------------------------------------------------------------------
input_rds  <- "C:/iLand/20230901_Bottoms_Up/plot_init/R/stsm_roma/ALS_clean_alive.rds"
output_dir <- "C:/P/DMP_CROSS_CASCADE/03_rawdata/02_process_storage/"
plot_dir   <- file.path(output_dir, "plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)

#------------------------------------------------------------------------------
# 1. Load data
#------------------------------------------------------------------------------
ALS_clean_alive <- readRDS(input_rds)

#------------------------------------------------------------------------------
# 2. Deciduous trick — replace "deciduous" with the most common broadleaf
#    species found in the SAME plot (your original logic, unchanged)
#------------------------------------------------------------------------------
broadleaf_species <- setdiff(unique(ALS_clean_alive$species), c("piab", "pisy", "deciduous"))

broadleaf_mode <- ALS_clean_alive %>%
  dplyr::filter(species %in% broadleaf_species) %>%
  dplyr::group_by(plotid, species) %>%
  dplyr::tally() %>%
  dplyr::group_by(plotid) %>%
  dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
  dplyr::select(plotid, species)

ALS_fixed <- ALS_clean_alive %>%
  dplyr::left_join(broadleaf_mode, by = "plotid", suffix = c("", "_mode")) %>%
  dplyr::mutate(species = ifelse(species == "deciduous" & !is.na(species_mode),
                                 species_mode, species)) %>%
  dplyr::select(-species_mode)

# Sanity check: plots where "deciduous" was the ONLY broadleaf present have no
# mode to fall back on and stay tagged "deciduous" — flag these for manual review.
still_deciduous <- ALS_fixed %>% dplyr::filter(species == "deciduous")
if (nrow(still_deciduous) > 0) {
  message(nrow(still_deciduous), " trees remain tagged 'deciduous' in ",
          dplyr::n_distinct(still_deciduous$plotid),
          " plot(s) with no other broadleaf species to borrow from: ",
          paste(unique(still_deciduous$plotid), collapse = ", "))
}

#------------------------------------------------------------------------------
# 3. Basal area per tree
#    dbh is in mm -> convert to cm -> BA in m2 = pi/4 * (dbh_cm/100)^2
#------------------------------------------------------------------------------
ALS_fixed <- ALS_fixed %>%
  dplyr::mutate(dbh_cm = dbh / 10,
                ba_m2  = pi / 4 * (dbh_cm / 100)^2)

#------------------------------------------------------------------------------
# 4. Detect plot clusters ("sites") automatically from spatial coordinates
#------------------------------------------------------------------------------
plot_coords <- ALS_fixed %>%
  dplyr::group_by(plotid) %>%
  dplyr::summarise(x = mean(x, na.rm = TRUE), y = mean(y, na.rm = TRUE), .groups = "drop")

# --- Diagnostic: k-nearest-neighbour distance plot to help choose eps ---
# Look for the "elbow" where distances jump: that's the within-cluster vs.
# between-cluster gap. Re-run with a different eps_dist if your dataset's
# elbow sits somewhere else than the value used below.
png(file.path(output_dir, "kNNdist_diagnostic.png"), width = 800, height = 600)
dbscan::kNNdistplot(as.matrix(plot_coords[, c("x", "y")]), k = 1)
abline(h = 800, col = "red", lty = 2)
dev.off()

eps_dist <- 800   # metres — set from the elbow in kNNdist_diagnostic.png
min_pts  <- 2

clust <- dbscan::dbscan(as.matrix(plot_coords[, c("x", "y")]), eps = eps_dist, minPts = min_pts)
plot_coords$cluster <- clust$cluster   # 0 = noise / isolated plot

n_clusters <- length(unique(clust$cluster[clust$cluster != 0]))
n_noise    <- sum(clust$cluster == 0)

if (n_clusters == 0) {
  # No spatial cluster structure detected -> whole landscape is one site
  message("No spatial clusters detected — using the whole landscape as a single site.")
  plot_coords$Site <- "WHOLE_LANDSCAPE"
} else {
  message(n_clusters, " plot cluster(s) detected as site groups; ",
          n_noise, " isolated plot(s) kept as their own singleton site.")
  plot_coords <- plot_coords %>%
    dplyr::mutate(Site = ifelse(cluster == 0,
                                paste0("SITE_solo_", plotid),
                                paste0("SITE_", sprintf("%02d", cluster))))
}

ALS_fixed <- ALS_fixed %>%
  dplyr::left_join(plot_coords %>% dplyr::select(plotid, Site), by = "plotid")

#------------------------------------------------------------------------------
# 5. Species proportions per PLOT (basal-area based)
#------------------------------------------------------------------------------
plot_species_ba <- ALS_fixed %>%
  dplyr::group_by(Site, plotid, species) %>%
  dplyr::summarise(TotalBA = sum(ba_m2, na.rm = TRUE), .groups = "drop_last") %>%
  dplyr::mutate(SpeciesProportion = TotalBA / sum(TotalBA),
                sp_per = SpeciesProportion * 100) %>%
  dplyr::ungroup()

#------------------------------------------------------------------------------
# 6. Species proportions per SITE / cluster (basal-area based)
#------------------------------------------------------------------------------
site_species_ba <- ALS_fixed %>%
  dplyr::group_by(Site, species) %>%
  dplyr::summarise(TotalBA = sum(ba_m2, na.rm = TRUE), .groups = "drop_last") %>%
  dplyr::mutate(SpeciesProportion = TotalBA / sum(TotalBA),
                sp_per = SpeciesProportion * 100) %>%
  dplyr::ungroup()

#------------------------------------------------------------------------------
# 7. Export summary tables
#------------------------------------------------------------------------------
write_xlsx(plot_species_ba, file.path(output_dir, "sp_prop_plot_ba.xlsx"))
write_xlsx(site_species_ba, file.path(output_dir, "sp_prop_site_ba.xlsx"))
write_xlsx(plot_coords %>% dplyr::select(plotid, Site),
           file.path(output_dir, "plot_to_site_lookup.xlsx"))

#------------------------------------------------------------------------------
# 8. Per-plot xlsx export (tree list, deciduous already resolved) — same
#    logic as your original loop, now driven off ALS_fixed
#------------------------------------------------------------------------------
unique_plots <- unique(ALS_fixed$plotid)
for (p in unique_plots) {
  subset_df <- ALS_fixed %>%
    dplyr::filter(plotid == p) %>%
    tidyr::drop_na(x, y)
  
  if (nrow(subset_df) == 0) next
  
  subset_df$x <- as.numeric(subset_df$x)
  subset_df$y <- as.numeric(subset_df$y)
  
  write_xlsx(subset_df, file.path(plot_dir, paste0("plot_", p, ".xlsx")))
}

#------------------------------------------------------------------------------
# 9. Visualization — pie charts of BA-based species proportion per site
#------------------------------------------------------------------------------
species.we.have <- unique(ALS_fixed$species)

# Palette keyed on the 4-letter species codes actually used in this dataset.
# Colours reused from your original latin-name palette where the species
# matches; algl (Alnus glutinosa) is new so it gets a similar Alnus-family grey.
cols.all <- c(
  "algl" = "#696969",  # Alnus glutinosa
  "acpl" = "red",  # Acer platanoides
  "bepe" = "#fadfad",  # Betula pendula
  "potr" = "#00468B",  # Populus tremula
  "frex" = "#fe9cb5",  # Fraxinus excelsior
  "tico" = "#128FC8",  # Tilia cordata
  "quro" = "#FF7F00",  # Quercus robur
  "cabe" = "#fe6181",  # Carpinus betulus
  "pisy" = "#A4DE02",  # Pinus sylvestris
  "piab" = "#006600",  # Picea abies
  "deciduous" = "#e0e0e0"
)

# Ecological ordering: pioneers/broadleaf first, shade-tolerant conifers last
new_order_gg.all <- c("algl", "acpl", "bepe", "potr", "frex", "tico", "quro",
                      "cabe", "pisy", "piab", "deciduous")

cols          <- cols.all[names(cols.all) %in% species.we.have]
new_order_gg  <- new_order_gg.all[new_order_gg.all %in% species.we.have]

n_sites <- dplyr::n_distinct(site_species_ba$Site)

x7wb <- ggplot(site_species_ba,
               aes(x = "", y = SpeciesProportion,
                   fill = factor(species, levels = new_order_gg))) +
  geom_bar(stat = "identity", width = 1, show.legend = TRUE) +
  scale_fill_manual(values = cols[new_order_gg], guide = guide_legend(reverse = TRUE)) +
  facet_wrap(~Site, ncol = 28) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = ifelse(sp_per > 5, round(sp_per, 1), "")),
            position = position_stack(vjust = 0.5), size = 2) +
  labs(x = NULL, y = NULL, fill = NULL) +
  ggtitle(paste0("Species proportions [%] based on site-level basal area (",
                 n_sites, " sites detected)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 6),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(file.path(output_dir, "species_proportion_per_site_BA.png"),
       x7wb, width = 16, height = 16, dpi = 200)

x7wb
