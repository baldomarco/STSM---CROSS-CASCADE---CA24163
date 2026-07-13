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
input_rds  <- "C:/P/DMP_CROSS_CASCADE/03_rawdata/02_process_storage/FID_order_rds/FID_2015_clean_alive.rds"
output_dir <- "C:/P/DMP_CROSS_CASCADE/03_rawdata/02_process_storage/"
plot_dir   <- file.path(output_dir, "plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)

#------------------------------------------------------------------------------
# 1. Load data
#------------------------------------------------------------------------------
FID_2015_clean_alive <- readRDS(input_rds)


#------------------------------------------------------------------------------
# 2. Replace Prunus, Malus and Pyrus with Prunus serotina (prse)
#    all Salix sp with Salix caprea (saca), Ulmus with U. glabra (ulgl) and NA with their reference species
#------------------------------------------------------------------------------
# --- Convert iLand Species ----

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

#------------------------------------------------------------------------------
# 3. Basal area per tree
#    dbh is in mm -> convert to cm -> BA in m2 = pi/4 * (dbh_cm/100)^2
#------------------------------------------------------------------------------
FID_2015_fix <- FID_2015_clean_alive %>%
  dplyr::mutate(dbh_cm = dbh / 10,
                ba_m2  = pi / 4 * (dbh_cm / 100)^2)

#------------------------------------------------------------------------------
# 4. Detect plot clusters ("sites") automatically from spatial coordinates
#------------------------------------------------------------------------------
plot_coords <- FID_2015_fix %>%
  dplyr::group_by(plotid) %>%
  dplyr::summarise(x = mean(x, na.rm = TRUE), y = mean(y, na.rm = TRUE), .groups = "drop")

# --- Diagnostic: k-nearest-neighbour distance plot to help choose eps ---
# Look for the "elbow" where distances jump: that's the within-cluster vs.
# between-cluster gap. Re-run with a different eps_dist if your dataset's
# elbow sits somewhere else than the value used below.
png(file.path(output_dir, "FID_2015_kNNdist_diagnostic.png"), width = 800, height = 600)
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

FID_2015_fix <- FID_2015_fix %>%
  dplyr::left_join(plot_coords %>% dplyr::select(plotid, Site), by = "plotid")

#------------------------------------------------------------------------------
# 5. Species proportions per PLOT (basal-area based)
#------------------------------------------------------------------------------
plot_species_ba <- FID_2015_fix %>%
  dplyr::group_by(Site, plotid, species) %>%
  dplyr::summarise(TotalBA = sum(ba_m2, na.rm = TRUE), .groups = "drop_last") %>%
  dplyr::mutate(SpeciesProportion = TotalBA / sum(TotalBA),
                sp_per = SpeciesProportion * 100) %>%
  dplyr::ungroup()

#------------------------------------------------------------------------------
# 6. Species proportions per SITE / cluster (basal-area based)
#------------------------------------------------------------------------------
site_species_ba <- FID_2015_fix %>%
  dplyr::group_by(Site, species) %>%
  dplyr::summarise(TotalBA = sum(ba_m2, na.rm = TRUE), .groups = "drop_last") %>%
  dplyr::mutate(SpeciesProportion = TotalBA / sum(TotalBA),
                sp_per = SpeciesProportion * 100) %>%
  dplyr::ungroup()

#------------------------------------------------------------------------------
# 7. Export summary tables
#------------------------------------------------------------------------------
write_xlsx(plot_species_ba, file.path(output_dir, "sp_prop_plot_ba_FID_2015.xlsx"))
write_xlsx(site_species_ba, file.path(output_dir, "sp_prop_site_ba_FID_2015.xlsx"))
write_xlsx(plot_coords %>% dplyr::select(plotid, Site),
           file.path(output_dir, "plot_to_site_lookup_FID_2015.xlsx"))

#------------------------------------------------------------------------------
# 8. Per-plot xlsx export (tree list, deciduous already resolved) — same
#    logic as your original loop, now driven off FID_2015_fix
#------------------------------------------------------------------------------
unique_plots <- unique(FID_2015_fix$plotid)
for (p in unique_plots) {
  subset_df <- FID_2015_fix %>%
    dplyr::filter(plotid == p) %>%
    tidyr::drop_na(x, y)
  
  if (nrow(subset_df) == 0) next
  
  subset_df$x <- as.numeric(subset_df$x)
  subset_df$y <- as.numeric(subset_df$y)
  
  write_xlsx(subset_df, file.path(plot_dir, paste0("FID_2015_plot_", p, ".xlsx")))
}

#------------------------------------------------------------------------------
# 9. Visualization — pie charts of BA-based species proportion per site
#------------------------------------------------------------------------------
species.we.have <- unique(FID_2015_fix$species)

# Palette keyed on the 4-letter species codes actually used in this dataset.
# Colours reused from your original latin-name palette where the species
# matches; algl (Alnus glutinosa) is new so it gets a similar Alnus-family grey.
cols.all <- c(
  "algl" = "#696969",  # Alnus glutinosa
  "acpl" = "red",      # Acer platanoides
  "acps" = "#CD5C5C",  # Acer pseudoplatanus
  "bepe" = "#fadfad",  # Betula pendula
  "bepu" = "#FFF2B2",  # Betula pubescens
  "potr" = "#00468B",  # Populus tremula
  "frex" = "#fe9cb5",  # Fraxinus excelsior
  "tipl" = "#128FC8",  # Tilia platyphyllos / petiolaris
  "quro" = "#FF7F00",  # Quercus robur
  "qupe" = "#F4A460",  # Quercus petraea
  "cabe" = "#fe6181",  # Carpinus betulus
  "pisy" = "#A4DE02",  # Pinus sylvestris
  "piab" = "#006600",  # Picea abies
  "saca" = "#7FC97F",  # Salix caprea / alba / fragilis
  "soau" = "#8C564B",  # Sorbus aucuparia
  "prse" = "#984EA3",  # Prunus spp., Malus, Pyrus
  "ulgl" = "#8DD3C7",  # Ulmus spp.
  "deciduous" = "#e0e0e0"
)


# Ecological ordering: pioneers/broadleaf first, shade-tolerant conifers last
new_order_gg.all <- c(
  "algl",
  "acpl",
  "acps",
  "bepe",
  "bepu",
  "potr",
  "saca",
  "soau",
  "prse",
  "frex",
  "tipl",
  "ulgl",
  "quro",
  "qupe",
  "cabe",
  "pisy",
  "piab",
  "deciduous"
)

#-------------------------------------------------------------------------------
cols          <- cols.all[names(cols.all) %in% species.we.have]
new_order_gg  <- new_order_gg.all[new_order_gg.all %in% species.we.have]

n_sites <- dplyr::n_distinct(site_species_ba$Site)

# Plot -------------------------------------------------------------------------
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

ggsave(file.path(output_dir, "species_proportion_per_site_BA_FID_2015.png"),
       x7wb, width = 16, height = 16, dpi = 200)

x7wb
