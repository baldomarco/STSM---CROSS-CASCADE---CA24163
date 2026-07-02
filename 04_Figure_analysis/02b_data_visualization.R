# Marco Baldo 21-06-2026 contact: baldo@fld.czu.cz
# Data visualization, descriptive statistics & diagnostic tables
# Article: Multi-temporal ALS-based initialization and validation of the iLand
#          forest dynamics model under drought and bark beetle disturbances
# Copyright GNU

# This script reads the cleaned/matched data frames saved at the end of the
# previous processing script and produces:
#   - Static ggplot figures (distributions, agreement, transitions)
#   - Richer static figures: ridgelines (ggridges), correlation matrix (corrplot),
#     pairwise relationships (GGally), species treemap (treemapify)
#   - Interactive figures: plotly (scatter/density/3D) and leaflet (spatial map)
#   - Descriptive statistics tables (CSV + gt-formatted HTML/PNG)

#-------------------------------------------------------------------------------
# 0 -- Setup
#-------------------------------------------------------------------------------
install.packages(c("ggridges", "corrplot", "GGally", "treemapify",
                   "plotly", "leaflet", "htmlwidgets", "gt", "broom"))


library(tidyverse)
library(dplyr)
library(glue)
library(scales)
library(patchwork)
library(viridis)

# Richer static
library(ggridges)
library(corrplot)
library(GGally)
library(treemapify)

# Interactive
library(plotly)
library(leaflet)
install.packages("htmlwidgets")
library(htmlwidgets)

# Tables
install.packages("gt")
library(gt)
library(broom)

rm(list = ls())
base_dir <- "C:/P/DMP_CROSS_CASCADE"

raw_process_data_path <- file.path(base_dir, "03_rawdata/02_process_storage")
figure_path            <- file.path(base_dir, "04_work/03_analysis/02_figure")
table_path             <- file.path(base_dir, "04_work/03_analysis/03_table")
interactive_path        <- file.path(figure_path, "interactive")

dir.create(figure_path,     recursive = TRUE, showWarnings = FALSE)
dir.create(table_path,      recursive = TRUE, showWarnings = FALSE)
dir.create(interactive_path, recursive = TRUE, showWarnings = FALSE)

# Consistent project palette / theme ------------------------------------------
pal_alive_dead <- c("alive" = "#2E7D32", "dead" = "#B23A48")
pal_source     <- c("FID 2015" = "#1B4F72", "FID 2019" = "#2E86C1",
                     "FID 2022" = "#85C1E9", "ALS"      = "#E67E22")

theme_proj <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey35", size = 10),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom"
  )
theme_set(theme_proj)

save_fig <- function(plot, name, width = 8, height = 6) {
  ggsave(file.path(figure_path, paste0(name, ".png")),
         plot = plot, width = width, height = height, dpi = 300, bg = "white")
}

save_html <- function(widget, name) {
  saveWidget(widget, file.path(interactive_path, paste0(name, ".html")),
             selfcontained = TRUE)
}

#-------------------------------------------------------------------------------
# 1 -- Load processed data
#-------------------------------------------------------------------------------

FID_2015_clean        <- readRDS(file.path(raw_process_data_path, "FID_2015_clean.rds"))
FID_2019_clean         <- readRDS(file.path(raw_process_data_path, "FID_2019_clean.rds"))
ALS_clean              <- readRDS(file.path(raw_process_data_path, "ALS_clean.rds"))

FID_2015_clean_alive   <- readRDS(file.path(raw_process_data_path, "FID_2015_clean_alive.rds"))
FID_2015_clean_dead    <- readRDS(file.path(raw_process_data_path, "FID_2015_clean_dead.rds"))
FID_2019_clean_alive   <- readRDS(file.path(raw_process_data_path, "FID_2019_clean_alive.rds"))
FID_2019_clean_dead    <- readRDS(file.path(raw_process_data_path, "FID_2019_clean_dead.rds"))
ALS_clean_alive        <- readRDS(file.path(raw_process_data_path, "ALS_clean_alive.rds"))
ALS_clean_dead         <- readRDS(file.path(raw_process_data_path, "ALS_clean_dead.rds"))

FID_2015_2019_long     <- readRDS(file.path(raw_process_data_path, "FID_2015_2019_long.rds"))
FID_2022                <- readRDS(file.path(raw_process_data_path, "FID_2022.rds"))
FID_2022_D              <- readRDS(file.path(raw_process_data_path, "FID_2022_D.rds"))
plot_count_check        <- readRDS(file.path(raw_process_data_path, "plot_count_check.rds"))
FID_2015_2019_Info_clean <- readRDS(file.path(raw_process_data_path, "FID_2015_2019_Info_clean.rds"))

#-------------------------------------------------------------------------------
# 1b -- Column sanity check (fail fast with a clear message, not a cryptic error)
#-------------------------------------------------------------------------------

required_cols <- list(
  FID_2015_clean = c("plotid", "treeid", "x", "y", "sp_name", "dbh", "h", "dead"),
  ALS_clean      = c("plotid", "treeid", "x", "y", "species", "dbh", "h", "dead")
)

check_cols <- function(df, df_name, cols) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf(
      "Column check failed for %s: missing %s. Available columns: %s",
      df_name, paste(missing, collapse = ", "), paste(names(df), collapse = ", ")
    ))
  }
}
check_cols(FID_2015_clean, "FID_2015_clean", required_cols$FID_2015_clean)
check_cols(ALS_clean,      "ALS_clean",      required_cols$ALS_clean)

struct_compare <- bind_rows(
  FID_2015_clean %>% mutate(source = "FID 2015") %>% select(plotid, dbh, h, dead, source),
  ALS_clean      %>% mutate(source = "ALS")      %>% select(plotid, dbh, h, dead, source)
) %>%
  mutate(status = if_else(dead == 1, "dead", "alive"))

#===============================================================================
# PART A -- Structural agreement FID vs ALS
#===============================================================================

## A1 -- ggridges: DBH distribution by source -----------------------------------
p_dbh_ridges <- struct_compare %>%
  filter(!is.na(dbh)) %>%
  ggplot(aes(x = dbh, y = source, fill = source)) +
  geom_density_ridges(alpha = 0.75, scale = 1.3, color = "white") +
  scale_fill_manual(values = pal_source) +
  labs(title = "DBH distribution by data source",
       subtitle = "Ridgeline comparison, matched plots Bia\u0142owie\u017Ca 2015",
       x = "DBH (cm)", y = NULL, fill = NULL) +
  theme(legend.position = "none")
save_fig(p_dbh_ridges, "A1_dbh_ridges")

## A2 -- ggridges: Height distribution by source --------------------------------
p_h_ridges <- struct_compare %>%
  filter(!is.na(h)) %>%
  ggplot(aes(x = h, y = source, fill = source)) +
  geom_density_ridges(alpha = 0.75, scale = 1.3, color = "white") +
  scale_fill_manual(values = pal_source) +
  labs(title = "Height distribution by data source",
       subtitle = "Ridgeline comparison, matched plots Bia\u0142owie\u017Ca 2015",
       x = "Height (m)", y = NULL, fill = NULL) +
  theme(legend.position = "none")
save_fig(p_h_ridges, "A2_height_ridges")

## A3 -- plotly interactive scatter: stem count FID vs ALS with hover ----------
density_compare <- plot_count_check %>% filter(!is.na(FID2015), !is.na(ALS))
r2_density <- cor(density_compare$FID2015, density_compare$ALS, use = "complete.obs")^2

p_density_scatter <- density_compare %>%
  ggplot(aes(x = FID2015, y = ALS, text = paste0("Plot: ", plotid,
                                                  "<br>FID 2015: ", FID2015,
                                                  "<br>ALS: ", ALS))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.2, alpha = 0.75, color = "#1B4F72") +
  geom_smooth(method = "lm", se = TRUE, color = "#E67E22", linewidth = 0.8) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = paste0("R\u00B2 = ", round(r2_density, 3)), size = 4.2, fontface = "italic") +
  labs(title = "Stem count per plot: field inventory vs ALS detection",
       x = "FID 2015 stem count", y = "ALS-detected stem count")
save_fig(p_density_scatter, "A3_stem_density_scatter")

p_density_scatter_ly <- ggplotly(p_density_scatter, tooltip = "text") %>%
  layout(title = list(text = "Stem count per plot: FID vs ALS (interactive)"))
save_html(p_density_scatter_ly, "A3_stem_density_scatter_interactive")

## A4 -- Bland-Altman agreement plot --------------------------------------------
ba_data <- density_compare %>%
  mutate(mean_count = (FID2015 + ALS) / 2, diff_count = ALS - FID2015)
ba_mean <- mean(ba_data$diff_count)
ba_sd   <- sd(ba_data$diff_count)

p_bland_altman <- ba_data %>%
  ggplot(aes(x = mean_count, y = diff_count, text = paste0("Plot: ", plotid))) +
  geom_hline(yintercept = ba_mean, color = "#1B4F72", linewidth = 0.9) +
  geom_hline(yintercept = ba_mean + 1.96 * ba_sd, linetype = "dashed", color = "#B23A48") +
  geom_hline(yintercept = ba_mean - 1.96 * ba_sd, linetype = "dashed", color = "#B23A48") +
  geom_point(size = 2.2, alpha = 0.75, color = "#2E7D32") +
  labs(title = "Agreement analysis: ALS minus FID stem counts per plot",
       subtitle = "Solid = mean bias, dashed = 95% limits of agreement",
       x = "Mean stem count (FID, ALS)", y = "Difference (ALS \u2212 FID)")
save_fig(p_bland_altman, "A4_bland_altman")

## A5 -- corrplot: correlation matrix of structural metrics ---------------------
corr_input <- struct_compare %>%
  mutate(dead = as.numeric(dead)) %>%
  select(dbh, h, dead) %>%
  filter(complete.cases(.))

corr_matrix <- cor(corr_input, use = "complete.obs")

png(file.path(figure_path, "A5_correlation_matrix.png"), width = 1600, height = 1600, res = 300)
corrplot(corr_matrix, method = "color", type = "upper", addCoef.col = "black",
         tl.col = "black", tl.srt = 45, col = colorRampPalette(c("#B23A48", "white", "#1B4F72"))(200),
         title = "Correlation matrix: DBH, Height, Mortality", mar = c(0, 0, 2, 0))
dev.off()

## A6 -- GGally pairwise plot: DBH, H, status (FID vs ALS) ----------------------
p_pairs_data <- struct_compare %>%
  filter(!is.na(dbh), !is.na(h)) %>%
  select(dbh, h, source, status)
p_pairs_data <- p_pairs_data %>% sample_n(min(2000, nrow(p_pairs_data)))

p_pairs <- p_pairs_data %>%
  ggpairs(columns = 1:2, mapping = aes(color = source, alpha = 0.5),
          upper = list(continuous = wrap("cor", size = 3)),
          lower = list(continuous = wrap("points", size = 0.6)),
          diag  = list(continuous = wrap("densityDiag", alpha = 0.4))) +
  scale_color_manual(values = pal_source) +
  scale_fill_manual(values = pal_source) +
  labs(title = "Pairwise structural relationships: DBH vs Height, FID vs ALS")
ggsave(file.path(figure_path, "A6_pairwise_dbh_height.png"), p_pairs,
       width = 8, height = 8, dpi = 300, bg = "white")

#===============================================================================
# PART B -- Survival / condition dynamics over time
#===============================================================================

status_summary <- bind_rows(
  FID_2015_clean %>% mutate(year = "2015"),
  FID_2019_clean %>% mutate(year = "2019"),
  FID_2022_D     %>% rename(dead = dead2022) %>% mutate(year = "2022")
) %>%
  mutate(status = if_else(dead == 1, "dead", "alive")) %>%
  count(year, status) %>%
  group_by(year) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

## B1 -- Stacked proportion bar --------------------------------------------------
p_status_bar <- status_summary %>%
  ggplot(aes(x = year, y = prop, fill = status)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = percent(prop, accuracy = 1)),
            position = position_stack(vjust = 0.5), color = "white", fontface = "bold") +
  scale_fill_manual(values = pal_alive_dead) +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Alive vs dead tree proportions over time",
       x = NULL, y = "Proportion of trees", fill = NULL)
save_fig(p_status_bar, "B1_alive_dead_proportions")

## B2 -- plotly interactive bar (frame by year) ---------------------------------
p_status_ly <- plot_ly(status_summary, x = ~year, y = ~prop, color = ~status,
                        colors = pal_alive_dead, type = "bar",
                        text = ~percent(prop, accuracy = 0.1), hoverinfo = "text") %>%
  layout(barmode = "stack", title = "Alive vs dead proportions over time (interactive)",
         yaxis = list(title = "Proportion", tickformat = ",.0%"))
save_html(p_status_ly, "B2_alive_dead_proportions_interactive")

## B3 -- Mortality rate violin + jitter ------------------------------------------
mortality_by_plot <- bind_rows(
  FID_2015_clean %>% mutate(year = "2015"),
  FID_2019_clean %>% mutate(year = "2019")
) %>%
  group_by(plotid, year) %>%
  summarise(mortality_rate = mean(dead == 1, na.rm = TRUE), n_trees = n(), .groups = "drop")

p_mortality_violin <- mortality_by_plot %>%
  ggplot(aes(x = year, y = mortality_rate, fill = year)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_jitter(width = 0.08, alpha = 0.4, size = 1.6) +
  scale_fill_manual(values = c("2015" = "#1B4F72", "2019" = "#2E86C1")) +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Distribution of plot-level mortality rates",
       x = NULL, y = "Mortality rate", fill = NULL)
save_fig(p_mortality_violin, "B3_mortality_rate_violin")

## B4 -- Tree-level transitions (2015 -> 2019) -----------------------------------
transition_data <- FID_2015_2019_long %>%
  select(plotid, treeid, year, cond) %>%
  pivot_wider(names_from = year, values_from = cond, names_prefix = "cond_") %>%
  mutate(
    transition = case_when(
      is.na(cond_2015) & !is.na(cond_2019) ~ "new (ingrowth)",
      !is.na(cond_2015) & is.na(cond_2019) ~ "removed/missing",
      TRUE ~ "present both years"
    )
  )

p_transitions <- transition_data %>%
  count(transition) %>%
  mutate(transition = fct_reorder(transition, n)) %>%
  ggplot(aes(x = transition, y = n, fill = transition)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = comma(n)), hjust = -0.15, fontface = "bold") +
  scale_fill_viridis_d(option = "D") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Tree-level presence transitions, 2015 \u2192 2019", x = NULL, y = "Number of trees")
save_fig(p_transitions, "B4_tree_transitions")

## B5 -- plotly 3D: DBH x Height x mortality (FID 2015) --------------------------
p_3d_data <- FID_2015_clean %>%
  filter(!is.na(dbh), !is.na(h)) %>%
  mutate(status = if_else(dead == 1, "dead", "alive"))
p_3d_data <- p_3d_data %>% sample_n(min(3000, nrow(p_3d_data)))

p_3d <- p_3d_data %>%
  plot_ly(x = ~dbh, y = ~h, z = ~as.numeric(factor(status)), color = ~status,
          colors = pal_alive_dead, type = "scatter3d", mode = "markers",
          marker = list(size = 2.5, opacity = 0.6)) %>%
  layout(title = "3D structure: DBH x Height x status (FID 2015)",
         scene = list(xaxis = list(title = "DBH (cm)"),
                       yaxis = list(title = "Height (m)"),
                       zaxis = list(title = "Status (coded)")))
save_html(p_3d, "B5_3d_dbh_height_status")

#===============================================================================
# PART C -- Species composition
#===============================================================================

sp_fid <- FID_2015_clean %>% count(sp_name) %>% mutate(source = "FID 2015")
sp_als <- ALS_clean      %>% count(species) %>% rename(sp_name = species) %>% mutate(source = "ALS")

sp_compare <- bind_rows(sp_fid, sp_als) %>%
  filter(!is.na(sp_name)) %>%
  group_by(source) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

top_species <- sp_compare %>%
  group_by(sp_name) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  slice_max(total, n = 12, with_ties = FALSE) %>%
  pull(sp_name)

## C1 -- Composition bars FID vs ALS ---------------------------------------------
p_species_compare <- sp_compare %>%
  filter(sp_name %in% top_species) %>%
  mutate(sp_name = fct_reorder(sp_name, prop, .fun = sum)) %>%
  ggplot(aes(x = sp_name, y = prop, fill = source)) +
  geom_col(position = position_dodge(0.75), width = 0.65) +
  scale_fill_manual(values = pal_source) +
  scale_y_continuous(labels = percent_format()) +
  coord_flip() +
  labs(title = "Species composition: field inventory vs ALS",
       subtitle = "Top species by combined frequency",
       x = NULL, y = "Proportion of stems", fill = NULL)
save_fig(p_species_compare, "C1_species_composition")

## C2 -- treemapify: species composition as treemap (FID 2015) ------------------
p_treemap <- sp_fid %>%
  filter(!is.na(sp_name)) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(area = n, fill = sp_name, label = paste0(sp_name, "\n", comma(n)))) +
  geom_treemap() +
  geom_treemap_text(color = "white", fontface = "bold", place = "centre", reflow = TRUE, size = 11) +
  scale_fill_viridis_d(option = "turbo") +
  labs(title = "Species composition treemap, FID 2015") +
  theme(legend.position = "none")
save_fig(p_treemap, "C2_species_treemap", width = 9, height = 7)

## C3 -- Species x mortality heatmap ---------------------------------------------
sp_status <- FID_2015_clean %>%
  filter(!is.na(sp_name)) %>%
  mutate(status = if_else(dead == 1, "dead", "alive")) %>%
  count(sp_name, status) %>%
  group_by(sp_name) %>%
  mutate(prop = n / sum(n), total = sum(n)) %>%
  ungroup() %>%
  filter(total >= 10)

p_species_heatmap <- sp_status %>%
  filter(status == "dead") %>%
  mutate(sp_name = fct_reorder(sp_name, prop)) %>%
  ggplot(aes(x = "Mortality", y = sp_name, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = percent(prop, accuracy = 1)), color = "white", fontface = "bold") +
  scale_fill_viridis_c(option = "C", labels = percent_format(), name = "Dead (%)") +
  labs(title = "Species-specific mortality rate", subtitle = "FID 2015, n \u2265 10 stems",
       x = NULL, y = NULL)
save_fig(p_species_heatmap, "C3_species_mortality_heatmap", width = 6, height = 7)

## C4 -- plotly interactive species composition treemap -------------------------
sp_fid_ly <- sp_fid %>% filter(!is.na(sp_name))
p_treemap_ly <- plot_ly(
  type = "treemap",
  labels = sp_fid_ly$sp_name,
  parents = rep("", nrow(sp_fid_ly)),
  values = sp_fid_ly$n,
  textinfo = "label+value+percent parent"
) %>%
  layout(title = "Species composition (interactive treemap), FID 2015")
save_html(p_treemap_ly, "C4_species_treemap_interactive")

#===============================================================================
# PART D -- Spatial distribution
#===============================================================================

top6_species <- ALS_clean %>%
  filter(!is.na(species)) %>%
  count(species) %>%
  slice_max(n, n = 6, with_ties = FALSE) %>%
  pull(species)

## D1 -- Static spatial map by species -------------------------------------------
p_spatial_species <- ALS_clean %>%
  filter(species %in% top6_species) %>%
  ggplot(aes(x = x, y = y, color = species)) +
  geom_point(size = 0.5, alpha = 0.6) +
  coord_equal() +
  scale_color_viridis_d(option = "turbo") +
  labs(title = "Spatial distribution of ALS-detected trees by species",
       subtitle = "Top 6 most frequent species",
       x = "Easting (m)", y = "Northing (m)", color = NULL) +
  theme(legend.position = "right")
save_fig(p_spatial_species, "D1_spatial_map_species", width = 9, height = 7)

## D2 -- plotly interactive spatial scatter (zoomable, hover info) --------------
p_spatial_ly <- ALS_clean %>%
  filter(species %in% top6_species) %>%
  mutate(status = if_else(dead == 1, "dead", "alive")) %>%
  plot_ly(x = ~x, y = ~y, color = ~species, colors = "Set1",
          text = ~paste0("Plot: ", plotid, "<br>Species: ", species,
                          "<br>DBH: ", dbh, "<br>Status: ", status),
          type = "scatter", mode = "markers",
          marker = list(size = 4, opacity = 0.6)) %>%
  layout(title = "ALS tree positions (interactive, zoomable)",
         xaxis = list(title = "Easting (m)", scaleanchor = "y"),
         yaxis = list(title = "Northing (m)"))
save_html(p_spatial_ly, "D2_spatial_map_interactive")

## D3 -- leaflet map: requires WGS84 lat/lon. If x/y are projected (e.g. local --
## grid or UTM), this block will only run correctly once coordinates are
## reprojected. Wrapped in tryCatch so the rest of the script keeps running.
tryCatch({
  library(sf)
  als_sf <- ALS_clean %>%
    filter(!is.na(x), !is.na(y)) %>%
    st_as_sf(coords = c("x", "y"), crs = 2180)  # Polish national grid (PUWG 1992) -- VERIFY this matches your actual CRS
  als_wgs84 <- st_transform(als_sf, crs = 4326)
  coords_ll <- st_coordinates(als_wgs84)

  als_leaflet_data <- als_sf %>%
    st_drop_geometry() %>%
    mutate(lon = coords_ll[, "X"], lat = coords_ll[, "Y"],
           status = if_else(dead == 1, "dead", "alive"))

  leaflet_map <- leaflet(als_leaflet_data) %>%
    addProviderTiles(providers$Esri.WorldImagery) %>%
    addCircleMarkers(lng = ~lon, lat = ~lat, radius = 2, stroke = FALSE,
                      fillOpacity = 0.6,
                      color = ~ifelse(status == "dead", "#B23A48", "#2E7D32"),
                      popup = ~paste0("Plot: ", plotid, "<br>Species: ", species,
                                       "<br>DBH: ", dbh, " cm<br>Status: ", status))
  save_html(leaflet_map, "D3_leaflet_spatial_map")
}, error = function(e) {
  message("Leaflet map skipped -- check/verify the CRS used for x/y in ALS_clean. Error: ", conditionMessage(e))
})

#===============================================================================
# PART E -- Descriptive statistics tables
#===============================================================================

## E1 -- Descriptive stats: DBH and Height by source and status -----------------
desc_stats <- struct_compare %>%
  pivot_longer(cols = c(dbh, h), names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(source, status, variable) %>%
  summarise(
    n      = n(),
    mean   = mean(value),
    sd     = sd(value),
    median = median(value),
    min    = min(value),
    max    = max(value),
    .groups = "drop"
  ) %>%
  arrange(variable, source, status)

write.csv(desc_stats, file.path(table_path, "E1_descriptive_stats_dbh_height.csv"), row.names = FALSE)

gt_desc_stats <- desc_stats %>%
  mutate(variable = recode(variable, dbh = "DBH (cm)", h = "Height (m)")) %>%
  gt(groupname_col = "variable") %>%
  tab_header(title = "Descriptive statistics: DBH and Height",
             subtitle = "By data source and tree status") %>%
  fmt_number(columns = c(mean, sd, median, min, max), decimals = 2) %>%
  cols_label(source = "Source", status = "Status", n = "N") %>%
  data_color(columns = mean, palette = c("#FFFFFF", "#1B4F72"))
gtsave(gt_desc_stats, file.path(table_path, "E1_descriptive_stats_dbh_height.html"))
gtsave(gt_desc_stats, file.path(table_path, "E1_descriptive_stats_dbh_height.png"))

## E2 -- Per-plot stem density summary statistics --------------------------------
density_stats <- density_compare %>%
  summarise(
    n_plots       = n(),
    mean_FID2015  = mean(FID2015), sd_FID2015 = sd(FID2015),
    mean_ALS      = mean(ALS),     sd_ALS     = sd(ALS),
    mean_diff     = mean(ALS - FID2015),
    sd_diff       = sd(ALS - FID2015),
    correlation_r = cor(FID2015, ALS),
    r_squared     = cor(FID2015, ALS)^2
  )

write.csv(density_stats, file.path(table_path, "E2_stem_density_agreement_stats.csv"), row.names = FALSE)

gt_density_stats <- density_stats %>%
  gt() %>%
  tab_header(title = "Stem density agreement: FID vs ALS") %>%
  fmt_number(columns = everything(), decimals = 2)
gtsave(gt_density_stats, file.path(table_path, "E2_stem_density_agreement_stats.html"))

## E3 -- Mortality rate summary by year ------------------------------------------
mortality_stats <- mortality_by_plot %>%
  group_by(year) %>%
  summarise(
    n_plots = n(),
    mean_mortality = mean(mortality_rate),
    sd_mortality   = sd(mortality_rate),
    median_mortality = median(mortality_rate),
    min_mortality  = min(mortality_rate),
    max_mortality  = max(mortality_rate),
    .groups = "drop"
  )

write.csv(mortality_stats, file.path(table_path, "E3_mortality_rate_summary.csv"), row.names = FALSE)

gt_mortality_stats <- mortality_stats %>%
  gt() %>%
  tab_header(title = "Plot-level mortality rate summary by year") %>%
  fmt_percent(columns = c(mean_mortality, sd_mortality, median_mortality,
                           min_mortality, max_mortality), decimals = 1) %>%
  cols_label(year = "Year", n_plots = "N plots")
gtsave(gt_mortality_stats, file.path(table_path, "E3_mortality_rate_summary.html"))
gtsave(gt_mortality_stats, file.path(table_path, "E3_mortality_rate_summary.png"))

## E4 -- Species frequency table FID vs ALS, with chi-squared test --------------
sp_wide <- sp_compare %>%
  select(sp_name, source, n) %>%
  pivot_wider(names_from = source, values_from = n, values_fill = 0) %>%
  filter(sp_name %in% top_species)

sp_contingency <- sp_wide %>% select(-sp_name) %>% as.matrix()
rownames(sp_contingency) <- sp_wide$sp_name

chisq_species <- suppressWarnings(chisq.test(sp_contingency))
chisq_summary <- tidy(chisq_species)

write.csv(sp_wide, file.path(table_path, "E4_species_frequency_FID_vs_ALS.csv"), row.names = FALSE)
write.csv(chisq_summary, file.path(table_path, "E4_species_chisq_test.csv"), row.names = FALSE)

gt_species_freq <- sp_wide %>%
  gt() %>%
  tab_header(title = "Species frequency: FID 2015 vs ALS",
             subtitle = paste0("Chi-squared test: \u03C7\u00B2 = ", round(chisq_summary$statistic, 2),
                                ", df = ", chisq_summary$parameter,
                                ", p ", ifelse(chisq_summary$p.value < 0.001, "< 0.001",
                                               paste0("= ", round(chisq_summary$p.value, 3))))) %>%
  cols_label(sp_name = "Species")
gtsave(gt_species_freq, file.path(table_path, "E4_species_frequency_FID_vs_ALS.html"))
gtsave(gt_species_freq, file.path(table_path, "E4_species_frequency_FID_vs_ALS.png"))

## E5 -- Plot-level richness statistics -------------------------------------------
richness_by_plot <- FID_2015_clean %>%
  filter(!is.na(sp_name)) %>%
  group_by(plotid) %>%
  summarise(richness = n_distinct(sp_name), n_trees = n(), .groups = "drop")

richness_stats <- richness_by_plot %>%
  summarise(
    n_plots = n(),
    mean_richness = mean(richness), sd_richness = sd(richness),
    mean_trees = mean(n_trees), sd_trees = sd(n_trees),
    correlation_richness_density = cor(richness, n_trees)
  )

write.csv(richness_stats, file.path(table_path, "E5_richness_summary.csv"), row.names = FALSE)

gt_richness <- richness_stats %>%
  gt() %>%
  tab_header(title = "Species richness and stand density summary") %>%
  fmt_number(columns = everything(), decimals = 2)
gtsave(gt_richness, file.path(table_path, "E5_richness_summary.html"))

## E5b -- Richness scatter plot ---------------------------------------------------
p_richness <- richness_by_plot %>%
  ggplot(aes(x = n_trees, y = richness)) +
  geom_point(alpha = 0.6, size = 2.2, color = "#2E7D32") +
  geom_smooth(method = "loess", se = TRUE, color = "#1B4F72", linewidth = 0.8) +
  labs(title = "Species richness vs stand density per plot", subtitle = "FID 2015",
       x = "Number of trees per plot", y = "Number of species per plot")
save_fig(p_richness, "C5_species_richness_vs_density")

## E6 -- Welch two-sample t-test: DBH, FID vs ALS (alive trees only) -------------
ttest_dbh <- t.test(dbh ~ source, data = struct_compare %>% filter(status == "alive", !is.na(dbh)))
ttest_h   <- t.test(h ~ source,   data = struct_compare %>% filter(status == "alive", !is.na(h)))

ttest_summary <- bind_rows(
  tidy(ttest_dbh) %>% mutate(variable = "DBH (cm)"),
  tidy(ttest_h)   %>% mutate(variable = "Height (m)")
) %>%
  select(variable, estimate, estimate1, estimate2, statistic, p.value, conf.low, conf.high)

write.csv(ttest_summary, file.path(table_path, "E6_ttest_FID_vs_ALS.csv"), row.names = FALSE)

gt_ttest <- ttest_summary %>%
  gt() %>%
  tab_header(title = "Welch two-sample t-test: FID vs ALS (alive trees)",
             subtitle = "estimate1 = FID mean, estimate2 = ALS mean") %>%
  fmt_number(columns = c(estimate, estimate1, estimate2, statistic, conf.low, conf.high), decimals = 2) %>%
  fmt_number(columns = p.value, decimals = 4)
gtsave(gt_ttest, file.path(table_path, "E6_ttest_FID_vs_ALS.html"))
gtsave(gt_ttest, file.path(table_path, "E6_ttest_FID_vs_ALS.png"))

#===============================================================================
# PART F -- Composite summary dashboard (patchwork, static figures only)
#===============================================================================

dashboard <- (p_dbh_ridges + p_h_ridges) /
             (p_status_bar + p_species_compare) /
             (p_spatial_species + p_richness) +
  plot_annotation(
    title = "Bia\u0142owie\u017Ca FID-ALS Initialization: Diagnostic Dashboard",
    subtitle = "Structural agreement, temporal dynamics, species composition, and spatial structure",
    theme = theme(plot.title = element_text(face = "bold", size = 18),
                  plot.subtitle = element_text(size = 12, color = "grey30"))
  )
save_fig(dashboard, "F1_diagnostic_dashboard", width = 14, height = 16)

#-------------------------------------------------------------------------------
# Done
#-------------------------------------------------------------------------------
message("Static figures saved to:    ", figure_path)
message("Interactive HTML saved to:  ", interactive_path)
message("Tables (CSV/gt) saved to:   ", table_path)
