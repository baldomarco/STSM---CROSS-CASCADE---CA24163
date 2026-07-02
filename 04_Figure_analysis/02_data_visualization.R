# Marco Baldo 21-06-2026 contact: baldo@fld.czu.cz
# Data visualization & diagnostic plots
# Article: Multi-temporal ALS-based initialization and validation of the iLand
#          forest dynamics model under drought and bark beetle disturbances
# Copyright GNU

# This script reads the cleaned/matched data frames saved at the end of the
# previous processing script and produces a comprehensive set of diagnostic
# and exploratory figures: structural agreement FID vs ALS, survival/condition
# over time (2015-2019-2022), species composition, and spatial distribution.

#-------------------------------------------------------------------------------
# 0 -- Setup
#-------------------------------------------------------------------------------

library(tidyverse)
library(scales)
library(patchwork)
library(ggrepel)
library(viridis)

rm(list = ls())
base_dir <- "C:/P/DMP_CROSS_CASCADE"

raw_process_data_path <- file.path(base_dir, "03_rawdata/02_process_storage")
figure_path            <- file.path(base_dir, "04_work/03_analysis/02_figure")

dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)

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

#-------------------------------------------------------------------------------
# 1 -- Load processed data
#-------------------------------------------------------------------------------

FID_2015_clean       <- readRDS(file.path(raw_process_data_path, "FID_2015_clean.rds"))
FID_2019_clean        <- readRDS(file.path(raw_process_data_path, "FID_2019_clean.rds"))
ALS_clean             <- readRDS(file.path(raw_process_data_path, "ALS_clean.rds"))

FID_2015_clean_alive  <- readRDS(file.path(raw_process_data_path, "FID_2015_clean_alive.rds"))
FID_2015_clean_dead   <- readRDS(file.path(raw_process_data_path, "FID_2015_clean_dead.rds"))
FID_2019_clean_alive  <- readRDS(file.path(raw_process_data_path, "FID_2019_clean_alive.rds"))
FID_2019_clean_dead   <- readRDS(file.path(raw_process_data_path, "FID_2019_clean_dead.rds"))
ALS_clean_alive       <- readRDS(file.path(raw_process_data_path, "ALS_clean_alive.rds"))
ALS_clean_dead        <- readRDS(file.path(raw_process_data_path, "ALS_clean_dead.rds"))

FID_2015_2019_long    <- readRDS(file.path(raw_process_data_path, "FID_2015_2019_long.rds"))
FID_2022               <- readRDS(file.path(raw_process_data_path, "FID_2022.rds"))
FID_2022_D             <- readRDS(file.path(raw_process_data_path, "FID_2022_D.rds"))
plot_count_check       <- readRDS(file.path(raw_process_data_path, "plot_count_check.rds"))
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

#-------------------------------------------------------------------------------
# 2 -- PART A: Structural agreement FID vs ALS (DBH, height, density)
#-------------------------------------------------------------------------------

## 2.1 -- Combine FID 2015 and ALS into one tidy comparison frame --------------
struct_compare <- bind_rows(
  FID_2015_clean %>% mutate(source = "FID 2015") %>% select(plotid, dbh, h, dead, source),
  ALS_clean      %>% mutate(source = "ALS")      %>% select(plotid, dbh, h, dead, source)
) %>%
  mutate(status = if_else(dead == 1, "dead", "alive"))

## 2.2 -- DBH distribution: FID vs ALS (density ridges via overlapping density) -
p_dbh_dist <- struct_compare %>%
  filter(!is.na(dbh)) %>%
  ggplot(aes(x = dbh, fill = source, color = source)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_fill_manual(values = pal_source) +
  scale_color_manual(values = pal_source) +
  labs(title = "DBH distribution: field inventory vs ALS",
       subtitle = "Matched plots, Bia\u0142owie\u017Ca 2015",
       x = "DBH (cm)", y = "Density", fill = NULL, color = NULL)
save_fig(p_dbh_dist, "A1_dbh_distribution_FID_vs_ALS")

## 2.3 -- Height distribution: FID vs ALS ---------------------------------------
p_h_dist <- struct_compare %>%
  filter(!is.na(h)) %>%
  ggplot(aes(x = h, fill = source, color = source)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_fill_manual(values = pal_source) +
  scale_color_manual(values = pal_source) +
  labs(title = "Tree height distribution: field inventory vs ALS",
       subtitle = "Matched plots, Bia\u0142owie\u017Ca 2015",
       x = "Height (m)", y = "Density", fill = NULL, color = NULL)
save_fig(p_h_dist, "A2_height_distribution_FID_vs_ALS")

## 2.4 -- Per-plot stem density comparison (scatter + 1:1 line) ----------------
density_compare <- plot_count_check %>%
  filter(!is.na(FID2015), !is.na(ALS))

r2_density <- cor(density_compare$FID2015, density_compare$ALS, use = "complete.obs")^2

p_density_scatter <- density_compare %>%
  ggplot(aes(x = FID2015, y = ALS)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5, alpha = 0.7, color = "#1B4F72") +
  geom_smooth(method = "lm", se = TRUE, color = "#E67E22", linewidth = 0.8) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = paste0("R\u00B2 = ", round(r2_density, 3)), size = 4.2, fontface = "italic") +
  labs(title = "Stem count per plot: field inventory vs ALS detection",
       subtitle = "Dashed line = perfect 1:1 agreement",
       x = "FID 2015 stem count", y = "ALS-detected stem count")
save_fig(p_density_scatter, "A3_stem_density_FID_vs_ALS_scatter")

## 2.5 -- Bland-Altman style agreement plot (mean vs difference, per plot) -----
ba_data <- density_compare %>%
  mutate(mean_count = (FID2015 + ALS) / 2,
         diff_count = ALS - FID2015)

ba_mean <- mean(ba_data$diff_count)
ba_sd   <- sd(ba_data$diff_count)

p_bland_altman <- ba_data %>%
  ggplot(aes(x = mean_count, y = diff_count)) +
  geom_hline(yintercept = ba_mean, color = "#1B4F72", linewidth = 0.9) +
  geom_hline(yintercept = ba_mean + 1.96 * ba_sd, linetype = "dashed", color = "#B23A48") +
  geom_hline(yintercept = ba_mean - 1.96 * ba_sd, linetype = "dashed", color = "#B23A48") +
  geom_point(size = 2.5, alpha = 0.7, color = "#2E7D32") +
  labs(title = "Agreement analysis: ALS minus FID stem counts per plot",
       subtitle = "Solid = mean bias, dashed = 95% limits of agreement",
       x = "Mean stem count (FID, ALS)", y = "Difference (ALS \u2212 FID)")
save_fig(p_bland_altman, "A4_bland_altman_stem_counts")

## 2.6 -- DBH agreement boxplot by alive/dead status ---------------------------
p_dbh_box <- struct_compare %>%
  filter(!is.na(dbh)) %>%
  ggplot(aes(x = source, y = dbh, fill = status)) +
  geom_boxplot(outlier.alpha = 0.3, position = position_dodge(0.8), width = 0.65) +
  scale_fill_manual(values = pal_alive_dead) +
  labs(title = "DBH by data source and tree status",
       x = NULL, y = "DBH (cm)", fill = NULL)
save_fig(p_dbh_box, "A5_dbh_boxplot_source_status")

#-------------------------------------------------------------------------------
# 3 -- PART B: Survival / condition dynamics over time
#-------------------------------------------------------------------------------

## 3.1 -- Overall alive/dead proportions per year (2015, 2019, 2022) -----------
status_summary <- bind_rows(
  FID_2015_clean %>% mutate(year = "2015"),
  FID_2019_clean %>% mutate(year = "2019"),
  FID_2022_D     %>% rename(dead = dead2022) %>% mutate(year = "2022")
) %>%
  mutate(status = if_else(dead == 1, "dead", "alive")) %>%
  count(year, status) %>%
  group_by(year) %>%
  mutate(prop = n / sum(n))

p_status_bar <- status_summary %>%
  ggplot(aes(x = year, y = prop, fill = status)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = percent(prop, accuracy = 1)),
            position = position_stack(vjust = 0.5), color = "white", fontface = "bold") +
  scale_fill_manual(values = pal_alive_dead) +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Alive vs dead tree proportions over time",
       subtitle = "Field Inventory Data, Bia\u0142owie\u017Ca",
       x = NULL, y = "Proportion of trees", fill = NULL)
save_fig(p_status_bar, "B1_alive_dead_proportions_over_time")

## 3.2 -- Per-plot mortality rate distribution (2015 vs 2019) ------------------
mortality_by_plot <- bind_rows(
  FID_2015_clean %>% mutate(year = "2015"),
  FID_2019_clean %>% mutate(year = "2019")
) %>%
  group_by(plotid, year) %>%
  summarise(mortality_rate = mean(dead == 1, na.rm = TRUE), .groups = "drop")

p_mortality_violin <- mortality_by_plot %>%
  ggplot(aes(x = year, y = mortality_rate, fill = year)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_jitter(width = 0.08, alpha = 0.4, size = 1.6) +
  scale_fill_manual(values = c("2015" = "#1B4F72", "2019" = "#2E86C1")) +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Distribution of plot-level mortality rates",
       x = NULL, y = "Mortality rate (% dead trees per plot)", fill = NULL)
save_fig(p_mortality_violin, "B2_mortality_rate_distribution")

## 3.3 -- Tree-level transitions: alive->alive, alive->dead, missing, etc. -----
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
  labs(title = "Tree-level presence transitions, 2015 \u2192 2019",
       x = NULL, y = "Number of trees")
save_fig(p_transitions, "B3_tree_transitions_2015_2019")

## 3.4 -- Plot-level tree count trajectory across all three years -------------
plot_counts_years <- bind_rows(
  FID_2015_clean %>% count(plotid) %>% mutate(year = 2015),
  FID_2019_clean %>% count(plotid) %>% mutate(year = 2019),
  FID_2022       %>% count(plotid) %>% mutate(year = 2022)
)

p_trajectory <- plot_counts_years %>%
  ggplot(aes(x = year, y = n, group = plotid)) +
  geom_line(alpha = 0.12, color = "#1B4F72") +
  stat_summary(fun = mean, geom = "line", color = "#E67E22", linewidth = 1.4, aes(group = 1)) +
  stat_summary(fun = mean, geom = "point", color = "#E67E22", size = 3, aes(group = 1)) +
  scale_x_continuous(breaks = c(2015, 2019, 2022)) +
  labs(title = "Tree count trajectory per plot, 2015\u20132022",
       subtitle = "Thin lines = individual plots, orange = mean trend",
       x = NULL, y = "Trees per plot")
save_fig(p_trajectory, "B4_plot_tree_count_trajectory")

#-------------------------------------------------------------------------------
# 4 -- PART C: Species composition
#-------------------------------------------------------------------------------

## 4.1 -- Species composition: FID 2015 vs ALS (side-by-side bars) ------------
sp_fid <- FID_2015_clean %>% count(sp_name) %>% mutate(source = "FID 2015")
sp_als <- ALS_clean      %>% count(species) %>% rename(sp_name = species) %>% mutate(source = "ALS")

sp_compare <- bind_rows(sp_fid, sp_als) %>%
  filter(!is.na(sp_name)) %>%
  group_by(source) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

top_species <- sp_compare %>%
  group_by(sp_name) %>%
  summarise(total = sum(n)) %>%
  slice_max(total, n = 12) %>%
  pull(sp_name)

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
save_fig(p_species_compare, "C1_species_composition_FID_vs_ALS")

## 4.2 -- Species x mortality status heatmap (FID 2015) ------------------------
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
  labs(title = "Species-specific mortality rate",
       subtitle = "FID 2015, species with n \u2265 10 stems",
       x = NULL, y = NULL)
save_fig(p_species_heatmap, "C2_species_mortality_heatmap", width = 6, height = 7)

## 4.3 -- Species richness per plot (diversity proxy) --------------------------
richness_by_plot <- FID_2015_clean %>%
  filter(!is.na(sp_name)) %>%
  group_by(plotid) %>%
  summarise(richness = n_distinct(sp_name), n_trees = n())

p_richness <- richness_by_plot %>%
  ggplot(aes(x = n_trees, y = richness)) +
  geom_point(alpha = 0.6, size = 2.2, color = "#2E7D32") +
  geom_smooth(method = "loess", se = TRUE, color = "#1B4F72", linewidth = 0.8) +
  labs(title = "Species richness vs stand density per plot",
       subtitle = "FID 2015",
       x = "Number of trees per plot", y = "Number of species per plot")
save_fig(p_richness, "C3_species_richness_vs_density")

#-------------------------------------------------------------------------------
# 5 -- PART D: Spatial distribution
#-------------------------------------------------------------------------------

## 5.1 -- Spatial map of ALS trees colored by species (top species only) ------
top6_species <- ALS_clean %>%
  filter(!is.na(species)) %>%
  count(species) %>%
  slice_max(n, n = 6) %>%
  pull(species)

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
save_fig(p_spatial_species, "D1_spatial_map_species_ALS", width = 9, height = 7)

## 5.2 -- Spatial map of ALS trees colored by alive/dead -----------------------
p_spatial_status <- ALS_clean %>%
  mutate(status = if_else(dead == 1, "dead", "alive")) %>%
  ggplot(aes(x = x, y = y, color = status)) +
  geom_point(size = 0.5, alpha = 0.6) +
  coord_equal() +
  scale_color_manual(values = pal_alive_dead) +
  labs(title = "Spatial distribution of ALS-detected trees by status",
       x = "Easting (m)", y = "Northing (m)", color = NULL) +
  theme(legend.position = "right")
save_fig(p_spatial_status, "D2_spatial_map_status_ALS", width = 9, height = 7)

## 5.3 -- Single example plot footprint: FID vs ALS tree positions ------------
example_plot <- plot_count_check %>%
  filter(!is.na(FID2015), FID2015 > 0) %>%
  slice_max(FID2015, n = 1) %>%
  pull(plotid) %>%
  first()

example_compare <- bind_rows(
  FID_2015_clean %>% filter(plotid == example_plot) %>%
    mutate(source = "FID 2015") %>% select(x, y, source),
  ALS_clean %>% filter(plotid == example_plot) %>%
    mutate(source = "ALS") %>% select(x, y, source)
)

p_example_plot <- example_compare %>%
  ggplot(aes(x = x, y = y, color = source, shape = source)) +
  geom_point(size = 3, alpha = 0.75) +
  coord_equal() +
  scale_color_manual(values = pal_source) +
  labs(title = paste0("Tree positions: FID vs ALS \u2014 plot ", example_plot),
       subtitle = "Example single-plot footprint comparison",
       x = "X (local)", y = "Y (local)", color = NULL, shape = NULL)
save_fig(p_example_plot, "D3_example_plot_footprint_FID_vs_ALS")

#-------------------------------------------------------------------------------
# 6 -- PART E: Composite summary dashboard (patchwork)
#-------------------------------------------------------------------------------

dashboard <- (p_dbh_dist + p_h_dist) /
             (p_status_bar + p_species_compare) /
             (p_spatial_status + p_richness) +
  plot_annotation(
    title = "Bia\u0142owie\u017Ca FID-ALS Initialization: Diagnostic Dashboard",
    subtitle = "Structural agreement, temporal dynamics, species composition, and spatial structure",
    theme = theme(plot.title = element_text(face = "bold", size = 18),
                  plot.subtitle = element_text(size = 12, color = "grey30"))
  )
save_fig(dashboard, "E1_diagnostic_dashboard", width = 14, height = 16)

#-------------------------------------------------------------------------------
# Done -- all figures written to figure_path
#-------------------------------------------------------------------------------
message("All figures saved to: ", figure_path)

