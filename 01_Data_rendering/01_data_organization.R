# Marco Baldo 17-06-2026 contact: baldo@fld.czu.cz
# Initial data organization from metadata
# Article: Multi-temporal ALS-based initialization and validation of the iLand forest dynamics model under drought and bark beetle disturbances
# To cite these functions, use:

# Copyright GNU 


# Upload the required libraries
library(readxl)
library(tidyverse)

rm(list = ls())
base_dir <- "C:/P/DMP_CROSS_CASCADE"

raw_metadata_path <- file.path(base_dir, "03_rawdata/01_metadata")
raw_process_data_path <- file.path(base_dir, "03_rawdata/02_process_storage")

figure_path <- file.path(base_dir, "04_work/03_analysis/02_figure")
table_path  <- file.path(base_dir, "04_work/03_analysis/03_table")



#1-------------------------------------------------------------------------------
# Load and prepare the Field Inventory Data 

# FID 2015 & 2019
FID_2015_2019 <- read_excel(
  file.path(raw_metadata_path, "bialowieza_field_2015-2019.xlsx"),
  sheet = "bialowieza_field_2015-2019")

head(FID_2015_2019)
unique(FID_2015_2019$plotid)

# FID 2022
FID_2022 <- read_excel(
  file.path(raw_metadata_path, "bialowieza_field_2022.xlsx"),
  sheet = "Sheet1")

head(FID_2022)
unique(FID_2022$plotid)

#-------------------------------------------------------------------------------
# Convert the species names in scientific names

species_ref <- tribble(
  ~sp_ibl, ~sp_name, ~species,
  "BRZ", "Betula pendula", "bepe",
  "BRZB", "Betula pendula", "bepe",
  "BRZO", "Betula pubescens", "bepu",
  "CZZ", "Prunus serotina", "prse",
  "DB", "Quercus species", "qupe",
  "DBB", "Quercus petraea", "qupe",
  "DBS", "Quercus robur", "quro",
  "GB", "Carpinus betulus", "cabe",
  "GR", "Pyrus communis", "pyco",
  "JB", "Malus silvestris", "masi",
  "JRZ", "Sorbus aucuparia", "soau",
  "JS", "Fraxinus excelsior", "frex",
  "JW", "Acer pseudoplatanus", "acps",
  "KLZ", "Acer platanoides", "acpl",
  "LESZCZ", "Corylus avellana", "coav",
  "LP", "Tilia cordata", "tico",
  "LPD", "Tilia petiolaris", "tipe",
  "OL", "Alnus glutinosa", "algl",
  "SL", "Prunus spinosa", "prse",
  "OS", "Populus tremula", "potr",
  "SO", "Pinus sylvestris", "pisy",
  "SW", "Picea abies", "piab",
  "WB", "Salix alba", "saal",
  "WBI", "Salix alba", "saal",
  "WBK", "Salix fragilis", "safr",
  "WBW", "Salix viminalis", "savi",
  "WIP", "Salix caprea", "saca",
  "WZ", "Ulmus", "ulgl",
  "WZG", "Ulmus", "ulgl",
  "WZP", "Ulmus carpinifolia", "ulca",
  "WZS", "Ulmus laevis", "ulla"
)

# Rename and order the column for species names and iland codes 2015 - 2019
FID_2015_2019 <- FID_2015_2019 %>%
  rename(sp_ibl = species) %>%
  left_join(species_ref, by = "sp_ibl") %>%
  relocate(sp_ibl, sp_name, species, .after = sp_ibl)

# Check the species codes for iLand
species <- unique(FID_2015_2019$sp_name)
species_iland <- unique(FID_2015_2019$species)

# Rename and order the column for species names and iland codes 2022
FID_2022 <- FID_2022 %>%
  rename(sp_ibl = species) %>%
  left_join(species_ref, by = "sp_ibl") %>%
  relocate(sp_ibl, sp_name, species, .after = sp_ibl)

# Check the species codes for iLand
species <- unique(FID_2022$sp_name)
species_iland <- unique(FID_2022$species)

#-------------------------------------------------------------------------------
# split the dead and alive trees in 0 - 1 dead/alive 2015_2019
FID_2015_2019_D <- FID_2015_2019 %>%
  mutate(dead2015 = ifelse(cond2015 %in% c("SC", "Z", "ZZ"), 0, 1),
         dead2019 = ifelse(cond2019 %in% c("SC", "Z", "ZZ"), 0, 1))

df_alive_2015 <- FID_2015_2019_D %>% filter(dead2015 == 0)
df_dead_2015  <- FID_2015_2019_D %>% filter(dead2015 == 1)

# Prepare the data for the splitting logically dividing the variables matching assigned by 2015 or 2019 inserting a year column
long_df <- FID_2015_2019_D %>%
  pivot_longer(
    cols = matches("(2015|2019)$"),
    names_to = c(".value", "year"),
    names_pattern = "(.*)(2015|2019)"
  )

#-------------------------------------------------------------------------------
# split the dead and alive trees in 0 - 1 dead/alive 2022
FID_2022_D <- FID_2022 %>%
  mutate(dead2022 = ifelse(cond %in% c("SC", "Z", "ZZ"), 0, 1))

df_alive_2022 <- FID_2022_D %>% filter(dead2022 == 0)
df_dead_2022  <- FID_2022_D %>% filter(dead2022 == 1)

#-------------------------------------------------------------------------------
# Create two separate df 2015 - 2019 & alive - dead
FID_2015 <- long_df %>% filter(year == "2015") %>% select(-year)
FID_2019 <- long_df %>% filter(year == "2019") %>% select(-year)

FID_2015_alive <- FID_2015 %>% filter(cond %in% c("SC","Z", "ZZ"))
FID_2015_dead  <- FID_2015 %>% filter(!(cond %in% c("SC","Z", "ZZ")))

FID_2019_alive <- FID_2019 %>% filter(cond %in% c("SC","Z", "ZZ"))
FID_2019_dead  <- FID_2019 %>% filter(!(cond %in% c("SC","Z", "ZZ")))

#-------------------------------------------------------------------------------
# 2022 alive / dead (same logic as above)
FID_2022_alive <- FID_2022_D %>% filter(cond %in% c("SC","Z", "ZZ"))
FID_2022_dead  <- FID_2022_D %>% filter(!(cond %in% c("SC","Z", "ZZ")))

#-------------------------------------------------------------------------------
# Remove NA rows in split dfs
FID_2015 <- FID_2015 %>% filter(!is.na(cond))
FID_2019 <- FID_2019 %>% filter(!is.na(cond))

#-------------------------------------------------------------------------------
# Test if is true that always a tree missing is present the year after or before
FID_2015_2019 %>%
  filter(is.na(cond2015) & is.na(cond2019)) %>%
  nrow()

#-------------------------------------------------------------------------------
# Check the whole n. of presence in 2015-2019 and the n. removal - regrowth
check_missing <- FID_2015_2019 %>%
  mutate(
    missing_2015 = is.na(cond2015),
    missing_2019 = is.na(cond2019)
  ) %>%
  count(missing_2015, missing_2019)

#-------------------------------------------------------------------------------
# Inspect cases of removal - regrowth
problem_cases <- FID_2015_2019 %>%
  filter(is.na(cond2015) | is.na(cond2019))

#-------------------------------------------------------------------------------
# time-series modeling - mixed models - growth/mortality transitions analysis use it
FID_2015_2019_long <- long_df

#-------------------------------------------------------------------------------
# FID 2015 vs 2019 differences
setdiff(unique(FID_2015$plotid), unique(FID_2019$plotid))
setdiff(unique(FID_2019$plotid), unique(FID_2015$plotid))

#-------------------------------------------------------------------------------
# ADDITIONAL: include 2022 in comparison (separate from 2015–2019 logic)

# FID 2015 vs 2022 differences
setdiff(unique(FID_2015$plotid), unique(FID_2022$plotid))
setdiff(unique(FID_2022$plotid), unique(FID_2015$plotid))

# FID 2019 vs 2022 differences
setdiff(unique(FID_2019$plotid), unique(FID_2022$plotid))
setdiff(unique(FID_2022$plotid), unique(FID_2019$plotid))



#2------------------------------------------------------------------------------
# Save the first part of the results - analysis on Forest Inventory Data

# SAVE IN CSV

# Species and sp code added 
write.csv(FID_2015_2019, file.path(raw_process_data_path, "FID_2015_2019.csv"), row.names = FALSE)
write.csv(FID_2022,      file.path(raw_process_data_path, "FID_2022.csv"), row.names = FALSE)

# Whole dataset included of dead alive 0-1 column
write.csv(FID_2015_2019_D, file.path(raw_process_data_path, "FID_2015_2019_D.csv"), row.names = FALSE)
write.csv(FID_2022_D,      file.path(raw_process_data_path, "FID_2022_D.csv"), row.names = FALSE)

# Separation per years and dead or alive for the 2015-2019
write.csv(FID_2015,        file.path(raw_process_data_path, "FID_2015.csv"), row.names = FALSE)
write.csv(FID_2019,        file.path(raw_process_data_path, "FID_2019.csv"), row.names = FALSE)

write.csv(FID_2015_alive,  file.path(raw_process_data_path, "FID_2015_alive.csv"), row.names = FALSE)
write.csv(FID_2015_dead,   file.path(raw_process_data_path, "FID_2015_dead.csv"), row.names = FALSE)

write.csv(FID_2019_alive,  file.path(raw_process_data_path, "FID_2019_alive.csv"), row.names = FALSE)
write.csv(FID_2019_dead,   file.path(raw_process_data_path, "FID_2019_dead.csv"), row.names = FALSE)

# The same but for 2022 
write.csv(FID_2019_alive,  file.path(raw_process_data_path, "FID_2019_alive.csv"), row.names = FALSE)
write.csv(FID_2019_dead,   file.path(raw_process_data_path, "FID_2019_dead.csv"), row.names = FALSE)

# More cohemprensive dataset for complex analysis
write.csv(FID_2015_2019_long, file.path(raw_process_data_path, "FID_2015_2019_long.csv"), row.names = FALSE)
write.csv(problem_cases,   file.path(raw_process_data_path, "problem_cases.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# SAVE IN RDS
saveRDS(FID_2015_2019, file.path(raw_process_data_path, "FID_2015_2019.rds"))
saveRDS(FID_2022,      file.path(raw_process_data_path, "FID_2022.rds"))

saveRDS(FID_2015_2019_D, file.path(raw_process_data_path, "FID_2015_2019_D.rds"))
saveRDS(FID_2022_D,      file.path(raw_process_data_path, "FID_2022_D.rds"))

saveRDS(FID_2015,        file.path(raw_process_data_path, "FID_2015.rds"))
saveRDS(FID_2019,        file.path(raw_process_data_path, "FID_2019.rds"))

saveRDS(FID_2015_alive,  file.path(raw_process_data_path, "FID_2015_alive.rds"))
saveRDS(FID_2015_dead,   file.path(raw_process_data_path, "FID_2015_dead.rds"))

saveRDS(FID_2019_alive,  file.path(raw_process_data_path, "FID_2019_alive.rds"))
saveRDS(FID_2019_dead,   file.path(raw_process_data_path, "FID_2019_dead.rds"))

saveRDS(FID_2022_alive,  file.path(raw_process_data_path, "FID_2022_alive.rds"))
saveRDS(FID_2022_dead,   file.path(raw_process_data_path, "FID_2022_dead.rds"))

saveRDS(FID_2015_2019_long, file.path(raw_process_data_path, "FID_2015_2019_long.rds"))
saveRDS(problem_cases,      file.path(raw_process_data_path, "problem_cases.rds"))


#3------------------------------------------------------------------------------
# Let's process the GIS data ALS

library(sf)

# Load SHP and extract attribute table
shp <- st_read("C:/P/DMP_CROSS_CASCADE/03_rawdata/01_metadata/rs_trees_2025/itd_pb_2025.shp")

attr_table <- st_drop_geometry(shp)

ALS_2015 <- attr_table

head(ALS_2015)
unique(ALS_2015$plot_id)
species <- unique(ALS_2015$species)

### TEST THE DATA SET AND SPECIES ####
# --- Check 1: is "other" all alive, or a mix of alive/dead? ---
ALS_2015 %>%
  filter(species == "other") %>%
  count(dead)

# Also worth seeing the raw values in case 'dead' isn't a clean 0/1
ALS_2015 %>%
  filter(species == "other") %>%
  pull(dead) %>%
  table()

# --- Check 2: are "decidious" and "snag" all dead, or mixed? ---
ALS_2015 %>%
  filter(species %in% c("decidious", "snag")) %>%
  count(species, dead)



###    CORRECTION ALS TABLE  ######
# Correct the errors in snag, and assigne the species as for FID where possible
ALS_2015 <- ALS_2015 %>%
  mutate(
    dead = case_when(
      species == "snag" ~ 1,
      TRUE ~ dead
    ),
    species = case_when(
      species == "other" ~ "deciduous",
      species == "decidious" ~ "snag_broadl",
      TRUE ~ species
    )
  )

ALS_2015 %>% count(species, dead)

#  LEt's assigne the proper names
ALS_2015 <- ALS_2015 %>%
  rename(sp_ibl = species) %>%
  mutate(
    sp_name = case_when(
      sp_ibl == "hornbeam" ~ "Carpinus betulus",
      sp_ibl == "maple"    ~ "Acer platanoides",
      sp_ibl == "spruce"   ~ "Picea abies",
      sp_ibl == "oak"      ~ "Quercus robur",
      sp_ibl == "lime"     ~ "Tilia cordata",
      sp_ibl == "pine"     ~ "Pinus sylvestris",
      sp_ibl == "ash"      ~ "Fraxinus excelsior",
      sp_ibl == "aspen"    ~ "Populus tremula",
      sp_ibl == "alder"    ~ "Alnus glutinosa",
      sp_ibl == "birch"    ~ "Betula pendula",
      TRUE ~ sp_ibl
    ),
    species = case_when(
      sp_ibl == "hornbeam" ~ "cabe",
      sp_ibl == "maple"    ~ "acpl",
      sp_ibl == "spruce"   ~ "piab",
      sp_ibl == "oak"      ~ "quro",
      sp_ibl == "lime"     ~ "tico",
      sp_ibl == "pine"     ~ "pisy",
      sp_ibl == "ash"      ~ "frex",
      sp_ibl == "aspen"    ~ "potr",
      sp_ibl == "alder"    ~ "algl",
      sp_ibl == "birch"    ~ "bepe",
      TRUE ~ sp_ibl
    )
  ) %>%
  relocate(sp_name, species, .after = sp_ibl)

species_ALS <- unique(ALS_2015$species)

length(unique(ALS_2015$plot_id))
length(unique(FID_2015_2019$plotid))

# fix ALS plot IDs to match FID format (KS001)
ALS_2015 <- ALS_2015 %>%
  mutate(plot_id = sprintf("KS%03d", as.numeric(plot_id)))

# contingency variables column names between FID and ALS using the FID
ALS_2015 <- ALS_2015 %>%
  rename(
    plotid = plot_id,
    treeid = tree_id
  ) %>%
  select(plotid, treeid, x, y, species, dbh, h, dead)

head(ALS_2015)
str(ALS_2015)

###    SAVE MATCHING CONTROL PLOTS FILE - NOT REPEAT IT CAREFUL      ###########

# check mismatch FID-ALS
fid_not_als <- setdiff(unique(FID_2015$plotid), unique(ALS_2015$plotid))
# check mismatch FID-ALS
als_not_fid <- setdiff(unique(ALS_2015$plotid), unique(FID_2015$plotid))

# matching plot IDs The most IMPORTANT
matching_ids <- intersect(unique(FID_2015$plotid), unique(ALS_2015$plotid))


#4------------------------------------------------------------------------------
# write to txt files for save the info
writeLines(matching_ids,    file.path(raw_process_data_path, "matching_plots.txt"))
writeLines(fid_not_als, file.path(raw_process_data_path, "FID_not_in_ALS.txt"))
writeLines(als_not_fid, file.path(raw_process_data_path, "ALS_not_in_FID.txt"))



#5------------------------------------------------------------------------------
# Let's create different data frame for alive and dead trees as in FID

# filter both datasets to matching the plots in 2015
FID_match <- FID_2015 %>% filter(plotid %in% matching_ids)
ALS_match <- ALS_2015 %>% filter(plotid %in% matching_ids)

#FID 2015 vs ALS differences
setdiff(unique(FID_2015$plotid), unique(ALS_match$plotid))
setdiff(unique(ALS_match$plotid), unique(FID_2015$plotid))


# define common plot set across ALL datasets ALS 2015 - FID 2015 and 2019
valid_plots <- Reduce(intersect, list(
  unique(FID_2015$plotid),
  unique(FID_2019$plotid),
  unique(ALS_match$plotid)
))

# filter ALL datasets to same plots
FID_2015_clean <- FID_2015 %>% filter(plotid %in% valid_plots)
FID_2019_clean <- FID_2019 %>% filter(plotid %in% valid_plots)
ALS_clean      <- ALS_match %>% filter(plotid %in% valid_plots)

# ALS split 
ALS_clean_alive <- ALS_clean %>% filter(dead == 0)
ALS_clean_dead  <- ALS_clean %>% filter(dead == 1)

# FID splits (2015)
FID_2015_clean_alive <- FID_2015_clean %>% filter(dead == 0)
FID_2015_clean_dead  <- FID_2015_clean %>% filter(dead == 1)

# FID splits (2019)
FID_2019_clean_alive <- FID_2019_clean %>% filter(dead == 0)
FID_2019_clean_dead  <- FID_2019_clean %>% filter(dead == 1)



#6------------------------------------------------------------------------------
# CSV exports
write.csv(ALS_2015, file.path(raw_process_data_path, "ALS_2015.csv"), row.names = FALSE)

write.csv(FID_2015_clean, file.path(raw_process_data_path, "FID_2015_clean.csv"), row.names = FALSE)
write.csv(FID_2019_clean, file.path(raw_process_data_path, "FID_2019_clean.csv"), row.names = FALSE)
write.csv(ALS_clean,      file.path(raw_process_data_path, "ALS_clean.csv"), row.names = FALSE)

write.csv(FID_2015_clean_alive, file.path(raw_process_data_path, "FID_2015_clean_alive.csv"), row.names = FALSE)
write.csv(FID_2015_clean_dead,  file.path(raw_process_data_path, "FID_2015_clean_dead.csv"), row.names = FALSE)

write.csv(FID_2019_clean_alive, file.path(raw_process_data_path, "FID_2019_clean_alive.csv"), row.names = FALSE)
write.csv(FID_2019_clean_dead,  file.path(raw_process_data_path, "FID_2019_clean_dead.csv"), row.names = FALSE)

write.csv(ALS_clean_alive, file.path(raw_process_data_path, "ALS_clean_alive.csv"), row.names = FALSE)
write.csv(ALS_clean_dead,  file.path(raw_process_data_path, "ALS_clean_dead.csv"), row.names = FALSE)

# RDS exports
saveRDS(ALS_2015, file.path(raw_process_data_path, "ALS_2015.rds"))

saveRDS(FID_2015_clean, file.path(raw_process_data_path, "FID_2015_clean.rds"))
saveRDS(FID_2019_clean, file.path(raw_process_data_path, "FID_2019_clean.rds"))
saveRDS(ALS_clean,      file.path(raw_process_data_path, "ALS_clean.rds"))

saveRDS(FID_2015_clean_alive, file.path(raw_process_data_path, "FID_2015_clean_alive.rds"))
saveRDS(FID_2015_clean_dead,  file.path(raw_process_data_path, "FID_2015_clean_dead.rds"))

saveRDS(FID_2019_clean_alive, file.path(raw_process_data_path, "FID_2019_clean_alive.rds"))
saveRDS(FID_2019_clean_dead,  file.path(raw_process_data_path, "FID_2019_clean_dead.rds"))

saveRDS(ALS_clean_alive, file.path(raw_process_data_path, "ALS_clean_alive.rds"))
saveRDS(ALS_clean_dead,  file.path(raw_process_data_path, "ALS_clean_dead.rds"))


#7------------------------------------------------------------------------------
# FINAL TEST TO SEE IF EVERYTHING IS FINE

# TEST 1 — identical plot sets across all data sets
setequal(unique(FID_2015_clean$plotid),
         unique(FID_2019_clean$plotid))

setequal(unique(FID_2015_clean$plotid),
         unique(ALS_clean$plotid))

setequal(unique(FID_2019_clean$plotid),
         unique(ALS_clean$plotid))


# TEST 2 — random sample check - same plots exist everywhere
set.seed(123)
sample_plots <- sample(unique(FID_2015_clean$plotid), 10)

sample_plots %in% ALS_clean$plotid
sample_plots %in% FID_2019_clean$plotid


# !! TEST 3 — row-level consistency - counts per plot
test_counts <- data.frame(
  plotid = sample_plots,
  FID2015 = sapply(sample_plots, function(x) sum(FID_2015_clean$plotid == x)),
  FID2019 = sapply(sample_plots, function(x) sum(FID_2019_clean$plotid == x)),
  ALS     = sapply(sample_plots, function(x) sum(ALS_clean$plotid == x))
)

test_counts




#-------------------------------------------------------------------------------
###   FINAL PART - MERGING DATA FRAMES      ##############
# Now let's upload also the plot information table from FID_2015_2019
# And let's make the final match also for the tree ID and species

#8------------------------------------------------------------------------------
# Load Field Inventory Data general information
FID_2015_2019_Info <- read_excel(
  file.path(raw_metadata_path, "bialowieza_field_2015-2019.xlsx"),
  sheet = "plots_info")

FID_2015_2019_Info_clean <- FID_2015_2019_Info %>% filter(plotid %in% valid_plots)

# Change the species names as above
FID_2015_2019_Info_clean <- FID_2015_2019_Info_clean %>%
  left_join(species_ref %>% rename(domspecies2015 = sp_ibl,
                                   sp_name_2015 = sp_name,
                                   species_2015 = species),
            by = "domspecies2015") %>%
  left_join(species_ref %>% rename(domspecies2019 = sp_ibl,
                                   sp_name_2019 = sp_name,
                                   species_2019 = species),
            by = "domspecies2019") %>%
  relocate(sp_name_2015, species_2015, .after = domspecies2015) %>%
  relocate(sp_name_2019, species_2019, .after = domspecies2019)


####    TEST OF PLOTS NUMBERS AND TREE COUNTS     #############

# 1) total unique plots per dataset
data.frame(
  dataset = c("FID_2015_clean", "FID_2019_clean", "ALS_clean",
              "FID_2015_clean_alive", "FID_2015_clean_dead",
              "FID_2019_clean_alive", "FID_2019_clean_dead",
              "ALS_clean_alive", "ALS_clean_dead",
              "FID_2015_2019_Info_clean"),
  n_unique_plots = c(
    length(unique(FID_2015_clean$plotid)),
    length(unique(FID_2019_clean$plotid)),
    length(unique(ALS_clean$plotid)),
    length(unique(FID_2015_clean_alive$plotid)),
    length(unique(FID_2015_clean_dead$plotid)),
    length(unique(FID_2019_clean_alive$plotid)),
    length(unique(FID_2019_clean_dead$plotid)),
    length(unique(ALS_clean_alive$plotid)),
    length(unique(ALS_clean_dead$plotid)),
    length(unique(FID_2015_2019_Info_clean$plotid))
  )
)

# 2) per-plot tree counts (example: FID 2015)
table(FID_2015_clean$plotid)

# 3) per-plot tree counts comparison
plot_count_check <- data.frame(
  plotid = sort(unique(FID_2015_clean$plotid)),
  FID2015 = as.integer(table(FID_2015_clean$plotid)),
  FID2019 = as.integer(table(FID_2019_clean$plotid)),
  ALS     = as.integer(table(ALS_clean$plotid))
)

plot_count_check

# SAVE THE NEW DATA FRAME ABOUT INFO

#-------------------------------------------------------------------------------
# save path
out <- raw_process_data_path

#-------------------------------------------------------------------------------
# save main dataframe
write.csv(FID_2015_2019_Info_clean,
          file.path(out, "FID_2015_2019_Info_clean.csv"),
          row.names = FALSE)

saveRDS(FID_2015_2019_Info_clean,
        file.path(out, "FID_2015_2019_Info_clean.rds"))

#-------------------------------------------------------------------------------
# save plot count summary table
plot_summary <- data.frame(
  dataset = c("FID_2015_clean", "FID_2019_clean", "ALS_clean",
              "FID_2015_clean_alive", "FID_2015_clean_dead",
              "FID_2019_clean_alive", "FID_2019_clean_dead",
              "ALS_clean_alive", "ALS_clean_dead",
              "FID_2015_2019_Info_clean"),
  n_unique_plots = c(
    length(unique(FID_2015_clean$plotid)),
    length(unique(FID_2019_clean$plotid)),
    length(unique(ALS_clean$plotid)),
    length(unique(FID_2015_clean_alive$plotid)),
    length(unique(FID_2015_clean_dead$plotid)),
    length(unique(FID_2019_clean_alive$plotid)),
    length(unique(FID_2019_clean_dead$plotid)),
    length(unique(ALS_clean_alive$plotid)),
    length(unique(ALS_clean_dead$plotid)),
    length(unique(FID_2015_2019_Info_clean$plotid))
  )
)

write.csv(plot_summary,
          file.path(out, "plot_summary.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# save per-plot counts table
write.csv(plot_count_check,
          file.path(out, "plot_count_check.csv"),
          row.names = FALSE)

saveRDS(plot_count_check,
        file.path(out, "plot_count_check.rds"))


#9imp---------------------------------------------------------------------------
# Now let's merge the ALS with the FID 2015 to match the trees IDs

NOT WORKING WELL


# Let's join the 



# Important information column structure for init table in iLand: x y species dbh height age