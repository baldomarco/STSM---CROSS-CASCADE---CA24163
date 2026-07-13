# Load required libraries
library(tidyr)
library(dplyr)
library(RSQLite)
library(vegan)
library(fields)

rm(list=ls()) 

#setwd("C:/iLand/2026/Bialowieza/project_files/ALS_V1/DB_V1/")
setwd("C:/iLand/2026/Bialowieza/output/")

# Path to the directory containing your SQLite databases
#dataroot <- "C:/iLand/2026/Bialowieza/project_files/ALS_V1/DB_V1/"
dataroot <- "C:/iLand/2026/Bialowieza/output/"

{# Create an empty list to store data frames
  dfs <- list() # not working for several subset, only one
}


# Create the list of dataframes I want to save
abeStand_scen <- c()

bb_scen <- c()

tree_scen <- c()

stand_scen <- c()

lnd_scen <-c()

#carbon_scen <- c()

damage.all <- c()

H_BA_heterogenity_scen <- c()

variables.all <- c()

plot_variables_all <- c()

#-------------------------------------------------------------------------------
# import the list of files within the folder in dataroot with .sqlite extension

database_files <- list.files(path = dataroot, pattern = ".sqlite", full.names = FALSE)                              # alternative way to select all the databases within a folder select all the file with sqlite format
database_files <- list.files(dataroot, ".sqlite")  

for (i in (1:length(database_files)))  {    # We read in the files in the loop. The "i" is for the x from 1 to i lenght of the dataset of files 
  
  #i<-1
  # i<-1 # to test but remember to don't run also the }
  
  
  # Name of the database
  file <-paste0(database_files[i])                                    # File to read here the case is always the actual case in the loop
  # "file"= name of the object, "paste0"+ function to create a NAME for a computer path of selection of data/objects
  
  # case<- strsplit(all_v[i],".s")[[1]][1]      # Why we used it in the paper analysis? # explanation of this function is needed
  case <- tools::file_path_sans_ext(database_files[i]) 
  
  
  # Control
  print(file)
  
  # connect to the database 
  sqlite.driver <- dbDriver("SQLite")
  db <- dbConnect(sqlite.driver, dbname = file)  # connect to the file
  tables.in.the.file <- dbListTables(db)         # explore the tables in the file
  print(tables.in.the.file)
  
  #-----------------------------------------------------------------------------
  # LOAD THE DATABASE # Read tables from the database
  abeStand <- dbReadTable(db, "abeStand")
  barkbeetle <- dbReadTable(db, "barkbeetle")
  #carbon <- dbReadTable(db, "carbon")
  landscape <- dbReadTable(db, "landscape")
  stand <- dbReadTable(db, "stand")
  tree <- dbReadTable(db, "tree")
  
  dbDisconnect(db)    # close the file
  
  # TO UNDERSTAND THE OPERATORS %>% AND %IN% HAVE TO STUDY THEM IN DATACAMP AND IN DPLYR CRAN PACKAGES
  
  #-----------------------------------------------------------------------------
  # CREATE SHANNON VARIABLE
  
  # CREATE A TABLE FOR THE ANNUAL VOL,BA,COUNT DATA
  annual.data<-landscape %>% 
    group_by(year) %>% 
    filter(year>=0) %>% 
    summarize(VOL.tot=sum(volume_m3), BA.tot=sum(basal_area_m2), count.tot=sum(count_ha))
  
  # SAME BUT PER SPECIES
  annual.spec.data<-landscape %>% 
    group_by(year, species) %>%  
    filter(year>=0) %>% 
    summarize(VOL=(volume_m3), BA=(basal_area_m2), count=sum(count_ha))
  
  # print(head(annual.data))
  # print(head(annual.spec.data))
                       
  
  # BECARFUL THE NUMBERS OF VALUES CAN BE LESS THEN THE NUMBER OF YEARS IN YOUNG STANDS.
  S <- landscape %>%
    filter(year >= 0) %>%
    group_by(year) %>%
    summarise(
      n = n_distinct(species[volume_m3 > 0]),
      .groups = "drop"
    )
  
  S <- tibble(year = 0:max(landscape$year)) %>%
    left_join(S, by = "year") %>%
    mutate(n = replace_na(n, 1))
  
  # number of species in each year  (added the filter to count non-zero volumes, now it is okay)
  
  #print(S$n)
  #length(S$n)
  
  t<-annual.spec.data %>% right_join(annual.data, by="year")  %>% 
    mutate(prop.VOL=VOL/VOL.tot, prop.BA=BA/BA.tot, prop.count=count/count.tot) %>% 
    filter(prop.VOL>0)  # here I also filtering them to have only records with m3>0
  
  #https://www.statology.org/shannon-diversity-index/
  #https://www.statology.org/shannon-diversity-index-calculator/
  
  # Shannon diversity index (SDI): this is already Shannon.... that extra step by dividing the by S$n is making equitability index based on the link above.
  # so maybe we can make shannon based on BA, VOL and number of trees.... (? can discuss or save all 3 and will see...)
  
  # Calculate H.BA
  H.BA <- t %>%
    group_by(year) %>%
    summarize(H = -1 * sum(prop.BA * log(prop.BA)))
  
  # Calculate H.VOL
  H.VOL <- t %>%
    group_by(year) %>%
    summarize(H = -1 * sum(prop.VOL * log(prop.VOL)))
  
  # Calculate H.count
  H.count <- t %>%
    group_by(year) %>%
    summarize(H = -1 * sum(prop.count * log(prop.count)))
  
  # Find the maximum year across all dataframes
  max_year <- max(
    max(H.BA$year, na.rm = TRUE),
    max(H.VOL$year, na.rm = TRUE),
    max(H.count$year, na.rm = TRUE)
  )
  
  # Create a reference dataframe with all years from 1 to the maximum year
  reference_df <- tibble(year = 0:max(H.BA$year, H.VOL$year, H.count$year))
  
  
  # Left join each dataframe with the reference dataframe and replace NA with 0
  H.BA <- reference_df %>%
    left_join(H.BA, by = "year") %>%
    replace(is.na(.), 0)
  
  H.VOL <- reference_df %>%
    left_join(H.VOL, by = "year") %>%
    replace(is.na(.), 0)
  
  H.count <- reference_df %>%
    left_join(H.count, by = "year") %>%
    replace(is.na(.), 0)
  
  # here I just put the proportion of number of trees
  set.panel(3,1)
  plot(H.BA$year, H.BA$H, type="l")
  plot(H.VOL$year, H.VOL$H, type="l")
  plot(H.count$year, H.count$H, type="l")
  
  #The higher the value of H, the higher the diversity of species in a particular community. 
  #The lower the value of H, the lower the diversity. 
  #A value of H = 0 indicates a community that only has one species.
  
  
  # Shannon Equitability Index (SEI).... maybe we do not need this at all.....
  SEI<-H.BA$H/log(S$n)
  SEI[which(is.na(SEI)==T)]<-0
  
  
  # Pioneer species proportion based on BA
  
  earlyspecs<-c("lade","acps","frex","cabe","bepe","alin","algl","acca","acpl","soau","soar","coav","alvi","potr","poni","ulgl","saca","rops")
  
  early.spec.prop<-t %>% filter(species %in% earlyspecs) %>% summarize(BA=sum(prop.BA))
  
  esp<-data.frame(year=c(0:max(landscape$year))) # FILL THE MISSING YEARS WITH ZERO
  esp<-esp %>% left_join(early.spec.prop, by="year")
  esp$BA[which(is.na(esp$BA)==T)]<-0
  
  #-----------------------------------------------------------------------------
  # Create Shannon on BA heterogeneity with VEGAN
  
  # Calculate the Shannon diversity index
  H_BA_heterogenity <- t %>%
    group_by(year) %>%
    summarize(shannon_ba_heterog = diversity(prop.BA, base = exp(1)))
  
  # Print the resulting dataframe
  print(H_BA_heterogenity)
  
  # PLOT THE EARLY SPECIES PROPORTION AND THE SHANNON EQUITABILITY INDEX
  set.panel(3,1)
  plot(SEI)
  plot(esp$year,esp$BA)
  plot(H_BA_heterogenity$year, H_BA_heterogenity$shannon_ba_heterog)
  
  
  
  #-----------------------------------------------------------------------------
  # CREATE THE DATA FRAME FOR ADD VARIABLES ABOUT CARBON IN THE FINAL DATA FRAME
  
  #ab.tot.c<- data.frame(carbon %>% 
   #                       group_by(year) %>% 
    #                      summarise(tot_carbon=sum(stem_c, branch_c, foliage_c, coarseRoot_c, fineRoot_c, regeneration_c, snags_c, snagsOther_c, downedWood_c, litter_c, soil_c)))
  
  
  #-----------------------------------------------------------------------------
  # CREATE A NEW DATA FRAME TO WORK ON THE SPECIFIC VARIABLES NEEDED 
  
  head(landscape)
  
  # TO UNDERSTAND THE OPERATORS %>% AND %IN% HAVE TO STUDY THEM IN DATACAMP AND IN DPLER CRAN PACKAGES
  
  ab.lnd.v<- data.frame(landscape %>% 
                          group_by(year) %>% 
                          filter(year>=0) %>%
                          summarise(tot_volume=sum(volume_m3),living_c=sum(total_carbon_kg), count_ha=sum(count_ha), tot_ba=sum(basal_area_m2),npp=sum(NPP_kg), LAI=sum(LAI), sapling=sum(cohort_count_ha),growth_m3=sum(gwl_m3)))
  
  
  #-------------------------------------------------------------------------------
  # Create dynamic stand statistics from individual tree output
  
  dynamicstand_1 <- tree %>%
    group_by(year) %>%
    summarise(
      height_mean = mean(height, na.rm = TRUE),
      dbh_mean    = mean(dbh, na.rm = TRUE),
      age_mean    = mean(age, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Fill missing years so the dataframe matches carbon outputs
  dynamicstand_1 <- tibble(
    year = H.VOL$year
  ) %>%
    left_join(dynamicstand_1, by = "year") %>%
    mutate(
      height_mean = tidyr::replace_na(height_mean, 0),
      dbh_mean    = tidyr::replace_na(dbh_mean, 0),
      age_mean    = tidyr::replace_na(age_mean, 0)
    )
  
  #-------------------------------------------------------------------------------
  # Create variables data frame
  
  variables <- data.frame(
    case      = case,
    year      = H.VOL$year,
    h         = dynamicstand_1$height_mean,
    dbh       = dynamicstand_1$dbh_mean,
    age       = dynamicstand_1$age_mean,
    H.count   = H.count$H,
    H.VOL     = H.VOL$H,
    H.BA      = H.BA$H,
    SEI       = SEI
  )
  
  #-------------------------------------------------------------------------------
  # Add landscape variables
  
  variables <- variables %>%
    inner_join(ab.lnd.v, by = "year") 
  
  head(variables)
  
  
  # ---- Let's create the BDV predictors variables ----
  
  # Preventing the row jump if 0
  all_years <- tibble(year = sort(unique(landscape$year)))
  
  #-------------------------------------------------------------------------------
  # Create the LAI time series
  
  LAI <- landscape %>%
    group_by(year) %>%
    summarise(LAI = sum(LAI), .groups = "drop") %>%
    right_join(all_years, by = "year") %>%
    mutate(LAI = tidyr::replace_na(LAI, 0))
  
  # EDIT IT FOR THE PREDICTORS BA BROADLEAF - TREE_10_40 - BROAD>40
  
  #-------------------------------------------------------------------------------
  # Calculate the number of trees between 10 and 40 cm
  
  # Specify the range for DBH
  dbh_min <- 10
  dbh_max <- 40
  
  # Calculate the number of trees in each year with DBH between 10 and 40 cm
  tree_10_40 <- tree %>%
    filter(dbh >= dbh_min, dbh <= dbh_max) %>%
    group_by(year) %>%
    summarise(tree_10_40 = n(), .groups = "drop") %>%
    right_join(all_years, by = "year") %>%
    mutate(tree_10_40 = tidyr::replace_na(tree_10_40, 0))
  
  
  #-------------------------------------------------------------------------------
  # To define the species to be removed
  unique_sp <- unique(landscape$species)  # alternative unique_plots <- unique(CZ_JH1[,"plotID"])
  
  
  species_to_remove <- c("piab", "pisy", "abal",
                         "lade", "psme", "pini", 
                         "pice")
  
  # Use subset to filter the dataframe
  filtered_broadl <- subset(landscape, !species %in% species_to_remove)
  
  # Sum BA for every broadleaf species in every year
  summed_broadl_ba <- filtered_broadl %>%
    group_by(year) %>%
    summarise(basal_area_m2 = sum(basal_area_m2), .groups = "drop") %>%
    right_join(all_years, by = "year") %>%
    mutate(basal_area_m2 = tidyr::replace_na(basal_area_m2, 0))
  
  #-------------------------------------------------------------------------------
  # Calculate the basal area only of the broadleave with a dbh > 40cm
  
  # Perform left join with the summarization result
  # Create the broadl_40 summary with continuous years and fill missing values with 0
  broadl_40 <- tree %>%
    filter(dbh > 40 & !species %in% species_to_remove) %>%
    group_by(year) %>%
    summarise(broadl_40 = n(), .groups = "drop") %>%
    right_join(all_years, by = "year") %>%
    mutate(broadl_40 = tidyr::replace_na(broadl_40, 0))
  
  #-------------------------------------------------------------------------------
  # Age - comes from the avarage age
  
  # Summarize the mean age for each year and round the ages
  stand_age <- tree %>%
    group_by(year) %>%
    summarise(age = round(mean(age, na.rm = TRUE)), .groups = "drop") %>%
    right_join(all_years, by = "year") %>%
    mutate(age = tidyr::replace_na(age, 0))
  
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # Merge the data frames Plot L1_10
  plot_L1_10_df_simul <- bind_cols(stand_age, 
                                   lai_sim = LAI$LAI,
                                   ba_broadl = summed_broadl_ba$basal_area_m2,
                                   trees_10_40 = tree_10_40$tree_10_40,
                                   broadl_40 = broadl_40$broadl_40)
  
  
  
  #-----------------------------------------------------------------------------
  # CREATE THE CALCULATION FOR DAMAGES LOOK LINE 370
  
  landscape.area<-landscape$area[1]                                             # CREATE THE VARIABLE FOR LANDSCAPE AREA          
  
  lnd_volume = landscape %>% group_by(year)  %>%                                # CREATE THE SUMMARIZATION OF THE SPECIES VOLUME PROPORTION TO CREATE A TOTAL LANDSCAPE VOLUME
    summarise(tot_vol = sum(volume_m3),
              .groups = 'drop')
  
  head(lnd_volume)                                                              
  
  # wind and barkbeetle merging:
  head(barkbeetle)
  
  # CREATE THE DATA FRAME FOR DAMAGE OF BARKBEETLE
  damage <- data.frame(year = barkbeetle$year, 
                       barkbeetle = barkbeetle$killedVolume, 
                       case = case)
  
  # ADD THE LANDSCAPE VOLUME IN THE DAMAGE DATA FRAME
  damage <- left_join(damage, lnd_volume, by = "year")
  
  # GIVE THE NAME TO EVERY VARIABLE
  colnames(damage) <- c("year", "barkbeetle", "case", "volume")
  
  head(damage)
  
  
  #------------------------------                                               # CLOSE THE litte FOR CYCLE
  
  # Collect abeUnit data FOR CREATE THE VARIABLE ABEUNIT FOR ALL THE RUNS
  #landscape_removed <-(landscape_removed %>% mutate(run=case))
  #landscape_removed_scen<-rbind(landscape_removed_scen, landscape_removed)
  
  # Collect abeStand data for create the variable abeStand for all the run
  abeStand <-(abeStand %>% mutate(run=case))
  abeStand_scen <-rbind(abeStand_scen, abeStand)
  
  # Collect landscape data FOR CREATE THE VARIABLE LND FOR ALL THE RUNS
  landscape <- (landscape %>% mutate(run=case))
  lnd_scen <- rbind(lnd_scen, landscape)
  
  # Collect dynamicstand data FOR CREATE THE VARIABLE WIND FOR ALL THE RUNS
  tree <-(tree %>% mutate(run=case))
  tree_scen <-rbind(tree_scen, tree)
  
  # Collect dynamicstand data FOR CREATE THE VARIABLE WIND FOR ALL THE RUNS
  stand <-(stand %>% mutate(run=case))
  stand_scen <-rbind(stand_scen, stand)

  # Collect barkbeetle data FOR CREATE THE VARIABLE BB FOR ALL THE RUNS
  barkbeetle <-(barkbeetle %>% mutate(run=case))
  bb_scen <-rbind(bb_scen, barkbeetle)

  # CREATE THE VARIABLE DAMAGE FOR ALL THE RUNS
  damage <-(damage %>% mutate(run=case))
  damage.all<-rbind(damage.all, damage)                                         # PUT ALL THE DAMAGE RUNS INTO A SINGLE DATAFRAME WITH DIFF CASES TO BE PLOT ALL TOGETHER IN LINE 370
  
  # CREATE THE VARIABLES DF FOR ALL THE RUNS
  variables <-(variables %>% mutate(run=case))
  variables.all<-rbind(variables.all, variables) 
  
  # CREATE THE VARIABLE H of BA per species FOR ALL THE RUNS
  H_BA_heterogenity <-(H_BA_heterogenity %>% mutate(run=case))
  H_BA_heterogenity_scen <- rbind(H_BA_heterogenity_scen, H_BA_heterogenity)
  
  # CREATE THE VARIABLES NEEDED IN THE BDV STUDY
  plot_L1_33_df_simul <-(plot_L1_10_df_simul %>% mutate(run=case))
  plot_variables_all <- rbind(plot_variables_all, plot_L1_33_df_simul)
  
}                                                                               # CLOSE THE WHOLE FOR CYCLE


# ---- SAVE THE DATA IN RDS ----
# Output directory
dataroot <- "C:/iLand/2026/Bialowieza/RDS/"

if (!dir.exists(dataroot)) {
  dir.create(dataroot, recursive = TRUE)
}

# Save all aggregated data frames
saveRDS(abeStand_scen, file.path(dataroot, "abeStand_scen_FID_V2_FULL.rds"))
saveRDS(lnd_scen, file.path(dataroot, "lnd_scen_FID_V2_FULL.rds"))
saveRDS(tree_scen, file.path(dataroot, "tree_scen_FID_V2_FULL.rds"))
saveRDS(stand_scen, file.path(dataroot, "stand_scen_FID_V2_FULL.rds"))
saveRDS(bb_scen, file.path(dataroot, "bb_scen_FID_V2_FULL.rds"))
saveRDS(damage.all, file.path(dataroot, "damage_all_FID_V2_FULL.rds"))
saveRDS(variables.all, file.path(dataroot, "variables_all_FID_V2_FULL.rds"))
saveRDS(H_BA_heterogenity_scen, file.path(dataroot, "H_BA_heterogenity_scen_FID_V2_FULL.rds"))
saveRDS(plot_variables_all, file.path(dataroot, "plot_variables_all_FID_V2_FULL.rds"))

cat("All .rds files saved to:", dataroot, "\n")


library(ggplot2)
library(GGally)
library(corrplot)


# Start with plots

#_______________________________________________
library(ggplot2)
library(gridExtra) # To arrange the graphs in a grid

# NEED TO OPEN A PDF WRITER AND GIVE IT THE ROOT, THE NAME, AND THE SIZE
dataroot <- "C:/iLand/2026/Bialowieza/Figure/"
pdf(paste0(dataroot, "20260708_FID_V1_16_plots.pdf"), height=8, width=12)


#-------------------------------------------------------------------------------
# This tells the colors:

species.we.have<-unique(lnd_scen$species)                                            # IT IS SAYING WHICH SPECIES WE HAVE IN THE DATABASE IN THIS CASE LANDSCAPE


# LIST OF ALL POSSIBLE SPECIES

cols.all=c( "rops"="#e0e0e0", "acpl"="#A9A9A9",   "alin"="#696969", "alvi"="#2e2e2e",
            "bepe"="#fadfad", "bepu"= "#5F9EA0", "prse" = "#8B2323",
            "casa"="#7eeadf", "coav"="#20c6b6",  
            "tipl"="#645394", "ulgl"="#311432" ,
            "saca"="#D8BFD8",  "soar"="#DDA0DD", "soau"="#BA55D3",
            "pice"="#D27D2D", "pini"="#a81c07",
            "algl"="#2ECBE9","tico"="#128FC8",  "potr"="#00468B","poni"="#5BAEB7",
            "frex"="#fe9cb5","cabe"="#fe6181","acps"="#fe223e",
            "lade"="#FFFE71","abal"="#FFD800", "pisy"="#A4DE02",
            "fasy"="#76BA1B", "piab"="#006600",
            "quro"="#FF7F00", "qupe"="#FF9900", "qupu"="#CC9900" 
)


# COLORATION ORDER FOR ALL THE POSSIBLE SPECIES

new_order_gg.all=c("alvi","alin", "acpl", "rops","bepe" ,"bepu", "prse", "coav", "casa", "ulgl", "tipl",  "soau", "soar", "saca",  "pini", "pice",
                   "poni", "algl", "tico", "potr",  "frex","cabe", "acps",  "lade", "abal",  "qupu", "qupe","quro","pisy", "fasy", "piab")


# This will show at the end only the species we really have on the landscape. 

cols<-cols.all[names(cols.all) %in% species.we.have]
new_order_gg<- new_order_gg.all[new_order_gg.all %in% species.we.have]

# STARTING PLOTS

# Make a plot with ggplot, volume, colored by species for the transitional period for Clear cut management system
#-------------------------------------------------------------------------------
# PLOT LANDSCAPE VOLUME PLOT FOR CASES (GEOM AREA)

g1 <- ggplot(lnd_scen, aes(year, volume_m3, fill=factor(species, levels=new_order_gg)))+
  geom_area() +
  scale_fill_manual(values=cols[new_order_gg], guide=guide_legend(reverse=TRUE))+
  ggtitle("Landscape Volume by species")+
  facet_wrap(~run, ncol=2)+
  labs(x = "Year",y="Volume [m3/ha]",fill = "Species")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylim(0,800)+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )


# Create a bar chart of volume by species
ggplot(lnd_scen, aes(x = year, y = volume_m3, fill = factor(species, levels = new_order_gg))) +
  geom_bar(stat = "identity") +  # Change from geom_area to geom_bar
  scale_fill_manual(values = cols[new_order_gg], guide = guide_legend(reverse = TRUE)) +
  ggtitle("Landscape Volume by Species") +
  facet_wrap(~run, ncol = 2) +
  labs(x = "Year", y = "Volume [m3/ha]", fill = "Species") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0, 400) +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )

#-------------------------------------------------------------------------------
# PLOT BASAL AREA GEOM_LINE AT LANDSCAPE LEVEL BY SPECIES SELECTED

# SPECIES specificaly BA:

species.to.keep<-c("piab", "algl","bepe", "pisy", "quro")


lnd_scen2 <- lnd_scen %>% filter(species %in% species.to.keep)

b1 <- ggplot(data=lnd_scen2, aes(x=year, y=basal_area_m2, color=species)) + 
  geom_line(size=1.2)+
  ggtitle("Basal area by species") +
  facet_wrap(~run, ncol=4)+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("Basal area [m2/ha]")+  
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )

#-------------------------------------------------------------------------------
# PLOT TOTAL AVG BASAL AREA AT LANDSCAPE LEVEL BY SPECIES

g3 <- ggplot(lnd_scen, aes(year, basal_area_m2, fill=factor(species, levels=new_order_gg)))+
  geom_area() +
  scale_fill_manual(values=cols[new_order_gg], guide=guide_legend(reverse=TRUE))+
  ggtitle("Total Basal Area")+
  facet_wrap(~run, ncol=2)+
  labs(x = "Year",y="Basal Area [m2/ha]",fill = "Species")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )


#-------------------------------------------------------------------------------
# PLOT DBH GEOM_LINE AT LANDSCAPE LEVEL BY SPECIES

species.to.keep<-c("piab", "algl","quro", "bepe")


lnd_scen2 <- lnd_scen %>% filter(species %in% species.to.keep)

b3 <- ggplot(data=lnd_scen2, aes(x=year, y=dbh_avg_cm, color=species)) + 
  geom_line(size=1.2)+
  ggtitle("Avarage DBH by species") +
  facet_wrap(~run, ncol=4)+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("DBH [cm]")+  
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )

#-------------------------------------------------------------------------------
# PLOT NUMBER OF STEMS GEOM_LINE AT LANDSCAPE LEVEL BY SPECIES

# lnd2 <- lnd %>% filter(species %in% species.to.keep)    --- ADD IF YOU WANT SELECT SPECIES TO KEEP

b4 <- ggplot(data=lnd_scen2, aes(x=year, y=count_ha, color=species)) + 
  geom_line(size=1.2)+
  ggtitle("N. individual stems by species") +
  facet_wrap(~run, ncol=4)+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("Individual stems")+  
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )


#-------------------------------------------------------------------------------
# PLOT NUMBER OF STEMS GEOM_AREA AT LANDSCAPE LEVEL BY SPECIES

g5 <- ggplot(lnd_scen, aes(x=year, y=count_ha, fill=factor(species, levels=new_order_gg)))+ 
  geom_area(size=1.2)+
  scale_fill_manual(values=cols[new_order_gg], guide=guide_legend(reverse=TRUE))+
  ggtitle("N. individual stems by species") +
  facet_wrap(~run, ncol=2)+
  labs(x = "Year",y="Individual Stems", fill = "Species")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )

# Create a bar chart of n of trees by species
ggplot(lnd_scen, aes(x = year, y = count_ha, fill = factor(species, levels = new_order_gg))) +
  geom_bar(stat = "identity") +  # Change from geom_area to geom_bar
  scale_fill_manual(values = cols[new_order_gg], guide = guide_legend(reverse = TRUE)) +
  ggtitle("N. individual stems by species") +
  facet_wrap(~run, ncol = 4) +
  labs(x = "Year", y = "Individual Stems", fill = "Species") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )
#-------------------------------------------------------------------------------
# PLOT HEIGHT GEOM_LINE AT LANDSCAPE LEVEL BY SPECIES

# lnd2 <- lnd %>% filter(species %in% species.to.keep)    ---   ADD IF YOU WANT SELECT SPECIES TO KEEP

b5 <- ggplot(data=lnd_scen2, aes(x=year, y=height_avg_m, color=species)) + 
  geom_line(size=1.2)+
  ggtitle("Avarage Height by species") +
  facet_wrap(~run, ncol=4)+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("Height [m]")+  
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )


#-------------------------------------------------------------------------------
# Total Carbon in Kg (total_carbon_kg	double	total carbon in living biomass (aboveground compartments and roots) of all living trees (including regeneration layer) (kg/ha))

g8 <- ggplot(lnd_scen, aes(year, total_carbon_kg, fill=factor(species, levels=new_order_gg)))+
  geom_area() +
  scale_fill_manual(values=cols[new_order_gg], guide=guide_legend(reverse=TRUE))+
  ggtitle("Total Carbon in Living Biomass")+
  facet_wrap(~run, ncol=2)+
  labs(x = "Year",y="[kg/ha]",fill = "Species")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )

#-------------------------------------------------------------------------------
# PLOT LAI AT LANDSCAPE LEVEL BY SPECIES

g9 <- ggplot(lnd_scen, aes(year, LAI, fill=factor(species, levels=new_order_gg)))+
  geom_area() +
  scale_fill_manual(values=cols[new_order_gg], guide=guide_legend(reverse=TRUE))+
  ggtitle("LAI index by species")+
  facet_wrap(~run, ncol=4)+
  labs(x = "Year",y="LAI index",fill = "Species")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )


#-------------------------------------------------------------------------------
# PLOT NPP AT LANDSCAPE LEVEL BY SPECIES

g10 <- ggplot(lnd_scen, aes(year, NPP_kg, fill=factor(species, levels=new_order_gg)))+
  geom_area() +
  scale_fill_manual(values=cols[new_order_gg], guide=guide_legend(reverse=TRUE))+
  ggtitle("Net Primary Productivity")+
  facet_wrap(~run, ncol=4)+
  labs(x = "Year",y="NPP [kg/ha]",fill = "Species")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )


#grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10, ncol=10)

#grid.arrange(b1,b2,b3,b4,b5,ncol=5)

#grid.arrange(g1,g5,g10, ncol=3)

#grid.arrange(g1,g5,g4, ncol=3)




# Shannon

H.count <- ggplot(variables.all, aes(x=year, y=H.count))+
  geom_line() +
  facet_wrap(~case, ncol=2)+
  ggtitle("Shannon Index on numbers of individuals per species")+
  labs(x = "Year",y="Shannon Index")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()

# Shannon

H.VOL <- ggplot(variables.all, aes(x=year, y=H.VOL))+
  geom_line() +
  facet_wrap(~case, ncol=2)+
  ggtitle("Shannon Index on volume proportion per species")+
  labs(x = "Year",y="Shannon Index")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()

# Shannon

H.BA <- ggplot(variables.all, aes(x=year, y=H.BA))+
  geom_line() +
  facet_wrap(~case, ncol=2)+
  ggtitle("Shannon Index on basal area per species")+
  labs(x = "Year",y="Shannon Index")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()

#-------------------------------------------------------------------------------
# Make new plots


plot_long <- plot_variables_all %>%
  pivot_longer(
    cols = c(lai_sim, ba_broadl, trees_10_40, broadl_40, age),
    names_to = "variable",
    values_to = "value"
  )

labs <- c(
  lai_sim = "LAI",
  ba_broadl = "Broadleaf BA (m²/ha)",
  trees_10_40 = "Trees (10–40 cm)",
  broadl_40 = "Broadleaf >40 cm",
  age = "Stand age"
)

plot_long$variable <- factor(
  plot_long$variable,
  levels = names(labs),
  labels = labs
)

ggplot(plot_long,
       aes(year, value, colour = run, group = run)) +
  geom_line(linewidth = 1) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  )


########################################################## CLOSE EVERY PLOTs

dev.off()
