library(tidyverse)
library(FCSplankton)

setwd("PATH/TO/PROJECT/")

project <- basename(getwd())
PSD_all <- read_csv(paste0("Influx_", project,"_PSD.csv"))
PSD_all[1:3,]

# Note for conversion from carbon per cell to equivalent spherical diameter
# Menden-Deuer, S. & Lessard, E. J. Carbon to volume relationships for dinoflagellates, diatoms, and other protist plankton. Limnol. Oceanogr. 45, 569–579 (2000).
d <- 0.261; e <- 0.860 # < 3000 µm3

# Calculate cell abundance, biomass, carbon quota, and cell diameter
influx_replicates <- PSD_all %>%
  group_by(time,station,cast,lat,lon,depth,pop,replicate) %>%
  dplyr::summarise(
    lat = unique(lat, na.rm = TRUE),
    lon = unique(lon, na.rm = TRUE),
    cell_abundance = sum(abundance_per_bin, na.rm = TRUE), # calculate cell abundance per population
    carbon_biomass = sum(biomass_per_bin, na.rm = TRUE)) %>% # calculate carbon biomass per population
  mutate(carbon_quota = carbon_biomass / cell_abundance, # calculate carbon per cell
         diameter = round(2*(3/(4*base::pi)*(carbon_quota/d)^(1/e))^(1/3),5)) %>% # convert carbon per cell to equivalent spherical diameter
  arrange(time)

# Average data of sample replicates
influx <- influx_replicates %>%
  group_by(time,station,cast,lat,lon,depth,pop) %>%
  dplyr::summarise(
    lat = mean(lat, na.rm = TRUE),
    lon = mean(lon, na.rm = TRUE),
    cell_abundance = mean(cell_abundance, na.rm = TRUE), # calculate cell abundance per population
    carbon_biomass = mean(carbon_biomass, na.rm = TRUE), # calculate carbon biomass per population
    carbon_quota = mean(carbon_quota), # calculate carbon per cell
    diameter = mean(diameter, na.rm = TRUE)) %>% # convert carbon per cell to equivalent spherical diameter
  arrange(time)

# Make sure time is in the proper format for CMAP
influx$time <- str_replace_all(as.character(influx$time), " ","T")

# Save the data in an excel spreadsheet following CMAP data submission guidlines
project <- basename(getwd())
cruise <- "" # Cruise ID (ex. KM1906); leave blank if samples were not collected during a cruise
cruise_keywords <- "" # Cruise keywords, alternate cruise names commonly referred to (ex. "Gradients 2, Gradients 2017, NPSG, Thomas G. Thompson"), any relevant geographical, personal, dataset specific keywords for querying data

cmap_convert(data = influx , cruise, cruise_keywords, project, version = "v1.0")
