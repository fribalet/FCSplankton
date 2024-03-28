library(tidyverse)
library(FCSplankton)

setwd("PATH/TO/PROJECT/")

unstained <- FALSE

### Load Data

## Load mie lookup table
mie <- as.data.frame(read.csv(system.file("scatter", paste0("calibrated-mieINFLUX.csv"),package="FCSplankton")))
mie[1:3,] ## NOTE: Leo and Penny are included in the same Mie lookup table.

## Load metadata
meta <- read_csv("metadata.txt")
meta[1:3,]

## Convert metadata
file <- paste0(meta$file,".fcs") # format  sample name to filename (.fcs)
time <- meta$time
sample <- meta$sample
station <- meta$station
cast <- meta$cast
lat <- meta$lat
lon <- meta$lon
depth <- meta$depth
replicate <- meta$replicate
volume <- meta$volume
stain <- meta$stain
flag <- meta$flag
comments <- meta$comments

# create new metadata tibble
metadata <- tibble(file, time, sample, station, cast, lat, lon, depth, replicate, volume, stain, flag, comments)
metadata[1:3,]


### Create and bin particle size distribution (PSD) data

## Set PSD bin widths
m <- 137
v_min <- 0.002
delta_v_inv <- 11 # define the width of each bin
breaks <- round(v_min*2^(((1:m)+1)*delta_v_inv^-1),4) # log 2 spaced binning
print(breaks)

## Load gated FCS data to generate PSD
if(unstained){folder <- "unstained"
}else{folder <- c("unstained","stained")}

distribution <- tibble()

for(folder_name in folder){
  # Read FCS
  file_list <- list.files(paste0(folder_name,"/raw"), pattern = ".fcs$", full.names=T)

  # load one file
  this_file <- file_list[1]
  fcs <- read_influx(this_file, transformation = TRUE) # create a dataframe WITH log-amplified data)

  # replace number by column indice of FSC, 692, 580 and 530 respectively
  id <- c(13,23,19,17)
  names(fcs)[id]
  names_pmt <- names(fcs)[id] # original names of FSC, 692, 580 and 530 respectively

  ## PSD
  for (this_file in file_list){

    # Read data
    print(paste(this_file))
    fcs <- read_influx(this_file, transformation = TRUE) # create a flowFrame WITH transformation (for log-amplified data)

    # change header
    id <- match(names_pmt, names(fcs))
    names(fcs)[id] <- c("scatter", "red", "orange", "green")

    # Load gating parameters
    previous <- sub("raw","gating",paste0(this_file, ".RData"))
    if(file.exists(previous)) load(previous)

    # Assign "beads" label to particles
    fcs <- classify_fcs(fcs, gates.log[][1])

    # Normalize to beads
    beads <- fcs[which(fcs$pop == "beads"),]
    fcs$norm.scatter <- fcs$scatter / mean(beads$scatter)
    fcs$norm.orange <- fcs$orange / mean(beads$orange)
    fcs$norm.red <- fcs$red / mean(beads$red)
    fcs$norm.green <- fcs$green / mean(beads$green)

    # Apply gates and label particles according to 'gates.log'
    fcs <- try(classify_fcs(fcs, gates.log))
    if(class(fcs) == "try-error") next

    # Convert scatter to Qc
    id <- findInterval(fcs$norm.scatter, mie$scatter)
    id[which(id == 0)] <- 1
    fcs$quotas <- mie[id, "Qc_Penny_lwr"]

    # Select wanted and remove unwanted populations
    if(folder == "unstained") fcs.pop <- dplyr::filter(fcs, pop != "unknown" & pop != "beads")
    if(folder == "stained") fcs.pop <- dplyr::filter(fcs, pop== "bacteria")

    # create PSD
    psd <- fcs.pop %>%
      group_by(pop, breaks = cut(quotas, breaks = c(breaks,Inf), right=FALSE), .drop=F) %>%
      count(breaks) %>%
      pivot_wider(names_from = breaks, values_from = n)
    psd <- psd %>% add_column(file = basename(this_file), .before = 1)
    distribution <- dplyr::bind_rows(distribution, psd)
  }

  # Merge sample metadata with PSD
  DIST <- as_tibble(merge(distribution, metadata, by="file", all.x=TRUE))

  # convert breaks to geometric mean values
  clmn <- grep("^\\[.*?)$", names(DIST))

  # calculate abundance in each size class bin
  DIST[,clmn] <-  DIST[,clmn] / DIST[["volume"]]


}


### Data Correction

# Get PSD bin column names
clmn <- grep("^\\[.*?)$", names(DIST))

## Remove pro population from bacteria + prochlorococcus gate PSD
PSD_all <- cyano_psd <- DIST %>%
  dplyr::filter(stain!=1) # Only include unstained populations. Bacteria samples will be added back during the correction step.

pro <- DIST%>%
  dplyr::filter(pop == "prochloro" & stain!=1)

# Bacteria population correction
if(unstained == FALSE){
  bact <- hetero <- DIST %>%
    dplyr::filter(pop == "bacteria")

  # find matching sample numbers
  id.pro <- match(bact$sample, pro$sample, nomatch=0)
  id.bact <- which(id.pro != 0)

  # sanity check
  bact$sample[id.bact] == pro$sample[id.pro]

  # Remove Pro population from bacteria PSD - NOTE: this Assumes staining does not impact scatter values
  clmn <- grep("^\\[.*?)$", names(bact))
  hetero[id.bact,clmn] <- bact[id.bact,clmn] - pro[id.pro, clmn]
  hetero <- hetero[id.bact, ]
  hetero[,clmn][hetero[,clmn] < 0] <- 0
  hetero <- dplyr::select(hetero,-sample)
  hetero[1:3,]

  PSD_all <- rbind(cyano_psd, hetero) # Combine unstained PSD with corrected bacteria PSD
}

PSD_all[1:3,]

# calculate the geometric mean Qc for each size class
clmn <- grep(")", names(PSD_all))
b <- strsplit(sub("\\[","",sub("\\)","",colnames(PSD_all)[clmn])),",")
Qc_geom_mean <- unlist(list(lapply(b, function(x) sqrt(mean(as.numeric(x))*max(as.numeric(x))))))

# convert dataframe to long format
PSD_long <- PSD_all %>%
  pivot_longer(
    cols = (clmn),
    names_to = "Qc",
    values_to = "abundance_per_bin",
    values_drop_na = TRUE
  )

# Remove infinite values from PSD bin breaks
PSD_long$Qc <- as.numeric(PSD_long$Qc)
PSD_long$Qc[sapply(PSD_long$Qc, is.infinite)] <- NA
PSD <- dplyr::filter(PSD_long,abundance_per_bin != 0)

# Calculate biomass (µgC/L) per bin
PSD$biomass_per_bin <- PSD$Qc * PSD$abundance_per_bin # (pgC/cell * cell/µL)

# Find closest equivilent spherical diameter matches to the carbon quota from the Mie lookup table for each particle
mie <- as.data.frame(read.csv(system.file("scatter", paste0("calibrated-mieINFLUX.csv"),package="FCSplankton")))
id <- findInterval(PSD$Qc, mie$"Qc_Penny_lwr", all.inside = TRUE)
PSD[,"esd"] <- mie[id,"diam_Penny_lwr"]
PSD[1:3,]

###  Save PSD data in a .csv file
PSD[1:3,]
project <- basename(getwd())
write_csv(PSD,file = paste0("Influx_", project,"_PSD.csv"))
