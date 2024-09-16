library(tidyverse)
library(FCSplankton)

setwd("PATH/TO/PROJECT/")

## Read in PSD data
project <- basename(getwd())
PSD_all <- read_csv(paste0("Influx_", project,"_PSD.csv"))
PSD_all[1:3,]

## Summarize PSD across sample replicates
PSD_all <- PSD_all %>%
  dplyr::select(-sample,-time,-file,-stain,-flag,-replicate,-comments,-volume) %>% # remove columns that would inhibit averaging across replicates
  dplyr::group_by(station,depth,Qc,pop) %>%
  summarize_all(function(x) mean(x, na.rm=T)) %>%
  arrange(station)

## Create station position labels
PSD_all$lat <- trunc(PSD_all$lat *10^2)/10^2
PSD_all$lon <- trunc(PSD_all$lon*10^2)/10^2
PSD_all$position <- paste0("St. ",PSD_all$station," (",PSD_all$lat,",",PSD_all$lon,")")
PSD_all$position <- as_factor(PSD_all$position)

## Set population group colors
group_colors <- c(bacteria = "lightsalmon1",
                  prochloro=viridis::viridis(4)[1],
                  synecho=viridis::viridis(4)[2],
                  picoeuk=viridis::viridis(4)[3])

## Plot all pop distributions with height allow biomass comparisons across population
PSD_all %>%
  group_by(position,depth) %>%
  ggplot() +
  ggridges::geom_density_ridges(aes(x = Qc, y = -depth, height = biomass_per_bin, fill =pop, group = interaction(depth,pop)), stat="identity", color="darkgrey", alpha=0.55,size=.15,panel_scaling=FALSE) + # panel_scaling=TRUE relative scaling is calculated separately for each panel. panel_scaling=FALSE, relative scaling is calculated globally
  scale_x_continuous(trans = "log10") +
  scale_fill_manual(name = 'Population', values = group_colors, breaks = c("bacteria",'prochloro','synecho',"picoeuk"), labels = c("Bacteria",'Pro','Syn',"Picoeuk")) +
  theme(legend.key.size = unit(.35, 'cm')) +
  annotation_logticks(sides = "b")  +
  theme_bw() +
  facet_wrap( . ~ position,ncol=4) +
  labs(x="Carbon Content Distribution (pgC)",
       y= "Depth (m)")
ggsave("biomass_distribution.png", path = "./plots",height=10,width=8)

## Plot all pop distributions with height allow abundance comparisons across population
PSD_all %>%
  group_by(position,depth) %>%
  ggplot() +
  ggridges::geom_density_ridges(aes(x = Qc, y = -depth, height = abundance_per_bin, fill =pop, group = interaction(depth,pop)), stat="identity", color="darkgrey", alpha=0.55,size=.15,panel_scaling=FALSE) + # panel_scaling=TRUE relative scaling is calculated separately for each panel. panel_scaling=FALSE, relative scaling is calculated globally
  scale_x_continuous(trans = "log10") +
  scale_fill_manual(name = 'Population', values = group_colors, breaks = c("bacteria",'prochloro','synecho',"picoeuk"), labels = c("Bacteria",'Pro','Syn',"Picoeuk")) +
  theme(legend.key.size = unit(.35, 'cm')) +
  annotation_logticks(sides = "b")  +
  theme_bw() +
  facet_wrap( . ~ position,ncol=4) +
  labs(x="Carbon Content Distribution (pgC)",
       y= "Depth (m)")
ggsave("abundance_distribution.png", path = "./plots",height=10,width=8)

## Create individual PSDs for each population (distribution heights not comparable across populations)
PSD_wide <- PSD_all %>%
  pivot_wider(names_from = pop, values_from = c("abundance_per_bin","biomass_per_bin"),values_fill=0) #issue is with cell_diameter column
PSD_wide[1:3,]

## Plot individual population PSDs (distribution heights not comparable across populations)
# Plot only to see the the distributions of each population, not quantitive or comparable across populations and panels
PSD_wide %>%
  group_by(position,depth) %>%
  ggplot() +
  ggridges::geom_density_ridges2(aes(x = Qc, y = -depth, height = abundance_per_bin_bacteria, fill = "bacteria", group = depth), stat="identity",  color="darkgrey", alpha=0.4,size=.25,panel_scaling=TRUE) + # panel_scaling=TRUE relative scaling is calculated separately for each panel. panel_scaling=FALSE, relative scaling is calculated globally
  ggridges::geom_density_ridges2(aes(x = Qc, y = -depth, height = abundance_per_bin_prochloro, fill ="prochloro", group = depth), stat="identity", color="darkgrey", alpha=0.55,size=.25,panel_scaling=TRUE) + # panel_scaling=TRUE relative scaling is calculated separately for each panel. panel_scaling=FALSE, relative scaling is calculated globally
  ggridges::geom_density_ridges2(aes(x = Qc, y = -depth, height = abundance_per_bin_synecho, fill="synecho", group = depth), stat="identity", color="darkgrey", alpha=0.55,size=.25,panel_scaling=TRUE) + # panel_scaling=TRUE relative scaling is calculated separately for each panel. panel_scaling=FALSE, relative scaling is calculated globally
  ggridges::geom_density_ridges2(aes(x = Qc, y = -depth, height = abundance_per_bin_picoeuk, fill='picoeuk', group = depth), stat="identity",  color="darkgrey", alpha=0.55,size=.25,panel_scaling=TRUE) + # panel_scaling=TRUE relative scaling is calculated separately for each panel. panel_scaling=FALSE, relative scaling is calculated globally
    scale_x_continuous(trans = "log10",breaks = c(.01,.1,1,1,10), minor_breaks = c(.005,0.05,.5,5)) +
  scale_fill_manual(name = 'Population', values = group_colors, breaks = c("bacteria",'prochloro','synecho',"picoeuk"), labels = c("Bacteria",'Pro','Syn',"Picoeuk")) +
  theme(legend.key.size = unit(.35, 'cm')) +
  annotation_logticks(sides = "b")  +
  theme_bw() +
  facet_wrap( . ~ position,ncol=4) +
  labs(x="Carbon Content Distribution (pgC)",
       y= "Depth (m)")
ggsave("individual_abundance_distribution.png", path = "./plots",height=10,width=8)
