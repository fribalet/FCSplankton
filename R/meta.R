#' Format Influx data into CAMP compatible format.
#'
#' @param data data from Influx.
#' @param cruise cruise name, if any (otherwise NA).
#' @param cruise_nickname cruise nickname name, if any (otherwise NA).
#' @param project cruise name, if any (otherwise NA).
#' @param unstained Sample stain status (TRUE/FALSE).
#' @param version Version of the dataset.
#' @return None
#' @examples
#' \dontrun{
#' csv_convert(db, meta, path)
#' }
#' @export
cmap_convert<- function(data, cruise, cruise_nickname, project, unstained, version = "v1.0") {

  ## FORMAT DATA

  ## Assign specific refractive indexes for each population for Mie theory conversion
  # Make dataframe with select populations that have a high index of refraction (low estimates)
  lwr <- data %>%
        dplyr::filter(population == "picoeuk" | population == "large-picoeuk" | population == "small-picoeuk" | population == "prochloro" | population == "synecho" | population == "bacteria") %>%
        dplyr::select(-diam_mid, -diam_upr, -Qc_mid, -Qc_upr) %>%
        dplyr::rename(cell_diameter = diam_lwr, carbon_content = Qc_lwr)

  # Make dataframe with select populations that have a mid index of refraction (middle estimates)
  mid <- data %>%
        dplyr::filter(population == "unknown" | population == "beads" | population == "croco") %>%
        dplyr::select(-diam_lwr, -diam_upr, -Qc_lwr, -Qc_upr) %>%
        dplyr::rename(cell_diameter = diam_mid, carbon_content = Qc_mid)

  # merge dataframes
  data <- merge(lwr, mid, all = TRUE)

  data$biomass <- data$carbon_content * data$abundance

  # Set core keywords
  core <- paste("discrete flow cytometry, BD Influx cell sorter, insitu, in-situ, biology, phytoplankton, picophytoplankton, Armbrust, UW, University of Washington", cruise, cruise_nickname, sep = ", ")

  # Pivot data from long format to wide; split population data for each variable (value) column
  data.pivot <- data %>%
        pivot_wider(names_from = population,
                    values_from = c(count, scatter, red, orange, green, abundance, cell_diameter, carbon_content, biomass))
if(unstained == TRUE){
    data.pivot <- select(data.pivot, !starts_with("green"))
  }

  ## Format population data
  # master list of possible population names; add to as needed
  long_pop <- c("prochlorococcus", "synechococcus", "picoeukaryote", "large-picoeukaryote", "small-picoeukaryote", "crocosphaera", "bacteria", "unknown", "beads")

  # get unique population names from dataset
  data_pop <- str_sort(unique(data$population))

  # anchor population name to subset correct population names based on matching the beginning of population names
  match_pop <- str_replace_all(data_pop, "^", "^")

  # loop through list of populations from data to get list of full population names
  pop <- c()
  for(i in 1:length(data_pop)){
    pop[[i]] <- str_subset(long_pop, match_pop[i])
  }

  ## FORMAT METADATA
  #set order of standard variable metadata
  var_cols <- colnames(data.pivot)

  var_data <- c("file" ,
            str_sort(str_subset(var_cols, "count_")),
            str_sort(str_subset(var_cols, "scatter_")),
            str_sort(str_subset(var_cols, "red_")),
            str_sort(str_subset(var_cols, "orange_")),
            str_sort(str_subset(var_cols, "green_")),
            "volume",
            str_sort(str_subset(var_cols, "abundance_")),
            str_sort(str_subset(var_cols, "cell_diameter_")),
            #str_sort(str_subset(var_cols, "diam_lwr_")),
            #str_sort(str_subset(var_cols, "diam_mid_")),
            #str_sort(str_subset(var_cols, "diam_upr_")),
            str_sort(str_subset(var_cols, "carbon_content_")),
            #str_sort(str_subset(var_cols, "Qc_lwr_")),
            #str_sort(str_subset(var_cols, "Qc_mid_")),
            #str_sort(str_subset(var_cols, "Qc_upr_")),
            str_sort(str_subset(var_cols, "biomass_")),
            "flag",
            "stain")

  # list of each variable (value) that got pivotted with each population for required CMAP variable metadata
  long_count <- c()
  long_scatter <- c()
  long_red <- c()
  long_orange <- c()
  long_green <- c()
  long_abundance <- c()
  long_cell_diameter <- c()
  #long_diameter_lwr <- c()
  #long_diameter_mid <- c()
  #long_diameter_upr <- c()
  long_carbon_content <- c()
  #long_carbon_lwr <- c()
  #long_carbon_mid <- c()
  #long_carbon_upr <- c()
  long_biomass <- c()
  comment_count <- c()
  comment_scatter <- c()
  comment_red <- c()
  comment_orange <- c()
  comment_green <- c()
  comment_abundance <- c()
  comment_cell_diameter <- c()
  #comment_diameter_lwr <- c()
  #comment_diameter_mid <- c()
  #comment_diameter_upr <- c()
  comment_carbon_content <- c()
  #comment_carbon_lwr <- c()
  #comment_carbon_mid <- c()
  #comment_carbon_upr <- c()
  comment_biomass <- c()
  keyword_count <- c()
  keyword_scatter <- c()
  keyword_red <- c()
  keyword_orange <- c()
  keyword_green <- c()
  keyword_abundance <- c()
  keyword_cell_diameter <- c()
  #keyword_diameter_lwr <- c()
  #keyword_diameter_mid <- c()
  #keyword_diameter_upr <- c()
  keyword_carbon_content <- c()
  #keyword_carbon_lwr <- c()
  #keyword_carbon_mid <- c()
  #keyword_carbon_upr <- c()
  keyword_biomass <- c()

  # run through loop to add all the population dependent variable metadata
  for(i in 1:length(pop)){
    if(pop[i] == "picoeukaryote" | pop[i] == "prochlorococcus" | pop[i] == "synechococcus" | pop[i] == "bacteria"){
      RI <- "high refractive index"
      IR <- "1.41"
      }else{
      RI <- "mid refractive index"
      IR <- "1.38"
      }

    long_count[[i]] <- paste0("number of ", pop[i],"-like particles counted by the instrument")
    long_scatter[[i]] <- paste0("forward angle light scatter of ", pop[i], "-like particles, normalized to 1 micron beads (proxy of cell diameter)")
    long_red[[i]] <- paste0("red fluorescence of ", pop[i], "-like particles, normalized to 1 micron beads (proxy of chlorophyll content)")
    long_orange[[i]] <- paste0("orange fluorescence of ", pop[i], "-like particles normalized to 1 micron beads (proxy of phycoerythrin content)")
    long_green[[i]] <- paste0("green fluorescence of ", pop[i], "-like particles normalized to 1 micron beads (nucleic acid stain)")
    long_abundance[[i]] <- paste0("cell abundance of ", pop[i], "-like particles")
    long_cell_diameter[[i]] <- paste0("equivalent spherical diameter of ", pop[i], "-like particles using ", RI)
    #long_diameter_lwr[[i]] <- paste(pop[i], "equivalent spherical diameter using high refractive index")
    #long_diameter_mid[[i]] <- paste(pop[i], "equivalent spherical diameter using mid refractive index")
    #long_diameter_upr[[i]] <- paste(pop[i], "equivalent spherical diameter using low refractive index")
    long_carbon_content[[i]] <- paste0("cellular carbon content of ", pop[i], "-like particles using ", RI)
    #long_carbon_lwr[[i]] <- paste(pop[i], "cellular carbon content using high refractive index")
    #long_carbon_mid[[i]] <- paste(pop[i], "cellular carbon content using mid refractive index")
    #long_carbon_upr[[i]] <- paste(pop[i], "cellular carbon content using low refractive index")
    long_biomass[[i]] <- paste0("Carbon biomass of ", pop[i], "population")
    comment_count[[i]] <- paste(pop[i], "count needs to be > 30 to be trusted")
    comment_scatter[[i]] <- paste(pop[i], "light scatter collected using a 457-50 bandpass filter")
    comment_red[[i]] <- paste(pop[i], "red fluorescence collected using a 692-40 bandpass filter")
    comment_orange[[i]] <- paste(pop[i], "orange fluorescence collected using a 572-27 bandpass filter")
    comment_green[[i]] <- paste(pop[i], "green fluorescence collected using a 530-40 bandpass filter")
    comment_abundance[[i]] <- paste(pop[i], "cell abundance, number of particles divided by the volume of sample, see https://github.com/fribalet/FCSplankton for more details")
    comment_cell_diameter[[i]] <- paste(pop[i], "cell diameter based on Mie theory using an index of refraction of", IR, "for phytoplankton and 1.337 for seawater, see https://github.com/seaflow-uw/fsc-size-calibration for more details")
    #comment_diameter_lwr[[i]] <- paste(pop[i], "cell diameter based on Mie theory using an index of refraction of 1.41 for phytoplankton and 1.337 for seawater, see https://github.com/seaflow-uw/fsc-size-calibration for more details")
    #comment_diameter_mid[[i]] <- paste(pop[i], "cell diameter based on Mie theory using an index of refraction of 1.38 for phytoplankton and 1.337 for seawater, see https://github.com/seaflow-uw/fsc-size-calibration for more details")
    #comment_diameter_upr[[i]] <- paste(pop[i], "cell diameter based on Mie theory using an index of refraction of 1.36 for phytoplankton and 1.337 for seawater, see https://github.com/seaflow-uw/fsc-size-calibration for more details")
    comment_carbon_content[[i]] <- paste(pop[i], "carbon content based on the equation fgC cell-1 = 0.261 x Volume^0.860, where Volume is calculated from cell diameter (Mie-based, using refractive index of", IR, "for phytoplankton) assuming spherical particle; see https://github.com/seaflow-uw/fsc-poc-calibration for more details")
    #comment_carbon_lwr[[i]] <- paste(pop[i], "carbon content based on the equation fgC cell-1 = 0.261 x Volume^0.860, where Volume is calculated from cell diameter (Mie-based, using refractive index of 1.41 for phytoplankton) assuming spherical particle; see https://github.com/seaflow-uw/fsc-poc-calibration for more details")
    #comment_carbon_mid[[i]] <- paste(pop[i], "carbon content based on the equation fgC cell-1 = 0.261 x Volume^0.860, where Volume is calculated from cell diameter (Mie-based, using refractive index of 1.38 for phytoplankton) assuming spherical particle; see https://github.com/seaflow-uw/fsc-poc-calibration for more details")
    #comment_carbon_upr[[i]] <-paste(pop[i], "carbon content based on the equation fgC cell-1 = 0.261 x Volume^0.860, where Volume is calculated from cell diameter (Mie-based, using refractive index of 1.36 for phytoplankton) assuming spherical particle; see https://github.com/seaflow-uw/fsc-poc-calibration for more details")
    comment_biomass[[i]] <- paste(pop[i], "carbon biomass = cell abundance x carbon content")
    keyword_count[[i]] <- paste(pop[i], "particle count,", core)
    keyword_scatter[[i]] <- paste(pop[i], "forward angle light scatter, FSC, FALS,", core)
    keyword_red[[i]] <- paste(pop[i], "red fluorescence, chlorophyll,", core)
    keyword_orange[[i]] <- paste(pop[i], "orange fluorescence, phycoerythrin,", core)
    keyword_green[[i]] <- paste(pop[i], "green fluorescence, SYBR Green I, nucleic acid stain,", core)
    keyword_abundance[[i]] <- paste(pop[i],"abundance, cell concentration, cell count, cell abundance,", core)
    keyword_cell_diameter[[i]] <- paste(pop[i],"size, diameter, ESD,", core)
    #keyword_diameter_lwr[[i]] <- paste(pop[i],"size, diameter, ESD,", core)
    #keyword_diameter_mid[[i]] <- paste(pop[i],"size, diameter, ESD,", core)
    #keyword_diameter_upr[[i]] <- paste(pop[i],"size, diameter, ESD,", core)
    keyword_carbon_content[[i]] <- paste(pop[i], "quotas, carbon, biomass, POC,", core)
    #keyword_carbon_lwr[[i]] <- paste(pop[i], "quotas, carbon, biomass, POC,", core)
    #keyword_carbon_mid[[i]] <- paste(pop[i], "quotas, carbon, biomass, POC,", core)
    #keyword_carbon_upr[[i]] <- paste(pop[i], "quotas, carbon, biomass, POC,", core)
    keyword_biomass[[i]] <- paste(pop[i], "carbon biomass, POC,", core)
  }

  var_long_name <- c("flow cytometry standard file (.fcs)",
                     long_count,
                     long_scatter,
                     long_red,
                     long_orange,
                     long_green,
                     "volume measured by pressure sensor",
                     long_abundance,
                     long_cell_diameter,
                     #long_diameter_lwr,
                     #long_diameter_mid,
                     #long_diameter_upr,
                     long_carbon_content,
                     #long_carbon_lwr,
                     #long_carbon_mid,
                     #long_carbon_upr,
                     long_biomass,
                     "sample status flag",
                     "sample stain flag")

  var_comment <- c("",
                   comment_count,
                   comment_scatter,
                   comment_red,
                   comment_orange,
                   comment_green,
                   "volume measured by pressure sensor",
                   comment_abundance,
                   comment_cell_diameter,
                   #comment_diameter_lwr,
                   #comment_diameter_mid,
                   #comment_diameter_upr,
                   comment_carbon_content,
                   #comment_carbon_lwr,
                   #comment_carbon_mid,
                   #comment_carbon_upr,
                   comment_biomass,
                   "outliers (0 = quality data; 1 = issue related to instrument performance; 2 = issue related to population classification)",
                   "DNA stain (0 = unstained sample; 1 = stained sample)")

  var_unit <- c("",
                rep("",length(pop)*5),
                "microliter",
                rep("cells per microliter",length(pop)),
                rep("micrometer", length(pop)),
                rep("picogram carbon per cell", length(pop)),
                rep("microgram carbon per liter", length(pop)),
                "",
                "")

  var_sensor <- rep("Flow cytometry",length(var_data))

  var_discipline <- c("",
                     rep("Biology",length(pop)*5),
                     "",
                     rep("Biology",length(pop)),
                     rep("Biology, Biogeochemistry", length(pop)*3),
                     "",
                     "")

  visualize <- c(0, rep(1,length(pop)*5), 0, rep(1,length(pop)*4), 0, 0)

  var_keywords <- c(paste("file,", core),
        keyword_count,
        keyword_scatter,
        keyword_red,
        keyword_orange,
        keyword_green,
        paste("volume, vol,", core),
        keyword_abundance,
        keyword_cell_diameter,
        #keyword_diameter_lwr,
        #keyword_diameter_mid,
        #keyword_diameter_upr,
        keyword_carbon_content,
        #keyword_carbon_lwr,
        #keyword_carbon_mid,
        #keyword_carbon_upr,
        keyword_biomass,
        paste("flag,", core),
        paste("DNA stain, SYBR stain,", core))

  if(unstained == TRUE){
    var_data <- var_data[lapply(var_data,function(x) length(grep("green",x,value=FALSE))) == 0]

    var_long_name <- var_long_name[lapply(var_long_name,function(x) length(grep("green",x,value=FALSE))) == 0]

    var_comment <- var_comment[lapply(var_comment,function(x) length(grep("green",x,value=FALSE))) == 0]

    var_unit <- c("",
                  rep("",length(pop)*4),
                  "microliter",
                  rep("cells per microliter",length(pop)),
                  rep("micrometer", length(pop)),
                  rep("picogram carbon per cell", length(pop)),
                  rep("microgram carbon per liter", length(pop)),
                  "",
                  "")

    var_sensor <- rep("Flow cytometry",length(var_data))

    var_discipline <- c("",
                       rep("Biology",length(pop)*4),
                       "",
                       rep("Biology",length(pop)),
                       rep("Biology, Biogeochemistry", length(pop)*3),
                       "",
                       "")

    visualize <- c(0, rep(1,length(pop)*4), 0, rep(1,length(pop)*4), 0, 0)

    var_keywords <- var_keywords[lapply(var_keywords,function(x) length(grep("green",x,value=FALSE))) == 0]
  }


  # add custom column to metadata
  coordinates <- c("time","lat","lon","depth")
  id0 <- match(coordinates,colnames(data.pivot)) # which column are standard
  id1 <- match(var_data,colnames(data.pivot)) # which column are standard
  id2 <- which(is.na(match(1:ncol(data.pivot), c(id0, id1)))) # which column are not standard
  id <- c(id0, id1, id2)

  # vars_metadata
  vars_metadata <- dplyr::tibble(
                          var_short_name = var_data,
                          var_long_name,
                          var_sensor,
                          var_unit,
                          var_spatial_res = "irregular",
                          var_temporal_res = "irregular",
                          var_discipline,
                          visualize,
                          var_keywords,
                          var_comment)

  # custom metadata
  custom_metadata <- dplyr::tibble(
                          var_short_name = colnames(data.pivot)[id2],
                          var_long_name = "",
                          var_sensor = "",
                          var_unit = "",
                          var_spatial_res = "irregular",
                          var_temporal_res = "irregular",
                          var_discipline = "",
                          visualize = rep(0, length(id2)),
                          var_keywords = core,
                          var_comment = "")

  allvars_metadata <- rbind(vars_metadata, custom_metadata)

  # dataset_metadata
  dataset_metadata <- dplyr::tibble(
                        dataset_short_name = paste0("Influx_",project),
                        dataset_long_name = paste0("Discrete flow cytometry of ",str_replace_all(project, "_", " "), " using a BD Influx Cell Sorter"),
                        dataset_version = version,
                        dataset_release_date = as.Date(Sys.time()),
                        dataset_make = "observation",
                        dataset_source = "Armbrust Lab, University of Washington",
                        dataset_acknowledgement = "",
                        dataset_description = "to be added by data owner",
                        dataset_references = "",
                        climatology = NULL,
                        cruise_names = cruise)

  # reorder order of column so they match metadata
  data.pivot <- data.pivot[,id]

  ## Save data
  openxlsx::write.xlsx(x=list(data.pivot, dataset_metadata, allvars_metadata),
                      file = paste0("Influx_", project, "_",as.Date(Sys.time()),"_" ,version,".xlsx"),
                      sheetName=c('data','dataset_meta_data','vars_meta_data'))

}
