#' Format Influx data into CAMP compatible format.
#'
#' @param data data from Influx.
#' @param cruise cruise name, if any (otherwise NA).
#' @param cruise_keywords cruise keywords, if any (otherwise NA).
#' @param project cruise name, if any (otherwise NA).
#' @param version Version of the dataset.
#' @return None
#' @examples
#' \dontrun{
#' csv_convert(db, meta, path)
#' }
#' @export
cmap_convert <- function(data, cruise, cruise_keywords, project, version = "v1.0") {

  ## FORMAT DATA

  # Set core keywords
  core <- paste("BD Influx Cell Sorter, biogeochemistry, biology, cruise, cyanobacteria, discrete flow cytometry, FACS, FCM, in situ, insitu, in-situ, observation, phytoplankton, picophytoplankton, Armbrust, UW, University of Washington, Kelsy Cain, François Ribalet", cruise, cruise_keywords, sep = ", ")

  # Pivot data from long format to wide; split population data for each variable (value) column
  data.pivot <- data %>%
    pivot_wider(names_from = pop,
                values_from = c(cell_abundance, diameter, carbon_quota, carbon_biomass))
  # if(unstained == TRUE){data.pivot <- select(data.pivot, !starts_with("green"))}

  ## Format population data
  # master list of possible population names; add to as needed
  long_pop <- c("Prochlorococcus", "Synechococcus", "Picoeukaryote", "Large-Picoeukaryote", "Small-Picoeukaryote", "Crocosphaera", "Bacteria", "Unknown", "Beads")

  # get unique population names from dataset
  data_pop <- str_sort(unique(data$pop))
  str_sub(data_pop, 1, 1) <- str_sub(data_pop, 1, 1) %>% str_to_upper()

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

  var_data <- c(str_sort(str_subset(var_cols, "cell_abundance_")),
                str_sort(str_subset(var_cols, "diameter_")),
                str_sort(str_subset(var_cols, "carbon_quota_")),
                str_sort(str_subset(var_cols, "carbon_biomass_")))

  # list of each variable (value) that got pivotted with each population for required CMAP variable metadata
  long_cell_abundance <- c()
  long_diameter <- c()
  long_carbon_quota <- c()
  long_carbon_biomass <- c()
  comment_cell_abundance <- c()
  comment_diameter <- c()
  comment_carbon_quota <- c()
  comment_carbon_biomass <- c()
  keyword_cell_abundance <- c()
  keyword_diameter <- c()
  keyword_carbon_quota <- c()
  keyword_carbon_biomass <- c()

  # run through loop to add all the population dependent variable metadata
  for(i in 1:length(pop)){
    if(pop[i] == "picoeukaryote" | pop[i] == "prochlorococcus" | pop[i] == "synechococcus" | pop[i] == "bacteria"){
      RI <- "high refractive index"
      IR <- "1.41"
    }else{
      RI <- "mid refractive index"
      IR <- "1.38"
    }

    long_cell_abundance[[i]] <- paste0("abundance of ", pop[i], "-like particles")
    long_diameter[[i]] <- paste0("diameter of ", pop[i], "-like particles")
    long_carbon_quota[[i]] <- paste0("carbon quota of ", pop[i], "-like particles")
    long_carbon_biomass[[i]] <- paste0("carbon biomass of ", pop[i], " population")
    comment_cell_abundance[[i]] <- paste(pop[i], "cell abundance, number of particles divided by the volume of sample, see https://github.com/fribalet/FCSplankton for more details")
    comment_diameter[[i]] <- paste(pop[i], "cell diameter based on Mie theory using an index of refraction of", IR, "for phytoplankton and 1.337 for seawater, see https://github.com/seaflow-uw/fsc-size-calibration for more details")
    comment_carbon_quota[[i]] <- paste(pop[i], "carbon content based on the equation pgC cell-1 = 0.261 x Volume^0.860, where Volume is calculated from cell diameter (Mie-based, using refractive index of", IR, "for phytoplankton) assuming spherical particle; see https://github.com/seaflow-uw/fsc-poc-calibration for more details")
    comment_carbon_biomass[[i]] <- paste(pop[i], "carbon biomass = cell abundance x carbon content")
    keyword_cell_abundance[[i]] <- paste(pop[i],"abundance, cell concentration, cell abundance,", core)
    keyword_diameter[[i]] <- paste(pop[i],"cell size, diameter, equivalent spherical diameter, ESD,", core)
    keyword_carbon_quota[[i]] <- paste(pop[i], "carbon content, carbon quota,", core)
    keyword_carbon_biomass[[i]] <- paste(pop[i], "carbon biomass,", core)
  }

  var_long_name <- c(long_cell_abundance,
                     long_diameter,
                     long_carbon_quota,
                     long_carbon_biomass)

  var_long_name <- stringr::str_to_title(var_long_name)

  var_comment <- c(comment_cell_abundance,
                   comment_diameter,
                   comment_carbon_quota,
                   comment_carbon_biomass)

  var_unit <- c(rep("cells per microliter",length(pop)),
                rep("micrometer", length(pop)),
                rep("picogram carbon per cell", length(pop)),
                rep("microgram carbon per liter", length(pop)))

  var_sensor <- rep("Flow Cytometer",length(var_data))

  var_discipline <- c(rep("Biology",length(pop)),
                      rep("Biology+Biogeochemistry", length(pop)*3))

  visualize <- c(rep(1,length(pop)*4))

  var_keywords <- c(keyword_cell_abundance,
                    keyword_diameter,
                    keyword_carbon_quota,
                    keyword_carbon_biomass)

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
    var_sensor = "Uncategorized",
    var_unit = "",
    var_spatial_res = "irregular",
    var_temporal_res = "irregular",
    var_discipline = "Uncategorized",
    visualize = rep(0, length(id2)),
    var_keywords = core,
    var_comment = "")

  allvars_metadata <- rbind(vars_metadata, custom_metadata)

  # dataset_metadata
  dataset_metadata <- dplyr::tibble(
    dataset_short_name = paste0("Influx_",project,version),
    dataset_long_name = paste0("Discrete Flow Cytometry of ",str_replace_all(project, "_", " "), " Using a BD Influx Cell Sorter"),
    dataset_version = version,
    dataset_release_date = as.Date(Sys.time()),
    dataset_make = "observation",
    dataset_source = "Your Institution/ Your Lab; Armbrust Lab, University of Washington",
    dataset_distributor = "to be added by data owner",
    dataset_acknowledgement = "Samples provided by: Your Name, Your Lab, Your Institution; Data analyzed by: Data provided by: Kelsy Cain and François Ribalet, Armbrust lab, University of Washington",
    dataset_history = "",
    dataset_description = "to be added by data owner",
    dataset_references = "",
    climatology = "",
    cruise_names = cruise)

  # reorder order of column so they match metadata
  data.pivot <- data.pivot[,id]

  ## Save data
  openxlsx::write.xlsx(x=list(data.pivot, dataset_metadata, allvars_metadata),
                       file = paste0("Influx_", project, "_",as.Date(Sys.time()),"_" ,version,".xlsx"),
                       sheetName=c('data','dataset_meta_data','vars_meta_data'),overwrite = TRUE)

}
