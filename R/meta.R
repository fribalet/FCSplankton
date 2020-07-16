#' Format Influx data into CAMP compatible format.
#'
#' @param data data from Influx.
#' @param cruise cruise name, if any (otherwise NA).
#' @param project cruise name, if any (otherwise NA).
#' @param version Version of the dataset.
#' @return None
#' @examples
#' \dontrun{
#' csv_convert(db, meta, path)
#' }
#' @export
cmap_convert<- function(data, cruise, project, version = "v1.0") {

  ## FORMAT DATA
  # Split population data into columns, and reorder columns in dataset
  data.pivot <- data %>%
          pivot_wider(names_from = population, values_from = c(count, scatter, red, orange, abundance))

  var_cols <- colnames(data.pivot)

  core <- paste("discrete flow cytometry, BD Influx cell sorter, insitu, in-situ, biology, phytoplankton, picophytoplankton, Armbrust, UW, University of Washington,", cruise)

  data_pop <- str_sort(unique(data$population)) # get unique population names

  # master list of possible population names; add to as needed
  long_pop <- c("prochlorococcus", "synechococcus", "picoeukaryote", "large-picoeukaryote", "small-picoeukaryote", "crocosphaera", "bacteria", "unknown", "beads")
  # anchor population name to subset correct population names based on matching the beginning of population names
  match_pop <- str_replace_all(data_pop, "^", "^")
  # loop through list of populations from data to get list of full population names
  pop <- c()
  for(i in 1:length(data_pop)){
    pop[[i]] <- str_subset(long_pop, match_pop[i])
  }

  ## FORMAT METADATA
  #set order of standard variable metadata
  var_data <- c("file" ,
                str_sort(str_subset(var_cols, "count_")),
                str_sort(str_subset(var_cols, "scatter_")),
                str_sort(str_subset(var_cols, "red_")),
                str_sort(str_subset(var_cols, "orange_")),
                "volume",
                str_sort(str_subset(var_cols, "abundance_")),
                "flag",
                "stain")


  long_count <- c()
  long_scatter <- c()
  long_red <- c()
  long_orange <- c()
  long_abundance <- c()
  comment_count <- c()
  comment_scatter <- c()
  comment_red <- c()
  comment_orange <- c()
  comment_abundance <- c()
  keyword_count <- c()
  keyword_scatter <- c()
  keyword_red <- c()
  keyword_orange <- c()
  keyword_abundance <- c()

  # run through loop to add all the population dependent metadata
  for(i in 1:length(pop)){
    long_count[[i]] <- paste("number of", pop[i], "particles counted by the instrument")
    long_scatter[[i]] <- paste("50% percentile of", pop[i], "forward angle light scatter normalized to 1 micron beads (proxy of cell diameter)")
    long_red[[i]] <- paste("50% percentile of", pop[i], "red fluorescence normalized to 1 micron beads (proxy of chlorophyll content)")
    long_orange[[i]] <- paste("50% percentile of", pop[i], "orange fluorescence normalized to 1 micron beads (proxy of phycoerythrin content)")
    long_abundance[[i]] <- paste(pop[i], "cell abundance")
    comment_count[[i]] <- paste(pop[i], "count needs to be > 30 to be trusted")
    comment_scatter[[i]] <- paste(pop[i], "light scatter collected using a 457-50 bandpass filter")
    comment_red[[i]] <- paste(pop[i], "red fluorescence collected using a 692-40 bandpass filter")
    comment_orange[[i]] <- paste(pop[i], "orange fluorescence collected using a 572-27 bandpass filter")
    comment_abundance[[i]] <- paste(pop[i], "cell abundance, number of particles divided by the volume of sample, see https://github.com/fribalet/FCSplankton for more details")
    keyword_count[[i]] <- paste(pop[i], "particle count,", core)
    keyword_scatter[[i]] <- paste(pop[i], "forward angle light scatter, FSC, FALS,", core)
    keyword_red[[i]] <- paste(pop[i], "red fluorescence, chlorophyll,", core)
    keyword_orange[[i]] <- paste(pop[i], "orange fluorescence, phycoerythrin,", core)
    keyword_abundance[[i]] <- paste(pop[i],"abundance, cell concentration, cell count, cell abundance,", core)
  }

  var_long_name <- c("flow cytometry standard file (.fcs)",
                     long_count,
                     long_scatter,
                     long_red,
                     long_orange,
                     "volume measured by pressure sensor",
                     long_abundance,
                     "sample status flag",
                     "sample stain flag")

  var_comment <- c("",
                   comment_count,
                   comment_scatter,
                   comment_red,
                   comment_orange,
                   "volume measured by pressure sensor",
                   comment_abundance,
                   "outliers (0 = quality data; 1 = issue related to instrument performance; 2 = issue related to population classification)",
                   "DNA stain (0 = unstained sample; 1 = stained sample)")

  var_unit <- c("",
                rep("",length(pop)*4),
                "microliter",
                rep("cells/μL",length(pop)),
                "",
                "")

  var_sensor <- rep("Flow cytometry",length(var_data))

  var_discipline <- c("",
                     rep("Biology",length(pop)*4),
                     "",
                     rep("Biology",length(pop)),
                     "",
                     "")

  visualize <- c(0, rep(1,length(pop)*4), 0, rep(1,length(pop)), 0, 0)

  var_keywords <- c(paste("file", core),
        keyword_count,
        keyword_scatter,
        keyword_red,
        keyword_orange,
        paste("volume, vol,", core),
        keyword_abundance,
        paste("flag,", core),
        paste("DNA stain, SYBR stain,", core))


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

  # Save data
  openxlsx::write.xlsx(x=list(data.pivot, dataset_metadata, allvars_metadata),
                      file = paste0("Influx_", project, "_",as.Date(Sys.time()),"_" ,version,".xlsx"),
                      sheetName=c('data','dataset_meta_data','vars_meta_data'))

}
