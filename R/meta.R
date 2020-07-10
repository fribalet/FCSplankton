#' Format Influx data into CAMP compatible format.
#'
#' @param data data from Influx.
#' @param cruise cruise name, if any (otherwise NA).
#' @param project cruise name, if any (otherwise NA).
#' @param unstained Whether samples were unstained (no fluorochrome added) or stained (e.g. SYBRGreen). Default is TRUE
#' @param version Version of the dataset.
#' @return None
#' @examples
#' \dontrun{
#' csv_convert(db, meta, path)
#' }
#' @export
cmap_convert<- function(data, cruise, project, unstained=TRUE, version = "v1.0") {

  ## FORMAT DATA
  # Split population data into columns, and reorder columns in dataset
  data <- data %>%
          pivot_wider(names_from = population, values_from = c(count, scatter, red, orange, abundance))


  ## FORMAT METADATA
  var_data <- c("file",
              "count_prochloro","count_picoeuk","count_synecho","count_beads","count_unknown","count_bacteria",
              "scatter_prochloro","scatter_picoeuk","scatter_synecho","scatter_beads","scatter_unknown","scatter_bacteria",
              "red_prochloro","red_picoeuk","red_synecho","red_beads","red_unknown","red_bacteria",
              "orange_prochloro","orange_picoeuk","orange_synecho","orange_beads","orange_unknown","orange_bacteria",
              "volume",
              "abundance_prochloro","abundance_picoeuk","abundance_synecho","abundance_beads","abundance_unknown","abundance_bacteria",
              "flag","stain")

  var_standard_name <- c("filename",
                    "prochlorococcus particle count","picoeukaryote particle count","synechococcus particle count","beads particle count","unknown particle count","bacteria particle count",
                    "median of prochlorococcus forward angle light scatter","median of picoeukaryote forward angle light scatter","median of synechococcus forward angle light scatter","median of beads forward angle light scatter","median of unknown forward angle light scatter","median of bacteria forward angle light scatter",
                    "median of prochlorococcus red fluorescence","median of picoeukaryote red fluorescence","median of synechococcus red fluorescence","median of beads red fluorescence","median of unknown red fluorescence","median of bacteria red fluorescence",
                    "median of prochlorococcus orange fluorescence","median of picoeukaryote orange fluorescence","median of synechococcus orange fluorescence","median of beads orange fluorescence","median of unknown orange fluorescence","median of bacteria orange fluorescence",
                    "sample volume",
                    "prochlorococcus cell concentration","picoeukaryote cell concentration","synechococcus cell concentration","beads cell concentration","unknown cell concentration","bacteria cell concentration",
                    "status flag",
                    "stain flag")

  var_long_name <- c("flow cytometry standard file (.fcs)",
                    "number of prochlorococcus particles counted by the instrument","number of picoeukaryote particles counted by the instrument","number of synechococcus particles counted by the instrument","number of bead particles counted by the instrument","number of unknown particles counted by the instrument","number of bacteria particles counted by the instrument",
                    "50% percentile of prochlorococcus forward angle light scatter (proxy of cell diameter)","50% percentile of picoeukaryote forward angle light scatter (proxy of cell diameter)","50% percentile of synechococcus forward angle light scatter (proxy of cell diameter)","50% percentile of bead forward angle light scatter (proxy of cell diameter)","50% percentile of unknown forward angle light scatter (proxy of cell diameter)","50% percentile of bacteria forward angle light scatter (proxy of cell diameter)",
                    "50% percentile of prochlorococcus red fluorescence (proxy of chlorophyll content)","50% percentile of picoeukaryote red fluorescence (proxy of chlorophyll content)","50% percentile of synechococcus red fluorescence (proxy of chlorophyll content)","50% percentile of bead red fluorescence (proxy of chlorophyll content)","50% percentile of unknown red fluorescence (proxy of chlorophyll content)","50% percentile of bacteria red fluorescence (proxy of chlorophyll content)",
                    "50% percentile of prochlocococcus orange fluorescence (proxy of phycoerythrin content)","50% percentile of picoeukaryote orange fluorescence (proxy of phycoerythrin content)","50% percentile of synechococcus orange fluorescence (proxy of phycoerythrin content)","50% percentile of beads orange fluorescence (proxy of phycoerythrin content)","50% percentile of unknown orange fluorescence (proxy of phycoerythrin content)","50% percentile of bacteria orange fluorescence (proxy of phycoerythrin content)",
                    "volume of sample analyzed",
                    "prochlorococcus cell abundance","picoeukaryote cell abundance","synechococcus cell abundance","beads cell abundance","unknown cell abundance","bacteria cell abundance",
                    "sample status flag","sample stain flag")

  var_comment <-  c("",
                  "prochlorococcus count needs to be > 30 to be trusted","picoeukaryote count needs to be > 30 to be trusted","synechococcus count needs to be > 30 to be trusted","bead count needs to be > 30 to be trusted","unknown count needs to be > 30 to be trusted","bacteria count needs to be > 30 to be trusted",
                  "prochlorococcus light scatter collected using a 457-50 bandpass filter","picoeukaryote light scatter collected using a 457-50 bandpass filter","synechococcus light scatter collected using a 457-50 bandpass filter"," bead light scatter collected using a 457-50 bandpass filter","unknown light scatter collected using a 457-50 bandpass filter","bacteria light scatter collected using a 457-50 bandpass filter",
                  "prochlorococcus red fluorescence collected using a 692-40 bandpass filter","picoeukaryote red fluorescence collected using a 692-40 bandpass filter","synechococcus red fluorescence collected using a 692-40 bandpass filter", "bead red fluorescence collected using a 692-40 bandpass filter","unknown red fluorescence collected using a 692-40 bandpass filter","bacteria red fluorescence collected using a 692-40 bandpass filter",
                  "prochlorococcus orange fluorescence collected using a 572-27 bandpass filter","picoeukaryote orange fluorescence collected using a 572-27 bandpass filter","synechococcus orange fluorescence collected using a 572-27 bandpass filter","bead orange fluorescence collected using a 572-27 bandpass filter","unknown orange fluorescence collected using a 572-27 bandpass filter","bacteria orange fluorescence collected using a 572-27 bandpass filter",
                  "volume measured by pressure sensor",
                  "prochlorococcus cell abundance, number of particles divided by the volume of sample, see https://github.com/fribalet/FCSplankton for more details","picoeukaryote cell abundance, number of particles divided by the volume of sample, see https://github.com/fribalet/FCSplankton for more details","synechococcus cell abundance, number of particles divided by the volume of sample, see https://github.com/fribalet/FCSplankton for more details",
                  "bead cell abundance, number of particles divided by the volume of sample, see https://github.com/fribalet/FCSplankton for more details","unknown cell abundance, number of particles divided by the volume of sample, see https://github.com/fribalet/FCSplankton for more details","bacteria cell abundance, number of particles divided by the volume of sample, see https://github.com/fribalet/FCSplankton for more details",
                  "outliers (0 = quality data; 1 = issue related to instrument performance; 2 = issue related to population classification)",
                  "DNA stain (0 = unstained sample; 1 = stained sample)")

  var_unit <- c(rep("",25),
                "microliter",
                rep("cells/μL",6),
                "",
                "")

  var_sensor <- rep("Flow cytometry",34)

  var_discipline <- c("",
                     rep("Biology",24),
                     "",
                     rep("Biology",6),
                     "",
                     "")

  visualize <- c(rep(0,7), rep(1,18), 0, rep(1,6), 0, 0)

  core <- paste("discrete flow cytometry, BD Influx cell sorter, insitu, in-situ, biology, Armbrust, UW, University of Washington", cruise)

  var_keywords <- c(paste("file", core),
        paste("prochlorococcus particle count, phytoplankton,",core),
        paste("picoeukaryote particle count, picophytoplankton,",core),
        paste("synechococcus particle count, phytoplankton,",core),
        paste("bead particle count,",core),
        paste("unknown, particle count,",core),
        paste("bacteria, heterotrophic bacteria, particle, bacteria particle count,", core),
        paste("prochlorococcus forward angle light scatter, FSC, FALS, phytoplankton,",core),
        paste("picoeukaryote forward angle light scatter, FSC, FALS, picophytoplankton,",core),
        paste("synechococcus forward angle light scatter, FSC, FALS, phytoplankton,",core),
        paste("bead forward angle light scatter, FSC, FALS, beads,",core),
        paste("unknown forward angle light scatter, FSC, FALS, unknown,",core),
        paste("bacteria forward angle light scatter, FSC, FALS,", core),
        paste("prochlorococcus red fluorescence, chlorophyll, phytoplankton,",core),
        paste("picoeukaryote red fluorescence, chlorophyll, picophytoplankton,",core),
        paste("synechococcus red fluorescence, chlorophyll, phytoplankton,",core),
        paste("bead red fluorescence, chlorophyll,",core),
        paste("unknown red fluorescence, chlorophyll,",core),
        paste("bacteria red fluorescence, chlorophyll, fluorescence,", core),
        paste("prochlorococcus orange fluorescence, phycoerythrin, phytoplankton,",core),
        paste("picoeukaryote orange fluorescence, phycoerythrin, picophytoplankton,",core),
        paste("synecochoccus orange fluorescence, phycoerythrin, phytoplankton,",core),
        paste("bead orange fluorescence, phycoerythrin,",core),
        paste("unknown orange fluorescence, phycoerythrin,",core),
        paste("bacteria orange fluorescence, phycoerythrin,", core),
        paste("volume, vol,", core),
        paste("prochlorococcus abundance, cell concentration, cell count, cell abundance, phytoplankton,",core),
        paste("pikoeukaryote abundance, cell concentration, cell count, cell abundance, picophytoplankton,",core),
        paste("synechococcus abundance, cell concentration, cell count, cell abundance, phytoplankton,",core),
        paste("bead abundance, cell concentration, cell count, cell abundance,",core),
        paste("unknown abundance, cell concentration, cell count, cell abundance,",core),
        paste("bacteria abundance, cell concentration, cell count, cell abundance,", core),
        paste("flag,", core),
        paste("DNA stain, SYBR stain", core))

  if(unstained){
    var_keywords <- var_keywords[!str_detect(var_keywords,pattern="bacteria")]
    var_data <- var_data[!str_detect(var_data,pattern="bacteria")]
    var_standard_name <- var_standard_name[!str_detect(var_standard_name,pattern="bacteria")]
    var_long_name <- var_long_name[!str_detect(var_long_name,pattern="bacteria")]
    var_comment <- var_comment[!str_detect(var_comment,pattern="bacteria")]
    var_unit <- c(rep("",21), "microliter", rep("cells/μL",5), "", "")
    var_sensor <- rep("Flow cytometry",29)
    var_discipline <- c("", rep("Biology",20), "", rep("Biology",5), "", "")
    visualize <- c(rep(0,6), rep(1,15), 0, rep(1,5), 0, 0)
  }


  # add custom column to metadata
  coordinates <- c("time","lat","lon","depth")
  id0 <- match(coordinates,colnames(data)) # which column are standard
  id1 <- match(var_data,colnames(data)) # which column are standard
  id2 <- which(is.na(match(1:ncol(data), c(id0, id1)))) # which column are not standard
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
                          var_short_name = colnames(data)[id2],
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
                        dataset_long_name = paste0("Influx_",project),
                        dataset_version = version,
                        dataset_release_date = as.Date(Sys.time()),
                        dataset_make = "observation",
                        dataset_source = "University of Washington / Armbrust lab",
                        "official_cruise_name(s)" = cruise,
                        dataset_acknowledgement = "",
                        contact_email = "kcain97@uw.edu, ribalet@uw.edu",
                        dataset_description = "to be added by data owner",
                        dataset_references = "",
                        climatology = NULL)

  # reorder order of column so they match metadata
  data <- data[,id]

  # Save data
  openxlsx::write.xlsx(x=list(data, dataset_metadata, allvars_metadata),
                      file = paste0("Influx_", project, "_",as.Date(Sys.time()),"_" ,version,".xlsx"),
                      sheetName=c('data','dataset_meta_data','vars_meta_data'))

}
