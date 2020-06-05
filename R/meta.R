var_data <- c("time", "lat", "lon", "depth",
              "file",
              "population",
              "count",
              "scatter","red","orange",
              "volume",
              "abundance",
              "flag",
              "stain")

var_standard_name <- c("time",
                  "latitude",
                  "longitude",
                  "depth",
                  "filename",
                  "population name",
                  "particle count",
                  "median of forward angle light scatter",
                  "median of red fluorescence",
                  "median of orange fluorescence",
                  "sample volume",
                  "cell concentration",
                  "status flag",
                  "stain flag")


var_long_name <- c("time of sample collection (UTC)",
                  "latitude",
                  "longitude",
                  "depth",
                  "flow cytometry standard file (.fcs)",
                  "name of cytometric population",
                  "number of particles counted by the instrument",
                  "50% percentile of forward angle light scatter (proxy of cell diameter)",
                  "50% percentile of red fluorescence (proxy of chlorophyll content)",
                  "50% percentile of orange fluorescence (proxy of phycoerythrin content)",
                  "volume of sample analyzed",
                  "cell abundance",
                  "sample status flag",
                  "sample stain flag")

var_comment <-  c(rep("",5),
                  "prochloro (Prochlorococcus) synecho (Synechococcus) picoeuk (large phytoplankton) beads (internal standard) croco (Crocosphaera-like particles) bacteria (heterotrophic bacteria) unknown (unclassified particles)",
                  "needs to be > 30 to be trusted",
                  "light scatter collected using a 457-50 bandpass filter",
                  "red fluorescence collected using a 692-40 bandpass filter",
                  "orange fluorescence collected using a 572-27 bandpass filter",
                  "volume measured by pressure sensor",
                  "cell abundance, number of particles divided by the volume of sample, see https://github.com/fribalet/FCSplankton for more details",
                  "outliers (0 = quality data; 1 = issue related to instrument performance; 2 = issue related to population classification)",
                  "DNA stain (0 = unstained sample; 1 = stained sample)")

var_unit <- c("%Y-%m-%dT%H:%M:%S",
              "decimal degree North",
              "decimal degree East",
              "meter",
              rep("",6),
              "microliter",
              "cells/μL",
              "",
              "")

var_keywords <- c("time, UTC, date",
                "latitude, lat",
                "longitude, lon",
                "depth",
                paste0("file, flow, cytometry, flow cytometry, discrete flow cytometry, influx, BD Influx cell sorter, insitu, in-situ, bio, biology, armbrust, UW, University of Washington, ", cruise,", ", cruise_common,", ",project),
                paste0("Prochlorococcus, pro, prochloro, Synechococcus, syn, synecho, Crocosphaera, croc, croco bacteria, picoeukaryotes, pico, picoeuks, phytoplankton, picophytoplankton, unknown, bacteria, heterotrophic bacteria, flow, cytometry, flow cytometry, discrete flow cytometry, influx, BD Influx cell sorter, insitu, in-situ, bio, biology, armbrust, UW, University of Washington, ",cruise,", ", cruise_common,", ",project),
                paste0("particle, particle count, flow, cytometry, flow cytometry, discrete flow cytometry, influx, BD Influx cell sorter, insitu, in-situ, bio, biology, armbrust, UW, University of Washington", cruise,", ", cruise_common,", ",project),
                paste0("forward angle light scatter, FSC, FALS, flow, cytometry, flow cytometry, discrete flow cytometry, influx, BD Influx cell sorter, insitu, in-situ, bio, biology, armbrust, UW, University of Washington", cruise,", ", cruise_common,", ",project),
                paste0("red fluorescence, chlorophyll, fluorescence, flow, cytometry, flow cytometry, discrete flow cytometry, influx, BD Influx cell sorter, insitu, in-situ, bio, biology, armbrust, UW, University of Washington", cruise,", ", cruise_common,", ",project),
                paste0("orange fluorescence, phycoerythrin, fluorescence, synechococcus, syn, synecho, flow, cytometry, flow cytometry, discrete flow cytometry, influx, BD Influx cell sorter, insitu, in-situ, bio, biology, armbrust, UW, University of Washington", cruise,", ", cruise_common,", ",project),
                paste0("volume, vol, flow, cytometry, flow cytometry, discrete flow cytometry, influx, BD Influx cell sorter, insitu, in-situ, bio, biology, armbrust, UW, University of Washington", cruise,", ", cruise_common,", ",project),
                paste0("abundance, cell concentration, concentration, cell count, cell abundance, prochloro abundance, Prochlorococcus Abundance, pro abundance, synecho abundance, synechococcus abundance, syn abundance, bacteria abundance, heterotrophic bacteria abundace, Crocosphaera abundance, croco abundance, picoeukaryote abundance, pico abundance, picoeuk abundance, flow, cytometry, flow cytometry, discrete flow cytometry, influx, BD Influx cell sorter, insitu, in-situ, bio, biology, armbrust, UW, University of Washington", cruise,", ", cruise_common,", ",project),
                paste0("flag, flow, cytometry, flow cytometry, discrete flow cytometry, influx, BD Influx cell sorter, insitu, in-situ, bio, biology, armbrust, UW, University of Washington", cruise,", ", cruise_common,", ",project),
                paste0("stain, DNA stain, SYBR, SYBR stain, bacteria, heterotrophic bacteria, bacteria abundance, flow, cytometry, flow cytometry, discrete flow cytometry, influx, BD Influx cell sorter, insitu, in-situ, bio, biology, armbrust, UW, University of Washington", cruise,", ", cruise_common,", ",project))

var_sensor <- rep("Flow cytometry",14)

var_discipline <- c(rep("", 5),
                   rep("Biology",5),
                   "",
                   "Biology",
                   "",
                   "")





#' Format Influx data into an Excel spreadsheet, along with metadata.
#'
#' @param data data from Influx.
#' @param cruise cruise name, if any (otherwise NA).
#' @param version Version of the dataset.
#' @return None
#' @examples
#' \dontrun{
#' csv_convert(db, meta, path)
#' }
#' @export
xls_convert<- function(data, cruise, version = "v1.0") {

    # add custom column to metadata
    id1 <- match(var_data,colnames(data)) # which column are standard
    id2 <- which(is.na(match(1:ncol(data), id1))) # which column are not standard
    id <- c(id1,id2)

    # vars_metadata
    vars_metadata <- dplyr::tibble(
                          var_short_name = var_data,
                          var_long_name,
                          var_sensor,
                          var_unit,
                          var_spatial_res = "irregular",
                          var_temporal_res = "irregular",
                          var_discipline,
                          var_keywords,
                          var_comment)

    # custom metadata
    custom_metadata <- dplyr::tibble(
                          var_short_name = colnames(data[,id2]),
                          var_long_name = "",
                          var_sensor = "",
                          var_unit = "",
                          var_spatial_res = "irregular",
                          var_temporal_res = "irregular",
                          var_discipline = "",
                          var_keywords = paste0("flow, cytometry, flow cytometry, discrete flow cytometry, influx, BD Influx cell sorter, insitu, in-situ, bio, biology, armbrust, UW, University of Washington, ", cruise,", ", cruise_nickname,", ",project),
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
                          official_cruise_name = cruise,
                          dataset_acknowledgement = "",
                          contact_email = "kcain97@uw.edu, ribalet@uw.edu",
                          dataset_description = "to be added by data owner",
                          dataset_references = "",
                          climatology = NULL)

    # reorder columns in dataset
    data <- data[,id]

    openxlsx::write.xlsx(x=list(data, dataset_metadata, allvars_metadata),
                         file = paste0("Influx_", project, "_dataset_", version,".xlsx"),
                         sheetName=c('data','dataset_meta_data','vars_meta_data'))

}
