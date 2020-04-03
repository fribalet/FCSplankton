var_data <- c("time", "lat", "lon",
              "file",
              "population",
              "count",
              "scatter","red","orange",
              "volume",
              "abundance",
              "flag")

var_standard_name <- c("time",
                  "latitude",
                  "longitude",
                  "filename",
                  "population name",
                  "particle count",
                  "median of forward angle light scatter",
                  "median of red fluorescence",
                  "median of orange fluorescence",
                  "sample volume",
                  "cell concentration",
                  "status flag")

var_long_name <- c("time of sample collection (UTC)",
                  "latitude",
                  "longitude",
                  "flow cytometry standard file (.fcs)",
                  "name of cytometric population",
                  "number of particles counted by the instrument",
                  "50% percentile of forward angle light scatter (proxy of cell diameter)",
                  "50% percentile of red fluorescence (proxy of chlorophyll content)",
                  "50% percentile of orange fluorescence (proxy of phycoerythrin content)",
                  "volume of sample analyzed",
                  "cell abundance",
                  "outliers")

var_comment <-  c(rep("none",4),
                  "prochloro (Prochlorococcus) synecho (Synechococcus) picoeuk (large phytoplankton) beads (internal standard) croco (Crocosphaera-like particles) unknown (unclassified particles)",
                  "needs to be > 30 to be trusted",
                  "light scatter collected using a 457-50 bandpass filter",
                  "red fluorescence collected using a 692-40 bandpass filter",
                  "orange fluorescence collected using a 572-27 bandpass filter",
                  "volume measured by pressure sensor",
                  "number of particles divided by the volume of sample",
                  'outliers (0 = Quality data; 1 = issue related to instrument performance; 2 = issue related to population classification)')

var_unit <- c("%Y-%m-%dT%H:%M:%S",
              "decimal degree North",
              "decimal degree East",
              rep("unitless",6),
              "microliter",
              "cells per microliter",
              "unitless")

var_keywords <- c("time+UTC+date",
                "latitude",
                "longitude",
                "",
                "Prochlorococcus+Synechococcus+Crocosphaera+picoeukaryotes+phytoplankton+picophytoplankton+unknown",
                "",
                "forward+angle+light+scatter+FSC+FALS",
                "red+fluorescence+chlorophyll",
                "orange+fluorescence+chlorophyll",
                "",
                "abundance+concentration+density",
                "")

var_sensor <- rep("Influx",12)

var_discipline <- c(rep("", 4),
                   "taxonomy+cytometry+Influx",
                   rep("optics+cytometry+Influx",4),
                   "",
                   "biology+cytometry+Influx",
                   "")





#' Format Influx data into an Excel spreadsheet, along with metadata.
#'
#' @param data data from Influx.
#' @param project Version of the dataset.
#' @param cruise cruise name, if any (otherwise NA).
#' @param version Version of the dataset.
#' @return None
#' @examples
#' \dontrun{
#' csv_convert(db, meta, path)
#' }
#' @export
xls_convert<- function(data, project, cruise, version = "v1.0") {

    # add custom column to metadata
    id1 <- match(var_data,colnames(data)) # which column are standard
    id2 <- which(is.na(match(1:ncol(data), id1))) # which column are not standard
    id <- c(id1,id2) 
    
    # vars_metadata
    vars_metadata <- dplyr::tibble(
                          var_short_name = var_data,
                          var_standard_name,
                          var_sensor,
                          var_unit,
                          var_spatial_res = "irregular",
                          var_temporal_res = "irregular",
                          var_missing_value = "NA",
                          var_discipline,
                          var_keywords,
                          var_comment)

    # custom metadata
    custom_metadata <- dplyr::tibble(
                          var_short_name = colnames(data[,id2]),
                          var_standard_name = "NA",
                          var_sensor = "NA",
                          var_unit = "NA",
                          var_spatial_res = "irregular",
                          var_temporal_res = "irregular",
                          var_missing_value = "NA",
                          var_discipline = "NA",
                          var_keywords = "NA",
                          var_comment = "NA")

    allvars_metadata <- rbind(vars_metadata, custom_metadata)

    # dataset_metadata
    dataset_metadata <- dplyr::tibble(
                          dataset_short_name = paste0("Influx_",project),
                          dataset_long_name = paste0("Influx_",project),
                          dataset_cruise = cruise,
                          dataset_project = project,
                          dataset_version = version,
                          dataset_release_date = as.Date(Sys.time()),
                          dataset_make = "observation",
                          dataset_source = "Influx data generated by University of Washington / Armbrust lab (kcain97@uw.edu, ribalet@uw.edu)",
                          dataset_doi = "to be added by the data owner",
                          dataset_history = "",
                          dataset_description = "to be added by data owner",
                          dataset_references = "")
    
    # reorder columns in dataset
    data <- data[,id]
                   
    openxlsx::write.xlsx(x=list(data, allvars_metadata, dataset_metadata), 
                         file = paste0("Influx_", project, "_dataset_", version,".xlsx"), 
                         sheetName=c('data','vars_metadata','dataset_metadata'))

}
