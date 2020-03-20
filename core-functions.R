#' Plot cytogram with only builtin R graphics.
#'
#' @param fcs FCS data frame.
#' @param para.x Channel to use as x axis.
#' @param para.y Channel to use as y axis.
#' @param ... Additional parameters for plot()
#' @return None
#' @usage plot.cytogram(fcs, para.x = "scatter", para.y = "red", ...)
#' @export plot.cytogram
plot.cytogram <- function(fcs, para.x = "scatter", para.y = "red", ...) {
  
  par(pty="s")
  plot(fcs[,c(para.x, para.y)], pch=16, cex=0.3, col = grDevices::densCols(log10(fcs[,c(para.x, para.y)]), colramp = viridis::viridis), log="xy", ...)
}

#' Plot cytogram with particles colored by population.
#'
#' @param fcs FCS data frame.
#' @param para.x Channel to use as x axis.
#' @param para.y Channel to use as y axis.
#' @param ... Additional parameters for plot()
#' @return None
#' @usage plot.vct.cytogram(fcs, para.x = "scatter", para.y = "red")
#' @export plot.vct.cytogram
plot.vct.cytogram <- function (fcs, para.x = "scatter", para.y = "red", ...){
    group.colors <- c(unknown="grey", beads="red3", 
                      prochloro=viridis::viridis(4)[1],
                      synecho=viridis::viridis(4)[2],
                      picoeuk=viridis::viridis(4)[3], 
                      croco=viridis::viridis(4)[4], 
                      small-picoeuk=viridis::viridis(4)[3], 
                      large-picoeuk=viridis::viridis(4)[4])
    fcs$pop <- factor(fcs$pop, levels = names(group.colors))
    caption <- group.colors[unique(fcs$pop)]

    par(pty = "s")
    plot(fcs[, c(para.x, para.y)], pch = 16, cex = 0.3, col = group.colors[fcs$pop], log="xy", ...)
    legend("topleft", legend = names(caption), col = caption,
        pch = 16, pt.cex = 0.6, bty = "n")
    abline(v=1, h=1, col="grey", lty=2)          

}

#' Define polygons for population gating.
#'
#' @param fcs FCS data frame.
#' @param popname Population name.
#' @param para.x Channel to use as x axis.
#' @param para.y Channel to use as y axis.
#' @param poly.log Named list of gating polygon definitions. If a definition for
#'   popname already exists it will be updated. If it doesn't exist it will be
#'   appended to the end to the list. If poly.log is NULL a new list will be
#'   created.
#' @return Version of poly.log with a new polygon defintion for popname.
#' @examples
#' \dontrun{
#' poly.log <- set.gating.params(opp, "beads", "fsc_small", "pe")
#' poly.log <- set.gating.params(opp, "prochloro", "fsc_small", "chl_small",
#'                               poly.log)
#' }
#' @export
set.gating.params <- function(fcs, popname, para.x, para.y, poly.log=NULL) {
  popname <- as.character(popname)
  para.x <- as.character(para.x)
  para.y <- as.character(para.y)

  par(mfrow=c(1,1))
  plot.cytogram(fcs, para.x, para.y)
  mtext(paste("Set Gate for:",popname), font=2)
  poly <- splancs::getpoly(quiet=TRUE) # Draw Gate
  colnames(poly) <- c(para.x, para.y)

  poly.l <- list(poly)
  names(poly.l) <- popname

  if (is.null(poly.log)) {
    # Start a new gating entry
    poly.log <- poly.l
  } else {
    # if gate parameters for the same population already exist, overwrite,
    # otherwise append gate parameters for new population
    poly.log[popname] <- poly.l
  }
  return(poly.log)
}

#' Classify particles based on manually defined population gates.
#'
#' @param fcs FCS data frame.
#' @param params Named list of gating parameters. Must contain a params$poly
#'   entry with polygon definitions.
#' @param popname Name of the population
#' @return List of per particle classifications.
#' @examples
#' \dontrun{
#' vct <- manual.classify(opp, gates.log, "beads")
#' }
#' @export
manual.classify <- function(fcs, params, popname) {
  if (is.null(fcs$pop)) {
    fcs$pop <- "unknown"
  }

  if (is.null(params)) {
    stop(paste0("No gate parameters found for ", popname))
  }

  poly <- params # Get gating polygon definition
  para <- colnames(poly)  # channels
  df <- fcs[fcs$pop=="unknown", para]

  if (nrow(df) > 0) {
    colnames(poly) <- colnames(df) <- c("x","y") # to stop stupid Warnings from splancs::inout()
    vct <- df[splancs::inout(df,poly=poly, bound=TRUE, quiet=TRUE), ] # subset particles based on Gate
    fcs[row.names(vct), "pop"] <- popname
  }

  return(fcs)
}

#' Classify particles from an FCS dataframe.
#'
#' Classify particles from an FCS dataframe using a gating scheme provided by gates.log.
#'
#' @param fcs FCS data frame.
#' @param gates.log A gating scheme from the function "add.manual.classification()" or "add.auto.classification()"
#' @return List of per particle classifications
#' @examples
#' \dontrun{
#' opp <- classify.fcs(fcs, gates.log)
#' }
#' @export
classify.fcs <- function(fcs, gates.log) {
  for (popname in names(gates.log)) {
    params <- gates.log[[popname]]
    fcs <- manual.classify(fcs, params, popname)
   }
  if (! is.null(fcs$pop)) {
    fcs$pop <- factor(fcs$pop)
  }
  return(fcs)
}

