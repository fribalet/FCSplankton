# Load packages
library(viridis)
library(tidyverse)
library(FCSplankton)


# set the path to the Project folder
setwd("PATH/TO/PROJECT")

# Set analysis for unstained (TRUE) or stained samples (FALSE)
unstained <- FALSE

# load list of FCS files
if(unstained){folder <- "unstained"
}else{folder <- "stained"}
file_list <- list.files(paste0(folder,"/raw"), pattern = ".fcs$", full.names=T)

### Initialization ###
######################
# load one file
this_file <- file_list[1]
fcs <- read_influx(this_file, transformation = TRUE) # create a dataframe WITH log-amplified data)
# Rename PMTs
print(names(fcs))
id <- c(2,3,4,5) ## replace number by column indice of FSC, 692, 580 and 530 respectively
names(fcs)[id]
names_pmt <- names(fcs)[id] # original names of FSC, 692, 580 and 530 respectively


### Gating populations of interest ###
######################################
gating <- TRUE
summary_table <- NULL
# summary_table <- read_csv(paste0(folder, "summary.csv")) # if you want to append results to an existing summary_table

# Path where to save Gating output
system(paste0("mkdir ",folder, "/gating"))

for (this_file in file_list){

    # Read data
    # this_file <- file_list[1]
    print(paste("gating", this_file))
    fcs <- read_influx(this_file, transformation = TRUE) # create a flowFrame WITH transformation (for log-amplified data)

    # change header
    id <- match(names_pmt, names(fcs))
    names(fcs)[id] <- c("scatter", "red", "orange", "green")

    # Load gating if exist
    if(!gating){ 
      previous <- sub("raw","gating",paste0(this_file, ".RData"))
      if(file.exists(previous)) load(previous)
      }

    # Gates Beads
    if(gating & folder == "unstained") gates.log <- set_gating_params(fcs, "beads", "scatter", "orange")
    if(gating & folder == "stained") gates.log <- set_gating_params(fcs, "beads", "scatter", "red")
    
    # Assign "beads" label to particles
    fcs <- classify_fcs(fcs, gates.log[][1])

    # Normalization
    beads <- fcs[which(fcs$pop == "beads"),]
    fcs$norm.scatter <- fcs$scatter / median(beads$scatter)
    fcs$norm.orange <- fcs$orange / median(beads$orange)
    fcs$norm.red <- fcs$red / median(beads$red)
    fcs$norm.green <- fcs$green / median(beads$green)

    # Gating population - WARNINGS: don't change population names!
    if(gating & folder == "unstained"){
        gates.log <- set_gating_params(fcs, "synecho", "norm.scatter", "norm.orange", gates.log)
        gates.log <- set_gating_params(fcs, "prochloro", "norm.scatter", "norm.red", gates.log)
        gates.log <- set_gating_params(fcs, "picoeuk", "norm.scatter", "norm.red", gates.log)
        #gates.log <- set_gating_params(fcs, "small-picoeuk", "norm_scatter", "norm_red", gates.log)
        #gates.log <- set_gating_params(fcs, "large-picoeuk", "norm_scatter", "norm_red", gates.log)
    }
    if(gating & folder == "stained"){
       gates.log <- set_gating_params(fcs, "bacteria", "norm.scatter", "norm.orange", gates.log)
    }
    
    # Apply gates and label particles according to 'gates.log'
    fcs <- classify_fcs(fcs, gates.log)

    # Save Gating
    if(gating){
        save(gates.log, file=sub("raw","gating",paste0(this_file, ".RData")))
    }

    # Save plot
    if(gating & folder == "stained"){
        png(sub("raw","gating",paste0(this_file, ".png")),width=12, height=6, unit="in", res=200)
        par(mfrow=c(1,2), pty="s", cex=1.2, oma=c(0,0,1,0), mar=c(5,5,1,1))
        plot_vct_cytogram(fcs, "norm.scatter","norm.red", ylab="red\n (normalized to beads)", xlab="scatter\n (normalized to beads)")
        plot_vct_cytogram(fcs, "norm.scatter","norm.orange", ylab="orange\n (normalized to beads)", xlab="scatter\n (normalized to beads)")
        dev.off()
    }
    if(gating & folder == "unstained"){
      png(sub("raw","gating",paste0(this_file, ".png")),width=12, height=6, unit="in", res=200)
      par(mfrow=c(1,2), pty="s", cex=1.2, oma=c(0,0,1,0), mar=c(5,5,1,1))
      plot_vct_cytogram(fcs, "norm.scatter","norm.red", ylab="red\n (normalized to beads)", xlab="scatter\n (normalized to beads)")
      plot_vct_cytogram(fcs, "norm.scatter","norm.orange", ylab="orange\n (normalized to beads)", xlab="scatter\n (normalized to beads)")
      dev.off()
    }
    
    # Aggregate statistics
    stat_table <- NULL
    for(population in unique(fcs$pop)){
        #print(i)

        p <- subset(fcs, pop == population)
        count <- nrow(p)

        if(count == 0) {
            scatter <- 0
            red <- 0
        }else{
            scatter <- round(median(p$norm.scatter),6)
            red <- round(median(p$norm.red),6)
            orange <- round(median(p$norm.orange),6)
            green <- round(median(p$norm.green),6)
        }
        var <- cbind(population,count,scatter,red,orange)
        stat_table <- rbind(stat_table, var)
    }
    table <- data.frame(cbind(stat_table, file=basename(this_file)))

    #  remove entries that already exist
    id <- which(!is.na(match(summary_table$file,unique(table$file))))
     if(length(id) > 0) summary_table <- summary_table[-id,]
    # add replace with new entries
    summary_table <- rbind(summary_table, table)

}


### save summary table ###
##########################
write_csv(summary_table,path=paste0(folder, "/summary.csv"))
