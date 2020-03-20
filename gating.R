### 1. Load libraries and other goodies
#######################################
library(viridis)
library(tidyverse)
library(FCS_analysis)

### 2. Gating populations of interest
#####################################
# set the path to the FCS files
setwd("/Volumes/GoogleDrive/Shared\ drives/Influx/Light_Dark_Experiment_2019")

# path to save analysis
system("mkdir gating")

# load list of FCS files
file.list <- dir(".", pattern = ".fcs$", recursive=F)

# # load metadata
# meta <- read.csv(file = paste0("metadata.csv"), header = T, sep = ",")
# ok_meta <- subset(meta, Flag == 0)      # only the unflagged files


# a. Initialization
this_file <- file.list[1]
fcs <- caroline::tab2df(flowCore::exprs(flowCore::read.FCS(this_file, transformation = T, emptyValue=F))) # create a flowFrame WITH transformation (for log-amplified data)

    # Rename PMTs
    print(colnames(fcs))
    id <- c(2,3,4) ## replace number by column indice of FSC, 692, 580 respectively
    colnames(fcs)[id]
    names.pmt <- colnames(fcs)[id] # original names of FSC, 692, 580 respectively

gating <- TRUE

summary.table <- NULL
# summary.table <- read_csv("gating/summary.csv") # if you want to append results to an existing summary.table.

for (this_file in file.list){

    ### Read data
    # this_file <- file.list[1]
    print(paste("gating", this_file))
    fcs <- caroline::tab2df(flowCore::exprs(flowCore::read.FCS(this_file, transformation = T, emptyValue=F))) # create a flowFrame WITH transformation (for log-amplified data)

    # change header
    id <- match(names.pmt, colnames(fcs))
    colnames(fcs)[id] <- c("scatter", "red", "orange")
    
    # Load gating if exist
    if(!gating) load(file=paste0("gating/",this_file, ".RData"))

    ### Beads Normalization
    # Gates Beads
    if(gating) gates.log <- set.gating.params(fcs, "beads", "scatter", "orange")
    fcs <- classify.fcs(fcs, gates.log)
    
    # Normalization
    beads <- fcs[which(fcs$pop == "beads"),]
    fcs$norm.scatter <- fcs$scatter / median(beads$scatter)
    fcs$norm.orange <- fcs$orange / median(beads$orange)
    fcs$norm.red <- fcs$red / median(beads$red)

    ### Gating population - WARNINGS: don't change population names!
    if(gating){
        gates.log <- set.gating.params(fcs, "synecho", "norm.scatter", "norm.orange", gates.log)
        #gates.log <- set.gating.params(fcs, "prochloro", "norm.scatter", "norm.red", gates.log)
        #gates.log <- set.gating.params(fcs, "picoeuk", "norm.scatter", "norm.red", gates.log)
        gates.log <- set.gating.params(fcs, "small-picoeuk", "norm.scatter", "norm.red", gates.log)
        gates.log <- set.gating.params(fcs, "large-picoeuk", "norm.scatter", "norm.red", gates.log)
   }
    fcs$pop <- "unknown"
    fcs <- classify.fcs(fcs, gates.log)

    ### Save Gating
    if(gating){ 
        save(gates.log, file=paste0("gating/",this_file,".RData"))
    }
   
    ### Save plot
    if(gating){
        png(paste0("gating/",this_file,".png"),width=12, height=6, unit="in", res=200)
        par(mfrow=c(1,2), pty="s", cex=1.2, oma=c(0,0,1,0), mar=c(5,5,1,1))
        plot.vct.cytogram(fcs, "norm.scatter","norm.red", ylab="red\n (normalized to beads)", xlab="scatter\n (normalized to beads)")
        plot.vct.cytogram(fcs, "norm.scatter","norm.orange", ylab="orange\n (normalized to beads)", xlab="scatter\n (normalized to beads)")
        mtext(paste(basename(this_file)), outer=TRUE)
        dev.off()
    }

    ### Aggregate statistics 
    stat.table <- NULL
    for(population in unique(fcs$pop)){
        #print(i)
   
        p <- subset(fcs, pop == population)
        n <- nrow(p)
    
        if(n ==0) {
            scatter <- 0
            red <- 0
        }else{
            scatter <- round(median(p$norm.scatter),6)
            red <- round(median(p$norm.red),6)
            orange <- round(median(p$norm.orange),6)
        }
        var <- cbind(population,n,scatter,red,orange)
        stat.table <- rbind(stat.table, var)
    }


    table <- data.frame(cbind(stat.table, file=basename(this_file)))
    summary.table <- rbind(summary.table, table)

}


write.csv(summary.table,file="gating/summary.csv", row.names=FALSE, quote=FALSE)
