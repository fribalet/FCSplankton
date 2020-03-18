## 1. Install dependencies
##########################
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("flowCore",
                        "splanc", 
                        "caroline", 
                        "viridis",
                        "tidyverse",
                        dependencies = TRUE)

### 2. Load libraries and other goodies
#######################################
library(flowCore)
library(splancs)
library(caroline)
library(viridis)
library(tidyverse)
library(grDevices)

group.colors <- c(unknown="grey", beads="red3", prochloro=viridis::viridis(4)[1],synecho=viridis::viridis(4)[2],picoeuk=viridis::viridis(4)[3], croco=viridis::viridis(4)[4])


plot.cytogram <- function (evtopp, para.x = "scatter", para.y = "red", ...){
    cols <- colorRampPalette(c("blue4", "royalblue4", "deepskyblue3",
        "seagreen3", "yellow", "orangered2", "darkred"))
    par(pty = "s")
    plot(evtopp[, c(para.x, para.y)], pch = 16, cex = 0.3,col = densCols(log10(evtopp[, c(para.x, para.y)]),
          colramp = viridis), log="xy",...)
}

plot.vct.cytogram <- function (opp, para.x = "scatter", para.y = "red", ...){
    opp$pop <- factor(opp$pop, levels = names(group.colors))
    caption <- group.colors[unique(opp$pop)]

    plot(opp[, c(para.x, para.y)], pch = 16, cex = 0.3, col = group.colors[opp$pop], log="xy", ...)
    legend("topleft", legend = names(caption), col = caption,
        pch = 16, pt.cex = 0.6, bty = "n")
    abline(v=1, h=1, col="grey", lty=2)          

}




### 3. Gating populations of interest
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
opp <- caroline::tab2df(exprs(read.FCS(this_file, transformation = T, emptyValue=F))) # create a flowFrame WITH transformation (for log-amplified data)

    # Rename PMTs
    print(colnames(opp))
    id <- c(2,3,4) ## replace number by column indice of FSC, 692, 580 respectively

gating <- TRUE

summary.table <- NULL

for (this_file in file.list){

    # this_file <- file.list[3]
    print(paste("gating", this_file))
    opp <- caroline::tab2df(exprs(read.FCS(this_file, transformation = T, emptyValue=F))) # read FCS file (Transformation = TRUE for log-amplified data)
    colnames(opp)[c(2,3,4)] <- c("scatter", "red", "orange")
    opp$pop <- "unknown" # add a colum "pop" to the table


    # # Load gating if exist
    # if(!gating){    
    #         load(file=paste0("gating/",this_file, ".RData")
    #         b.gates <- gates[[1]]
    #         syn.gates <- gates[[2]]
    #         pro.gates <- gates[[3]]
    #         pico.gates <- gates[[4]]
    #        }

    ### Beads Normalization
    # a. Gate beads
    params1 <- c("scatter","orange")
    if(gating){
            plot.cytogram(opp, params1[1], params1[2], main="Gate Beads")
            b.gates <- splancs::getpoly(); colnames(b.gates) <- params1 # draw gate
            polygon(b.gates, lwd=2,  border="red3")
            }
    # b. Filter beads
    beads <- opp[inout(opp[,params1], b.gates),] 
    opp[row.names(beads),"pop"] <- "beads"

    # c. Normalization
        opp$norm.scatter <- opp$scatter / median(beads$scatter)
        opp$norm.orange <- opp$orange / median(beads$orange)
        opp$norm.red <- opp$red / median(beads$red)

    # d. Gate Synecho
    params2 <- c("norm.scatter","norm.orange")
    x <- subset(opp, pop=="unknown")
    if(gating){
            par(mfrow=c(1,1), pty="s")
            plot.cytogram(x, params2[1], params2[2], main="Gate Synechococcus")
            syn.gates <- splancs::getpoly(); colnames(syn.gates) <- params2  # draw gate
            polygon(syn.gates, lwd=2,  border="red3")
            }
        syn <- x[inout(x[,params2], syn.gates),]  
        opp[row.names(syn),"pop"] <- "synecho"

    # e. Gate Prochlorococcus
    params3 <- c("norm.scatter","norm.red")
    x <- subset(opp, pop=="unknown")
    if(gating){
            par(mfrow=c(1,1), pty="s")
            plot.cytogram(x, params3[1], params3[2], main="Gate Prochlorococcus")
            pro.gates <- splancs::getpoly(); colnames(pro.gates) <- params3 # draw gate
            polygon(pro.gates, lwd=2,  border="red3")
            }
        pro <- x[inout(x[,params3], pro.gates),] 
        opp[row.names(pro),"pop"] <- "prochloro"

    # e. Gate Picoeukaryotes
    params4 <- c("norm.scatter","norm.red")
    x <- subset(opp, pop=="unknown")
    if(gating){
            par(mfrow=c(1,1), pty="s")
            plot.cytogram(x, params4[1], params4[2], main="Gate Picoeukaryotes")
            pico.gates <- splancs::getpoly(); colnames(pico.gates) <- params4 # draw gate
            polygon(pico.gates, lwd=2,  border="red3")
            }
        pico <- x[inout(x[,params4], pico.gates),] 
        opp[row.names(pico),"pop"] <- "picoeuk"

    # f. Save gating
    gates <- list(b.gates, syn.gates, pro.gates, pico.gates) 
    save(gates, file=paste0("gating/",this_file,".RData"))
   
    ### 4. SAVE PLOT
    ################
    png(paste0("gating/",this_file,".png"),width=12, height=6, unit="in", res=200)
        par(mfrow=c(1,2), pty="s", cex=1.2, oma=c(0,0,1,0), mar=c(5,5,1,1))
        plot.vct.cytogram(opp, "norm.scatter","norm.red", ylab="red\n (normalized to beads)", xlab="scatter\n (normalized to beads)")
        plot.vct.cytogram(opp, "norm.scatter","norm.orange", ylab="orange\n (normalized to beads)", xlab="scatter\n (normalized to beads)")
        mtext(paste(basename(this_file)), outer=TRUE)
    dev.off()

    ### 5. SUMMARY 
    ##############

    stat.table <- NULL
    for(population in unique(opp$pop)){
        #print(i)
   
        p <- subset(opp, pop == population)
        n <- nrow(p)
    
        if(n ==0) {
            scatter <- 0
            red <- 0
        }else{
            scatter <- round(median(p$norm.scatter))
            red <- round(median(p$norm.red))
            orange <- round(median(p$norm.orange))
        }
        var <- cbind(population,n,scatter,red,orange)
        stat.table <- rbind(stat.table, var)
    }


    table <- data.frame(cbind(stat.table, file=basename(this_file)))
    summary.table <- rbind(summary.table, table)

}


write.csv(summary.table,file="gating/summary.csv", row.names=FALSE)
