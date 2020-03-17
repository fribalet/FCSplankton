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

plot.cytogram <- function (evtopp, para.x = "scatter", para.y = "red", ...){
    par(pty = "s")
    plot(evtopp[, c(para.x, para.y)], pch = 16, cex = 0.3,col = densCols(log10(evtopp[, c(para.x, para.y)]),
          colramp = viridis::viridis), log='xy',...)
    }

plot.vct.cytogram <- function (opp, para.x = "scatter", para.y = "red", ...){
    group.colors <- c(unknown="grey", beads="red3", prochloro=viridis::viridis(4)[1],synecho=viridis::viridis(4)[2],picoeuk=viridis::viridis(4)[3], croco=viridis::viridis(4)[4])
    opp$pop <- factor(opp$pop, levels = names(group.colors))

    par(pty = "s")        
    plot(opp[, c(para.x, para.y)], pch = 16, cex = 0.3, col = as.numeric(as.factor(opp$pop)), log='xy', ...)
    legend("topleft", legend = (unique(opp$pop)), col = unique(as.numeric(as.factor(opp$pop))), pch = 16, pt.cex = 0.6, bty = "n")
    }


### 3. Gating populations of interest
#####################################
# set the path to the FCS files
setwd("/Volumes/GoogleDrive/My\ Drive/_PROJECT_")

# load list of FCS files
file.list <- dir(".", pattern = ".fcs$", recursive=F)

# load metadata
meta <- read.csv(file = paste0('meta.csv'), header = T, sep = ',')
ok_meta <- subset(meta, Flag == 0)      # only the unflagged files


# a. Initialization
this_file <- file.list[1]
opp <- caroline::tab2df(exprs(read.FCS(this_file, transformation = T, emptyValue=F))) # create a flowFrame WITH transformation (for log-amplified data)

    # Rename PMTs
    print(colnames(opp))
    colnames(opp)[i] <- "scatter" # replace i by the column indice of FSC
    colnames(opp)[i] <- "red" # replace i by the column indice of 692
    colnames(opp)[i] <- "orange" # replace i by the column indice of 580
    # colnames(opp)[i] <- "green" # replace i by the column indice of 530


gating <- TRUE

for (this_file in file.list){

    # this_file <- file.list[3]
    print(paste("gating", this_file)
    opp <- caroline::tab2df(exprs(read.FCS(this_file, transformation = T, emptyValue=F))) # read FCS file (Transformation = TRUE for log-amplified data)
    opp$pop <- 0 # add a colum 'pop' to the table


    # Load gating if exist
    if(!gating){    
            load(file=paste0("gating/",this_file, ".RData")
            b.gates <- gates[[1]]
            syn.gates <- gates[[2]]
            pro.gates <- gates[[3]]
            pico.gates <- gates[[4]]
           }

    ### Beads Normalization
    # a. Gate beads
    params1 <- c("scatter","orange")
    if(gating){
            plot.cytogram(opp[,params1], main="Gate Beads")
            b.gates <- splancs::getpoly(); colnames(b.gates) <- params1 # draw gate
            polygon(b.gates, lwd=2,  border='red3')
            }
    # b. Filter beads
    beads <- opp[inout(opp[,params1], b.gates),] 

    # c. Normalization
        opp$norm.scatter <- opp$scatter / median(beads$scatter)
        opp$norm.orange <- opp$orange / median(beads$orange)
        opp$norm.red <- opp$red / median(beads$red)
        norm.opp <- opp[!inout(opp[,params1], b.gates),] # exclude beads

    # d. Gate Synecho
    params2 <- c("norm.scatter","norm.orange")
    x <- subset(norm.opp, pop==0)
    if(gating){
            plot.cytogram(x[,params2], main="Gate Synechococcus")
            syn.gates <- splancs::getpoly(); colnames(syn.gates) <- params2  # draw gate
            polygon(syn.gates, lwd=2,  border='red3')
            }
        syn <- x[inout(x[,params2], syn.gates),]  
        norm.opp[row.names(syn),'pop'] <- "synecho"

    # e. Gate Prochlorococcus
    params3 <- c("norm.scatter","norm.red")
    x <- subset(norm.opp, pop==0)
    if(gating){
            plot.cytogram(x[,params3], main="Gate Prochlorococcus")
            pro.gates <- splancs::getpoly(); colnames(pro.gates) <- params3 # draw gate
            polygon(pro.gates, lwd=2,  border='red3')
            }
        pro <- x[inout(x[,params3], pro.gates),] 
        norm.opp[row.names(pro),'pop'] <- "prochloro"

    # e. Gate Picoeukaryotes
    params4 <- c("norm.scatter","norm.red")
    x <- subset(norm.opp, pop==0)
    if(gating){
            plot.cytogram(x[,params3], main="Gate Picoeukaryotes")
            pico.gates <- splancs::getpoly(); colnames(pico.gates) <- params4  # draw gate
            polygon(pico.gates, lwd=2,  border='red3')
            }
        pico <- x[inout(x[,params4], pico.gates),] 
        norm.opp[row.names(pico),'pop'] <- "picoeuk"

    # f. Save gating
    gates <- list(b.gates, syn.gates, pro.gates, pico.gates) 
    save(gates, file=paste0("gating/",this_file,".RData"))
   
    ### 4. SAVE PLOT
    ################
    png(paste0("gating/",this_file,".png"),width=12, height=9, unit='in', res=100)
        par(mfrow=c(1,2))
        plot.vct.cytogram(norm.opp, "scatter","red")
        plot.vct.cytogram(norm.opp, "scatter","orange")
    dev.off()

    ### 5. SUMMARY 
    ##############

    stat.table <- NULL
    for(population in unique(norm.opp$pop)){
        #print(i)
        if(population == 0) next
    
        p <- subset(norm.opp, pop == population)
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


    table <- data.frame(cbind(stat.table, file=basename(file)))
    summary.table <- rbind(summary.table, table)

}


write.csv(summary.table,file="gating/summary.csv", row.names=FALSE)
