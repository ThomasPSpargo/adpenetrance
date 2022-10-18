#####
# Script for running a simulation study based on the simADPenetrance function
# This function and associated scripts were written to validate the ADPenetrance (https://github.com/ThomasPSpargo/adpenetrance/) approach to calculate genetic penetrance
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk) 
# Please get in touch if you have any issues.
####

# This specific simulation tests a scenario where the people with and without variant f do not have equal onset variability.
# onsetRateDiff <1 and therefore the disease onset window is more compressed for people without variant f which may affect penetrance estimation
# [this simulation is approx. the inverse of time simulation 3]

#Load necessary packages, and ADPenetrance functions
library(plyr)
library(ggplot2)
library(reshape2)
source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/adpenetrance_function.R") #Main function and subfunctions
source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/simADPenetrance.R") #main simulation function
#Load subfunctions functions for simulating families, via lapply
simFunctions_root<- ("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/subfunctions/")
simFunctions <- c("affAtAge.R","bindRows.R","genFamily.R","varChars.R")
invisible(lapply(paste0(simFunctions_root,simFunctions),source))
rm(simFunctions_root,simFunctions) #Clean temp objects from environment

start <- Sys.time()
cat("Start time:",as.character(start),"\n")

simResult <- simADPenetrance(onsetRateDiff = 0.77,
                             f = c(0.25,0.50,0.75,1),
                             g=0,
                             numfamilies_var = 100000,
                             numsteps = 10,
                             f_cohort_only=FALSE)

outfile <- file.path(".","TimeSimulation4_0pt77onset.pdf") #Set outpath
plotHeight <- 5*40  #Set plot height; 5 state combinations modelled

#Save the output plot, and respective data
ggsave(filename=outfile,plot=simResult$ggfigure,device='pdf',units="mm",width=210,height=plotHeight,dpi=300)
write.table(simResult$results,file=gsub(".pdf",".tsv",outfile),sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE) 


warnings()

end <- Sys.time()
cat("Endtime:",as.character(end),"\n")
cat("Time taken: ",as.character(end-start),"\n")
