#####
# Script for running a simulation study based on the simADPenetrance function
# This function and associated scripts were written to validate the ADPenetrance (https://github.com/ThomasPSpargo/adpenetrance/) approach to calculate genetic penetrance
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk) 
# Please get in touch if you have any issues.
####

# This specific simulation tests a scenario where the only families in which the variant f occurs have been sampled.
# Therefore disease state rates are expected to vary over time, affecting accuracy in estimating lifetime penetrance according to age of sampling.

#Load necessary packages, and ADPenetrance functions
library(plyr)
library(ggplot2)
library(reshape2)
source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/adpenetrance_function.R") #Main function and subfunctions
#Load subfunctions functions for simulating families, via lapply
simFunctions_root<- ("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/subfunctions/")
simFunctions <- c("affAtAge.R","bindRows.R","genFamily.R","varChars.R")
invisible(lapply(paste0(root,simFunctions),source))
rm(simFunctions_root,simFunctions) #Clean temp objects from environment

start <- Sys.time()
cat("Start time:",as.character(start),"\n")

simResult <- simADPenetrance(onsetRateDiff = 1,
                             f = c(0.25,0.50,0.75,1),
                             g=0,
                             numfamilies_var = 100000,
                             numsteps = 10,
                             f_cohort_only=TRUE)

outfile = file.path(".","TimeSimulation1_onesample.pdf") #Set outpath
plotHeight<-5*40  #Set plot height; 5 state combinations modelled

#Save the output plot, and respective data
write.table(simResult$results,file=gsub(".pdf",".tsv",outfile),sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE) 
ggsave(filename=outfile,plot=simResult$ggfigure,device='pdf',units="mm",width=210,height=plotHeight,dpi=300)


warnings()

end <- Sys.time()
cat("Endtime:",as.character(end),"\n")
cat("Time taken: ",as.character(end-start),"\n")
