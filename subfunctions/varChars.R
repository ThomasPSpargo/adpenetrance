#####
# Specify varChars function, which defines a matrix indicating disease onset characteristics of people with and without a variant
# Output is to be passed to the affAtAge/genFamily functions
# Written to prepare simulated data for use in validating the ADPenetrance (https://github.com/ThomasPSpargo/adpenetrance/) approach to calculate genetic penetrance
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk) 
# Please get in touch if you have any issues.
####

## INPUT:
#f              = max (lifetime) penetrance for people with variant
#g              = max (lifetime) risk for people without variant
#numsteps       = final time interval to consider (between time 0:numsteps)
#onsetRateDiff  = time-scaling factor for group f vs group g

## OUTPUT:
# A matrix of class var.Char.matrix. The matrix contains 5 columns:
# step              #Sampling timepoint between 0 and numsteps where 0 is the last age before the disease becomes penetrant in any person
# f_propAff         #Proportion of all people harbouring a variant of lifetime penetrance f who will be affected at time point step
# g_propAff         #Proportion of all people with lifetime disease risk g (since varint f is absent) who will be affected at time point step
# f_step            #Probability of being affected by time step if harbouring f
# g_step            #Probability of being affected by time step if having lifeitme risk g (i.e. not harbouring f)

# Example of use at the bottom of script

varChars <- function(f,g,numsteps,onsetRateDiff){
  
  #Scale propor differently according to which population has compressed onset
  base_propor<- c(0,(1:numsteps/numsteps)) #Time (assuming years) from first age affected
  if(onsetRateDiff<1){
    flipped <- TRUE
    onsetRateDiff <- 1/onsetRateDiff
  } else {
    flipped <- FALSE
  }
  
  faster_propor<- base_propor*onsetRateDiff
  if(any(faster_propor>1)){
    faster_propor[which(faster_propor>=1)][1] <- 1 #Restrict the first value to the max
    faster_propor[which(faster_propor>1)]    <- 0  #set remainder to no risk, since there is none additional
  }
  
  if(onsetRateDiff>1 && flipped == FALSE){
    #Assuming f has more compressed onset density
    g_propor<- base_propor
    f_propor<- faster_propor
  } else if(flipped){
    #Assuming g has more compressed onset density
    g_propor <- faster_propor
    f_propor <- base_propor
  } else {
    #Equal onset proportions
    g_propor <- f_propor <- base_propor
  }
  
  
  f_prob<-f*f_propor            #Probability of being affected by given time, according to proportion of total penetrance
  g_prob<-g*g_propor            #repeat for people with risk g
  
  var.Char<- cbind(step=0:numsteps,    #sampling timepoint
                   f_propAff=f_propor, #Proportion of all people with f who will be affected at time point step
                   g_propAff=g_propor, #Proportion of all people with g who will be affected at time point step
                   f_step=f_prob,      #probability of being affected by age step if harbouring f
                   g_step=g_prob)      #probability of being affected by age step if harbouring g
  
  #### Set class for var.Char
  class(var.Char) <- c(class(var.Char),"var.Char.matrix")
  
  return(var.Char)
}

############ END FUNCTION SPECIFICATION ############

# #Quick functionality test
# 
# #Generate var.Char
# f=0.75
# g=0.00
# numsteps=10
# onsetRateDiff=1
# var.Char<- varChars(f,g,numsteps,onsetRateDiff)
