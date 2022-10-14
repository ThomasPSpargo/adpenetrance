#####
#Helper function for ADPenetrance (https://github.com/ThomasPSpargo/adpenetrance).
#Purpose is to calculate risk of disease among people not harbouring a given variant (denoted M)
#Calculation based on lifetime risk of disease in population (PA), and variant frequency among people affected (MA).
#####

## INPUT:
#PA must be given directly
#MA can be specified directly or alternatively calculated as a weighted sum of variant frequencies MF and MS, with weighting factor PF.
#MF = variant frequency in familial state
#MS = variant frequency in sporadic state
#PF = familial frequency in population, where 1-PF is the sporadic frequency

## OUTPUT:
# A numeric indicating the lifetime risk of disease for people not harbouring the variant 'M'

getResidualRisk <- function(PA, MA=NA, MF=NA, MS=NA, PF=NA){
  
  #Either derive PMgA if required inputs are given, or proceed as specified in the input
  if(!is.na(MF) && !is.na(MS) && !is.na(PF)){
    message("Deriving total frequency of variant in people affected (MA) based on all of MF, MS, PF")
    #PSgA = 1-PF #compute probability sporadic given no variant
    MA = MF*PF+MS*(1-PF) #Prob variant given affected
  } else if (!is.na(MA)) {
    message("Proceeding with MA as defined by user") #do nothing
  } else {
    stop("Please specify MA directly, or all of MF, MS, PF so that this can be derived")
  }
  g= PA*(1-MA)
  
  return(g)
}