#####
#Helper function for ADPenetrance (https://github.com/ThomasPSpargo/adpenetrance).
#Purpose is to calculate risk of disease among people not harbouring a given variant (denoted M)
#Calculation based on lifetime risk of disease in population (PA), and variant frequency among people affected (MA).
#####

## INPUT:
#PA must be given directly
#MA can be specified directly or alternatively calculated as a weighted sum of variant frequencies MF and MS, with weighting factor PF.
#PM can be specified directly or alternatively calculated as a weighted sum of variant frequencies MA and MU, with weighting factor PA.
#MF = variant frequency in familial state
#MS = variant frequency in sporadic state
#MU = variant frequency in unaffected state
#PF = familial frequency in population, where 1-PF is the sporadic frequency

## OUTPUT:
# A numeric indicating the lifetime risk of disease for people not harbouring the variant 'M'

getResidualRisk <- function(PA, MA=NA, PM=NA, MF=NA, MS=NA, MU=NA, PF=NA){
  
  ### PA input
  if (!is.numeric(PA) || !(PA<=1 && PA>=0)) { #Return error if PA is impossible
    stop("Please specify PA as a numeric between 0 and 1.")
  }
  
  ### MA input
  if (is.na(MA) || !is.numeric(MA)) { #Unless MA specified directly, try to derive PMgA from alternative inputs
    if(!is.na(MF) && !is.na(MS) && !is.na(PF)){
      message("Deriving total frequency of variant in people affected (MA) based on all of MF, MS, PF")
      #PSgA = 1-PF #compute probability sporadic given no variant
      MA = MF*PF+MS*(1-PF) #Prob variant given affected
    } 
  }
  if (!is.numeric(MA) || !(MA<=1 && MA>=0)) { #Return error if MA is impossible 
    stop("MA is not a numeric between 0 and 1. Please specify MA directly, or all of MF, MS, PF so that this can be derived.")
  }
  
  ### PM input
  if (is.na(PM) || !is.numeric(PM)) { #Unless PM specified directly, try to derive PM from alternative inputs
    if(!is.na(MU)){
      message("Deriving total probability of having variant (PM) based on all of MA, MU, PA")
      PM=MA*PA+MU*(1-PA)
    } 
  }
  if (!is.numeric(PM) || !(PM<=1 && PM>=0)) { #Return error if PM is impossible 
    stop("PM is not a numeric between 0 and 1. Please specify PM directly, or all of MA, MU, PA so that this can be derived.")
  }
  
  ### Derive g
  g = (PA*(1-MA))/(1-PM)
  
  return(g)
}