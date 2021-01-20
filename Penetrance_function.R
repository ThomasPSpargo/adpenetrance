#Version updated 19/01/2021


#NOTE:
#Before using this function, please first consult the README.md documentation in the associated github repository:
#https://github.com/thomaspspargo/penetrance-calculator


#--------------------------------------------#
#       Penetrance calculator function       #
#--------------------------------------------#


#Specify function
Pen_calculator <- function(N, MF=0, MS=0, MA=0, MU=0, PA=0, PF=0, MF_SE=0, MS_SE=0, MU_SE=0, MA_SE=0, Zout=1.96){
  
  #Construct penetrance value sequence for lookup table. Values generated are between 0 and 1 at increments of .0001
  f = seq(0,1, by=.0001)
  
  #--------------------------------------------#
  #   IF ELSE statements for input variables   #
  #--------------------------------------------#
  #These statements evaluate which disease states are specified in funcition to determine how each input variable is stored and calculation is performed.
  #All acceptable combinations are defined unacceptable selections will be stored as NA and the output will not be valid
  
  
  #Mutation frequency 1 - stored as state X
  
  if(MF>0 && MS>0 && MU>0 && MA==0) { #If familial sporadic, and unaffected specified, and affected not specified
    X <- MF #Familial state
  } else if(MF>0 && MS>0 && MU==0 && MA==0) { #If familial and sporadic specified, and unaffected and affected not specified
    X <- MF #Familial state
  } else if(MF>0 && MS==0 && MU>0 && MA==0) { #If familial and population specified, and sporadic and affected not specified
    X <- MF #Familial state
  } else if(MF==0 && MS>0 && MU>0 && MA==0) { #If sporadic and unaffected specified, and familial and affected not specified
    X <- MS  #Sporadic state
  } else if(MF==0 && MS==0 && MU>0 && MA>0) { #If unaffected and affected specified, familial and sporadic not specified
    X <- MA #Affected state
  } else {
    X <- NA #If none of the above, value is NA 
  }
  
  
  #Weighting variable for mutation frequency 1
  if(MF>0 && MS>0 && MU>0 && MA==0) {
    WeiX <- PF*PA
  }else if(MF>0 && MS>0 && MU==0 && MA==0) {
    WeiX <- PF
  } else if(MF>0 && MS==0 && MU>0 && MA==0) {
    WeiX <- PF*PA
  } else if(MF==0 && MS>0 && MU>0 && MA==0) {
    WeiX <- (1-PF)*PA
  } else if(MF==0 && MS==0 && MU>0 && MA>0) {
    WeiX <- PA
  } else {
    WeiX <- NA
  }
  
  #Mutation frequency 2 - stored as state Y
  
  if(MF>0 && MS>0 && MU>0 && MA==0) {
    Y <- MS #Sporadic state
  } else if(MF>0 && MS>0 && MU==0 && MA==0) {
    Y <- MS #Sporadic state
  } else if(MF>0 && MS==0 && MU>0 && MA==0) {
    Y <- MU #Unaffected state
  } else if(MF==0 && MS>0 && MU>0 && MA==0) {
    Y <- MU #Unaffected state
  } else if(MF==0 && MS==0 && MU>0 && MA>0) {
    Y <- MU #Unaffected state
  } else {
    Y <- NA
  }
  
  #Weighting variable for mutation frequency 2
  if(MF>0 && MS>0 && MU>0 && MA==0) {
    WeiY = (1-PF)*PA
  } else if(MF>0 && MS>0 && MU==0 && MA==0) {
    WeiY <- (1-PF)
  } else if(MF>0 && MS==0 && MU>0 && MA==0) {
    WeiY <- 1-PA
  } else if(MF==0 && MS>0 && MU>0 && MA==0) {
    WeiY <- 1-PA
  } else if(MF==0 && MS==0 && MU>0 && MA>0) {
    WeiY <- 1-PA
  } else {
    WeiY <- NA
  }
  
  
  #Mutation frequency 3 - stored as state Z
  #Will only be MU
  
  if(MF>0 && MS>0 && MU>0 && MA==0) { #If familial sporadic, and unaffected specified, and affected not specified
    Z <- MU #Unaffected
  } else {
    Z <- NA
  }
  
  #Weighting variable for mutation frequency 3
  if(MF>0 && MS>0 && MU>0 && MA==0) {
    WeiZ <- 1-PA
  } else {
    WeiZ <- NA
  }
  
  
  #Additional IF ELSE statements for error propagation variables
  
  
  #Standard error for mutation frequency 1
  if(MF>0 && MS>0 && MU>0 && MA==0 && MF_SE>0) {
    XSE <- MF_SE
  } else if(MF>0 && MS>0 && MU==0 && MA==0 && MF_SE>0) {
    XSE <- MF_SE
  } else if(MF>0 && MS==0 && MU>0 && MA==0 && MF_SE>0) {
    XSE <- MF_SE
  } else if(MF==0 && MS>0 && MU>0 && MA==0 && MS_SE>0) {
    XSE <- MS_SE
  } else if(MF==0 && MS==0 && MU>0 && MA>0 && MA_SE>0) {
    XSE <- MA_SE
  } else {
    XSE <- NA
  }
  
  #Standard error for mutation frequency 2
  if(MF>0 && MS>0 && MU>0 && MA==0 && MS_SE>0) {
    YSE <- MS_SE
  } else if(MF>0 && MS>0 && MU==0 && MA==0 && MS_SE>0) {
    YSE <- MS_SE
  } else if(MF>0 && MS==0 && MU>0 && MA==0 && MU_SE>0) {
    YSE <- MU_SE
  } else if(MF==0 && MS>0 && MU>0 && MA==0 && MU_SE>0) {
    YSE <- MU_SE
  } else if(MF==0 && MS==0 && MU>0 && MA>0 && MU_SE>0) {
    YSE <- MU_SE
  } else {
    YSE <- NA
  }
  
  #Standard error for mutation frequency 3
  if(MF>0 && MS>0 && MU>0 && MA==0 && MU_SE>0) {
    ZSE <- MU_SE
  } else {
    ZSE <- NA
  }
  
  
  
  
  
  
  
  
  ###Construct a string labels for disease state of mutation frequency 1
  #This is used in constructing the output tables
  if(MF>0 && MS>0 && MU>0 && MA==0) {
    LabelX = "familial"
  } else if(MF>0 && MS>0 && MU==0 && MA==0) {
    LabelX = "familial"
  } else if(MF>0 && MS==0 && MU>0 && MA==0) {
    LabelX = "familial"
  } else if(MF==0 && MS>0 && MU>0 && MA==0) {
    LabelX = "sporadic"
  } else if(MF==0 && MS==0 && MU>0 && MA>0) {
    LabelX = "affected"
  } else {
    LabelX <- NA
  }
  
  
  #-------------------------------#
  #       Main calculations       #
  #-------------------------------#
  
  #Calculate the observed probability of disease state X (ObsProbX) from mutation frequency data and weighting factors
  #This is a weighted proportion calculation
  
  if(MF>0 && MS>0 && MU>0 && MA==0) { #Calculation if data specified for familial, sporadic, and unaffected
    ObsProbX = (X*WeiX)/((X*WeiX)+(Y*WeiY)+(Z*WeiZ)) #Weighted proportion of disease states
  } else {  #Calculation for all other valid combinations of data 
    ObsProbX = (X*WeiX)/((X*WeiX)+(Y*WeiY)) #Weighted proportion of disease states
  }
  
  
  
  #Optionally calculate standard error of ObsProbX - perform error propagation
  #Applies the calculus approach to error propagation from Hughes and Hase (2010)
  #Step 1: calculate partial derivatives of the ObsProbX equation according to each variable
  #Step 2: Apply Hughes & Hase equation
  
  if(!is.na(XSE) && !is.na(YSE)) { #Do this step if error terms are given
    
    if(MF>0 && MS>0 && MU>0 && MA==0 && !is.na(ZSE)) { #Calculate with respect to variables X Z and Z if data specified for familial, sporadic and unaffected
      
      DifX = (WeiX*(Y*WeiY+Z*WeiZ))/((X*WeiX+Y*WeiY+Z*WeiZ)^2) #Partial derivative of the ObsProbX calculation with respect to X
      
      DifY = -((X*WeiX*WeiY)/((X*WeiX+Y*WeiY+Z*WeiZ)^2)) #Partial derivative of the ObsProbX calculation with respect to Y
      
      DifZ = -((X*WeiX*WeiZ)/((X*WeiX+Y*WeiY+Z*WeiZ)^2)) #Partial derivative of the ObsProbX calculation with respect to Z
      
      ObsProbXSE <- sqrt(DifX^2*XSE^2 + DifY^2*YSE^2 + DifZ^2*ZSE^2) #Use partial derivatives and Std errors to propagate error for observed probability of disease state X
      
      
    } else { #In all other instances, calculate with respect to only X and Y
      
      DifX <- (WeiX*WeiY*Y)/((WeiX*X+WeiY*Y)^2) #Partial derivative of the ObsProbX calculation with respect to X
      DifY <- -((WeiX*WeiY*X)/((WeiY*Y+WeiX*X)^2)) #Partial derivative of the ObsProbX calculation with respect to Y
      
      
      ObsProbXSE <- sqrt(DifX^2*XSE^2 + DifY^2*YSE^2) #Partial derivatives*Std errors+...= SE in ObsProbX
      
    }
    
  } else {ObsProbXSE <- as.numeric(NA)} #If error terms not given, skip this operation
  
  
  #If error has been propagated: Convert ObsProbXSE into confidence intervals through z-score conversion
  if (!is.na(ObsProbXSE)) {
    lowerCI = ObsProbX-(Zout*ObsProbXSE)   #Lower interval
    upperCI = ObsProbX+(Zout*ObsProbXSE)   #Upper interval
    
  } else {
    lowerCI=as.numeric(NA)
    upperCI=as.numeric(NA)
  }
  
  
  
  #Begin "operation 2", 
  #Apply Al-Chalabi & Lewis (2011) equations to calculate P(unaffected), P(sporadic), P(familial) 
  #Do this at the defined specified sibship size for a sequence of f (penetrance) values
  #Store output of each in data.frame with respective f value
  dis_states <- data.frame(
    f,
    Punasc = (1-f)*(1-f/2)^N, #Probability unaffected
    Pspor  = f*(1-f/2)^N+N*(f/2)*(1-f/2)^(N-1)*(1-f), #Probability sporadic/simplex
    Pfam   = 1-((1-f/2)^N+N*(f/2)*((1-f/2)^(N-1))*(1-f)) #Probability familial (= 1 - Punasc - Pspor)
  )
  
  #Build lookup table recording EXPECTED value of ObsProbX at each f value.
  #A series of LookupX values, each respective to a value of f, are calculated as a proportion of those states defined in Operation 1.
  #Calculation: Probability of state X / sum(probabilities of states defined)
  #Formating of if else statements as before
  
  if(MF>0 && MS>0 && MU>0 && MA==0) {
    LookupX = dis_states[,"Pfam"]/(dis_states[,"Pfam"]+dis_states[,"Pspor"]+dis_states[,"Punasc"]) #F/(F+S+U) = F/1
  } else if(MF>0 && MS>0 && MU==0 && MA==0) { 
    LookupX = dis_states[,"Pfam"]/(dis_states[,"Pfam"]+dis_states[,"Pspor"]) #F/(F+S)
  } else if(MF>0 && MS==0 && MU>0 && MA==0) {
    LookupX = dis_states[,"Pfam"]/(dis_states[,"Pfam"]+dis_states[,"Punasc"]) #F/(F+U)
  } else if(MF==0 && MS>0 && MU>0 && MA==0) {
    LookupX = dis_states[,"Pspor"]/(dis_states[,"Pspor"]+dis_states[,"Punasc"]) #S/(S+U)
  } else if(MF==0 && MS==0 && MU>0 && MA>0) {
    LookupX = (dis_states[,"Pfam"]+dis_states[,"Pspor"])/(dis_states[,"Pfam"]+dis_states[,"Pspor"]+dis_states[,"Punasc"]) #(F+S)/(F+S+U) = (F+S)/1 = A/1
  } else {
    LookupX = as.numeric(NA)
  }
  
  #Build a two column Lookup table of LookupX values with respective f values
  LookupTable <- cbind(LookupX, f)
  
  
  #Begin operation 3
  #Compare ObsProbX to LookupTable - retrieve loci of LookupX with nearest value to ObsProbX (and intervals) 
  #The values should be comparable - large disparity suggests ObsProbX exceeds or is lesser than rate expected at penetrance= 1 or 0
  
  #Prepare matrix to store the locus of the nearest LookupX value
  loci <- matrix(NA,ncol=3,nrow=1)
  colnames(loci) <- c("Lower CI", "Estimate", "Upper CI")
  
  #Matching ObsProbX to closest LookupX
  #Store locus retrieved in 'loci' object defined above
  if(!is.na(ObsProbX)) {
    loci[,"Estimate"] <- which(abs(LookupTable[,1]-as.vector(ObsProbX))==min(abs(LookupTable[,1]-as.vector(ObsProbX)),na.rm=T)) #Locus for the estimate
  }
  
  #If error propagation has been included:
  #Repeat process at disease state rates observed at confidence interval bounds
  if(!is.na(ObsProbXSE)) {
    loci[,"Lower CI"] <- which(abs(LookupTable[,1]-as.vector(lowerCI))==min(abs(LookupTable[,1]-as.vector(lowerCI)),na.rm=T)) #Locus for the lower bound
    loci[,"Upper CI"] <- which(abs(LookupTable[,1]-as.vector(upperCI))==min(abs(LookupTable[,1]-as.vector(upperCI)),na.rm=T)) #Locus for the upper bound
  }
  
  
  #Build matrix 1 for output:
  #Row 1: ObsProbX (with or without CI and standard error)
  #Row 2: Value of LookupX closest to ObsProbX (with or without CI) - identified at row stored in 'loci' object
  #Row 3: Penetrance estimate associated with LookupX (with or without CI)
  #This will be table 1
  
  if(!is.na(ObsProbXSE)) { #Include confidence intervals if propagation performed
    
    Output = matrix(c(lowerCI, ObsProbX, upperCI, ObsProbXSE, #Estimate for X rate, with CI and SE given
                      LookupTable[,1][loci[,"Lower CI"]], LookupTable[,1][loci[,"Estimate"]], LookupTable[,1][loci[,"Upper CI"]], NA, #LookupX nearest estimate and CI
                      LookupTable[,2][loci[,"Lower CI"]], LookupTable[,2][loci[,"Estimate"]], LookupTable[,2][loci[,"Upper CI"]], NA), #Corresponding penetrance estimate
                    byrow =T, ncol=4, nrow=3)
    
    colnames(Output) <- c("Lower CI","Estimate", "Upper CI", "Standard error") #Columnames
    
  } else {   #Create without confidence intervals
    
    Output = matrix(c(ObsProbX, #Estimate for X rate
                      LookupTable[,1][loci[,"Estimate"]], #LookupX nearest estimate
                      LookupTable[,2][loci[,"Estimate"]]), #Corresponding penetrance estimate
                    ncol=1, nrow=3)
    
    colnames(Output) <- c("Estimate") #Columnames
  }
  
  #Construct Row names for the first two rows
  #LabelX is the name of the state modelled as state X
  R1 <- paste("Observed", LabelX, "rate")
  R2 <- paste("Expected", LabelX, "rate")
  
  #Assign row names
  rownames(Output) <- c(R1, R2, "Penetrance")
  
  
  #Build matrix 2 for output:
  #Calculate probabilities of being familial, sporadic or unaffected expected at the penetrance value obtained
  graphloci<- which(abs(dis_states[,"f"]-Output["Penetrance","Estimate"])==min(abs(dis_states[,"f"]-Output["Penetrance","Estimate"]))) #Find penetrance locus for estimate
  Famprob <- dis_states[graphloci[1],"Pfam"] #Find P(familial) at locus
  Sporprob <- dis_states[graphloci[1],"Pspor"]  #Find P(sporadic) at locus
  Unascprob <- dis_states[graphloci[1],"Punasc"]  #Find P(unaffected) at locus
  
  if(!is.na(ObsProbXSE)) { #If error terms given, repeat for penetrance values at CI bounds
    #Repeat for confidence intervals
    lowloci<- which(abs(dis_states[,"f"]-Output["Penetrance","Lower CI"])==min(abs(dis_states[,"f"]-Output["Penetrance","Lower CI"]))) #Find penetrance locus for lower interval
    Famlow <- dis_states[lowloci[1],"Pfam"]
    Sporlow <- dis_states[lowloci[1],"Pspor"]
    Unasclow <- dis_states[lowloci[1],"Punasc"]
    
    uploci<- which(abs(dis_states[,"f"]-Output["Penetrance","Upper CI"])==min(abs(dis_states[,"f"]-Output["Penetrance","Upper CI"]))) #Find penetrance locus for upper interval
    Famup <- dis_states[uploci[1],"Pfam"]
    Sporup <- dis_states[uploci[1],"Pspor"]
    Unascup <- dis_states[uploci[1],"Punasc"]
    
    #Store in 3x3 matrix
    Probs<- matrix(c(Famlow, Famprob, Famup,  #Familial estimates
                     Sporlow, Sporprob, Sporup, #Sporadic estimates
                     Unasclow, Unascprob, Unascup), #Unaffected estimates
                   nrow=3,ncol=3, byrow=T)
    
    colnames(Probs) <- c("Lower CI","Estimate", "Upper CI") #Assign column names
    
    
  } else { #if error terms not given, construct table just at estimates
    #Store in matrix
    Probs<- matrix(c(Famprob,Sporprob,Unascprob),nrow=3,ncol=1)
    
    colnames(Probs) <- c("Estimate") #Assign column names
  }
  
  rownames(Probs) <- c("Familial", "Sporadic","Unaffected") #Assign row names
  
  
  
  return(list(Calculated=Output, Probabilities=Probs))
  
}

