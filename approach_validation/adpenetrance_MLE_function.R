#This script is associated with the following GitHub repository: https://github.com/thomaspspargo/adpenetrance.
#It generates the adpenetrance function for calculating genetic penetrance in autosomal dominant traits,
  #including an additional step to allow Maximum Likelihood Estimation of unadjusted penetrance estimates.
    #This is used for validation of the lookup table approach primarily employed by the method,
    #and is to be called in the adpenetrance_MLE_validation.R script contained within the associated github repository.

  #The function relies upon the sub-functions adpenetrance.errorfit and adpenetrance.unadjusted, which are downladed from GitHub when running this script

#--------------------------------------------#
#       Penetrance calculator function       #
#--------------------------------------------#

#First, if not already loaded in the environment, download the adpenetrance.unadjusted and adpenetrance.errorfit subfunctions from GitHub
if(!"adpenetrance.errorfit" %in% ls()){
  source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/subfunctions/adpenetrance_errorfit_function.R")
}

#If the function is not defined in the environment, download from github
if(!"adpenetrance.unadjusted" %in% ls()){
  source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/subfunctions/adpenetrance_unadjusted_function.R")
}


#Second, specify main adpenetrance function
adpenetrance.MLE <- function(N, MF=0, MS=0, MA=0, MU=0, PA=0, PF=0, MF_SE=0, MS_SE=0, MU_SE=0, MA_SE=0, Zout=1.96,
                             RX=0, RX_SE=0, states="none",define_sibstructure=NULL,include_MLE=TRUE){
  
  
  #--------------------------------------------#
  #   IF ELSE statements for input variables   #
  #--------------------------------------------#
  #First, define the states modelled
  
  #Perform all these calculations only if RX is not given directly:
  
  #Generate value for states term within the function - necessary to pass through states to adpenetrance.errorfit
  
  if(states=="none"){ 
    if(MF>0 && MS>0 && MU>0 && MA==0) {
      states="fsu"
    } else if(MF>0 && MS>0 && MU==0 && MA==0) {
      states="fs"
    } else if(MF>0 && MS==0 && MU>0 && MA==0) {
      states="fu"
    } else if(MF==0 && MS>0 && MU>0 && MA==0) {
      states="su"
    } else if(MF==0 && MS==0 && MU>0 && MA>0) {
      states="au"
    } else {
      stop("No valid disease state combination has been defined. Variant frequency or RX estimates should be defined for any two or three of the familial, sporadic, and unaffected states or the affected and unaffected states. Please check that variant frequency estimates have been defined for a valid combination of states or that the 'RX' and 'states' arguments are properly defined.")
    }
  } else if(!any(c('fsu','fu','fs','au','su') %in% states)){ #Change states value to "incorrect" if nothing has been defined
    stop("No valid disease state combination has been defined. Variant frequency or RX estimates should be defined for any two or three of the familial, sporadic, and unaffected states or the affected and unaffected states. Please check which states have been indicated in the 'states' object")
  }
  
  #The following statements evaluate how data are specified in the function to determine how to calculate the observed rate of state X.
  #All valid combinations are defined and unacceptable selections will be stored as NA and the output will not be valid
  if(RX==0){ 
    if(MF>0 && MS>0 && MU>0 && MA==0) { #If familial, sporadic, unaffected specified, and affected not specified
      X <- MF           #Familial state is state X
      WeiX <- PF*PA     #Weighting variable
      
      Y <- MS           #Sporadic state is state Y
      WeiY = (1-PF)*PA  #Weighting variable
      
      Z <- MU           #Unaffected is state Z (only used in this formulation - after this only X and Y need specification)
      WeiZ <- 1-PA      #Weighting variable
      
      if(MF_SE>0 && MS_SE>0 && MU_SE>0) { #If error terms are given for all appropriate variables, define error terms
        XSE <- MF_SE
        YSE <- MS_SE
        ZSE <- MU_SE
      } else if (MF_SE>0 || MS_SE>0 || MU_SE>0) { #If error terms are given for some but not all variables, do not define error terms and return warning message
        XSE <- NA
        YSE <- NA
        ZSE <- NA
        
        warning("Error terms have only been given for a subset of the variant frequencies defined. Please ensure that these are given for all included states")
      } else { 
        XSE <- NA
        YSE <- NA
        ZSE <- NA
      }
      
    } else if(MF>0 && MS>0 && MU==0 && MA==0) { #If familial and sporadic specified, and unaffected and affected not specified
      X <- MF           #Familial state is state X
      WeiX <- PF        
      
      Y <- MS           #Sporadic state is state Y
      WeiY <- 1-PF      
      
      if(MF_SE>0 && MS_SE>0) { #If error terms are given for all appropriate variables, define error terms
        XSE <- MF_SE
        YSE <- MS_SE
      } else if (MF_SE>0 || MS_SE>0) { #If error terms are given for some but not all variables, do not define error terms and return warning message
        XSE <- NA
        YSE <- NA
        
        warning("Error terms have only been given for a subset of the variant frequencies defined. Please ensure that these are given for all included states")
      } else {
        XSE <- NA
        YSE <- NA
      }
      
      
    } else if(MF>0 && MS==0 && MU>0 && MA==0) { #If familial and unaffected specified, and sporadic and affected not specified
      X <- MF         #Familial state is X
      WeiX <- PF*PA 
      
      Y <- MU         #Unaffected state is Y
      WeiY <- 1-PA
      
      if(MF_SE>0 && MU_SE>0) { #If error terms are given for all appropriate variables, define error terms
        XSE <- MF_SE
        YSE <- MU_SE
      } else if (MF_SE>0 || MU_SE>0) { #If error terms are given for some but not all variables, do not define error terms and return warning message
        XSE <- NA
        YSE <- NA
        
        warning("Error terms have only been given for a subset of the variant frequencies defined. Please ensure that these are given for all included states")
      } else {
        XSE <- NA
        YSE <- NA
      }
      
    } else if(MF==0 && MS>0 && MU>0 && MA==0) { #If sporadic and unaffected specified, and familial and affected not specified
      X <- MS         #Sporadic state is X
      WeiX <- (1-PF)*PA
      
      Y <- MU         #Unaffected state is Y
      WeiY <- 1-PA
      
      if(MS_SE>0 && MU_SE) { #If error terms are given for all appropriate variables, define error terms
        XSE <- MS_SE
        YSE <- MU_SE
      } else if (MS_SE>0 || MU_SE>0) { #If error terms are given for some but not all variables, do not define error terms and return warning message
        XSE <- NA
        YSE <- NA
        
        warning("Error terms have only been given for a subset of the variant frequencies defined. Please ensure that these are given for all included states")
      } else {
        XSE <- NA
        YSE <- NA
      }
      
    } else if(MF==0 && MS==0 && MU>0 && MA>0) { #If unaffected and affected specified, familial and sporadic not specified
      X <- MA         #Affected state is X
      WeiX <- PA
      
      Y <- MU         #Unaffected state is Y
      WeiY <- 1-PA
      
      if(MA_SE>0 && MU_SE>0) { #If error terms are given for all appropriate variables, define error terms
        XSE <- MA_SE
        YSE <- MU_SE
      } else if (MA_SE>0 || MU_SE>0) { #If error terms are given for some but not all variables, do not define error terms and return warning message
        XSE <- NA
        YSE <- NA
        
        warning("Error terms have only been given for a subset of the variant frequencies defined. Please ensure that these are given for all included states")
      } else {
        XSE <- NA
        YSE <- NA
      }
      
    } else { #If no valid conditions are satisfied
      
      stop("No valid disease state combination has been defined. Variant frequency or RX estimates should be defined for any two or three of the familial, sporadic, and unaffected states or the affected and unaffected states. Please check that variant frequency estimates have been defined for a valid combination of states or that the 'RX' and 'states' arguments are properly defined.")
      
    }
    
    #Calculate the observed probability of disease state X (ObsProbX) from state data given
    #This is a weighted proportion calculation
    if(states=="fsu") {                                 #Calculation if data specified for familial, sporadic, and unaffected
      ObsProbX = (X*WeiX)/((X*WeiX)+(Y*WeiY)+(Z*WeiZ))
    } else {                                            #Calculation for all other valid combinations of data 
      ObsProbX = (X*WeiX)/((X*WeiX)+(Y*WeiY))
    }
    
    #Perform error propagation if error terms are provided
    #Applies the calculus approach to error propagation from Hughes and Hase (2010)
    #Step 1: Calculate partial derivatives of the ObsProbX equation according to each variable
    #Step 2: Apply Hughes & Hase equation
    if(!is.na(XSE) && !is.na(YSE)) { #Do this step only if error terms are given
      
      if(states=="fsu" && !is.na(ZSE)) { #Calculate with respect to variables X, Y, and Z if data specified for familial, sporadic, and unaffected
        
        DifX = (WeiX*(Y*WeiY+Z*WeiZ))/((X*WeiX+Y*WeiY+Z*WeiZ)^2) #Partial derivative of the ObsProbX calculation with respect to X
        DifY = -((X*WeiX*WeiY)/((X*WeiX+Y*WeiY+Z*WeiZ)^2))       #Partial derivative of the ObsProbX calculation with respect to Y
        DifZ = -((X*WeiX*WeiZ)/((X*WeiX+Y*WeiY+Z*WeiZ)^2))       #Partial derivative of the ObsProbX calculation with respect to Z
        
        ObsProbXSE <- sqrt(DifX^2*XSE^2 + DifY^2*YSE^2 + DifZ^2*ZSE^2) #Use partial derivatives and Std errors to propagate error for observed probability of disease state X
        
      } else { #In all other instances, calculate with respect to only X and Y
        
        DifX <- (WeiX*WeiY*Y)/((WeiX*X+WeiY*Y)^2) #Partial derivative of the ObsProbX calculation with respect to X
        DifY <- -((WeiX*WeiY*X)/((WeiY*Y+WeiX*X)^2)) #Partial derivative of the ObsProbX calculation with respect to Y
        
        ObsProbXSE <- sqrt(DifX^2*XSE^2 + DifY^2*YSE^2) #Partial derivatives*Std errors+...= SE in ObsProbX
        
      }
    } else {ObsProbXSE <- as.numeric(NA)} #If error terms not given, skip this operation
    
    
    
    
  } else if(RX>0){ #IF RX is specified directly, the above is skipped and ObsProbX is defined directly by the user
    if(MF==0 && MS==0 && MU==0 && MA==0){ #Do not do this if the variant frequencies have also been defined
      ObsProbX = RX #Define ObsProbX directly from RX
      
      if(RX_SE>0){
        ObsProbXSE = RX_SE #Define SE in ObsProbX from the error specified
      } else {
        ObsProbXSE <- as.numeric(NA) #Define as NA if error is unspecified
      }
    } else {
      stop("Input data have been provided for at least 1 variant frequency as well as RX. Please specify either the 'RX' and the 'states' arguments, or specify variant frequency estimates for each represented disease state")
    }
  }
  
  
  #Calculate confidence intervals for ObsProbX using z-score conversion, if error is known.
  #No if statement according to whether RX is given or variant frequencies
  if (!is.na(ObsProbXSE)) {
    lowerCI = ObsProbX-(Zout*ObsProbXSE)   #Lower interval
    upperCI = ObsProbX+(Zout*ObsProbXSE)   #Upper interval
    
  } else {
    lowerCI=as.numeric(NA)
    upperCI=as.numeric(NA)
  }
  
  #Step 2:
  #Construct penetrance value sequence for lookup table. Values generated are between 0 and 1 at increments of 0.0001
  f = seq(0,1, by=0.0001)
  
  #Apply Al-Chalabi & Lewis (2011) equations to calculate P(unaffected), P(sporadic), P(familial) 
  #Do this at the defined specified sibship size for all values of f- store in matrix of length(f)
  dis_states <- matrix(c(
    f,
    Punasc = (1-f)*(1-f/2)^N, #Probability unaffected
    Pspor  = f*(1-f/2)^N+N*(f/2)*(1-f/2)^(N-1)*(1-f), #Probability sporadic/simplex
    Pfam   = 1-((1-f/2)^N+N*(f/2)*((1-f/2)^(N-1))*(1-f)) #Probability familial (= 1 - Punasc - Pspor)
  ),nrow=length(f))
  colnames(dis_states) <- c("f","Punasc","Pspor","Pfam")
  
  #Build lookup table recording EXPECTED value of ObsProbX at each f value.
  #A series of LookupX values, each respective to a value of f, are calculated as a proportion of those states defined in Operation 1.
  #Calculation: Probability of state X / sum(probabilities of states defined)
  if(states=="fsu") {
    LookupX = dis_states[,"Pfam"] #F/(F+S+U) = F/1
    LabelX = "familial"
  } else if(states=="fs") { 
    LookupX = dis_states[,"Pfam"]/(dis_states[,"Pfam"]+dis_states[,"Pspor"]) #F/(F+S)
    LookupX[which(is.na(LookupX))] <- 0 #Replace any NaNs with 0 - NaN is expected for first value if penetrance is 0
    LabelX = "familial"
  } else if(states=="fu") {
    LookupX = dis_states[,"Pfam"]/(dis_states[,"Pfam"]+dis_states[,"Punasc"]) #F/(F+U)
    LabelX = "familial"
  } else if(states=="su") {
    LookupX = dis_states[,"Pspor"]/(dis_states[,"Pspor"]+dis_states[,"Punasc"]) #S/(S+U)
    LabelX = "sporadic"
  } else if(states=="au") {
    LookupX = (dis_states[,"Pfam"]+dis_states[,"Pspor"])/(dis_states[,"Pfam"]+dis_states[,"Pspor"]+dis_states[,"Punasc"]) #(F+S)/(F+S+U) = (F+S)/1 = A/1
    LabelX = "affected"
  } else {
    LookupX = as.numeric(NA)
    LabelX <- NA
  }
  
  #Build a two column Lookup table of LookupX values with respective f values
  LookupTable <- matrix(c(LookupX, f),nrow=length(f))
  
  
  #Steps 3 and 4:
  #Prepare for Step 4 estimate adjustment:
  #Generate error curve for sample, according to the specified N, and whether upload_sibships or define_sibships have been defined
  errorfit <- adpenetrance.errorfit(states=states, setmean=N,samp_size=90000,seed=24,define_sibstructure=define_sibstructure)
  
  #Compare ObsProbX to LookupTable - retrieve loci of LookupX with nearest value to ObsProbX (and intervals) 
  #The values should be comparable - large disparity suggests ObsProbX exceeds or is lesser than rate expected at penetrance= 1 or 0
  
  #Prepare matrices to store the locus of the nearest LookupX values and a matrix to store adjusted f values
  Corr_f <- loci <- matrix(NA,ncol=3,nrow=1)
  colnames(Corr_f) <- colnames(loci) <- c("Lower CI", "Estimate", "Upper CI")
  
  #Matching ObsProbX to closest LookupX - store the locus retrieved in 'loci' object defined above
  if(!is.na(ObsProbX)) {
    loci[,"Estimate"] <- which(abs(LookupTable[,1]-as.vector(ObsProbX))==min(abs(LookupTable[,1]-as.vector(ObsProbX)),na.rm=T)) #Locus for the estimate
    
    #Predict the error in f for the estimate, based on the value of f obtained, using the correction curve
    correctionVal<- unname(predict(errorfit$fitbest,data.frame(x=LookupTable[,2][loci[,"Estimate"]])))
    #Adjust penetrance estimate obtained by the value obtained. Store in Corr_f matrix
    Corr_f[,"Estimate"] <- LookupTable[,2][loci[,"Estimate"]]+correctionVal
  }
  
  #If error propagation has been included:
  #Repeat process at disease state rates observed at confidence interval bounds
  if(!is.na(ObsProbXSE)) {
    loci[,"Lower CI"] <- which(abs(LookupTable[,1]-as.vector(lowerCI))==min(abs(LookupTable[,1]-as.vector(lowerCI)),na.rm=T)) #Locus for the lower bound
    correctionVal<- unname(predict(errorfit$fitbest,data.frame(x=LookupTable[,2][loci[,"Lower CI"]])))
    Corr_f[,"Lower CI"] <- LookupTable[,2][loci[,"Lower CI"]]+correctionVal
    
    loci[,"Upper CI"] <- which(abs(LookupTable[,1]-as.vector(upperCI))==min(abs(LookupTable[,1]-as.vector(upperCI)),na.rm=T)) #Locus for the upper bound
    correctionVal<- unname(predict(errorfit$fitbest,data.frame(x=LookupTable[,2][loci[,"Upper CI"]])))
    Corr_f[,"Upper CI"] <- LookupTable[,2][loci[,"Upper CI"]]+correctionVal
    
  }
  
  
  #Build output matrix:
  #Row 1: ObsProbX (with or without CI and standard error)
  #Row 2: Value of LookupX closest to ObsProbX (with or without CI) - identified at row stored in 'loci' object
  #Row 3: Penetrance estimate associated with LookupX (with or without CI)
  #Row 4: Adjusted penetrance estimate
  if(!is.na(ObsProbXSE)) { #Include confidence intervals if propagation performed
    #Truncate adjusted estimates to fit within [0,1]
    Corr_f[which(Corr_f<0)] <- 0
    Corr_f[which(Corr_f>1)] <- 1
    
    Output = matrix(c(lowerCI, ObsProbX, upperCI, ObsProbXSE, #Estimate for X rate, with CI and SE given
                      LookupTable[,1][loci[,"Lower CI"]], LookupTable[,1][loci[,"Estimate"]], LookupTable[,1][loci[,"Upper CI"]], NA, #LookupX nearest estimate and CI
                      LookupTable[,2][loci[,"Lower CI"]], LookupTable[,2][loci[,"Estimate"]], LookupTable[,2][loci[,"Upper CI"]], NA, #Corresponding penetrance estimate
                      Corr_f,NA), #Adjusted penetrance estimate
                    byrow =T, ncol=4, nrow=4)
    
    colnames(Output) <- c("Lower CI","Estimate", "Upper CI", "Standard error")
    
  } else {   #Create without confidence intervals
    #Truncate corrected estimates to fit within [0,1]
    Corr_f[which(Corr_f<0)] <- 0
    Corr_f[which(Corr_f>1)] <- 1
    
    Output = matrix(c(ObsProbX, #Estimate for X rate
                      LookupTable[,1][loci[,"Estimate"]], #LookupX nearest estimate
                      LookupTable[,2][loci[,"Estimate"]], #Corresponding penetrance estimate
                      Corr_f[,"Estimate"]), #Adjusted penetrance estimate
                    ncol=1, nrow=4)
    
    
    colnames(Output) <- c("Estimate")
  }
  
  #Construct Row names for the first two rows
  #LabelX is the name of the state modelled as state X
  R1 <- paste("Observed", LabelX, "rate")
  R2 <- paste("Expected", LabelX, "rate")
  
  #Assign row names
  rownames(Output) <- c(R1, R2, "Unadjusted Penetrance","Adjusted Penetrance")
  
  if(isTRUE(include_MLE)){
    #ADD MLE STEP
    #Define likelihood function, f (penetrance) to be maximised
    #MLE to be performed using nlm, so negML function returns dbinom()*-1 
    negML <- function(f){
      #Al-Chalabi & Lewis (2011) model equations
      Punasc = (1-f)*(1-f/2)^N                                #Probability unaffected
      Pspor  = f*(1-f/2)^N+N*(f/2)*(1-f/2)^(N-1)*(1-f)        #Probability sporadic/simplex
      Pfam   = 1-((1-f/2)^N+N*(f/2)*((1-f/2)^(N-1))*(1-f))    #Probability familial (= 1 - Punasc - Pspor)
      
      #Define r for comparison to ObsProbX - derive r according to represented states 
      if(states=="fsu") {
        r = Pfam #F/(F+S+U) = F/1
      } else if(states=="fs") { 
        r=Pfam/(Pfam+Pspor)     #F/(F+S)
      } else if(states=="fu") {
        r=Pfam/(Pfam+Punasc)     #F/(F+U)
      } else if(states=="su") {
        r=Pfam/(Pspor+Punasc)     #S/(S+U)
      } else if(states=="au") {
        r=(Pfam+Pspor)/(Pfam+Pspor+Punasc)     #(F+S)/(F+S+U) = (F+S)/1 = A/1
      }
      
      trials=10000000 #Define number of trials - used convert RX into integer necessary for parsing by dbinom 
      RXml = round(ObsProbX,7)*trials #Adjust the ObsProbX proportion by the number of trials - rounding to guarantee integer
      
      seed = 24
      set.seed(seed)
      dbinom(RXml,trials,r)*-1       #Negative likelihood estimate: dbinom(number of events,number of trials, probabilities)
    }
    
    #Prepare results table for MLE estimate
    MLresults <-matrix(NA,ncol=3,nrow=3)
    colnames(MLresults) <- c("Lower CI","Estimate", "Upper CI") #Assign column names
    rownames(MLresults) <- c("Unadjusted Penetrance", "Adjusted Penetrance","NLM convergence") #Assign row names
    
    #Use nlm function to minimise negML, finding ML penetrance estimate
    #Starting value is set by the lookup table unadjusted estimate, close to true penetrance value...
    MLestimate <- nlm(negML,Output["Unadjusted Penetrance","Estimate"],stepmax=0.5)
    
    MLresults[1,2] <- MLestimate$estimate #Uncorrected ML estimate
    correctionVal<- unname(predict(errorfit$fitbest,data.frame(x=MLresults[1,2]))) #Predict the error in f for the estimate, based on the value of f obtained, using the correction curve
    MLresults[2,2] <- MLresults[1,2]+correctionVal #Corrected ML estimate
    MLresults[3,2] <- MLestimate$code #Check for convergence of MLE model
    
    if(!is.na(ObsProbXSE)) { #If error terms given, at CI bounds of penetrance estimate
      MLlower <- nlm(negML,Output["Unadjusted Penetrance","Lower CI"],stepmax=0.5)
      MLupper <- nlm(negML,Output["Unadjusted Penetrance","Upper CI"],stepmax=0.5)
      
      MLresults[1,1] <- MLlower$estimate #Uncorrected ML estimate
      correctionVal<- unname(predict(errorfit$fitbest,data.frame(x=MLresults[1,1]))) #Predict the error in f for the estimate, based on the value of f obtained, using the correction curve
      MLresults[2,1] <- MLresults[1,1]+correctionVal #Corrected ML estimate
      MLresults[3,1] <- MLestimate$code #Check for convergence of MLE model
      
      MLresults[1,3] <- MLupper$estimate #Uncorrected ML estimate
      correctionVal<- unname(predict(errorfit$fitbest,data.frame(x=MLresults[1,3]))) #Predict the error in f for the estimate, based on the value of f obtained, using the correction curve
      MLresults[2,3] <- MLresults[1,3]+correctionVal #Corrected ML estimate
      MLresults[3,3] <- MLestimate$code #Check for convergence of MLE model
    }
    
    return(list(output=Output, MLE=MLresults,errorfit=summary(errorfit$fitbest)))
  } else {
    
    return(list(output=Output,errorfit=summary(errorfit$fitbest)))
  }
  
}