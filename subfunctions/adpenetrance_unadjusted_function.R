#This script is a reduced version of the adpenetrance function which follows only Steps 1-3 of the approach
  #and therefore penetrance estimates are not adjusted as in Step 4
  #The function is to be called by adpenetrance.errorfit function for determining Step 3 penetrance estimates
    #expected across known ground truth penetrance values estimated for a population of mean sib-size N.
    #It should be used only within the main adpenetrance function, which itself calls adpenetrance.errorfit
    #The main adpenetrance function can be loaded from the "adpenetrance_function.R" script available within our github repository:
    #https://github.com/thomaspspargo/adpenetrance

adpenetrance.unadjusted <- function(N, MF=NA, MS=NA, MA=NA, MU=NA, PA=NA, PF=NA, MF_SE=NA, MS_SE=NA, MU_SE=NA, MA_SE=NA, Zout=1.96,
                                    useG=0,
                                    RX=NA, RX_SE=NA, states="none"){
  
  #--------------------------------------------#
  #   IF ELSE statements for input variables   #
  #--------------------------------------------#
  #First, define the states modelled
  
  #Perform these calculations only if RX is not given directly:
  
  #Generate value for states term within the function - necessary to pass through states to adpenetrance.errorfit
  if(states=="none"){ 
    if(!is.na(MF) && !is.na(MS) && !is.na(MU) && is.na(MA)) {
      states="fsu"
    } else if(!is.na(MF) && !is.na(MS) && is.na(MU) && is.na(MA)) {
      states="fs"
    } else if(!is.na(MF) && is.na(MS) && !is.na(MU) && is.na(MA)) {
      states="fu"
    } else if(is.na(MF) && !is.na(MS) && !is.na(MU) && is.na(MA)) {
      states="su"
    } else if(is.na(MF) && is.na(MS) && !is.na(MU) && !is.na(MA)) {
      states="au"
    } else {
      stop("No valid disease state combination has been defined. Variant frequency or RX estimates should be defined for any two or three of the familial, sporadic, and unaffected states or the affected and unaffected states. Please check that variant frequency estimates have been defined for a valid combination of states or that the 'RX' and 'states' arguments are properly defined.")
    }
    
    #Logic check weighting factors
    if(grepl("f|s",states)){
      if(is.na(PF)) {
        stop("PF has not been defined; This is required as a weighting factor in calculations involving variant frequencies MF or MS")
      }
    }
    if(grepl("a|u",states)){
      if(is.na(PA)) {
        stop("PA has not been defined. This is required as a weighting factor in calculations involving variant frequencies MA or MU")
      }
    }
  } else if(!any(c('fsu','fu','fs','au','su') %in% states)){ #Change states value to "incorrect" if nothing has been defined
      stop("No valid disease state combination has been defined. Variant frequency or RX estimates should be defined for any two or three of the familial, sporadic, and unaffected states or the affected and unaffected states. Please check which states have been indicated in the 'states' object")
  }
  
  
  #The following statements evaluate how data are specified in the function to determine how to calculate the observed rate of state X.
  #All valid combinations are defined and unacceptable selections will be stored as NA and the output will not be valid
  if(is.na(RX)){
    
    if(states=="fsu") { #If familial, sporadic, unaffected specified, and affected not specified
      X <- MF           #Familial state is state X
      WeiX <- PF*PA     #Weighting variable
      
      Y <- MS           #Sporadic state is state Y
      WeiY = (1-PF)*PA  #Weighting variable
      
      Z <- MU           #Unaffected is state Z (only used in this formulation - after this only X and Y need specification)
      WeiZ <- 1-PA      #Weighting variable
      
      if(!is.na(MF_SE) && !is.na(MS_SE) && !is.na(MU_SE)) { #If error terms are given for all appropriate variables, define error terms
        XSE <- MF_SE
        YSE <- MS_SE
        ZSE <- MU_SE
      } else if (!is.na(MF_SE) || !is.na(MS_SE) || !is.na(MU_SE)) { #If error terms are given for some but not all variables, do not define error terms and return warning message
        XSE <- NA
        YSE <- NA
        ZSE <- NA
        
        warning("Error terms have only been given for a subset of the variant frequencies defined. Please ensure that these are given for all included states")
        } else { 
        XSE <- NA
        YSE <- NA
        ZSE <- NA
      }
      
    } else if(states=="fs") { #If familial and sporadic specified, and unaffected and affected not specified
      X <- MF           #Familial state is state X
      WeiX <- PF        
      
      Y <- MS           #Sporadic state is state Y
      WeiY <- 1-PF      
      
      if(!is.na(MF_SE) && !is.na(MS_SE)) { #If error terms are given for all appropriate variables, define error terms
        XSE <- MF_SE
        YSE <- MS_SE
      } else if (!is.na(MF_SE) || !is.na(MS_SE)) { #If error terms are given for some but not all variables, do not define error terms and return warning message
        XSE <- NA
        YSE <- NA
        
        warning("Error terms have only been given for a subset of the variant frequencies defined. Please ensure that these are given for all included states")
      } else {
        XSE <- NA
        YSE <- NA
      }
      
      
    } else if(states=="fu") { #If familial and unaffected specified, and sporadic and affected not specified
      X <- MF         #Familial state is X
      WeiX <- PF*PA 
      
      Y <- MU         #Unaffected state is Y
      WeiY <- 1-PA
      
      if(!is.na(MF_SE) && !is.na(MU_SE)) { #If error terms are given for all appropriate variables, define error terms
        XSE <- MF_SE
        YSE <- MU_SE
      } else if (!is.na(MF_SE) || !is.na(MU_SE)) { #If error terms are given for some but not all variables, do not define error terms and return warning message
        XSE <- NA
        YSE <- NA
        
        warning("Error terms have only been given for a subset of the variant frequencies defined. Please ensure that these are given for all included states")
      } else {
        XSE <- NA
        YSE <- NA
      }
      
    } else if(states=="su") { #If sporadic and unaffected specified, and familial and affected not specified
      X <- MS         #Sporadic state is X
      WeiX <- (1-PF)*PA
      
      Y <- MU         #Unaffected state is Y
      WeiY <- 1-PA
      
      if(!is.na(MS_SE) && !is.na(MU_SE)) { #If error terms are given for all appropriate variables, define error terms
        XSE <- MS_SE
        YSE <- MU_SE
      } else if (!is.na(MS_SE) || !is.na(MU_SE)) { #If error terms are given for some but not all variables, do not define error terms and return warning message
        XSE <- NA
        YSE <- NA
        
        warning("Error terms have only been given for a subset of the variant frequencies defined. Please ensure that these are given for all included states")
      } else {
        XSE <- NA
        YSE <- NA
      }
      
    } else if(states=="au") { #If unaffected and affected specified, familial and sporadic not specified
      X <- MA         #Affected state is X
      WeiX <- PA
      
      Y <- MU         #Unaffected state is Y
      WeiY <- 1-PA
      
      if(!is.na(MA_SE) && !is.na(MU_SE)) { #If error terms are given for all appropriate variables, define error terms
        XSE <- MA_SE
        YSE <- MU_SE
      } else if (!is.na(MA_SE) || !is.na(MU_SE)) { #If error terms are given for some but not all variables, do not define error terms and return warning message
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
    } else {ObsProbXSE <- NA_real_} #If error terms not given, skip this operation
    
    
    
    
  } else if(!is.na(RX)){ #IF RX is specified directly, the above is skipped and ObsProbX is defined directly by the user
    
    #Pass error if any variant frequencies have also been defined
    if(any(!is.na(MF),!is.na(MS),!is.na(MU),!is.na(MA))){
      stop("Input data have been provided for at least 1 variant frequency as well as RX. Please specify either the 'RX' and the 'states' terms, or specify variant frequency estimates and required weighting factors for each represented disease state")
    }
    
    ObsProbX = RX #Define ObsProbX directly from RX
    
    if(!is.na(RX_SE)){
      ObsProbXSE = RX_SE #Define SE in ObsProbX from the error specified
    } else {
      ObsProbXSE <- NA_real_ #Define as NA if error is unspecified
    }
    
  }
  

  #Calculate confidence intervals for ObsProbX using z-score conversion, if error is known.
  #No if statement according to whether RX is given or variant frequencies
  if (!is.na(ObsProbXSE)) {
    lowerCI = ObsProbX-(Zout*ObsProbXSE)   #Lower interval
    upperCI = ObsProbX+(Zout*ObsProbXSE)   #Upper interval
    
  } else {
    lowerCI=NA_real_
    upperCI=NA_real_
  }
  
  #Step 2:
  #Construct penetrance value sequence for lookup table. Values generated are between 0 and 1 at increments of 0.0001
  f = seq(0,1, by=0.0001)
  
  #Apply extended Al-Chalabi & Lewis (2011) equations to calculate P(unaffected), P(sporadic), P(familial) 
  #Do this at the defined specified sibship size for all values of f- store in matrix of length(f)
  #when g=0, the default, these eqns are synonymous with the original Al-Chalabi & Lewis formulae
  g=useG
  
  #Calculate proportions under the extended model including residual disease risk g. If g=0, synonymous with original equations
  Punasc = (1-f)*(1-f/2-g/2)^N*(1-g) #Probability unaffected
  Pspor  = f*(1-f/2-g/2)^N*(1-g)+
    (1-f)*N*(f/2+g/2)*(1-f/2-g/2)^(N-1)*(1-g)+
    (1-f)*(1-f/2-g/2)^N*g #Probability sporadic/simplex
  Pfam   = 1-(Punasc+Pspor) #Probability familial (= 1 - Punasc - Pspor)
  message("Disease state proportions are calculated expecting residual disease risk (g) of ",signif(g,3)," among family members not harbouring the variant")
  
  #Store state probabilities in a matrix
  dis_states <- matrix(c(
    f,
    Punasc = Punasc,
    Pspor  = Pspor,
    Pfam   = Pfam
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
    LookupX = NA_real_
    LabelX <- NA
  }
  
  #Build a two column Lookup table of LookupX values with respective f values
  LookupTable <- matrix(c(LookupX, f),nrow=length(f))

  
  #Step 3:
  #Compare ObsProbX to LookupTable - retrieve loci of LookupX with nearest value to ObsProbX (and intervals) 
    #The values should be comparable - large disparity suggests ObsProbX exceeds or is lesser than rate expected at penetrance= 1 or 0
  
  #Prepare matrices to store the locus of the nearest LookupX values and a matrix to store adjusted f values
  loci <- matrix(NA,ncol=3,nrow=1)
  colnames(loci) <- c("Lower CI", "Estimate", "Upper CI")
  
  #Matching ObsProbX to closest LookupX - store the locus retrieved in 'loci' object defined above
  if(!is.na(ObsProbX)) {
    tryCatch(
      loci[,"Estimate"] <- which(abs(LookupTable[,1]-as.vector(ObsProbX))==min(abs(LookupTable[,1]-as.vector(ObsProbX)),na.rm=T)) #Locus for the estimate
      ,error = function(x){warning(x)
        
        #Extract matched indices and sanitise for error reporting
        match<- which(abs(LookupTable[,1]-as.vector(ObsProbX))==min(abs(LookupTable[,1]-as.vector(ObsProbX)),na.rm=T))
        matchtab <- LookupTable[match,]
        colnames(matchtab) <- c("Expected RX_i","corresponding penetrance")
        
        cat("------------------------------\n")
        cat("ADPenetrance error catch when indexing lookup table to identify expected disease state rate corresponding to RX.\n")
        cat("Multiple rows in the lookup table may have been matched.\nThis could reflect equal distance between the Observed disease state rate and two positions in the lookup table\n current parameter details are as follows:\n")
        
        cat("Observed disease state rate (RX):\n")
        print(ObsProbX)
        cat("Estimated average sibship size (N):\n")
        print(N)
        cat("Estimated residual disease risk (g):\n")
        print(useG)
        cat("States modelled:\n")
        print(states)
        
        cat("Matched row indices in lookup table (only 1 match is expected):\n")
        print(match)
        cat("Corresponding rows in lookup table:\n")
        print(matchtab)
        
        cat("Penetrance estimation will proceed using the first index match\n")
        loci[,"Estimate"] <<- match[1] #Locus for the estimate
        cat("-----------------------------\n")
      }) #End trycatch
  }
  
  #If error propagation has been included:
    #Repeat process at disease state rates observed at confidence interval bounds
  if(!is.na(ObsProbXSE)) {
    
    #Perform indexing; use trycatch to give more informative error handling
    tryCatch(
      loci[,"Lower CI"] <- which(abs(LookupTable[,1]-as.vector(lowerCI))==min(abs(LookupTable[,1]-as.vector(lowerCI)),na.rm=T)) #Locus for the lower bound
      ,error = function(x){warning(x)
        
        #Extract matched indices and sanitise for error reporting
        match<- which(abs(LookupTable[,1]-as.vector(lowerCI))==min(abs(LookupTable[,1]-as.vector(lowerCI)),na.rm=T))
        matchtab <- LookupTable[match,]
        colnames(matchtab) <- c("Expected RX_i","corresponding penetrance")
        
        cat("------------------------------\n")
        cat("ADPenetrance error catch when indexing lookup table to identify expected disease state rate corresponding to lower confidence interval of RX.\n")
        cat("Multiple rows in the lookup table may have been matched.\nThis could reflect equal distance between the observed disease state rate and two positions in the lookup table\n current parameter details are as follows:\n")
        
        cat("Observed disease state rate (RX) lower CI value:\n")
        print(lowerCI)
        cat("Estimated average sibship size (N):\n")
        print(N)
        cat("Estimated residual disease risk (g):\n")
        print(useG)
        cat("States modelled:\n")
        print(states)
        
        cat("Matched row indices in lookup table (only 1 match is expected):\n")
        print(match)
        cat("Corresponding rows in lookup table:\n")
        print(matchtab)
        
        cat("Penetrance estimation will proceed using the first index match\n")
        loci[,"Lower CI"] <<- match[1] #Locus for the estimate
        cat("-----------------------------\n")
      }) #End trycatch
    
    #Perform indexing; use trycatch to give more informative error handling
    tryCatch(
      loci[,"Upper CI"] <- which(abs(LookupTable[,1]-as.vector(upperCI))==min(abs(LookupTable[,1]-as.vector(upperCI)),na.rm=T)) #Locus for the upper bound
      ,error = function(x){warning(x)
        
        #Extract matched indices and sanitise for error reporting
        match<- which(abs(LookupTable[,1]-as.vector(upperCI))==min(abs(LookupTable[,1]-as.vector(upperCI)),na.rm=T))
        matchtab <- LookupTable[match,]
        colnames(matchtab) <- c("Expected RX_i","corresponding penetrance")
        
        cat("------------------------------\n")
        cat("ADPenetrance error catch when indexing lookup table to identify expected disease state rate corresponding to upper confidence interval of RX.\n")
        cat("Multiple rows in the lookup table may have been matched.\nThis could reflect equal distance between the observed disease state rate and two positions in the lookup table.\nCurrent parameter details are as follows:\n")
        
        cat("Observed disease state rate (RX) upper CI value:\n")
        print(upperCI)
        cat("Estimated average sibship size (N):\n")
        print(N)
        cat("Estimated residual disease risk (g):\n")
        print(useG)
        cat("States modelled:\n")
        print(states)
        
        cat("Matched row indices in lookup table (only 1 match is expected):\n")
        print(match)
        cat("Corresponding rows in lookup table:\n")
        print(matchtab)
        
        cat("Penetrance estimation will proceed using the first index match\n")
        loci[,"Upper CI"] <<- match[1]
        cat("-----------------------------\n")
      }) #End trycatch 
      
  }
  
  
  #Build output matrix:
  #Row 1: ObsProbX (with or without CI and standard error)
  #Row 2: Value of LookupX closest to ObsProbX (with or without CI) - identified at row stored in 'loci' object
  #Row 3: Penetrance estimate associated with LookupX (with or without CI)
  if(!is.na(ObsProbXSE)) { #Include confidence intervals if propagation performed
   
    
    Output = matrix(c(lowerCI, ObsProbX, upperCI, ObsProbXSE, #Estimate for X rate, with CI and SE given
                      LookupTable[,1][loci[,"Lower CI"]], LookupTable[,1][loci[,"Estimate"]], LookupTable[,1][loci[,"Upper CI"]], NA, #LookupX nearest estimate and CI
                      LookupTable[,2][loci[,"Lower CI"]], LookupTable[,2][loci[,"Estimate"]], LookupTable[,2][loci[,"Upper CI"]], NA #Corresponding penetrance estimate
                      ),
                    byrow =T, ncol=4, nrow=3)
    
    colnames(Output) <- c("Lower CI","Estimate", "Upper CI", "Standard error")
    
  } else {   #Create without confidence intervals
    
    Output = matrix(c(ObsProbX, #Estimate for X rate
                      LookupTable[,1][loci[,"Estimate"]], #LookupX nearest estimate
                      LookupTable[,2][loci[,"Estimate"]] #Corresponding penetrance estimate
                      ), #Adjusted penetrance estimate
                    ncol=1, nrow=3)

    
    colnames(Output) <- c("Estimate")
  }
  
  #Construct Row names for the first two rows
  #LabelX is the name of the state modelled as state X
  R1 <- paste("Observed", LabelX, "rate")
  R2 <- paste("Expected", LabelX, "rate")
  
  #Assign row names
  rownames(Output) <- c(R1, R2, "Unadjusted Penetrance")
  
  return(list(output=Output,states=states,ResidualRiskG=useG))
  
}