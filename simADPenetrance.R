#####
# Specifications for a function to simulate error in lifetime Penetrance estimation according to age of sampling under various scenarios: simADPenetrance
# Examples of use are provided in the approach_validation portion of the GitHub repository where this is function is hosted (https://github.com/ThomasPSpargo/adpenetrance/)
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk) 
# Please get in touch if you have any issues.
####

#Load adpenetrance function and subfunctions
source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/adpenetrance_function.R")

#Load subfunctions functions for simulating families, via lapply
simFunctions_root<- ("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/subfunctions/")
simFunctions <- c("affAtAge.R","bindRows.R","genFamily.R","varChars.R")
invisible(lapply(paste0(simFunctions_root,simFunctions),source))
rm(simFunctions_root,simFunctions) #Clean temp objects from environment


######
### Set simADPenetrance function
######

# ## Defaults:
# simADPenetrance(onsetRateDiff = 1, f = c(0.25,0.50,0.75,1), f_compare = c(0.2,0.4,0.6,0.8,1), amplify_f = 3,g=0,states=c("fsu","fs","fu","su","au"),
#                 numfamilies_var = 50000, scale_novar = 1, sibstructureCustom = NULL, which_sibstructure = c(1,2), nameSibstructureCustom = "Custom sibstructure",numsteps = 10,                      
#                 f_cohort_only=FALSE, eldestAt0=FALSE, stepHazard=FALSE, benchmark = FALSE)

# ## INPUTS:
# f                       = Numeric vector (all values must be between 0 and 1) defining the lifetime penetrance values to test (the f vector can have only one element).
# f_compare               = Numeric vector (all values must be between 0 and 1) which can have only one element, specifying alternative variants occurring in families without the each tested f variant;
#                           each element of f_compare is assigned to subset of the simulated families not harbouring f (see: scale_novar argument).
# onsetRateDiff           = Define whether differences exist in rate of onset between variant and non-variant group.
# amplify_f               = Set number of loops for each defined f. Simulation results will be averaged across these to reduce any effects of pseudo-randomisation.
# g                       = Residual disease risk in people without either f or f_compare variant, defaults to 0.
# states                  = Define the disease state combinations to model Defaults to: c("fsu","fs","fu","su","au"), but any subset of these 5 options can be modelled.
# numfamilies_var         = Number of families to simulate which harbour each variant f.
# scale_novar             = Numeric to indicate scaling of families relative to the value of numfamilies_var. Generates numfamilies_var*scale_novar families harbouring one of the competing variants.
#                           defined in f_compare (e.g. if numfamilies_var*scale_novar=50,000, and f_compare is a vector defining 5 penetrance values for competing variants,
#                           then numfamilies_var/length(f_compare)=50,000/5=10,000 families with each f_compare variant).
# sibstructureCustom      = Specify a custom sibship population.
# which_sibstructure      = Specify which sibship distributions (sibstructure) to test; 1=UK,2=NS,3=Custom. Defaults to c(1,2).
#                           Further details about the preset sibship structures are provided in Figure S3 of the manuscript associate with this repository.
# nameSibstructureCustom  = If providing a custom sibstructure (see sibstructureCustom and which_sibstructure), specify a string to name the custom structue in the output plot.
# numsteps                = Number of time points (ages) across which a person may become affected after the 0 time.
# f_cohort_only           = Logical, defaults to FALSE. Run simulation with f cohort only, generating RX based on proportion of states within families
#                           Certain arguments are ignored if TRUE. These are: f_compare, scale_novar.
# eldestAt0               = Logical, defaults to FALSE and is passed to the genFamily.
#                           If TRUE, adjusts the family ages such that the eldest parent is assigned age 0 at the first time of sampling,
#                           where 0 is the last time point where no family members could be affected by disease as indicated by risks stored within the var.Char matrix.
#                           If FALSE, the youngest family member is at age 0, and thus at the first time of sampling all other family members have some probability of being affected by disease, according to their variant status and risk indicated in var.Char.
# stepHazard              = Logical, defaults to FALSE and is passed to the genFamily subfunction.
#                           If FALSE, disease risks at each time point are determined using the affAtAge function according to the person's age and variant status (i.e. "var", "novar", or "possvar", where "possvar" are sibs who have a 50% chance of inheriting. a variant from a "var" parent)
#                           If TRUE, disease risk at the first time of sampling are determined using the affAtAge function according to person's age, variant status ("var", "novar", or "possvar").
#                           Thereafter, cumulative risk across each age step is determined based on additional risk between each step of var.Chars across the relevant possible ages of onset and variant status ("var", "novar", or "possvar").
#                           Setting TRUE is not recommended because higher (familial/sporadic) disease state proportions are systematically underrepresented under the 'stepwise' approach.
#benchmark = FALSE        = Provide benchmark to indicate how long the analysis may take given current settings. Benchmarks are currently very approximate



## OUTPUT:
#1. A ggplot indicating divergence sampling time on the x-axis and divergence from true penetrance estimates on the y axis.
#2. A dataframe storing data visualised in the plot
#Column headers for the data frame file are:
 #Population              = The modelled population as defined in which_sibstructure
 #Tailoring               = Whether adjusted penetrance estimates were adjusted using:
 #                          Poisson: the default sibhship distribution generated by default within the adpenetrance function.
 #                          Tailored: The an approximation of the sibship distribution for the population modelled
 #                          No correction: 'Step-3' penetrance estimates before applying any correction
 #TruePenetrance          = The true penetrance, each row will have one of the values specified in the argument 'f'
 #States.Modelled         = One of the 5 valid disease state combinations, each row will be one of the options set in the states argument
 #modif                   = Time from first sampling, starting from 0 until the youngest member of each simulated family has passed the timepoint indicated in the 'numsteps' argument
 #Estimate.difference     = Difference between true and estimated penetrance, averaged across number of repetitions for each value of 'f' (as specified in 'amplify_f' argument).


# Examples of use are provided in the GitHub repository hosting this function (https://github.com/ThomasPSpargo/adpenetrance/).
 
simADPenetrance <- function(
  onsetRateDiff = 1,                   
  f = c(0.25,0.50,0.75,1),             
  f_compare = c(0.2,0.4,0.6,0.8,1),   
  amplify_f = 3,                      
  g=0,                                
  states=c("fsu","fs","fu","su","au"), 
  numfamilies_var = 50000,            
  scale_novar = 1,
  sibstructureCustom = NULL,          
  which_sibstructure = c(1,2),        
  nameSibstructureCustom = "Custom sibstructure",
  numsteps = 10,                      
  f_cohort_only=FALSE,                
  eldestAt0=FALSE,
  stepHazard=FALSE,
  benchmark = FALSE
){
  
  ######
  ### Return Pre-analysis warnings/errors, as required
  ######
  if(!is.null(sibstructureCustom) && !(3 %in% which_sibstructure)){
    warning("A custom sibstructure is defined but value '3' is not included in the which_sibstructure argument and so this has been ignored.")
  }
  
  if(any(f>1) || any(f<0)){
    stop("Impossible values have been passed to the 'f' argument. Please make sure all values are between 0 and 1")
  }
  
  if(any(f_compare>1) || any(f_compare<0)){
    stop("Impossible values have been passed to the 'f_compare' argument. Please make sure all values are between 0 and 1")
  }
  
  ######
  ### Set cohort sizes and generate benchmark estimate of time to run if benchmark=TRUE
  ######
  if(f_cohort_only==FALSE){
    #Indicate numer of families with comparator variants
    numfamilies_novar = numfamilies_var*scale_novar  #people with different risk variants
    
    #Identify the total number of families to generate
    numfamilies <- numfamilies_var+numfamilies_novar
  } else if (f_cohort_only){
    numfamilies <- numfamilies_var
  }
  if(benchmark){ #If benchmarking, run only this section
    
    #Approx time in seconds for major steps
    gen100000family11times <- 60                            #to generate 100,000 families across 11 time points (in seconds)
    gen1family1time <- (gen100000family11times/100000)/11   #to generate 1 family across 1 time points (in seconds)
    onerun <- 1.5                                           #to run ADPenetrance once
    
    #Estimate number of families and time needed
    numToGen <- numfamilies*
      length(f)*amplify_f*       #Number of l loops
      length(which_sibstructure) #number of p loops
    timeToGen <- gen1family1time*numToGen*numsteps+1
    
    #Estimate number of adpenetrance calls and time needed
    timesCalled <- length(f)*amplify_f* #Number of l loops
      length(which_sibstructure)* #number of p loops
      numsteps+1* #number of time points to estimate
      length(states)*2 #Number of 's' and 't' loops, one with poisson and one with a tailored Sibstructure
    estRuntime <- onerun*timesCalled
    
    message("Approximate benchmarking for the most time-consuming steps:\n
    A total of ", signif(numToGen,2), " families will be simulated to obtain disease state assignments at ",numsteps+1," time points across all loops.
    This requires approximately ", signif(timeToGen,2), " seconds.\n
    Adpenetrance will be called ", signif(timesCalled,2)," times during this simulation. This alone requires approximately ",signif(estRuntime,2)," seconds.\n
    We recommend allowing over ", signif((estRuntime+timeToGen)/60,2)," minutes to run this simulation.
    Note that these benchmarks are a very approximate lower bound for time required.\n
    Call function again with benchmarking=FALSE to perform simulation.")
    
    option <- options(show.error.messages=FALSE)
    on.exit(options(option))
    stop()
  }
  
  ######
  ### set true f (penetrance) values to test prepare output objects
  ######
  f <- sort(rep(f,amplify_f)) #Amplify to test each value amplify_f times, sort into ascending order
  sims <- length(f)   #identify number of top-level loops to run
  
  #Prepare empty list/vectors to store output summary at each states modelled
  output_colnames<- c("Population","Tailoring","Seed","mean_N","MF","MS","MU","PF","PA","RX","modif","States.Modelled","TruePenetrance","UncorrectedPen", "CorrectedPen")
  statetemp <- vector(mode="list",length=2) 
  result <- NULL 
  
  sim_seed <- sample(1:30000,sims) #Generate random seeds for use in each loop
  
  ######
  ### Set population structure(s)
  ######
  
  if(1 %in% which_sibstructure){
    #Define UK population weights:
    #This structure is to be supplied to define_sibstructure in some simulations, when Step-4 error correction is tailored to a known population distribution
    #Column 1 stores sibship levels, Column 2 stores probabilities of each level
    sibstructureUK<- cbind(0:4,c(0.18,0.18,0.37,0.16,0.11))
    
  } else {
    sibstructureUK <- NULL
  }
  
  if(2 %in% which_sibstructure){
    #Determine Next Steps population weights (see final two rows of Table 2 in publication):
    NextStepN<- 3647+10588                                                    #Total sample size
    Sibs <- c(393+662, 1231+3926,952+3059,556+1591,248+648,119+344,148+358)   #Sib levels 1:7, vector contains number in each sib level
    #Note that 'only child' represents sibship size 1
    Psibs <- Sibs/NextStepN                                                   #Psibs stores the probabilities of each sibship size
    
    #Build matrix
    #This structure is to be supplied to define_sibstructure in some simulations, when Step-4 error correction is tailored to a known population distribution
    #Column 1 stores sibship levels, Column 2 stores probabilities of each level
    sibstructureNS<- cbind(1:7, Psibs)
    
  } else{
    sibstructureNS <- NULL
  }
  
  if(3 %in% which_sibstructure){
    if (is.vector(sibstructureCustom,mode="integer")){
      #Extract sibstructure as proportions
      sibstructureCustom <- as.matrix(prop.table(table(sibstructureCustom)))
      sibstructureCustom <- cbind(as.numeric(rownames(sibstructureCustom)), sibstructureCustom[,1])
    }

    if(ncol(sibstructureCustom)!=2){ #Check that the format is correct
      stop("Information passed to the sibstructureCustom argument does not correspond with an expected format")
    }
    
  } else {
    sibstructureCustom <- NULL
  }
  
  #Store the various sibstructures in a list.
  strucs<- list(sibstructureUK,sibstructureNS,sibstructureCustom)
  
  message(paste0("A total of ", sims, " simulation(s) will be performed. Indications of progress across the ",sims," top-level loop will be returned."))
  #Loop across penetrance values
  for(l in 1:sims){
    ######
    ### Derive families per sibstructure
    ######
    if(1 %in% which_sibstructure){
      #Define seed for the randomisation
      #Generate population of families resembling UK structure
      set.seed(sim_seed[l])
      sibshipsUK <-sample(sibstructureUK[,1], numfamilies,replace=TRUE,prob=sibstructureUK[,2])
    } else {
      sibshipsUK <- NULL
    }
    
    if(2 %in% which_sibstructure){
      #Generate population of families resembling NS structure
      set.seed(sim_seed[l])
      sibshipsNS <-sample(sibstructureNS[,1], numfamilies,replace=TRUE,prob=sibstructureNS[,2])
    } else{
      sibshipsNS <- NULL
    }
    
    if(3 %in% which_sibstructure){
      #Generate a custom sibstructure
      set.seed(sim_seed[l])
      sibshipsCustom <-sample(sibstructureCustom[,1], numfamilies,replace=TRUE,prob=sibstructureCustom[,2])
    } else {
      sibshipsCustom <- NULL
    }
    
    #Store populations generated in column matrix
    sibships <- matrix(c(sibshipsUK,sibshipsNS,sibshipsCustom),
                       nrow=numfamilies)
    
    
    #Loop across simulated population structure
    for(p in 1:ncol(sibships)){
      #Model the poisson-tuned correction and correction tailored to specified sibstructure from strucs
      sibstructure <- list(NULL,strucs[[which_sibstructure[p]]])
      pop <- which_sibstructure[p] #label the population structure being tested
      
      currPop<- sibships[,p]
      
      obsN <- sum(currPop)/numfamilies #Identify mean N value in sample

      ######
      ### Generate and extract variant and non-variant families, 
      ######
      
      #Generate variant cohort: varChars for variant characteristics, genFamily for family wise results
      if(onsetRateDiff<1){
        #If onsetRateDiff<1 then onset time for comparator group is more compressed, flip the diff and supply this to var.Char_comparison calculation 
        var.Char_f <- varChars(f=f[l],g=g,numsteps=numsteps,onsetRateDiff=1)
      } else {
        #Otherwise, pass diff to the var.Char_f calculation. if onsetRateDiff==1, has no influence
        var.Char_f <- varChars(f=f[l],g=g,numsteps=numsteps,onsetRateDiff=onsetRateDiff)
      }
      
      cohort_var <- currPop[1:numfamilies_var] #Select families
      cohort_var <- lapply(cohort_var,genFamily,group="var",var.Char=var.Char_f,eldestAt0=eldestAt0,stepHazard=stepHazard)
        cohort_var<- cbind("var", bindRows(cohort_var))
        
      if(f_cohort_only==FALSE){
        #Repeat for novar group, generating a variety of variants as indicated by f_compare
        cohort_inp <- currPop[(numfamilies_var+1):(numfamilies_var+numfamilies_novar)]
        chunksize <- length(cohort_inp)/length(f_compare)
        cohort_novar <- NULL
        for(i in 1:length(f_compare)){
          if(onsetRateDiff<1){ #Apply varChars to compute onset characteristics for each f_compare[i]
            #If the comparator group is more compressed, flip the diff and supply this to var.Char_comparison calculation 
            var.Char_comparison <- varChars(f=f_compare[i],g=g,numsteps=numsteps,onsetRateDiff=1/onsetRateDiff)
          } else {
            var.Char_comparison <- varChars(f=f_compare[i],g=g,numsteps=numsteps,onsetRateDiff=1)
          }
          #Obtain subset of no variant cohort (who have various different variants), and generate assignments
          chunk <- cohort_inp[1:chunksize]
          temp <- lapply(chunk,genFamily,group="var",var.Char=var.Char_comparison,eldestAt0=eldestAt0,stepHazard=stepHazard) #Note that group is also denoted var. They have a comparison variant
          temp <- cbind("novar",bindRows(temp))
          
          cohort_novar<- bindRows(list(cohort_novar,temp),isMatrix = TRUE,silent=TRUE) #Save the result. #Use silent because first novar is expectedly NULL
          cohort_inp<- cohort_inp[-(1:chunksize)] #Drop the chunk used from input
        }
        
        
      #Unify the cohorts by passing a list to bindRows; this will drop any if they are null
      famchars <- bindRows(list(cohort_var,cohort_novar),isMatrix = TRUE)
      } else if(f_cohort_only){
        famchars<- cohort_var
      }
      colnames(famchars) <- c("set_group",paste0("t_",2:ncol(famchars)-2)) 
      
      if(f_cohort_only==FALSE){
        ######
        ### Unless specified, set weighting factors for all ADPenetrance runs, based on disease characteristics at final time
        ######
        
        #Number affected
        NAff <- length(which(famchars[,ncol(famchars)] %in% c("Familial","Sporadic")))
        PA=NAff/nrow(famchars)
        
        #Number familial
        NFam <- length(which(famchars[,ncol(famchars)] %in% c("Familial")))
        PF= NFam/NAff

      } else if (f_cohort_only) {
        #Required for output
        PF=NA
        PA=NA
      }
      
      #Loop across famchars columns (time points)
      for(n in 2:ncol(famchars)){
        # #Retain for 
        # cat("Table famchars:\n")
        # print(table(famchars[,n],famchars[,"set_group"]))
        
        if(f_cohort_only==FALSE){          
          
        #Determine the proportion of sample with familial sporadic and unaffected assignments, store in vector
        obsed <- c(
          MF=length(which(famchars[,n]=="Familial" & famchars[,"set_group"]=="var"))/length(which(famchars[,n]=="Familial")),
          MS=length(which(famchars[,n]=="Sporadic" & famchars[,"set_group"]=="var"))/length(which(famchars[,n]=="Sporadic")),
          MU=length(which(famchars[,n]=="Unaffected" & famchars[,"set_group"]=="var"))/length(which(famchars[,n]=="Unaffected"))
        )
        
        } else if(f_cohort_only){
          #Determine the proportion of sample with familial sporadic and unaffected assignments, store in vector
          obsed <- c(
            length(which(famchars[,n]=="Familial"))/numfamilies,
            length(which(famchars[,n]=="Sporadic"))/numfamilies,
            length(which(famchars[,n]=="Unaffected"))/numfamilies
          )
        }
        
        #If a state is now 0 or NaN, substitute nominally small freq to avoid any division by 0 within adpenetrance
        #NaN may occur if division by 0 when a state isn't present for any cohort
        #0 may occur if the state isn't present in the f cohort
        if(any(is.nan(obsed))){
          obsed[which(is.nan(obsed))] <- 1e-10
        }
        if(any(obsed==0)){
          obsed[which(obsed==0)] <- 1e-10
        }
        
        #### Population has been defined - now estimate penetrance and produce summary
        #seed     =     seed used in randomisation
        #f        =     true penetrance
        #famchars = df containing full population
        #obsed stores disease state rates across population
        
        #Run the adpenetrance analysis for all disease state combinations
        for(s in 1:length(states)){
          # #Retain for debugging
          # print(paste0("Current states loop is: ", s))
          #Run Adpenetrance for the states modelled at the given f value
          #loop across sibstructure which contains NULL in element 1, for a poisson distribution, and then the correct, tailored dist
          adpen_estimate<- list(poisson=NULL,tailored=NULL)
          if(f_cohort_only==FALSE){
            if(states[s]=="fsu") {
              stateval <- 1
              for(t in 1:2){
                suppressMessages(
                  adpen_estimate[[t]] <- adpenetrance(N=obsN,
                                                      MF=obsed["MF"], 
                                                      MS=obsed["MS"], 
                                                      MU=obsed["MU"], 
                                                      PA=PA, PF=PF,
                                                      define_sibstructure = sibstructure[[t]],
                                                      useG = g
                  )
                )
              }
            } else if(states[s]=="fs") {
              
              #If both MF and MS are 1e-10 then the actual variant frequency is NaN or 0 at this point and these nominally small frequencies have been 'substituted in'.
              #When sampling F and S only, RX=F*W_F/(F*W_F+S*W_S). Thus, even the nominally small values will cause inappropriate penetrance estimates.
              #Therefore, for this state combination only, return MF to 0
              MFsubto0 <- FALSE #Flag variable to indicate substitution is performed
              if(obsed["MF"]==1e-10 && obsed["MS"]==1e-10){
                obsed["MF"] <- 0
                MFsubto0 <- TRUE
              }
              
              stateval <- 2
              for(t in 1:2){
                suppressMessages(
                  adpen_estimate[[t]] <- adpenetrance(N=obsN,
                                                      MF=obsed["MF"], 
                                                      MS=obsed["MS"], 
                                                      PF=PF,
                                                      define_sibstructure = sibstructure[[t]],
                                                      useG = g
                  )
                )
              }
              
              if(MFsubto0){ #If substitution has been triggered, return to 1e-10 for subsequent state combinations
                obsed["MF"] <- 1e-10
              }
              
            } else if(states[s]=="fu") {
              stateval <- 3
              for(t in 1:2){
                suppressMessages(
                  adpen_estimate[[t]] <- adpenetrance(N=obsN,
                                                      MF=obsed["MF"],
                                                      MU=obsed["MU"], 
                                                      PA=PA, PF=PF,
                                                      define_sibstructure = sibstructure[[t]],
                                                      useG = g
                  )
                )
              }
            } else if(states[s]=="su") {
              stateval <- 4
              for(t in 1:2){
                suppressMessages(
                  adpen_estimate[[t]] <- adpenetrance(N=obsN,
                                                      MS=obsed["MS"], 
                                                      MU=obsed["MU"], 
                                                      PA=PA,
                                                      PF=PF,
                                                      define_sibstructure = sibstructure[[t]],
                                                      useG = g
                  )
                )
              }
            } else if(states[s]=="au") {
              stateval <- 5
              for(t in 1:2){
                suppressMessages(
                  adpen_estimate[[t]] <- adpenetrance(N=obsN,
                                                      MA=(obsed["MF"]*PF)+(obsed["MS"]*(1-PF)),
                                                      MU=obsed["MU"], 
                                                      PA=PA,
                                                      define_sibstructure = sibstructure[[t]],
                                                      useG = g
                  )
                )
              }
            }
          } else if(f_cohort_only){
            #Calculate RX according to states represented in 'states' object:
            #Obsed[1] = familial
            #Obsed[2] = sporadic
            #Obsed[3] = unaffected
            if(states[s]=="fsu") {
              stateval <- 1
              RX<- obsed[1] #No division needed because all states are represented
            } else if(states[s]=="fs") {
              stateval <- 2
              
              #If both obsed 1 and 2 1e-10 then the actual state frequency is 0 at this point and these nominally small frequencies have been 'substituted in'.
              #When sampling F and S only, RX=F/(F+S). Thus, even the nominally small values will cause inappropriate penetrance estimates.
              #Therefore, for this state combination only, return obsed[1] to 0
              MFsubto0 <- FALSE #Flag variable to indicate when substitution is performed
              if(obsed[1]==1e-10 && obsed[2]==1e-10){
                obsed[1] <- 0
                MFsubto0 <- TRUE
              }
              
              #Calculate RX
              RX<- obsed[1]/(obsed[1]+obsed[2])
              
              if(MFsubto0){ #If substitution has been triggered, return to 1e-10 for subsequent state combinations
                obsed[1] <- 1e-10
              }
            } else if(states[s]=="fu") {
              stateval <- 3
              RX<- obsed[1]/(obsed[1]+obsed[3])
            } else if(states[s]=="su") {
              stateval <- 4
              RX<- obsed[2]/(obsed[2]+obsed[3])
            } else if(states[s]=="au") {
              stateval <- 5
              RX<- (obsed[1]+obsed[2])/(obsed[1]+obsed[2]+obsed[3])
            }

            for(t in 1:2){
              suppressMessages(
                adpen_estimate[[t]] <- adpenetrance(N=obsN,
                                                    RX=RX,states=states[s],
                                                    define_sibstructure = sibstructure[[t]],
                                                    useG = g
                )
              )
            }
            
          } #End onecohortonly
          
          #Summary structure for data:
          #c("Population","Tailoring","Seed","mean_N","MF","MS","MU","PF","PA","RX","modif","States.Modelled","TruePenetrance","UncorrectedPen", "CorrectedPen")
          
          #Store results for poisson and tailored in the summary list - define population name and that step-4 simulated sibs follow Poisson distribution
          for(t in 1:2){
            #Store results for poisson and tailored in the summary matrices - define population name and that step-4 simulated sibs follow Poisson distribution
            statetemp[[t]] <- c(pop,t,sim_seed[l],obsN,obsed,PF,PA,
                                adpen_estimate[[t]]$output[1,"Estimate"],n-2,stateval,f[l], adpen_estimate[[t]]$output[3,], adpen_estimate[[t]]$output[4,])
            
          }
          
        result <- rbind(result,statetemp[[1]],statetemp[[2]]) #Combine temporary vectors with other results
        } #End states loop (s)
      } #End time loop (n)
    } #End population structure loop (p)
    print(l) #print values of l for indications of progress
  } #End loop across all test f values (l)
  
  #Dataset is generated and stored in the 'results' object - now assign colnames and extract useful information 
  colnames(result) <- output_colnames
  
  #Calculate the difference between true and estimated penetrance
  CorrMinusTrue <- result[,"CorrectedPen"]-result[,"TruePenetrance"]
  
  #Calculate the difference between true and uncorrected penetrance
  UncorrMinusTrue <- result[,"UncorrectedPen"]-result[,"TruePenetrance"]
  
  #Combine with difference column and convert into data.frame
  sum_diff<- as.data.frame(cbind(result,CorrMinusTrue,UncorrMinusTrue))
  
  #Average the repeated seeds to get the mean difference in True penetrance minus estimated across states modelled at each value of true penetrance
  SimFunc <- ddply(sum_diff,.(Population,Tailoring,TruePenetrance,States.Modelled,modif),summarise,CorrMinusTrue=mean(CorrMinusTrue),UncorrMinusTrue=mean(UncorrMinusTrue))
  
  #Convert states modelled column into factors with appropriate names
  SimFunc$States.Modelled <- factor(SimFunc$States.Modelled,levels= c(1:5),labels= c("F, S, U",
                                                                                     "F, S",
                                                                                     "F, U",
                                                                                     "S, U",
                                                                                     "A, U"))
  
  SimFunc$Population <- factor(SimFunc$Population,levels= c(1:3),labels= c("UK (1974)",
                                                                           "Next Steps",
                                                                           nameSibstructureCustom))
  
  SimFunc$Tailoring <- factor(SimFunc$Tailoring,levels= c(1:3),labels= sim_labs <- c("Poisson",
                                                                                     "Tailored",
                                                                                     "No correction"))
  
  #Melt so that the corrected and uncorrected penentrance estiamtes are in the same column
  mSimFunc<- reshape2::melt(SimFunc,measure.vars=c("CorrMinusTrue","UncorrMinusTrue"))
  
  #Unadjusted estimates are duplicated (i.e. Poisson and Tailored show the same, therefore remove any 'tailored' uncorrected rows)
  mSimFunc <- mSimFunc[-which(mSimFunc$variable=="UncorrMinusTrue" & mSimFunc$Tailoring=="Tailored"),] 
  
  #Then relabel uncorrected into a third 'tailoring' category
  mSimFunc$Tailoring[which(mSimFunc$variable=="UncorrMinusTrue")] <- "No correction"
  
  fig <- ggplot(data=mSimFunc, aes(modif,value,colour=factor(TruePenetrance)))+
    geom_line()+
    facet_grid(States.Modelled~Population+Tailoring)+
    ylab("Difference from lifetime penetrance (Estimate - True lifetime penetrance)")+
    xlab("Time from first sampling")+
    geom_hline(yintercept=0,lty=2,alpha=0.75)+
    scale_colour_manual(name = 'True Lifetime Penetrance',
                        values =c(#"#e3d16e", #0.1
                          "#ff7943", #0.25
                          "#be3044", #0.5
                          "#796c5d", #0.75
                          "#202547"), #1
    )+  #https://jaredhuling.org/jcolors/ colour palette pal5
    scale_x_continuous(breaks=function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+ #Force integer x axis
    theme(legend.position = "bottom",
          panel.background = element_rect(fill="white"),
          panel.grid = element_line(colour = "gray90"),
          strip.background = element_rect(fill="gray94"),
          panel.border = element_rect(fill=NA, colour="gray94")
    )
  
  
  
  #Drop the 'variable' column since it is meaningless, a product of using melt which is superseded by the information provided in the "Tailoring" column
  #Rename 'value' column to something more interpretable
  mSimFunc <- mSimFunc[,-6]
  colnames(mSimFunc)[6] <- "Estimate.difference"
  
  ## Remove 'ggsaving' functionality from within function; instead return output as a list containing the figure and input data
    # plotHeight<-length(states)*40            #Scale height of the plot by number of vertical panels expected
    # if(plotHeight<150){plotHeight <- 150}    #However, also set a minimum plot height
  
    # Save the output plot, and input data
    # ggsave(filename=outfile,plot=fig,device='pdf',units="mm",width=210,height=plotHeight,dpi=300)
    # write.table(mSimFunc,file=gsub(".pdf",".tsv",outfile),sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE) 
  
  #Return the summary data.frame that is passed to ggplot
  return(list(ggfigure=fig,results=mSimFunc))
               
}

### Print message if dependent packages are not installed [Move within function??]
req_packages <- c("plyr", "ggplot2", "reshape2")
pkg_missing <- req_packages[!req_packages %in% installed.packages()[ , "Package"]]
if(length(pkg_missing)>0){ 
  warning("The simADPenetrance simulation function is dependent upon on the several packages which were not identified when calling installed.packages().\nPlease install and load the following before using the function: ",paste0(pkg_missing,collapse=", "),"\n")
}
rm(pkg_missing,req_packages)
