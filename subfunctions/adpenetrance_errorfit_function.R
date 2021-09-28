#This script produces a function for generating a simulated population representative of the sample used in penetrance estimation
  #and predicting the error in penetrance estimates of a given value, as determined under a polynomial regression model fitted to the simulated population.

#By default, the simulated sibship distribution follows the Poisson distribution, with lambda defined by mean sibship size in the real sample data.
  #The user of the main adpenetrance function can opt to supply information about the sibship structure of their sample
  #and the simulated sibship distribution produced be generated to match this structure.

#The adpenetrance.unadjusted function must be loaded in order for this script to function
  #This function can be loaded from the "penetrance_unadjusted_function.R" script available within our github repository


#Important terms for this function
  #states = the disease states represented in the main adpenetrance function calculation, options are:
    #'fsu' = familial, sporadic, and unaffected    
    #'fs'= familial and sporadic
    #'fu' = familial and unaffected
    #'su' = sporadic and unaffected
    #'au' = affected and unaffected
  #setmean = the value of N used in the main adpenetrance function.
  #define_sibstructure = optionally specify the structure of sibships sampled
    #to allow tailoring of sibship structures generated in the simulated population to match the real sample datat
    #This allows more tailored adjustment of errors in penetrance values and giving a more precise penetrance estimate


#Define function
adpenetrance.errorfit <- function(states,setmean,samp_size=90000,seed=24,define_sibstructure=NULL){

#Define simulated sibships by Poison distribution if no data are given in 'define_sibstructure'
  #If data are given as a 2 column object, col1 represents siblevels and col2 represents proportions of population at each level  
  #If data do not have 2 columns, expect a vector of sibship sizes from which the proportions of each sib-size level can be determined
  if(!is.null(define_sibstructure)){
    if(is.matrix(define_sibstructure) && ncol(define_sibstructure)==2){
      #define_sibstructure[,1] - sibship levels in population
      #define_sibstructure[,2] - sibship weightings (e.g. proportions of each sib size number of each sibships)
      
      set.seed(seed)
      sibships <- sample(as.numeric(define_sibstructure[,1]), 90000,replace=TRUE,prob=as.numeric(define_sibstructure[,2]))
    } else {
      #Extract sibship proportions specified in data
      size_probs<- prop.table(table(define_sibstructure))
      #Define sibship levels by:  as.numeric(names(size_probs))
      #Define probabilities by:   unname(size_probs)
      
      #Generate simulation sample based on data 
      set.seed(seed)
      sibships <- sample(as.numeric(names(size_probs)), 90000,replace=TRUE,prob=unname(size_probs))
    }
  } else {
    #Generate simulated population with Poisson distribution
    set.seed(seed)
    sibships <- rpois(n=samp_size,lambda = setmean)
  }

#Identify and sort unique sibship levels in sample
  sim_N<- sort(unique(sibships))
  
##Identify mean N value in sample
  obsN <- sum(sibships)/length(sibships)

#Generate sequence of 25 penetrance values between 0 and 1 - curves are stable across sample sizes
  sim_f=seq(from=0.01,to=1,length.out = 25)

#Prepare empty matrix for use in loop, storing penetrance outputs
  Preds <- matrix(NA,ncol=4,nrow=length(sim_f))
  colnames(Preds) <- c("SeedNo.","TrueF", "PredF","diffs")
  
#Prepare empty matrix for storing the sibship size and assigned disease state
  famchars <- matrix(c(sibships,rep(NA,length(sibships))),ncol=2,nrow=length(sibships))
  
#Name states to be assigned in $states column
  all_states <- c("Familial", "Sporadic", "Unaffected") 
  
#Loop across all values of penetrance
for(m in 1:length(sim_f)){
  
  ###Calculate the probability of each disease state expected each value of sim_N with penetrance f
  Punasc = (1-sim_f[m])*(1-sim_f[m]/2)^sim_N                                                            #Probability unaffected 
  Pspor  = sim_f[m]*(1-sim_f[m]/2)^sim_N+sim_N*(sim_f[m]/2)*(1-sim_f[m]/2)^(sim_N-1)*(1-sim_f[m])       #Probability sporadic/simplex
  Pfam   = 1-((1-sim_f[m]/2)^sim_N+sim_N*(sim_f[m]/2)*((1-sim_f[m]/2)^(sim_N-1))*(1-sim_f[m]))          #Probability familial (= 1 - Punasc - Pspor)
  
  #Combine the calculated probabilities into a matrix, with row names defined per the value of sim_N
  Probabilities <- matrix(c(Pfam, Pspor, Punasc),
                             nrow=length(sim_N))

  #Loop across sibship sizes in sample
  for(i in 1:nrow(Probabilities)){
    #Identify which famchars sib_size has the value of sibship at index i 
    index <- which(famchars[,1]==sim_N[i])
    
    #Pseudo-randomly assign disease states to families according to probabilities expected at each sib size
    set.seed(seed)
    famchars[index,2] <- sample(all_states,size=length(index), replace = TRUE,prob=Probabilities[i,])
  }
  
  #Determine familial, sporadic, and unaffected rates across the whole population
    obsed <- c(
      length(which(famchars[,2]=="Familial"))/length(famchars[,2]),
      length(which(famchars[,2]=="Sporadic"))/length(famchars[,2]),
      length(which(famchars[,2]=="Unaffected"))/length(famchars[,2])
    )
  
  #Calculate observed rate of state X (obsRX) according to states represented in 'states' object:
    #Obsed[1] = familial
    #Obsed[2] = sporadic
    #Obsed[3] = unaffected
  if(states=="fsu") {
    obsRX<- obsed[1] #No division needed because all states are represented
  } else if(states=="fs") {
    obsRX<- obsed[1]/(obsed[1]+obsed[2])
  } else if(states=="fu") {
    obsRX<- obsed[1]/(obsed[1]+obsed[3])
  } else if(states=="su") {
    obsRX<- obsed[2]/(obsed[2]+obsed[3])
  } else if(states=="au") {
    obsRX<- (obsed[1]+obsed[2])/(obsed[1]+obsed[2]+obsed[3])
  }
    
  #Run adpenetrance.unadjusted
  adpen_out <- adpenetrance.unadjusted(N=obsN,RX=obsRX,states=states)

  #Store important information in the output object Preds
  Preds[m,] <- c(seed,sim_f[m], adpen_out$output[3],sim_f[m]-adpen_out$output[3])
  
} #End loop across all values of sim_f


###############################
#  Fit curve of difference    #
###############################
x <- Preds[,3]     #RHS of eqn - specify predicted values as x axis of curve
y <- Preds[,4]     #LHS of eqn - specify difference between ground truth and each penetrance estimate

#Find best fitting polynomial for a the difference curve
polyfit <- function(j) x <- AIC(lm(y~poly(x,j)))                  #Define function to optimise
best<- as.integer(optimize(polyfit,interval = c(1,5))$minimum)    #optimise function

#Fit the best model
fitbest <- lm(y~poly(x,best))

#Return the best fitting object
return(list(fitbest=fitbest,Preds=Preds))

} #End function

########################################################################