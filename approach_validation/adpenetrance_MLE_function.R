#This script is associated with the following GitHub repository: https://github.com/thomaspspargo/adpenetrance.
#It generates the adpenetrance function for calculating genetic penetrance in autosomal dominant traits,
  #including an additional step to allow Maximum Likelihood Estimation of unadjusted penetrance estimates.
    #This is used for validation of the lookup table approach primarily employed by the method,
    #and is to be called in the adpenetrance_MLE_validation.R script contained within the associated github repository.

  #The function relies upon the sub-functions adpenetrance.errorfit and adpenetrance.unadjusted, which are downladed from GitHub when running this script

#--------------------------------------------#
#       Penetrance calculator function       #
#--------------------------------------------#

#First, if not already loaded in the environment, download the adpenetrance functions from GitHub
if(!"adpenetrance" %in% ls()){
  source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/adpenetrance_function.R")
}



#Second, specify main adpenetrance function
adpenetrance.MLE <- function(N, MF = NA, MS = NA, MA = NA, MU = NA, PA = NA, PF = NA, MF_SE = NA, MS_SE = NA, MU_SE = NA, MA_SE = NA, Zout = 1.96,
                             useG = 0,
                             RX = NA, RX_SE = NA, states = "none", define_sibstructure = NULL, include_MLE=TRUE){
  
  
  #------------------------------------------------#
  #   Pass all standard arguments to adpenetrance  #
  #------------------------------------------------#
  adpenetrance.out <- adpenetrance(N=N, MF = MF, MS = MS, MA = MA, MU = MU, PA = PA, PF = PF, MF_SE = MF_SE, MS_SE = MS_SE, MU_SE = MU_SE, MA_SE = MA_SE, Zout = Zout,                
    useG = useG, RX = RX, RX_SE = RX_SE, states = states, define_sibstructure = define_sibstructure) 
  
  #------------------------------------------------#
  #               Prep and run MLE step            #
  #------------------------------------------------#
  
  if(isTRUE(include_MLE)){
    
    states=adpenetrance.out$states
    RXinp=adpenetrance.out$output[grepl("Observed",rownames(adpenetrance.out$output)),"Estimate"]
    
    #ADD MLE STEP
    #Define likelihood function, f (penetrance) to be maximised
    #MLE to be performed using nlm, so negML function returns dbinom()*-1 
    negML <- function(f,useG,RXinp){
      g=useG
      #Extended Al-Chalabi & Lewis (2011) model equations (If g=0, synonymous with original equations)
      Punasc = (1-f)*(1-f/2-g/2)^N*(1-g) #Probability unaffected
      Pspor  = f*(1-f/2-g/2)^N*(1-g)+
        (1-f)*N*(f/2+g/2)*(1-f/2-g/2)^(N-1)*(1-g)+
        (1-f)*(1-f/2-g/2)^N*g #Probability sporadic/simplex
      Pfam   = 1-(Punasc+Pspor) #Probability familial (= 1 - Punasc - Pspor)
      
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
      RXml = round(RXinp,7)*trials #Adjust the ObsProbX proportion by the number of trials - rounding to guarantee integer
      
      set.seed(24)
      dbinom(RXml,trials,r)*-1       #Negative likelihood estimate: dbinom(number of events,number of trials, probabilities)
    }
    
    #Prepare results table for MLE estimate
    MLresults <-matrix(NA,ncol=3,nrow=3)
    colnames(MLresults) <- c("Lower CI","Estimate", "Upper CI") #Assign column names
    rownames(MLresults) <- c("Unadjusted Penetrance", "Adjusted Penetrance","NLM convergence") #Assign row names
    
    #Use nlm function to minimise negML, finding ML penetrance estimate
    #Starting value is set by the lookup table unadjusted estimate, close to true penetrance value...
    MLestimate <- nlm(negML,adpenetrance.out$output["Unadjusted Penetrance","Estimate"],useG=useG,RXinp=RXinp,stepmax=0.5)
    
    MLresults[1,2] <- MLestimate$estimate #Uncorrected ML estimate
    correctionVal<- unname(predict(adpenetrance.out$errorfit,data.frame(x=MLresults[1,2]))) #Predict the error in f for the MLR estimate, based on the value of f obtained, using the correction curve
    MLresults[2,2] <- MLresults[1,2]+correctionVal #Corrected ML estimate
    MLresults[3,2] <- MLestimate$code #Check for convergence of MLE model
    
    if(ncol(adpenetrance.out$output)==4) { #If error terms given, at CI bounds of penetrance estimate
      MLlower <- nlm(negML,adpenetrance.out$output["Unadjusted Penetrance","Lower CI"],useG=useG,RXinp=RXinp,stepmax=0.5)
      MLupper <- nlm(negML,adpenetrance.out$output["Unadjusted Penetrance","Upper CI"],useG=useG,RXinp=RXinp,stepmax=0.5)
      
      #Lower results
      MLresults[1,1] <- MLlower$estimate #Uncorrected ML estimate
      correctionVal<- unname(predict(adpenetrance.out$errorfit,data.frame(x=MLresults[1,1]))) 
      MLresults[2,1] <- MLresults[1,1]+correctionVal #Corrected ML estimate
      MLresults[3,1] <- MLestimate$code #Check for convergence of MLE model
      
      #Upper results
      MLresults[1,3] <- MLupper$estimate #Uncorrected ML estimate
      correctionVal<- unname(predict(adpenetrance.out$errorfit,data.frame(x=MLresults[1,3])))
      MLresults[2,3] <- MLresults[1,3]+correctionVal #Corrected ML estimate
      MLresults[3,3] <- MLestimate$code #Check for convergence of MLE model
    }
    
    return(list(default=adpenetrance.out,MLE.extension=MLresults))
  } else {
    
    return(adpenetrance.out)
  }
  
}
