#This script generates the adpenetrance function for calculating genetic penetrance in autosomal dominant traits.
#The function relies upon the sub-functions adpenetrance.errorfit and adpenetrance.unadjusted,
  #which are downladed from GitHub when running this script (See https://github.com/thomaspspargo/adpenetrance)

#The README.md documentation in the associated GitHub repository outlines how the adpenetrance function should be used
  #https://github.com/thomaspspargo/adpenetrance

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
adpenetrance <- function(N, MF=0, MS=0, MA=0, MU=0, PA=0, PF=0, MF_SE=0, MS_SE=0, MU_SE=0, MA_SE=0, Zout=1.96,
                         useG=NA,
                         RX=0, RX_SE=0, states="none",define_sibstructure=NULL){
  
  ######
  ### Run ADPenetrance unadjusted function for steps 1-3, passing across relevant inputs
  ######
  unadjusted.out<- adpenetrance.unadjusted(N=N, MF=MF, MS=MS, MA=MA, MU=MU, PA=PA, PF=PF, MF_SE=MF_SE, MS_SE=MS_SE, MU_SE=MU_SE, MA_SE=MA_SE, Zout=Zout,
                                           useG=useG,
                                           RX=RX, RX_SE=RX_SE, states=states)
  
  ######
  ### Perform step 4 Adjustment
  ######
  
  #Generate error curve for sample, according to the specified N, the states combination returned by the unadjusted function, and whether upload_sibships or define_sibships have been defined
  errorfit <- adpenetrance.errorfit(states=unadjusted.out$states, setmean=N,samp_size=90000,seed=24,define_sibstructure=define_sibstructure,useG = useG)
  
  #Prepare matrix to store adjusted f values
  Corr_f <- matrix(NA,ncol=4,nrow=1)
  colnames(Corr_f) <- c("Lower CI", "Estimate", "Upper CI","Standard error")
  rownames(Corr_f) <- c("Adjusted Penetrance")
  
  #Predict the error in f for the point estimate, based on the value of f obtained, using the correction curve
  correctionVal<- unname(
    predict(errorfit$fitbest,data.frame(x=unadjusted.out$output["Unadjusted Penetrance","Estimate"]))
  )
  #Adjust penetrance estimate obtained by the value obtained. Store in Corr_f matrix
  Corr_f[,"Estimate"] <- unadjusted.out$output["Unadjusted Penetrance","Estimate"]+correctionVal

  #If error propagation has been included, above process at confidence interval bounds
  if(ncol(unadjusted.out$output)==4) {
    #As above, predict error in f as above for the Lower CI bound then adjust by value obtained
    correctionVal<- unname(
      predict(errorfit$fitbest,data.frame(x=unadjusted.out$output["Unadjusted Penetrance","Lower CI"]))
    )
    Corr_f[,"Lower CI"] <- unadjusted.out$output["Unadjusted Penetrance","Lower CI"]+correctionVal
    
    #As above, predict error in f as above for the Upper CI bound then adjust by value obtained
    correctionVal<- unname(
      predict(errorfit$fitbest,data.frame(x=unadjusted.out$output["Unadjusted Penetrance","Upper CI"]))
    )
    Corr_f[,"Upper CI"] <- unadjusted.out$output["Unadjusted Penetrance","Upper CI"]+correctionVal
  }
  
  #Truncate adjusted estimates to fit within [0,1]
  Corr_f[which(Corr_f<0)] <- 0
  Corr_f[which(Corr_f>1)] <- 1
  
  #Append adjusted elements to unadjusted.out, subset columns if only point estimate is being made
  if(ncol(unadjusted.out$output)==4) {
    output <- rbind(unadjusted.out$output,Corr_f)
  } else {
    output <- rbind(unadjusted.out$output,Corr_f[,"Estimate"])
  }
  rownames(output)[4] <- "Adjusted Penetrance" #add new rowname for adjusted penetrance
  
  return(list(output=output,states=states,errorfit=summary(errorfit$fitbest)))
  
}
