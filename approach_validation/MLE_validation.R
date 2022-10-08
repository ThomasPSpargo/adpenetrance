#This script calls the adapted adpenetrance function which makes additional Maximum Likelihood Estimates (MLE) for unadjusted penetrance
  #The adapted function is applied within this script to the case studies presented in Table 2 of the manuscript associated with this GitHub repository
  #The unadjusted penetrance estimates obtained via the lookup and MLE approaches are compared and results are returned as a .csv file
  #The MLE is checked for convergence in all instances, as defined within the Non-Linear Minimisation approach employed (see ?nlm)
    #Values of 1 or 2 indicate convergence and were observed in all studies

#First, if not already loaded in the environment, download the adpenetrance functions from GitHub
if(!"adpenetrance" %in% ls()){
  source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/adpenetrance_function.R")
}
#Load penetrance function with additional MLE step and subfunctions
source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/approach_validation/adpenetrance_MLE_function.R")

#Extra functions
#Source function for calculating residual risk g
source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/getResidualRisk.R")

### Quick function for calculating standard error for proportion based on variant counts
#count is the number of variants in the sample
#size is the sample size
se_p <- function(count,size){
  p <- count/size
  n <- size
  
  sqrt(p*(1-p)/n)
}

#Import dataset
Cases <- read.csv("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/case_studies/case_data.csv",header=TRUE)

#Select the populations presented in Table 2, stratify by data supply method (variant counts or variant frequency estimates with CIs)
Cases_T2_counts <- Cases[c(41,29:32),]
Cases_T2_freqs  <- Cases[c(36:39),]


#Loop across Case Studies 1 and 2 (in which variant counts were provided)
Cases_T2_counts_out <- NULL
for(k in 1:length(Cases_T2_counts$Population_sampled)){
  
  #Compute the remaining disease for a population member who does not harbour the variant
  residualRiskG<- getResidualRisk(MF=Cases_T2_counts[k,"Var_fam"]/Cases_T2_counts[k,"N_familial"],
                                  MS=Cases_T2_counts[k,"Var_spor"]/Cases_T2_counts[k,"N_sporadic"],
                                  PF=Cases_T2_counts[k,"PF"],
                                  PA=1/Cases_T2_counts[k,"PA"])
  
  #FS
  #Every sample has data for the familial and sporadic states, run for every sample
  #ObsX is the familial/(familial+sporadic) rate
  cases_out <- adpenetrance.MLE(N=Cases_T2_counts[k,"TFR"],
                            MF=Cases_T2_counts[k,"Var_fam"]/Cases_T2_counts[k,"N_familial"],
                            MF_SE=se_p(Cases_T2_counts[k,"Var_fam"],Cases_T2_counts[k,"N_familial"]),
                            MS=Cases_T2_counts[k,"Var_spor"]/Cases_T2_counts[k,"N_sporadic"],
                            MS_SE=se_p(Cases_T2_counts[k,"Var_spor"],Cases_T2_counts[k,"N_sporadic"]),
                            #MU=Cases_T2_counts[k,"Var_unaf"]/Cases_T2_counts[k,"N_unaffected"],
                            #MU_SE=se_p(Cases_T2_counts[k,"Var_unaf"],Cases_T2_counts[k,"N_unaffected"]),
                            PF=Cases_T2_counts[k,"PF"],
                            #PA=1/Cases_T2_counts[k,"PA"]
                            useG = residualRiskG,
                            include_MLE = TRUE)
  
  #Matrix of output
  est_out_fs <- matrix(c(2,
                         cases_out$default$output[1,1], cases_out$default$output[1,2], cases_out$default$output[1,3],
                         cases_out$default$output[3,1], cases_out$default$output[3,2], cases_out$default$output[3,3],
                         cases_out$MLE[1,1],    cases_out$MLE[1,2],    cases_out$MLE[1,3],
                         cases_out$MLE[3,1],    cases_out$MLE[3,2],    cases_out$MLE[3,3],
                         cases_out$default$output[4,1], cases_out$default$output[4,2], cases_out$default$output[4,3]),nrow=1)
  
  
  
  #Run additional analyses for the Case 1 data
  if(k==1){
    
    #FSU
    #ObsX is the familial/(familial+sporadic+unaffected) rate
    cases_out <- adpenetrance.MLE(N=Cases_T2_counts[k,"TFR"],
                                  MF=Cases_T2_counts[k,"Var_fam"]/Cases_T2_counts[k,"N_familial"],
                                  MF_SE=se_p(Cases_T2_counts[k,"Var_fam"],Cases_T2_counts[k,"N_familial"]),
                                  MS=Cases_T2_counts[k,"Var_spor"]/Cases_T2_counts[k,"N_sporadic"],
                                  MS_SE=se_p(Cases_T2_counts[k,"Var_spor"],Cases_T2_counts[k,"N_sporadic"]),
                                  MU=Cases_T2_counts[k,"Var_unaf"]/Cases_T2_counts[k,"N_unaffected"],
                                  MU_SE=se_p(Cases_T2_counts[k,"Var_unaf"],Cases_T2_counts[k,"N_unaffected"]),
                                  PF=Cases_T2_counts[k,"PF"],
                                  PA=1/Cases_T2_counts[k,"PA"],
                                  useG = residualRiskG,
                                  include_MLE = TRUE)
    
    #Matrix of output
    est_out_fsu <- matrix(c(1,
                            cases_out$default$output[1,1], cases_out$default$output[1,2], cases_out$default$output[1,3],
                            cases_out$default$output[3,1], cases_out$default$output[3,2], cases_out$default$output[3,3],
                            cases_out$MLE[1,1],    cases_out$MLE[1,2],    cases_out$MLE[1,3],
                            cases_out$MLE[3,1],    cases_out$MLE[3,2],    cases_out$MLE[3,3],
                            cases_out$default$output[4,1], cases_out$default$output[4,2], cases_out$default$output[4,3]),nrow=1)
    
    
    
    
    #FU
    #ObsX is the familial/(familial+unaffected) rate
    cases_out <- adpenetrance.MLE(N=Cases_T2_counts[k,"TFR"],
                                  MF=Cases_T2_counts[k,"Var_fam"]/Cases_T2_counts[k,"N_familial"],
                                  MF_SE=se_p(Cases_T2_counts[k,"Var_fam"],Cases_T2_counts[k,"N_familial"]),
                                  #MS=Cases_T2_counts[k,"Var_spor"]/Cases_T2_counts[k,"N_sporadic"],
                                  #MS_SE=se_p(Cases_T2_counts[k,"Var_spor"],Cases_T2_counts[k,"N_sporadic"]),
                                  MU=Cases_T2_counts[k,"Var_unaf"]/Cases_T2_counts[k,"N_unaffected"],
                                  MU_SE=se_p(Cases_T2_counts[k,"Var_unaf"],Cases_T2_counts[k,"N_unaffected"]),
                                  PF=Cases_T2_counts[k,"PF"],
                                  PA=1/Cases_T2_counts[k,"PA"],
                                  useG = residualRiskG,
                                  include_MLE = TRUE)
    
    
    #Matrix of output
    est_out_fu <- matrix(c(3,
                           cases_out$default$output[1,1], cases_out$default$output[1,2], cases_out$default$output[1,3],
                           cases_out$default$output[3,1], cases_out$default$output[3,2], cases_out$default$output[3,3],
                           cases_out$MLE[1,1],    cases_out$MLE[1,2],    cases_out$MLE[1,3],
                           cases_out$MLE[3,1],    cases_out$MLE[3,2],    cases_out$MLE[3,3],
                           cases_out$default$output[4,1], cases_out$default$output[4,2], cases_out$default$output[4,3]),nrow=1)
    
    
    
    #SU
    #ObsX is the sporadic/(sporadic+unaffected) rate
    cases_out <- adpenetrance.MLE(N=Cases_T2_counts[k,"TFR"],
                                  #MF=Cases_T2_counts[k,"Var_fam"]/Cases_T2_counts[k,"N_familial"],
                                  #MF_SE=se_p(Cases_T2_counts[k,"Var_fam"],Cases_T2_counts[k,"N_familial"]),
                                  MS=Cases_T2_counts[k,"Var_spor"]/Cases_T2_counts[k,"N_sporadic"],
                                  MS_SE=se_p(Cases_T2_counts[k,"Var_spor"],Cases_T2_counts[k,"N_sporadic"]),
                                  MU=Cases_T2_counts[k,"Var_unaf"]/Cases_T2_counts[k,"N_unaffected"],
                                  MU_SE=se_p(Cases_T2_counts[k,"Var_unaf"],Cases_T2_counts[k,"N_unaffected"]),
                                  PF=Cases_T2_counts[k,"PF"],
                                  PA=1/Cases_T2_counts[k,"PA"],
                                  useG = residualRiskG,
                                  include_MLE = TRUE)
    
    #Matrix of output
    est_out_su <- matrix(c(4,
                           cases_out$default$output[1,1], cases_out$default$output[1,2], cases_out$default$output[1,3],
                           cases_out$default$output[3,1], cases_out$default$output[3,2], cases_out$default$output[3,3],
                           cases_out$MLE[1,1],    cases_out$MLE[1,2],    cases_out$MLE[1,3],
                           cases_out$MLE[3,1],    cases_out$MLE[3,2],    cases_out$MLE[3,3],
                           cases_out$default$output[4,1], cases_out$default$output[4,2], cases_out$default$output[4,3]),nrow=1)
    
    #Combine all the out objects, where all states have been modelled
    est_out <- rbind(est_out_fsu,est_out_fs,est_out_fu,est_out_su)
    
  } else {
    #If only fs modelled, take just that object
    est_out <- matrix(est_out_fs,nrow=1)
  }
  
  rownames(est_out) <- rep(Cases_T2_counts$Population_sampled[k],nrow(est_out))
  #Combine into an overall output object
  Cases_T2_counts_out <- rbind(Cases_T2_counts_out,est_out)
  
}







#Loop across Case Studies 3 and 4 (in which variant frequencies and confidence intervals were provided)
Cases_T2_freqs_out <- matrix(as.numeric(NA),ncol=16,nrow=4) #Prepare empty matrix
rownames(Cases_T2_freqs_out)<- Cases_T2_freqs$Population_sampled
for(k in 1:length(Cases_T2_freqs$Population_sampled)){
  
  #Compute the remaining disease for a population member who does not harbour the variant
  residualRiskG<- getResidualRisk(MF=Cases_T2_freqs[k,"Freq_fam"],
                                  MS=Cases_T2_freqs[k,"Freq_spor"],
                                  PF=Cases_T2_freqs[k,"PF"],
                                  PA=1/Cases_T2_freqs[k,"PA"])
  
  #ObsX is the familial/(familial+sporadic) rate
  cases_out <- adpenetrance.MLE(N=Cases_T2_freqs[k,"TFR"],
                            MF=Cases_T2_freqs[k,"Freq_fam"], MF_SE=(Cases_T2_freqs[k,"Freq_fam"]-Cases_T2_freqs[k,"Freq_fam_LCI"])/1.96,
                            MS=Cases_T2_freqs[k,"Freq_spor"], MS_SE=(Cases_T2_freqs[k,"Freq_spor"]-Cases_T2_freqs[k,"Freq_spor_LCI"])/1.96,
                            PF=Cases_T2_freqs[k,"PF"], Zout=1.96,
                            useG = residualRiskG,
                            include_MLE = TRUE)
  
  #Matrix of output
  Cases_T2_freqs_out[k,] <- c(2,
                               cases_out$default$output[1,1], cases_out$default$output[1,2], cases_out$default$output[1,3],
                               cases_out$default$output[3,1], cases_out$default$output[3,2], cases_out$default$output[3,3],
                               cases_out$MLE[1,1],    cases_out$MLE[1,2],    cases_out$MLE[1,3],
                               cases_out$MLE[3,1],    cases_out$MLE[3,2],    cases_out$MLE[3,3],
                               cases_out$default$output[4,1], cases_out$default$output[4,2], cases_out$default$output[4,3])
}

#Combine the results from Case Studies 1-4
Cases_T2_out <- rbind(Cases_T2_counts_out,Cases_T2_freqs_out)

colnames(Cases_T2_out) <- c("states_modelled","RX_lower","RX_est","RX_upper",
                          "Un_pen_lower","Un_pen_est","Un_pen_upper",
                          "MLE_pen_lower","MLE_pen_est","MLE_pen_upper",
                          "MLE_convergence_lower","MLE_convergence_est","MLE_convergence_upper",
                          "Adj_pen_lower","Adj_pen_est","Adj_pen_upper")


#Convert to data.frame and truncate the values to 3dp
sum_Table <- round(as.data.frame(Cases_T2_out),digits=3)

#Assign factor levels states modelled
sum_Table$states_modelled <- factor(sum_Table$states_modelled,levels=1:4,labels=c("fsu","fs","fu","su"))


#Summarise results in "estimate (LowerCI, UpperCI)" format
sums <- cbind(paste0(rownames(sum_Table)," (",sum_Table$states_modelled,")"),
              paste0(sum_Table$RX_est," (",sum_Table$RX_lower,", ",sum_Table$RX_upper,")"),
              paste0(sum_Table$Un_pen_est," (",sum_Table$Un_pen_lower,", ",sum_Table$Un_pen_upper,")"),
              paste0(sum_Table$MLE_pen_est," (",sum_Table$MLE_pen_lower,", ",sum_Table$MLE_pen_upper,")"),
              paste0(sum_Table$MLE_convergence_est," (",sum_Table$MLE_convergence_lower,", ",sum_Table$MLE_convergence_upper,")"),
              paste0(sum_Table$Adj_pen_est," (",sum_Table$Adj_pen_lower,", ",sum_Table$Adj_pen_upper,")")
)

#Assign colnames
colnames(sums) <- c("Analysis (states modelled)", "RXobs","lookup unadjusted pen", "MLE unadjusted penetrance", "MLE convergence (see: ?nlm)", "adjusted penetrance")

#Quick check to compare results:
  #identical(sums[,3],sums[,4]) #Check if all results are identical when rounded to 3dp
    #[1] TRUE                  
  
  #All results correspond and MLE convergence values are all 1 or 2, indicating correct solution found [See: ?nlm]


#Return results in CSV file
write.csv(sums,"adpenetrance_TS5_MLE_validation.csv")

