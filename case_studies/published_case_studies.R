#NOTE:
#This script is associated with the following GitHub repository: https://github.com/thomaspspargo/adpenetrance

#It details the calculations performed in the case studies upon which the adpenetrance approach was tested

#See the publication for full description of these examples
  #Calculations can also be performed without the error terms
  #The results will be returned as .csv files in the working directory
  #Separate results tables are generated according to tables presented in the article and supplementary materials

#Load penetrance function, subfunctions, and the extra function for calculating residual risk g
source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/adpenetrance_function.R")
source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/getResidualRisk.R")

### Define function for calculating standard error for proportion based on variant counts
  #count is the number of variants in the sample
  #size is the sample size
se_p <- function(count,size){
  p <- count/size
  n <- size
  
  sqrt(p*(1-p)/n)
}

#Import dataset
Cases <- read.csv("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/case_studies/case_data.csv",header=TRUE)


#Take data for case studies 1, 2, and some of 3
  #filter case study 1 to include only sample populations with >5 variant observations in familial and sporadic states
  #filter case study 3 to include only the variant-specific results
Cases_wcount <- Cases[c(which(Cases$Case_study<3 & Cases$Var_spor>5 & Cases$Var_fam>5),
                        which(Cases$Population_sampled %in% c("SOD1-A5V","SOD1-D91A","SOD1-I114T")))
                      ,]


#Identify the estimates aggregated across European regions and then total across the world.
aggregate<-c(
  grep("Europe",Cases_wcount$Population_sampled),
  grep("Total",Cases_wcount$Population_sampled)
)

#Sort into ascending order
aggregate<- sort(aggregate)

#Reorder the dataset
Cases_wcount <- rbind(
  Cases_wcount[-aggregate,][Cases_wcount[-aggregate,]$Case_study==1,], #Extract the individual estimates from case 1
  Cases_wcount[aggregate,],                                            #Extract the joint-populations sampled from case 1 (i.e. those aggregated)
  Cases_wcount[Cases_wcount$Case_study==2,],                           #Extract case 2
  Cases_wcount[Cases_wcount$Case_study==3,]                            #Extract case 3
)

#Loop across case Studies 1 and 2 to derive penetrance estimates
for(k in 1:length(Cases_wcount$Population_sampled)){
  
  #If the variant frequency is attained by 0/0 division, NaN produced.
  #This represents having no sampled familial state, given that in this state the variant frequency of 0
  #Accordingly, substitute NaN with 0
  MF_get <- Cases_wcount[k,"Var_fam"]/Cases_wcount[k,"N_familial"]
  if(is.nan(MF_get)){MF_get <- 0}
  
  MS_get <- Cases_wcount[k,"Var_spor"]/Cases_wcount[k,"N_sporadic"]
  if(is.nan(MS_get)){MS_get <- 0}
  
  #Compute the remaining disease for a population member who does not harbour the variant
  residualRiskG<- getResidualRisk(MF=MF_get,
                                  MS=MS_get,
                                  PF=Cases_wcount[k,"PF"],
                                  PA=1/Cases_wcount[k,"PA"])
  
  
  residualRiskG<- c(0, residualRiskG) #Convert to 2-element vector, for running function in loop across both base and extended model, useG=NA runs non-extended model
  
  if(Cases_wcount[k,"Var_fam"]>0 & Cases_wcount[k,"Var_spor"]>0){
    
    cases_out <- list(basemodel=NULL,extendedmodel=NULL)
    for(j in 1:length(residualRiskG)){
      #FS
      #Every sample has data for the familial and sporadic states, run for every sample
      #ObsX is the familial/(familial+sporadic) rate
      cases_out[[j]] <- adpenetrance(N=Cases_wcount[k,"TFR"],
                                MF=Cases_wcount[k,"Var_fam"]/Cases_wcount[k,"N_familial"],
                                MF_SE = se_p(Cases_wcount[k,"Var_fam"],Cases_wcount[k,"N_familial"]),
                                
                                MS=Cases_wcount[k,"Var_spor"]/Cases_wcount[k,"N_sporadic"],
                                MS_SE = se_p(Cases_wcount[k,"Var_spor"],Cases_wcount[k,"N_sporadic"]),
                                
                                # MU=Cases_wcount[k,"Var_unaf"]/Cases_wcount[k,"N_unaffected"],
                                # MU_SE = se_p(Cases_wcount[k,"Var_unaf"],Cases_wcount[k,"N_unaffected"]),
                                PF=Cases_wcount[k,"PF"],
                                #PA=1/Cases_wcount[k,"PA"]
                                useG=residualRiskG[j]
      )
    }
    
    #Vector of output
    est_out_fs <- c(2,
                    cases_out$basemodel$output[1,1], cases_out$basemodel$output[1,2], cases_out$basemodel$output[1,3],
                    cases_out$basemodel$output[3,1], cases_out$basemodel$output[3,2], cases_out$basemodel$output[3,3],
                    cases_out$basemodel$output[4,1], cases_out$basemodel$output[4,2], cases_out$basemodel$output[4,3],
                    cases_out$extendedmodel$output[3,1], cases_out$extendedmodel$output[3,2], cases_out$extendedmodel$output[3,3],
                    cases_out$extendedmodel$output[4,1], cases_out$extendedmodel$output[4,2], cases_out$extendedmodel$output[4,3],
                    cases_out$extendedmodel$ResidualRiskG
                    )
    
  }
    
  #Run additional checks and analyses for those for which unaffected state data are provided
  if(Cases_wcount[k,"Var_unaf"]>1){
    
    if(Cases_wcount[k,"Var_fam"]>0 & Cases_wcount[k,"Var_spor"]>0 & Cases_wcount[k,"Var_unaf"]>0 ){
      
      cases_out <- list(basemodel=NULL,extendedmodel=NULL)
      for(j in 1:length(residualRiskG)){
      #FSU
      #Run penetrance for every population with fsu states sampled
      #ObsX is the familial/(familial+sporadic+unaffected) rate
      cases_out[[j]] <- adpenetrance(N=Cases_wcount[k,"TFR"],
                                MF=Cases_wcount[k,"Var_fam"]/Cases_wcount[k,"N_familial"],
                                MF_SE = se_p(Cases_wcount[k,"Var_fam"],Cases_wcount[k,"N_familial"]),
                                
                                MS=Cases_wcount[k,"Var_spor"]/Cases_wcount[k,"N_sporadic"],
                                MS_SE = se_p(Cases_wcount[k,"Var_spor"],Cases_wcount[k,"N_sporadic"]),
                                
                                MU=Cases_wcount[k,"Var_unaf"]/Cases_wcount[k,"N_unaffected"],
                                MU_SE = se_p(Cases_wcount[k,"Var_unaf"],Cases_wcount[k,"N_unaffected"]),
                                PF=Cases_wcount[k,"PF"],
                                PA=1/Cases_wcount[k,"PA"],
                                useG=residualRiskG[j]
      )
      }
      
      #Vector of output
      est_out_fsu <- c(1,
                       cases_out$basemodel$output[1,1], cases_out$basemodel$output[1,2], cases_out$basemodel$output[1,3],
                       cases_out$basemodel$output[3,1], cases_out$basemodel$output[3,2], cases_out$basemodel$output[3,3],
                       cases_out$basemodel$output[4,1], cases_out$basemodel$output[4,2], cases_out$basemodel$output[4,3],
                       cases_out$extendedmodel$output[3,1], cases_out$extendedmodel$output[3,2], cases_out$extendedmodel$output[3,3],
                       cases_out$extendedmodel$output[4,1], cases_out$extendedmodel$output[4,2], cases_out$extendedmodel$output[4,3],
                       cases_out$extendedmodel$ResidualRiskG
      )
      
    }
    

    if(Cases_wcount[k,"Var_fam"]>0 & Cases_wcount[k,"Var_unaf"]>0){
      
      cases_out <- list(basemodel=NULL,extendedmodel=NULL)
      for(j in 1:length(residualRiskG)){
    #FU
    #Run penetrance for every population with fu states sampled
      #ObsX is the familial/(familial+unaffected) rate
    cases_out[[j]] <- adpenetrance(N=Cases_wcount[k,"TFR"],
                              MF=Cases_wcount[k,"Var_fam"]/Cases_wcount[k,"N_familial"],
                              MF_SE = se_p(Cases_wcount[k,"Var_fam"],Cases_wcount[k,"N_familial"]),
                              
                              # MS=Cases_wcount[k,"Var_spor"]/Cases_wcount[k,"N_sporadic"],
                              # MS_SE = se_p(Cases_wcount[k,"Var_spor"],Cases_wcount[k,"N_sporadic"]),
                              
                              MU=Cases_wcount[k,"Var_unaf"]/Cases_wcount[k,"N_unaffected"],
                              MU_SE = se_p(Cases_wcount[k,"Var_unaf"],Cases_wcount[k,"N_unaffected"]),
                              PF=Cases_wcount[k,"PF"],
                              PA=1/Cases_wcount[k,"PA"],
                              useG=residualRiskG[j]
    )
      }
    
    #Vector of output
    est_out_fu <- c(3,
                      cases_out$basemodel$output[1,1], cases_out$basemodel$output[1,2], cases_out$basemodel$output[1,3],
                      cases_out$basemodel$output[3,1], cases_out$basemodel$output[3,2], cases_out$basemodel$output[3,3],
                      cases_out$basemodel$output[4,1], cases_out$basemodel$output[4,2], cases_out$basemodel$output[4,3],
                      cases_out$extendedmodel$output[3,1], cases_out$extendedmodel$output[3,2], cases_out$extendedmodel$output[3,3],
                      cases_out$extendedmodel$output[4,1], cases_out$extendedmodel$output[4,2], cases_out$extendedmodel$output[4,3],
                    cases_out$extendedmodel$ResidualRiskG
                    )
    
    }
    
    if(Cases_wcount[k,"Var_spor"]>0 & Cases_wcount[k,"Var_unaf"]>0 ){
      
      cases_out <- list(basemodel=NULL,extendedmodel=NULL)
      for(j in 1:length(residualRiskG)){
      #SU
      #Run penetrance for every population with su states sampled
      #ObsX is the sporadic/(sporadic+unaffected) rate
      cases_out[[j]] <- adpenetrance(N=Cases_wcount[k,"TFR"],
                                # MF=Cases_wcount[k,"Var_fam"]/Cases_wcount[k,"N_familial"],
                                # MF_SE = se_p(Cases_wcount[k,"Var_fam"],Cases_wcount[k,"N_familial"]),
                                
                                MS=Cases_wcount[k,"Var_spor"]/Cases_wcount[k,"N_sporadic"],
                                MS_SE = se_p(Cases_wcount[k,"Var_spor"],Cases_wcount[k,"N_sporadic"]),
                                
                                MU=Cases_wcount[k,"Var_unaf"]/Cases_wcount[k,"N_unaffected"],
                                MU_SE = se_p(Cases_wcount[k,"Var_unaf"],Cases_wcount[k,"N_unaffected"]),
                                PF=Cases_wcount[k,"PF"],
                                PA=1/Cases_wcount[k,"PA"],
                                useG=residualRiskG[j]
      )
      }
      
      #Vector of output
      est_out_su <- c(4,
                      cases_out$basemodel$output[1,1], cases_out$basemodel$output[1,2], cases_out$basemodel$output[1,3],
                      cases_out$basemodel$output[3,1], cases_out$basemodel$output[3,2], cases_out$basemodel$output[3,3],
                      cases_out$basemodel$output[4,1], cases_out$basemodel$output[4,2], cases_out$basemodel$output[4,3],
                      cases_out$extendedmodel$output[3,1], cases_out$extendedmodel$output[3,2], cases_out$extendedmodel$output[3,3],
                      cases_out$extendedmodel$output[4,1], cases_out$extendedmodel$output[4,2], cases_out$extendedmodel$output[4,3],
                      cases_out$extendedmodel$ResidualRiskG
      )
      
    }
    
  }
    
    
  #Combine each output into summary matrix according to whether they're present in global environment
  est_out <- rbind(
    if("est_out_fsu" %in% ls()){est_out_fsu},
    if("est_out_fs" %in% ls()){est_out_fs},
    if("est_out_fu" %in% ls()){est_out_fu},
    if("est_out_su" %in% ls()){est_out_su}
  )
  
  #Append the case study number and the name of the population sampled; convert to data frame
  est_out <- data.frame(study = Cases_wcount$Case_study[k],Pop_Samp = Cases_wcount$Population_sampled[k],est_out)
  
  #Combine into an overall output object
  if(k==1){ #For the first loop, make the df
    est_out_wcount <- est_out
  } else {
    est_out_wcount <- rbind(est_out_wcount,est_out)
  }
  
  #Any present output objects from global environment
  suppressWarnings(rm(list=c("est_out_fsu","est_out_fs","est_out_fu","est_out_su")))
  
  
}
  

#Take main analyses for case study 3 and analyses for case 4, which are based on variant frequencies and lower confidence intervals
Cases_wfreq   <- Cases[which(Cases$Case_study>=3 & Cases$Freq_fam>0),]


est_out_wfreq <- matrix(NA_real_,ncol=17,nrow=length(Cases_wfreq$Population_sampled)) #Prepare output matrix of nrows=4
for(k in 1:length(Cases_wfreq$Population_sampled)){ #Loop across each row of input data

  #Compute the remaining disease for a population member who does not harbour the variant
  residualRiskG<- getResidualRisk(MF=Cases_wfreq[k,"Freq_fam"],
                                  MS=Cases_wfreq[k,"Freq_spor"],
                                  PF=Cases_wfreq[k,"PF"],
                                  PA=1/Cases_wfreq[k,"PA"])
  
  residualRiskG<- c(0, residualRiskG) #Convert to 2-element vector, for running function in loop across both base and extended model, useG=NA runs non-extended model
  
  #Run model, looping across base and extended model
  cases_out <- list(basemodel=NULL,extendedmodel=NULL)
  for(j in 1:length(residualRiskG)){
    #ObsX is the familial/(familial+sporadic) rate
    cases_out[[j]] <- adpenetrance(N=Cases_wfreq[k,"TFR"],
                                   MF=Cases_wfreq[k,"Freq_fam"], MF_SE=(Cases_wfreq[k,"Freq_fam"]-Cases_wfreq[k,"Freq_fam_LCI"])/1.96,
                                   MS=Cases_wfreq[k,"Freq_spor"], MS_SE =(Cases_wfreq[k,"Freq_spor"]-Cases_wfreq[k,"Freq_spor_LCI"])/1.96,
                                   PF=Cases_wfreq[k,"PF"],Zout=1.96,
                                   useG=residualRiskG[j]
    )
  }
  #Add vector to output matrix
  est_out_wfreq[k,] <- c(2,
                         cases_out$basemodel$output[1,1], cases_out$basemodel$output[1,2], cases_out$basemodel$output[1,3],
                         cases_out$basemodel$output[3,1], cases_out$basemodel$output[3,2], cases_out$basemodel$output[3,3],
                         cases_out$basemodel$output[4,1], cases_out$basemodel$output[4,2], cases_out$basemodel$output[4,3],
                         cases_out$extendedmodel$output[3,1], cases_out$extendedmodel$output[3,2], cases_out$extendedmodel$output[3,3],
                         cases_out$extendedmodel$output[4,1], cases_out$extendedmodel$output[4,2], cases_out$extendedmodel$output[4,3],
                         cases_out$extendedmodel$ResidualRiskG
  )
}

#Append the case study number and name of the population sampled; convert to data frame
est_out_wfreq <- data.frame(study = Cases_wfreq$Case_study,Pop_Samp = Cases_wfreq$Population_sampled,est_out_wfreq)



#Convert analyses into format for write-up using custom function
  
  #Define custom function -x is the only argument, the orignal format results table
  reformatTab<- function(x){
    #Assign column names, to be called when summarising
    colnames(x) <- c("study", "Pop_sampled","states_modelled","RX_lower","RX_est","RX_upper",
                              "BaseUn_pen_lower","BaseUn_pen_est","BaseUn_pen_upper",
                              "BaseAdj_pen_lower","BaseAdj_pen_est","BaseAdj_pen_upper",
                     "ExtUn_pen_lower","ExtUn_pen_est","ExtUn_pen_upper",
                     "ExtAdj_pen_lower","ExtAdj_pen_est","ExtAdj_pen_upper","ExtResidRiskG")
    
    x[3:ncol(x)] <- signif(x[3:ncol(x)],digits=3) #Round numeric columns
    
    x$states_modelled <- factor(x$states_modelled,  #Convert states modelled into a factor variable
                                levels = c(1,2,3,4), 
                                labels = c("fsu","fs","fu","su")
                                ) 
    
    
    #Summarise results in "estimate (LowerCI, UpperCI)" format
    sums <- cbind(paste0(x$Pop_sampled," (",x$states_modelled,")"),
                  paste0(x$ExtResidRiskG),
                  paste0(x$RX_est," (",x$RX_lower,", ",x$RX_upper,")"),
                  paste0(x$BaseUn_pen_est," (",x$BaseUn_pen_lower,", ",x$BaseUn_pen_upper,")"),
                  paste0(x$BaseAdj_pen_est," (",x$BaseAdj_pen_lower,", ",x$BaseAdj_pen_upper,")"),
                  paste0(x$ExtUn_pen_est," (",x$ExtUn_pen_lower,", ",x$ExtUn_pen_upper,")"),
                  paste0(x$ExtAdj_pen_est," (",x$ExtAdj_pen_lower,", ",x$ExtAdj_pen_upper,")")
    )
    
    
    colnames(sums) <- c("Analysis","ResidRiskG", "RXobs","unadjusted penetrance Gzero", "adjusted penetrance Gzero","unadjusted penetrance useG", "adjusted penetrance useG") #Assign colnames
    
    return(sums)
  }
  
  
#Prepare Table 2 results
    #Rbind: the European gnomad estimate, all BMPR2 estimates, main SOD1 and C9orf72 analyses
  est_out_T2<- rbind(
    est_out_wcount[grep("Europe total_gnomad",est_out_wcount$Pop_Samp),],
    est_out_wcount[est_out_wcount$study==2,],
    est_out_wfreq)
  
  #Reformat with function and store in Table 2 results object
  T2 <- reformatTab(est_out_T2)

  
  
#Prepare Table S2 results
  est_out_TS2 <- est_out_wcount[est_out_wcount$study==1,]
  
  #Reformat with function and store in Table S2 results object
  TS2 <- reformatTab(est_out_TS2)

  
#Prepare Table S3 results
  est_out_TS3 <- est_out_wcount[est_out_wcount$study==3,]
                          
  #Store in Table S3 output object
  TS3 <- reformatTab(est_out_TS3)
  
  #Write results into CSV files:
  write.csv(T2,"adpenetrance_T2.csv")
  write.csv(TS2,"adpenetrance_TS2.csv")
  write.csv(TS3,"adpenetrance_TS3.csv")
  
   