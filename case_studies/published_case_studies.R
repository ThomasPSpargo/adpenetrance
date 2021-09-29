#NOTE:
#This script is associated with the following GitHub repository: https://github.com/thomaspspargo/adpenetrance

#It details the calculations performed in the case studies upon which the adpenetrance approach was tested

#See the publication for full description of these examples
  #Calculations can also be performed without the error terms
  #The results will be returned as .csv files in the working directory
  #Separate results tables are generated according to tables presented in the article and supplementary materials

#Load penetrance function and subfunctions
source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/adpenetrance_function.R")

#Import dataset
Cases <- read.csv("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/case_studies/case_data.csv",header=TRUE)


#Take data for case studies 1 and 2, filtering case study 1 to include only sample populations with >5 variant observations in familial and sporadic states
Cases_c1_2 <- Cases[which(Cases$Case_study<3 & Cases$Var_spor>5 & Cases$Var_fam>5),]

#Reorder the dataset
Cases_c1_2 <- rbind(
  Cases_c1_2[-c(2,6,7,10:14),], #Extract the individual estimates from case 1
  Cases_c1_2[c(2,6,7,10),],     #Extract the joint-populations sampled from case 1
  Cases_c1_2[11:14,]            #Extract case 2
)

#Loop across case Studies 1 and 2 to derive penetrance estimates
for(k in 1:length(Cases_c1_2$Population_sampled)){
  
  #FS
  #Every sample has data for the familial and sporadic states, run for every sample
    #ObsX is the familial/(familial+sporadic) rate
  cases_out <- adpenetrance(N=Cases_c1_2[k,"TFR"],
                           MF=Cases_c1_2[k,"Var_fam"]/Cases_c1_2[k,"N_familial"],
                            MF_SE=sqrt(Cases_c1_2[k,"Var_fam"]/Cases_c1_2[k,"N_familial"]*(1-Cases_c1_2[k,"Var_fam"]/Cases_c1_2[k,"N_familial"])/Cases_c1_2[k,"N_familial"]),
                           MS=Cases_c1_2[k,"Var_spor"]/Cases_c1_2[k,"N_sporadic"],
                            MS_SE=sqrt(Cases_c1_2[k,"Var_spor"]/Cases_c1_2[k,"N_sporadic"]*(1-Cases_c1_2[k,"Var_spor"]/Cases_c1_2[k,"N_sporadic"])/Cases_c1_2[k,"N_sporadic"]),
                           #MU=Cases_c1_2[k,"Var_unaf"]/Cases_c1_2[k,"N_unaffected"],
                            #MU_SE=sqrt(Cases_c1_2[k,"Var_unaf"]/Cases_c1_2[k,"N_unaffected"]*(1-Cases_c1_2[k,"Var_unaf"]/Cases_c1_2[k,"N_unaffected"])/Cases_c1_2[k,"N_unaffected"]),
                           PF=Cases_c1_2[k,"PF"]#,
                           #PA=Cases_c1_2[k,"PA"]
                           )$output
  
  #Matrix of output
  est_out_fs <- matrix(c(Cases_c1_2$Population_sampled[k],2,
    cases_out[1,1], cases_out[1,2], cases_out[1,3],
    cases_out[3,1], cases_out[3,2], cases_out[3,3],
    cases_out[4,1], cases_out[4,2], cases_out[4,3]),nrow=1)
    
  
  
  #Run additional analyses for those in which all state combinations are being modelled
  if(Cases_c1_2[k,"Var_unaf"]>1){

    #FSU
    #Run penetrance for every population with fsu states sampled
      #ObsX is the familial/(familial+sporadic+unaffected) rate
    cases_out <- adpenetrance(N=Cases_c1_2[k,"TFR"],
                              MF=Cases_c1_2[k,"Var_fam"]/Cases_c1_2[k,"N_familial"],
                              MF_SE=sqrt(Cases_c1_2[k,"Var_fam"]/Cases_c1_2[k,"N_familial"]*(1-Cases_c1_2[k,"Var_fam"]/Cases_c1_2[k,"N_familial"])/Cases_c1_2[k,"N_familial"]),
                              MS=Cases_c1_2[k,"Var_spor"]/Cases_c1_2[k,"N_sporadic"],
                              MS_SE=sqrt(Cases_c1_2[k,"Var_spor"]/Cases_c1_2[k,"N_sporadic"]*(1-Cases_c1_2[k,"Var_spor"]/Cases_c1_2[k,"N_sporadic"])/Cases_c1_2[k,"N_sporadic"]),
                              MU=Cases_c1_2[k,"Var_unaf"]/Cases_c1_2[k,"N_unaffected"],
                              MU_SE=sqrt(Cases_c1_2[k,"Var_unaf"]/Cases_c1_2[k,"N_unaffected"]*(1-Cases_c1_2[k,"Var_unaf"]/Cases_c1_2[k,"N_unaffected"])/Cases_c1_2[k,"N_unaffected"]),
                              PF=Cases_c1_2[k,"PF"],
                              PA=1/Cases_c1_2[k,"PA"]
                              )$output
    
    #Matrix of output
    est_out_fsu <- matrix(c(Cases_c1_2$Population_sampled[k],1,
                        cases_out[1,1], cases_out[1,2], cases_out[1,3],
                        cases_out[3,1], cases_out[3,2], cases_out[3,3],
                        cases_out[4,1], cases_out[4,2], cases_out[4,3]),nrow=1)
    
    
    

    #FU
    #Run penetrance for every population with fu states sampled
      #ObsX is the familial/(familial+unaffected) rate
    cases_out <- adpenetrance(N=Cases_c1_2[k,"TFR"],
                              MF=Cases_c1_2[k,"Var_fam"]/Cases_c1_2[k,"N_familial"],
                              MF_SE=sqrt(Cases_c1_2[k,"Var_fam"]/Cases_c1_2[k,"N_familial"]*(1-Cases_c1_2[k,"Var_fam"]/Cases_c1_2[k,"N_familial"])/Cases_c1_2[k,"N_familial"]),
                              #MS=Cases_c1_2[k,"Var_spor"]/Cases_c1_2[k,"N_sporadic"],
                              #MS_SE=sqrt(Cases_c1_2[k,"Var_spor"]/Cases_c1_2[k,"N_sporadic"]*(1-Cases_c1_2[k,"Var_spor"]/Cases_c1_2[k,"N_sporadic"])/Cases_c1_2[k,"N_sporadic"]),
                              MU=Cases_c1_2[k,"Var_unaf"]/Cases_c1_2[k,"N_unaffected"],
                              MU_SE=sqrt(Cases_c1_2[k,"Var_unaf"]/Cases_c1_2[k,"N_unaffected"]*(1-Cases_c1_2[k,"Var_unaf"]/Cases_c1_2[k,"N_unaffected"])/Cases_c1_2[k,"N_unaffected"]),
                              PF=Cases_c1_2[k,"PF"],
                              PA=1/Cases_c1_2[k,"PA"]
                              )$output
    
    #Matrix of output
    est_out_fu <- matrix(c(Cases_c1_2$Population_sampled[k],3,
                            cases_out[1,1], cases_out[1,2], cases_out[1,3],
                            cases_out[3,1], cases_out[3,2], cases_out[3,3],
                            cases_out[4,1], cases_out[4,2], cases_out[4,3]),nrow=1)
    
    
    
    #SU
    #Run penetrance for every population with su states sampled
      #ObsX is the sporadic/(sporadic+unaffected) rate
    cases_out <- adpenetrance(N=Cases_c1_2[k,"TFR"],
                              #MF=Cases_c1_2[k,"Var_fam"]/Cases_c1_2[k,"N_familial"],
                              #MF_SE=sqrt(Cases_c1_2[k,"Var_fam"]/Cases_c1_2[k,"N_familial"]*(1-Cases_c1_2[k,"Var_fam"]/Cases_c1_2[k,"N_familial"])/Cases_c1_2[k,"N_familial"]),
                              MS=Cases_c1_2[k,"Var_spor"]/Cases_c1_2[k,"N_sporadic"],
                              MS_SE=sqrt(Cases_c1_2[k,"Var_spor"]/Cases_c1_2[k,"N_sporadic"]*(1-Cases_c1_2[k,"Var_spor"]/Cases_c1_2[k,"N_sporadic"])/Cases_c1_2[k,"N_sporadic"]),
                              MU=Cases_c1_2[k,"Var_unaf"]/Cases_c1_2[k,"N_unaffected"],
                              MU_SE=sqrt(Cases_c1_2[k,"Var_unaf"]/Cases_c1_2[k,"N_unaffected"]*(1-Cases_c1_2[k,"Var_unaf"]/Cases_c1_2[k,"N_unaffected"])/Cases_c1_2[k,"N_unaffected"]),
                              PF=Cases_c1_2[k,"PF"],
                              PA=1/Cases_c1_2[k,"PA"]
                              )$output
    
    #Matrix of output
    est_out_su <- matrix(c(Cases_c1_2$Population_sampled[k],4,
                            cases_out[1,1], cases_out[1,2], cases_out[1,3],
                            cases_out[3,1], cases_out[3,2], cases_out[3,3],
                            cases_out[4,1], cases_out[4,2], cases_out[4,3]),nrow=1)
    
    #Combine all the out objects, where all states have been modelled
    est_out <- rbind(est_out_fsu,est_out_fs,est_out_fu,est_out_su)
    
  } else {
    #If only fs modelled, take just that object
    est_out <- est_out_fs
  }
  
  #Combine into an overall output object
  if(k==1){
    est_out_c1_2 <- est_out
  } else {
    est_out_c1_2 <- rbind(est_out_c1_2,est_out)
  }
  
}
  
  #Analyse supplementary analyses for case study 3, which are based upon variant counts
    
  #A5V
  #ObsX is the familial/(familial+sporadic) rate
  k <- which(Cases$Population_sampled=="SOD1-A5V")
  cases_out <- adpenetrance(N=Cases[k,"TFR"],
                            MF=Cases[k,"Var_fam"]/Cases[k,"N_familial"],
                            MF_SE=sqrt(Cases[k,"Var_fam"]/Cases[k,"N_familial"]*(1-Cases[k,"Var_fam"]/Cases[k,"N_familial"])/Cases[k,"N_familial"]),
                            MS=Cases[k,"Var_spor"]/Cases[k,"N_sporadic"],
                            MS_SE=sqrt(Cases[k,"Var_spor"]/Cases[k,"N_sporadic"]*(1-Cases[k,"Var_spor"]/Cases[k,"N_sporadic"])/Cases[k,"N_sporadic"]),
                            PF=Cases[k,"PF"]
  )$output
  
  #Matrix of output
  est_out_A5V <- matrix(c(Cases$Population_sampled[k],2,
                         cases_out[1,1], cases_out[1,2], cases_out[1,3],
                         cases_out[3,1], cases_out[3,2], cases_out[3,3],
                         cases_out[4,1], cases_out[4,2], cases_out[4,3]),nrow=1)
  
  #D91A
  #ObsX is the sporadic/(sporadic+unaffected) rate
  k <- which(Cases$Population_sampled=="SOD1-D91A")
  cases_out <- adpenetrance(N=Cases[k,"TFR"],
                            MS=Cases[k,"Var_spor"]/Cases[k,"N_sporadic"],
                            MS_SE=sqrt(Cases[k,"Var_spor"]/Cases[k,"N_sporadic"]*(1-Cases[k,"Var_spor"]/Cases[k,"N_sporadic"])/Cases[k,"N_sporadic"]),
                            MU=Cases[k,"Var_unaf"]/Cases[k,"N_unaffected"],
                            MU_SE=sqrt(Cases[k,"Var_unaf"]/Cases[k,"N_unaffected"]*(1-Cases[k,"Var_unaf"]/Cases[k,"N_unaffected"])/Cases[k,"N_unaffected"]),
                            PF=Cases[k,"PF"],
                            PA=1/Cases[k,"PA"]
  )$output
  
  #Matrix of output
  est_out_D91A <- matrix(c(Cases$Population_sampled[k],4,
                         cases_out[1,1], cases_out[1,2], cases_out[1,3],
                         cases_out[3,1], cases_out[3,2], cases_out[3,3],
                         cases_out[4,1], cases_out[4,2], cases_out[4,3]),nrow=1)
  
  
  
  #I114T
  #ObsX is the familial/(familial+sporadic) rate
  k <- which(Cases$Population_sampled=="SOD1-I114T")
  cases_out <- adpenetrance(N=Cases[k,"TFR"],
                            MF=Cases[k,"Var_fam"]/Cases[k,"N_familial"],
                            MF_SE=sqrt(Cases[k,"Var_fam"]/Cases[k,"N_familial"]*(1-Cases[k,"Var_fam"]/Cases[k,"N_familial"])/Cases[k,"N_familial"]),
                            MS=Cases[k,"Var_spor"]/Cases[k,"N_sporadic"],
                            MS_SE=sqrt(Cases[k,"Var_spor"]/Cases[k,"N_sporadic"]*(1-Cases[k,"Var_spor"]/Cases[k,"N_sporadic"])/Cases[k,"N_sporadic"]),
                            PF=Cases[k,"PF"]
  )$output
  
  #Matrix of output
  est_out_I114T <- matrix(c(Cases$Population_sampled[k],2,
                         cases_out[1,1], cases_out[1,2], cases_out[1,3],
                         cases_out[3,1], cases_out[3,2], cases_out[3,3],
                         cases_out[4,1], cases_out[4,2], cases_out[4,3]),nrow=1)
  
  
  #Take main analyses for case study 3 and analyses for case 4, which are based on variant frequencies and lower confidence intervals
  Cases_c3_4   <- Cases[which(Cases$Case_study>=3 & Cases$Freq_fam>0),]
    
  #Prepare output matrix of nrows=4
  est_out_c3_4 <- matrix(as.numeric(NA),ncol=11,nrow=length(Cases_c3_4$Population_sampled))
  
    #Loop across each row of input data
    for(k in 1:length(Cases_c3_4$Population_sampled)){
      
      #ObsX is the familial/(familial+sporadic) rate
      cases_out <- adpenetrance(N=Cases_c3_4[k,"TFR"],
               MF=Cases_c3_4[k,"Freq_fam"], MF_SE=(Cases_c3_4[k,"Freq_fam"]-Cases_c3_4[k,"Freq_fam_LCI"])/1.96,
               MS=Cases_c3_4[k,"Freq_spor"], MS_SE =(Cases_c3_4[k,"Freq_spor"]-Cases_c3_4[k,"Freq_spor_LCI"])/1.96,
               PF=Cases_c3_4[k,"PF"], Zout=1.96)$output
    
      #Matrix of output
      est_out_c3_4[k,] <- c(Cases_c3_4$Population_sampled[k],2,
                                cases_out[1,1], cases_out[1,2], cases_out[1,3],
                                cases_out[3,1], cases_out[3,2], cases_out[3,3],
                                cases_out[4,1], cases_out[4,2], cases_out[4,3])
      
    }
    
    
#Convert analyses into format for write-up
  
#Prepare Table 2 results
  est_out_T2<- rbind(est_out_c1_2[19:26,], est_out_c3_4)
  
  #Assign column names
  colnames(est_out_T2) <- c("Pop_sampled","states_modelled","RX_lower","RX_est","RX_upper",
                             "Un_pen_lower","Un_pen_est","Un_pen_upper",
                             "Adj_pen_lower","Adj_pen_est","Adj_pen_upper")
  
  #Convert to data.frame and truncate the values to 3dp
  sum_Table <- round(as.data.frame(est_out_T2),digits=3)
  
  #Assign factor levels to pop sampled and states modelled
  sum_Table$Pop_sampled <- as.factor(sum_Table$Pop_sampled)
    levels(sum_Table$Pop_sampled) <- levels(Cases$Population_sampled)[as.numeric(levels(sum_Table$Pop_sampled))]
  sum_Table$states_modelled <- as.factor(sum_Table$states_modelled)
    levels(sum_Table$states_modelled)  <- c("fsu","fs","fu","su")
  
  
  #Summarise results in "estimate (LowerCI, UpperCI)" format
  sums <- cbind(paste0(sum_Table$Pop_sampled," (",sum_Table$states_modelled,")"),
                paste0(sum_Table$RX_est," (",sum_Table$RX_lower,", ",sum_Table$RX_upper,")"),
                paste0(sum_Table$Un_pen_est," (",sum_Table$Un_pen_lower,", ",sum_Table$Un_pen_upper,")"),
                paste0(sum_Table$Adj_pen_est," (",sum_Table$Adj_pen_lower,", ",sum_Table$Adj_pen_upper,")")
  )
  
  #Assign colnames
  colnames(sums) <- c("Analysis", "RXobs","unadjusted penetrance", "adjusted penetrance")
  
  #Store in Table 2 results object
  T2 <- sums
  
  
  

  
  

#Prepare Table S2 results
  est_out_TS2 <- est_out_c1_2[1:22,]

  #Assign column names
  colnames(est_out_TS2) <- c("Pop_sampled","states_modelled","RX_lower","RX_est","RX_upper",
                            "Un_pen_lower","Un_pen_est","Un_pen_upper",
                            "Adj_pen_lower","Adj_pen_est","Adj_pen_upper")
  
  #Convert to data.frame and truncate the values to 3dp
  sum_Table <- round(as.data.frame(est_out_TS2),digits=3)
  
  #Assign factor levels to pop sampled and states modelled
  sum_Table$Pop_sampled <- as.factor(sum_Table$Pop_sampled)
    levels(sum_Table$Pop_sampled) <- levels(Cases$Population_sampled)[as.numeric(levels(sum_Table$Pop_sampled))]
  sum_Table$states_modelled <- as.factor(sum_Table$states_modelled)
    levels(sum_Table$states_modelled)  <- c("fsu","fs","fu","su")
  
  
  #Summarise results in "estimate (LowerCI, UpperCI)" format
  sums <- cbind(paste0(sum_Table$Pop_sampled," (",sum_Table$states_modelled,")"),
                paste0(sum_Table$RX_est," (",sum_Table$RX_lower,", ",sum_Table$RX_upper,")"),
                paste0(sum_Table$Un_pen_est," (",sum_Table$Un_pen_lower,", ",sum_Table$Un_pen_upper,")"),
                paste0(sum_Table$Adj_pen_est," (",sum_Table$Adj_pen_lower,", ",sum_Table$Adj_pen_upper,")")
  )
  #Assign column names
  colnames(sums) <- c("Analysis", "RXobs","unadjusted penetrance", "adjusted penetrance")
  
  #Store in Table S2 results object
  TS2 <- sums
  
  
#Prepare Table S3 results
  est_out_TS3 <- rbind(est_out_A5V,est_out_D91A,est_out_I114T)
  
  #Assign column names
  colnames(est_out_TS3) <- c("Pop_sampled","states_modelled","RX_lower","RX_est","RX_upper",
                             "Un_pen_lower","Un_pen_est","Un_pen_upper",
                             "Adj_pen_lower","Adj_pen_est","Adj_pen_upper")
  
  #Convert to data.frame and truncate the values.
    #Do not round - manual rounding, because of small D91A estimate
  sum_Table <- as.data.frame(est_out_TS3)
  
  #Assign factor levels to pop sampled and states modelled
  sum_Table$Pop_sampled <- as.factor(sum_Table$Pop_sampled)
  levels(sum_Table$Pop_sampled) <- levels(Cases$Population_sampled)[as.numeric(levels(sum_Table$Pop_sampled))]
  sum_Table$states_modelled <- as.factor(sum_Table$states_modelled)
  levels(sum_Table$states_modelled)  <- c("fsu","fs","fu","su")
  
  
  #Summarise results in "estimate (LowerCI, UpperCI)" format
  sums <- cbind(paste0(sum_Table$Pop_sampled," (",sum_Table$states_modelled,")"),
                paste0(sum_Table$RX_est," (",sum_Table$RX_lower,", ",sum_Table$RX_upper,")"),
                paste0(sum_Table$Un_pen_est," (",sum_Table$Un_pen_lower,", ",sum_Table$Un_pen_upper,")"),
                paste0(sum_Table$Adj_pen_est," (",sum_Table$Adj_pen_lower,", ",sum_Table$Adj_pen_upper,")")
  )
  #Assign column names
  colnames(sums) <- c("Analysis", "RXobs","unadjusted penetrance", "adjusted penetrance")
  
  #Store in Table S3 output object
  TS3 <- sums
  
  

  #Write results into CSV files:
  write.csv(T2,"adpenetrance_T2.csv")
  write.csv(TS2,"adpenetrance_TS2.csv")
  write.csv(TS3,"adpenetrance_TS3.csv")
  
   