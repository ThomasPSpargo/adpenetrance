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
  
  #Vectorise output
  est_out_fs <- c(2,
                  cases_out[1,1], cases_out[1,2], cases_out[1,3],
                  cases_out[3,1], cases_out[3,2], cases_out[3,3],
                  cases_out[4,1], cases_out[4,2], cases_out[4,3])
    
  
  
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
    
    #Vectorise output
    est_out_fsu <- c(1,
                     cases_out[1,1], cases_out[1,2], cases_out[1,3],
                     cases_out[3,1], cases_out[3,2], cases_out[3,3],
                     cases_out[4,1], cases_out[4,2], cases_out[4,3])
    
    
    

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
    
    #Vectorise output
    est_out_fu <- c(3,
                    cases_out[1,1], cases_out[1,2], cases_out[1,3],
                    cases_out[3,1], cases_out[3,2], cases_out[3,3],
                    cases_out[4,1], cases_out[4,2], cases_out[4,3])
    
    
    
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
    
    #Vectorise output
    est_out_su <- c(4,
                    cases_out[1,1], cases_out[1,2], cases_out[1,3],
                    cases_out[3,1], cases_out[3,2], cases_out[3,3],
                    cases_out[4,1], cases_out[4,2], cases_out[4,3])
    
    #Combine all the out objects, where all states have been modelled
    est_out <- rbind(est_out_fsu,est_out_fs,est_out_fu,est_out_su)
    
  } else {
    #If only fs modelled, take just that object
    est_out <- matrix(est_out_fs,nrow=1)
  }
  
  #Append the name of the population sampled and convert to data frame
  est_out <- data.frame(Pop_Samp = Cases_c1_2$Population_sampled[k],est_out)
  
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
  
  #Vector of output
  est_out_A5V <- c(2,
                   cases_out[1,1], cases_out[1,2], cases_out[1,3],
                   cases_out[3,1], cases_out[3,2], cases_out[3,3],
                   cases_out[4,1], cases_out[4,2], cases_out[4,3])
  
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
  
  #Vector of output
  est_out_D91A <- c(4,
                    cases_out[1,1], cases_out[1,2], cases_out[1,3],
                    cases_out[3,1], cases_out[3,2], cases_out[3,3],
                    cases_out[4,1], cases_out[4,2], cases_out[4,3])
  
  
  
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
  
  #Vector of output
  est_out_I114T <- c(2,
                     cases_out[1,1], cases_out[1,2], cases_out[1,3],
                     cases_out[3,1], cases_out[3,2], cases_out[3,3],
                     cases_out[4,1], cases_out[4,2], cases_out[4,3])
  
  
  #Take main analyses for case study 3 and analyses for case 4, which are based on variant frequencies and lower confidence intervals
  Cases_c3_4   <- Cases[which(Cases$Case_study>=3 & Cases$Freq_fam>0),]
    
  #Prepare output matrix of nrows=4
  est_out_c3_4 <- matrix(NA_real_,ncol=10,nrow=length(Cases_c3_4$Population_sampled))
  
    #Loop across each row of input data
    for(k in 1:length(Cases_c3_4$Population_sampled)){
      
      #ObsX is the familial/(familial+sporadic) rate
      cases_out <- adpenetrance(N=Cases_c3_4[k,"TFR"],
               MF=Cases_c3_4[k,"Freq_fam"], MF_SE=(Cases_c3_4[k,"Freq_fam"]-Cases_c3_4[k,"Freq_fam_LCI"])/1.96,
               MS=Cases_c3_4[k,"Freq_spor"], MS_SE =(Cases_c3_4[k,"Freq_spor"]-Cases_c3_4[k,"Freq_spor_LCI"])/1.96,
               PF=Cases_c3_4[k,"PF"], Zout=1.96)$output
    
      #Add vector to output matrix
      est_out_c3_4[k,] <- c(2,
                            cases_out[1,1], cases_out[1,2], cases_out[1,3],
                            cases_out[3,1], cases_out[3,2], cases_out[3,3],
                            cases_out[4,1], cases_out[4,2], cases_out[4,3])
      
    }
  
  #Append the name of the population sampled and convert to data frame
  est_out_c3_4 <- data.frame(Pop_Samp = Cases_c3_4$Population_sampled,est_out_c3_4)
  
    
    
#Convert analyses into format for write-up using custom function
  
  #Define custom function -x is the only argument, the orignal format results table
  reformatTab<- function(x){
    #Assign column names, to be called when summarising
    colnames(x) <- c("Pop_sampled","states_modelled","RX_lower","RX_est","RX_upper",
                              "Un_pen_lower","Un_pen_est","Un_pen_upper",
                              "Adj_pen_lower","Adj_pen_est","Adj_pen_upper")
    
    x[3:ncol(x)] <- round(x[3:ncol(x)],digits=3) #Round numeric columns
    
    x$states_modelled <- factor(x$states_modelled,  #Convert states modelled into a factor variable
                                levels = c(1,2,3,4), 
                                labels = c("fsu","fs","fu","su")
                                ) 
    
    
    #Summarise results in "estimate (LowerCI, UpperCI)" format
    sums <- cbind(paste0(x$Pop_sampled," (",x$states_modelled,")"),
                  paste0(x$RX_est," (",x$RX_lower,", ",x$RX_upper,")"),
                  paste0(x$Un_pen_est," (",x$Un_pen_lower,", ",x$Un_pen_upper,")"),
                  paste0(x$Adj_pen_est," (",x$Adj_pen_lower,", ",x$Adj_pen_upper,")")
    )
    
    
    colnames(sums) <- c("Analysis", "RXobs","unadjusted penetrance", "adjusted penetrance") #Assign colnames
    
    return(sums)
  }
  
  
#Prepare Table 2 results
  est_out_T2<- rbind(est_out_c1_2[19:26,], est_out_c3_4)
  
  #Reformat with function and store in Table 2 results object
  T2 <- reformatTab(est_out_T2)

  
  
#Prepare Table S2 results
  est_out_TS2 <- est_out_c1_2[1:22,]
  
  #Reformat with function and store in Table 2 results object
  TS2 <- reformatTab(est_out_TS2)

  
#Prepare Table S3 results
  est_out_TS3 <- data.frame(Pop_Samp = c("SOD1-A5V","SOD1-D91A","SOD1-I114T")
                            ,rbind(est_out_A5V,est_out_D91A,est_out_I114T))
    
    
  
  
  #Store in Table S3 output object
  TS3 <- reformatTab(est_out_TS3)
  
  

  #Write results into CSV files:
  write.csv(T2,"adpenetrance_T2.csv")
  write.csv(TS2,"adpenetrance_TS2.csv")
  write.csv(TS3,"adpenetrance_TS3.csv")
  
   