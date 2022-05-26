#This script is associated with the following GitHub repository: https://github.com/ThomasPSpargo/adpenetrance
#It is used for validation of the adpenetrance approach to penetrance estimation using in simulated datasets
  #This simulation generates two example populations (defined based on real sample populations)
    #and estimates the accuracy of adpenetrance estimation in these populations across a series of ground truth penetrance values
    #when input parameters are correctly defined for the sample population
  #A more detailed description is available in the associated documentation

  #This script prints values to indicate progress along loops of the l object - 100 loops are to be performed

  #The simulation data will be written into a .csv file and a ggplot will be constructed, written into a .png file

  #Regarding populations simulated:
    #Structure 1 is representative of sibship distribution in the UK population 1974 birth cohort, see: https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/conceptionandfertilityrates/datasets/childbearingforwomenbornindifferentyearsreferencetable
    #Structure 2 represents the sample of English families represented within the 'Next Steps' cohort study, see: https://doi.org/10.1017/ehs.2020.54

#load necessary packages
  library(plyr)
  library(ggplot2)

#Load adpenetrance function and adpenetrance.errorfit and adpenetrance.unadjusted subfunctions:
  source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/adpenetrance_function.R")
  
#Define simulation populations:
    
    #Define UK population weights:
      #This structure is to be supplied to define_sibstructure in some simulations, when Step-4 error correction is tailored to a known population distribution
      #Column 1 stores sibship levels, Column 2 stores probabilities of each level
    sibstructureUK<- cbind(0:4,c(0.18,0.18,0.37,0.16,0.11))
    
    #Determine Next Steps population weights (see final two rows of Table 2 in publication):
    NextStepN<- 3647+10588                                                    #Total sample size
    Sibs <- c(393+662, 1231+3926,952+3059,556+1591,248+648,119+344,148+358)   #Sib levels 1:7, vector contains number in each sib level
                                                                                #Note that 'only child' represents sibship size 1
    Psibs <- Sibs/NextStepN                                                   #Psibs stores the probabilities of each sibship size
    
    #Build matrix
      #This structure is to be supplied to define_sibstructure in some simulations, when Step-4 error correction is tailored to a known population distribution
      #Column 1 stores sibship levels, Column 2 stores probabilities of each level
    sibstructureNS<- cbind(1:7, Psibs)
    
#Ground truth penetrance values (f):
  f <- seq(0.05,1,length=20) #Set the penetrance values to test
  f <- sort(rep(f,5)) #Amplify to test each value 5 times, sort into ascending order

  sims <- length(f)   #Extract number of penetrance loops to run
  
  #Define required chr vectors:
  all_states <- c("Familial", "Sporadic", "Unaffected") #Disease state levels to be assigned to families
  states=c("fsu","fs","fu","su","au") #Define the disease states to model
  
  #Prepare empty matrices to store output summary at each states modelled
  statetemp <- statetempTailored <- matrix(as.numeric(NA),ncol=10,nrow=5)
  result <- matrix(as.numeric(NA),ncol=10)
  
  colnames(result) <- colnames(statetemp) <- colnames(statetempTailored) <- c("Population","Tailoring","Seed","mean_N","RX","modif","States.Modelled","TruePenetrance","UncorrectedPen", "CorrectedPen")
  
  
  #Factor levels used throughout:
  #Population:
    #1=UK
    #2=NS
  
  #Tailoring: 
    #1=Poisson
    #2=Tailored
  
  #States:
    #1=fsu
    #2=fs
    #3=fu 
    #4=su
    #5=au
  
  sim_seed <- sample(1:30000,sims) #Generate random seeds for use in each loop
  
  #Loop across penetrance values
  for(l in 1:sims){
  #Define seed for the randomisation
  #Generate population of 90000 families resembling UK structure
  set.seed(sim_seed[l])
  sibshipsUK <-sample(sibstructureUK[,1], 90000,replace=TRUE,prob=sibstructureUK[,2])
  
  #Generate population of 90000 families resembling NS structure
  set.seed(sim_seed[l])
  sibshipsNS <-sample(sibstructureNS[,1], 90000,replace=TRUE,prob=sibstructureNS[,2])

  #Store populations in 2 column matrix
  sibships <- matrix(c(sibshipsUK,sibshipsNS),
                     nrow=90000)
  
  #Loop across simulated population
  for(p in 1:ncol(sibships)){
    #If modelling for the UK sibship sample, define sibstructure by UK distribution and define population
    if(p==1){sibstructure <- sibstructureUK
            pop <- 1}
    
    #If modelling for the NS sibship sample define sibstructure by NS distribution and define population
    if(p==2){sibstructure <- sibstructureNS
            pop <- 2}
    
    #Identify mean N value in sample
    obsN <- sum(sibships[,p])/length(sibships[,p])
    
    #Return vector of unique sib_sizes represented in sample
    N <- sort(unique(sibships[,p]))
  
  
  ###Calculate the probability of each disease state expected each value of N at penetrance f
  Punasc = (1-f[l])*(1-f[l]/2)^N                                                  #Probability unaffected 
  Pspor  = f[l]*(1-f[l]/2)^N+N*(f[l]/2)*(1-f[l]/2)^(N-1)*(1-f[l])                 #Probability sporadic/simplex
  Pfam   = 1-((1-f[l]/2)^N+N*(f[l]/2)*((1-f[l]/2)^(N-1))*(1-f[l]))                #Probability familial (= 1 - Punasc - Pspor)

  #Combine the calculated probabilities into a matrix, with row names defined by the value of N
  Probabilities <- matrix(c(Pfam, Pspor, Punasc),
                          nrow=length(N))
  rownames(Probabilities) <- N
  
  #Prepare empty matrix for storing the sibship size and assigned disease state
  famchars <- matrix(c(sibships[,p],rep(NA,length(sibships[,p]))),ncol=2,nrow=length(sibships[,p]))
  
  #Loop across sibship sizes in sample
  for(i in 1:nrow(Probabilities)){
    #Identify which famchars sib_size has the value of sibship at index i 
    index <- which(famchars[,1]==N[i])
    
    #Pseudo-randomly assign disease states to families according to probabilities expected at each sib size
    set.seed(sim_seed)
    famchars[index,2] <- sample(all_states,size=length(index), replace = TRUE,prob=Probabilities[i,])
  }
  
  #Determine the proportion of sample with familial sporadic and unaffected assignments, store in vector
  obsed <- c(
    length(which(famchars[,2]=="Familial"))/length(famchars[,2]),
    length(which(famchars[,2]=="Sporadic"))/length(famchars[,2]),
    length(which(famchars[,2]=="Unaffected"))/length(famchars[,2])
  )
  
#### Population has been defined - now estimate penetrance and produce summary
    #seed     =     seed used in randomisation
    #f        =     true penetrance
    #famchars = df containing full population
    #obsed stores disease state rates across population
    
    #Run the adpenetrance analysis for all disease state combinations
    for(s in 1:length(states)){
      
      #Calculate RX according to states represented in 'states' object:
        #Obsed[1] = familial
        #Obsed[2] = sporadic
        #Obsed[3] = unaffected
      if(states[s]=="fsu") {
        RX<- obsed[1] #No division needed because all states are represented
      } else if(states[s]=="fs") {
        RX<- obsed[1]/(obsed[1]+obsed[2])
      } else if(states[s]=="fu") {
        RX<- obsed[1]/(obsed[1]+obsed[3])
      } else if(states[s]=="su") {
        RX<- obsed[2]/(obsed[2]+obsed[3])
      } else if(states[s]=="au") {
        RX<- (obsed[1]+obsed[2])/(obsed[1]+obsed[2]+obsed[3])
      }
      
      #Summary structure for data:
        #c("Population","Tailoring","Seed","mean_N","RX","modif","States.Modelled","TruePenetrance","UncorrectedPen", "CorrectedPen")
      
      #Run Adpenetrance for the states modelled at the given f value
        #with Step 4 correction by error predicted across the default Poisson distribution sample generated in adpenetrance.errorfit
      adpen_estimate <- adpenetrance(N=obsN,RX=RX,states=states[s])
        
      #Store results in the summary matrix - define population name and that step-4 simulated sibs follow Poisson distribution
        statetemp[s,] <- c(pop,1,sim_seed[l],obsN,RX,0,s,f[l], adpen_estimate$output[3,], adpen_estimate$output[4,])
      
     
      #Run Adpenetrance for the states modelled at the given f value
        #with Step 4 correction by error predicted across the Tailored sample distribution,
        #generated in adpenetrance.errorfit according to the structure defined in define_sibstructure
      adpen_estimateTailored <- adpenetrance(N=obsN,RX=RX,states=states[s],define_sibstructure = sibstructure)
      
      #Store results in the summary matrix - define population name and that step-4 simulated sibs follow tailored distribution
      statetempTailored[s,] <- c(pop,2,sim_seed[l],obsN,RX,0,s,f[l], adpen_estimateTailored$output[3,], adpen_estimateTailored$output[4,])
      
    } #End states loop
  
    #Take the temporary matrices across states modelled and add to previous outputs in these loops
    if(is.na(result[1,1])){
      result <- plyr::rbind.fill.matrix(statetemp,statetempTailored)
    } else {
      result <- plyr::rbind.fill.matrix(result,statetemp,statetempTailored)
    }
    
    } #End loop across values of p
  print(l) #print values of l for indications of progress
  } #End loop across all test f values
    
  #Dataset is generated and stored in the 'results' object - now obtain useful information from this
  
  #Calculate the difference between true and estimated penetrance
  TrueMinusCorr <- result[,"TruePenetrance"]-result[,"CorrectedPen"]

  #Combine with difference column and convert into data.frame
  sum_diff<- as.data.frame(cbind(result,TrueMinusCorr)) 
  
  #Average the repeated seeds to get the mean difference in True penetrance minus estimated across states modelled at each value of true penetrance
  sum_diff_ave<- ddply(sum_diff,.(Population,Tailoring,TruePenetrance,States.Modelled),summarise,TrueMinusCorr=mean(TrueMinusCorr))
  
  #Switch the valence to obtain estimated penetrance minus true
  sum_diff_ave$CorrMinusTrue <- sum_diff_ave$TrueMinusCorr*-1
  
  #Convert to new object for unique plot call
  Sim1 <- sum_diff_ave
  
  #Write the data into file and define name by system time - commented to avoid unnecessary writing
  ##write.csv(Sim1,file=paste0("Simulation1_facet_grid_",Sys.Date(),".csv"))
  
  #Convert states modelled column into factors
  Sim1$States.Modelled <- as.factor(Sim1$States.Modelled)

  #Compute labels for each facet type
  sim_labs <- c(
    `1` = "Poisson",
    `2` ="Tailored"
  )
  
  pop_labs <- c(
    `1` = "UK (1974)",
    `2`  ="Next Steps"
  )
  
    ggplot(data=Sim1 ,aes(TruePenetrance,CorrMinusTrue))+
      geom_line(aes(colour=States.Modelled))+
      facet_grid(Population~Tailoring,labeller = labeller(.rows=pop_labs,.cols=sim_labs))+
      xlim(c(0,1))+
      ylab("Estimate Error (Estimate - True penetrance)")+
      xlab("True penetrance value")+
      geom_hline(yintercept=0,lty=2,alpha=0.75)+
      scale_colour_manual(name = 'States modelled',
                          values =c("deepskyblue3", #fsu
                                    "firebrick", #fs
                                    "black", #fu
                                    "forestgreen", #su
                                    "orange"), #au
                          labels =c('F, S, U ','F, S','F, U','S, U','A, U'))+
      theme(legend.position = "bottom",
            panel.background = element_rect(fill="white"),
            panel.grid = element_line(colour = "gray90"),
            strip.background = element_rect(fill="gray94"),
            panel.border = element_rect(fill=NA, colour="gray94")
      )

    ggsave(path="Supp_figures",filename=paste0("Figure_S4.png"),units="mm",width=210,height=150,dpi=300)
    #ggsave(path="Supp_figures",device='pdf',filename=paste0("Figure_S4.pdf"),units="mm",width=210,height=150,dpi=300)
    