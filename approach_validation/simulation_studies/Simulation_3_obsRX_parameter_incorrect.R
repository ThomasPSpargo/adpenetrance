#####
# Script for running a simulation study as part validation for the ADPenetrance (https://github.com/ThomasPSpargo/adpenetrance/) approach to calculate genetic penetrance
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk) 
# Please get in touch if you have any issues.
####

# This simulation tests the role of misspecifying the input parameter RX^obs upon penetrance estimation accuracy
# This analysis is performed five times, for variants with true penetrance of 0.1, 0.25, 0.50, 0.75, 1.0
# Each analysis is repeated 3 times per penetrance value and results are averaged across these triplicates

#The data from the simulation will be written into a .tsv file and a ggplot will be constructed, written into a .pdf file

#load necessary packages and adpenetrance function/subfunctions
library(plyr)
library(ggplot2)
library(reshape2)
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

#Store the sibstructures in a list
strucs<- list(sibstructureUK,sibstructureNS)

#Ground truth penetrance values (f):
f <- c(0.1,0.25,0.5,0.75,1) #Set the penetrance values to test
f <- sort(rep(f,3)) #Amplify to test each value 3 times, sort into ascending order

sims <- length(f)   #Extract number of penetrance loops to run

g=0 #Fix residual disease risk to 0

#Define required chr vectors:
all_states <- c("Familial", "Sporadic", "Unaffected") #Disease state levels to be assigned to families
states=c("fsu","fs","fu","su","au") #Define the disease states to model

#Define levels of modification to RX
modif= seq(-0.15,0.15,length=7) 

# #Prepare empty list to store output summary at each states modelled, also set colnames
output_colnames <- c("Population","Tailoring","Seed","mean_N","RX","modif","States.Modelled","TruePenetrance","UncorrectedPen", "CorrectedPen")
statetemp <- vector(mode="list",length=2) 
result <- NULL 

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
    #Model the poisson-tuned correction and correction tailored to specified sibstructure from strucs
    sibstructure <- list(NULL,strucs[[p]])
    pop <- p #label the population structure being tested
    
    #Identify mean N value in sample
    obsN <- sum(sibships[,p])/length(sibships[,p])
    
    #Return vector of unique sib_sizes represented in sample
    N <- sort(unique(sibships[,p]))
    
    
    ###Calculate the probability of each disease state expected each value of N at penetrance f
    Punasc = (1-f[l])*(1-f[l]/2-g/2)^N*(1-g)                                              #Probability unaffected
    Pspor  = f[l]*(1-f[l]/2-g/2)^N*(1-g)+
      (1-f[l])*N*(f[l]/2+g/2)*(1-f[l]/2-g/2)^(N-1)*(1-g)+
      (1-f[l])*(1-f[l]/2-g/2)^N*g                                                         #Probability sporadic/simplex
    Pfam   = 1-(Punasc+Pspor)                                                         #Probability familial (= 1 - Punasc - Pspor)
    
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
      
      #Loop across values of modification defined in modif - i.e. manipulate misspecification of RX
      for(n in 1:length(modif)){
        
        #Derive the modifiedRX value to be used
        modifiedRX <-RX+modif[n]
        #Truncate the rate to fit within the feasible range of [0,1]
        if(modifiedRX<0){
          modifiedRX <- 1e-10
        } else if (modifiedRX>1) {
          modifiedRX <- 1
        }
        
        #Run ADPenetrance, loop across the poisson and tailored penetrance correction
        for(t in 1:2){
          suppressMessages(
            adpen_estimate <- adpenetrance(N=obsN,
                                           RX=modifiedRX,
                                           states=states[s],
                                           define_sibstructure = sibstructure[[t]]
            )
          )
          
          #Store results for poisson and tailored in the summary list
          #c("Population","Tailoring","Seed","mean_N","RX","modif","States.Modelled","TruePenetrance","UncorrectedPen", "CorrectedPen")
          statetemp[[t]] <- c(pop,t,sim_seed[l],obsN,RX,modif[n],s,f[l], adpen_estimate$output[3,], adpen_estimate$output[4,])
          
        }
        result <- rbind(result,statetemp[[1]],statetemp[[2]]) #Combine temporary vectors with other results
      } #End the modification loop
    } #End states loop
  } #End loop across values of p
  print(l) #print values of l for indications of progress
} #End loop across all test f values

#Dataset is generated and stored in the 'results' object - now obtain useful information from this
colnames(result) <- output_colnames

#Calculate the difference between true and estimated penetrance (corrected and uncorrected), then bind to main dataframe
CorrMinusTrue <- result[,"CorrectedPen"]-result[,"TruePenetrance"]
UncorrMinusTrue <- result[,"UncorrectedPen"]-result[,"TruePenetrance"]
sum_diff<- as.data.frame(cbind(result,CorrMinusTrue,UncorrMinusTrue)) 

#Average the repeated seeds to get the mean difference in True penetrance minus estimated across states modelled at each value of true penetrance
Sim3 <- ddply(sum_diff,.(Population,Tailoring,TruePenetrance,States.Modelled,modif),summarise,CorrMinusTrue=mean(CorrMinusTrue),UncorrMinusTrue=mean(UncorrMinusTrue))

#Compute labels for each facet type ('no correction' level will be established by melting the df)
#Convert states modelled column into factors with appropriate names
Sim3$States.Modelled <- factor(Sim3$States.Modelled,levels= c(1:5),labels= c("F, S, U",
                                                                             "F, S",
                                                                             "F, U",
                                                                             "S, U",
                                                                             "A, U"))

Sim3$Population <- factor(Sim3$Population,levels= c(1:2),labels= c("UK (1974)",
                                                                   "Next Steps"))

Sim3$Tailoring <- factor(Sim3$Tailoring,levels= c(1:3),labels= sim_labs <- c("Poisson",
                                                                             "Tailored",
                                                                             "No correction"))

#Melt so that the corrected and uncorrected penentrance estiamtes are in the same column
mSim3<- reshape2::melt(Sim3,measure.vars=c("CorrMinusTrue","UncorrMinusTrue"))

#Unadjusted estimates are duplicated (i.e. Poisson and Tailored show the same, therefore remove any 'tailored' uncorrected rows)
mSim3 <- mSim3[-which(mSim3$variable=="UncorrMinusTrue" & mSim3$Tailoring=="Tailored"),] 

#Then relabel uncorrected into a third 'tailoring' category
mSim3$Tailoring[which(mSim3$variable=="UncorrMinusTrue")] <- "No correction"

fig <- ggplot(data=mSim3 ,aes(modif,value,colour=factor(TruePenetrance)))+
  geom_line()+
  facet_grid(States.Modelled~Population+Tailoring)+
  ylab("Estimate Error (Estimate - True penetrance)")+
  xlab(bquote("Degree of error in estimate of R(X)"^"obs"))+
  geom_hline(yintercept=0,lty=2,alpha=0.75)+
  scale_colour_manual(name = 'True Penetrance',
                      values =c("#e3d16e", #0.1
                                "#ff7943", #0.25
                                "#be3044", #0.5
                                "#796c5d", #0.75
                                "#202547"), #1
  )+  #https://jaredhuling.org/jcolors/ colour palette pal5
  theme(legend.position = "bottom",
        panel.background = element_rect(fill="white"),
        panel.grid = element_line(colour = "gray90"),
        strip.background = element_rect(fill="gray94"),
        panel.border = element_rect(fill=NA, colour="gray94")
  )

ggsave(filename="Simulation_3.pdf",plot=fig,device='pdf',units="mm",width=210,height=210,dpi=300)
write.table(mSim3,file="Simulation_3.tsv",sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE) 

