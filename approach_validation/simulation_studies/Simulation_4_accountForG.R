#####
# Script for running a simulation study as part validation for the ADPenetrance (https://github.com/ThomasPSpargo/adpenetrance/) approach to calculate genetic penetrance
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk) 
# Please get in touch if you have any issues.
####

# This simulation tests the role of residual risk g upon penetrance estimation according to the relative rarity of g.
# Families are generated using a model whereby people with no variant have g= (0, ..., 0.2) disease risk
# Analysis with ADPenetrance is run taking into account values of g and when g is assumed to =0
# The result plot examines differences in penetrance estimation accuracy according to the true value of g
# This is run four times, for variants with true penetrance of 0.25, 0.50, 0.75, 1.0
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
f <- seq(0.25,1,length=4) #Set the penetrance values to test
f <- sort(rep(f,3)) #Amplify to test each value 3 times, sort into ascending order
sims <- length(f)   #Extract number of penetrance loops to run

#Indicate tested values of g #Error is so small that it is pointless going smaller than 0.001
useG <- c(0,0.001,0.025,0.050,0.075,0.1,0.125,0.150,0.175,0.2)

#Define required chr vectors:
all_states <- c("Familial", "Sporadic", "Unaffected") #Disease state levels to be assigned to families
states=c("fsu","fs","fu","su","au") #Define the disease states to model

# #Prepare empty list to store output summary at each states modelled and set colnames
output_colnames<- c("Population","Tailoring","Seed","mean_N","RX","modif","States.Modelled","TruePenetrance","BaseUncorrectedPen", "BaseCorrectedPen","ExtendedUncorrectedPen","ExtendedCorrectedPen")
statetemp <- vector(mode="list",length=2) 
result <- NULL

sim_seed <- sample(1:30000,sims) #Generate random seeds for use in each loop

#Indicate the number of families to generate
numfamilies <- 90000

#Loop across penetrance values
for(l in 1:sims){
  #Define seed for the randomisation
  #Generate population of 90000 families resembling UK structure
  set.seed(sim_seed[l])
  sibshipsUK <-sample(sibstructureUK[,1], numfamilies,replace=TRUE,prob=sibstructureUK[,2])
  
  #Generate population of 90000 families resembling NS structure
  set.seed(sim_seed[l])
  sibshipsNS <-sample(sibstructureNS[,1], numfamilies,replace=TRUE,prob=sibstructureNS[,2])
  
  #Store populations in 2 column matrix
  sibships <- matrix(c(sibshipsUK,sibshipsNS),
                     nrow=numfamilies)
  
  #Loop across simulated population
  for(p in 1:ncol(sibships)){
    #Model the poisson-tuned correction and correction tailored to specified sibstructure from strucs
    sibstructure <- list(NULL,strucs[[p]])
    pop <- p #label the population structure being tested
    
    #Identify mean N value in sample
    obsN <- sum(sibships[,p])/numfamilies
    
    #Return vector of unique sib_sizes represented in sample
    N <- sort(unique(sibships[,p]))
    
    #Generate family probabilities according to G
    for(j in 1:length(useG)){
      g <- useG[j] #Specify the value of g to use
      
      #Calculate disease state proportions under the extended model, considering N and penetrance f, with residual disease risk g
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
      famchars <- matrix(c(sibships[,p],rep(NA,numfamilies)),ncol=2,nrow=numfamilies)
      
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
        length(which(famchars[,2]=="Familial"))/numfamilies,
        length(which(famchars[,2]=="Sporadic"))/numfamilies,
        length(which(famchars[,2]=="Unaffected"))/numfamilies
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
        
        #Set a nominally small value for RX if it's determined to be 0
        if((RX)==0){
          RX <- 1e-10
        }
        
        #Run ADPenetrance, with and without consideration for g, loop across the poisson and tailored penetrance correction
        considerG <- list(noG=NULL,withG=NULL)
        for(t in 1:2){
          
          suppressMessages(
            considerG[[1]] <- adpenetrance(N=obsN,
                                           RX=RX,states=states[s],
                                           define_sibstructure = sibstructure[[t]],
                                           useG = 0
            )
          )
          suppressMessages(
            considerG[[2]] <- adpenetrance(N=obsN,
                                           RX=RX,states=states[s],
                                           define_sibstructure = sibstructure[[t]],
                                           useG = g
            )
          )
          
          
          #Store results for poisson and tailored in the summary list
          c("Population","Tailoring","Seed","mean_N","RX","modif","States.Modelled","TruePenetrance","BaseUncorrectedPen", "BaseCorrectedPen","ExtendedUncorrectedPen","ExtendedCorrectedPen")
          statetemp[[t]] <- c(pop,t,sim_seed[l],obsN,RX,g,s,f[l],
                              considerG$noG$output[3,], considerG$noG$output[4,],
                              considerG$withG$output[3,], considerG$withG$output[4,])
        }
        result <- rbind(result,statetemp[[1]],statetemp[[2]]) #Combine temporary vectors with other results
      } #End states loop
    } #end loop across values of g      
  } #End loop across values of p
print(l) #print values of l for indications of progress
} #End loop across all test f values

#Dataset is generated and stored in the 'results' object - now obtain useful information from this
colnames(result) <- output_colnames

#Estimate values
#"TruePenetrance","BaseUncorrectedPen", "BaseCorrectedPen","ExtendedUncorrectedPen","ExtendedCorrectedPen"

#In base model, calculate difference between uncor and corr estimates and true pen
BaseUncorrMinusTrue <- result[,"BaseUncorrectedPen"]-result[,"TruePenetrance"]
BaseCorrMinusTrue <- result[,"BaseCorrectedPen"]-result[,"TruePenetrance"]

#Repeat in extended model
ExtUncorrMinusTrue <- result[,"ExtendedUncorrectedPen"]-result[,"TruePenetrance"]
ExtCorrMinusTrue <- result[,"ExtendedCorrectedPen"]-result[,"TruePenetrance"]

#Combine with difference column and convert into data.frame
sum_diff<- as.data.frame(cbind(result,
                               BaseUncorrMinusTrue,BaseCorrMinusTrue,ExtUncorrMinusTrue,ExtCorrMinusTrue
                               )) 

#Average the repeated seeds to get the mean difference in True penetrance minus estimated across states modelled at each value of true penetrance
Sim4 <- ddply(sum_diff,.(Population,Tailoring,TruePenetrance,States.Modelled,modif),summarise,
              BaseUncorrMinusTrue=mean(BaseUncorrMinusTrue),BaseCorrMinusTrue=mean(BaseCorrMinusTrue),ExtUncorrMinusTrue=mean(ExtUncorrMinusTrue),ExtCorrMinusTrue=mean(ExtCorrMinusTrue)
)

#Compute labels for each facet type ('no correction' level will be established by melting the df)
Sim4$States.Modelled <- factor(Sim4$States.Modelled,levels= c(1:5),labels= c("F, S, U",
                                                                             "F, S",
                                                                             "F, U",
                                                                             "S, U",
                                                                             "A, U"))

Sim4$Population <- factor(Sim4$Population,levels= c(1:2),labels= c("UK (1974)",
                                                                   "Next Steps"))

Sim4$Tailoring <- factor(Sim4$Tailoring,levels= c(1:3),labels= sim_labs <- c("Poisson",
                                                                             "Tailored",
                                                                             "No correction"))

#Melt so that the corrected and uncorrected penentrance estiamtes are in the same column
#mSim4<- reshape2::melt(Sim4,measure.vars=c("CorrMinusTrue","UncorrMinusTrue"))
mSim4<- reshape2::melt(Sim4,measure.vars=c("BaseCorrMinusTrue","BaseUncorrMinusTrue","ExtUncorrMinusTrue","ExtCorrMinusTrue"))

#Define whether penetrance estimate is made using a model assuming no g or including g
mSim4$withG<- NULL
mSim4$withG[grepl("Ext",mSim4$variable)] <- "g = x"
mSim4$withG[grepl("Base",mSim4$variable)] <- "g = 0"
  
#Unadjusted estimates are duplicated (i.e. Poisson and Tailored show the same, therefore remove any 'tailored' uncorrected rows)
mSim4 <- mSim4[-which(grepl("UncorrMinusTrue",mSim4$variable) & mSim4$Tailoring=="Tailored"),] 

#Then relabel uncorrected into a third and fourth 'tailoring' category
mSim4$Tailoring[grepl("UncorrMinusTrue",mSim4$variable)] <- "No correction"

fig <- ggplot(data=mSim4, aes(modif,value,colour=factor(TruePenetrance)))+
  geom_line()+
  facet_grid(withG+States.Modelled~Population+Tailoring)+
  ylab("Estimate Error (Estimate - True penetrance)")+
  xlab("Residual disease risk in population")+
  geom_hline(yintercept=0,lty=2,alpha=0.75)+
  scale_colour_manual(name = 'True Penetrance',
                      values =c(#"#e3d16e", #0.1
                                "#ff7943", #0.25
                                "#be3044", #0.5
                                "#796c5d", #0.75
                                "#202547"), #1
  )+  #https://jaredhuling.org/jcolors/ colour palette pal5
  theme(legend.position = "bottom",
        panel.background = element_rect(fill="white"),
        panel.grid = element_line(colour = "gray90"),
        strip.background = element_rect(fill="gray94"),
        panel.border = element_rect(fill=NA, colour="gray94"),
        axis.text.x = element_text(angle=90,vjust = 0.5)
  )

ggsave(filename="Simulation_4.pdf",plot=fig,device='pdf',units="mm",width=210,height=250,dpi=300)
write.table(mSim4,file="Simulation_4.tsv",sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE) 
