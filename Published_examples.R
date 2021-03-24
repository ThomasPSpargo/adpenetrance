#Version updated 19/01/2021

#NOTE:
#Before using this function, please first consult the README.md documentation in the associated github repository:
#https://github.com/thomaspspargo/adpenetrance


#---------------------------------------------------------#
#       ADPenetrance: Penetrance calculator function      #
#---------------------------------------------------------#
#Set your working directory

#setwd("/Users/tom/OneDrive - King's College London/PhD/PhD project/Sibship_penetrance/Git_files")
#getwd()
#Download and run the script containing the adpenetrance tool 
source("Penetrance_function.R") 

#--------------------------------------------#
#     Examples modelled in publication       #
#--------------------------------------------#
#See the publication for full description of these examples
#Only the data used in each instance is shown here
#Calculations can also be performed without the error terms 

#Case 1 - Parkinsons and LRRK2 (p.G2019S)
#Various combinations of disease state data used for this case

#States modelled: familial, sporadic, unaffected
adpenetrance(N=1.646, 
               MF=201/5123,MF_SE=sqrt(201/5123*(1-201/5123)/5123),
               MS=179/14253,MS_SE=sqrt(179/14253*(1-179/14253)/14253),
               MU=11/14886,MU_SE=sqrt(11/14886*(1-11/14886)/14866),
               PF=.105, PA=1/37
               )

#P(fam) = 0.0982 (0.0593, 0.1371)
#f      = 0.3360 (0.2581, 0.4011)


#States modelled: familial, sporadic
adpenetrance(N=1.646, #Weighted aggregation of named regions in 2018
               MF=201/5123,MF_SE=sqrt(201/5123*(1-201/5123)/5123),
               MS=179/14253,MS_SE=sqrt(179/14253*(1-179/14253)/14253),
               PF=.105
               )

#P(fam)=  0.2682 (0.2292, 0.3072)
#f=       0.4532 (0.3935, 0.5111)


#States modelled: familial, unaffected
adpenetrance(N=1.646,
               MF=201/5123,MF_SE=sqrt(201/5123*(1-201/5123)/5123),
               MU=11/14886,MU_SE=sqrt(11/14886*(1-11/14886)/14866),
               PF=.105, PA=1/37
               )

#P(fam) = 0.1341 (0.0637, 0.2045)
#f      = 0.3055 (0.2205, 0.3676)


#States modelled: sporadic, unaffected
adpenetrance(N=1.646,
               MS=179/14253,MS_SE=sqrt(179/14253*(1-179/14253)/14253),
               MU=11/14886,MU_SE=sqrt(11/14886*(1-11/14886)/14866),
               PF=.105, PA=1/37
                 )

#P(spor) = 0.2970 (0.1699, 0.4241)
#f       = 0.1960 (0.1032, 0.3055)



#Case 2 - PAH and BMRP2
#in all cases, states modelled: familial, sporadic
#Different datasets used in each version of this calculation

#Evans et al - BMRP2
adpenetrance(N=1.543,
               MF=202/247, MF_SE=sqrt(202/247*(1-202/247)/247),
               MS=200/1174, MS_SE=sqrt(200/1174*(1-200/1174)/1174),
               PF=.055
               )

#P(fam) = 0.2184 (0.1946, 0.2422)
#f      = 0.3951 (0.3560, 0.4333)


#Aldred et al - BMRP2 (any variant)
adpenetrance(N=1.543,
               MF=40/58, MF_SE=sqrt(40/58*(1-40/58)/58),
               MS=26/126, MS_SE=sqrt(26/126*(1-26/126)/126),
               PF=.055
               )

#P(fam) = 0.1628 (0.1106, 0.2151)
#f      = 0.3025 (0.2108, 0.3898)


#Aldred et al - BMRP2 (small variants)
adpenetrance(N=1.543,
               MF=33/58, MF_SE=sqrt(33/58*(1-33/58)/58),
               MS=20/126, MS_SE=sqrt(20/126*(1-20/126)/126),
               PF=.055
               )

#P(fam) = 0.1726 (0.1069, 0.2383)
#f      = 0.3191 (0.2042, 0.4272)


#Aldred et al - BMRP2 (structural variants)
adpenetrance(N=1.543,
               MF=7/58, MF_SE=sqrt(7/58*(1-7/58)/58),
               MS=6/126, MS_SE=sqrt(6/126*(1-6/126)/126),
               PF=.055
               )

#P(fam) = 0.1285 (0.0115, 0.2456)
#f      = 0.2429 (0.0230, 0.4388)


#Case 3 - ALS and SOD1
#In all cases, states modelled: familial, sporadic
#Modelled independently for European and Asian populations
#Error terms for input derived from confidence intervals reported within meta-analysis from which values are drawn

#SOD1 (Asian) -
adpenetrance(N=1.823,
               MF=.3, MF_SE=(.3-.251)/1.96,
               MS=.015, MS_SE=(.015-.010)/1.96,
               PF=.05
               )

#P(fam) = 0.5128 (0.4201, 0.6056)
#f      = 0.7486 (0.6291, 0.8640)


#SOD1 (European)
adpenetrance(N=1.543,
               MF=.148, MF_SE=(.148-.115)/1.96,
               MS=.012, MS_SE =(.012-.007)/1.96,
               PF=.05)

#P(fam) = 0.3936 (0.2808, 0.5064)
#f      = 0.6597 (0.4938, 0.8123)


#Case 4: ALS and C9orf72
#In all cases, states modelled: familial, sporadic
#Modelled independently for European and Asian populations
#Error terms for input derived from confidence intervals reported within meta-analysis from which values are drawn

#C9orf72 (Asian)
adpenetrance(N=1.823,
               MF=.04, MF_SE= (.04-.02)/1.96,
               MS=.01, MS_SE =(.01-.00)/1.96,
               PF=.05)

#P(fam) = 0.1739 (0.0133, 0.3345)
#f      = 0.2824 (0.0230, 0.5142)


#C9orf72 (European)
adpenetrance(N=1.543,
               MF=.32, MF_SE=(.32-.28)/1.96,
               MS=.05, MS_SE =(.05-.04)/1.96,
               PF=.05, Zout=1.96)

#P(fam) = 0.2520 (0.2075, 0.2964)
#f      = 0.4489 (0.3773, 0.5176)





#Supplementary analyses:
#ALS and named SOD1 variants
#A5V
adpenetrance(N=1.543,
             MF=7/1125, MF_SE=sqrt(7/1125*(1-7/1125)/1125),
             MS=1/4366, MS_SE=sqrt(1/4366*(1-1/4366)/4366),
             PF=.05
             )

#P(fam) = 0.588 (0.081, 1.096)
#f      = 0.916 (0.157, 1.000*)
#*truncated as described in publication


#D91A (including gnomad data)
adpenetrance(N=1.543,
             MS=4/4366, MS_SE=sqrt(4/4366*(1-4/4366)/4366),
             MU=33/24143, MU_SE=sqrt(33/24143*(1-33/24143)/24123),
             PF=.05, PA=1/400
             )

#P(fam) = 1.59x10^-3 (0.000, 3.24x10^-3)
#f      = 0.0009, (0.0000, 0.0018)


#I114T
adpenetrance(N=1.543,
             MF=17/1140, MF_SE=sqrt(17/1140*(1-17/1140)/1140),
             MS=6/4366, MS_SE=sqrt(6/4366*(1-6/4366)/4366),
             PF=.05
             )

#P(fam) = 0.364 (0.149, 0.578)
#f      = 0.617 (0.278, 0.904)

