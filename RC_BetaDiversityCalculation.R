###########################################################################################
#Script for Iacarella et al. ............                                                 #
#                                                                                         #
#All code and explanation is from Chase et al. (2011) Ecosphere                           #
#EXCEPT: modified code Lines 40-51 to account for different # of sites/disturbance        #
#category that make up the null pool                                                      #
#Accounting for site # discrepancy follows method by Velland et al. (2007) J. of Ecology  #
###########################################################################################

##Dataframe has sites as rows, species as columns
Null.fish<-read.csv("Oregonian.spp.comb.csv")
attach(Null.fish)

#Call Raup Crick Function
source("~/XX/RaupCrickFunction_Chase_etal2011.R")

##Remove any additonal columns defining regions, etc.
#Leave only site names in first column, followed by species counts 
Null.fish2<-Null.fish[,c(-1:-3)]

##Turn dataframe into matrix, with site names column converted to row names
Null.fish2<-as.matrix(Null.fish2)
row.names(Null.fish2) <- Null.fish2[,1]
Null.fish2 <- Null.fish2[,-1]
class(Null.fish2)<-"numeric"

##Obtain Raup Crick dissimilarity matrix (corrected for alpha div)
#Note: added in function calculation of mean occurrence frequencies by disturbance category (sensu Vellend et al. 2007)
#Dist_bysite is added function call to define sites by disturbance category
Beta_RC_class<-raup_crick(Null.fish2,Dist_bysite,plot_names_in_col1=FALSE, classic_metric=TRUE, 
                          split_ties=TRUE,reps=999, set_all_species_equal=FALSE, 
                          as.distance.matrix=TRUE, report_similarity=FALSE)

beta_RC <- as.matrix(as.dist(Beta_RC_class)) 
diag(beta_RC) <- NA

##Re-attached columns with region names, etc.
beta_RC2<-cbind(Null.fish[,c(1,3)],beta_RC)
#Save
write.csv(beta_RC2,"XX")

##Follow procedure in 'BC_BetaDeviationCalculation.R' to set up columns defining site-site comparisons 
#Used for later bootstrapping of within-disturbance level comparisons defined by within- and across-region comparisons
#In 'Bootstrap&Graph_Fig2.R'
