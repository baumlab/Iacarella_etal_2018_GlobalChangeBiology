################################################################################
#Script for Iacarella et al. ............                                      #
#                                                                              #
#Null model assembly, calculation of BC beta-deviation,                        #
#and set-up of dataframes for within-disturbance category site comparisons     #
#defined by within- or across-region spatial scales                            #
#Null model assembly code adapted from Tucker et al. 2016 Oikos                #
################################################################################

library(vegan)
library(dplyr)

##Load dataset with all sites representing a provincial pool
#One column has site names, all following columns represent species counts (i.e., species are columns, rows are sites)
#E.g.
Null.fish<-read.csv("Oregonian.spp.comb.csv")
attach(Null.fish)

Null.fish<-read.csv("Aleutian.spp.comb.csv")
attach(Null.fish)

##Remove any additonal columns defining regions, etc.
#Leave only site names in first column, followed by species counts 
Null.fish2<-Null.fish[,c(-1:-3)]

##Turn dataframe into matrix, with site names column converted to row names
Null.fish2<-as.matrix(Null.fish2)
row.names(Null.fish2) <- Null.fish2[,1]
Null.fish2 <- Null.fish2[,-1]
class(Null.fish2)<-"numeric"

##Null model assembly and calculation of Bray-Curtis for null models, following Tucker et al. (2016)
#Site abundance and richness is randomized 
#Total count per species is fixed (i.e. rare species remain rare),but # of times a spp. occurs across sites is random

rand<-999
all_null_sim2<-array(0,dim=c(nrow(Null.fish2),nrow(Null.fish2),rand))
bbs.sp.site <- Null.fish2
#Generate null patches - fixes abundance for each SPECIES, but not for SITES
for (randomize in 1:rand) {  
  null.dist = Null.fish2
  for (species in 1:ncol(null.dist)) {
    tot.abund = sum(null.dist[,species])
    null.dist[,species] = 0
    for (individual in 1:tot.abund) {
      sampled.site = sample(c(1:nrow(bbs.sp.site)), 1)
      null.dist[sampled.site, species] = null.dist[sampled.site, species] + 1
    }}
  #Log transfrom matrix (sensu Anderson & Willis, 2003, Ecology)
  null.dist.ln <- log(null.dist+1)
  #Calculate Bray-Curtis dissimilarity for null model
  bucket_bray <- as.matrix(vegdist(null.dist.ln, "bray"))  #Expected beta-diversity, for one null matrix at a time
  diag(bucket_bray) <- NA
  #For pairwise Beta expected across all null matrices (to take ave and SD later)
  all_null_sim2[,,randomize]<- bucket_bray
} ## end randomize loop

##Calculate observed beta-diversity
Null.fish.ln<-log(Null.fish2+1)
beta_comm_abund <- vegdist(Null.fish.ln, "bray")
beta_obs <- as.matrix(as.dist(beta_comm_abund)) #beta(observed)
diag(beta_obs) <- NA

##Calculate BC beta-deviation
#Use mean and SD across all pairwise comparisons and null matrices (sensu Karp et al., 2012, Ecology Letters)
mean_null<-mean(all_null_sim2,na.rm = TRUE)
var_null<-var(all_null_sim2,na.rm = TRUE)
beta_species<-(beta_obs-mean_null)/(sqrt(var_null))

#Save BC beta-deviation measures
write.csv(beta_species,"OregZscore_ALL.csv")

##Set up columns defining site-site comparisons 
#Used for later bootstrapping of within-disturbance level comparisons defined by within- and across-region comparisons

library(reshape2)
#Change matrix of beta-deviation values to a column list of site-site comparisons
Beta.col<-melt(beta_species)[melt(lower.tri(beta_species))$value,]

#Set Region values to define within and across region comparisons
Region<-Null.fish[,1]

#Create matrix of Region values 
mat<-matrix(nrow=length(Region),ncol=length(Region),dimnames=list(Region,t(Region)))
x=matrix(data=NA,nrow=length(Region),ncol=length(Region))

for(i in 1:length(Region)) {
  for(j in 1:length(Region)) {
    x[i,j] <- paste(rownames(mat)[i],colnames(mat)[j])
  } }

#Convert matrix to column list of site-site comparisons, to match beta-diversity comparisons
Reg.col<-melt(x)[melt(lower.tri(x))$value,]
Reg.col<-as.data.frame(Reg.col[,'value'])

#Do the same for disturbance categories
Dist<-Null.fish[,3]
mat<-matrix(nrow=length(Dist),ncol=length(Dist),dimnames=list(Dist,t(Dist)))
y=matrix(data=NA,nrow=length(Dist),ncol=length(Dist))

for(i in 1:length(Dist)) {
  for(j in 1:length(Dist)) {
    y[i,j] <- paste(rownames(mat)[i],colnames(mat)[j])
  } }

Dist.col<-melt(y)[melt(lower.tri(x))$value,]
Dist.col<-as.data.frame(Dist.col[,'value'])

#Bind region and disturbance category columns
Reg.Dist<-cbind(Reg.col,Dist.col)
colnames(Reg.Dist)<-c("Reg.comp","Dist.comp")
#Bind with beta-diversity values
Beta_all<-cbind(Beta.col,Reg.Dist)

##Remove all comparisons across disturbance levels (e.g. high vs. low)
#Final dataframe is ALL site comparisons within same dist level
Beta_dist_all<-subset(Beta_all,Dist.comp!="high low" & Dist.comp!="low high" & Dist.comp!="low medium" & Dist.comp!="medium high" &  Dist.comp!="medium low" & Dist.comp!="high medium")
Beta_dist_all$Dist.comp<-as.character(Beta_dist_all$Dist.comp)
#Rename category comparisons
Beta_dist_all$Dist.comp[Beta_dist_all$Dist.comp=="low low"] <-"low"
Beta_dist_all$Dist.comp[Beta_dist_all$Dist.comp=="medium medium"] <-"medium"
Beta_dist_all$Dist.comp[Beta_dist_all$Dist.comp=="high high"] <-"high"

##Remove all comparisons within regions 
#Final dataframe is ACROSS region comparisons within same dist level
#E.g.
Beta_dist_cross<-subset(Beta_dist_all, !(Reg.comp %in% c("Victoria Victoria",
                                                         "Qualicum Qualicum","Comox Comox","Clayoquot Clayoquot",
                                                         "Fraser Fraser","GulfI GulfI","Barkley Barkley")))
#Save dataframe for later bootstrapping
#E.g.
write.csv(Beta_dist_cross,"Oreg_BetabyImpact_acrossRegions_Zscore.csv")

##Remove all comparisons across regions 
#final dataframe is WITHIN region comparisons within same dist level
#E.g.
Beta_dist_within<-subset(Beta_dist_all, Reg.comp %in% c("Victoria Victoria",
                                                        "Qualicum Qualicum","Comox Comox","Clayoquot Clayoquot",
                                                        "Fraser Fraser","GulfI GulfI","Barkley Barkley"))
#E.g.
write.csv(Beta_dist_within,"Oreg_BetabyImpact_withinRegions_Zscore.csv")