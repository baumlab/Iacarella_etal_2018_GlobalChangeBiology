############################################################
#Script for Iacarella et al. ............                  #
#                                                          #
#Gower distances for physical and biotic variables         #
############################################################

##Distances for retained variables across best-fit models from dbRDA
#Load all data from both zoogeographic provinces
Env.dist<-read.csv("Region.env.comb.csv")
attach(Env.dist)

##Select biotic variables retained in best-fit models from dataframe
Biot.dist<-Env.dist[,c(1:3,13,17:19)]
#remove missing rows with missing data - missing most eelgrass characteristics for some sites
Biot.data<-Biot.dist[!apply(is.na(Biot.dist), 1, any),]
Biot.data2<-Biot.data[,c(-1:-3)]

##Calculate Gower distances
Biot.dist2<-daisy(Biot.data2,"gower")

#Try removing tidal range
Biot.data3<-Biot.data2[,-4]
Biot.dist3<-daisy(Biot.data3,"gower")

##Select physical variables retained in best-fit models from dataframe
Phys.data<-Env.dist[,c(1:3,26:28,31)] #leave in missing data, only missing data for salinity

##Calculate Gower distances
Phys.dist2<-daisy(Phys.data,"gower")

##Set up dataframe for bootstrapping ('Bootstrap&Graph_Fig2.R')
Dist2<-as.data.frame(Biot.dist2)
Dist2<-cbind(Biot.data[,2],Dist2)
Dist3<-as.matrix(Dist2)
row.names(Dist3) <- Dist3[,1]
Dist3 <- Dist3[,-1]
class(Dist3)<-"numeric"

##Save and repeat physical variable distances
write.csv(Dist3,"Biot.Gower.csv")

##Add region and disturbance category information to select within-disturbance category and 
#within/across region site comparisons
#Same steps as with BC and RC site comparisons - 'BCBetaDeviationCalculation.R'
library(reshape2)
Phys_col<-melt(Phys.dist2)[melt(lower.tri(Phys.dist2))$value,]
Biot_col<-melt(Biot.dist2)[melt(lower.tri(Biot.dist2))$value,]

Region<-Biot.data[,1]
mat<-matrix(nrow=length(Region),ncol=length(Region),dimnames=list(Region,t(Region)))
x=matrix(data=NA,nrow=length(Region),ncol=length(Region))

for(i in 1:length(Region)) {
  for(j in 1:length(Region)) {
    x[i,j] <- paste(rownames(mat)[i],colnames(mat)[j])
  } }

Reg.col<-melt(x)[melt(lower.tri(x))$value,]
Reg.col<-as.data.frame(Reg.col[,'value'])

Dist<-Biot.data[,3]
mat<-matrix(nrow=length(Dist),ncol=length(Dist),dimnames=list(Dist,t(Dist)))
y=matrix(data=NA,nrow=length(Dist),ncol=length(Dist))

for(i in 1:length(Dist)) {
  for(j in 1:length(Dist)) {
    y[i,j] <- paste(rownames(mat)[i],colnames(mat)[j])
  } }

Dist.col<-melt(y)[melt(lower.tri(x))$value,]
Dist.col<-as.data.frame(Dist.col[,'value'])
Reg.Dist<-cbind(Reg.col,Dist.col)
colnames(Reg.Dist)<-c("Reg.comp","Dist.comp")

Phys_all<-cbind(Phys_col,Reg.Dist)
Biot_all<-cbind(Biot_col,Reg.Dist)

##Remove all comparisons across disturbance levels (e.g. high low)
#Final dataframe is ALL site comparisons within same dist level
Physdist_all<-subset(Phys_all,Dist.comp!="high low" & Dist.comp!="low high" & Dist.comp!="low medium" & Dist.comp!="medium high" &  Dist.comp!="medium low" & Dist.comp!="high medium")
Physdist_all$Dist.comp<-as.character(Physdist_all$Dist.comp)
Physdist_all$Dist.comp[Physdist_all$Dist.comp=="low low"] <-"low"
Physdist_all$Dist.comp[Physdist_all$Dist.comp=="medium medium"] <-"medium"
Physdist_all$Dist.comp[Physdist_all$Dist.comp=="high high"] <-"high"

##Remove all comparisons within regions 
#Final dataframe is ACROSS region comparisons within same disturbance category
#E.g.
PhysDist_cross<-subset(Physdist_all, !(Reg.comp %in% c("Victoria Victoria","Qualicum Qualicum","Comox Comox","Clayoquot Clayoquot",
                                                       "Fraser Fraser","GulfIslands GulfIslands","Barkley Barkley","Skeena Skeena","CCoast CCoast","Skeena Victoria","Skeena Qualicum","Skeena Comox","Skeena Clayoquot",
                                                       "Skeena Fraser","Skeena GulfIslands","Skeena Barkley","CCoast Victoria","CCoast Qualicum","CCoast Comox","CCoast Clayoquot",
                                                       "CCoast Fraser","CCoast GulfIslands","CCoast Barkley","Victoria CCoast","Qualicum CCoast","Comox CCoast","Clayoquot CCoast",
                                                       "Fraser CCoast","GulfIslands CCoast","Barkley CCoast","Victoria Skeena","Qualicum Skeena","Comox Skeena","Clayoquot Skeena",
                                                       "Fraser Skeena","GulfIslands Skeena","Barkley Skeena")))

##Remove all comparisons across regions & across disturbance levels
#Final dataframe is WITHIN region comparisons within same disturbance category
#E.g.
PhysDist_within<-subset(Physdist_all, Reg.comp %in% c("Victoria Victoria",
                                                      "Qualicum Qualicum","Comox Comox","Clayoquot Clayoquot",
                                                      "Fraser Fraser","GulfIslands GulfIslands","Barkley Barkley","Skeena Skeena","CCoast CCoast"))

##Save and load into 'Bootstrap&Graph_Fig2.R'
write.csv(PhysDist_within,"PhysDist_ALLReg_within.csv")
write.csv(PhysDist_cross,"PhysDist_ALLReg_across.csv")
