#####################################################
#Script for Iacarella et al. ............           #
#                                                   #
#Indicator species and fourth-corner trait analysis #
#Code for Table 2 and Figure 5                      #
#####################################################

##Load file with sites as rows and species as columns
#E.g.
AllReg.sp<-read.csv("Region.spp.comb.csv")

##Separate vector for disturbance
H.Dist<-AllReg.sp[,3]

##Remove columns with extra info, transform to matrix with sites as row names
All.Sp.m<-AllReg.sp[,c(-1:-3,-69)]
All.Sp.m<-as.matrix(All.Sp.m)
row.names(All.Sp.m) <- All.Sp.m[,1]
All.Sp.m <- All.Sp.m[,-1]
class(All.Sp.m)<-"numeric"

##Remove rare counts by species and by site

##Remove rare occurrences acros species (<10 individuals)
#Calculate counts for each species across sites
sp.sum<-apply(All.Sp.m,2,sum)
#Find no. of species with < 10 individuals 
sum(sp.sum < 10)  
#Order columns by most abundant spp
spe.sorted<-All.Sp.m[,order(sp.sum,decreasing = TRUE)] 
#Keep first 42 species, removes species that were counted less than 10 times
spe.abund<-spe.sorted[,c(1:42)] 

##Remove rare occurrences across sites (<5 sites)
#Calculate no. of sites each species occurred
spe.pres<-apply(spe.abund > 0, 2, sum)
#Find no. of species at less than 5 sites
sum(spe.pres < 5) 
#Order columns by spp with the most sites
spe.sorted2<-spe.abund[,order(spe.pres,decreasing = TRUE)]
#Keep first 36 species, removes species that were counted at fewer than 5 sites
spe.common<-spe.sorted2[,c(1:36)] 

########################################################################################\
##Indicator species analysis##

##Calculate indicator values by disturbance category
library(indicspecies)

##Set up dataframes for abundance and incidence measures
spe.common<-as.data.frame(spe.common)
spe.common.pa<-as.data.frame(ifelse(spe.common>0,1,0))

##Indicator analysis within disturbance categories, Table 2
#By abundance
indval<-multipatt(spe.common,H.Dist,duleg=TRUE,func="IndVal.g",control=how(nperm=999))
summary(indval, indvalcomp = TRUE) 
#By presence/absence
indval<-multipatt(spe.common.pa,H.Dist,duleg=TRUE,func="IndVal.g",control=how(nperm=999))
summary(indval, indvalcomp = TRUE) 

########################################################################################
##Trait analysis##

##Make sure species and sites are in the same order in all three datasets:
#species counts (prepared above), species traits, site characteristics (i.e. disturbance scores)
#species counts from above used here, with four species removed owing to missing trait information

##Load fish trait matrix - species are rows, traits are columns
Fish.traits<-read.csv("FishTraitsSimple.common.csv")
attach(Fish.traits)

#Remove uneeded columns, set species names to row names
Fish.traits2<-Fish.traits[,-c(2:6,12)]
row.names(Fish.traits2) <- Fish.traits2[,1]
Fish.traits2 <- Fish.traits2[,-1]
#Set any cateogorical triats to factors, categorical traits are dummy coded
Fish.traits2$Res<-as.factor(Fish.traits2$Res)
Fish.traits2$Rep.Guild<-as.factor(Fish.traits2$Rep.Guild)

##Load site disturbance scores - sites are rows
Site.dist<-read.csv("Dist.allsites.csv")
attach(Site.dist)
#Remove any uneeded columns
Site.dist2<-Site.dist[,c(2:7)]

##Run analysis
library(mvabund)

##Use method='glm1path' method for fourth-corner coeff fit with LASSO penalties
trait.test<-traitglm(spe.common, Site.dist2[,c(2)], Fish.traits2[,-c(4,5,7:10)],method="glm1path")
trait.test$fourth.corner

##Use method='manyglm' to test significance of the trait:environment (i.e. disturbance) interaction
trait.test2<-traitglm(spe.common, Site.dist2[,c(2)], Fish.traits2[,-c(4,5,7:10)],method="manyglm")
#Bootstrap significance test
trait.sig<-anova(trait.test2,nBoot=999)
trait.sig

##Graph fourth-corner coefficients (from 'glm1path'), Figure 5
library(lattice)
a        = max( abs(trait.test$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th = levelplot(t(as.matrix(trait.test$fourth.corner)), xlab="Anthropogenic disturbance",
                     ylab="Species traits",
                     col.regions=colort(100), at=seq(-.4, .4, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)



