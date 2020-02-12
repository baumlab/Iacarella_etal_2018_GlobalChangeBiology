#################################################################
#Script for Iacarella et al. ............                       #
#                                                               #
#DbRDA, variance partitioning, Table 1, Figure 4, and LMM code  #
#################################################################

##Load dataframe for analysis
#Dataframe is matrix of site-site beta-diversity comparisons
Oreg.z<-read.csv("OregZscore_RDA_ALL.csv")

#Turn into matrix with site names as row names
Null.fish2<-as.matrix(Oreg.z)
row.names(Null.fish2) <- Null.fish2[,1]
Null.fish2 <- Null.fish2[,-1]
class(Null.fish2)<-"numeric"

#Convert to distance matrix
Oreg.fish.Z<-as.dist(Null.fish2)

##Load biotic covariates
Oreg.eelg<-read.csv("Oreg_Meadow_Char.csv")
#For analysis of sepaarate disturbance categories, extract sites needed for analysis
Env.data.sep<-subset(Oreg.eelg, H.Disturbance %in% c("low"))
#For all sites
Env.data.sep<-Oreg.eelg

##Load physical covariates
Oreg.phys<-read.csv("Oreg_Physical_Char.csv")
#For analysis of sepaarate disturbance categories, extract sites needed for analysis
Phys.data.sep<-subset(Oreg.phys, H.Disturbance %in% c("low"))
#For overall region
Phys.data.sep<-Oreg.phys

##Standardize continuous variables
Temp.insitu.z<-scale(Phys.data$Temp.insitu,center=TRUE,scale=TRUE)
Sal.insitu.z<-scale(Phys.data$Sal.insitu,center=TRUE,scale=TRUE)
SST.July.z<-scale(Phys.data$SST.July,center=TRUE,scale=TRUE)
SST.annual.z<-scale(Phys.data$SST.annual,center=TRUE,scale=TRUE)
Tide.z<-scale(Phys.data$Tide_RMS,center=TRUE,scale=TRUE)
River.z<-scale(Phys.data$River_calc,center=TRUE,scale=TRUE)

Shoot.z<-scale(Env.data$Shoot.ave,center=TRUE,scale=TRUE)
Biom.z<-scale(Env.data$Eel.ave,center=TRUE,scale=TRUE)
Epi.z<-scale(Env.data$Epi.lia,center=TRUE,scale=TRUE)

Pop.z<-scale(Phys.data$Wtshd.pop,center=TRUE,scale=TRUE)
Struct.z<-scale(Phys.data$Struct.cover,center=TRUE,scale=TRUE)
Shore.z<-scale(Phys.data$Shore.mod,center=TRUE,scale=TRUE)
H.add.z<-Pop.z+Struct.z+Shore.z

##Check correlations among variables
pairs.test<-as.data.frame(cbind(Shoot.z,Epi.z,Biom.z,Sal.insitu.z,Temp.insitu.z,SST.July.z,SST.annual.z,Tide.z,River.z,Pop.z,Struct.z,Shore.z,H.add.z))
colnames(pairs.test)<-c("Shoot.z","Epi.z","Biom.z","Sal.insitu.z","Temp.insitu.z","SST.July.z","SST.annual.z","Tide.z","River.z","Pop.z","Struct.z","Shore.z","H.add.z")
corr.test<-as.data.frame(cor(pairs.test,use="pairwise.complete.obs"))

##Bind standardized variables to original dataframes
Env.data.z<-cbind(Env.data,Shoot.z,Biom.z,Epi.z)
Phys.data.z<-cbind(Phys.data,Sal.insitu.z,SST.July.z,SST.annual.z,Tide.z,River.z)
HDist.data.z<-as.data.frame(cbind(Pop.z,Struct.z,Shore.z,H.add.z))
colnames(HDist.data.z)<-c("Pop.z","Struct.z","Shore.z","H.add.z")

library(vegan)

##########Distance-based redundancy analysis, Table 1##########
#Run dbRDA process on disturbance ("all" data only), physical, and biotic variable components separately
#Run for BC and RC beta-diversity metrics separately

##Anthropogenic disturbance model
##Run full model
Reg.comp.all<-capscale(Oreg.fish.Z~Pop.z + Shore.z + Struct.z, HDist.data.z, add=TRUE) 
anova(Reg.comp.all,step=999,perm.max=999) #should be significant before running selection process
RsquareAdj(Reg.comp.all)$adj.r.squared #this r-squared should not be exceeded

##Check for collinearity - avoid VIF>4
#Remove variable with highest VIF if >4, re-run and check, repeat if needed
vif.cca(Reg.comp.all)

##Forward selection process
step.forward<-ordistep(capscale(Oreg.fish.Z~1,HDist.data.z,add=TRUE),
                       scope=formula(Reg.comp.all),direction="forward",pstep=1000)

##Check r-squared of each model
RsquareAdj(capscale(Oreg.fish.Z~ Shore.z,HDist.data.z,add=TRUE))$adj.r.squared #check if exceeds global r-sq

#Run best-fit model
Reg.dist<-capscale(Oreg.fish.Z~ Shore.z + Struct.z,HDist.data.z,add=TRUE)
anova(Reg.dist,step=999,perm.max=999)

##Make dataframe with selected variables
hdist.pars<-cbind(HDist.data.z$Shore.z,HDist.data.z$Struct.z)
colnames(hdist.pars)<-c("Shore.z","Struct.z")

##Biotic meadow characteristics model
##Run full model
Reg.comp.all<-capscale(Oreg.fish.Z~Bed.form+Bed.distribution+Tide.range.code+Epi.z+Shoot.z+Biom.z,Env.data.z,add=TRUE)
anova(Reg.comp.all,step=999,perm.max=999) #should be significant before running selection process
RsquareAdj(Reg.comp.all)$adj.r.squared #this r-squared should not be exceeded

##Check for collinearity - avoid VIF>4
#Remove variable with highest VIF if >4, re-run and check, repeat if needed
vif.cca(Reg.comp.all)

##Forward selection process
step.forward<-ordistep(capscale(Oreg.fish.Z~1,Env.data.z,add=TRUE),
                       scope=formula(Reg.comp.all),direction="forward",pstep=999)

##Check r-squared of each model
RsquareAdj(capscale(Oreg.fish.Z~Tide.range.code + Epi.z ,
                    Env.data.z,add=TRUE))$adj.r.squared #check if exceeds global r-sq

##Run best-fit model
bio.pars<-capscale(Oreg.fish.Z~ Tide.range.code + Epi.z + Shoot.z + Biom.z,Env.data.z,add=TRUE)
anova(bio.pars,step=999,perm.max=999)

##Make dataframe with selected variables
meadow.pars<-cbind(Env.data.z$Tide.range.code, Env.data.z$Epi.z, Env.data.z$Shoot.z, Env.data.z$Biom.z) 
colnames(meadow.pars)<-c("Tide.range","Epi.z","Shoot.z","Biom.z")

##Physical characteristics model
##Run full model
Reg.comp.all<-capscale(Oreg.fish.Z~Exp.code + SST.July.z + Tide.z + River.z + Type_simple + SST.annual.z +
                         Sal.insitu.z,Phys.data.z,add=TRUE) 
anova(Reg.comp.all,step=999,perm.max=999) #should be significant before running selection process
RsquareAdj(Reg.comp.all)$adj.r.squared #this r-squared should not be exceeded

##Check for collinearity - avoid VIF>4
#Remove variable with highest VIF if >4, re-run and check, repeat if needed
vif.cca(Reg.comp.all)

##Forward selection process
step.forward<-ordistep(capscale(Oreg.fish.Z~1,Phys.data.z,add=TRUE),
                       scope=formula(Reg.comp.all),direction="forward",pstep=999)

##Check r-squared of each model
RsquareAdj(capscale(Oreg.fish.Z~SST.July.z + Tide.z,
                    Phys.data.z,add=TRUE))$adj.r.squared #check if exceeds global r-sq

##Run best-fit model
anova.cca(capscale(Oreg.fish.Z~ SST.July.z + Tide.z + Sal.insitu.z + SST.annual.z, 
                   Phys.data.z,add=TRUE))

##Make dataframe with selected variables
physical.pars<-cbind(Phys.data.z$Sal.insitu.z,Phys.data.z$SST.July.z,Phys.data.z$SST.annual.z,Phys.data.z$Tide.z)
colnames(physical.pars)<-c("Sal.insitu.z" ,"SST.July.z", "SST.annual.z", "Tide.z")

##########Variance partitioning##########
#Uses selected variable components saved from process above
spe.part.all<- varpart(Oreg.fish.Z,hdist.pars,physical.pars,meadow.pars,sqrt.dist=FALSE,add=TRUE)  

#Save explained variance for 'all' data
indfract<-as.data.frame(spe.part.all$part[2])
part<-c("a","b","c","d","e","f","g","resid")
indfract.adjr<-as.data.frame(cbind(part,indfract[,3]))
write.csv(indfract.adjr,"VarParts_ALL_Zscore.csv")

#Save explained variance for each disturbance category
indfract<-as.data.frame(spe.part.all$part[3])
part<-c("physical","both","meadow","resid")
indfract.adjr<-as.data.frame(cbind(part,indfract[,3]))
write.csv(indfract.adjr,"VarParts_HIGH_Zscore.csv")

##########Plot variance explained using stacked bars, Figure 4##########
library(ggplot2)
library(grid)

##Load dataframe that has all variance components combined
Var2<-read.csv("VarParts_Zscore_comb.csv")
attach(Var2)
#Set up dataframe for graphing
Var2$dist<-factor(Var2$dist,levels=c("A","B","C","D"))
Var2$part<-factor(Var2$part,levels=c("Physical","Biotic","Physical & Biotic","Disturbance","Physical & Disturbance", "Physical, Biotic, Disturbance"))
Var2<-Var2[order(Var2$dist,Var2$part),]

##Plot
Var.plot<-ggplot(aes(x=dist,y=perc),data=Var2) +
  geom_bar(aes(fill=part),stat="identity",position="stack",col="black",width=0.6) + 
  scale_y_continuous(limits=c(0,50),breaks=seq(0,50,by=10)) + 
  scale_fill_manual(values=c("gray10","gray40","gray85","gray20","gray50","gray99"),
                    breaks=c("Physical","Biotic","Physical & Biotic","Disturbance","Physical & Disturbance", "Physical, Biotic, Disturbance"),
                    labels=c("P","B","P & B","D","P & D", "P, B, D")) +
  scale_x_discrete(breaks=c("A","B","C","D"),labels=c("All","Low","Medium","High")) +
  theme_bw()+
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(angle=90,size=18), 
        axis.text.y=element_text(size=16), 
        legend.position=c(.01,1),legend.justification=c(.01,1), 
        legend.background = element_rect(fill=alpha("white", 0)),
        legend.text=element_text(size=16),
        legend.title=element_blank(), legend.key=element_blank(),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  xlab("Human disturbance") + ylab("Variance explained (%)") 
Var.plot

#Add epiphyte inset graph (see lmer test below)
vp <- viewport(width = 0.4, height = 0.4, x = unit(12,"lines"),
               y = unit(30, "lines"), just = c("left", "top"))

full <- function() {
  print(Var.plot)
  print(epi.plot, vp = vp)
}
full()

##########LMERs on physical and biotic variables across disturbance categories##########
library(lme4)

##Run for all best-fit physical and biotic variables
epi.t<-lmer(Epi.lia~H.Disturbance + (1|Region) -1,Eel,REML=TRUE, control = lmerControl(optimizer="nloptwrap"))
summary(epi.t)
confint(epi.t)


