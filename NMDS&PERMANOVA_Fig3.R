#######################################################
#Script for Iacarella et al. ............             #
#                                                     #
#NMDS, Figure 3, and PERMANOVA code                   #
#######################################################

##Load dataframes with ALL BC beta-deviation and RC beta-diversity values (Oregonian region only)
#(i.e. without removal of across disturbance category comparisons or within/across region comparisons, as used by bootstrapping test)
#Dataframes are matrices of site-site comparisons
Oreg.z<-read.csv("OregZscore_ALL.csv")
Dist<-Oreg.z$H.Disturbance #Distubance category values for each site
Oreg.z<-Oreg.z[,c(-1,-2)]

#Set up matrix of beta-diversity values with sites as row names
Null.fish2<-as.matrix(Oreg.z)
row.names(Null.fish2) <- Null.fish2[,1]
Null.fish2 <- Null.fish2[,-1]
class(Null.fish2)<-"numeric"
Oreg.fish.Z<-as.dist(Null.fish2)

#Do the same set up above for both BC and RC values

library(car)
Dist.ABC<-recode(Dist,'"high"="c";"medium"="b"; "low"="a"')

##NMDS
library(vegan)
#Run for BC values
spe.nmds.bc<- metaMDS(Oreg.fish.Z,labels=Dist.ABC,autotransform = FALSE)
spe.nmds.bc
#Run for RC values
spe.nmds<- metaMDS(Oreg.fish.RC,labels=Dist.ABC,autotransform = FALSE,trymax=100)
spe.nmds
stressplot(spe.nmds)

##Plot NMDS
NMDS.bc<-data.frame(NMDS1.bc=spe.nmds.bc$points[,1],NMDS2.bc=spe.nmds.bc$points[,2],group=Dist.ABC)
NMDS.mean.bc<-aggregate(. ~ group,NMDS.bc,mean)
ord<-ordiellipse(spe.nmds.bc,Dist.ABC,display="sites",kind="se",conf=0.95,label=T)

NMDS<-data.frame(NMDS1=spe.nmds$points[,1],NMDS2=spe.nmds$points[,2],group=Dist.ABC)
NMDS.mean<-aggregate(. ~ group,NMDS,mean)
ord.rc<-ordiellipse(spe.nmds,Dist.ABC,display="sites",kind="se",conf=0.95,label=T)

##Run ellipse functions
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell.bc <- data.frame()
for(g in levels(NMDS.bc$group)){
  df_ell.bc <- rbind(df_ell.bc, cbind(as.data.frame(with(NMDS.bc[NMDS.bc$group==g,],
                                                         veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale))),group=g))
}

df_ell.rc <- data.frame()
for(g in levels(NMDS$group)){
  df_ell.rc <- rbind(df_ell.rc, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                         veganCovEllipse(ord.rc[[g]]$cov,ord.rc[[g]]$center,ord.rc[[g]]$scale))),group=g))
}

#plot with ellipses of 95% CI
library(ggplot2)
library(gridExtra)

z.nmds.plot<- 
  ggplot(NMDS.bc, aes(NMDS1.bc, NMDS2.bc)) +
  geom_point(stat="identity",aes(shape=group,fill=group),size=2,color="black") +
  scale_shape_manual(values=c(22,21,24),breaks=c("a","b","c"),labels=c("Low","Medium","High"), name="Human disturbance") +
  scale_fill_manual(values=c("white","gray47","gray21"),breaks=c("a","b","c"),labels=c("Low","Medium","High"), name="Human disturbance") +
  scale_y_continuous(limits=c(-10,20),breaks=seq(-10,20,by=5),name="NMDS2, Bray-Curtis Beta-deviation") +
  scale_x_continuous(limits=c(-15,20),breaks=seq(-15,20,by=5),name="NMDS1, Bray-Curtis Beta-deviation") + 
  theme_bw()+
  theme(axis.text.x=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(angle=90,size=8), 
        axis.text.y=element_text(size=8), 
        legend.position=c(.01,0.99),legend.justification=c(.1,.95), 
        legend.background = element_rect(fill="white", size=.2, linetype="solid",color="black"),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8), legend.key=element_blank(),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  annotate("text",x=15,y=19,label="stress = 0.20",size=3)
z.nmds.plot

rc.nmds.plot<- 
  ggplot(NMDS, aes(NMDS1, NMDS2)) +
  geom_point(stat="identity",aes(shape=group,fill=group),size=2,color="black", show.legend = FALSE) +
  scale_shape_manual(values=c(22,21,24),breaks=c("a","b","c"),labels=c("Low","Medium","High"), name="Human disturbance") +
  scale_fill_manual(values=c("white","gray47","gray21"),breaks=c("a","b","c"),labels=c("Low","Medium","High"), name="Human disturbance") +
  scale_y_continuous(limits=c(-.6,.6),breaks=seq(-.6,.6,by=0.3),name="NMDS2, Raup-Crick Beta-diversity") +
  scale_x_continuous(limits=c(-.6,.6),breaks=seq(-.6,.6,by=0.3),name="NMDS1, Raup-Crick Beta-diversity") + 
  theme_bw()+
  theme(axis.text.x=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(angle=90,size=8), 
        axis.text.y=element_text(size=8), 
        panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  annotate("text",x=0.45,y=0.55,label="stress = 0.24",size=3)
rc.nmds.plot 

#Set plots next to each other, Figure 3
grid.arrange(z.nmds.plot,rc.nmds.plot,ncol=2)

##PERMANOVA
Oreg.anova<-adonis(Oreg.fish.Z ~ Dist)
Oreg.anova



