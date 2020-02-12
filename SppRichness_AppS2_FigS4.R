################################################################
#Script for Iacarella et al. ............                      #
#                                                              #
#Sample based species richness rarefaction curves              #
################################################################

library(iNEXT)
library(ggplot2)

##Sample-size based rarefaction using incidence based data

##Load dataframe of species counts for all sites (sites are rows, species are columns)
Reg.fish<-read.csv("Region.spp.comb.csv")
attach(Reg.fish)

#Remove uneeded columns and sort by disturbance category
Reg.simple<-Reg.fish[,c(-1,-2)]
Reg.simple<-Reg.simple[order(Reg.simple$Dist),]

##Subset into disturbance categories, transpose for analysis (spp are rows, sites are columns), 
#and make prescence/absence matrix
High<-subset(Reg.simple,Dist == "high")
High<-High[,-1]
High.t<-t(High)
colnames(High.t) = High.t[1, ] # the first row will be the header
High.t = High.t[-1, ]
High.n <- apply(High.t,2,as.numeric)
High.pres <- ifelse(High.n > 0, 1, 0) #make matrix pres/abs
High.pres <- apply(High.pres,2,as.integer)
rownames(High.pres)<-rownames(High.t)

Med<-subset(Reg.simple,Dist == "medium")
Med<-Med[,-1]
Med.t<-t(Med)
colnames(Med.t) = Med.t[1, ] # the first row will be the header
Med.t = Med.t[-1, ]
Med.n <- apply(Med.t,2,as.numeric)
Med.pres <- ifelse(Med.n > 0, 1, 0) #make matrix pres/abs
Med.pres <- apply(Med.pres,2,as.integer)
rownames(Med.pres)<-rownames(Med.t)

Low<-subset(Reg.simple,Dist == "low")
Low<-Low[,-1]
Low.t<-t(Low)
colnames(Low.t) = Low.t[1, ] # the first row will be the header
Low.t = Low.t[-1, ]
Low.n <- apply(Low.t,2,as.numeric)
Low.pres <- ifelse(Low.n > 0, 1, 0) #make matrix pres/abs
Low.pres <- apply(Low.pres,2,as.integer)
rownames(Low.pres)<-rownames(Low.t)

##Make list of matrices (one for each disturbance category)
my.list<-sapply(list(set1,set2), identity, simplify="array")
my.list[[1]] <- Low.pres  
my.list[[2]] <- Med.pres  
my.list[[3]] <- High.pres  

my.list<-list(Low.pres,Med.pres,High.pres)
names(my.list)<-c("Low","Med","High")
str(my.list)

##Run rarefaction on matrix list
Reg.rar<-iNEXT(my.list,datatype = "incidence_raw",endpoint=150,nboot = 999)

ggiNEXT(Reg.rar)

##Save results for graph
df<-fortify(Reg.rar,type=1)
df.point <- df[which(df$method=="observed"),]
df.line <- df[which(df$method!="observed"),]
df.line$method <- factor(df.line$method, 
                         c("interpolated", "extrapolated"),
                         c("interpolation", "extrapolation"))

##Plot Appendix S2, FigS4
gg.rich<- ggplot(df, aes(x=x,y=y,shape=site)) + 
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr), fill="gray68", alpha=0.7) +
  geom_point(aes(shape=site), size=5, data=df.point, colour="black") +
  scale_shape_manual(values=c(15,16,17),breaks=c("Low","Med","High"),labels=c("Low","Medium","High"), name="Anthropogenic disturbance") +
  geom_line(aes(linetype=method), lwd=1.5, data=df.line, colour="black") +
  scale_linetype_manual(values=c("solid","dashed"),breaks=c("interpolation","extrapolation"),labels=c("Interpolation","Extrapolation"),name="Method",guide=FALSE) +
  labs(x="Number of sites", y="Species richness") +
  theme_bw() + 
  scale_y_continuous(limits=c(0,80),breaks=seq(0,80,by=20)) +
  scale_x_continuous(limits=c(0,150),breaks=seq(0,150,by=50)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(angle=90,size=18), 
        axis.text.y=element_text(size=16), 
        legend.position=c(.8,0.2),legend.justification=c(.8,.2),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16), legend.key=element_blank(),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank())
gg.rich