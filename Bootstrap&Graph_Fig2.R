#####################################################
#Script for Iacarella et al. ............           #
#                                                   #
#Bootstrapping and graphing code for Figure 2       #
#####################################################

##Load dataframes representing beta-diversity metric within or across regions by disturbance category
#Same steps taken for beta-diversity metrics and Gower distances of physical and biotic variables

#E.g.
Oreg_col<-read.csv("Oreg_BetabyImpact_withinRegions_Zscore.csv")
attach(Oreg_col)
Aleut_col<-read.csv("Aleut_BetabyImpact_withinRegions_Zscore.csv")
attach(Aleut_col)

##Bind dataframes across provinces
Reg_across<-rbind(Oreg_col,Aleut_col) 
High.dist<-subset(Reg_across,Dist.comp %in% "high")
#Extract beta-diversity measure for bootstrapping
High.dist<-High.dist[,5] 
Med.dist<-subset(Reg_across,Dist.comp %in% "medium")
Med.dist<-Med.dist[,5]
Low.dist<-subset(Reg_across,Dist.comp %in% "low")
Low.dist<-Low.dist[,5]

##Bootstrap mean and confidence interval for each disturbance category
library(boot)
boot.mean <- function(x,i){boot.mean <- mean(x[i])}

#95% CIs using BCa method (bias corrected & accelerated)
H <- boot(High.dist, boot.mean, R = 2000)
conf.H<-boot.ci(H,type="bca") 

M <- boot(Med.dist, boot.mean, R = 2000)
conf.M<-boot.ci(M,type="bca")

L <- boot(Low.dist, boot.mean, R = 2000)
conf.L<-boot.ci(L,type="bca") 

##Assemble boot mean values
mean.H<-H$t0
mean.M<-M$t0
mean.L<-L$t0
Vals.mean<-rbind(mean.H,mean.M,mean.L)

##Assemble upper confidence interval values
ciU.H<-conf.H$bca[5]
ciU.M<-conf.M$bca[5]
ciU.L<-conf.L$bca[5]

##Assemble lower confidence interval values
ciL.H<-conf.H$bca[4]
ciL.M<-conf.M$bca[4]
ciL.L<-conf.L$bca[4]

##Bind disturbance category values
Vals.ciU<-rbind(ciU.H,ciU.M,ciU.L)
Vals.ciL<-rbind(ciL.H,ciL.M,ciL.L)
Dist<-c("C","B","A") #Label disturbance cateogories to fix order

#Bind all values
Vals.all<-cbind(Vals.mean,Vals.ciL,Vals.ciU)
Vals.all<-as.data.frame(Vals.all)
colnames(Vals.all)<-c("mean","low.ci","upp.ci")
Vals.all<-cbind(Dist,Vals.all)

##Save values for later graphing
#E.g.
write.csv(Vals.all,"Vals.all.RCclassic.within.csv")
########

library(ggplot2)
library(gridExtra)

##load all graph files for combined plot (Figure 2)
#E.g.
RC.within<-read.csv("Vals.all.RCclassic.within.csv")
RC.across<-read.csv("Vals.all.RCclassic.across.csv")

##Define and save each plot seperately
plot.P.W<-ggplot(aes(x=Dist,y=mean),data=Phys.within) +
  geom_point(stat="identity",aes(shape=factor(Dist)),size=4,fill="gray47",color="black",show.legend = FALSE) +
  scale_shape_manual(values=c(22,21,24)) +
  geom_errorbar(stat="identity",aes(ymin=low.ci, ymax=upp.ci),width=0.1,colour="black") +
  scale_y_continuous(limits=c(0,.5),breaks=seq(0,.5,by=.1)) +
  scale_x_discrete(breaks=c("A","B","C"),labels=c("Low","Medium","High")) +
  theme_bw()+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(angle=90,size=16), plot.title=element_text(size=16),
        axis.text.y=element_text(size=14),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  xlab("") + ylab("Gower distance") + ggtitle("") +
  annotate(geom="text",x=c(1,2,3),y=c(Phys.within[3,5]+0.05,Phys.within[2,5]+0.05,Phys.within[1,5]+0.05),label=c("a","a","a"),size=5) 
plot.P.W

plot.P.A<-ggplot(aes(x=Dist,y=mean),data=Phys.across) +
  geom_point(stat="identity",aes(shape=factor(Dist)),size=4,fill="gray47",color="black",show.legend = FALSE) +
  scale_shape_manual(values=c(22,21,24)) +
  geom_errorbar(stat="identity",aes(ymin=low.ci, ymax=upp.ci),width=0.1,colour="black") +
  scale_y_continuous(limits=c(0,.5),breaks=seq(0,.5,by=.1)) +
  scale_x_discrete(breaks=c("A","B","C"),labels=c("Low","Medium","High")) +
  theme_bw()+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        plot.title=element_text(size=16),
        axis.text.y=element_blank(), 
        panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  xlab("") + ylab("") + ggtitle("") +
  annotate(geom="text",x=c(1,2,3),y=c(Phys.across[3,5]+.05,Phys.across[2,5]+.05,Phys.across[1,5]+.05),label=c("a","b","c"),size=5)
plot.P.A

plot.B.W<-ggplot(aes(x=Dist,y=mean),data=Biot.within) +
  geom_point(stat="identity",aes(shape=factor(Dist)),size=4,fill="gray47",color="black",show.legend = FALSE) +
  scale_shape_manual(values=c(22,21,24)) +
  geom_errorbar(stat="identity",aes(ymin=low.ci, ymax=upp.ci),width=0.1,colour="black") +
  scale_y_continuous(limits=c(0,0.4),breaks=seq(0,0.4,by=.1)) +
  scale_x_discrete(breaks=c("A","B","C"),labels=c("Low","Medium","High")) +
  theme_bw()+
  theme(axis.text.x=element_text(size=14), 
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(angle=90,size=16), 
        axis.text.y=element_text(size=14),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  xlab("Anthropogenic disturbance") + ylab("Gower distance")  +
  annotate(geom="text",x=c(1,2,3),y=c(Biot.within[3,5]+.05,Biot.within[2,5]+.05,Biot.within[1,5]+.05),label=c("a","b","b"),size=5)
plot.B.W

plot.B.A<-ggplot(aes(x=Dist,y=mean),data=Biot.across) +
  geom_point(stat="identity",aes(shape=factor(Dist)),size=4,fill="gray47",color="black",show.legend = FALSE) +
  scale_shape_manual(values=c(22,21,24)) +
  geom_errorbar(stat="identity",aes(ymin=low.ci, ymax=upp.ci),width=0.1,colour="black") +
  scale_y_continuous(limits=c(0,0.4),breaks=seq(0,0.4,by=.1)) +
  scale_x_discrete(breaks=c("A","B","C"),labels=c("Low","Medium","High")) +
  theme_bw()+
  theme(axis.text.x=element_text(size=14), 
        axis.title.x=element_text(size=16),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(), 
        panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  xlab("") + ylab("")  +
  annotate(geom="text",x=c(1,2,3),y=c(Biot.across[3,5]+.05,Biot.across[2,5]+.05,Biot.across[1,5]+.05),label=c("a","b","b"),size=5)
plot.B.A

plot.Z.W<-ggplot(aes(x=Dist,y=mean),data=Z.within) + 
  geom_point(stat="identity",aes(shape=factor(Dist)),size=4,fill="gray47",color="black",show.legend = FALSE) +
  scale_shape_manual(values=c(22,21,24)) +
  geom_errorbar(stat="identity",aes(ymin=low.ci, ymax=upp.ci),width=0.1,colour="black") +
  scale_y_continuous(limits=c(0,25),breaks=seq(0,25,by=5)) +
  scale_x_discrete(breaks=c("A","B","C"),labels=c("Low","Medium","High")) +
  theme_bw()+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(angle=90,size=16), plot.title=element_text(size=16),
        axis.text.y=element_text(size=14),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  xlab("") + ylab("Bray-Curtis B-dev.") + ggtitle("Within region") +
  annotate(geom="text",x=c(1,2,3),y=c(Z.within[3,5]+2,Z.within[2,5]+2,Z.within[1,5]+2),label=c("a","b","b"),size=5) 
plot.Z.W

plot.Z.A<-ggplot(aes(x=Dist,y=mean),data=Z.across) +
  geom_point(stat="identity",aes(shape=factor(Dist)),size=4,fill="gray47",color="black",show.legend = FALSE) +
  scale_shape_manual(values=c(22,21,24)) +
  geom_errorbar(stat="identity",aes(ymin=low.ci, ymax=upp.ci),width=0.1,colour="black") +
  scale_y_continuous(limits=c(0,25),breaks=seq(0,25,by=5)) +
  scale_x_discrete(breaks=c("A","B","C"),labels=c("Low","Medium","High")) +
  theme_bw()+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        plot.title=element_text(size=16),
        axis.text.y=element_blank(), 
        panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  xlab("") + ylab("") + ggtitle("Across region") +
  annotate(geom="text",x=c(1,2,3),y=c(Z.across[3,5]+2,Z.across[2,5]+2,Z.across[1,5]+2),label=c("a","b","ab"),size=5)
plot.Z.A

plot.RC.W<-ggplot(aes(x=Dist,y=mean),data=RC.within) +
  geom_point(stat="identity",aes(shape=factor(Dist)),size=4,fill="gray47",color="black",show.legend = FALSE) +
  scale_shape_manual(values=c(22,21,24)) +
  geom_errorbar(stat="identity",aes(ymin=low.ci, ymax=upp.ci),width=0.1,colour="black") +
  scale_y_continuous(limits=c(0,0.6),breaks=seq(0,0.6,by=.1)) +
  scale_x_discrete(breaks=c("A","B","C"),labels=c("Low","Medium","High")) +
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(angle=90,size=16), 
        axis.text.y=element_text(size=14),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  xlab("") + ylab("Raup-Crick B-div.")  +
  annotate(geom="text",x=c(1,2,3),y=c(RC.within[3,5]+.05,RC.within[2,5]+.05,RC.within[1,5]+.05),label=c("a","b","b"),size=5)
plot.RC.W

plot.RC.A<-ggplot(aes(x=Dist,y=mean),data=RC.across) +
  geom_point(stat="identity",aes(shape=factor(Dist)),size=4,fill="gray47",color="black",show.legend = FALSE) +
  scale_shape_manual(values=c(22,21,24)) +
  geom_errorbar(stat="identity",aes(ymin=low.ci, ymax=upp.ci),width=0.1,colour="black") +
  scale_y_continuous(limits=c(0,0.6),breaks=seq(0,0.6,by=.2)) +
  scale_x_discrete(breaks=c("A","B","C"),labels=c("Low","Medium","High")) +
  theme_bw()+
  theme(axis.text.x=element_blank(), 
        axis.title.x=element_text(size=16),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(), 
        panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  xlab("") + ylab("")  +
  annotate(geom="text",x=c(1,2,3),y=c(RC.across[3,5]+.05,RC.across[2,5]+.05,RC.across[1,5]+.05),label=c("a","a","b"),size=5)
plot.RC.A


library(gtable)
library(grid)

##Combine plots
#Figure 2a
g1 <- ggplotGrob(plot.Z.W)
index <- subset(g1$layout, name == "axis-b") 
names <- g1$layout$name[g1$layout$t<=index$t]
g1 <- gtable_filter(g1, paste(names, sep="", collapse="|"))
# set height of remaining, empty rows to 0
for (i in (index$t+1):length(g1$heights))
{g1$heights[[i]] <- unit(0, "cm")}

#Figure 2c
g2 <- ggplotGrob(plot.RC.W)
index <- subset(g2$layout, name == "panel") 
# need to work with b here instead of t, to prevent deletion of background
names <- g2$layout$name[g2$layout$b>=index$b]
g2 <- gtable_filter(g2, paste(names, sep="", collapse="|"))
# set height of remaining, empty rows to 0
for (i in 1:(index$b-1))
{g2$heights[[i]] <- unit(0, "cm")}

index <- subset(g2$layout, name == "axis-b") 
names <- g2$layout$name[g2$layout$t<=index$t]
g2 <- gtable_filter(g2, paste(names, sep="", collapse="|"))
# set height of remaining, empty rows to 0
for (i in (index$t+1):length(g2$heights))
{g2$heights[[i]] <- unit(0, "cm")}

#Figure 2e
g3 <- ggplotGrob(plot.P.W)
index <- subset(g3$layout, name == "panel") 
# need to work with b here instead of t, to prevent deletion of background
names <- g3$layout$name[g3$layout$b>=index$b]
g3 <- gtable_filter(g3, paste(names, sep="", collapse="|"))
# set height of remaining, empty rows to 0
for (i in 1:(index$b-1))
{g3$heights[[i]] <- unit(0, "cm")}

index <- subset(g3$layout, name == "axis-b") 
names <- g3$layout$name[g3$layout$t<=index$t]
g3 <- gtable_filter(g3, paste(names, sep="", collapse="|"))
# set height of remaining, empty rows to 0
for (i in (index$t+1):length(g3$heights))
{g3$heights[[i]] <- unit(0, "cm")}

#Figure 2g
g4 <- ggplotGrob(plot.B.W)
index <- subset(g4$layout, name == "panel") 
# need to work with b here instead of t, to prevent deletion of background
names <- g4$layout$name[g4$layout$b>=index$b]
g4 <- gtable_filter(g4, paste(names, sep="", collapse="|"))
# set height of remaining, empty rows to 0
for (i in 1:(index$b-1))
{g4$heights[[i]] <- unit(0, "cm")}

#Figure 2b
g1b <- ggplotGrob(plot.Z.A)
index <- subset(g1b$layout, name == "axis-b") 
names <- g1b$layout$name[g1b$layout$t<=index$t]
g1b <- gtable_filter(g1b, paste(names, sep="", collapse="|"))
# set height of remaining, empty rows to 0
for (i in (index$t+1):length(g1b$heights))
{g1b$heights[[i]] <- unit(0, "cm")}

#Figure 2d
g2b <- ggplotGrob(plot.RC.A)
index <- subset(g2b$layout, name == "panel") 
# need to work with b here instead of t, to prevent deletion of background
names <- g2b$layout$name[g2b$layout$b>=index$b]
g2b <- gtable_filter(g2b, paste(names, sep="", collapse="|"))
# set height of remaining, empty rows to 0
for (i in 1:(index$b-1))
{g2b$heights[[i]] <- unit(0, "cm")}

index <- subset(g2b$layout, name == "axis-b") 
names <- g2b$layout$name[g2b$layout$t<=index$t]
g2b <- gtable_filter(g2b, paste(names, sep="", collapse="|"))
# set height of remaining, empty rows to 0
for (i in (index$t+1):length(g2b$heights))
{g2b$heights[[i]] <- unit(0, "cm")}

#Figure 2f
g2c <- ggplotGrob(plot.P.A)
index <- subset(g2c$layout, name == "panel") 
# need to work with b here instead of t, to prevent deletion of background
names <- g2c$layout$name[g2c$layout$b>=index$b]
g2c <- gtable_filter(g2c, paste(names, sep="", collapse="|"))
# set height of remaining, empty rows to 0
for (i in 1:(index$b-1))
{g2c$heights[[i]] <- unit(0, "cm")}

index <- subset(g2c$layout, name == "axis-b") 
names <- g2c$layout$name[g2c$layout$t<=index$t]
g2c <- gtable_filter(g2c, paste(names, sep="", collapse="|"))
# set height of remaining, empty rows to 0
for (i in (index$t+1):length(g2c$heights))
{g2c$heights[[i]] <- unit(0, "cm")}

#Figure 2h
g3b <- ggplotGrob(plot.B.A)
index <- subset(g3b$layout, name == "panel") 
# need to work with b here instead of t, to prevent deletion of background
names <- g3b$layout$name[g3b$layout$b>=index$b]
g3b <- gtable_filter(g3b, paste(names, sep="", collapse="|"))
# set height of remaining, empty rows to 0
for (i in 1:(index$b-1))
{g3b$heights[[i]] <- unit(0, "cm")}

##Bind the plots together
g.main <- rbind(g1, g2, g3, g4, size="last") #left panel
g.main2 <- rbind(g1b, g2b, g2c, g3b, size="last") #right panel
#Combined panels
grid.arrange(g.main,g.main2,ncol=2,widths=c(2,1.8))

plot(g.main)
