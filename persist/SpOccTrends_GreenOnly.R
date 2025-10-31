#============================SECTION A: Net Species Changes to Native and Introduced species (if native anywhere in study area, restricted to native range)
#==================================================================================
#---------------------Net Species Changes to Only Native Populations
NET=read.csv("PersistenceMetrics.csv")
nogreen=read.csv("NOGREENmetrics.csv")
NOGREENsites=nogreen$RepetID


removes=read.csv("REMOVES.csv")
removes=removes[,c(1,3)]
#removes=removes%>%rename("RepetID"="RepeatID")
NET=left_join(NET,removes,by="RepeatID")
NET=subset(NET,NET$GearMiss!="Y")
NET=subset(NET, NET$F_MAUG_<1000)
NET=subset(NET, NET$RepeatID!=234 & NET$RepeatID!=74& NET$RepeatID!=237& NET$RepeatID!=238)
'%!in%' <- function(x,y)!('%in%'(x,y))
NET=subset(NET,NET$RepeatID%!in%NOGREENsites)

NET=NET%>%pivot_longer(cols = c(28:121),names_to = "Species",values_to = "change")

NET=subset(NET,!is.na(NET$change))
#NET$Species[which(NET$Species=="RMCOT" | NET$Species=="COLCOT")]="MOTCOT"
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
NET=left_join(NET,traits,by="Species")
NET$Status[which(NET$Species=="RDSH")]="Native"
NET$Status[which(NET$Species=="RDSH" & NET$RepeatID %in% rdshnn)]="Introduced"
NET$Status[which(NET$Species=="LKCH" & NET$RepeatID %in% lkchnn)]="Introduced"
NET$Status[which(NET$Species=="FHMN" & NET$RepeatID %in% FHMNnn)]="Introduced"
NET$Status[which(NET$Species=="LNDC" & NET$RepeatID %in% LNDCnn)]="Introduced"
NET$Status[which(NET$Species=="CRCH" & NET$RepeatID %in% CRCHnn)]="Introduced"
NET$Status[which(NET$Species=="BLBH")]="Native"
NET$Status[which(NET$Species=="BLBH" & NET$RepeatID %in% BLBHnn)]="Introduced"
NET$Status[which(NET$Species=="CCAT" & NET$RepeatID %in% CCATnn)]="Introduced"
NET$Status[which(NET$Species=="LING" & NET$RepeatID %in% LINGnn)]="Introduced"
NET$Status[which(NET$Species=="WSU" & NET$RepeatID %in% WSUnn)]="Introduced"
NET$Status[which(NET$Species=="RMCT" & NET$RepeatID %in% RMCTnn)]="Introduced"
NET$Status[which(NET$Species=="DRUM")]="Native"
NET$Status[which(NET$Species=="DRUM" & NET$RepeatID %in% drumnn)]="Introduced"
NET$Status[which(NET$Species=="PKF")]="Native"
NET$Status[which(NET$Species=="PKF" & NET$RepeatID %in% pkfnn)]="Introduced"
NET$Status[which(NET$Species=="PTMN")]="Native"
NET$Status[which(NET$Species=="PTMN" & NET$RepeatID %in% ptmnnn)]="Introduced"
NET$Status[which(NET$Species=="BRSB")]="Native"
NET$Status[which(NET$Species=="BRSB" & NET$RepeatID %in% brsbnn)]="Introduced"
invasive20s=c("LL","CARP","GSUN","RB","RSSH","EB","SMB","NP")
NET=subset(NET,NET$Status=="Native"|NET$Species%in%invasive20s)
nativeSites=NET%>%
  group_by(Species)%>%
  summarise(sites=length(RepeatID))
#write.csv(nativeSites,"nativesites.csv")

NET2=NET%>%group_by(Species)%>%summarise(net=mean(change), std.E = sd(change)/sqrt(length(change)), n=length(change))
NET2$CI90 = NA
NET2$CI90 = NET2$std.E*1.645
NETsub = subset(NET2, NET2$n>=20)
NETsub30 = subset(NET2, NET2$n>=50)


#---------------COMBINED GUILDS------------------------
traits$ThermalFlow=NA
traits$ThermalFlow=paste(traits$SimpleThermal,traits$LowFlowSens,sep = "")

traits$ThermalLifHis=NA
traits$ThermalLifHis=paste(traits$SimpleThermal,traits$LifeHist,sep = "")

traits$FlowLifHis=NA
traits$FlowLifHis=paste(traits$LowFlowSens,traits$LifeHist,sep = "")

#-----------------------------------------------------


traits2=traits[,c(1,3)]
NET=left_join(NETsub,traits,by="Species")
NETsub=subset(NET,NET$Species!="ONC" & NET$Species!="FMxWSU")
mean(NETsub$net)
#write.csv(NETsub, "table1materialWITHconfidence.csv")

NETsub$CommonName[which(NETsub$CommonName=="Northern Redbelly Dace")]="Nor. Redbelly Dace"
NETsub$CommonName[which(NETsub$CommonName=="Rocky Mountain Cutthroat Trout")]="R.M. Cutthroat Trout"


NETsub$Direc=NA
NETsub$Direc[which(NETsub$net>=1)]="Pos"
NETsub$Direc[which(NETsub$net<1)]="Neg"
NETsub$Status=NA
NETsub$Status="Native"
NETsub$Status[which(NETsub$Species%in%invasive20s)]="Introduced"

library(ggbreak)
library(ggh4x)
changenative=NETsub%>%
  ggplot(aes(x=reorder(CommonName,-net),y=net, shape = Status, color=Status))+
  geom_hline(yintercept = 0, color=alpha("black", alpha = 0.3), linewidth=1.1)+
  geom_segment( aes(x=CommonName, xend=CommonName, y=net-CI90, yend=net+CI90), color="black") +
  geom_point(size=3)+
  theme_light()+
  ylab(label = "Site Occupancy Trend")+
  xlab(label="")+
  scale_color_manual(values = c("#ff9999","#005555"))+
  scale_shape_manual(values=c(15,17))+
  scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1, by=0.2))+
  coord_flip()+
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face = 2),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_text(size=12, color = "black"),
    axis.text.x = element_text(size=10, color = "black"),
    legend.position = "none")+
  theme(ggh4x.axis.ticks.length.minor=rel(1))
changenative


#ggsave(filename="NativeSpChange_CI.tiff",dpi = 400, width = 10, height = 8, units = "in")

NETsub$LifeHist[which(NETsub$LifeHist=="Equil")]="Equilibrium"
NETsub$LifeHist[which(NETsub$LifeHist=="Opp")]="Opportunistic"
NETsub$LifeHist[which(NETsub$LifeHist=="Per")]="Periodic"

graphA=NETsub%>%
  ggplot(aes(x = Status, y= net))+
  geom_hline(yintercept = 0, color=alpha("black", alpha = 0.3), linewidth=1.1)+
  geom_boxplot()+
  theme_light()+
  ylab("Site Occupancy Trend")+
  xlab("Native vs. Introduced")+
  ggtitle("A. Native v. Introduced")+
  annotate("text", x = 2.3, y = -0.45, label = "p = 0.22")+
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    axis.text.x = element_text(size=12, color = "black"),
    axis.title.x = element_blank(),
    legend.position = "none")+
  theme(ggh4x.axis.ticks.length.minor=rel(1))


graphB=NETsub%>%filter(SimpleThermal=="Cold"|SimpleThermal=="Cool"|SimpleThermal=="Warm")%>%
  ggplot(aes(x = SimpleThermal, y= net))+
  geom_hline(yintercept = 0, color=alpha("black", alpha = 0.3), linewidth=1.1)+
  geom_boxplot()+
  theme_light()+
  ylab("Site Occupancy Trend")+
  xlab("Thermal Guild")+
  ggtitle("B. Thermal Guild")+
  annotate("text", x = 3, y = -0.48, label = "p = 0.18")+
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    axis.text.x = element_text(size=12, color = "black"),
    axis.title = element_blank(),
    legend.position = "none")+
  theme(ggh4x.axis.ticks.length.minor=rel(1))

graphC=NETsub%>%filter(LowFlowSens=="Low"|LowFlowSens=="Moderate"|LowFlowSens=="High")%>%
  ggplot(aes(x = LowFlowSens, y= net))+
  geom_hline(yintercept = 0, color=alpha("black", alpha = 0.3), linewidth=1.1)+
  geom_boxplot()+
  theme_light()+
  ylab("Site Occupancy Trend")+
  xlab("Low Flow Sensitivity Guild")+
  ggtitle("C. Flow Preference")+
  annotate("text", x = 2.3, y = -0.45, label = "p = 0.06")+
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    axis.text.x = element_text(size=12, color = "black"),
    axis.title = element_blank(),
    legend.position = "none")+
  theme(ggh4x.axis.ticks.length.minor=rel(1))

graphD=NETsub%>%
  ggplot(aes(x = LifeHist, y= net))+
  geom_hline(yintercept = 0, color=alpha("black", alpha = 0.3), linewidth=1.1)+
  geom_boxplot()+
  theme_light()+
  ylab("Site Occupancy Trend")+
  ggtitle("D. Life History Strategy")+
  annotate("text", x = 3, y = -0.45, label = "p = 0.01")+
  annotate("text", x = 1.1, y = 0.4, label = "a")+
  annotate("text", x = 2.1, y = 0.15, label = "b")+
  annotate("text", x = 3.2, y = 0.3, label = "ab")+
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    axis.text.x = element_text(size=12, color = "black"),
    axis.title = element_blank(),
    legend.position = "none")+
  theme(ggh4x.axis.ticks.length.minor=rel(1))
graphD


library(ggpubr)


guild.net.graph=ggarrange(graphA, graphB, graphC, graphD,
                          nrow = 1)
guild.net.graph

#ggsave(filename = "guildNETgraph.tiff", dpi=400, height = 4, width = 15, units = "in")






#---------------ANOVAs & t.tests---------------
#NATIVE v INTRODUCED
NETsubnative=subset(NETsub, NETsub$Status=="Native")
NETsubintroduced=subset(NETsub, NETsub$Status!="Native")
t.test(NETsubnative$net,NETsubintroduced$net)
wilcox.test(NETsubnative$net,NETsubintroduced$net)
median(NETsubnative$net)
median(NETsubintroduced$net)

#THERMAL
NETsub2=NETsub%>%filter(SimpleThermal=="Cold"|SimpleThermal=="Cool"|SimpleThermal=="Warm")
therm.aov=aov(net~SimpleThermal, data=NETsub2)
summary(therm.aov) #p=0.248
#plot(therm.aov)
#kruskal.test(net~SimpleThermal, data=NETsub)
NETsub2%>%group_by(SimpleThermal)%>%summarise(mean(net), median(net))

#Flow
NETsubLow=subset(NETsub, NETsub$LowFlowSens=="Low")
NETsubMod=subset(NETsub, NETsub$LowFlowSens=="Moderate")
t.test(NETsubLow$net,NETsubMod$net)
mean(NETsubLow$net)
mean(NETsubMod$net)

#LifeHist
lifhist.aov=aov(net~LifeHist, data=NETsub)
summary(lifhist.aov) #p=0.014
TukeyHSD(lifhist.aov)
#plot(lifhist.aov)
NETsub%>%group_by(LifeHist)%>%summarise(mean(net), median(net))








#===================================================================================
#============================SECTION B: Green Basin ONLY ANALYSES & GRAPHS
#==================================================================================
#####################################################################################################
library(SSNbler)
library(SSN2)
library(caret)
library(tidyverse)
library(sf)
library(spmodel)


#===========================================================
#============================SECTION 1: Dataset preparation
#===========================================================
######Creation of Persistence Metrics Dataset from Repeated Sites Data
##Load repeated site data
wc=read.csv("S2DR_v_1_1_WILDCARD.csv")

wc=wc[-240,]#remove sample that was added twice
#wcincluded=subset(wc, wc$RepeatID%in%included)
wc_nn=wc
wc_nn=wc_nn[,c(29,137)]
wc_nn=wc_nn%>%group_by(RepeatID)%>%summarise(pu=first(PU))
wc_nn=wc_nn%>%rename("PU"="pu")
write.csv(wc_nn,file = "processingunits.csv")

############Correct for samples that differentially classified Hybognathus sp.
wchybogchecks=subset(wc, wc$HYBOG==1 | wc$WS_PLMN==1)
wchybogchecks2=as.list(wchybogchecks$RepeatID)
wchybogchecks=subset(wc,wc$RepeatID%in%wchybogchecks2)
wchybogchecks=wchybogchecks[,c(29,30,74:77,134,136)]
wc$PLMN[which(wc$RepeatID==23 & wc$TIME=="LATE")]=0
wc$WS_PLMN[which(wc$RepeatID==23 & wc$TIME=="LATE")]=1
wc$PLMN[which(wc$RepeatID==167 & wc$TIME=="LATE")]=0
wc$WS_PLMN[which(wc$RepeatID==167 & wc$TIME=="LATE")]=1
wc$PLMN[which(wc$RepeatID==266 & wc$TIME=="LATE")]=0
wc$WS_PLMN[which(wc$RepeatID==266 & wc$TIME=="LATE")]=1
wc$PLMN[which(wc$RepeatID==341 & wc$TIME=="EARLY")]=0
wc$WSMN[which(wc$RepeatID==341 & wc$TIME=="EARLY")]=0
wc$HYBOG[which(wc$RepeatID==341 & wc$TIME=="EARLY")]=1
wc$numSp[which(wc$RepeatID==341 & wc$TIME=="EARLY")]=6
wc$PLMN[which(wc$RepeatID==387 & wc$TIME=="EARLY")]=0
wc$HYBOG[which(wc$RepeatID==387 & wc$TIME=="EARLY")]=1


wc_early=wc%>%filter(TIME=="EARLY")
wc_late=wc%>%filter(TIME=="LATE")

wc_meta=wc_early[,c(2,8:10,13:20,22,26:29,33:37,39,41)]

wc_earlylong=wc_early%>%pivot_longer(cols = c(42:135), names_to = "Species", values_to = "Present")
wc_earlylong=wc_earlylong[,c(29,42,44,45)]
wc_earlylong=wc_earlylong%>%rename("numSpE"="numSp","PresentE"="Present")

wc_latelong=wc_late%>%pivot_longer(cols = c(42:135), names_to = "Species", values_to = "Present")
wc_latelong=wc_latelong[,c(29,42,44,45)]
wc_latelong=wc_latelong%>%rename("numSpL"="numSp","PresentL"="Present")

wc_long=left_join(wc_earlylong,wc_latelong,by=c("RepeatID","Species"))

wc_long$Change=NA
wc_long$richChange=NA
wc_long$Change[which(wc_long$PresentE==1 & wc_long$PresentL==1)]=0
wc_long$Change[which(wc_long$PresentE==1 & wc_long$PresentL==0)]=-1
wc_long$Change[which(wc_long$PresentE==0 & wc_long$PresentL==0)]=NA
wc_long$Change[which(wc_long$PresentE==0 & wc_long$PresentL==1)]=1
wc_long$richChange=wc_long$numSpL-wc_long$numSpE
#TURNOVER METRICS
wc_long2=wc_long
wcT=wc_long2%>%filter(Change==1 | Change==-1)%>%group_by(RepeatID)%>%summarise(numchan=length(Species))
wcA=wc_long2%>%filter(!is.na(Change))%>%group_by(RepeatID)%>%summarise(numall=length(Species))
wcT=left_join(wcA,wcT,by="RepeatID")
wcT$numchan[which(is.na(wcT$numchan))]=0
wcT$Turnover=NA
wcT$Turnover=wcT$numchan/wcT$numall
wcT$numall=NULL
wcT$numchan=NULL
write.csv(wcT,file = "TurnoverMetric.csv")

#wc_long$wc_NULL#wc_long$wc_long#wc_long$numSpE=NULL
wc_long$numSpL=NULL
wc_long$PresentE=NULL
wc_long$PresentL=NULL
wcc=wc_long%>%pivot_wider(names_from = "Species", values_from = "Change")

wc=left_join(wc_meta,wcc,by="RepeatID")
wc_rumham=wc

#Now incorporate information on native/introduced status and postglacial pioneers
traits=read.csv("traits.csv")
traits=traits[,c(3,7,14)]
traits=traits%>%rename("Species"="Code")
wctraits=left_join(wc_long,traits,by="Species")
###############For species that are both native and introduced to some parts of study area
#############list of repeat ID's where they are nonnative
brsbnn=c(2,29,32,35,37,38,112,157,158,172,185,191,219,220,246,283,342,361,383,415,433,460)
ptmnnn=c(46,95,295,298,327,424,425,426)
pkfnn=c(22,32,36,37,92,93,94,95,96,97,165,225,226,229,230,263,264,330,346,349,352,358,378,381,425,426,429,430,433,434,436)
drumnn=c(185,184)
rdshnn=c(60,61,63,64)

wctraits=left_join(wctraits,wc_nn,by="RepeatID")
lkchnn=subset(wctraits,wctraits$Species=="LKCH" & wctraits$PU=="GREEN")
lkchnn=lkchnn$RepeatID
LNDCnn=subset(wctraits,wctraits$Species=="LNDC" & wctraits$PU=="GREEN")
LNDCnn=LNDCnn$RepeatID
CRCHnn=subset(wctraits,wctraits$Species=="CRCH" & wctraits$PU=="GREEN")
CRCHnn=CRCHnn$RepeatID
CCATnn=subset(wctraits,wctraits$Species=="CCAT" & wctraits$PU=="GREEN")
CCATnn=CCATnn$RepeatID
LINGnn=subset(wctraits,wctraits$Species=="LING" & wctraits$PU=="GREEN")
LINGnn=LINGnn$RepeatID
WSUnn=subset(wctraits,wctraits$Species=="WSU" & wctraits$PU=="GREEN")
WSUnn=WSUnn$RepeatID
FHMNnn=subset(wctraits,wctraits$Species=="FHMN" & wctraits$PU=="GREEN")
FHMNnn=FHMNnn$RepeatID
BLBHnn=subset(wctraits,wctraits$Species=="BLBH")
BLBHnn=subset(BLBHnn,BLBHnn$PU=="UPMO" | BLBHnn$PU=="YELL"|BLBHnn$PU=="LTMO"|BLBHnn$PU=="GREEN")
BLBHnn=BLBHnn$RepeatID
RMCTnn=subset(wctraits,wctraits$Species=="RMCT")
RMCTnn=subset(RMCTnn,RMCTnn$PU=="CHEY" | RMCTnn$PU=="PLAT")
RMCTnn=RMCTnn$RepeatID




wctraits$Status[which(wctraits$Status=="Mixed")]="Native"
wctraits$Status[which(wctraits$Species=="BRSB" & wctraits$RepeatID %in% brsbnn)]="Introduced"
wctraits$Status[which(wctraits$Species=="PTMN" & wctraits$RepeatID %in% ptmnnn)]="Introduced"
wctraits$Status[which(wctraits$Species=="PKF" & wctraits$RepeatID %in% pkfnn)]="Introduced"
wctraits$Status[which(wctraits$Species=="DRUM" & wctraits$RepeatID %in% drumnn)]="Introduced"
wctraits$Status[which(wctraits$Species=="RDSH" & wctraits$RepeatID %in% rdshnn)]="Introduced"

wctraits$Status[which(wctraits$Species=="LKCH" & wctraits$RepeatID %in% lkchnn)]="Introduced"
wctraits$Status[which(wctraits$Species=="FHMN" & wctraits$RepeatID %in% FHMNnn)]="Introduced"
wctraits$Status[which(wctraits$Species=="LNDC" & wctraits$RepeatID %in% LNDCnn)]="Introduced"
wctraits$Status[which(wctraits$Species=="CRCH" & wctraits$RepeatID %in% CRCHnn)]="Introduced"
wctraits$Status[which(wctraits$Species=="BLBH" & wctraits$RepeatID %in% BLBHnn)]="Introduced"
wctraits$Status[which(wctraits$Species=="CCAT" & wctraits$RepeatID %in% CCATnn)]="Introduced"
wctraits$Status[which(wctraits$Species=="LING" & wctraits$RepeatID %in% LINGnn)]="Introduced"
wctraits$Status[which(wctraits$Species=="WSU" & wctraits$RepeatID %in% WSUnn)]="Introduced"
wctraits$Status[which(wctraits$Species=="RMCT" & wctraits$RepeatID %in% RMCTnn)]="Introduced"



wctraits$Glacial[which(wctraits$Status=="Introduced")]=NA


wctraits$richChange=NULL
wctraits$Species=NULL
wctraits=subset(wctraits,wctraits$Status=="Native")

################Native Species Net Change
wcNATIVE=wctraits
wcNATIVE$Glacial=NULL
wcNATIVE=subset(wcNATIVE,!is.na(wcNATIVE$Change))
wcNATIVE=wcNATIVE%>%group_by(RepeatID)%>%summarise(NativeNet=sum(Change))
wc=left_join(wc,wcNATIVE,by="RepeatID")



################Native Species Persistence
wcNATIVE=wctraits
wcNATIVE$Glacial=NULL
wcNATIVE=subset(wcNATIVE,!is.na(wcNATIVE$Change))
wcNATIVE=subset(wcNATIVE,wcNATIVE$Change!=1)#remove colonizations
wcNATIVE2=wcNATIVE%>%group_by(RepeatID)%>%summarise(NativeE=length(numSpE))
wcNATIVE3=wcNATIVE%>%filter(Change==-1)%>%group_by(RepeatID)%>%summarise(lost=length(numSpE))
wcNATIVE=left_join(wcNATIVE2,wcNATIVE3,by="RepeatID")
wcNATIVE$lost[which(is.na(wcNATIVE$lost))]=0
wcNATIVE$propLost=NA
wcNATIVE$propLost=wcNATIVE$lost/wcNATIVE$NativeE
wcNATIVE$persNative=NA
wcNATIVE$persNative=1-wcNATIVE$propLost
wcNATIVE$NativeE=NULL
wcNATIVE$lost=NULL
wcNATIVE$propLost=NULL
wc=left_join(wc,wcNATIVE,by="RepeatID")



##############Community Persistence
wc_long2=subset(wc_long,wc_long$Change!=1)
wc_long2$Species=NULL
wc_long2=wc_long2%>%group_by(RepeatID)%>%
  summarise(numsp=mean(numSpE),change=sum(Change))
wc_long2$oppopersComm=wc_long2$change/wc_long2$numsp
wc_long2$persComm=NA
wc_long2$persComm=1+wc_long2$oppopersComm

wc_long2=wc_long2[,c(1,5)]
wc=left_join(wc,wc_long2,by="RepeatID")


######################Postglacial Pioneer Analysis
wcGLACIAL=wctraits
wcGLACIAL$Status=NULL
wcGLACIAL=subset(wcGLACIAL,wcGLACIAL$Glacial=="L" | wcGLACIAL$Glacial=="E")
wcGLACIAL2=subset(wcGLACIAL,wcGLACIAL$Change==-1 | wcGLACIAL$Change==0)
wcGLACIAL3=wcGLACIAL2%>%filter(Glacial=="E")%>%group_by(RepeatID)%>%summarise(n=length(Change), numlost=-1*sum(Change))
wcGLACIAL3$propLost=NA
wcGLACIAL3$propLost=wcGLACIAL3$numlost/wcGLACIAL3$n
wcGLACIAL3$persGlacE=NA
wcGLACIAL3$persGlacE=1-wcGLACIAL3$propLost
wcGLACIAL3=wcGLACIAL3[,c(1,5)]
wc=left_join(wc,wcGLACIAL3,by="RepeatID")

##############Non Pioneers for comparison
wcGLACIAL=wctraits
wcGLACIAL$Status=NULL
wcGLACIAL=subset(wcGLACIAL,wcGLACIAL$Glacial=="L" | wcGLACIAL$Glacial=="E")
wcGLACIAL2=subset(wcGLACIAL,wcGLACIAL$Change==-1 | wcGLACIAL$Change==0)
wcGLACIAL3=wcGLACIAL2%>%filter(Glacial=="L")%>%group_by(RepeatID)%>%summarise(n=length(Change), numlost=-1*sum(Change))
wcGLACIAL3$propLost=NA
wcGLACIAL3$propLost=wcGLACIAL3$numlost/wcGLACIAL3$n
wcGLACIAL3$persGlacL=NA
wcGLACIAL3$persGlacL=1-wcGLACIAL3$propLost
wcGLACIAL3=wcGLACIAL3[,c(1,5)]
wc=left_join(wc,wcGLACIAL3,by="RepeatID")

#############Predatory Fish Covariate
preds=c("NP","LL","GSUN","RKB","SMB","LMB", "TGT","SPK","BLCR","WHCR","WBS","WE")
wcL=wc_latelong
wcL=subset(wcL,wcL$PresentL==1)
wcL=subset(wcL,wcL$Species%in%preds)
wcL=wcL%>%group_by(RepeatID)%>%summarise(nnPreds=length(PresentL))
wcL$nnPreds[which(wcL$nnPreds>=1)]=1
wc=left_join(wc,wcL,by="RepeatID")
wc$nnPreds[which(is.na(wc$nnPreds))]=0


#############Number Colonizing Species Covariate
col=wc_long
col=subset(col,col$Change==1)
col=col%>%group_by(RepeatID)%>%summarise(nColonizSp=sum(Change))
wc=left_join(wc,col,by="RepeatID")
wc$nColonizSp[which(is.na(wc$nColonizSp))]=0

write.csv(wc,file = "PersistenceMetrics.csv")
#st_write(wc, "PersistenceMetrics.shp")





####################Split into Colonization and Persistence Datasets
#Persistence
wc2=wc%>%pivot_longer(cols=c(27:120), names_to = "Species", values_to = "V")
wc2$V[which(wc2$V==1)]=NA#REmove colonizations
wc2$V[which(wc2$V==0)]=1#fish persistence assigned value 1
wc2$V[which(wc2$V==-1)]=0#fish extirpation assigned value 0
wc2=subset(wc2,!is.na(wc2$V))
wc2=wc2%>%pivot_wider(names_from = "Species", values_from = "V")
write.csv(wc2, file = "PersistenceSubset.csv")
snap=read_sf("WILDCARD_SNAPPED.shp")#pull CRS from shapefile
wgs84=st_crs(snap)
wc2sf <- st_as_sf(wc2, coords = c("LONGITUDE", "LATITUDE"), 
                  crs = wgs84)#convert to sf object with WGS84 projection
#Test that object plots correctly
ggplot() +
  geom_sf(data = wc2sf) +
  coord_sf(datum = st_crs(wc2sf))
#Write SHAPEFILE
st_write(wc2sf,
         "Persistence_Extirpation.shp", driver = "ESRI Shapefile",append = F)



#Colonization
wc3=wc%>%pivot_longer(cols=c(27:120), names_to = "Species", values_to = "V")
wc3$V[which(wc3$V==-1)]=NA#REmove extirpations
wc3=subset(wc3,!is.na(wc3$V))
wc3=wc3%>%pivot_wider(names_from = "Species", values_from = "V")
write.csv(wc3, file = "ColonizationSubset.csv")
wc3sf <- st_as_sf(wc3, coords = c("LONGITUDE", "LATITUDE"), 
                  crs = wgs84)#convert colonizations to sf object with WGS84 projection
#Test that object plots correctly
ggplot() +
  geom_sf(data = wc3sf) +
  coord_sf(datum = st_crs(wc3sf))
#Write SHAPEFILE
st_write(wc3sf,
         "Colonization.shp", driver = "ESRI Shapefile",append = F)



##############Create Refugia (persistence or colonization =1) Extirpation Dataset
wc4=wc%>%pivot_longer(cols=c(27:120), names_to = "Species", values_to = "V")
wc4$V[which(wc4$V==0)]=1#fish persistence & colonization assigned value 1
wc4$V[which(wc4$V==-1)]=0#fish extirpation assigned value 0
wc4=subset(wc4,!is.na(wc4$V))
wc4=wc4%>%pivot_wider(names_from = "Species", values_from = "V")
write.csv(wc4, file = "RefugiaSubset.csv")
snap=read_sf("WILDCARD_SNAPPED.shp")#pull CRS from shapefile
wgs84=st_crs(snap)
wc4sf <- st_as_sf(wc4, coords = c("LONGITUDE", "LATITUDE"), 
                  crs = wgs84)#convert to sf object with WGS84 projection
#Test that object plots correctly
ggplot() +
  geom_sf(data = wc4sf) +
  coord_sf(datum = st_crs(wc4sf))
#Write SHAPEFILE
st_write(wc4sf,
         "RefugiaSites.shp", driver = "ESRI Shapefile",append = F)












#======================================================
#============================SECTION 2: SSN Preparation
#======================================================

#------------------Importing Required Data---------------------

## import the streams, observation sites
streams <- st_read("NSI_fix6.shp")
obs <- st_read("newmetrics.shp")
removes=read.csv("REMOVES.csv")
removes=removes[,c(1,3)]
removes=removes%>%rename("RepetID"="RepeatID")
obs=left_join(obs,removes,by="RepetID")
obs=subset(obs,obs$GearMiss!="Y")
obs=subset(obs, obs$F_MAUG_<1000)
obs=subset(obs, obs$RepetID!=234 & obs$RepetID!=74& obs$RepetID!=237& obs$RepetID!=238)

streams=st_cast(streams, to="LINESTRING")

streams <- st_transform(streams, crs = 5070)
obs <- st_transform(obs, crs = 5070)


## Plot the data using ggplot2
ggplot() +
  geom_sf(data = streams) +
  geom_sf(data = obs, color = "blue", size = 2) +
  coord_sf(datum = st_crs(streams))

# --------------- Build the LSN -------------------------------
## Set the lsn.path variable
lsn.path <- "persist_LSN/"

## Build the LSN
edges <- lines_to_lsn(
  streams = streams,
  lsn_path = lsn.path,
  check_topology = TRUE,
  snap_tolerance = 0.05,
  topo_tolerance = 20,
  overwrite = TRUE
)


# ----------- Incorporate sites into LSN -----------------------
## Incorporate observations. Pay attention to output messages in R
obs <- sites_to_lsn(
  sites = obs,
  edges = edges,
  lsn_path = lsn.path,
  file_name = "obs",
  snap_tolerance = 200,
  save_local = TRUE,
  overwrite = TRUE
)


# --------- Calculate upstream distance for edges ------------
edges$GNIS_ID=as.numeric(edges$GNIS_ID)

edges <- updist_edges(
  edges = edges,
  save_local = TRUE,
  lsn_path = lsn.path,
  calc_length = TRUE
)


# ---------- Calculate upstream distance for sites -----------
site.list <- updist_sites(
  sites = list(
    obs = obs),
  edges = edges,
  length_col = "Length",
  save_local = TRUE,
  lsn_path = lsn.path
)


## Plot the upstream distances for edges and obs
ggplot() +
  geom_sf(data = edges, aes(color = upDist)) +
  geom_sf(data = site.list$obs, aes(color = upDist)) +
  coord_sf(datum = st_crs(streams)) +
  scale_color_viridis_c()

# ------- Calculate AFV for edges --------------------------------
## Summarize h2oAreaKm2 and check for zeros
summary(edges$TotDASqKM) 

edges <- afv_edges(
  edges = edges,
  infl_col = "TotDASqKM",
  segpi_col = "areaPI",
  afv_col = "afvArea",
  lsn_path = lsn.path
)


# ------- Calculate AFV for sites -------------------------------
site.list <- afv_sites(
  sites = site.list,    ## Input is a named list
  edges = edges,
  afv_col = "afvArea",
  save_local = TRUE,
  lsn_path = lsn.path
)


# ---------- Assemble the SSN Object -----------------------------
nssn <- ssn_assemble(
  edges = edges,
  lsn_path = lsn.path,
  obs_sites = site.list$obs,
  ssn_path = "nssn.ssn",
  import = TRUE,
  check = TRUE,
  afv_col = "afvArea",
  overwrite = TRUE
)




## Look at the observations in the SSN object. Notice the new columns
## added when the SSN was assembled
class(nssn$obs)
names(nssn$obs)

## Plot ssn
ggplot() +
  geom_sf(
    data = nssn$edges,
    color = "medium blue",
    aes(linewidth = log10(TotDASqKM))
  ) +
  scale_linewidth(range = c(0.1, 2)) +
  geom_sf(
    data = nssn$obs,
    size = 1.7
  )




# ---- Create Distance Matrices ----------------------------------
## Generate hydrologic distance matrices
ssn_create_distmat(nssn)

## Output a list of distance matrices for obs
obs.distmat<- ssn_get_stream_distmat(nssn,
                                     name = "obs")
names(obs.distmat)
colnames(obs.distmat[[1]])

## Create symmetric hydrologic distance matrix
obs.distmat2 <- obs.distmat[[1]] + t(obs.distmat[[1]])



view(nssn$obs)