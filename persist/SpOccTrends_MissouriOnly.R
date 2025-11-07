library(tidyverse)


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

NET=subset(NET,NET$RepeatID%in%NOGREENsites)

NET=NET%>%pivot_longer(cols = c(28:121),names_to = "Species",values_to = "change")

####FROM HERE, MUST RUN "DatasetPrep.R" for section
NET=subset(NET,!is.na(NET$change))
NET$Species[which(NET$Species=="RMCOT" | NET$Species=="COLCOT")]="MOTCOT"
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
write.csv(NETsub, "table1materialWITHconfidence.csv")

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


ggsave(filename="NativeSpChange_CI_MISSOURIONLY.tiff",dpi = 400, width = 10, height = 8, units = "in")

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

ggsave(filename = "guildNETgraph.tiff", dpi=400, height = 4, width = 15, units = "in")






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
#============================SECTION B: Missouri Basin ONLY ANALYSES & GRAPHS
#==================================================================================
#Manually removed all green river basin data from nssn$obs
nogreen=read.csv("NOGREENmetrics.csv")
library(lme4)
library(cv)
library(MuMIn)


###-----TURNOVER-----
ngFULL=glm(Turnover~scale(temp)*scale(barrier)*scale(length)*(Yrange), family="binomial", data=nogreen)
summary(ngFULL)

ngFULL=glm(Turnover~scale(temp)*scale(barrier)*scale(length)*scale(prosper), family="binomial", data=nogreen,  na.action = na.fail)
summary(ngFULL)
dredge.results=dredge(ngFULL)
dredge.results=subset(dredge.results,dredge.results$delta<=2)

ngnointeractions=glm(Turnover~scale(temp)+scale(barrier)+scale(length)+scale(prosper), family="binomial", data=nogreen,  na.action = na.fail)
summary(ngnointeractions)

ngTEMP=glm(Turnover~temp, data=nogreen, family = "binomial")
summary(ngTEMP)

ngTEMPPROSP=glm(Turnover~scale(temp)*scale(prosper), data=nogreen, family = "binomial")
summary(ngTEMPPROSP)

ngINTERCEPT=glm(Turnover~1, data = nogreen, family = "binomial")
summary(ngINTERCEPT)

cv(ngTEMP, k="loo")
cv(ngTEMPPROSP, k="loo")
cv(ngnointeractions,k="loo")
cv(ngFULL, k="loo")

#VAR IMP
coeff=as.data.frame(ngnointeractions$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractions$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Missouri Basin Community Turnover")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))

####
#BASIN SPECIFIC MODELS
nogreen$PU[which(nogreen$PU=="LTMO")]="CHEY"
ngUPMO=subset(nogreen,nogreen$PU=="UPMO")
ngPLAT=subset(nogreen,nogreen$PU=="PLAT")
ngCHEY=subset(nogreen,nogreen$PU=="CHEY")
ngYELL=subset(nogreen,nogreen$PU=="YELL")







#Upper Missouri
ngnointeractionsUPMO=glm(Turnover~scale(temp)+scale(barrier)+scale(length)+scale(prosper), family="binomial", data=ngUPMO,  na.action = na.fail)
summary(ngnointeractionsUPMO)

#VAR IMP
coeff=as.data.frame(ngnointeractionsUPMO$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractionsUPMO$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Upper Missouri Community Turnover")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))



#Yellowstone
ngnointeractionsYELL=glm(Turnover~scale(temp)+scale(barrier)+scale(length)+scale(prosper), family="binomial", data=ngYELL,  na.action = na.fail)
summary(ngnointeractionsYELL)

#VAR IMP
coeff=as.data.frame(ngnointeractionsYELL$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractionsYELL$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Yellowstone Community Turnover")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))









#Cheyenne-Little Missouri
ngnointeractionsCHEY=glm(Turnover~scale(temp)+scale(barrier)+scale(length)+scale(prosper), family="binomial", data=ngCHEY,  na.action = na.fail)
summary(ngnointeractionsCHEY)

#VAR IMP
coeff=as.data.frame(ngnointeractionsCHEY$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractionsCHEY$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Cheyenne-Little Missouri Community Turnover")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))







#Platte
ngnointeractionsPLAT=glm(Turnover~scale(temp)+scale(barrier)+scale(length)+scale(prosper), family="binomial", data=ngPLAT,  na.action = na.fail)
summary(ngnointeractionsPLAT)

#VAR IMP
coeff=as.data.frame(ngnointeractionsPLAT$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractionsPLAT$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Platte Community Turnover")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))






nogreen%>%group_by(PU)%>%
  summarise(avgT=mean(Turnover))




###-----NATIVE PERSISTENCE-----
ngFULL=glm(pNat~scale(temp)*scale(barrier)*scale(length)*(Yrange), family="binomial", data=nogreen)
summary(ngFULL)

nogreenNATIVE=nogreen%>%filter(!is.na(pNat))
ngFULL=glm(pNat~scale(temp)*scale(barrier)*scale(length)*scale(prosper)*scale(pisc), family="binomial", data=nogreenNATIVE,  na.action = na.fail)
summary(ngFULL)
dredge.results=dredge(ngFULL)
dredge.results=subset(dredge.results,dredge.results$delta<=2)

ngnointeractions=glm(pNat~scale(temp)+scale(barrier)+scale(length)+scale(prosper)+scale(pisc), family="binomial", data=nogreenNATIVE,  na.action = na.fail)
summary(ngnointeractions)

ngTEMP=glm(pNat~temp, data=nogreenNATIVE, family = "binomial")
summary(ngTEMP)

ngTEMPPROSP=glm(pNat~scale(temp)*scale(prosper), data=nogreenNATIVE, family = "binomial")
summary(ngTEMPPROSP)

ngINTERCEPT=glm(pNat~1, data = nogreen, family = "binomial")
summary(ngINTERCEPT)

cv(ngTEMP, k="loo")
cv(ngTEMPPROSP, k="loo")
cv(ngnointeractions,k="loo")
cv(ngFULL, k="loo")

#VAR IMP
coeff=as.data.frame(ngnointeractions$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractions$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length*"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Native Species Persistence")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))

NOGREENsites=nogreen$RepetID






ngUPMO=subset(nogreenNATIVE,nogreenNATIVE$PU=="UPMO")
ngPLAT=subset(nogreenNATIVE,nogreenNATIVE$PU=="PLAT")
ngCHEY=subset(nogreenNATIVE,nogreenNATIVE$PU=="CHEY")
ngYELL=subset(nogreenNATIVE,nogreenNATIVE$PU=="YELL")


nguppmo=glm(pNat~scale(temp)*scale(length)*scale(prosper)*scale(pisc), family="binomial", data=ngUPMO,  na.action = na.fail)
summary(nguppmo)
dredge.results.uppmo=dredge(nguppmo)
dredge.results.uppmo=subset(dredge.results.uppmo,dredge.results.uppmo$delta<2)

#####Upper Missouri basin
ngnointeractions=glm(pNat~scale(temp)+scale(barrier)+scale(length)+scale(prosper)+scale(pisc), family="binomial", data=ngUPMO,  na.action = na.fail)
summary(ngnointeractions)


#VAR IMP
coeff=as.data.frame(ngnointeractions$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractions$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Upper Missouri Native Species Persistence")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))








#####Yellowstone basin
ngnointeractions=glm(pNat~scale(temp)+scale(barrier)+scale(length)+scale(prosper)+scale(pisc), family="binomial", data=ngYELL,  na.action = na.fail)
summary(ngnointeractions)


#VAR IMP
coeff=as.data.frame(ngnointeractions$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractions$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Yellowstone Native Species Persistence")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))










#####Little Mo Chey basin
ngnointeractions=glm(pNat~scale(temp)+scale(barrier)+scale(length)+scale(prosper)+scale(pisc), family="binomial", data=ngCHEY,  na.action = na.fail)
summary(ngnointeractions)


#VAR IMP
coeff=as.data.frame(ngnointeractions$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractions$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Cheyenne-Little Missouri Native Species Persistence")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))








#####Platte basin
ngnointeractions=glm(pNat~scale(temp)+scale(barrier)+scale(length)+scale(prosper)+scale(pisc), family="binomial", data=ngPLAT,  na.action = na.fail)
summary(ngnointeractions)


#VAR IMP
coeff=as.data.frame(ngnointeractions$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractions$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Platte Native Species Persistence")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))









NOGREENsites=nogreen$RepetID

###-----NATIVE PERSISTENCE-----
ngFULL=glm(pNat~scale(temp)*scale(barrier)*scale(length)*(Yrange), family="binomial", data=nogreen)
summary(ngFULL)

nogreenNATIVE=nogreen%>%filter(!is.na(pNat))
ngFULL=glm(pNat~scale(temp)*scale(barrier)*scale(length)*scale(prosper)*scale(pisc), family="binomial", data=nogreenNATIVE,  na.action = na.fail)
summary(ngFULL)
dredge.results=dredge(ngFULL)
dredge.results=subset(dredge.results,dredge.results$delta<=2)

ngnointeractions=glm(pNat~scale(temp)+scale(barrier)+scale(length)+scale(prosper)+scale(pisc), family="binomial", data=nogreenNATIVE,  na.action = na.fail)
summary(ngnointeractions)

ngTEMP=glm(pNat~temp, data=nogreenNATIVE, family = "binomial")
summary(ngTEMP)

ngTEMPPROSP=glm(pNat~scale(temp)*scale(prosper), data=nogreenNATIVE, family = "binomial")
summary(ngTEMPPROSP)

ngINTERCEPT=glm(pNat~1, data = nogreen, family = "binomial")
summary(ngINTERCEPT)

cv(ngTEMP, k="loo")
cv(ngTEMPPROSP, k="loo")
cv(ngnointeractions,k="loo")
cv(ngFULL, k="loo")

#VAR IMP
coeff=as.data.frame(ngnointeractions$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractions$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Native Species Persistence")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")




#####-----INDIVIDUAL SPECIES PERSISTENCE MODELS-----
library(dplyr)
library(ggplot2)

species_cols <- NETsub30$Species
predictors_all <- c("temp", "barrier", "length", "prosper", "pisc")

# species that should NOT include "pisc"
no_pisc_species <- c("LL", "GSUN", "SMB", "NP")

models <- setNames(vector("list", length(species_cols)), species_cols)

for (i in species_cols) {
  cat("\n--------------------------------------------\n")
  cat("Fitting:", i, "\n")
  
  # Choose predictors
  predictors <- if (i %in% no_pisc_species) {
    setdiff(predictors_all, "pisc")
  } else {
    predictors_all
  }
  
  # Keep rows with no NAs in the response or predictors
  keep <- complete.cases(nogreen[, c(i, predictors)])
  i.sub <- nogreen[keep, , drop = FALSE]
  
  # Build formula and fit model
  terms_scaled <- paste0("scale(", predictors, ")")
  fml <- reformulate(terms_scaled, response = i)
  i.model <- glm(fml, family = binomial, data = i.sub)
  models[[i]] <- i.model
  
  # Print model summary
  print(summary(i.model))
  
  # --- Variable importance plot ---
  coeff <- as.data.frame(i.model$coefficients)
  coeff$Variable <- rownames(coeff)
  names(coeff)[1] <- "Score"
  coeff <- coeff[-1, , drop = FALSE]  # remove intercept
  coeff$Score2 <- abs(coeff$Score)
  coeff$Sign <- ifelse(coeff$Score > 0, "Pos", "Neg")
  
  coeff$Var2 <- NA
  coeff$Var2[coeff$Variable == "scale(temp)"]    <- "Temperature"
  coeff$Var2[coeff$Variable == "scale(barrier)"] <- "Barriers"
  coeff$Var2[coeff$Variable == "scale(length)"]  <- "Fragment Length"
  coeff$Var2[coeff$Variable == "scale(prosper)"] <- "Flow Permanence"
  coeff$Var2[coeff$Variable == "scale(pisc)"]    <- "Piscivores"
  
  coeff <- coeff %>%
    arrange(Score2) %>%
    mutate(Var2 = factor(Var2, levels = Var2))
  
  p <- ggplot(coeff, aes(x = Var2, y = Score2, colour = Sign)) +
    geom_segment(aes(xend = Var2, y = 0, yend = Score2)) +
    geom_point(size = 4) +
    theme_light() +
    scale_color_manual(values = c("#ff9999", "#018080")) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(label = paste(i, "- Variable Importance")) +
    xlab("") +
    coord_flip() +
    ylab("Absolute Coefficient")
  
  print(p)
  
  # Optional: save plot
  # ggsave(filename = paste0("varimp_", i, ".png"), plot = p, width = 6, height = 4)
}





# --- Build a species x predictors table (signed version) ---

predictors_all <- c("temp", "barrier", "length", "prosper", "pisc")

term_name <- function(var) paste0("scale(", var, ")")

# choose scaling method: "sum" or "max"
scale_method <- "sum"  # or "max"

make_species_row <- function(sp, fit) {
  sm <- summary(fit)
  coefs <- sm$coefficients
  
  # Pull coefficients and p-values for our predictors
  est <- sapply(predictors_all, function(v) {
    rn <- term_name(v)
    if (rn %in% rownames(coefs)) coefs[rn, "Estimate"] else NA_real_
  })
  pvals <- sapply(predictors_all, function(v) {
    rn <- term_name(v)
    if (rn %in% rownames(coefs)) coefs[rn, "Pr(>|z|)"] else NA_real_
  })
  
  # Scale coefficients by sum or max of absolute values
  denom <- if (scale_method == "sum") sum(abs(est), na.rm = TRUE) else max(abs(est), na.rm = TRUE)
  if (is.finite(denom) && denom != 0) {
    scaled_est <- est / denom
  } else {
    scaled_est <- est
  }
  
  # Format with asterisk for p < 0.1 (but keep sign)
  fmt <- function(val, p) {
    if (is.na(val)) return(NA_character_)
    paste0(format(round(val, 3), nsmall = 3),
           ifelse(!is.na(p) && p < 0.1, "*", ""))
  }
  cells <- mapply(fmt, scaled_est, pvals, USE.NAMES = FALSE)
  names(cells) <- predictors_all
  
  # AIC
  aic_val <- AIC(fit)
  
  # Return one row
  data.frame(
    species = sp,
    as.list(cells),
    AIC = round(aic_val, 2),
    check.names = FALSE
  )
}

# Build the full table
species_vec <- names(models)
rows <- lapply(species_vec, function(sp) make_species_row(sp, models[[sp]]))
varimp_table_signed <- do.call(rbind, rows)
varimp_table_signed <- varimp_table_signed[, c("species", predictors_all, "AIC")]

# Preview
print(varimp_table_signed)

write.csv(varimp_table_signed, "species_var_importance_signed.csv", row.names = FALSE)











#####################################################################################################
#============================SECTION C: ORDINATION
#####################################################################################################
#============================Dataset preparation

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
wc=wc%>%rename("temp"="S1_93_11", "size"="F_MAUG_HIS", "barrier"="DSbarrier","length"="Length_km")


removes=read.csv("REMOVES.csv")
removes=removes[,c(1,3)]
#removes=removes%>%rename("RepetID"="RepeatID")
wc=left_join(wc,removes,by="RepeatID")
wc=subset(wc,wc$GearMiss!="Y")
wc=subset(wc, wc$size<1000)
wc=subset(wc, wc$RepeatID!=234 & wc$RepeatID!=74& wc$RepeatID!=237& wc$RepeatID!=238)
#write.csv(NET, "PersistenceMetrics_filtered.csv")
wc=subset(wc,wc$PU!="GREEN")
wc=subset(wc,wc$RepeatID!=368)
wc=subset(wc,wc$RepeatID!=431)
wc=subset(wc,wc$RepeatID!=358)
wc=subset(wc,wc$RepeatID!=356)
wc=subset(wc,wc$RepeatID!=397)
wc=subset(wc,wc$RepeatID!=329)

wc$Site_ID=NA
wc$Site_ID=paste(wc$RepeatID, wc$TIME, sep = "_")
wc.comm=wc[,c(139,42:135)]
colSums(wc.comm[2:95])
wc.comm <- wc.comm[, c(1, which(colSums(wc.comm[, 2:95]) > 0) + 1)]

wc.meta=wc[,c(29,30,137,26,36,39,41)]


#####################################################
#---------------------NMDS--------------------------
####################################################
library(vegan)



rowSums(wc.comm[,-1])
# Make first column row names and convert to matrix
wc.mat <- as.matrix(wc.comm[,-1])
rownames(wc.mat) <- wc.comm[,1]
# Example:
nmds <- metaMDS(wc.mat, distance = "jaccard", binary=T, trymax = 20,autotransform = F)

sc <- as.data.frame(scores(nmds, display = "sites"))
sc=cbind(sc, wc.meta)
sc$NMDS1=as.numeric(sc$NMDS1)
sc$NMDS2=as.numeric(sc$NMDS2)
# Ensure TIME and PU are proper factors
sc$TIME <- factor(sc$TIME, levels = c("EARLY", "LATE"))

sc$temp_cat <- NA
sc$temp_cat[sc$temp >= 20] <- "WARM"
sc$temp_cat[sc$temp <  20] <- "COOL"
sc$temp_cat <- factor(sc$temp_cat, levels = c("COOL","WARM"))  # nice legend order

# Colors for temperature categories
temp_cols <- setNames(c("steelblue", "tomato"), levels(sc$temp_cat))

# Base empty plot (no points)
plot(sc$NMDS1, sc$NMDS2, type = "n", xlab = "NMDS1", ylab = "NMDS2")

# Arrows per site, colored by temp category at that site; ensure EARLY -> LATE
by(sc, sc$RepeatID, function(df) {
  if (nrow(df) == 2) {
    df <- df[order(df$TIME), ]  # EARLY (row 1) -> LATE (row 2)
    col_here <- temp_cols[as.character(df$temp_cat[1])]  # same site temp_cat
    arrows(df$NMDS1[1], df$NMDS2[1],
           df$NMDS1[2], df$NMDS2[2],
           length = 0.06, angle = 20, lwd = 2,
           col = col_here)
  }
})

# Legend for temperature categories
legend("bottomleft", lwd = 2, col = temp_cols,
       legend = names(temp_cols), title = "Temperature", bty = "n")





library(dplyr)

# Make sure TIME is ordered
sc$TIME <- factor(sc$TIME, levels = c("EARLY", "LATE"))

# Calculate arrow lengths (distance between EARLY and LATE)
arrow_lengths <- sc %>%
  arrange(RepeatID, TIME) %>%
  group_by(RepeatID) %>%
  filter(n() == 2) %>%                    # only include sites with both samples
  summarize(
    PU        = first(PU),
    temp_cat  = first(temp_cat),
    length    = sqrt((NMDS1[2] - NMDS1[1])^2 + (NMDS2[2] - NMDS2[1])^2)
  )

# Inspect results
head(arrow_lengths)

# Average arrow length by temperature category
arrow_lengths %>%
  group_by(temp_cat) %>%
  summarize(
    mean_length = mean(length, na.rm = TRUE),
    sd_length   = sd(length, na.rm = TRUE),
    se_length = sd(length)/sqrt(length(length)),
    n           = n()
  )




#BIOTIC HOMOGENIZATION TEST
sc$PU   <- factor(sc$PU)

# 2. Build your community matrix and distance matrix
# (Assuming wc.mat or comm is your sites × species presence/absence matrix)
dist_jac <- vegdist(wc.mat, method = "jaccard", binary = TRUE)

# 3. Combine metadata (must match row order of wc.mat)
env <- sc[, c("RepeatID", "PU", "TIME", "temp_cat")]

# 4. PERMDISP setup
# Test whether within-basin dispersions differ between TIME periods
# (First, group by PU × TIME)
group <- interaction(env$PU, env$TIME)

disp <- betadisper(dist_jac, group = group)

# 5. Test for differences in dispersion between time periods within each basin
anova(disp)              # overall test, p-value <0.01
perm_res=permutest(disp, pairwise = TRUE)  # permutation test (pairwise comparisons)

perm_res$pairwise$permuted


# 6. Extract mean distances to centroid for plotting
centroid_df <- as.data.frame(disp$group)
centroid_df$distance <- disp$distances
centroid_df <- cbind(centroid_df, env)

aggregate(distance ~ PU + TIME, data = centroid_df, FUN = mean)





#BIOTIC HOMOGENIZATION TEST WITH TEMPERATURE
sc$PU   <- factor(sc$PU)
sc$temp_cat = factor(sc$temp_cat)
# 2. Build your community matrix and distance matrix
# (Assuming wc.mat or comm is your sites × species presence/absence matrix)
dist_jac <- vegdist(wc.mat, method = "jaccard", binary = TRUE)

# 3. Combine metadata (must match row order of wc.mat)
env <- sc[, c("RepeatID", "PU", "TIME", "temp_cat")]

# 4. PERMDISP setup
# Test whether within-basin dispersions differ between TIME periods
# (First, group by PU × TIME)
group <- interaction(env$PU, env$temp_cat,env$TIME)

disp <- betadisper(dist_jac, group = group)

# 5. Test for differences in dispersion between time periods within each basin
anova(disp)              # overall test, p-value <0.01
perm_res=permutest(disp, pairwise = TRUE)  # permutation test (pairwise comparisons)

perm_res$pairwise$permuted


# 6. Extract mean distances to centroid for plotting
centroid_df <- as.data.frame(disp$group)
centroid_df$distance <- disp$distances
centroid_df <- cbind(centroid_df, env)

aggregate(distance ~ PU + temp_cat+TIME, data = centroid_df, FUN = mean)



