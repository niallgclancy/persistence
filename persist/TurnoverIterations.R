###################################################
########################################
######Creation of Persistence Metrics Dataset from Repeated Sites Data
library(tidyverse)
library(sf)
##Load repeated site data
wc=read.csv("S2DR_v_1_1_WILDCARD.csv")

wc=wc[-240,]#remove sample that was added twice
removes=read.csv("REMOVES.csv")
removes=removes[,c(1,3)]

wc=left_join(wc,removes,by="RepeatID")
wc=subset(wc,wc$GearMiss!="Y")
wc=subset(wc, wc$F_MAUG_HIS<1000)
wc=subset(wc, wc$RepeatID!=234 & wc$RepeatID!=74& wc$RepeatID!=237& wc$RepeatID!=238)


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


#First Resample Iteration of Turnover calculation
EARLofTed=sample_n(wc_early,286)
LORDofBill=sample_n(wc_late,286)

EARLofTed.sum = EARLofTed%>%
  pivot_longer(42:135, names_to = "Species")%>%filter(value==1)%>%
  group_by(Species)%>%
  summarise(n=sum(value))
EARLofTed.sum$n[which(EARLofTed.sum$n>0)]=1  

LORDofBill.sum = LORDofBill%>%
  pivot_longer(42:135, names_to = "Species")%>%filter(value==1)%>%
  group_by(Species)%>%
  summarise(n=sum(value))
LORDofBill.sum$n[which(LORDofBill.sum$n>0)]=1  

Duke=full_join(EARLofTed.sum,LORDofBill.sum,by="Species")
Duke$n.x[which(is.na(Duke$n.x))]=0
Duke$n.y[which(is.na(Duke$n.y))]=0
Duke$dif=NA
Duke$dif=Duke$n.x-Duke$n.y
Duke$dif=abs(Duke$dif)
Duke=Duke%>%
  summarise(turn=sum(dif)/length(dif))
Duke$Grab=NA
Duke$Grab=0



#Loop iterations
iterations=seq(from=1, to=999, by=1)
for (i in iterations) {
  

EARLofTed=sample_n(wc_early,286)
LORDofBill=sample_n(wc_late,286)

EARLofTed.sum = EARLofTed%>%
  pivot_longer(42:135, names_to = "Species")%>%filter(value==1)%>%
  group_by(Species)%>%
  summarise(n=sum(value))
EARLofTed.sum$n[which(EARLofTed.sum$n>0)]=1  

LORDofBill.sum = LORDofBill%>%
  pivot_longer(42:135, names_to = "Species")%>%filter(value==1)%>%
  group_by(Species)%>%
  summarise(n=sum(value))
LORDofBill.sum$n[which(LORDofBill.sum$n>0)]=1  

Duke2=full_join(EARLofTed.sum,LORDofBill.sum,by="Species")
Duke2$n.x[which(is.na(Duke2$n.x))]=0
Duke2$n.y[which(is.na(Duke2$n.y))]=0
Duke2$dif=NA
Duke2$dif=Duke2$n.x-Duke2$n.y
Duke2$dif=abs(Duke2$dif)
Duke2=Duke2%>%
  summarise(turn=sum(dif)/length(dif))
Duke2$Grab=NA
Duke2$Grab=i

print(i)
Duke=rbind(Duke,Duke2)

}

mean(Duke$turn)#15.5%
