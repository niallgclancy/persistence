###################################################
########################################
######Creation of Persistence Metrics Dataset from Repeated Sites Data
library(tidyverse)
library(sf)
##Load repeated site data
wc=read.csv("S2DR_v_1_1_WILDCARD.csv")

wc=wc[-240,]#remove sample that was added twice

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

######################################################################
############Now incorporate information on native/introduced status and glacial relicts
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


######################Glacial Relicts Analysis
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

##############LATE colonists for comparison
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





#--------------Split into Colonization and Persistence Datasets------------------------
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



#--------------Create Refugia (persistence or colonization =1) Extirpation Dataset------------------------
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
