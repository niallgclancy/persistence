#####################################################################################################
# R script for Clancy Dissertation Chapter on Empirical Characteristics of Refugia & Species Trends

#Version 1.0 "Train Rocks"
#12-20-2024

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
wcincluded=subset(wc, wc$RepeatID%in%included)
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
#==================================================================================
#============================SECTION 3: SSN Models for Community & Functional Level
#==================================================================================
#Identify predictors with multicolinearity
library(PerformanceAnalytics)
mdata=nssn$obs[,c(16,21,23,24,32)]
mdata=as.data.frame(mdata)
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)
#None of the selected predictors (for data exploration) are highly correlated (max r = 0.28)

#rename variables
nssn$obs=nssn$obs%>%rename("temp"="S1_93_1", "size"="F_MAUG_", "barrier"="DSbarrr","length"="Lngth_k","pisc"="nnPreds")
nssn$obs=nssn$obs%>%rename("pComm"="persCmm", "pNat"="persNtv", "rGlac"="prsSmGl","rPela"="pelagcs","rPeri"="peridcs","rMont"="montane")


#-----------PROSPER VARIABLE INCOPORATION AND IMPUTATION OF MISSING VALUES
nssn$obs$prosper=NULL
nssn$obs$prosper.x=NULL
nssn$obs$prosper.y=NULL
nssn$obs$X=NULL
nssn$obs$X.x=NULL
nssn$obs$X.y=NULL
nssn$obs$Imputed.x=NULL
nssn$obs$Imputed.y=NULL
nssn$obs$Imputed=NULL
nssn$obs$prosperREG=NULL
#View(nssn$obNULL#VNULL#View(nssn$obNULL#View(nssn$obs)
drought=read.csv("PROSPERfromSando.csv")
drought$prosper=as.numeric(drought$prosper)
#Add values to sites with PROSPER value
nssn$obs=left_join(nssn$obs, drought, by="RepetID")

forimpute=as.data.frame(nssn$obs)
forimpute=forimpute[,c(8,19,21,141,142)]
sizeproLM=lm(prosper~log(size)*ELEV, data = forimpute)
sizeproLM
summary(sizeproLM)
forimpute2=subset(forimpute,is.na(forimpute$prosper))
forimpute3=subset(forimpute,forimpute$Imputed=="fromnear")
library(modelr)
forimpute2=forimpute2[,c(1:3)]
forimpute2=as.data.frame(forimpute2)
forimpute2=forimpute2 %>% add_predictions(sizeproLM, var = "prosperREG")
forimpute=forimpute2[,c(2,4)]

nssn$obs=left_join(nssn$obs,forimpute,by="RepetID")
nssn$obs$prosper[which(is.na(nssn$obs$prosper))]=0
nssn$obs$prosperREG[which(is.na(nssn$obs$prosperREG))]=0
nssn$obs$prosper2=NA
nssn$obs$prosper2=nssn$obs$prosper+nssn$obs$prosperREG
nssn$obs$prosper=NULL
nssn$obs$prosperREG=NULL
nssn$obs$Imputed=NULL
nssn$obs=nssn$obs%>%rename(prosper=prosper2)

nssn$obs$size2=NA
nssn$obs$size2=nssn$obs$size*0.02831683199881
nssn$obs$size=NA
nssn$obs$size=nssn$obs$size2
nssn$obs$size2=NULL

#Colinearity check
mdata=nssn$obs[,c(16,21,23,24,32,141)]
mdata=as.data.frame(mdata)
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)




#-----------Full Community Models-------------------------------------------------------------------------

quantile(nssn$obs$Yrange)
nssn$obs$YrangeCat=NA
nssn$obs$YrangeCat[which(nssn$obs$Yrange<16)]="S"
nssn$obs$YrangeCat[which(nssn$obs$Yrange>=16 & nssn$obs$Yrange<22)]="M"
nssn$obs$YrangeCat[which(nssn$obs$Yrange>=22 & nssn$obs$Yrange<30)]="L"
nssn$obs$YrangeCat[which(nssn$obs$Yrange>=30)]="XL"

range(nssn$obs$Yrange)

#full model -- max likelihood
ssn_mod <- ssn_glm(
  formula =   pComm~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc)+scale(prosper),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)
summary(ssn_mod)

#null model -- max likelihood
ssn_null <- ssn_glm(
  formula =   pComm~ 1,
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)

glances(ssn_mod,ssn_null)
loocv(ssn_mod)#0.258
loocv(ssn_null)#0.260

#MODEL FOR VAR IMP
coeff=as.data.frame(ssn_mod$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ssn_mod$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size*"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length*"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores*"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
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
    legend.position = "none",
    plot.title = element_text(hjust = -0.2)
  ) +
  ggtitle(label="Community Persistence")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")
ggsave(filename = "VarImp_pComm.tiff",height = 7, width = 4, units = "in", dpi=400)




#-----------Native Species Persistence Models-------------------------------------------------------------------------
mdata=subset(nssn$obs,!is.na(nssn$obs$pNat))
mdata=mdata[,c(16,21,23,24,32)]
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)
disp_init <- dispersion_initial(family = "beta", dispersion = 1, known = "dispersion")

nssn$obs$pNat[which(nssn$obs$pNat==1)]=0.999
nssn$obs$pNat[which(nssn$obs$pNat==0)]=0.001

#full model -- max likelihood
ssn_modN <- ssn_glm(
  formula =   pNat~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc)+scale(prosper),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)
summary(ssn_modN)
varcomp(ssn_modN)
#significant only -- max likelihood
ssn_modN_sig <- ssn_glm(
  formula =   pNat~ scale(temp)*scale(length),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  nugget_type = "none"
)
summary(ssn_modN_sig)
#significant only -- max likelihood
ssn_modN_sigNI <- ssn_glm(
  formula =   pNat~ scale(temp)+scale(length),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  nugget_type = "none"
)
loocv(ssn_modN_sig)
loocv(ssn_modN_sigNI)

#null model -- max likelihood
ssn_null <- ssn_glm(
  formula =   pNat~ 1,
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)

glances(ssn_modN,ssn_null)
loocv(ssn_modN)#0.287
loocv(ssn_null)#0.302
loocv(ssn_modN_sig)#0.287 (with interaction)

#-------------------Native Sp without "m" gear sites
nssn2=nssn
nssn2$obs$pNat[which(nssn2$obs$GearMiss=="M")]=NA

#full model -- max likelihood
ssn_modN_hiconf <- ssn_glm(
  formula =   pNat~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc)+scale(prosper),
  family = "beta",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)
summary(ssn_modN_hiconf) #order: 
summary(ssn_modN)


#MODEL FOR VAR IMP
coeff=as.data.frame(ssn_modN$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ssn_modN$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature**"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length**"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"

coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"
VarImpNative=coeff%>%arrange(Score2) %>%
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
    axis.title.x = element_blank(),
    axis.text = element_text(size=10, color = "black"),
    legend.position = "none"
  ) +
  ggtitle(label="A. Native Species Persistence")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")
VarImpNative
ggsave(filename = "VarImp_pNative.tiff",height = 4, width = 5, units = "in", dpi=400)


#-----------Postglacial Pioneer Refugia Models-------------------------------------------------------------------------
mdata=subset(nssn$obs,!is.na(nssn$obs$rGlac))
mdata=mdata[,c(16,21,23,24,32)]
mdata=as.data.frame(mdata)
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)

#full model -- max likelihood
disp_init <- dispersion_initial(family = "beta", dispersion = 1, known = "dispersion")


ssn_modG <- ssn_glm(
  formula =   rGlac~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)
summary(ssn_modG)

#null model -- max likelihood
ssn_null <- ssn_glm(
  formula =   rGlac~ 1,
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)

#Significant only
ssn_modG_sig <- ssn_glm(
  formula =   rGlac~ scale(length),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  nugget_type = "none"
)
summary(ssn_modG)


glances(ssn_modG,ssn_null)
loocv(ssn_modG)#0.313
loocv(ssn_null)#0.320
loocv(ssn_modG_sig)#0.315


#VAR IMP
coeff=as.data.frame(ssn_modG$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ssn_modG$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length*"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

VarImpGlacial=coeff%>%arrange(Score2) %>%
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
    legend.position = "none",
    plot.title = element_text(hjust = -0.2)
  ) +
  ggtitle(label="C. Postglacial Pioneer Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")
VarImpGlacial

ggsave(filename = "VarImp_rGlacial.tiff",height = 4, width = 5, units = "in", dpi=400)

#-----------Postglacial Pioneer Persistence-------------------------
#---------------------create persistence metric-----------------------
pers_sub=read.csv("PersistenceSubset.csv")
obs2=pers_sub
relict=obs2%>%pivot_longer(cols=35:118,names_to = "Species", values_to = "change")
####FROM HERE, MUST RUN "DatasetPrep.R" for section
relict=subset(relict,!is.na(relict$change))
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
relict=left_join(relict,traits,by="Species")
relict$Status[which(relict$Species=="RDSH")]="Native"
relict$Status[which(relict$Species=="RDSH" & relict$RepeatID %in% rdshnn)]="Introduced"
relict$Status[which(relict$Species=="LKCH" & relict$RepeatID %in% lkchnn)]="Introduced"
relict$Status[which(relict$Species=="FHMN" & relict$RepeatID %in% FHMNnn)]="Introduced"
relict$Status[which(relict$Species=="LNDC" & relict$RepeatID %in% LNDCnn)]="Introduced"
relict$Status[which(relict$Species=="CRCH" & relict$RepeatID %in% CRCHnn)]="Introduced"
relict$Status[which(relict$Species=="BLBH")]="Native"
relict$Status[which(relict$Species=="BLBH" & relict$RepeatID %in% BLBHnn)]="Introduced"
relict$Status[which(relict$Species=="CCAT" & relict$RepeatID %in% CCATnn)]="Introduced"
relict$Status[which(relict$Species=="LING" & relict$RepeatID %in% LINGnn)]="Introduced"
relict$Status[which(relict$Species=="WSU" & relict$RepeatID %in% WSUnn)]="Introduced"
relict$Status[which(relict$Species=="RMCT" & relict$RepeatID %in% RMCTnn)]="Introduced"
relict$Status[which(relict$Species=="DRUM")]="Native"
relict$Status[which(relict$Species=="DRUM" & relict$RepeatID %in% drumnn)]="Introduced"
relict$Status[which(relict$Species=="PKF")]="Native"
relict$Status[which(relict$Species=="PKF" & relict$RepeatID %in% pkfnn)]="Introduced"
relict$Status[which(relict$Species=="PTMN")]="Native"
relict$Status[which(relict$Species=="PTMN" & relict$RepeatID %in% ptmnnn)]="Introduced"
relict$Status[which(relict$Species=="BRSB")]="Native"
relict$Status[which(relict$Species=="BRSB" & relict$RepeatID %in% brsbnn)]="Introduced"
relict=subset(relict,relict$Status=="Native")
relict=relict[,-c(39:52)]

flowavgs=read.csv("FlowAverages.csv")
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
traits=left_join(traits,flowavgs,by="Species")
smallrelicts=subset(traits, traits$Glacial=="E")
smallrelicts=smallrelicts$Species
relict$X=NULL
obsrelict=relict
obsrelict=subset(obsrelict,obsrelict$Species%in%smallrelicts)
obsrelict=subset(obsrelict,!is.na(obsrelict$change))
obsrelcol=obsrelict%>%group_by(RepeatID)%>%summarise(persSmGlacial=mean(change))
obsrelcol$geometry=NULL
obsrelcol=obsrelcol%>%rename("pGlac"="persSmGlacial")
obsrelcol2=obsrelcol
obsrelcol2=subset(obsrelcol2,!is.na(obsrelcol2$pGlac))

#Add to nssn
obsrelcol=obsrelcol%>%rename("RepetID"="RepeatID")
nssn$obs$pGlac=NULL
nssn$obs=left_join(nssn$obs,obsrelcol,by="RepetID")

PostglacialsXSite=obsrelict%>%
  group_by(RepeatID)%>%
  summarise(numPostglac=length(unique(Species)), pGlac = mean(change), LAT=mean(LATITUDE),LONG=mean(LONGITUDE))
PostNoComm=obsrelict%>%filter(Species!="FHMN"&Species!="WSU")%>%
  group_by(RepeatID)%>%
  summarise(numPostglac=length(unique(Species)), pGlac = mean(change), LAT=mean(LATITUDE),LONG=mean(LONGITUDE))
PostNoComm=PostNoComm%>%filter(RepeatID%in%finalRepeatlist)

PostComm=obsrelict%>%filter(Species=="FHMN"|Species=="WSU")%>%
  group_by(RepeatID)%>%
  summarise(numCommon=length(unique(Species)))

finalRepeatlist=unique(nssn$obs$RepetID)
PostglacialsXSite=PostglacialsXSite%>%filter(RepeatID%in%finalRepeatlist)
post=PostglacialsXSite
post=left_join(post,PostComm,by="RepeatID")
post$pCommon=NA
post$pCommon=post$numCommon/post$numPostglac
post$pCommon[which(is.na(post$pCommon))]=0

post%>%
  ggplot(aes(x=numPostglac, y=pGlac, color=pCommon))+
  geom_jitter()+
  scale_colour_viridis_b()+
  stat_smooth(method = "glm", method.args = list(family=binomial), se=T, linewidth=1.5, color="black")+
  theme_bw()+
  ylab(label = "Postglacial Pioneer Persistence at Site")+
  xlab("Postglacial Pioneer Species Richness at Site")
mean(post$pGlac)

PostNoComm%>%
  ggplot(aes(x=numPostglac, y=pGlac))+
  geom_jitter(color="darkgrey")+
  stat_smooth(method = "glm", method.args = list(family=binomial), se=T, linewidth=1.5, color="black")+
  theme_classic()+
  ggtitle(label = "D. Postglacial Pioneer Persistence by Species Richness")+
  ylab(label = "Postglacial Pioneer Persistence at Site")+
  xlab("Postglacial Pioneer Species Richness at Site")+
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=14))

ggsave(filename = "PPsizeRich.tiff",dpi=400,height = 4, width = 8)



mean(PostNoComm$pGlac)
sd(PostNoComm$pGlac)
length(PostNoComm$pGlac)


post.lm=glm(post$pGlac~post$numPostglac,family = "binomial")
summary(post.lm)
library(DescTools)
PseudoR2(post.lm, which = "Nagelkerke")


write.csv(PostglacialsXSite,"PostglacialsXSite.csv")


#Do same for non pioneers
pers_sub=read.csv("PersistenceSubset.csv")
obs2=pers_sub
relict=obs2%>%pivot_longer(cols=35:118,names_to = "Species", values_to = "change")
####FROM HERE, MUST RUN "DatasetPrep.R" for section
relict=subset(relict,!is.na(relict$change))
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
relict=left_join(relict,traits,by="Species")
relict$Status[which(relict$Species=="RDSH")]="Native"
relict$Status[which(relict$Species=="RDSH" & relict$RepeatID %in% rdshnn)]="Introduced"
relict$Status[which(relict$Species=="LKCH" & relict$RepeatID %in% lkchnn)]="Introduced"
relict$Status[which(relict$Species=="FHMN" & relict$RepeatID %in% FHMNnn)]="Introduced"
relict$Status[which(relict$Species=="LNDC" & relict$RepeatID %in% LNDCnn)]="Introduced"
relict$Status[which(relict$Species=="CRCH" & relict$RepeatID %in% CRCHnn)]="Introduced"
relict$Status[which(relict$Species=="BLBH")]="Native"
relict$Status[which(relict$Species=="BLBH" & relict$RepeatID %in% BLBHnn)]="Introduced"
relict$Status[which(relict$Species=="CCAT" & relict$RepeatID %in% CCATnn)]="Introduced"
relict$Status[which(relict$Species=="LING" & relict$RepeatID %in% LINGnn)]="Introduced"
relict$Status[which(relict$Species=="WSU" & relict$RepeatID %in% WSUnn)]="Introduced"
relict$Status[which(relict$Species=="RMCT" & relict$RepeatID %in% RMCTnn)]="Introduced"
relict$Status[which(relict$Species=="DRUM")]="Native"
relict$Status[which(relict$Species=="DRUM" & relict$RepeatID %in% drumnn)]="Introduced"
relict$Status[which(relict$Species=="PKF")]="Native"
relict$Status[which(relict$Species=="PKF" & relict$RepeatID %in% pkfnn)]="Introduced"
relict$Status[which(relict$Species=="PTMN")]="Native"
relict$Status[which(relict$Species=="PTMN" & relict$RepeatID %in% ptmnnn)]="Introduced"
relict$Status[which(relict$Species=="BRSB")]="Native"
relict$Status[which(relict$Species=="BRSB" & relict$RepeatID %in% brsbnn)]="Introduced"
relict=subset(relict,relict$Status=="Native")
relict=relict[,-c(39:52)]

flowavgs=read.csv("FlowAverages.csv")
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
traits=left_join(traits,flowavgs,by="Species")
smallrelicts=subset(traits, traits$Glacial=="L")
smallrelicts=smallrelicts$Species
relict$X=NULL
obsrelict=relict
obsrelict=subset(obsrelict,obsrelict$Species%in%smallrelicts)
obsrelict=subset(obsrelict,!is.na(obsrelict$change))
obsrelcol=obsrelict%>%group_by(RepeatID)%>%summarise(persSmGlacial=mean(change))
obsrelcol$geometry=NULL
obsrelcol=obsrelcol%>%rename("pNotGlac"="persSmGlacial")
obsrelcol2=obsrelcol
obsrelcol2=subset(obsrelcol2,!is.na(obsrelcol2$pNotGlac))

#Add to nssn
obsrelcol=obsrelcol%>%rename("RepetID"="RepeatID")
nssn$obs$pNotGlac=NULL
nssn$obs=left_join(nssn$obs,obsrelcol,by="RepetID")




#Do same for True Glacial Relicts
pers_sub=read.csv("PersistenceSubset.csv")
obs2=pers_sub
relict=obs2%>%pivot_longer(cols=35:118,names_to = "Species", values_to = "change")
####FROM HERE, MUST RUN "DatasetPrep.R" for section
relict=subset(relict,!is.na(relict$change))
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
relict=left_join(relict,traits,by="Species")
relict$Status[which(relict$Species=="RDSH")]="Native"
relict$Status[which(relict$Species=="RDSH" & relict$RepeatID %in% rdshnn)]="Introduced"
relict$Status[which(relict$Species=="LKCH" & relict$RepeatID %in% lkchnn)]="Introduced"
relict$Status[which(relict$Species=="FHMN" & relict$RepeatID %in% FHMNnn)]="Introduced"
relict$Status[which(relict$Species=="LNDC" & relict$RepeatID %in% LNDCnn)]="Introduced"
relict$Status[which(relict$Species=="CRCH" & relict$RepeatID %in% CRCHnn)]="Introduced"
relict$Status[which(relict$Species=="BLBH")]="Native"
relict$Status[which(relict$Species=="BLBH" & relict$RepeatID %in% BLBHnn)]="Introduced"
relict$Status[which(relict$Species=="CCAT" & relict$RepeatID %in% CCATnn)]="Introduced"
relict$Status[which(relict$Species=="LING" & relict$RepeatID %in% LINGnn)]="Introduced"
relict$Status[which(relict$Species=="WSU" & relict$RepeatID %in% WSUnn)]="Introduced"
relict$Status[which(relict$Species=="RMCT" & relict$RepeatID %in% RMCTnn)]="Introduced"
relict$Status[which(relict$Species=="DRUM")]="Native"
relict$Status[which(relict$Species=="DRUM" & relict$RepeatID %in% drumnn)]="Introduced"
relict$Status[which(relict$Species=="PKF")]="Native"
relict$Status[which(relict$Species=="PKF" & relict$RepeatID %in% pkfnn)]="Introduced"
relict$Status[which(relict$Species=="PTMN")]="Native"
relict$Status[which(relict$Species=="PTMN" & relict$RepeatID %in% ptmnnn)]="Introduced"
relict$Status[which(relict$Species=="BRSB")]="Native"
relict$Status[which(relict$Species=="BRSB" & relict$RepeatID %in% brsbnn)]="Introduced"
relict=subset(relict,relict$Status=="Native")
relict=relict[,-c(39:52)]
mean(relict$change)

flowavgs=read.csv("FlowAverages.csv")
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
traits=left_join(traits,flowavgs,by="Species")
smallrelicts=c("NRBDC","FSDC","LKCH","BRSB","PLSU","GR","NPDC","HHCH","IOWA")
relict$X=NULL
obsrelict=relict
obsrelict=subset(obsrelict,obsrelict$Species%in%smallrelicts)
obsrelict=subset(obsrelict,!is.na(obsrelict$change))
obsrelcol=obsrelict%>%group_by(RepeatID)%>%summarise(persSmGlacial=mean(change))
obsrelcol$geometry=NULL
obsrelcol=obsrelcol%>%rename("pGlacRel"="persSmGlacial")
obsrelcol2=obsrelcol
obsrelcol2=subset(obsrelcol2,!is.na(obsrelcol2$pGlacRel))

#Add to nssn
obsrelcol=obsrelcol%>%rename("RepetID"="RepeatID")
nssn$obs$pGlacRel=NULL
nssn$obs=left_join(nssn$obs,obsrelcol,by="RepetID")



#Do same for Non Relicts
pers_sub=read.csv("PersistenceSubset.csv")
obs2=pers_sub
relict=obs2%>%pivot_longer(cols=35:118,names_to = "Species", values_to = "change")
####FROM HERE, MUST RUN "DatasetPrep.R" for section
relict=subset(relict,!is.na(relict$change))
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
relict=left_join(relict,traits,by="Species")
relict$Status[which(relict$Species=="RDSH")]="Native"
relict$Status[which(relict$Species=="RDSH" & relict$RepeatID %in% rdshnn)]="Introduced"
relict$Status[which(relict$Species=="LKCH" & relict$RepeatID %in% lkchnn)]="Introduced"
relict$Status[which(relict$Species=="FHMN" & relict$RepeatID %in% FHMNnn)]="Introduced"
relict$Status[which(relict$Species=="LNDC" & relict$RepeatID %in% LNDCnn)]="Introduced"
relict$Status[which(relict$Species=="CRCH" & relict$RepeatID %in% CRCHnn)]="Introduced"
relict$Status[which(relict$Species=="BLBH")]="Native"
relict$Status[which(relict$Species=="BLBH" & relict$RepeatID %in% BLBHnn)]="Introduced"
relict$Status[which(relict$Species=="CCAT" & relict$RepeatID %in% CCATnn)]="Introduced"
relict$Status[which(relict$Species=="LING" & relict$RepeatID %in% LINGnn)]="Introduced"
relict$Status[which(relict$Species=="WSU" & relict$RepeatID %in% WSUnn)]="Introduced"
relict$Status[which(relict$Species=="RMCT" & relict$RepeatID %in% RMCTnn)]="Introduced"
relict$Status[which(relict$Species=="DRUM")]="Native"
relict$Status[which(relict$Species=="DRUM" & relict$RepeatID %in% drumnn)]="Introduced"
relict$Status[which(relict$Species=="PKF")]="Native"
relict$Status[which(relict$Species=="PKF" & relict$RepeatID %in% pkfnn)]="Introduced"
relict$Status[which(relict$Species=="PTMN")]="Native"
relict$Status[which(relict$Species=="PTMN" & relict$RepeatID %in% ptmnnn)]="Introduced"
relict$Status[which(relict$Species=="BRSB")]="Native"
relict$Status[which(relict$Species=="BRSB" & relict$RepeatID %in% brsbnn)]="Introduced"
relict=subset(relict,relict$Status=="Native")
relict=relict[,-c(39:52)]
mean(relict$change)

flowavgs=read.csv("FlowAverages.csv")
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
traits=left_join(traits,flowavgs,by="Species")
smallrelicts=c("NRBDC","FSDC","LKCH","BRSB","PLSU","GR","NPDC","HHCH","IOWA")
relict$X=NULL
obsrelict=relict
'%!in%' <- function(x,y)!('%in%'(x,y))
obsrelict=subset(obsrelict,obsrelict$Species%!in%smallrelicts)
obsrelict=subset(obsrelict,!is.na(obsrelict$change))
obsrelcol=obsrelict%>%group_by(RepeatID)%>%summarise(persSmGlacial=mean(change))
obsrelcol$geometry=NULL
obsrelcol=obsrelcol%>%rename("pNotGlacRel"="persSmGlacial")
obsrelcol2=obsrelcol
obsrelcol2=subset(obsrelcol2,!is.na(obsrelcol2$pNotGlacRel))

#Add to nssn
obsrelcol=obsrelcol%>%rename("RepetID"="RepeatID")
nssn$obs$pNotGlacRel=NULL
nssn$obs=left_join(nssn$obs,obsrelcol,by="RepetID")







#Do same for Introduced Species
pers_sub=read.csv("PersistenceSubset.csv")
obs2=pers_sub
relict=obs2%>%pivot_longer(cols=35:118,names_to = "Species", values_to = "change")
####FROM HERE, MUST RUN "DatasetPrep.R" for section
relict=subset(relict,!is.na(relict$change))
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
relict=left_join(relict,traits,by="Species")
relict$Status[which(relict$Species=="RDSH")]="Native"
relict$Status[which(relict$Species=="RDSH" & relict$RepeatID %in% rdshnn)]="Introduced"
relict$Status[which(relict$Species=="LKCH" & relict$RepeatID %in% lkchnn)]="Introduced"
relict$Status[which(relict$Species=="FHMN" & relict$RepeatID %in% FHMNnn)]="Introduced"
relict$Status[which(relict$Species=="LNDC" & relict$RepeatID %in% LNDCnn)]="Introduced"
relict$Status[which(relict$Species=="CRCH" & relict$RepeatID %in% CRCHnn)]="Introduced"
relict$Status[which(relict$Species=="BLBH")]="Native"
relict$Status[which(relict$Species=="BLBH" & relict$RepeatID %in% BLBHnn)]="Introduced"
relict$Status[which(relict$Species=="CCAT" & relict$RepeatID %in% CCATnn)]="Introduced"
relict$Status[which(relict$Species=="LING" & relict$RepeatID %in% LINGnn)]="Introduced"
relict$Status[which(relict$Species=="WSU" & relict$RepeatID %in% WSUnn)]="Introduced"
relict$Status[which(relict$Species=="RMCT" & relict$RepeatID %in% RMCTnn)]="Introduced"
relict$Status[which(relict$Species=="DRUM")]="Native"
relict$Status[which(relict$Species=="DRUM" & relict$RepeatID %in% drumnn)]="Introduced"
relict$Status[which(relict$Species=="PKF")]="Native"
relict$Status[which(relict$Species=="PKF" & relict$RepeatID %in% pkfnn)]="Introduced"
relict$Status[which(relict$Species=="PTMN")]="Native"
relict$Status[which(relict$Species=="PTMN" & relict$RepeatID %in% ptmnnn)]="Introduced"
relict$Status[which(relict$Species=="BRSB")]="Native"
relict$Status[which(relict$Species=="BRSB" & relict$RepeatID %in% brsbnn)]="Introduced"
relict=subset(relict,relict$Status=="Introduced")
relict=relict[,-c(39:52)]

flowavgs=read.csv("FlowAverages.csv")
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
traits=left_join(traits,flowavgs,by="Species")
relict$X=NULL
obsrelict=relict
#obsrelict=subset(obsrelict,obsrelict$Status=="Introduced")
obsrelict=subset(obsrelict,!is.na(obsrelict$change))
obsrelcol=obsrelict%>%group_by(RepeatID)%>%summarise(persSmGlacial=mean(change))
obsrelcol$geometry=NULL
obsrelcol=obsrelcol%>%rename("pIntroduced"="persSmGlacial")
obsrelcol2=obsrelcol
obsrelcol2=subset(obsrelcol2,!is.na(obsrelcol2$pIntroduced))

#Add to nssn
obsrelcol=obsrelcol%>%rename("RepetID"="RepeatID")
nssn$obs=left_join(nssn$obs,obsrelcol,by="RepetID")


#------------------------------build models------------------------------------------------------
nssn$obs$pGlacRel[which(nssn$obs$pGlacRel==1)]=0.999
nssn$obs$pGlacRel[which(nssn$obs$pGlacRel==0)]=0.001

#full model -- max likelihood
disp_init <- dispersion_initial(family = "beta", dispersion = 1, known = "dispersion")
hist(nssn$obs$pGlacRel)
#nssn2$obs%>%filter(PU=="GREEN")%>%summarise(length(RepeatID))

#-- will only converge without euclidean distance
ssn_modpG <- ssn_glm(
  formula =   pGlacRel~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc)+scale(prosper),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  #euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  #nugget_type = "none"
)
summary(ssn_modpG)

#significant only-- will only converge without euclidean distance
ssn_modpG_sig <- ssn_glm(
  formula =   pGlacRel~ scale(temp),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  #taildown_type = "exponential",
  #euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  nugget_type = "none"
)
summary(ssn_modpG_sig)

#null model -- will only converge without euclidean distance
ssn_null <- ssn_glm(
  formula =   pGlacRel~ 1,
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  #taildown_type = "exponential",
  #euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  nugget_type = "none"
)


loocv(ssn_modpG)#0.434
varcomp(ssn_modpG)
loocv(ssn_modpG_sig)#0.446
loocv(ssn_null)#0.454



#VAR IMP
coeff=as.data.frame(ssn_modpG$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ssn_modpG$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

VarImpPGlacial=coeff%>%arrange(Score2) %>%
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
    axis.title.x = element_blank(),
    axis.text = element_text(size=10, color = "black"),
    legend.position = "none"
  ) +
  ggtitle(label="B. Glacial Relict Persistence")+
  xlab("") +
  coord_flip() 
  #ylab("Variable Importance")
VarImpPGlacial

ggsave(filename = "VarImp_pGlacialRelict.tiff",height = 4, width = 5, units = "in", dpi=400)



#-----------Species Turnover--------------------------
#cant use piscivore metrics since it is included in turnover response
turn=read.csv("TurnoverMetric.csv")
turn=turn%>%rename("RepetID"="RepeatID")
nssn$obs$Turnover.x=NULL
nssn$obs$Turnover.y=NULL
nssn$obs$Turnover=NULL
nssn$obs=left_join(nssn$obs, turn, by="RepetID")

nssn$obs$X=NULL

nssn$obs$Turnover[which(nssn$obs$Turnover==0)]=0.001
nssn$obs$Turnover[which(nssn$obs$Turnover==1)]=0.999

ssn_modT <- ssn_glm(
  formula =   Turnover~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(prosper),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  #euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)
obs$rchC
T.torg=Torgegram(pNat~ 1,ssn.object = nssn,bins = 30)
plot(T.torg)
summary(ssn_modT)
varcomp(ssn_modT)
#sig only
ssn_modT_sig <- ssn_glm(
  formula =   Turnover~ scale(temp)+scale(length)+scale(prosper),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)
summary(ssn_modT_sig)

#null model -- max likelihood
ssn_null <- ssn_glm(
  formula =   Turnover~ 1,
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)

loocv(ssn_modT)#0.207
loocv(ssn_modT_sig)#0.204
loocv(ssn_null)#0.211



#MODEL FOR VAR IMP
coeff=as.data.frame(ssn_modT$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ssn_modT$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature**"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length**"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence**"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

VarImpTurn=coeff%>%arrange(Score2) %>%
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
    axis.text = element_text(size=10, color = "black"),
    legend.position = "none"
  ) +
  ggtitle(label="C. Community Turnover")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")
VarImpTurn

ggsave(filename = "VarImp_Turn.tiff",height = 4, width = 5, units = "in", dpi=400)


#st_write(nssn$obs, dsn = "metrics.shp")


#---------------------Graph Different Scenarios-----------

hist(nssn$obs$length)
quantile(nssn$obs$length)
nssn$obs$FragCat=NA
nssn$obs$FragCat[which(nssn$obs$length<62)]="Smaller (<62 km)"#180
nssn$obs$FragCat[which(nssn$obs$length>=62)]="Large (>62 km)" #177
nssn$obs%>%group_by(FragCat)%>%summarise(n=length(temp))

quantile(nssn$obs$temp)
nssn$obs$TempCat=NA
nssn$obs$TempCat[which(nssn$obs$temp<20)]="Cold (<20 C)" #261
nssn$obs$TempCat[which(nssn$obs$temp>=20)]="Warm (>20 C)" #96
nssn$obs%>%group_by(TempCat)%>%summarise(n=length(temp))

###
tempgraph=nssn$obs%>%
  group_by(length,FragCat,TempCat, GNIS_NA)%>%
  summarise(pNAT=mean(pNat), TEMP=mean(temp))%>%
  ggplot(aes(x=TEMP, y=pNAT, color = FragCat, label = FragCat))+
  geom_jitter(size=2)+
  stat_smooth(method = "glm", method.args = list(family=binomial), se=F, linewidth=1.5)+
  scale_color_manual(name="Fragment Length",values=c("black", "#B0B0B0"))+
  geom_text(aes(x=12,y=0.66),label="Large", color="black")+
  geom_text(aes(x=11,y=0.765),label="Small", color="#B0B0B0")+
  ylab(label = "Native Species Persistence")+
  xlab(label = "Stream Temperature (C)")+
  ggtitle(label = "Persistence by Fragment Length & Temperature")+
  theme_classic()+
  theme(legend.position=c(0.27,0.25),
        legend.box.background = element_rect(),
        legend.box.margin = margin(5,5,5,5),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10, color = "black"))
tempgraph

sizegraph=nssn$obs%>%
  group_by(length,FragCat,TempCat, GNIS_NA)%>%
  summarise(pNAT=mean(pNat))%>%
  ggplot(aes(x=log(length), y=pNAT, color = TempCat))+
  geom_jitter(size=2)+
  stat_smooth(method = "glm", method.args = list(family=binomial), se=F, linewidth=1.5)+
  scale_color_manual(name="Stream Temperature", values = c("dodgerblue","darkorange"))+
  #geom_text(aes(x=20,y=0.64),label="Cold", color="lightblue")+
  #geom_text(aes(x=17,y=0.43),label="Cool-Warm", color="firebrick3")+
  ylab(label = "Native Species Persistence")+
  xlab(label = "log(Fragment Length (km))")+
  ggtitle(label = "D. Native Sp. Pers. by Fragment Length & Temp.")+
  theme_classic()+
  theme(
        legend.position=c(0.19,0.17),
        legend.background = element_blank(),
        legend.box.margin = margin(2,2,2,2),
        legend.key.size = unit(0.15, 'in'),
        legend.title = element_text(face="bold"),
        legend.key.spacing.y = unit(0.5,"pt"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10, color = "black"))
sizegraph
ggsave(filename = "Temp&Length.tiff", dpi=400, height = 4,, width = 5, units = "in")

nssn$obs$prosper

permgraph=nssn$obs%>%
  group_by(prosper,TempCat, GNIS_NA)%>%
  summarise(pNAT=mean(pNat))%>%
  ggplot(aes(x=prosper, y=pNAT, color = TempCat))+
  geom_jitter(size=2)+
  stat_smooth(method = "glm", method.args = list(family=binomial), se=F, linewidth=1.5)+
  scale_color_manual(name="Stream Temperature", values = c("dodgerblue","darkorange"))+
  #geom_text(aes(x=20,y=0.64),label="Cold", color="lightblue")+
  #geom_text(aes(x=17,y=0.43),label="Cool-Warm", color="firebrick3")+
  ylab(label = "Native Species Persistence")+
  xlab(label = "Probability of Flow Permanence")+
  ggtitle(label = "E. Native Sp. Pers. by Flow Permanence & Temp.")+
  theme_classic()+
  theme(
    legend.position=c(0.19,0.17),
    legend.background = element_blank(),
    legend.box.margin = margin(2,2,2,2),
    legend.key.size = unit(0.15, 'in'),
    legend.title = element_text(face="bold"),
    legend.key.spacing.y = unit(0.5,"pt"),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10, color = "black"))
permgraph
ggsave(filename = "Temp&PROSPER.tiff", dpi=400, height = 4,, width = 5, units = "in")
















###
tempgraph2=nssn$obs%>%
  group_by(temp, GNIS_NA)%>%
  summarise(pNAT=mean(pNat))%>%
  ggplot(aes(x=temp, y=pNAT))+
  geom_point()+
  stat_smooth(method = "glm", method.args = list(family=binomial), se=T, linewidth=1.5, color="black")+
  ylab(label = "Native Species Persistence")+
  xlab(label = "Stream Temperature (C)")+
  ggtitle(label = "D. Persistence by Temperature")+
  theme_classic()+
  theme(legend.position="none",
        axis.title = element_text(size=12),
        axis.text = element_text(size=10, color = "black"))
tempgraph2

sizegraph2=nssn$obs%>%
  group_by(length, GNIS_NA)%>%
  summarise(pNAT=mean(pNat))%>%
  ggplot(aes(x=length, y=pNAT))+
  #geom_point()+
  stat_smooth(method = "glm", method.args = list(family=binomial), linewidth=1.5, color="black")+
  ylab(label = "Native Species Persistence")+
  xlab(label = "Fragment Length (km)")+
  ggtitle(label = "Persistence by Fragment Length")+
  #ylim(limits=c(0.3,1))+
  theme_classic()+
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=10, color = "black"))
sizegraph2



library(ggpubr)

VarImpComb=ggarrange(VarImpNative,VarImpPGlacial,VarImpTurn,
          nrow = 3)
VarImpComb

ModComb=ggarrange(sizegraph,permgraph,nrow = 2)
ModComb

ggarrange(VarImpComb,ModComb, ncol = 2)
ggsave(filename = "VariableGraph2.tiff",dpi=400, height = 10, width = 10, units = "in")


#==================================================================================
#============================SECTION 4: SSN Models for Individual Species
#==================================================================================

obs3=obs
obs3$geometry=NULL
obs3=obs3%>%pivot_longer(cols=28:121, names_to = "Species")%>%
  filter(!is.na(value))%>%
  group_by(Species)%>%
  summarise(sites=length(RepeatID))

traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
nssn2=nssn
nssn2$obs=nssn2$obs%>%rename("RepeatID"="RepetID")
nssn2$obs=nssn2$obs%>%pivot_longer(cols = c(34:126), names_to = "Species")
nssn2$obs=left_join(nssn2$obs,traits,by="Species")
nssn2$obs$Status[which(nssn2$obs$Species=="RDSH")]="Native"
nssn2$obs$Status[which(nssn2$obs$Species=="RDSH" & nssn2$obs$RepeatID %in% rdshnn)]="Introduced"
nssn2$obs$Status[which(nssn2$obs$Species=="LKCH" & nssn2$obs$RepeatID %in% lkchnn)]="Introduced"
nssn2$obs$Status[which(nssn2$obs$Species=="FHMN" & nssn2$obs$RepeatID %in% FHMNnn)]="Introduced"
nssn2$obs$Status[which(nssn2$obs$Species=="LNDC" & nssn2$obs$RepeatID %in% LNDCnn)]="Introduced"
nssn2$obs$Status[which(nssn2$obs$Species=="CRCH" & nssn2$obs$RepeatID %in% CRCHnn)]="Introduced"
nssn2$obs$Status[which(nssn2$obs$Species=="BLBH")]="Native"
nssn2$obs$Status[which(nssn2$obs$Species=="BLBH" & nssn2$obs$RepeatID %in% BLBHnn)]="Introduced"
nssn2$obs$Status[which(nssn2$obs$Species=="CCAT" & nssn2$obs$RepeatID %in% CCATnn)]="Introduced"
nssn2$obs$Status[which(nssn2$obs$Species=="LING" & nssn2$obs$RepeatID %in% LINGnn)]="Introduced"
nssn2$obs$Status[which(nssn2$obs$Species=="WSU" & nssn2$obs$RepeatID %in% WSUnn)]="Introduced"
nssn2$obs$Status[which(nssn2$obs$Species=="RMCT" & nssn2$obs$RepeatID %in% RMCTnn)]="Introduced"
nssn2$obs$Status[which(nssn2$obs$Species=="DRUM")]="Native"
nssn2$obs$Status[which(nssn2$obs$Species=="DRUM" & nssn2$obs$RepeatID %in% drumnn)]="Introduced"
nssn2$obs$Status[which(nssn2$obs$Species=="PKF")]="Native"
nssn2$obs$Status[which(nssn2$obs$Species=="PKF" & nssn2$obs$RepeatID %in% pkfnn)]="Introduced"
nssn2$obs$Status[which(nssn2$obs$Species=="PTMN")]="Native"
nssn2$obs$Status[which(nssn2$obs$Species=="PTMN" & nssn2$obs$RepeatID %in% ptmnnn)]="Introduced"
nssn2$obs$Status[which(nssn2$obs$Species=="BRSB")]="Native"
nssn2$obs$Status[which(nssn2$obs$Species=="BRSB" & nssn2$obs$RepeatID %in% brsbnn)]="Introduced"
nssn2$obs=subset(nssn2$obs,nssn2$obs$Status=="Native")


###Limit to native species with over 30 sites present (if ever present at a site)
nativeSites=read.csv("nativesites.csv")
nativeSites$X=NULL

nssn2$obs=nssn2$obs%>%pivot_wider(names_from = "Species", values_from = "value")

#BRMN
BRMN.global <- ssn_glm(
  formula =   BRMN~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  additive = "afvArea", 
  random = ~as.factor(YrangeCat),estmethod = "ml")
summary(BRMN.global)
BRMN.null <- ssn_glm(
  formula =   BRMN~ 1,
  family = "binomial",
  ssn.object = nssn22,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
BRMN.sig <- ssn_glm(
  formula =   BRMN~ scale(barrier)+scale(length),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(BRMN.sig)
loocv(BRMN.global)#0.423
loocv(BRMN.null)#0.341
loocv(BRMN.sig)#0.423 --null model is best
plot(BRMN~temp,data = nssn2$obs, main="BRMN Refugia vs. Temperature")
plot(BRMN~size,data = nssn2$obs, main="BRMN Refugia vs. Stream Size")
plot(jitter(BRMN,.2)~barrier, data = nssn2$obs, main="BRMN Refugia vs. Barriers (Jittered)")
plot(BRMN~length,data = nssn2$obs, main="BRMN Refugia vs. Fragment Length")
plot(jitter(BRMN,.2)~pisc,data = nssn2$obs, main="BRMN Refugia vs. Piscivores (Jittered)")

#MOTCOT
nssn2$obs$MOTCOT = NA
nssn2$obs$MOTCOT=0
nssn2$obs$MOTCOT[which(nssn2$obs$COLCOT==1 | nssn2$obs$RMCOT==1)]=1
MOTCOT.global <- ssn_glm(
  formula =   MOTCOT~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(MOTCOT.global)

MOTCOT.null <- ssn_glm(
  formula =   MOTCOT~ 1,
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
loocv(MOTCOT.global)#0.150
loocv(MOTCOT.null)#0.156
plot(MOTCOT~temp,data = nssn2$obs, main="MOTCOT Refugia vs. Temperature")
plot(MOTCOT~size,data = nssn2$obs, main="MOTCOT Refugia vs. Stream Size")
plot(jitter(MOTCOT,.2)~barrier, data = nssn2$obs, main="MOTCOT Refugia vs. Barriers (Jittered)")
plot(MOTCOT~length,data = nssn2$obs, main="MOTCOT Refugia vs. Fragment Length")
plot(jitter(MOTCOT,.2)~pisc,data = nssn2$obs, main="MOTCOT Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(MOTCOT.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="MOTCOT.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length*"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores*"
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
  ggtitle(label="Mottled Sculpin Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")

ggsave(filename = "VarImp_rMOTCOT.tiff",height = 7, width = 4, units = "in", dpi=400)



#FHCH
FHCH.global <- ssn_glm(
  formula =   FHCH~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(FHCH.global)

FHCH.null <- ssn_glm(
  formula =   FHCH~ 1,
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(FHCH.global)#0.334
loocv(FHCH.null)#0.345
plot(FHCH~temp,data = nssn2$obs, main="FHCH Refugia vs. Temperature")
plot(FHCH~size,data = nssn2$obs, main="FHCH Refugia vs. Stream Size")
plot(jitter(FHCH,.2)~barrier, data = nssn2$obs, main="FHCH Refugia vs. Barriers (Jittered)")
plot(FHCH~length,data = nssn2$obs, main="FHCH Refugia vs. Fragment Length")
plot(jitter(FHCH,.2)~pisc,data = nssn2$obs, main="FHCH Refugia vs. Piscivores (Jittered)")


#VAR IMP
coeff=as.data.frame(FHCH.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="FHCH.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers*"
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
  ggtitle(label="Flathead Chub Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")

ggsave(filename = "VarImp_rFHCH.tiff",height = 7, width = 4, units = "in", dpi=400)


#FHMN
FHMN.global <- ssn_glm(
  formula =   FHMN~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(FHMN.global)
FHMN.null <- ssn_glm(
  formula =   FHMN~ 1,
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
FHMN.sig <- ssn_glm(
  formula =   FHMN~ scale(size)+scale(length),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(FHMN.global)#0.466
loocv(FHMN.null)#0.418 - null is best
loocv(FHMN.sig)#0.466
plot(FHMN~temp,data = nssn2$obs, main="FHMN Refugia vs. Temperature")
plot(FHMN~size,data = nssn2$obs, main="FHMN Refugia vs. Stream Size")
plot(jitter(FHMN,.2)~barrier, data = nssn2$obs, main="FHMN Refugia vs. Barriers (Jittered)")
plot(FHMN~length,data = nssn2$obs, main="FHMN Refugia vs. Fragment Length")
plot(jitter(FHMN,.2)~pisc,data = nssn2$obs, main="FHMN Refugia vs. Piscivores (Jittered)")


#FMSU
FMSU.global <- ssn_glm(
  formula =   FMSU~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(FMSU.global)

FMSU.null <- ssn_glm(
  formula =   FMSU~ 1,
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

FMSU.sig <- ssn_glm(
  formula =   FMSU~ scale(length),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
loocv(FMSU.global) #0.450
loocv(FMSU.null) #0.365
loocv(FMSU.sig) #0.410

plot(FMSU~temp,data = nssn2$obs, main="FMSU Refugia vs. Temperature")
plot(FMSU~size,data = nssn2$obs, main="FMSU Refugia vs. Stream Size")
plot(jitter(FMSU,.2)~barrier, data = nssn2$obs, main="FMSU Refugia vs. Barriers (Jittered)")
plot(FMSU~length,data = nssn2$obs, main="FMSU Refugia vs. Fragment Length")
plot(jitter(FMSU,.2)~pisc,data = nssn2$obs, main="FMSU Refugia vs. Piscivores (Jittered)")



#LKCH
LKCH.global <- ssn_glm(
  formula =   LKCH~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(LKCH.global)

LKCH.null <- ssn_glm(
  formula =   LKCH~ 1,
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
loocv(LKCH.global) #0.465
loocv(LKCH.null) #0.499

print(loocv_mod$RMSPE)
plot(LKCH~temp,data = nssn2$obs, main="LKCH Refugia vs. Temperature")
plot(LKCH~size,data = nssn2$obs, main="LKCH Refugia vs. Stream Size")
plot(jitter(LKCH,.2)~barrier, data = nssn2$obs, main="LKCH Refugia vs. Barriers (Jittered)")
plot(LKCH~length,data = nssn2$obs, main="LKCH Refugia vs. Fragment Length")
plot(jitter(LKCH,.2)~pisc,data = nssn2$obs, main="LKCH Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(LKCH.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="LKCH.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores*"
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
  ggtitle(label="Lake Chub Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")

ggsave(filename = "VarImp_rLKCH.tiff",height = 7, width = 4, units = "in", dpi=400)

#LNDC
LNDC.global <- ssn_glm(
  formula =   LNDC~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(LNDC.global)

LNDC.null <- ssn_glm(
  formula =   LNDC~ 1,
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

LNDC.sig <- ssn_glm(
  formula =   LNDC~ scale(temp)+scale(size),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(LNDC.global)#0.435
loocv(LNDC.null)#0.402
loocv(LNDC.sig)#0.435
plot(LNDC~temp,data = nssn2$obs, main="LNDC Refugia vs. Temperature")
plot(LNDC~size,data = nssn2$obs, main="LNDC Refugia vs. Stream Size")
plot(jitter(LNDC,.2)~barrier, data = nssn2$obs, main="LNDC Refugia vs. Barriers (Jittered)")
plot(LNDC~length,data = nssn2$obs, main="LNDC Refugia vs. Fragment Length")
plot(jitter(LNDC,.2)~pisc,data = nssn2$obs, main="LNDC Refugia vs. Piscivores (Jittered)")


#LNSU
LNSU.global <- ssn_glm(
  formula =   LNSU~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(LNSU.global)

LNSU.null <- ssn_glm(
  formula =   LNSU~ 1,
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

LNSU.sig <- ssn_glm(
  formula =   LNSU~ scale(temp)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(LNSU.global)#0.467
loocv(LNSU.null)#0.420
loocv(LNSU.sig) #0.473

plot(LNSU~temp,data = nssn2$obs, main="LNSU Refugia vs. Temperature")
plot(LNSU~size,data = nssn2$obs, main="LNSU Refugia vs. Stream Size")
plot(jitter(LNSU,.2)~barrier, data = nssn2$obs, main="LNSU Refugia vs. Barriers (Jittered)")
plot(LNSU~length,data = nssn2$obs, main="LNSU Refugia vs. Fragment Length")
plot(jitter(LNSU,.2)~pisc,data = nssn2$obs, main="LNSU Refugia vs. Piscivores (Jittered)")


#MTSU
MTSU.global <- ssn_glm(
  formula =   MTSU~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(MTSU.global)
MTSU.null <- ssn_glm(
  formula =   MTSU~ 1,
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
MTSU.sig <- ssn_glm(
  formula =   MTSU~ scale(size)+scale(barrier),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
loocv(MTSU.global)#0.429
loocv(MTSU.null)#0.391
loocv(MTSU.sig)#0.445


plot(MTSU~temp,data = nssn2$obs, main="MTSU Refugia vs. Temperature")
plot(MTSU~size,data = nssn2$obs, main="MTSU Refugia vs. Stream Size")
plot(jitter(MTSU,.2)~barrier, data = nssn2$obs, main="MTSU Refugia vs. Barriers (Jittered)")
plot(MTSU~length,data = nssn2$obs, main="MTSU Refugia vs. Fragment Length")
plot(jitter(MTSU,.2)~pisc,data = nssn2$obs, main="MTSU Refugia vs. Piscivores (Jittered)")


#PLMN
PLMN.global <- ssn_glm(
  formula =   PLMN~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(PLMN.global)

PLMN.null <- ssn_glm(
  formula =   PLMN~ 1, 
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(PLMN.global)#0.372
loocv(PLMN.null)#0.414
plot(PLMN~temp,data = nssn2$obs, main="PLMN Refugia vs. Temperature")
plot(PLMN~size,data = nssn2$obs, main="PLMN Refugia vs. Stream Size")
plot(jitter(PLMN,.2)~barrier, data = nssn2$obs, main="PLMN Refugia vs. Barriers (Jittered)")
plot(PLMN~length,data = nssn2$obs, main="PLMN Refugia vs. Fragment Length")
plot(jitter(PLMN,.2)~pisc,data = nssn2$obs, main="PLMN Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(PLMN.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="PLMN.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
  ggtitle(label="Plains Minnow Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")

ggsave(filename = "VarImp_rPLMN.tiff",height = 7, width = 4, units = "in", dpi=400)




#PLSU
PLSU.global <- ssn_glm(
  formula =   PLSU~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(PLSU.global)

PLSU.null <- ssn_glm(
  formula =   PLSU~ 1,
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(PLSU.global)#0.481
loocv(PLSU.null)#0.5

plot(PLSU~temp,data = nssn2$obs, main="PLSU Refugia vs. Temperature")
plot(PLSU~size,data = nssn2$obs, main="PLSU Refugia vs. Stream Size")
plot(jitter(PLSU,.2)~barrier, data = nssn2$obs, main="PLSU Refugia vs. Barriers (Jittered)")
plot(PLSU~length,data = nssn2$obs, main="PLSU Refugia vs. Fragment Length")
plot(jitter(PLSU,.2)~pisc,data = nssn2$obs, main="PLSU Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(PLSU.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="PLSU.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
  ggtitle(label="Plains Sucker Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")

ggsave(filename = "VarImp_rPLSU.tiff",height = 7, width = 4, units = "in", dpi=400)



#SPDC
SPDC.global <- ssn_glm(
  formula =   SPDC~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(SPDC.global)

SPDC.null <- ssn_glm(
  formula =   SPDC~ 1,
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")


loocv(SPDC.global)#0.438
loocv(SPDC.null) #0.468
plot(SPDC~temp,data = nssn2$obs, main="SPDC Refugia vs. Temperature")
plot(SPDC~size,data = nssn2$obs, main="SPDC Refugia vs. Stream Size")
plot(jitter(SPDC,.2)~barrier, data = nssn2$obs, main="SPDC Refugia vs. Barriers (Jittered)")
plot(SPDC~length,data = nssn2$obs, main="SPDC Refugia vs. Fragment Length")
plot(jitter(SPDC,.2)~pisc,data = nssn2$obs, main="SPDC Refugia vs. Piscivores (Jittered)")


#VAR IMP
coeff=as.data.frame(SPDC.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="SPDC.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size*"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
  ggtitle(label="Speckled Dace Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")

ggsave(filename = "VarImp_rSPDC.tiff",height = 7, width = 4, units = "in", dpi=400)



#WSU
WSU.global <- ssn_glm(
  formula =   WSU~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(WSU.global)

WSU.null <- ssn_glm(
  formula =   WSU~ 1,
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")


loocv(WSU.global)#0.376
loocv(WSU.null) #0.378

plot(WSU~temp,data = nssn2$obs, main="WSU Refugia vs. Temperature")
plot(WSU~size,data = nssn2$obs, main="WSU Refugia vs. Stream Size")
plot(jitter(WSU,.2)~barrier, data = nssn2$obs, main="WSU Refugia vs. Barriers (Jittered)")
plot(WSU~length,data = nssn2$obs, main="WSU Refugia vs. Fragment Length")
plot(jitter(WSU,.2)~pisc,data = nssn2$obs, main="WSU Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(WSU.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="WSU.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size*"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers*"
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
  ggtitle(label="White Sucker Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")
ggsave(filename = "VarImp_rWSU.tiff",height = 7, width = 4, units = "in", dpi=400)


#BLBH
BLBH.global <- ssn_glm(
  formula =   BLBH~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(BLBH.global)

BLBH.null <- ssn_glm(
  formula =   BLBH~ 1,
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(BLBH.global)#0.426
loocv(BLBH.null) #0.497

plot(BLBH~temp,data = nssn2$obs, main="BLBH Refugia vs. Temperature")
plot(BLBH~size,data = nssn2$obs, main="BLBH Refugia vs. Stream Size")
plot(jitter(BLBH,.2)~barrier, data = nssn2$obs, main="BLBH Refugia vs. Barriers (Jittered)")
plot(BLBH~length,data = nssn2$obs, main="BLBH Refugia vs. Fragment Length")
plot(jitter(BLBH,.2)~pisc,data = nssn2$obs, main="BLBH Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(BLBH.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="BLBH.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers*"
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
  ggtitle(label="Black Bullhead Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")

ggsave(filename = "VarImp_rBLBH.tiff",height = 7, width = 4, units = "in", dpi=400)



########INTRODUCED SPECIES
#EB
EB.global <- ssn_glm(
  formula =   EB~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(EB.global)

EB.null <- ssn_glm(
  formula =   EB~ 1,
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(EB.global)#0.408
loocv(EB.null)#0.448

plot(EB~temp,data = nssn2$obs, main="EB Refugia vs. Temperature")
plot(EB~size,data = nssn2$obs, main="EB Refugia vs. Stream Size")
plot(jitter(EB,.2)~barrier, data = nssn2$obs, main="EB Refugia vs. Barriers (Jittered)")
plot(EB~length,data = nssn2$obs, main="EB Refugia vs. Fragment Length")
plot(jitter(EB,.2)~pisc,data = nssn2$obs, main="EB Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(EB.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="EB.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
  ggtitle(label="Brook Trout Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")

ggsave(filename = "VarImp_rEB.tiff",height = 7, width = 4, units = "in", dpi=400)


#RSSH
RSSH.global <- ssn_glm(
  formula =   RSSH~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(RSSH.global)

RSSH.null <- ssn_glm(
  formula =   RSSH~ 1,
  family = "binomial",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")


loocv(RSSH.global)#0.358
loocv(RSSH.null)#0.385

plot(RSSH~temp,data = nssn2$obs, main="RSSH Refugia vs. Temperature")
plot(RSSH~size,data = nssn2$obs, main="RSSH Refugia vs. Stream Size")
plot(jitter(RSSH,.2)~barrier, data = nssn2$obs, main="RSSH Refugia vs. Barriers (Jittered)")
plot(RSSH~length,data = nssn2$obs, main="RSSH Refugia vs. Fragment Length")
plot(jitter(RSSH,.2)~pisc,data = nssn2$obs, main="RSSH Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(RSSH.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="RSSH.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
  ggtitle(label="Redside Shiner Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")

ggsave(filename = "VarImp_rRSSH.tiff",height = 7, width = 4, units = "in", dpi=400)


#=====================================================================
#============================SECTION 5: Examine Net Change for Species
#=====================================================================
#---------------------Net Species Change------------------------------------------
NET=read.csv("PersistenceMetrics.csv")
removes=read.csv("REMOVES.csv")
removes=removes[,c(1,3)]
#removes=removes%>%rename("RepetID"="RepeatID")
NET=left_join(NET,removes,by="RepeatID")
NET=subset(NET,NET$GearMiss!="Y")
NET=subset(NET, NET$F_MAUG_<1000)
NET=subset(NET, NET$RepeatID!=234 & NET$RepeatID!=74& NET$RepeatID!=237& NET$RepeatID!=238)


NET=NET%>%pivot_longer(cols = c(28:121),names_to = "Species",values_to = "change")
NET=subset(NET,!is.na(NET$change))
NET$Species[which(NET$Species=="COLCOT"|NET$Species=="RMCOT")]="MOTCOT"
NET2=NET%>%group_by(Species)%>%summarise(net=sum(change), sites=length(unique(RepeatID)))
NET3=NET%>%filter(change!=-1)%>%group_by(Species)%>%summarise(LateNum=length(change))
NET=NET%>%filter(change!=1)%>%group_by(Species)%>%summarise(EarlyNum=length(change))


NET=left_join(NET2,NET,by="Species")
NET=left_join(NET,NET3,by="Species")
NET$EarlyNum[which(is.na(NET$EarlyNum))]=0
NET$LateNum[which(is.na(NET$LateNum))]=0
NET$propChange=NA
NET$propChange=NET$LateNum/NET$EarlyNum
NET=full_join(NET,traits,by="Species")
write.csv(NET,"NET2.csv")
NETsub=subset(NET,NET$sites>=30)
NETsub=subset(NETsub,NETsub$Species!="ONC" & NETsub$Species!="FMxWSU")
mean(NETsub$propChange)


traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
traits$Species[which(traits$Species=="COLCOT")]="MOTCOT"
NETsub=left_join(NETsub,traits,by="Species")
NETsub$CommonName[which(NETsub$CommonName=="Northern Redbelly Dace")]="Nor. Redbelly Dace"
NETsub$CommonName[which(NETsub$CommonName=="Rocky Mountain Cutthroat Trout")]="R.M. Cutthroat Trout"

NETsub$Direc=NA
NETsub$Direc[which(NETsub$propChange>=1)]="Pos"
NETsub$Direc[which(NETsub$propChange<1)]="Neg"

NETsubALL=subset(NETsub,NETsub$Status=="Introduced")

changenative=NETsub%>%
  ggplot(aes(x=reorder(CommonName,-propChange),y=propChange, colour = Direc))+
  geom_segment( aes(x=CommonName, xend=CommonName, y=1, yend=propChange), color="black") +
  geom_point(size=3)+
  theme_light()+
  ylab(label = "Percent Change in Sites Occupied")+
  scale_y_continuous(breaks=seq(0.25,2.5,by=0.25), limits = c(0.25,2.5), labels = c("-75%","-50%","-25%","0%","+25%","+50%","+75%","+100%","+125%","+150%"))+
  xlab(label="")+
  scale_color_manual(values = c("#ff9999", "#018081"))+
  coord_flip()+
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face = 2),
        legend.title = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size=12, color = "black"),
        axis.text.x = element_text(size=10, color = "black"),
        axis.ticks.x = element_blank(),
        legend.position = "none")
changenative
ggsave(filename="NativeSpChange.tiff",dpi = 400, width = 8, height = 6, units = "in")





#---------------------Net Species Changes to Only Native Populations------------------------------------------
NET=read.csv("PersistenceMetrics.csv")
removes=read.csv("REMOVES.csv")
removes=removes[,c(1,3)]
#removes=removes%>%rename("RepetID"="RepeatID")
NET=left_join(NET,removes,by="RepeatID")
NET=subset(NET,NET$GearMiss!="Y")
NET=subset(NET, NET$F_MAUG_<1000)
NET=subset(NET, NET$RepeatID!=234 & NET$RepeatID!=74& NET$RepeatID!=237& NET$RepeatID!=238)

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

NET2=NET%>%group_by(Species)%>%summarise(net=sum(change))
NET3=NET%>%filter(change!=-1)%>%group_by(Species)%>%summarise(LateNum=length(change))
NET=NET%>%filter(change!=1)%>%group_by(Species)%>%summarise(EarlyNum=length(change))
NET=left_join(NET2,NET,by="Species")
NET=left_join(NET,NET3,by="Species")
NET$EarlyNum[which(is.na(NET$EarlyNum))]=0
NET$LateNum[which(is.na(NET$LateNum))]=0
NET$propChange=NA
NET$propChange=NET$LateNum/NET$EarlyNum
NET=left_join(NET,nativeSites,by="Species")
NETsub=subset(NET,NET$sites>=20)
NETsub=subset(NETsub,NETsub$Species!="ONC" & NETsub$Species!="FMxWSU")
mean(NETsub$propChange)
write.csv(NETsub, "table1material.csv")

#View(NETsub)

traits2=traits[,c(1,3)]
NETsub=left_join(NETsub,traits2,by="Species")
NETsub$CommonName[which(NETsub$CommonName=="Northern Redbelly Dace")]="Nor. Redbelly Dace"
NETsub$CommonName[which(NETsub$CommonName=="Rocky Mountain Cutthroat Trout")]="R.M. Cutthroat Trout"


NETsub$Direc=NA
NETsub$Direc[which(NETsub$propChange>=1)]="Pos"
NETsub$Direc[which(NETsub$propChange<1)]="Neg"
NETsub$Status=NA
NETsub$Status="Native"
NETsub$Status[which(NETsub$Species%in%invasive20s)]="Introduced"
NETsub$percChange=NA
NETsub$percChange=NETsub$propChange*100
NETsub$Status2=NA
NETsub$Status2="PGP"
NETsub$Status2[which(NETsub$Species=="LKCH"|NETsub$Species=="BRSB"|NETsub$Species=="PLSU"|NETsub$Species=="NRBDC"|NETsub$Species=="IOWA")]="GR"
NETsub$Status2[which(NETsub$Status=="Introduced")]="NN"
NETsub$Status2[which(NETsub$Species=="STR"|NETsub$Species=="CCAT"|NETsub$Species=="CMSH"|NETsub$Species=="CRCH"|NETsub$Species=="FMSU"|NETsub$Species=="JDT"|NETsub$Species=="PLMN"|NETsub$Species=="RCSU"|NETsub$Species=="RTCH"|NETsub$Species=="SDSH")]="OtherN"

library(ggbreak)
library(ggh4x)
changenative=NETsub%>%
  ggplot(aes(x=reorder(CommonName,-propChange),y=propChange, shape = Status2, color=Status2))+
  geom_hline(yintercept = 1, color=alpha("black", alpha = 0.3), linewidth=1.1)+
  geom_segment( aes(x=CommonName, xend=CommonName, y=1, yend=propChange), color="black") +
  geom_point(size=3)+
  theme_light()+
  ylab(label = "Occupancy Change")+
  xlab(label="")+
  scale_color_manual(values = c("lightblue", "#ff9999","lightgreen","dodgerblue4"))+
  scale_shape_manual(values=c(15,17, 16, 15))+
  scale_y_continuous(limits = c(0,6), breaks = seq(0,5.4, by=0.2),minor_breaks = seq(0,5.4,by=0.1),guide = "axis_minor")+
  coord_flip()+
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face = 2),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_text(size=12, color = "black"),
    axis.text.x = element_text(size=10, color = "black"),
    legend.position = "none")+
  theme(ggh4x.axis.ticks.length.minor=rel(1))+
  scale_y_break(c(2.5, 5))
changenative


ggsave(filename="NativeSpChange.tiff",dpi = 400, width = 10, height = 8, units = "in")





#=====================================================================
#============================SECTION 6: Group Changes Stats
#=====================================================================
traits$Glacial[which(traits$Species=="BLBH")]="L"
NETsub2=left_join(NET,traits,by="Species")
thresh=read.csv("Thresholds.csv")
NETsub2=left_join(NETsub2,thresh,by = "Species")
flowavgs=read.csv("FlowAverages.csv")
NETsub2=left_join(NETsub2,flowavgs,by="Species")
NETsub2=filter(NETsub2, NETsub2$sites>=20)

#All pioneers including big rivers
NETsubR=NETsub2%>%filter(Glacial=="E")
NETsubNR=NETsub2%>%filter(Glacial=="L")
t.test(NETsubR$propChange,NETsubNR$propChange) # dif = 0.285, p=0.159
wilcox.test(NETsubR$propChange,NETsubNR$propChange)
median(NETsubR$propChange)
median(NETsubNR$propChange)

#Native v. introduced
NETsub2$Status[which(NETsub2$Species=="BRSB")]="Native"
NETsubnative=subset(NETsub2, NETsub2$Status=="Native")
NETsubintroduced=subset(NETsub2, NETsub2$Status!="Native")
t.test(NETsubnative$propChange,NETsubintroduced$propChange)
wilcox.test(NETsubnative$propChange,NETsubintroduced$propChange)
median(NETsubnative$propChange)
median(NETsubintroduced$propChange)

#Glacial Relicts
smallrelicts=c("NRBDC","FSDC","LKCH","BRSB","PLSU","GR","NPDC","HHCH","IOWA")
NETsub2$PGtype=NA
NETsub2$PGtype[which(NETsub2$Glacial=="E")]="Group Members"
NETsub2$PGtype[which(NETsub2$Glacial=="L")]="Others"
'%!in%' <- function(x,y)!('%in%'(x,y))
NETsub2$PGtype[which(NETsub2$Species%in%smallrelicts)]="Group Members"
NETsub2$PGtype[which(NETsub2$Species%!in%smallrelicts)]="Others"
NETsub2$PGtype[which(NETsub2$Status=="Introduced")]="Introduced"
NETsubGR=subset(NETsub2, NETsub2$PGtype=="Group Members")
NETsubNGR=subset(NETsub2, NETsub2$PGtype=="Others")
t.test(NETsubGR$propChange,NETsubNGR$propChange)
wilcox.test(NETsubGR$propChange,NETsubNGR$propChange)
median(NETsubGR$propChange)
median(NETsubNGR$propChange)

#Combine for Graph
NETsubR=NETsubR[,c(5,19)]
NETsubNR=NETsubNR[,c(5,19)]
NETsubcom1=rbind(NETsubR,NETsubNR)
NETsubcom1=NETsubcom1%>%rename("Status"="Glacial")
NETsubcom1$Status[which(NETsubcom1$Status=="E")]="Group Members"
NETsubcom1$Status[which(NETsubcom1$Status=="L")]="Others"
NETsubcom1$Group=NA
NETsubcom1$Group="Postglacial-Pioneer Species"

NETsubnative=NETsubnative[,c(5,12)]
NETsubintroduced=NETsubintroduced[,c(5,12)]
NETsubnative$Group=NA
NETsubnative$Group="Native Species"
NETsubnative$Status="Group Members"
NETsubcom1=rbind(NETsubcom1,NETsubnative)
NETsubintroduced$Group=NA
NETsubintroduced$Status="Others"
NETsubintroduced$Group="Native Species"
NETsubcom1=rbind(NETsubcom1,NETsubintroduced)

NETsubGR=NETsubGR[,c(5,30)]
NETsubNGR=NETsubNGR[,c(5,30)]
NETsubGR=NETsubGR%>%rename("Status"="PGtype")
NETsubNGR=NETsubNGR%>%rename("Status"="PGtype")
NETsubGR$Group=NA
NETsubGR$Group="Glacial Relicts"
NETsubNGR$Group=NA
NETsubNGR$Group="Glacial Relicts"
NETsubcom1=rbind(NETsubcom1,NETsubGR)
NETsubcom1=rbind(NETsubcom1,NETsubNGR)
NETsubcom1$GroupStatus=NA
NETsubcom1$GroupStatus=paste(NETsubcom1$Group,NETsubcom1$Status,sep = " ")

NETsubcom1$Group[which(NETsubcom1$Group=="Glacial Relicts")]="A. Glacial Relicts"
NETsubcom1$Group[which(NETsubcom1$Group=="Postglacial-Pioneer Species")]="B. Postglacial-Pioneer Sp"
NETsubcom1$Group[which(NETsubcom1$Group=="Native Species")]="C. Native Species"

graph.glac=NETsubcom1%>%
  ggplot(aes(x=Status, y=propChange, fill=GroupStatus))+
  geom_hline(yintercept = 1, color="grey", linewidth=1.1)+
  geom_boxplot()+
  ylab(label="Occupancy Change")+
  scale_fill_manual(values=c("lightblue", "grey","lightgreen","#ff9999","dodgerblue4","grey"))+
  #scale_y_continuous(limits = c(-75,450))+
  theme_classic()+
  ggtitle(label = "Occupancy Change by Group")+
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.text.x = element_blank(),
        legend.position = "none")+
  #scale_y_break(c(150, 400), ticklabels = c(-75,-50,-25,0,25,50,100,150,400,450))+
  facet_wrap(~factor(Group, c("A. Glacial Relicts", "B. Postglacial-Pioneer Sp", "C. Native Species")), scales = "free")+
  theme(strip.text = element_text(face="bold", size=11),
        strip.background = element_rect(fill="lightgrey",linewidth=1))

graph.glac
ggsave(filename = "NEtChange.tiff", dpi=400, height = 6, width = 8)




#Anova
NETsub2$propChange
glac.lm=lm(propChange~PGtype, data=NETsub2)
glac.aov=aov(glac.lm)
summary(glac.aov)
TukeyHSD(glac.aov)

#GR vs non-GR
NETsub4=subset(NETsub2, NETsub2$PGtype=="Group Members")
NETsub5=subset(NETsub2, NETsub2$PGtype!="Group Members")

t.test(NETsub4$propChange, NETsub5$propChange)
wilcox.test(NETsub4$propChange, NETsub5$propChange)

#Postglacial vs NonPG
NETsub4=subset(NETsub2, NETsub2$Glacial=="E")
NETsub5=subset(NETsub2, NETsub2$Glacial=="L")

t.test(NETsub4$propChange, NETsub5$propChange)
wilcox.test(NETsub4$propChange, NETsub5$propChange)


#Postglacial vs NonPG
NETsub4=subset(NETsub2, NETsub2$Status=="Native")
NETsub5=subset(NETsub2, NETsub2$Status=="Introduced")

t.test(NETsub4$propChange, NETsub5$propChange)
wilcox.test(NETsub4$propChange, NETsub5$propChange)


#==================================================================================
#============================SECTION 7: Site-based metrics
#==================================================================================
#Native Species
nssn$obs$pNat[which(nssn$obs$pNat==0.999)]=1
nssn$obs$pNat[which(nssn$obs$pNat==0.001)]=0
Native=nssn$obs%>%filter(!is.na(pNat))%>%
  summarise(Type="Native",Prop = mean(pNat), se = sd(pNat)/sqrt(length(RepetID)))  #0.593, se=0.017
Native$geometry=NULL
Native

#Turnover
nssn$obs$Turnover[which(nssn$obs$Turnover==0.999)]=1
nssn$obs$Turnover[which(nssn$obs$Turnover==0.001)]=0
nssn$obs%>%filter(!is.na(Turnover))%>%
  summarise(mTurn = mean(Turnover),se = sd(Turnover)/sqrt(length(RepetID))) #0.564, se=0.013

#Postglacial Pioneers
nssn$obs$pGlac[which(nssn$obs$pGlac==0.999)]=1
nssn$obs$pGlac[which(nssn$obs$pGlac==0.001)]=0
#pioneers
Postglacial=nssn$obs%>%filter(!is.na(pGlac))%>%
  summarise(Type="APostglacial-Pioneer Species",Prop = mean(pGlac), se = sd(pGlac)/sqrt(length(RepetID))) #0.613, se=0.02
Postglacial$geometry=NULL
Postglacial

#nonpioneers
nssn$obs$pNotGlac[which(nssn$obs$pNotGlac==0.999)]=1
nssn$obs$pNotGlac[which(nssn$obs$pNotGlac==0.001)]=0
NotPP=nssn$obs%>%filter(!is.na(pNotGlac))%>%
  summarise(Type="Other",Prop = mean(pNotGlac), se = sd(pNotGlac)/sqrt(length(RepetID)))
NotPP$geometry=NULL
NotPP

#true glacial 
nssn$obs$pGlacRel[which(nssn$obs$pGlacRel==0.999)]=1
nssn$obs$pGlacRel[which(nssn$obs$pGlacRel==0.001)]=0
Relict=nssn$obs%>%filter(!is.na(pGlacRel))%>%
  summarise(Type="Glacial Relicts",Prop = mean(pGlacRel), se = sd(pGlacRel)/sqrt(length(RepetID))) #0.417, se=0.05
Relict$geometry=NULL
Relict

#nonrelicts
nssn$obs$pNotGlacRel[which(nssn$obs$pGlacRel==0.999)]=1
nssn$obs$pNotGlacRel[which(nssn$obs$pNotGlacRel==0.001)]=0
NotRel=nssn$obs%>%filter(!is.na(pNotGlacRel))%>%
  summarise(Type="Other",Prop = mean(pNotGlacRel), se = sd(pNotGlacRel)/sqrt(length(RepetID)))
NotRel$geometry=NULL
NotRel

#introduced Species
nssn$obs$pIntroduced[which(nssn$obs$pIntroduced==0.999)]=1
nssn$obs$pIntroduced[which(nssn$obs$pIntroduced==0.001)]=0
Intro=nssn$obs%>%filter(!is.na(pIntroduced))%>%
  summarise(Type="Non-Native",Prop = mean(pIntroduced), se = sd(pIntroduced)/sqrt(length(RepetID)))
Intro$geometry=NULL
Intro


#Native v. Introduced
t.test(nssn$obs$pNat, nssn$obs$pIntroduced) #p=0.19
wilcox.test(nssn$obs$pNat, nssn$obs$pIntroduced) #0.49

#Postglacial vs. Non PP
t.test(nssn$obs$pGlac, nssn$obs$pNotGlac) #p=0.37
wilcox.test(nssn$obs$pGlac, nssn$obs$pNotGlac) #p =0.61 


#Glacial Relict vs. Non-Relict
t.test(nssn$obs$pGlacRel, nssn$obs$pNotGlacRel) #p<0.01
wilcox.test(nssn$obs$pGlacRel, nssn$obs$pNotGlacRel) #p<0.01




A=rbind(Native,Intro)
A$Group="C. Native Species"
B=rbind(Postglacial,NotPP)
B$Group="B. Postglacial-Pioneer Sp"
C=rbind(Relict, NotRel)
C$Group="A. Glacial Relicts"
A=rbind(A,B)
A=rbind(A,C)



graph.glac=A%>%
  ggplot(aes(x=Type, y=Prop, color=Type))+
  geom_errorbar(aes(ymin=Prop-se,ymax=Prop+se), width=0.1, color="darkgrey")+
  geom_point(size=5)+
  ylab(label="Mean Site-Level Persistence")+
  scale_color_manual(values=c("dodgerblue4","lightblue","lightgreen","#ff9999","grey"))+
  theme_classic()+
  ggtitle(label = "Mean Site-Level Persistence by Group")+
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.text.x = element_blank(),
        legend.position = "none")+
  #scale_y_break(c(150, 400), ticklabels = c(-75,-50,-25,0,25,50,100,150,400,450))+
  facet_wrap(~factor(Group, c("A. Glacial Relicts", "B. Postglacial-Pioneer Sp", "C. Native Species")),scales = "free")+
  theme(strip.text = element_text(face="bold", size=11),
        strip.background = element_rect(fill="lightgrey",linewidth=1),axis.line.y=element_line(colour = "black"))+
  scale_y_continuous(limits = c(0.3,0.7))

graph.glac
ggsave(filename = "SiteChange.tiff", dpi=400, height = 6, width = 8)








###########Survey years, etc.
min(nssn$obs$Yrange)
max(nssn$obs$Yrange)
median(nssn$obs$Yrange)
mean(nssn$obs$Yrange)
sd(nssn$obs$Yrange)
View(nssn$obs)

mean(wc_early$Year)
median(wc_early$Year)
mean(wc_late$Year)
median(wc_late$Year)
#==================================================================================
#============================SECTION 8: Extra Graphs
#==================================================================================
#Turnover
tempgraph.turn=nssn$obs%>%
  ggplot(aes(x=temp, y=Turnover,))+
  stat_smooth(method = "glm", method.args = list(family=binomial), se=T, linewidth=1.5, color="black")+
  geom_jitter()+
  ylab(label = "Community Turnover")+
  xlab(label = "Stream Temperature (C)")+
  ggtitle(label = "Turnover by Temperature")+
  #ylim(limits=c(0.3,1))+
  theme_classic()+
  theme(legend.position="none",
        axis.title = element_text(size=12),
        axis.text = element_text(size=10, color = "black"))
tempgraph.turn

sizegraph.turn=nssn$obs%>%
  ggplot(aes(x=length, y=Turnover))+
  stat_smooth(method = "glm", method.args = list(family=binomial), linewidth=1.5, color="black")+
  geom_jitter()+
  ylab(label = "Community Turnover")+
  xlab(label = "Fragment Length (km)")+
  ggtitle(label = "Turnover by Fragment Length")+
  #ylim(limits=c(0.3,1))+
  theme_classic()+
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=10, color = "black"))
sizegraph.turn


#Glacial Relicts
tempgraph.relict=nssn$obs%>%
  ggplot(aes(x=temp, y=pGlacRel,))+
  stat_smooth(method = "glm", method.args = list(family=binomial), se=T, linewidth=1.5, color="black")+
  geom_jitter()+
  ylab(label = "Persistence")+
  xlab(label = "Stream Temperature (C)")+
  ggtitle(label = "Glacial Relict Persistence by Temperature")+
  #ylim(limits=c(0.3,1))+
  theme_classic()+
  theme(legend.position="none",
        axis.title = element_text(size=12),
        axis.text = element_text(size=10, color = "black"))
tempgraph.relict

barriergraph.relict=nssn$obs%>%
  ggplot(aes(x=barrier, y=pGlacRel))+
  stat_smooth(method = "glm", method.args = list(family=binomial), linewidth=1.5, color="black")+
  geom_jitter()+
  ylab(label = "Persistence")+
  xlab(label = "Barrier Presence")+
  ggtitle(label = "Glacial Relict Persistence by Barrier Presence")+
  #ylim(limits=c(0.3,1))+
  theme_classic()+
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=10, color = "black"))
barriergraph.relict

piscivoregraph.relict=nssn$obs%>%
  ggplot(aes(x=pisc, y=pGlacRel))+
  stat_smooth(method = "glm", method.args = list(family=binomial), linewidth=1.5, color="black")+
  geom_jitter()+
  ylab(label = "Persistence")+
  xlab(label = "Piscivore Presence")+
  ggtitle(label = "Glacial Relict Persistence by Piscivore Presence")+
  #ylim(limits=c(0.3,1))+
  theme_classic()+
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=10, color = "black"))
piscivoregraph.relict


#==================================================================================
#============================SECTION 9: Sites Changes by Species
#==================================================================================
rumham=wc_rumham
removes=read.csv("REMOVES.csv")
removes=removes[,c(1,3)]
#removes=removes%>%rename("RepetID"="RepeatID")
rumham=left_join(rumham,removes,by="RepeatID")
rumham=subset(rumham,rumham$GearMiss!="Y")
rumham=subset(rumham, rumham$F_MAUG_<1000)
rumham=subset(rumham, rumham$RepeatID!=234 & rumham$RepeatID!=74& rumham$RepeatID!=237& rumham$RepeatID!=238)

rumham=rumham%>%pivot_longer(cols = c(27:120),names_to = "Species",values_to = "change")

####FROM HERE, MUST RUN "DatasetPrep.R" for section
rumham=subset(rumham,!is.na(rumham$change))
rumham$Species[which(rumham$Species=="RMCOT" | rumham$Species=="COLCOT")]="MOTCOT"
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
rumham=left_join(rumham,traits,by="Species")
rumham$Status[which(rumham$Species=="RDSH")]="Native"
rumham$Status[which(rumham$Species=="RDSH" & rumham$RepeatID %in% rdshnn)]="Introduced"
rumham$Status[which(rumham$Species=="LKCH" & rumham$RepeatID %in% lkchnn)]="Introduced"
rumham$Status[which(rumham$Species=="FHMN" & rumham$RepeatID %in% FHMNnn)]="Introduced"
rumham$Status[which(rumham$Species=="LNDC" & rumham$RepeatID %in% LNDCnn)]="Introduced"
rumham$Status[which(rumham$Species=="CRCH" & rumham$RepeatID %in% CRCHnn)]="Introduced"
rumham$Status[which(rumham$Species=="BLBH")]="Native"
rumham$Status[which(rumham$Species=="BLBH" & rumham$RepeatID %in% BLBHnn)]="Introduced"
rumham$Status[which(rumham$Species=="CCAT" & rumham$RepeatID %in% CCATnn)]="Introduced"
rumham$Status[which(rumham$Species=="LING" & rumham$RepeatID %in% LINGnn)]="Introduced"
rumham$Status[which(rumham$Species=="WSU" & rumham$RepeatID %in% WSUnn)]="Introduced"
rumham$Status[which(rumham$Species=="RMCT" & rumham$RepeatID %in% RMCTnn)]="Introduced"
rumham$Status[which(rumham$Species=="DRUM")]="Native"
rumham$Status[which(rumham$Species=="DRUM" & rumham$RepeatID %in% drumnn)]="Introduced"
rumham$Status[which(rumham$Species=="PKF")]="Native"
rumham$Status[which(rumham$Species=="PKF" & rumham$RepeatID %in% pkfnn)]="Introduced"
rumham$Status[which(rumham$Species=="PTMN")]="Native"
rumham$Status[which(rumham$Species=="PTMN" & rumham$RepeatID %in% ptmnnn)]="Introduced"
rumham$Status[which(rumham$Species=="BRSB")]="Native"
rumham$Status[which(rumham$Species=="BRSB" & rumham$RepeatID %in% brsbnn)]="Introduced"


rumham2=rumham%>%filter(Status=="Native")%>%group_by(Species, change)%>%
  summarise(n=length(RepeatID))
rumham3=rumham%>%filter(Status=="Native")%>%group_by(Species)%>%
  summarise(n.total=length(RepeatID))
rumham=left_join(rumham2,rumham3, by="Species")
rumham=subset(rumham, rumham$n.total>=20)

rumham$change[which(rumham$change==-1)]="Extirpated"
rumham$change[which(rumham$change==0)]="Persisted"
rumham$change[which(rumham$change==1)]="Colonized"
rumham=left_join(rumham, traits, by="Species")
rumham%>%
ggplot(aes(y=n, x=reorder(CommonName,-order(CommonName)), fill=factor(change,levels=c("Colonized","Persisted","Extirpated"))))+ 
  geom_bar(position="fill", stat="identity")+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.title.y = element_blank())+
  ylab(label = "Proportion of Sites Where Species Ever Present")+
  scale_fill_manual(values = c("#009E73","#F0E442","#D55E00"))+
  coord_flip()


rumham%>%filter(CommonName=="Brook Stickleback"|CommonName=="Lake Chub"|CommonName=="Plains Sucker"|CommonName=="Northern Redbelly Dace")%>%
  ggplot(aes(y=n, x=reorder(CommonName,-order(CommonName)), fill=factor(change,levels=c("Colonized","Persisted","Extirpated"))))+ 
  geom_bar(position="fill", stat="identity")+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.title.y = element_blank())+
  ylab(label = "Proportion of Sites Where Species Ever Present")+
  scale_fill_manual(values = c("#005555","#cccccc","firebrick"))+
  coord_flip()


rumham%>%filter(Glacial=="E")%>%
  ggplot(aes(y=n, x=reorder(CommonName,-order(CommonName)), fill=factor(change,levels=c("Colonized","Persisted","Extirpated"))))+ 
  geom_bar(position="fill", stat="identity")+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.title.y = element_blank())+
  ylab(label = "Proportion of Sites Where Species Ever Present")+
  scale_fill_manual(values = c("#005555","#cccccc","firebrick"))+
  coord_flip()



#==================================================================================
#============================SECTION 10: Old models and analyses
#==================================================================================
#-----------Pelagic Broadcast Refugia Models-------------------------------------------------------------------------
mdata=subset(nssn$obs,!is.na(nssn$obs$rPela))
mdata=mdata[,c(16,21,23,24,32)]
mdata=as.data.frame(mdata)
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)
##barrier & size are correlated (0.46), size and pisc are somewhat correlated (0.39)--drop size


#full model -- max likelihood

ssn_mod <- ssn_glm(
  formula =   rPela~ scale(temp)+scale(barrier)+scale(length)+scale(pisc),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)
summary(ssn_mod)

#null model -- max likelihood
ssn_null <- ssn_glm(
  formula =   rPela~ 1,
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)

glances(ssn_mod,ssn_null)
loocv(ssn_mod)#0.445
loocv(ssn_null)#0.448

#VAR IMP
coeff=as.data.frame(ssn_mod$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ssn_mod$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
  ggtitle(label="Pelagic Spawner Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")

#-----------Periodic Strategists Refugia Models-------------------------------------------------------------------------
mdata=subset(nssn$obs,!is.na(nssn$obs$rPeri))
mdata=mdata[,c(16,21,23,24,32)]
mdata=as.data.frame(mdata)
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)

#full model -- max likelihood

ssn_mod <- ssn_glm(
  formula =   rPeri~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)
summary(ssn_mod)

#null model -- max likelihood
ssn_null <- ssn_glm(
  formula =   rPeri~ 1,
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  nugget_type = "none"
)

glances(ssn_mod,ssn_null)
loocv(ssn_mod)#0.392
loocv(ssn_null)#0.388

#Just Significant Variables bc full model worse than Null
ssn_sig <- ssn_glm(
  formula =   rPeri~ scale(size)+scale(length)+scale(pisc),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  nugget_type = "none"
)

loocv(ssn_sig)#0.382

#VAR IMP
coeff=as.data.frame(ssn_mod$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ssn_mod$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size*"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length*"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores*"
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
  ggtitle(label="Periodic Strategist Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")



#-----------Small Montane Species Refugia Models-------------------------------------------------------------------------
mdata=subset(nssn$obs,!is.na(nssn$obs$rMont))
mdata=mdata[,c(16,21,23,24,32)]
mdata=as.data.frame(mdata)
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)
#length & size correlated (0.58) #remove length

disp_init <- dispersion_initial(family = "beta", dispersion = 1, known = "dispersion")
#full model -- max likelihood
ssn_mod <- ssn_glm(
  formula =   rMont~ scale(temp)+scale(size)+scale(barrier)+scale(pisc),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  nugget_type = "none"
)
summary(ssn_mod)

#null model -- max likelihood
ssn_null <- ssn_glm(
  formula =   rMont~ 1,
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  nugget_type = "none"
)

glances(ssn_mod,ssn_null)
loocv(ssn_mod)#0.438
loocv(ssn_null)#0.464

#VAR IMP
coeff=as.data.frame(ssn_mod$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ssn_mod$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size*"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
  ggtitle(label="Small Montane Species Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")







###Limit to native species with over 30 sites present (if ever present at a site)
sp30=c("FHMN","LNDC","CRCH","LNSU","WSU","GSUN","SDSH","CARP","LL","RB","BLBH","LKCH","BRMN","PLMN","SHRH","FHCH","SCAT","RCSU","EB","PLSU","CCAT","SPDC","MTSU","COLCOT","RSSH","FMSU")
nssn$obs$FHMN=as.numeric(nssn$obs$FHMN)
nssn$obs$LNDC=as.numeric(nssn$obs$LNDC)
nssn$obs$CRCH=as.numeric(nssn$obs$CRCH)
nssn$obs$LNSU=as.numeric(nssn$obs$LNSU)
nssn$obs$WSU=as.numeric(nssn$obs$WSU)
nssn$obs$GSUN=as.numeric(nssn$obs$GSUN)
nssn$obs$SDSH=as.numeric(nssn$obs$SDSH)
nssn$obs$CARP=as.numeric(nssn$obs$CARP)
nssn$obs$LL=as.numeric(nssn$obs$LL)
nssn$obs$RB=as.numeric(nssn$obs$RB)
nssn$obs$BLBH=as.numeric(nssn$obs$BLBH)
nssn$obs$LKCH=as.numeric(nssn$obs$LKCH)
nssn$obs$BRMN=as.numeric(nssn$obs$BRMN)
nssn$obs$PLMN=as.numeric(nssn$obs$PLMN)
nssn$obs$SHRH=as.numeric(nssn$obs$SHRH)
nssn$obs$FHCH=as.numeric(nssn$obs$FHCH)
nssn$obs$SCAT=as.numeric(nssn$obs$SCAT)
nssn$obs$RCSU=as.numeric(nssn$obs$RCSU)
nssn$obs$EB=as.numeric(nssn$obs$EB)
nssn$obs$PLSU=as.numeric(nssn$obs$PLSU)
nssn$obs$CCAT=as.numeric(nssn$obs$CCAT)
nssn$obs$SPDC=as.numeric(nssn$obs$SPDC)
nssn$obs$MTSU=as.numeric(nssn$obs$MTSU)
nssn$obs$COLCOT=as.numeric(nssn$obs$COLCOT)
nssn$obs$RSSH=as.numeric(nssn$obs$RSSH)
nssn$obs$FMSU=as.numeric(nssn$obs$FMSU)



net30=subset(NET,NET$Species%in%sp30)
net30 = subset(net30,net30$net<0)

sp30ex = subset(net30, net30$Species!="RSSH" & net30$Species!="EB" &net30$Species!="BLBH")

sp30ex=as.list(sp30ex$Species)
sp30ex



#BRMN
BRMN.global <- ssn_glm(
  formula =   BRMN~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(BRMN.global)
BRMN.null <- ssn_glm(
  formula =   BRMN~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
BRMN.sig <- ssn_glm(
  formula =   BRMN~ scale(barrier)+scale(length),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(BRMN.sig)
loocv(BRMN.global)#0.423
loocv(BRMN.null)#0.341
loocv(BRMN.sig)#0.423 --null model is best
plot(BRMN~temp,data = nssn$obs, main="BRMN Refugia vs. Temperature")
plot(BRMN~size,data = nssn$obs, main="BRMN Refugia vs. Stream Size")
plot(jitter(BRMN,.2)~barrier, data = nssn$obs, main="BRMN Refugia vs. Barriers (Jittered)")
plot(BRMN~length,data = nssn$obs, main="BRMN Refugia vs. Fragment Length")
plot(jitter(BRMN,.2)~pisc,data = nssn$obs, main="BRMN Refugia vs. Piscivores (Jittered)")

#MOTCOT
nssn$obs$MOTCOT = NA
nssn$obs$MOTCOT=0
nssn$obs$MOTCOT[which(nssn$obs$COLCOT==1 | nssn$obs$RMCOT==1)]=1
MOTCOT.global <- ssn_glm(
  formula =   MOTCOT~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(MOTCOT.global)

MOTCOT.null <- ssn_glm(
  formula =   MOTCOT~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
loocv(MOTCOT.global)#0.150
loocv(MOTCOT.null)#0.156
plot(MOTCOT~temp,data = nssn$obs, main="MOTCOT Refugia vs. Temperature")
plot(MOTCOT~size,data = nssn$obs, main="MOTCOT Refugia vs. Stream Size")
plot(jitter(MOTCOT,.2)~barrier, data = nssn$obs, main="MOTCOT Refugia vs. Barriers (Jittered)")
plot(MOTCOT~length,data = nssn$obs, main="MOTCOT Refugia vs. Fragment Length")
plot(jitter(MOTCOT,.2)~pisc,data = nssn$obs, main="MOTCOT Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(MOTCOT.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="MOTCOT.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length*"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores*"
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
  ggtitle(label="Mottled Sculpin Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")




#FHCH
FHCH.global <- ssn_glm(
  formula =   FHCH~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(FHCH.global)

FHCH.null <- ssn_glm(
  formula =   FHCH~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(FHCH.global)#0.334
loocv(FHCH.null)#0.345
plot(FHCH~temp,data = nssn$obs, main="FHCH Refugia vs. Temperature")
plot(FHCH~size,data = nssn$obs, main="FHCH Refugia vs. Stream Size")
plot(jitter(FHCH,.2)~barrier, data = nssn$obs, main="FHCH Refugia vs. Barriers (Jittered)")
plot(FHCH~length,data = nssn$obs, main="FHCH Refugia vs. Fragment Length")
plot(jitter(FHCH,.2)~pisc,data = nssn$obs, main="FHCH Refugia vs. Piscivores (Jittered)")


#VAR IMP
coeff=as.data.frame(FHCH.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="FHCH.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers*"
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
  ggtitle(label="Flathead Chub Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")



#FHMN
FHMN.global <- ssn_glm(
  formula =   FHMN~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(FHMN.global)
FHMN.null <- ssn_glm(
  formula =   FHMN~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
FHMN.sig <- ssn_glm(
  formula =   FHMN~ scale(size)+scale(length),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(FHMN.global)#0.466
loocv(FHMN.null)#0.418 - null is best
loocv(FHMN.sig)#0.466
plot(FHMN~temp,data = nssn$obs, main="FHMN Refugia vs. Temperature")
plot(FHMN~size,data = nssn$obs, main="FHMN Refugia vs. Stream Size")
plot(jitter(FHMN,.2)~barrier, data = nssn$obs, main="FHMN Refugia vs. Barriers (Jittered)")
plot(FHMN~length,data = nssn$obs, main="FHMN Refugia vs. Fragment Length")
plot(jitter(FHMN,.2)~pisc,data = nssn$obs, main="FHMN Refugia vs. Piscivores (Jittered)")


#FMSU
FMSU.global <- ssn_glm(
  formula =   FMSU~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(FMSU.global)

FMSU.null <- ssn_glm(
  formula =   FMSU~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

FMSU.sig <- ssn_glm(
  formula =   FMSU~ scale(length),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
loocv(FMSU.global) #0.450
loocv(FMSU.null) #0.365
loocv(FMSU.sig) #0.410

plot(FMSU~temp,data = nssn$obs, main="FMSU Refugia vs. Temperature")
plot(FMSU~size,data = nssn$obs, main="FMSU Refugia vs. Stream Size")
plot(jitter(FMSU,.2)~barrier, data = nssn$obs, main="FMSU Refugia vs. Barriers (Jittered)")
plot(FMSU~length,data = nssn$obs, main="FMSU Refugia vs. Fragment Length")
plot(jitter(FMSU,.2)~pisc,data = nssn$obs, main="FMSU Refugia vs. Piscivores (Jittered)")



#LKCH
LKCH.global <- ssn_glm(
  formula =   LKCH~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(LKCH.global)

LKCH.null <- ssn_glm(
  formula =   LKCH~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
loocv(LKCH.global) #0.465
loocv(LKCH.null) #0.499

print(loocv_mod$RMSPE)
plot(LKCH~temp,data = nssn$obs, main="LKCH Refugia vs. Temperature")
plot(LKCH~size,data = nssn$obs, main="LKCH Refugia vs. Stream Size")
plot(jitter(LKCH,.2)~barrier, data = nssn$obs, main="LKCH Refugia vs. Barriers (Jittered)")
plot(LKCH~length,data = nssn$obs, main="LKCH Refugia vs. Fragment Length")
plot(jitter(LKCH,.2)~pisc,data = nssn$obs, main="LKCH Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(LKCH.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="LKCH.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores*"
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
  ggtitle(label="Lake Chub Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")


#LNDC
LNDC.global <- ssn_glm(
  formula =   LNDC~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(LNDC.global)

LNDC.null <- ssn_glm(
  formula =   LNDC~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

LNDC.sig <- ssn_glm(
  formula =   LNDC~ scale(temp)+scale(size),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(LNDC.global)#0.435
loocv(LNDC.null)#0.402
loocv(LNDC.sig)#0.435
plot(LNDC~temp,data = nssn$obs, main="LNDC Refugia vs. Temperature")
plot(LNDC~size,data = nssn$obs, main="LNDC Refugia vs. Stream Size")
plot(jitter(LNDC,.2)~barrier, data = nssn$obs, main="LNDC Refugia vs. Barriers (Jittered)")
plot(LNDC~length,data = nssn$obs, main="LNDC Refugia vs. Fragment Length")
plot(jitter(LNDC,.2)~pisc,data = nssn$obs, main="LNDC Refugia vs. Piscivores (Jittered)")


#LNSU
LNSU.global <- ssn_glm(
  formula =   LNSU~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(LNSU.global)

LNSU.null <- ssn_glm(
  formula =   LNSU~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

LNSU.sig <- ssn_glm(
  formula =   LNSU~ scale(temp)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(LNSU.global)#0.467
loocv(LNSU.null)#0.420
loocv(LNSU.sig) #0.473

plot(LNSU~temp,data = nssn$obs, main="LNSU Refugia vs. Temperature")
plot(LNSU~size,data = nssn$obs, main="LNSU Refugia vs. Stream Size")
plot(jitter(LNSU,.2)~barrier, data = nssn$obs, main="LNSU Refugia vs. Barriers (Jittered)")
plot(LNSU~length,data = nssn$obs, main="LNSU Refugia vs. Fragment Length")
plot(jitter(LNSU,.2)~pisc,data = nssn$obs, main="LNSU Refugia vs. Piscivores (Jittered)")


#MTSU
MTSU.global <- ssn_glm(
  formula =   MTSU~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(MTSU.global)
MTSU.null <- ssn_glm(
  formula =   MTSU~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
MTSU.sig <- ssn_glm(
  formula =   MTSU~ scale(size)+scale(barrier),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
loocv(MTSU.global)#0.429
loocv(MTSU.null)#0.391
loocv(MTSU.sig)#0.445


plot(MTSU~temp,data = nssn$obs, main="MTSU Refugia vs. Temperature")
plot(MTSU~size,data = nssn$obs, main="MTSU Refugia vs. Stream Size")
plot(jitter(MTSU,.2)~barrier, data = nssn$obs, main="MTSU Refugia vs. Barriers (Jittered)")
plot(MTSU~length,data = nssn$obs, main="MTSU Refugia vs. Fragment Length")
plot(jitter(MTSU,.2)~pisc,data = nssn$obs, main="MTSU Refugia vs. Piscivores (Jittered)")


#PLMN
PLMN.global <- ssn_glm(
  formula =   PLMN~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(PLMN.global)

PLMN.null <- ssn_glm(
  formula =   PLMN~ 1, 
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(PLMN.global)#0.372
loocv(PLMN.null)#0.414
plot(PLMN~temp,data = nssn$obs, main="PLMN Refugia vs. Temperature")
plot(PLMN~size,data = nssn$obs, main="PLMN Refugia vs. Stream Size")
plot(jitter(PLMN,.2)~barrier, data = nssn$obs, main="PLMN Refugia vs. Barriers (Jittered)")
plot(PLMN~length,data = nssn$obs, main="PLMN Refugia vs. Fragment Length")
plot(jitter(PLMN,.2)~pisc,data = nssn$obs, main="PLMN Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(PLMN.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="PLMN.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
  ggtitle(label="Plains Minnow Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")


#PLSU
PLSU.global <- ssn_glm(
  formula =   PLSU~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(PLSU.global)

PLSU.null <- ssn_glm(
  formula =   PLSU~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(PLSU.global)#0.481
loocv(PLSU.null)#0.5

plot(PLSU~temp,data = nssn$obs, main="PLSU Refugia vs. Temperature")
plot(PLSU~size,data = nssn$obs, main="PLSU Refugia vs. Stream Size")
plot(jitter(PLSU,.2)~barrier, data = nssn$obs, main="PLSU Refugia vs. Barriers (Jittered)")
plot(PLSU~length,data = nssn$obs, main="PLSU Refugia vs. Fragment Length")
plot(jitter(PLSU,.2)~pisc,data = nssn$obs, main="PLSU Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(PLSU.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="PLSU.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
  ggtitle(label="Plains Sucker Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")


#SPDC
SPDC.global <- ssn_glm(
  formula =   SPDC~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(SPDC.global)

SPDC.null <- ssn_glm(
  formula =   SPDC~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")


loocv(SPDC.global)#0.438
loocv(SPDC.null) #0.468
plot(SPDC~temp,data = nssn$obs, main="SPDC Refugia vs. Temperature")
plot(SPDC~size,data = nssn$obs, main="SPDC Refugia vs. Stream Size")
plot(jitter(SPDC,.2)~barrier, data = nssn$obs, main="SPDC Refugia vs. Barriers (Jittered)")
plot(SPDC~length,data = nssn$obs, main="SPDC Refugia vs. Fragment Length")
plot(jitter(SPDC,.2)~pisc,data = nssn$obs, main="SPDC Refugia vs. Piscivores (Jittered)")


#VAR IMP
coeff=as.data.frame(SPDC.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="SPDC.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size*"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
  ggtitle(label="Speckled Dace Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")



#WSU
WSU.global <- ssn_glm(
  formula =   WSU~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(WSU.global)

WSU.null <- ssn_glm(
  formula =   WSU~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")


loocv(WSU.global)#0.376
loocv(WSU.null) #0.378

plot(WSU~temp,data = nssn$obs, main="WSU Refugia vs. Temperature")
plot(WSU~size,data = nssn$obs, main="WSU Refugia vs. Stream Size")
plot(jitter(WSU,.2)~barrier, data = nssn$obs, main="WSU Refugia vs. Barriers (Jittered)")
plot(WSU~length,data = nssn$obs, main="WSU Refugia vs. Fragment Length")
plot(jitter(WSU,.2)~pisc,data = nssn$obs, main="WSU Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(WSU.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="WSU.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size*"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers*"
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
  ggtitle(label="White Sucker Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")


#BLBH
BLBH.global <- ssn_glm(
  formula =   BLBH~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(BLBH.global)

BLBH.null <- ssn_glm(
  formula =   BLBH~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(BLBH.global)#0.426
loocv(BLBH.null) #0.497

plot(BLBH~temp,data = nssn$obs, main="BLBH Refugia vs. Temperature")
plot(BLBH~size,data = nssn$obs, main="BLBH Refugia vs. Stream Size")
plot(jitter(BLBH,.2)~barrier, data = nssn$obs, main="BLBH Refugia vs. Barriers (Jittered)")
plot(BLBH~length,data = nssn$obs, main="BLBH Refugia vs. Fragment Length")
plot(jitter(BLBH,.2)~pisc,data = nssn$obs, main="BLBH Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(BLBH.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="BLBH.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers*"
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
  ggtitle(label="Black Bullhead Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")



########INTRODUCED SPECIES
#EB
EB.global <- ssn_glm(
  formula =   EB~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(EB.global)

EB.null <- ssn_glm(
  formula =   EB~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")

loocv(EB.global)#0.408
loocv(EB.null)#0.448

plot(EB~temp,data = nssn$obs, main="EB Refugia vs. Temperature")
plot(EB~size,data = nssn$obs, main="EB Refugia vs. Stream Size")
plot(jitter(EB,.2)~barrier, data = nssn$obs, main="EB Refugia vs. Barriers (Jittered)")
plot(EB~length,data = nssn$obs, main="EB Refugia vs. Fragment Length")
plot(jitter(EB,.2)~pisc,data = nssn$obs, main="EB Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(EB.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="EB.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
  ggtitle(label="Brook Trout Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")



#RSSH
RSSH.global <- ssn_glm(
  formula =   RSSH~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
summary(RSSH.global)

RSSH.null <- ssn_glm(
  formula =   RSSH~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")


loocv(RSSH.global)#0.358
loocv(RSSH.null)#0.385

plot(RSSH~temp,data = nssn$obs, main="RSSH Refugia vs. Temperature")
plot(RSSH~size,data = nssn$obs, main="RSSH Refugia vs. Stream Size")
plot(jitter(RSSH,.2)~barrier, data = nssn$obs, main="RSSH Refugia vs. Barriers (Jittered)")
plot(RSSH~length,data = nssn$obs, main="RSSH Refugia vs. Fragment Length")
plot(jitter(RSSH,.2)~pisc,data = nssn$obs, main="RSSH Refugia vs. Piscivores (Jittered)")

#VAR IMP
coeff=as.data.frame(RSSH.global$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="RSSH.global$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
  ggtitle(label="Redside Shiner Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")

#-----------ANOVAs for basin-------
pu=read.csv("units_ecoregions.csv")
pu=pu%>%rename("RepetID"="RepeatID")
basinSSN=nssn
basinSSN$obs=basinSSN$obs%>%rename("RepetID"="RepeatID")
basinSSN$obs=left_join(basinSSN$obs, pu,by="RepetID")
obsrelcol=obsrelcol%>%rename("RepetID"="RepeatID")
basinSSN$obs=left_join(basinSSN$obs,obsrelcol, by="RepetID")
basinSSN$obs$pGlac[which(basinSSN$obs$pGlac==1)]=0.999
basinSSN$obs$pGlac[which(basinSSN$obs$pGlac==0)]=0.001

#basinSSN$obs$Pubasin=as.factor(basinSSN$obs$Pubasin)
class(basinSSN$obs$Pubasin)
quantile(basinSSN$obs$size)
basinSSN$obs$sizeCat=NA
basinSSN$obs$sizeCat[which(nssn$obs$size<4.06)]="s"
basinSSN$obs$sizeCat[which(nssn$obs$size>=4.06 & nssn$obs$size<14.88)]="m"
basinSSN$obs$sizeCat[which(nssn$obs$size>=14.88 & nssn$obs$size<68.99)]="l"
basinSSN$obs$sizeCat[which(nssn$obs$size>=68.99)]="xl"

basinSSN$obs$cheybasin=NA
basinSSN$obs$cheybasin[which(basinSSN$obs$Pubasin=="CHEY")]="zCHEY"
basinSSN$obs$cheybasin[which(basinSSN$obs$Pubasin!="CHEY")]="Nope"
basinSSN$obs$cheybasin

#pNat by Basin
ssn_mod <- ssn_glm(
  formula =   pNat~Pubasin,
  family = "beta",
  ssn.object = basinSSN,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat)+as.factor(sizeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  estmethod = "reml")

ssn_mod_sum=summary(ssn_mod) #UPMO sig lower, Chey sig higher
ssn_mod_sum
anova(ssn_mod)
tidy(anova(ssn_mod))

ssn_mod <- ssn_glm(
  formula =   pNat~cheybasin,
  family = "beta",
  ssn.object = basinSSN,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat)+as.factor(sizeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  estmethod = "reml")
ssn_mod_sum=summary(ssn_mod) #UPMO sig lower, Chey sig higher
ssn_mod_sum

basinSSN$obs%>%filter(!is.na(pNat))%>%
  group_by(Pubasin)%>%
  summarise(mpNat = mean(pNat), mpTurn = mean(Turnover), mpGlacial =mean(pGlac))


#Turnover by Basin
ssn_mod <- ssn_glm(
  formula =   Turnover~Pubasin,
  family = "beta",
  ssn.object = basinSSN,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat)+as.factor(sizeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  estmethod = "reml")

ssn_mod_sum=summary(ssn_mod)
ssn_mod_sum

ssn_mod <- ssn_glm(
  formula =   Turnover~cheybasin,
  family = "beta",
  ssn.object = basinSSN,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat)+as.factor(sizeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  estmethod = "reml")
ssn_mod_sum=summary(ssn_mod)
ssn_mod_sum

basinSSN$obs%>%filter(!is.na(Turnover))%>%
  group_by(Pubasin)%>%
  summarise(mTurn = mean(Turnover))




#pGlac by Basin
ssn_mod <- ssn_glm(
  formula =   pGlac~Pubasin,
  family = "beta",
  ssn.object = basinSSN,
  #tailup_type = "exponential",
  taildown_type = "exponential",
  #euclid_type = "exponential",
  random = ~as.factor(YrangeCat)+as.factor(sizeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  estmethod = "reml",
  nugget_type = "none")

ssn_mod_sum=summary(ssn_mod)
ssn_mod_sum

ssn_mod <- ssn_glm(
  formula =   pGlac~cheybasin,
  family = "beta",
  ssn.object = basinSSN,
  #tailup_type = "exponential",
  taildown_type = "exponential",
  #euclid_type = "exponential",
  random = ~as.factor(YrangeCat)+as.factor(sizeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  estmethod = "reml",
  nugget_type = "none")
ssn_mod_sum=summary(ssn_mod) #UPMO sig lower, Chey sig higher
ssn_mod_sum

basinSSN$obs%>%filter(!is.na(pGlac))%>%
  group_by(Pubasin)%>%
  summarise(mpGlac = mean(pGlac))




basinSSN$obs$Pubasin[which(basinSSN$obs$Pubasin=="CHEY")]="zCHEY"

#pNat by Basin
ssn_mod <- ssn_glm(
  formula =   pNat~Pubasin,
  family = "beta",
  ssn.object = basinSSN,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat)+as.factor(sizeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  estmethod = "reml")

ssn_mod_sum=summary(ssn_mod) #UPMO sig lower, Chey sig higher
ssn_mod_sum

#Turnover by Basin
ssn_mod <- ssn_glm(
  formula =   Turnover~Pubasin,
  family = "beta",
  ssn.object = basinSSN,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat)+as.factor(sizeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  estmethod = "reml")

ssn_mod_sum=summary(ssn_mod)
ssn_mod_sum

#pGlac by Basin
ssn_mod <- ssn_glm(
  formula =   pGlac~Pubasin,
  family = "beta",
  ssn.object = basinSSN,
  #tailup_type = "exponential",
  taildown_type = "exponential",
  #euclid_type = "exponential",
  random = ~as.factor(YrangeCat)+as.factor(sizeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  estmethod = "reml",
  nugget_type = "none")

ssn_mod_sum=summary(ssn_mod)
ssn_mod_sum



#---------- Residuals Graph---------------
nssn$obs$pNat[which(nssn$obs$pNat==1)]=0.999
nssn$obs$pNat[which(nssn$obs$pNat==0)]=0.001

#significant only -- max likelihood
temp.mod <- ssn_glm(
  formula =   pNat~ scale(temp),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  #nugget_type = "none"
)

temp.resid=as.data.frame(temp.mod$residuals)
temp.resid <- tibble::rownames_to_column(temp.resid, "RepetID")
temp.resid$RepetID=as.numeric(temp.resid$RepetID)
temp.resid$res2=NA
temp.resid$res2=abs(temp.resid$response)
temp.resid=temp.resid[,c(1,6)]
obviate=nssn$obs
obviate=left_join(obviate,temp.resid,by="RepetID")

obviate%>%
  ggplot(aes(x=length, y=res2))+
  geom_point()+
  geom_smooth(method = "lm") #basically no relationship

