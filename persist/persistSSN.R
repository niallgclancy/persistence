################################################################
# R script for Clancy Dissertation Chapter on Empirical Characteristics of Refugia


#######LIMITING ANALYSIS TO PERSISTENCE/EXTIRPATION SITES (COLONIZATION POINTS REMOVED) 11-19-2024


################################################################
library(SSNbler)
library(SSN2)
library(caret)
library(tidyverse)
library(sf)





#======================================================
#============================SECTION 1: SSN Preparation
#======================================================

#------------------Importing Required Data---------------------

## import the streams, observation sites
streams <- st_read("NSI_fix6.shp")
obs <- st_read("RefugiaSites2.shp")

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




#========================================================================
#============================SECTION 2: Models for Community-Wide Metrics
#========================================================================
#----------------Identify predictors with multicolinearity---------------
library(PerformanceAnalytics)
mdata=obs[,c(14,19,20,21,22,23,30,31)]
mdata=as.data.frame(mdata)
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)
#None of the selected predictors (for data exploration) are highly correlated (max r = 0.32)

# ---- Fit spatial stream-network model -------------------------------
#------Beta regression can't use 0 or 1...convert to near 0 and 1------
nssn2=nssn

#Full Community Persistence
nssn2$obs$persComm=as.numeric(nssn2$obs$persComm)
nssn2$obs$persComm[which(nssn2$obs$persComm==0)]=0.000000001
nssn2$obs$persComm[which(nssn2$obs$persComm==1)]=0.999999999
anyNA(nssn2$obs$persComm)

#Native Sp. Persistence
nssn2$obs$persNative=as.numeric(nssn2$obs$persNative)
nssn2$obs$persNative=as.numeric(nssn2$obs$persNative)
nssn2$obs$persNative[which(nssn2$obs$persNative==0)]=0.000000001
nssn2$obs$persNative[which(nssn2$obs$persNative==1)]=0.999999999

#Early Postglacial Colonists Persistence
nssn2$obs$persGlacE=as.numeric(nssn2$obs$persGlacE)
nssn2$obs$persGlacE[which(nssn2$obs$persGlacE==0)]=0.000000001
nssn2$obs$persGlacE[which(nssn2$obs$persGlacE==1)]=0.999999999

#Late Postglacial Colonists Persistence
nssn2$obs$persGlacL=as.numeric(nssn2$obs$persGlacL)
nssn2$obs$persGlacL[which(nssn2$obs$persGlacL==0)]=0.000000001
nssn2$obs$persGlacL[which(nssn2$obs$persGlacL==1)]=0.999999999

################################
###############################Beta Regression for proportion metrics not working 11-19-2024
###############################
##Scale covariates to determine variable importance
ssn_mod <- ssn_glm(
  formula =   persComm~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  #spcov_initial=1,
  #dispersion_initial = 1,
  spcov_type = "exponential",
  additive = "afvArea", estmethod = "reml"
)
############################None of this working, attempted messing with spcov_initial, dispersion_initial, and spcov_type...none of which worked





## ---- Fit spatial (not stream-network) model -------------------------------
######Attempt with spmodel using euclidean only
library(spmodel)
# ---- Fit model for Community Persistence ----------------------------
obs2=obs
obs2$persComm[which(obs2$persComm==0)]=0.0000001
obs2$persComm[which(obs2$persComm==1)]=0.9999999
obs2$persNative[which(obs2$persNative==0)]=0.0000001
obs2$persNative[which(obs2$persNative==1)]=0.9999999

obsMO=subset(obs2, obs2$PUname!="Upper Green")
obsGR=subset(obs2, obs2$PUname=="Upper Green")

#
#------------------Missouri Basin Models-------------------------------------------------
#
test.mod<- spglm(
  formula = persComm~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml")
summary(test.mod)
loocv(test.mod) #RMSPE=0.254

####Compare to intercept only (null)   #NOT WORKING EITHER...same error as with SSN2
test.mod.null<- spglm(
  formula = persComm~ 1,
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml")
summary(test.mod.null)
loocv(test.mod.null) #0.264

#compare to only significant variables
test.mod.sig<- spglm(
  formula = persComm~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) #RMSPE=0.277


#Non-Spatial Model
test.mod.nonspatial<- glm(
  formula = persComm~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds)+(1|Yrange),
  family = "binomial",
  data = obsMO)

test.mod.nonspatial<- spglm(
  formula = persComm~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml", spcov_type = "none")
summary(test.mod.nonspatial)


glances(test.mod,test.mod.null,
        test.mod.sig,test.mod.nonspatial)

#
#------------------Green Basin Models-------------------------------------------------
#
#
test.mod<- spglm(
  formula = persComm~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml")
summary(test.mod)
loocv(test.mod) 

####Compare to intercept only (null)   -no convergence
test.mod.null<- spglm(
  formula = scale(persComm)~ 1,
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml")
summary(test.mod.null)
loocv(test.mod.null) 

#compare to only significant variables--no convergence
test.mod.sig<- spglm(
  formula = persComm~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) #RMSPE=0.277

#Non-Spatial Model
test.mod.nonspatial<- spglm(
  formula = persComm~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml", spcov_type = "none")
summary(test.mod.nonspatial)
loocv(test.mod.nonspatial)

###NS Intercept
test.mod.nonspatial.null<- spglm(
  formula = persComm~ 1,
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml", spcov_type = "none")
summary(test.mod.nonspatial.null)
loocv(test.mod.nonspatial.null)

glances(test.mod,test.mod.nonspatial,test.mod.nonspatial.null)



# ---- Fit model for Native Species Persistence----------------------------
#------------------Missouri Basin Models-------------------------------------------------
#
test.mod<- spglm(
  formula = persNative~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml")
summary(test.mod)
loocv(test.mod) #RMSPE=0.254

####Compare to intercept only (null)   #NOT WORKING EITHER...same error as with SSN2
test.mod.null<- spglm(
  formula = persNative~ 1,
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml")
summary(test.mod.null)
loocv(test.mod.null) #0.264

#compare to only significant variables
test.mod.sig<- spglm(
  formula = persNative~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(Length_km),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) #RMSPE=0.277


#Non-Spatial Model
test.mod.nonspatial<- glm(
  formula = persNative~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds)+(1|Yrange),
  family = "binomial",
  data = obsMO)

test.mod.nonspatial<- spglm(
  formula = persNative~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml", spcov_type = "none")
summary(test.mod.nonspatial)
loocv(test.mod.nonspatial)

glances(test.mod,test.mod.null,
        test.mod.sig,test.mod.nonspatial)

#
#------------------Green Basin Models-------------------------------------------------
#
#
test.mod<- spglm(
  formula = persNative~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml")
summary(test.mod)
loocv(test.mod) 

####Compare to intercept only (null)   -no convergence
test.mod.null<- spglm(
  formula = persNative~ 1,
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml")
summary(test.mod.null)
loocv(test.mod.null) 

#compare to only significant variables--no convergence
test.mod.sig<- spglm(
  formula = persNative~ scale(F_MAUG_HIS),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) #RMSPE=0.277

#Non-Spatial Model
test.mod.nonspatial<- spglm(
  formula = persNative~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml", spcov_type = "none")
summary(test.mod.nonspatial)
loocv(test.mod.nonspatial)

glances(test.mod,test.mod.null,test.mod.sig,test.mod.nonspatial)






# ---- Fit model for Glacial-Relict Persistence----------------------------
obs2$persGlacE[which(obs2$persGlacE==0)]=0.000000001
obs2$persGlacE[which(obs2$persGlacE==1)]=0.999999999

relict.mod<- spglm(
  formula = persGlacE~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2
)
summary(relict.mod)
loocv(relict.mod)

# ---- Compare to model fit to non-relicts
obs2$persGlacL[which(obs2$persGlacL==0)]=0.000000001
obs2$persGlacL[which(obs2$persGlacL==1)]=0.999999999

nonrelict.mod<- spglm(
  formula = persGlacE~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2
)
summary(nonrelict.mod)
loocv(nonrelict.mod)



# ---- Fit model for Glacial-Relict Persistence----------------------------



#-------------------------------Mean persistence, uncorrected for autocorrelation or colonization------------------------------- 
obs3=obs2
obs3=as.data.frame(obs3)
obs3%>%filter(!is.na(persGlacE))%>%summarise(avg=mean(persGlacE), sd(persGlacE)) #58% (sd=34%)
obs3%>%filter(!is.na(persGlacL))%>%summarise(avg=mean(persGlacL), sd(persGlacL)) #53% (sd=43%)
obs3%>%filter(!is.na(persComm))%>%summarise(avg=mean(persComm), sd(persComm)) #56% (sd=30%)
obs3%>%filter(!is.na(persNative))%>%summarise(avg=mean(persNative),sd(persNative)) #56% (sd=32%)


#=====================================================================
#============================SECTION 3: Examine Net Change for Species
#=====================================================================

#---------------------Net Species Change------------------------------------------
NET=read.csv("PersistenceMetrics.csv")
NET=NET%>%pivot_longer(cols = c(28:121),names_to = "Species",values_to = "change")

NET=subset(NET,!is.na(NET$change))
NET2=NET%>%group_by(Species)%>%summarise(net=sum(change))
NET3=NET%>%filter(change!=-1)%>%group_by(Species)%>%summarise(LateNum=length(change))
NET=NET%>%filter(change!=1)%>%group_by(Species)%>%summarise(EarlyNum=length(change))
NET=left_join(NET2,NET,by="Species")
NET=left_join(NET,NET3,by="Species")
NET$EarlyNum[which(is.na(NET$EarlyNum))]=0
NET$LateNum[which(is.na(NET$LateNum))]=0
NET$propChange=NA
NET$propChange=NET$LateNum/NET$EarlyNum
NETsub=subset(NET,NET$EarlyNum>=10|NET$LateNum>=10)
NETsub=subset(NETsub,NETsub$Species!="ONC" & NETsub$Species!="FMxWSU")

traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
NETsub=left_join(NETsub,traits,by="Species")
NETsub$CommonName[which(NETsub$CommonName=="Northern Redbelly Dace")]="Nor. Redbelly Dace"
NETsub$CommonName[which(NETsub$CommonName=="Rocky Mountain Cutthroat Trout")]="R.M. Cutthroat Trout"

declining=NETsub%>%filter(propChange<1)%>%
  ggplot(aes(x=reorder(CommonName,propChange),y=propChange))+
  geom_segment( aes(x=CommonName, xend=CommonName, y=1, yend=propChange), color="grey") +
  geom_point(aes(colour = as.factor(Status)), size=3)+
  scale_color_manual(values = c("#ff9999","#005555"))+
  scale_x_discrete(position = "top")+
  theme_light()+
  ylab(label = "Net Change in Proportional Occurrence")+
  xlab(label="")+
  ggtitle(label = "Declining Species")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0,face = 2),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size=12),
        axis.ticks.x = element_blank())
declining
ggsave(filename="DecliningSp.tiff",dpi = 400, width = 8, height = 6, units = "in")
increasing=NETsub%>%filter(propChange>=1)%>%
  ggplot(aes(x=reorder(CommonName,-propChange),y=propChange))+
  geom_segment( aes(x=CommonName, xend=CommonName, y=1, yend=propChange), color="grey") +
  geom_point(aes(colour = as.factor(Status)), size=3)+
  scale_color_manual(values = c("#ff9999","#cccccc","#005555"))+
  theme_light()+
  ylab(label = "Net Change in Proportional Occurrence")+
  xlab(label="")+
  ggtitle(label = "Increasing or Stable Species")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face = 2),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size=12),
        axis.ticks.x = element_blank())
increasing
ggsave(filename="IncreasingSp.tiff",dpi = 400, width = 8, height = 6, units = "in")




#----------------------Check against Ch1 thresholds and Elizabeth's thresholds---------------------
pers=read.csv("PersistenceMetricsFlow.csv")
pers=pers%>%pivot_longer(cols = c(28:121),names_to = "Species",values_to = "change")


#FLOW
thresh=read.csv("Thresholds.csv")
persthresh=left_join(pers,thresh, by="Species")
#for S1S30 avg
#persthresh$Savg=NA
#persthresh$Savg=persthresh$S30_2040D+persthresh$S1_93_11
#persthresh$Savg=persthresh$Savg/2
persthresh$flowdif =NA
persthresh$flowdif = persthresh$F_MAUG_HIS-persthresh$FlowTolerance
persthresh$tempdif =NA
persthresh$tempdif = persthresh$S1_93_11-persthresh$Upper_therm

persthresh$OverThresh=NA
persthresh$OverThresh[which(persthresh$flowdif>=0 | persthresh$tempdif>=0)]="Y"
persthresh$OverThresh[which(persthresh$flowdif<0 & persthresh$OverThresh!="Y")]="N"
persthresh$OverThresh[which(persthresh$tempdif<0 & persthresh$OverThresh!="Y")]="N"
persthresh$CorrectOver=NA
persthresh$CorrectOver="N"
persthresh$CorrectOver[which(persthresh$change==-1 & persthresh$OverThresh=="Y")]="Y"

cfsthresh=persthresh%>%filter(!is.na(change) & OverThresh=="Y")
cfsovercorr=cfsthresh%>%group_by(Species, CorrectOver)%>%
  summarise(n=length(GNIS_NAME))
cfsovertot=cfsthresh%>%group_by(Species)%>%
  summarise(total=length(GNIS_NAME))
cfsthresh=full_join(cfsovercorr,cfsovertot,by="Species")
cfsthresh=cfsthresh%>%pivot_wider(names_from = "CorrectOver", values_from = "n")
cfsthresh$N[which(is.na(cfsthresh$N))]=0
cfsthresh$Y[which(is.na(cfsthresh$Y))]=0
cfsthresh$prop=NA
cfsthresh$prop=cfsthresh$Y/cfsthresh$total
cfsthresh$total=NULL
cfsthresh$N=NULL
cfsthresh$Y=NULL
cfsthresh=cfsthresh%>%rename("ExtirpOverThresh" = "prop")


NETsub=left_join(NETsub,cfsthresh,by="Species")
NETsub$ExtirpGeneral=NA
NETsub$ExtirpGeneral=1-NETsub$propChange
write.csv(NETsub, "NETsub.csv")


#===========================================================================
#============================SECTION 4: Individual Species Models
#================================================================

#:::::::::::::::::::::::::::::::::::::
# Individual Species analysis###
#:::::::::::::::::::::::::::::::::::::
#Scale covariates to determine variable importance

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
  formula =   BRMN~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(BRMN.global)

loocv_mod <- loocv(BRMN.global)
print(loocv_mod$RMSPE)


#COLCOT
COLCOT.global <- ssn_glm(
  formula =   COLCOT~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(COLCOT.global)

loocv_mod <- loocv(COLCOT.global)
print(loocv_mod$RMSPE)

#FHCH
FHCH.global <- ssn_glm(
  formula =   FHCH~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(FHCH.global)

loocv_mod <- loocv(FHCH.global)
print(loocv_mod$RMSPE)


#FHMN
FHMN.global <- ssn_glm(
  formula =   FHMN~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(FHMN.global)

loocv_mod <- loocv(FHMN.global)
print(loocv_mod$RMSPE)


#FMSU
FMSU.global <- ssn_glm(
  formula =   FMSU~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(FMSU.global)

loocv_mod <- loocv(FMSU.global)
print(loocv_mod$RMSPE)


#LKCH
LKCH.global <- ssn_glm(
  formula =   LKCH~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(LKCH.global)

loocv_mod <- loocv(LKCH.global)
print(loocv_mod$RMSPE)


#LNDC
LNDC.global <- ssn_glm(
  formula =   LNDC~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(LNDC.global)

loocv_mod <- loocv(LNDC.global)
print(loocv_mod$RMSPE)



#LNSU
LNSU.global <- ssn_glm(
  formula =   LNSU~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(LNSU.global)

loocv_mod <- loocv(LNSU.global)
print(loocv_mod$RMSPE)



#MTSU
MTSU.global <- ssn_glm(
  formula =   MTSU~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(MTSU.global)

loocv_mod <- loocv(MTSU.global)
print(loocv_mod$RMSPE)


#PLMN
PLMN.global <- ssn_glm(
  formula =   PLMN~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(Length_km)+scale(nnPreds), #all sites must have a DS barrier (or not)---wouldn't run with it
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(PLMN.global)

loocv_mod <- loocv(PLMN.global)
print(loocv_mod$RMSPE)


#PLSU
PLSU.global <- ssn_glm(
  formula =   PLSU~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(PLSU.global)

loocv_mod <- loocv(PLSU.global)
print(loocv_mod$RMSPE)


#SPDC
SPDC.global <- ssn_glm(
  formula =   SPDC~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(SPDC.global)

loocv_mod <- loocv(SPDC.global)
print(loocv_mod$RMSPE)


#WSU
WSU.global <- ssn_glm(
  formula =   WSU~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(WSU.global)

loocv_mod <- loocv(WSU.global)
print(loocv_mod$RMSPE)





#-------------------Introduced Species------------------------



#BLBH
BLBH.global <- ssn_glm(
  formula =   BLBH~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(BLBH.global)

loocv_mod <- loocv(BLBH.global)
print(loocv_mod$RMSPE)


#EB
EB.global <- ssn_glm(
  formula =   EB~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(EB.global)

loocv_mod <- loocv(EB.global)
print(loocv_mod$RMSPE)


#RSSH
RSSH.global <- ssn_glm(
  formula =   RSSH~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

summary(RSSH.global)

loocv_mod <- loocv(RSSH.global)
print(loocv_mod$RMSPE)






#####Intercept RMSPE


#BRMN
BRMN.global <- ssn_glm(
  formula =   BRMN~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

loocv_mod <- loocv(BRMN.global)
print(loocv_mod$RMSPE)


#COLCOT
COLCOT.global <- ssn_glm(
  formula =   COLCOT~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")


loocv_mod <- loocv(COLCOT.global)
print(loocv_mod$RMSPE)

#FHCH
FHCH.global <- ssn_glm(
  formula =   FHCH~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")


loocv_mod <- loocv(FHCH.global)
print(loocv_mod$RMSPE)


#FHMN
FHMN.global <- ssn_glm(
  formula =   FHMN~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")


loocv_mod <- loocv(FHMN.global)
print(loocv_mod$RMSPE)


#FMSU
FMSU.global <- ssn_glm(
  formula =   FMSU~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")


loocv_mod <- loocv(FMSU.global)
print(loocv_mod$RMSPE)


#LKCH
LKCH.global <- ssn_glm(
  formula =   LKCH~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

loocv_mod <- loocv(LKCH.global)
print(loocv_mod$RMSPE)


#LNDC
LNDC.global <- ssn_glm(
  formula =   LNDC~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

loocv_mod <- loocv(LNDC.global)
print(loocv_mod$RMSPE)



#LNSU
LNSU.global <- ssn_glm(
  formula =   LNSU~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

loocv_mod <- loocv(LNSU.global)
print(loocv_mod$RMSPE)



#MTSU
MTSU.global <- ssn_glm(
  formula =   MTSU~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

loocv_mod <- loocv(MTSU.global)
print(loocv_mod$RMSPE)


#PLMN
PLMN.global <- ssn_glm(
  formula =   PLMN~ 1, #all sites must have a DS barrier (or not)---wouldn't run with it
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

loocv_mod <- loocv(PLMN.global)
print(loocv_mod$RMSPE)


#PLSU
PLSU.global <- ssn_glm(
  formula =   PLSU~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

loocv_mod <- loocv(PLSU.global)
print(loocv_mod$RMSPE)


#SPDC
SPDC.global <- ssn_glm(
  formula =   SPDC~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

loocv_mod <- loocv(SPDC.global)
print(loocv_mod$RMSPE)


#WSU
WSU.global <- ssn_glm(
  formula =   WSU~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

loocv_mod <- loocv(WSU.global)
print(loocv_mod$RMSPE)





#-------------------Introduced Species------------------------



#BLBH
BLBH.global <- ssn_glm(
  formula =   BLBH~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

loocv_mod <- loocv(BLBH.global)
print(loocv_mod$RMSPE)


#EB
EB.global <- ssn_glm(
  formula =   EB~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

loocv_mod <- loocv(EB.global)
print(loocv_mod$RMSPE)


#RSSH
RSSH.global <- ssn_glm(
  formula =   RSSH~ 1,
  family = "binomial",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(Yrange),
  additive = "afvArea", estmethod = "ml")

loocv_mod <- loocv(RSSH.global)
print(loocv_mod$RMSPE)

