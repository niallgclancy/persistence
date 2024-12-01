################################################################
# R script for Clancy Dissertation Chapter on Empirical Characteristics of Refugia


#######LIMITING ANALYSIS TO PERSISTENCE/EXTIRPATION SITES (COLONIZATION POINTS REMOVED) 11-19-2024


################################################################
library(SSNbler)
library(SSN2)
library(caret)
library(tidyverse)
library(sf)
library(spmodel)




#======================================================
#============================SECTION 1: SSN Preparation
#======================================================

#------------------Importing Required Data---------------------

## import the streams, observation sites
streams <- st_read("NSI_fix6.shp")
obs <- st_read("RefugiaSites.shp")

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
#Identify predictors with multicolinearity
library(PerformanceAnalytics)
mdata=obs[,c(14,19,20,21,22,23,30,31)]
mdata=as.data.frame(mdata)
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)
#None of the selected predictors (for data exploration) are highly correlated (max r = 0.32)

# ---- Fit spatial stream-network model -------------------------------
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

###############################Beta Regression for proportion metrics not working 11-19-2024
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





##Fit spatial (not stream-network) model
######Attempt with spmodel using euclidean only
# ---- Fit model for Community Persistence ----------------------------
obs2=obs
obs2$persComm[which(obs2$persComm==0)]=0.0000001
obs2$persComm[which(obs2$persComm==1)]=0.9999999
obs2$persNative[which(obs2$persNative==0)]=0.0000001
obs2$persNative[which(obs2$persNative==1)]=0.9999999
obs2$persGlacE[which(obs2$persGlacE==0)]=0.000000001
obs2$persGlacE[which(obs2$persGlacE==1)]=0.999999999
obs2$persGlacL[which(obs2$persGlacL==0)]=0.000000001
obs2$persGlacL[which(obs2$persGlacL==1)]=0.999999999
pu=read.csv("processingunits.csv")
pu$X=NULL
obs2=left_join(obs2,pu,by="RepeatID")

obs2%>%pivot_longer(cols=32:124, names_to = "Species", values_to = "change")%>%filter(!is.na(change))%>%summarise(mean(change))

obsMO=subset(obs2, obs2$PU!="GREEN")
obsGR=subset(obs2, obs2$PU=="GREEN")

#
#Missouri Basin Models
test.mod<- spglm(
  formula = persComm~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml")
summary(test.mod)
loocv(test.mod) #RMSPE=0.254
plot(persComm~S1_93_11,data = obsMO, main="MO Basin Community Persistence vs. Temperature")
plot(persComm~F_MAUG_HIS,data = obsMO, main="MO Basin Community Persistence vs. Stream Size")
plot(persComm~DSbarrier,data = obsMO, main="MO Basin Community Persistence vs. Barriers")
plot(persComm~Length_km,data = obsMO, main="MO Basin Community Persistence vs. Fragment Length")
plot(persComm~nnPreds,data = obsMO, main="MO Basin Community Persistence vs. Piscivores")

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
test.mod.nonspatial<- spglm(
  formula = persComm~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml", spcov_type = "none")
summary(test.mod.nonspatial)

#Compare model performance
glances(test.mod,test.mod.null,
        test.mod.sig,test.mod.nonspatial)

#
#Green Basin Models
test.mod<- spglm(
  formula = persComm~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml")
summary(test.mod)
loocv(test.mod) 
plot(persComm~S1_93_11,data = obsGR, main="GR Basin Community Persistence vs. Temperature")
plot(persComm~F_MAUG_HIS,data = obsGR, main="GR Basin Community Persistence vs. Stream Size")
plot(persComm~DSbarrier,data = obsGR, main="GR Basin Community Persistence vs. Barriers")
plot(persComm~Length_km,data = obsGR, main="GR Basin Community Persistence vs. Fragment Length")
plot(persComm~nnPreds,data = obsGR, main="GR Basin Community Persistence vs. Piscivores")
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
#Missouri Basin Models
test.mod<- spglm(
  formula = persNative~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml")
summary(test.mod)
loocv(test.mod) #RMSPE=0.254

plot(persNative~S1_93_11,data = obsMO, main="MO Basin Native Sp. Persistence vs. Temperature")
plot(persNative~F_MAUG_HIS,data = obsMO, main="MO Basin Native Sp. Persistence vs. Stream Size")
plot(persNative~DSbarrier,data = obsMO, main="MO Basin Native Sp. Persistence vs. Barriers")
plot(persNative~Length_km,data = obsMO, main="MO Basin Native Sp. Persistence vs. Fragment Length")
plot(persNative~nnPreds,data = obsMO, main="MO Basin Native Sp. Persistence vs. Piscivores")



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

#Temp
test.mod.sig<- spglm(
  formula = persNative~ scale(S1_93_11),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#Size
test.mod.sig<- spglm(
  formula = persNative~ scale(F_MAUG_HIS),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#Barrier
test.mod.sig<- spglm(
  formula = persNative~ scale(DSbarrier),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#ReachL
test.mod.sig<- spglm(
  formula = persNative~ scale(Length_km),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#Pisciv
test.mod.sig<- spglm(
  formula = persNative~ scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) #RMSPE=0.277

#Non-Spatial Model

test.mod.nonspatial<- spglm(
  formula = persNative~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml", spcov_type = "none")
summary(test.mod.nonspatial)
loocv(test.mod.nonspatial)



#Interaction attempts---nothing better than reach length alone
test.mod.sig<- spglm(
  formula = persNative~ scale(S1_93_11)*scale(F_MAUG_HIS)*scale(Length_km),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsMO, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 


#
#Green Basin Models
#
test.mod<- spglm(
  formula = persNative~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml")
summary(test.mod)
loocv(test.mod) 
plot(persNative~S1_93_11,data = obsGR, main="GR Basin Native Sp. Persistence vs. Temperature")
plot(persNative~F_MAUG_HIS,data = obsGR, main="GR Basin Native Sp. Persistence vs. Stream Size")
plot(persNative~DSbarrier,data = obsGR, main="GR Basin Native Sp. Persistence vs. Barriers")
plot(persNative~Length_km,data = obsGR, main="GR Basin Native Sp. Persistence vs. Fragment Length")
plot(persNative~nnPreds,data = obsGR, main="GR Basin Native Sp. Persistence vs. Piscivores")
####Compare to intercept only (null)  
test.mod.null<- spglm(
  formula = persNative~ 1,
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml")
summary(test.mod.null)
loocv(test.mod.null) 

#Temp
test.mod.sig<- spglm(
  formula = persNative~ scale(S1_93_11),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#Size
test.mod.sig<- spglm(
  formula = persNative~ scale(F_MAUG_HIS),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#Barrier
test.mod.sig<- spglm(
  formula = persNative~ scale(DSbarrier),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#ReachL
test.mod.sig<- spglm(
  formula = persNative~ scale(Length_km),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obsGR, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#Pisciv
test.mod.sig<- spglm(
  formula = persNative~ scale(nnPreds),
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

#glances(test.mod,test.mod.null,test.mod.sig,test.mod.nonspatial)






# ---- Fit model for Tributary Glacial-Relict Persistence----------------------------
obs2=obs
relict=obs2%>%pivot_longer(cols=32:124,names_to = "Species", values_to = "change")
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
relict=relict[,-c(39:53)]
#relict=relict%>%pivot_wider(names_from = "Species",values_from = "change")
mean(relict$change)

      #Create Small Glaical relict persistence column
flowavgs=read.csv("FlowAverages.csv")
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
traits=left_join(traits,flowavgs,by="Species")
smallrelicts=subset(traits, traits$Glacial=="E" & traits$FlowMedian<50)
smallrelicts=smallrelicts$Species
relict$X=NULL
obsrelict=relict
obsrelict=subset(obsrelict,obsrelict$Species%in%smallrelicts)
obsrelict=subset(obsrelict,!is.na(obsrelict$change))
obsrelcol=obsrelict%>%group_by(RepeatID)%>%summarise(persSmGlacial=mean(change))
obsrelcol$geometry=NULL
obs2=left_join(obs2,obsrelcol,by="RepeatID")
  #Convert 0 and 1 for beta regression
obs2$persSmGlacial[which(obs2$persSmGlacial==0)]=0.0000001
obs2$persSmGlacial[which(obs2$persSmGlacial==1)]=0.9999999


obs2%>%filter(!is.na(persSmGlacial))%>%summarise(mean(persSmGlacial))

#Full Model
test.mod<- spglm(
  formula = persSmGlacial~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod)
loocv(test.mod) 
plot(persSmGlacial~S1_93_11,data = obs2, main="Glacial Relict Persistence vs. Temperature")
plot(persSmGlacial~F_MAUG_HIS,data = obs2, main="Glacial Relict Persistence vs. Stream Size")
plot(persSmGlacial~DSbarrier,data = obs2, main="Glacial Relict Persistence vs. Barriers")
plot(persSmGlacial~Length_km,data = obs2, main="Glacial Relict Persistence vs. Fragment Length")
plot(persSmGlacial~nnPreds,data = obs2, main="Glacial Relict vs. Piscivores")

####Compare to intercept only (null)  
test.mod.null<- spglm(
  formula = persSmGlacial~ 1,
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.null)
loocv(test.mod.null) 

#Temp
test.mod.sig<- spglm(
  formula = persSmGlacial~ scale(S1_93_11),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#Size
test.mod.sig<- spglm(
  formula = persSmGlacial~ scale(F_MAUG_HIS),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) #0.330

#Barrier
test.mod.sig<- spglm(
  formula = persSmGlacial~ scale(DSbarrier),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 
obs_smallfrag=subset(obs2,obs2$Length_km<200)

#ReachL
test.mod.sig<- spglm(
  formula = persSmGlacial~ scale(Length_km),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#Pisciv
test.mod.sig<- spglm(
  formula = persSmGlacial~ scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#Non-Spatial Model
test.mod.nonspatial<- spglm(
  formula = persSmGlacial~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml", spcov_type = "none")
summary(test.mod.nonspatial)
loocv(test.mod.nonspatial)#0.328

#Full Model with interactions
test.mod<- spglm(
  formula = persSmGlacial~ scale(S1_93_11)*scale(F_MAUG_HIS)*scale(DSbarrier)*scale(Length_km)*scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod)
loocv(test.mod) #0.322 best model

#Refit with REML
test.mod<- spglm(
  formula = persSmGlacial~ scale(S1_93_11)*scale(F_MAUG_HIS)*scale(DSbarrier)*scale(Length_km)*scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "reml")
summary(test.mod)
loocv(test.mod) #0.313 

#No interactions with REML
test.mod<- spglm(
  formula = persSmGlacial~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "reml")
summary(test.mod)
loocv(test.mod) 

#VARIABLE IMPORTANCE PLOT
coeff=as.data.frame(test.mod$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="test.mod$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(S1_93_11)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(F_MAUG_HIS)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(DSbarrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(Length_km)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(nnPreds)")]="Piscivores"
coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point( color="purple", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ggtitle(label="Glacial Relicts")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")


#Non-Spatial REML Model with interactions#VariableNon-Spatial REML Model with interactions
nsGlacial<- spglm(
  formula = persSmGlacial~ scale(S1_93_11)*scale(F_MAUG_HIS)*scale(DSbarrier)*scale(Length_km)*scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "reml", spcov_type = "none")
summary(nsGlacial)
loocv(nsGlacial)#0.322

st_write(obs2, "SmGlacial.shp",append = F)



#----- Fit model for Periodic Life History-------------------------------------------
obs2=obs
relict=obs2%>%pivot_longer(cols=32:124,names_to = "Species", values_to = "change")
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
relict=relict[,-c(39:53)]

#Create Periodic LH persistence column
flowavgs=read.csv("FlowAverages.csv")
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
traits=left_join(traits,flowavgs,by="Species")
periodics=subset(traits, traits$LifeHist=="Per" #& traits$FlowMedian<150
                 )
#periodics=subset(periodics, periodics$Species=="WSU" | periodics$Species=="RTCH"|periodics$Species=="FMSU"|periodics$Species=="LNSU"|periodics$Species=="BHSU"|periodics$Species=="GILA"|periodics$Species=="GR")
periodics=periodics$Species
relict$X=NULL
peri=relict
peri=subset(peri,peri$Species%in%periodics)
peri=subset(peri,!is.na(peri$change))
obsrelcol=peri%>%group_by(RepeatID)%>%summarise(periodics=mean(change))
obsrelcol$geometry=NULL
obs2=left_join(obs2,obsrelcol,by="RepeatID")
#Convert 0 and 1 for beta regression
obs2$periodics[which(obs2$periodics==0)]=0.0000001
obs2$periodics[which(obs2$periodics==1)]=0.9999999
#obsMO=obs2%>%filter(PU!="GR")

obs2%>%filter(!is.na(periodics))%>%summarise(mean(periodics))



library(PerformanceAnalytics)
obsper=obs2%>%filter(!is.na(periodics))
mdata=obsper[,c(14,19,20,21,22,23,30,31)]
mdata=as.data.frame(mdata)
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)



#Full Model--Cant run with nnPreds
test.mod<- spglm(
  formula = periodics~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
full.perio.mod=test.mod
summary(test.mod)
loocv(test.mod) 
plot(periodics~S1_93_11,data = obs2, main="Periodics Persistence vs. Temperature")
plot(periodics~F_MAUG_HIS,data = obs2, main="Periodics Persistence vs. Stream Size")
plot(periodics~DSbarrier,data = obs2, main="Periodics Persistence vs. Barriers")
plot(periodics~Length_km,data = obs2, main="Periodics Persistence vs. Fragment Length")
plot(periodics~nnPreds,data = obs2, main="Periodics vs. Piscivores")

####Compare to intercept only (null)  
test.mod.null<- spglm(
  formula = periodics~ 1,
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.null)
loocv(test.mod.null) 

#Temp
test.mod.sig<- spglm(
  formula = periodics~ scale(S1_93_11),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#Size
test.mod.sig<- spglm(
  formula = periodics~ scale(F_MAUG_HIS),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) #0.330

#Barrier
test.mod.sig<- spglm(
  formula = periodics~ scale(DSbarrier),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#ReachL
test.mod.sig<- spglm(
  formula = periodics~ scale(Length_km),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) #0.382--best

#Pisciv
test.mod.sig<- spglm(
  formula = periodics~ scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#Non-Spatial Model
test.mod.nonspatial<- spglm(
  formula = periodics~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml", spcov_type = "none")
summary(test.mod.nonspatial)
loocv(test.mod.nonspatial)

#Refit with REML
test.mod<- spglm(
  formula = periodics~ scale(Length_km),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "reml")
summary(test.mod)
loocv(test.mod) #0.313

st_write(obs2, "Periodics.shp",append = F)

#PLOT FOR VAR IMP
test.mod<- spglm(
  formula = periodics~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "reml")


#VARIABLE IMPORTANCE PLOT
coeff=as.data.frame(test.mod$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="test.mod$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(S1_93_11)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(F_MAUG_HIS)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(DSbarrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(Length_km)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(nnPreds)")]="Piscivores"
coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point( color="purple", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ggtitle(label="Periodic Species")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")


#----- Fit model for Pelagic Broadcasters-------------------------------------------
obs2=obs
relict=obs2%>%pivot_longer(cols=32:124,names_to = "Species", values_to = "change")
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
relict=relict[,-c(39:53)]

#Create persistence column
flowavgs=read.csv("FlowAverages.csv")
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
traits=left_join(traits,flowavgs,by="Species")
pelagics=subset(traits, traits$Species=="STCH"|traits$Species=="SFCH"|traits$Species=="FHCH"|traits$Species=="PLMN")
pelagics=pelagics$Species
relict$X=NULL
pela=relict
pela=subset(pela,pela$Species%in%pelagics)
pela=subset(pela,!is.na(pela$change))
obsrelcol=pela%>%group_by(RepeatID)%>%summarise(pelagics=mean(change))
obsrelcol$geometry=NULL
obs2=left_join(obs2,obsrelcol,by="RepeatID")
#Convert 0 and 1 for beta regression
obs2$pelagics[which(obs2$pelagics==0)]=0.0000001
obs2$pelagics[which(obs2$pelagics==1)]=0.9999999

obs2%>%filter(!is.na(pelagics))%>%summarise(mean(pelagics))
#obsMO=obs2%>%filter(PU!="GR")
st_write(obs2, "Pelagics.shp",append = F)

library(PerformanceAnalytics)
obsper=obs2%>%filter(!is.na(pelagics))
mdata=obsper[,c(14,19,20,21,22,23,30,31)]
mdata=as.data.frame(mdata)
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)#CANT USE both Size and nnPreds or Size and DSbarrier


#Full Model--Cant run with nnPreds
test.mod<- spglm(
  formula = pelagics~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod)
loocv(test.mod) 
plot(pelagics~S1_93_11,data = obs2, main="pelagics Persistence vs. Temperature")
plot(pelagics~F_MAUG_HIS,data = obs2, main="pelagics Persistence vs. Stream Size")
plot(pelagics~DSbarrier,data = obs2, main="pelagics Persistence vs. Barriers")
plot(pelagics~Length_km,data = obs2, main="pelagics Persistence vs. Fragment Length")
plot(pelagics~nnPreds,data = obs2, main="pelagics vs. Piscivores")

####Compare to intercept only (null)  
test.mod.null<- spglm(
  formula = pelagics~ 1,
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.null)
loocv(test.mod.null) #0.442

#Temp
test.mod.sig<- spglm(
  formula = pelagics~ scale(S1_93_11),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#Size
test.mod.size<- spglm(
  formula = pelagics~ scale(F_MAUG_HIS),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.size)
loocv(test.mod.size) 

#Barrier -->no convergence
test.mod.sig<- spglm(
  formula = pelagics~ scale(DSbarrier),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
#obs_smallfrag=subset(obs2,obs2$Length_km<200)
#ReachL
test.mod.sig<- spglm(
  formula = pelagics~ scale(Length_km),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) ##best model so far 0.439

#Pisciv
test.mod.sig<- spglm(
  formula = pelagics~ scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml")
summary(test.mod.sig)
loocv(test.mod.sig) 

#Non-Spatial Model
test.mod.nonspatial<- spglm(
  formula = pelagics~ scale(Length_km),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "ml", spcov_type = "none")
summary(test.mod.nonspatial)
loocv(test.mod.nonspatial)#0.439

#Refit with REML
test.mod<- spglm(
  formula = pelagics~ scale(Length_km),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "reml")
summary(test.mod)
loocv(test.mod) #0.439 

#MODEL FOR VAR IMP
test.mod<- spglm(
  formula = pelagics~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "reml")

coeff=as.data.frame(test.mod$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="test.mod$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(S1_93_11)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(F_MAUG_HIS)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(DSbarrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(Length_km)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(nnPreds)")]="Piscivores"
coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point( color="purple", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ggtitle(label="Pelagic Broadcasters")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")

#----- Fit model for Small-Montane Sp-------------------------------------------
obs2=obs
relict=obs2%>%pivot_longer(cols=32:124,names_to = "Species", values_to = "change")
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
relict=relict[,-c(39:53)]
relict$Species[which(relict$Species=="RMCOT" | relict$Species=="COLCOT")]="MOTCOT"

#Identify small, montane species
slopes=read.csv("SlopeAverages.csv")
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
traits=left_join(traits,slopes,by="Species")
traits$X=NULL
lengths=read.csv("maxlengths.csv")
traits=left_join(traits,lengths,by="Species")
traits=subset(traits, traits$Species!="RMCOT" & traits$Species!="COLCOT")
traits%>%
  ggplot(aes(x=MaxL_cm, y=SlopeAvg, label=Species))+
  geom_text()+
  geom_point()+
  geom_hline(yintercept = 0.01)+
  geom_vline(xintercept = 26)

montane=subset(traits, traits$Species=="LSCH"|traits$Species=="MOTCOT"|traits$Species=="FSDC"|traits$Species=="HHCH")
montane=montane$Species
relict$X=NULL
mont=relict
mont=subset(mont,mont$Species%in%montane)
mont=subset(mont,!is.na(mont$change))
obsrelcol=mont%>%group_by(RepeatID)%>%summarise(montane=mean(change))
obsrelcol$geometry=NULL
obs2=left_join(obs2,obsrelcol,by="RepeatID")
#Convert 0 and 1 for beta regression
obs2$montane[which(obs2$montane==0)]=0.0000001
obs2$montane[which(obs2$montane==1)]=0.9999999
#obsMO=obs2%>%filter(PU!="GR")
st_write(obs2, "montane.shp",append = F)



library(PerformanceAnalytics)
obsper=obs2%>%filter(!is.na(montane))
mdata=obsper[,c(14,19,20,21,22,23,30,31)]
mdata=as.data.frame(mdata)
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)#CANT USE length and size

mean(obsglm$montane)
mean(obs2$persComm)

#SPATIAL MODELS ARE OFTEN NOT CONVERGING...USING NONSPATIAL MODEL
##Full
full.mont=glm(montane~scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(nnPreds), data = obs2,
    family = "binomial")
summary(full.mont)
###CROSS VALIDATION (LOOCV) 
obsglm=obs2%>%filter(!is.na(montane))
#specify the cross-validation method
ctrl <- trainControl(method = "LOOCV")
#fit a regression model and use LOOCV to evaluate performance
model <- train(montane~scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(nnPreds), data = obsglm, method = "glm", trControl = ctrl)
#view summary of LOOCV               
print(model)#0.452

##NULL
null.mont=glm(montane~1, data = obsglm, family = "binomial")
summary(null.mont)
    ###CROSS VALIDATION (LOOCV) 
obsglm$int=NA
obsglm$int=1
model <- train(montane~int, data = obsglm, method = "glm", trControl = ctrl)
print(model)#0.467

##SIGNIFICANT COVARIATES ONLY NO INTERACTIONS
sig.mont=glm(montane~scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier), data = obsglm, family = "binomial")
summary(sig.mont)
model <- train(montane~scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier), data = obsglm, method = "glm", trControl = ctrl)
print(model)#0.444--best model

##SIGNIFICANT COVARIATES ALONE--TEMP
temp.mont=glm(montane~scale(S1_93_11), data = obsglm, family = "binomial")
summary(temp.mont)
model <- train(montane~scale(S1_93_11), data = obsglm, method = "glm", trControl = ctrl)
print(model)#0.463

##SIGNIFICANT COVARIATES ALONE--SIZE
size.mont=glm(montane~scale(F_MAUG_HIS), data = obsglm, family = "binomial")
summary(size.mont)
model <- train(montane~scale(F_MAUG_HIS), data = obsglm, method = "glm", trControl = ctrl)
print(model)#0.463

##SIGNIFICANT COVARIATES ALONE--barrier
barrier.mont=glm(montane~scale(DSbarrier), data = obsglm, family = "binomial")
summary(barrier.mont)
model <- train(montane~scale(DSbarrier), data = obsglm, method = "glm", trControl = ctrl)
print(model)#0.465

#BEST MODEL WITH SPATIAL
sigSPAT=spglm(montane~scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier), 
              family = "beta",
              random = ~as.factor(Yrange),
              data = obs2, estmethod = "ml")
summary(sigSPAT)
loocv(sigSPAT)#0.439
sigSPAT.best=sigSPAT

#BEST MODEL WITH SPATIAL, refit with REML
sigSPAT=spglm(montane~scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier), 
              family = "beta",
              random = ~as.factor(Yrange),
              data = obs2, estmethod = "reml")
summary(sigSPAT)
loocv(sigSPAT)#0.440--not better

#MODEL FOR VAR IMP
test.mod<- spglm(
  formula = montane~scale(S1_93_11)+scale(F_MAUG_HIS)+scale(DSbarrier)+scale(nnPreds),
  family = "beta",
  random = ~as.factor(Yrange),
  data = obs2, estmethod = "reml")

coeff=as.data.frame(test.mod$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="test.mod$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(S1_93_11)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(F_MAUG_HIS)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(DSbarrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(Length_km)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(nnPreds)")]="Piscivores"
coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point( color="purple", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ggtitle(label="Small, Montane Species")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")

#---------------------Graph Predictions for Different Glacial-Relict Scenarios-----------
#TEMPS
obsrelict2=subset(obsrelict,obsrelict$F_MAUG_HIS<50)
median(obsrelict2$F_MAUG_HIS)#SMAll stream median = 7.238985
obsrelict2=subset(obsrelict,obsrelict$F_MAUG_HIS>50 & obsrelict$F_MAUG_HIS<200)
median(obsrelict2$F_MAUG_HIS)#Medium stream median = 97.94914
obsrelict2=subset(obsrelict,obsrelict$F_MAUG_HIS>=200)
median(obsrelict2$F_MAUG_HIS)#Large stream median = 858.207
#Split into small medium and large size streams

###################Barrier
######SMALL STREAMS
S1_93_11=seq(from=10, to=40, by=0.1)
temps=as.data.frame(S1_93_11)
temps$F_MAUG_HIS=NA
temps$F_MAUG_HIS=7.238985
temps$DSbarrier=NA
temps$DSbarrier=0
temps$Length_km=NA
temps$Length_km=median(obsrelict$Length_km)
temps$nnPreds=NA
temps$nnPreds=0
temps$Yrange=NA
temps$Yrange=max(obsrelict$Yrange)
temps$persSmGlacial=NA
temps$persSmGlacial=predict(nsGlacial, temps)
temps$Size=NA
temps$Size="Small"
temps1=temps
######MEDIUM STREAMS
S1_93_11=seq(from=10, to=40, by=0.1)
temps2=as.data.frame(S1_93_11)
temps2$F_MAUG_HIS=NA
temps2$F_MAUG_HIS=97.94914
temps2$DSbarrier=NA
temps2$DSbarrier=0
temps2$Length_km=NA
temps2$Length_km=median(obsrelict$Length_km)
temps2$nnPreds=NA
temps2$nnPreds=0
temps2$Yrange=NA
temps2$Yrange=max(obsrelict$Yrange)
temps2$persSmGlacial=NA
temps2$persSmGlacial=predict(nsGlacial, temps2)
temps2$Size=NA
temps2$Size="Medium"
######LARGE STREAMS
S1_93_11=seq(from=10, to=40, by=0.1)
temps3=as.data.frame(S1_93_11)
temps3$F_MAUG_HIS=NA
temps3$F_MAUG_HIS=858.207
temps3$DSbarrier=NA
temps3$DSbarrier=0
temps3$Length_km=NA
temps3$Length_km=median(obsrelict$Length_km)
temps3$nnPreds=NA
temps3$nnPreds=0
temps3$Yrange=NA
temps3$Yrange=max(obsrelict$Yrange)
temps3$persSmGlacial=NA
temps3$persSmGlacial=predict(nsGlacial, temps3)
temps3$Size=NA
temps3$Size="Large"
######SMALL STREAMS
S1_93_11=seq(from=10, to=40, by=0.1)
temps1.barriers=as.data.frame(S1_93_11)
temps1.barriers$F_MAUG_HIS=NA
temps1.barriers$F_MAUG_HIS=7.238985
temps1.barriers$DSbarrier=NA
temps1.barriers$DSbarrier=1
temps1.barriers$Length_km=NA
temps1.barriers$Length_km=median(obsrelict$Length_km)
temps1.barriers$nnPreds=NA
temps1.barriers$nnPreds=0
temps1.barriers$Yrange=NA
temps1.barriers$Yrange=max(obsrelict$Yrange)
temps1.barriers$persSmGlacial=NA
temps1.barriers$persSmGlacial=predict(nsGlacial, temps1.barriers)
temps1.barriers$Size=NA
temps1.barriers$Size="Small"
temps1.barriers=temps1.barriers
######MEDIUM STREAMS
S1_93_11=seq(from=10, to=40, by=0.1)
temps2.barriers=as.data.frame(S1_93_11)
temps2.barriers$F_MAUG_HIS=NA
temps2.barriers$F_MAUG_HIS=97.94914
temps2.barriers$DSbarrier=NA
temps2.barriers$DSbarrier=1
temps2.barriers$Length_km=NA
temps2.barriers$Length_km=median(obsrelict$Length_km)
temps2.barriers$nnPreds=NA
temps2.barriers$nnPreds=0
temps2.barriers$Yrange=NA
temps2.barriers$Yrange=max(obsrelict$Yrange)
temps2.barriers$persSmGlacial=NA
temps2.barriers$persSmGlacial=predict(nsGlacial, temps2.barriers)
temps2.barriers$Size=NA
temps2.barriers$Size="Medium"
######LARGE STREAMS
S1_93_11=seq(from=10, to=40, by=0.1)
temps3.barriers=as.data.frame(S1_93_11)
temps3.barriers$F_MAUG_HIS=NA
temps3.barriers$F_MAUG_HIS=858.207
temps3.barriers$DSbarrier=NA
temps3.barriers$DSbarrier=1
temps3.barriers$Length_km=NA
temps3.barriers$Length_km=median(obsrelict$Length_km)
temps3.barriers$nnPreds=NA
temps3.barriers$nnPreds=0
temps3.barriers$Yrange=NA
temps3.barriers$Yrange=max(obsrelict$Yrange)
temps3.barriers$persSmGlacial=NA
temps3.barriers$persSmGlacial=predict(nsGlacial, temps3.barriers)
temps3.barriers$Size=NA
temps3.barriers$Size="Large"

temps=rbind(temps1,temps2)
temps=rbind(temps, temps3)
temps=rbind(temps,temps1.barriers)
temps=rbind(temps,temps2.barriers)
temps=rbind(temps, temps3.barriers)
temps$size_f = factor(temps$Size, levels=c("Small","Medium","Large"))

temps$DSbarrier[which(temps$DSbarrier==1)]="Yes"
temps$DSbarrier[which(temps$DSbarrier==0)]="No"
temps.barriers=temps%>%
  ggplot(aes(x=S1_93_11,y=persSmGlacial, linetype = as.factor(DSbarrier)))+
  geom_smooth(method = "lm", se=F, color="black")+
  scale_y_continuous(limits = c(0,1))+
  theme_classic()+
  facet_wrap(~size_f)+
  ylab(label="Proportion of Species Persisting")+
  xlab(label = "Mean August Stream Temperature")+
  labs(linetype = "Downstream Barrier")+
  theme(legend.title = element_text(size = 8))
temps.barriers

###################FragmentsOverTemp
plot(density(log(obsrelict$Length_km)))
median(log(obsrelict$Length_km))
obsrelict2=subset(obsrelict,obsrelict$Length_km<68.21101)
median(obsrelict2$Length_km)#30.423
obsrelict2=subset(obsrelict,obsrelict$Length_km>=68.21101)
median(obsrelict2$Length_km)#128.81

######SMALL STREAMS
S1_93_11=seq(from=10, to=40, by=0.1)
temps=as.data.frame(S1_93_11)
temps$F_MAUG_HIS=NA
temps$F_MAUG_HIS=7.238985
temps$DSbarrier=NA
temps$DSbarrier=0
temps$Length_km=NA
temps$Length_km=30.423
temps$nnPreds=NA
temps$nnPreds=0
temps$Yrange=NA
temps$Yrange=max(obsrelict$Yrange)
temps$persSmGlacial=NA
temps$persSmGlacial=predict(nsGlacial, temps)
temps$Size=NA
temps$Size="Small"
temps$Fragment=NA
temps$Fragment="Short"
temps1=temps
######MEDIUM STREAMS
S1_93_11=seq(from=10, to=40, by=0.1)
temps2=as.data.frame(S1_93_11)
temps2$F_MAUG_HIS=NA
temps2$F_MAUG_HIS=97.94914
temps2$DSbarrier=NA
temps2$DSbarrier=0
temps2$Length_km=NA
temps2$Length_km=30.423
temps2$nnPreds=NA
temps2$nnPreds=0
temps2$Yrange=NA
temps2$Yrange=max(obsrelict$Yrange)
temps2$persSmGlacial=NA
temps2$persSmGlacial=predict(nsGlacial, temps2)
temps2$Size=NA
temps2$Size="Medium"
temps2$Fragment=NA
temps2$Fragment="Short"
######LARGE STREAMS
S1_93_11=seq(from=10, to=40, by=0.1)
temps3=as.data.frame(S1_93_11)
temps3$F_MAUG_HIS=NA
temps3$F_MAUG_HIS=858.207
temps3$DSbarrier=NA
temps3$DSbarrier=0
temps3$Length_km=NA
temps3$Length_km=30.423
temps3$nnPreds=NA
temps3$nnPreds=0
temps3$Yrange=NA
temps3$Yrange=max(obsrelict$Yrange)
temps3$persSmGlacial=NA
temps3$persSmGlacial=predict(nsGlacial, temps3)
temps3$Size=NA
temps3$Size="Large"
temps3$Fragment=NA
temps3$Fragment="Short"
######SMALL STREAMS
S1_93_11=seq(from=10, to=40, by=0.1)
temps1.barriers=as.data.frame(S1_93_11)
temps1.barriers$F_MAUG_HIS=NA
temps1.barriers$F_MAUG_HIS=7.238985
temps1.barriers$DSbarrier=NA
temps1.barriers$DSbarrier=0
temps1.barriers$Length_km=NA
temps1.barriers$Length_km=128.81
temps1.barriers$nnPreds=NA
temps1.barriers$nnPreds=0
temps1.barriers$Yrange=NA
temps1.barriers$Yrange=max(obsrelict$Yrange)
temps1.barriers$persSmGlacial=NA
temps1.barriers$persSmGlacial=predict(nsGlacial, temps1.barriers)
temps1.barriers$Size=NA
temps1.barriers$Size="Small"
temps1.barriers=temps1.barriers
temps1.barriers$Fragment=NA
temps1.barriers$Fragment="Long"
######MEDIUM STREAMS
S1_93_11=seq(from=10, to=40, by=0.1)
temps2.barriers=as.data.frame(S1_93_11)
temps2.barriers$F_MAUG_HIS=NA
temps2.barriers$F_MAUG_HIS=97.94914
temps2.barriers$DSbarrier=NA
temps2.barriers$DSbarrier=0
temps2.barriers$Length_km=NA
temps2.barriers$Length_km=128.81
temps2.barriers$nnPreds=NA
temps2.barriers$nnPreds=0
temps2.barriers$Yrange=NA
temps2.barriers$Yrange=max(obsrelict$Yrange)
temps2.barriers$persSmGlacial=NA
temps2.barriers$persSmGlacial=predict(nsGlacial, temps2.barriers)
temps2.barriers$Size=NA
temps2.barriers$Size="Medium"
temps2.barriers$Fragment=NA
temps2.barriers$Fragment="Long"
######LARGE STREAMS
S1_93_11=seq(from=10, to=40, by=0.1)
temps3.barriers=as.data.frame(S1_93_11)
temps3.barriers$F_MAUG_HIS=NA
temps3.barriers$F_MAUG_HIS=858.207
temps3.barriers$DSbarrier=NA
temps3.barriers$DSbarrier=0
temps3.barriers$Length_km=NA
temps3.barriers$Length_km=128.81
temps3.barriers$nnPreds=NA
temps3.barriers$nnPreds=0
temps3.barriers$Yrange=NA
temps3.barriers$Yrange=max(obsrelict$Yrange)
temps3.barriers$persSmGlacial=NA
temps3.barriers$persSmGlacial=predict(nsGlacial, temps3.barriers)
temps3.barriers$Size=NA
temps3.barriers$Size="Large"
temps3.barriers$Fragment=NA
temps3.barriers$Fragment="Long"
#Merge
temps=rbind(temps1,temps2)
temps=rbind(temps, temps3)
temps=rbind(temps,temps1.barriers)
temps=rbind(temps,temps2.barriers)
temps=rbind(temps, temps3.barriers)
temps$size_f = factor(temps$Size, levels=c("Small","Medium","Large"))

temps.fragments=temps%>%
  ggplot(aes(x=S1_93_11,y=persSmGlacial, linetype = as.factor(Fragment)))+
  geom_smooth(method = "lm", se=F, color="black")+
  scale_y_continuous(limits = c(0,1))+
  theme_classic()+
  facet_wrap(~size_f)+
  ylab(label="Proportion of Species Persisting")+
  xlab(label = "Mean August Stream Temperature")+
  labs(linetype = "Fragment  Length  ")+
  theme(legend.title = element_text(size = 8))

temps.fragments

library(ggpubr)
relicts.predicted=ggarrange(temps.barriers,temps.fragments,ncol=1)
annotate_figure(relicts.predicted,top = text_grob("Glacial-Relict Persistence", face = "bold", size = 14))
ggsave(filename = "GlacialRelictsPredicted.tiff",dpi=400, width = 10, height = 6, units = "in")


#=====================================================================
#============================SECTION 3: Examine Net Change for Species
#=====================================================================
#---------------------Net Species Change------------------------------------------
obs3=obs2
obs3=as.data.frame(obs3)
obs3%>%filter(!is.na(persGlacE))%>%summarise(avg=mean(persGlacE), sd(persGlacE)) #58% (sd=34%)
obs3%>%filter(!is.na(persGlacL))%>%summarise(avg=mean(persGlacL), sd(persGlacL)) #53% (sd=43%)
obs3%>%filter(!is.na(persComm))%>%summarise(avg=mean(persComm), sd(persComm)) #56% (sd=30%)
obs3%>%filter(!is.na(persNative))%>%summarise(avg=mean(persNative),sd(persNative)) #56% (sd=32%)

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
NETsub=subset(NET,NET$EarlyNum>=15|NET$LateNum>=15)
NETsub=subset(NETsub,NETsub$Species!="ONC" & NETsub$Species!="FMxWSU")
mean(NETsub$propChange)


traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
traits$Species[which(traits$Species=="COLCOT")]="MOTCOT"
NETsub=left_join(NETsub,traits,by="Species")
NETsub$CommonName[which(NETsub$CommonName=="Northern Redbelly Dace")]="Nor. Redbelly Dace"
NETsub$CommonName[which(NETsub$CommonName=="Rocky Mountain Cutthroat Trout")]="R.M. Cutthroat Trout"

declining=NETsub%>%filter(propChange<1)%>%
  ggplot(aes(x=reorder(CommonName,propChange),y=propChange))+
  geom_segment( aes(x=CommonName, xend=CommonName, y=1, yend=propChange), color="grey") +
  geom_point(aes(colour = as.factor(Status)), size=3)+
  scale_color_manual(values = c("#ff9999","#005555"))+
  scale_x_discrete(position = "top")+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0.4,1))+
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
  scale_y_continuous(breaks=seq(1,4,by=1), limits = c(1,4), labels = c("1.0","2.0","3.0","4.0"))+
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





#---------------------Net Species Changes to Only Native Populations------------------------------------------
NET=read.csv("PersistenceMetrics.csv")
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
NET=subset(NET,NET$Status=="Native")
NET2=NET%>%group_by(Species)%>%summarise(net=sum(change))
NET3=NET%>%filter(change!=-1)%>%group_by(Species)%>%summarise(LateNum=length(change))
NET=NET%>%filter(change!=1)%>%group_by(Species)%>%summarise(EarlyNum=length(change))
NET=left_join(NET2,NET,by="Species")
NET=left_join(NET,NET3,by="Species")
NET$EarlyNum[which(is.na(NET$EarlyNum))]=0
NET$LateNum[which(is.na(NET$LateNum))]=0
NET$propChange=NA
NET$propChange=NET$LateNum/NET$EarlyNum
NETsub=subset(NET,NET$EarlyNum>=15|NET$LateNum>=15)
NETsub=subset(NETsub,NETsub$Species!="ONC" & NETsub$Species!="FMxWSU")
mean(NETsub$propChange)

traits2=traits[,c(1,3)]
NETsub=left_join(NETsub,traits2,by="Species")
NETsub$CommonName[which(NETsub$CommonName=="Northern Redbelly Dace")]="Nor. Redbelly Dace"
NETsub$CommonName[which(NETsub$CommonName=="Rocky Mountain Cutthroat Trout")]="R.M. Cutthroat Trout"


declining=NETsub%>%filter(propChange<1)%>%
  ggplot(aes(x=reorder(CommonName,propChange),y=propChange))+
  geom_segment( aes(x=CommonName, xend=CommonName, y=1, yend=propChange), color="grey") +
  geom_point(size=3, color="#005555")+
  scale_x_discrete(position = "top")+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0.4,1))+
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
ggsave(filename="DecliningSpNATIVE.tiff",dpi = 400, width = 8, height = 6, units = "in")
increasing=NETsub%>%filter(propChange>=1)%>%
  ggplot(aes(x=reorder(CommonName,-propChange),y=propChange))+
  geom_segment( aes(x=CommonName, xend=CommonName, y=1, yend=propChange), color="grey") +
  geom_point(size=3, color="#005555")+
  theme_light()+
  ylab(label = "Net Change in Proportional Occurrence")+
  scale_y_continuous(breaks=seq(1,2,by=0.5), limits = c(1,2))+
  xlab(label="")+
  ggtitle(label = "Increasing or Stable Species")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face = 2),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size=12),
        axis.ticks.x = element_blank())
increasing
#ggsave(filename="IncreasingSpNative.tiff",dpi = 400, width = 8, height = 6, units = "in")



#===========================================================================
#============================SECTION 4: Functional Groups
#===========================================================================
###Uncorrected Prop Change by Functional Groups
traits$Glacial[which(traits$Species=="BLBH")]="L"
NETsub2=left_join(NETsub,traits,by="Species")
thresh=read.csv("Thresholds.csv")
NETsub2=left_join(NETsub2,thresh,by = "Species")
flowavgs=read.csv("FlowAverages.csv")
NETsub2=left_join(NETsub2,flowavgs,by="Species")
#Glacial Relicts
NETsubR=NETsub2%>%filter(Glacial=="E")
NETsubNR=NETsub2%>%filter(Glacial=="L")
t.test(NETsubR$propChange,NETsubNR$propChange) # dif = 0.14, p=0.40
     ##CHI SQ of Number of Glacial Relicts declining vs non-relicts declining is p=0.138
#Small-Medium STream Glacial Relicts
NETsmallGR=NETsubR%>%filter(FlowMedian<50)
smallglacials=NETsmallGR$Species
NETother=NETsub2%>%filter(!(Species%in%smallglacials))
mean(NETsmallGR$propChange)
t.test(NETsmallGR$propChange,NETother$propChange)

#Pelagic Broadcasters #limit to MO basin
NETsubPB=NETsub2%>%filter(ReproductiveGuild=="PB")
NETsubPB=NETsubPB%>%filter(Species!="GE")
unique(NETsubPB$Species)
NETsubNPB=NETsub2%>%filter(ReproductiveGuild!="PB" & NativeBasin!="GR")
t.test(NETsubPB$propChange,NETsubNPB$propChange) #p=0.94
     ##Greater proportion of non-pelagics are declining compared to pelagics

#Thermal Guilds
NETcold=NETsub2%>%filter(SimpleThermal=="Cold")#only 1 fish...Cutthroat, don't examine
NETcool=NETsub2%>%filter(SimpleThermal=="Cool")
NETwarm=NETsub2%>%filter(SimpleThermal=="Warm")
mean(NETcool$propChange)
mean(NETwarm$propChange)
t.test(NETcool$propChange,NETwarm$propChange)#dif = 0.16, p=0.22


#Life History Strategies
NETeq=NETsub2%>%filter(LifeHist=="Equil")
NETopp=NETsub2%>%filter(LifeHist=="Opp")
NETper=NETsub2%>%filter(LifeHist=="Per")
LH.anova=aov(NETsub2$propChange~NETsub2$LifeHist)
summary(LH.anova)
TukeyHSD(LH.anova)
ggplot(NETsub2, aes(x=LifeHist,y=propChange,label=Species))+
  geom_violin(trim=F)+
  geom_text()

mean(NETeq$propChange)#1.37
mean(NETopp$propChange)#0.91
mean(NETper$propChange)#1.04

NETsub2%>%
  ggplot(aes(x=log(FlowMedian), y=propChange, color = LifeHist,label=Species))+
  #geom_point()+
  geom_vline(xintercept = 5.010635)+ #corresponds to mean august flow of 150 cfs
  geom_hline(yintercept = 1, color="grey", linetype="dashed")+
  scale_color_manual(values = c("red","blue","gold"))+
  geom_text(size=3)+
  geom_smooth(method="lm",se=F)+
  theme_classic()



######GLACIAL RELICT LOSS BY STREAM SIZE PREFERENCE
obs4=obs2%>%pivot_longer(cols = 32:124,names_to = "Species")
obs4=obs4%>%filter(!is.na(value))%>%group_by(Species)%>%summarise(meanStreamSize=mean(F_MAUG_HIS))
obs4$geometry=NULL
NETsub2=left_join(NETsub2,obs4,by="Species")

NETsub2%>%
  ggplot(aes(x=log(FlowMedian), y=propChange, color = Glacial,label=Species))+
  #geom_point()+
  geom_vline(xintercept = 3.912)+ #corresponds to mean august flow of 50 cfs
  geom_hline(yintercept = 1, color="grey", linetype="dashed")+
  scale_color_manual(values = c("#005555","#ff9999"))+
  geom_text(size=3)+
  geom_smooth(method="lm",se=F)+
  theme_classic()

NETGR=NETsub2%>%filter(Glacial=="L")%>%filter(Species!="MWF"&Species!="RMCT")
SizeRelict.lm=lm(NETGR$propChange~NETGR$meanStreamSize)
summary(SizeRelict.lm)









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
#============================SECTION 5: Individual Species Models
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
plot(BRMN~S1_93_11,data = nssn$obs, main="BRMN Refugia vs. Temperature")
plot(BRMN~F_MAUG_HIS,data = nssn$obs, main="BRMN Refugia vs. Stream Size")
plot(jitter(BRMN,.2)~DSbarrier, data = nssn$obs, main="BRMN Refugia vs. Barriers (Jittered)")
plot(BRMN~Length_km,data = nssn$obs, main="BRMN Refugia vs. Fragment Length")
plot(jitter(BRMN,.2)~nnPreds,data = nssn$obs, main="BRMN Refugia vs. Piscivores (Jittered)")

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
plot(COLCOT~S1_93_11,data = nssn$obs, main="COLCOT Refugia vs. Temperature")
plot(COLCOT~F_MAUG_HIS,data = nssn$obs, main="COLCOT Refugia vs. Stream Size")
plot(jitter(COLCOT,.2)~DSbarrier, data = nssn$obs, main="COLCOT Refugia vs. Barriers (Jittered)")
plot(COLCOT~Length_km,data = nssn$obs, main="COLCOT Refugia vs. Fragment Length")
plot(jitter(COLCOT,.2)~nnPreds,data = nssn$obs, main="COLCOT Refugia vs. Piscivores (Jittered)")


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
plot(FHCH~S1_93_11,data = nssn$obs, main="FHCH Refugia vs. Temperature")
plot(FHCH~F_MAUG_HIS,data = nssn$obs, main="FHCH Refugia vs. Stream Size")
plot(jitter(FHCH,.2)~DSbarrier, data = nssn$obs, main="FHCH Refugia vs. Barriers (Jittered)")
plot(FHCH~Length_km,data = nssn$obs, main="FHCH Refugia vs. Fragment Length")
plot(jitter(FHCH,.2)~nnPreds,data = nssn$obs, main="FHCH Refugia vs. Piscivores (Jittered)")


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
plot(FHMN~S1_93_11,data = nssn$obs, main="FHMN Refugia vs. Temperature")
plot(FHMN~F_MAUG_HIS,data = nssn$obs, main="FHMN Refugia vs. Stream Size")
plot(jitter(FHMN,.2)~DSbarrier, data = nssn$obs, main="FHMN Refugia vs. Barriers (Jittered)")
plot(FHMN~Length_km,data = nssn$obs, main="FHMN Refugia vs. Fragment Length")
plot(jitter(FHMN,.2)~nnPreds,data = nssn$obs, main="FHMN Refugia vs. Piscivores (Jittered)")


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
plot(FMSU~S1_93_11,data = nssn$obs, main="FMSU Refugia vs. Temperature")
plot(FMSU~F_MAUG_HIS,data = nssn$obs, main="FMSU Refugia vs. Stream Size")
plot(jitter(FMSU,.2)~DSbarrier, data = nssn$obs, main="FMSU Refugia vs. Barriers (Jittered)")
plot(FMSU~Length_km,data = nssn$obs, main="FMSU Refugia vs. Fragment Length")
plot(jitter(FMSU,.2)~nnPreds,data = nssn$obs, main="FMSU Refugia vs. Piscivores (Jittered)")



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
plot(LKCH~S1_93_11,data = nssn$obs, main="LKCH Refugia vs. Temperature")
plot(LKCH~F_MAUG_HIS,data = nssn$obs, main="LKCH Refugia vs. Stream Size")
plot(jitter(LKCH,.2)~DSbarrier, data = nssn$obs, main="LKCH Refugia vs. Barriers (Jittered)")
plot(LKCH~Length_km,data = nssn$obs, main="LKCH Refugia vs. Fragment Length")
plot(jitter(LKCH,.2)~nnPreds,data = nssn$obs, main="LKCH Refugia vs. Piscivores (Jittered)")



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
plot(LNDC~S1_93_11,data = nssn$obs, main="LNDC Refugia vs. Temperature")
plot(LNDC~F_MAUG_HIS,data = nssn$obs, main="LNDC Refugia vs. Stream Size")
plot(jitter(LNDC,.2)~DSbarrier, data = nssn$obs, main="LNDC Refugia vs. Barriers (Jittered)")
plot(LNDC~Length_km,data = nssn$obs, main="LNDC Refugia vs. Fragment Length")
plot(jitter(LNDC,.2)~nnPreds,data = nssn$obs, main="LNDC Refugia vs. Piscivores (Jittered)")



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
plot(LNSU~S1_93_11,data = nssn$obs, main="LNSU Refugia vs. Temperature")
plot(LNSU~F_MAUG_HIS,data = nssn$obs, main="LNSU Refugia vs. Stream Size")
plot(jitter(LNSU,.2)~DSbarrier, data = nssn$obs, main="LNSU Refugia vs. Barriers (Jittered)")
plot(LNSU~Length_km,data = nssn$obs, main="LNSU Refugia vs. Fragment Length")
plot(jitter(LNSU,.2)~nnPreds,data = nssn$obs, main="LNSU Refugia vs. Piscivores (Jittered)")



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
plot(MTSU~S1_93_11,data = nssn$obs, main="MTSU Refugia vs. Temperature")
plot(MTSU~F_MAUG_HIS,data = nssn$obs, main="MTSU Refugia vs. Stream Size")
plot(jitter(MTSU,.2)~DSbarrier, data = nssn$obs, main="MTSU Refugia vs. Barriers (Jittered)")
plot(MTSU~Length_km,data = nssn$obs, main="MTSU Refugia vs. Fragment Length")
plot(jitter(MTSU,.2)~nnPreds,data = nssn$obs, main="MTSU Refugia vs. Piscivores (Jittered)")


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
plot(PLMN~S1_93_11,data = nssn$obs, main="PLMN Refugia vs. Temperature")
plot(PLMN~F_MAUG_HIS,data = nssn$obs, main="PLMN Refugia vs. Stream Size")
plot(jitter(PLMN,.2)~DSbarrier, data = nssn$obs, main="PLMN Refugia vs. Barriers (Jittered)")
plot(PLMN~Length_km,data = nssn$obs, main="PLMN Refugia vs. Fragment Length")
plot(jitter(PLMN,.2)~nnPreds,data = nssn$obs, main="PLMN Refugia vs. Piscivores (Jittered)")


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
plot(PLSU~S1_93_11,data = nssn$obs, main="PLSU Refugia vs. Temperature")
plot(PLSU~F_MAUG_HIS,data = nssn$obs, main="PLSU Refugia vs. Stream Size")
plot(jitter(PLSU,.2)~DSbarrier, data = nssn$obs, main="PLSU Refugia vs. Barriers (Jittered)")
plot(PLSU~Length_km,data = nssn$obs, main="PLSU Refugia vs. Fragment Length")
plot(jitter(PLSU,.2)~nnPreds,data = nssn$obs, main="PLSU Refugia vs. Piscivores (Jittered)")



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
plot(SPDC~S1_93_11,data = nssn$obs, main="SPDC Refugia vs. Temperature")
plot(SPDC~F_MAUG_HIS,data = nssn$obs, main="SPDC Refugia vs. Stream Size")
plot(jitter(SPDC,.2)~DSbarrier, data = nssn$obs, main="SPDC Refugia vs. Barriers (Jittered)")
plot(SPDC~Length_km,data = nssn$obs, main="SPDC Refugia vs. Fragment Length")
plot(jitter(SPDC,.2)~nnPreds,data = nssn$obs, main="SPDC Refugia vs. Piscivores (Jittered)")



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
plot(WSU~S1_93_11,data = nssn$obs, main="WSU Refugia vs. Temperature")
plot(WSU~F_MAUG_HIS,data = nssn$obs, main="WSU Refugia vs. Stream Size")
plot(jitter(WSU,.2)~DSbarrier, data = nssn$obs, main="WSU Refugia vs. Barriers (Jittered)")
plot(WSU~Length_km,data = nssn$obs, main="WSU Refugia vs. Fragment Length")
plot(jitter(WSU,.2)~nnPreds,data = nssn$obs, main="WSU Refugia vs. Piscivores (Jittered)")





########Introduced Species
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
plot(BLBH~S1_93_11,data = nssn$obs, main="BLBH Refugia vs. Temperature")
plot(BLBH~F_MAUG_HIS,data = nssn$obs, main="BLBH Refugia vs. Stream Size")
plot(jitter(BLBH,.2)~DSbarrier, data = nssn$obs, main="BLBH Refugia vs. Barriers (Jittered)")
plot(BLBH~Length_km,data = nssn$obs, main="BLBH Refugia vs. Fragment Length")
plot(jitter(BLBH,.2)~nnPreds,data = nssn$obs, main="BLBH Refugia vs. Piscivores (Jittered)")

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
plot(EB~S1_93_11,data = nssn$obs, main="EB Refugia vs. Temperature")
plot(EB~F_MAUG_HIS,data = nssn$obs, main="EB Refugia vs. Stream Size")
plot(jitter(EB,.2)~DSbarrier, data = nssn$obs, main="EB Refugia vs. Barriers (Jittered)")
plot(EB~Length_km,data = nssn$obs, main="EB Refugia vs. Fragment Length")
plot(jitter(EB,.2)~nnPreds,data = nssn$obs, main="EB Refugia vs. Piscivores (Jittered)")

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
plot(RSSH~S1_93_11,data = nssn$obs, main="RSSH Refugia vs. Temperature")
plot(RSSH~F_MAUG_HIS,data = nssn$obs, main="RSSH Refugia vs. Stream Size")
plot(jitter(RSSH,.2)~DSbarrier, data = nssn$obs, main="RSSH Refugia vs. Barriers (Jittered)")
plot(RSSH~Length_km,data = nssn$obs, main="RSSH Refugia vs. Fragment Length")
plot(jitter(RSSH,.2)~nnPreds,data = nssn$obs, main="RSSH Refugia vs. Piscivores (Jittered)")



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


#Introduced Species
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

