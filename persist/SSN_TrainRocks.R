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
obs <- st_read("newmetrics.shp")
removes=read.csv("REMOVES.csv")
removes=removes[,c(1,3)]
removes=removes%>%rename("RepetID"="RepeatID")
obs=left_join(obs,removes,by="RepetID")
obs=subset(obs,obs$GearMiss!="Y")

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
#============================SECTION 2: SSN Models for Community & Functional Level
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




#-----------Full Community Models-------------------------------------------------------------------------

quantile(nssn$obs$Yrange)
nssn$obs$YrangeCat=NA
nssn$obs$YrangeCat[which(nssn$obs$Yrange<16)]="S"
nssn$obs$YrangeCat[which(nssn$obs$Yrange>=16 & nssn$obs$Yrange<22)]="M"
nssn$obs$YrangeCat[which(nssn$obs$Yrange>=22 & nssn$obs$Yrange<30)]="L"
nssn$obs$YrangeCat[which(nssn$obs$Yrange>=30)]="XL"


#full model -- max likelihood
ssn_mod <- ssn_glm(
  formula =   pComm~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
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
View(nssn$obs)

#full model -- max likelihood
ssn_modN <- ssn_glm(
  formula =   pNat~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
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
loocv(ssn_modN)#0.294
loocv(ssn_null)#0.299

#-------------------Native Sp without "m" gear sites
nssn2=nssn
nssn2$obs$pNat[which(nssn2$obs$GearMiss=="M")]=NA

#full nssn#full model -- max likelihood
ssn_modN <- ssn_glm(
  formula =   pNat~ scale(temp)+scale(size)+scale(barrier)+scale(length)+scale(pisc),
  family = "beta",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)
summary(ssn_modN)

#null model -- max likelihood
ssn_null <- ssn_glm(
  formula =   pNat~ 1,
  family = "beta",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)

glances(ssn_modN,ssn_null)
loocv(ssn_modN)#0.289
loocv(ssn_null)#0.293


#MODEL FOR VAR IMP
coeff=as.data.frame(ssn_modN$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ssn_modN$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size*"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
    legend.position = "none",
    plot.title = element_text(hjust = -0.2)
  ) +
  ggtitle(label="A. Native Species Persistence")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")
ggsave(filename = "VarImp_pNative.tiff",height = 4, width = 4, units = "in", dpi=400)


#-----------Glacial Relict Persistence Models
#-----------Glacial Relict Refugia Models-------------------------------------------------------------------------
mdata=subset(nssn$obs,!is.na(nssn$obs$rGlac))
mdata=mdata[,c(16,21,23,24,32)]
mdata=as.data.frame(mdata)
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)

#full model -- max likelihood
disp_init <- dispersion_initial(family = "beta", dispersion = 1, known = "dispersion")


ssn_modG <- ssn_glm(
  formula =   pNat~ scale(temp)+scale(log(size))+scale(barrier)+scale(length)+scale(pisc),
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

glances(ssn_modG,ssn_null)
loocv(ssn_modG)#0.325
loocv(ssn_null)#0.335



#VAR IMP
coeff=as.data.frame(ssn_modG$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ssn_modG$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size*"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
    legend.position = "none",
    plot.title = element_text(hjust = -0.2)
  ) +
  ggtitle(label="C. Postglacial Pioneer Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")
ggsave(filename = "VarImp_rGlacial.tiff",height = 4, width = 4, units = "in", dpi=400)

#-----------Glacial Relict Persistence-------------------------
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
mean(relict$change)

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
obsrelcol=obsrelcol%>%rename("pGlac"="persSmGlacial")
#Add to nssn
nssn$obs=left_join(nssn$obs,obsrelcol,by="RepeatID")


#------------------------------build models------------------------------------------------------
nssn$obs$pGlac[which(nssn$obs$pGlac==1)]=0.999
nssn$obs$pGlac[which(nssn$obs$pGlac==0)]=0.001
#st_write(nssn$obs,"testpGlac.shp")

#full model not working...try for just GREEN basin
nssn2=nssn
nssn2$obs$pGlac[which(nssn2$obs$PU!="GREEN")]=NA
#full nssn#full model -- max likelihood
disp_init <- dispersion_initial(family = "beta", dispersion = 1, known = "dispersion")
hist(nssn2$obs$pGlac)
nssn2$obs%>%filter(PU=="GREEN")%>%summarise(length(RepeatID))

ssn_modpG <- ssn_glm(
  formula =   pGlac~ scale(temp)+scale(log(size))+scale(length)+scale(pisc),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  #euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  nugget_type = "none"
)
summary(ssn_modpG)

#null model -- will only converge without euclidean distance
ssn_null <- ssn_glm(
  formula =   pGlac~ 1,
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
 # euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  nugget_type = "none"
)

glances(ssn_modG,ssn_null)
loocv(ssn_modpG)#0.356
loocv(ssn_null)#0.368



#VAR IMP
coeff=as.data.frame(ssn_modG$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ssn_modG$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size*"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
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
    legend.position = "none",
    plot.title = element_text(hjust = -0.2)
  ) +
  ggtitle(label="C. Postglacial Pioneer Refugia")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")
ggsave(filename = "VarImp_rGlacial.tiff",height = 4, width = 4, units = "in", dpi=400)



#-----------Non-Glacial Relict Persistence Models-------------
#-----------Species Turnover--------------------------
#cant use piscivore metrics since it is included in turnover response
turn=read.csv("TurnoverMetric.csv")
#nssn$obs=left_join(nssn$obs, turn, by="RepeatID")
nssn$obs$X=NULL

nssn$obs$Turnover[which(nssn$obs$Turnover==0)]=0.001
nssn$obs$Turnover[which(nssn$obs$Turnover==1)]=0.999

ssn_modT <- ssn_glm(
  formula =   Turnover~ scale(temp)+scale(size)+scale(barrier)+scale(length),
  family = "beta",
  ssn.object = nssn,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)
summary(ssn_modT)

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

glances(ssn_modT,ssn_null)
loocv(ssn_modT)#0.203
loocv(ssn_null)#0.204



#MODEL FOR VAR IMP
coeff=as.data.frame(ssn_modT$coefficients$fixed)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ssn_modT$coefficients$fixed")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
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
  ggtitle(label="B. Community Turnover")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")
ggsave(filename = "VarImp_Turn.tiff",height = 4, width = 4, units = "in", dpi=400)


#st_write(nssn$obs, dsn = "metrics.shp")


#---------------------Graph Predictions for Different Scenarios-----------

#TEMPS
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
  ylab(label="rValue")+
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
  ylab(label="rValue")+
  xlab(label = "Mean August Stream Temperature")+
  labs(linetype = "Fragment  Length  ")+
  theme(legend.title = element_text(size = 8))

temps.fragments

library(ggpubr)
relicts.predicted=ggarrange(temps.barriers,temps.fragments,ncol=1)
annotate_figure(relicts.predicted,top = text_grob("Glacial-Relict Refugia", face = "bold", size = 14))
ggsave(filename = "GlacialRelictsPredicted.tiff",dpi=400, width = 10, height = 6, units = "in")



ggplot(nssn$obs, aes(x=temp, y=pNat, size=log(size))) +
  geom_point(alpha=0.5)+
  scale_size(range = c(1, 7), name="Native Sp Persistence")+
  scale_color_viridis_b()








#==================================================================================
#============================SECTION 3: SSN Models for Individual Species
#==================================================================================

obs3=obs
obs3$geometry=NULL
obs3=obs3%>%pivot_longer(cols=32:124, names_to = "Species")%>%
  filter(!is.na(value))%>%
  group_by(Species)%>%
  summarise(sites=length(RepetID))

traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
nssn2=nssn
nssn2$obs=nssn2$obs%>%pivot_longer(cols = c(34:126,140), names_to = "Species")
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
  random = ~as.factor(YrangeCat),
  additive = "afvArea", estmethod = "ml")
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


#==================================================================================
#============================SECTION 4: Old models for other (unused) functional groups
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

#==================================================================================
#============================SECTION 4: SSN ANOVAs for basin
#==================================================================================
Bstats=read.csv("basinstats.csv")
Bstats$huc6[which(Bstats$huc6==100200|Bstats$huc6==100301|Bstats$huc6==100302)]="100200"
Bmeans=Bstats%>%
  group_by(huc6)%>%
  summarise(mpNat=mean(pNat),mTurn=mean(Turnovr), mColo=mean(nClnzSp),mrGlac=mean(rGlac))
write.csv(Bmeans, file = "BasinMeans.csv")

basinSSN=nssn
Bstats=Bstats[,c(2,8)]
basinSSN$obs=left_join(basinSSN$obs, Bstats,by="RepeatID")
basinSSN$obs$huc6=as.factor(basinSSN$obs$huc6)
class(basinSSN$obs$huc6)
quantile(basinSSN$obs$size)
basinSSN$obs$sizeCat=NA
basinSSN$obs$sizeCat[which(nssn$obs$size<5.32)]="s"
basinSSN$obs$sizeCat[which(nssn$obs$size>=5.32 & nssn$obs$size<22.3)]="m"
basinSSN$obs$sizeCat[which(nssn$obs$size>=22.3 & nssn$obs$size<182)]="l"
basinSSN$obs$sizeCat[which(nssn$obs$size>=182)]="xl"


#pNat
ssn_mod <- ssn_glm(
  formula =   pNat~huc6,
  family = "beta",
  ssn.object = basinSSN,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat)+as.factor(sizeCat),
  spcov_type = "exponential",
  additive = "afvArea",
  nugget_type = "none"
)
summary(ssn_mod) #no basins sig dif
anova(ssn_mod)

#Turnover
ssn_modT <- ssn_glm(
  formula =   Turnover~huc6,
  family = "beta",
  ssn.object = basinSSN,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat)+as.factor(sizeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)
summary(ssn_modT) #no basins sig dif
anova(ssn_modT)

#Glacial
ssn_modG <- ssn_glm(
  formula =   rGlac~huc6,
  family = "beta",
  ssn.object = basinSSN,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  random = ~as.factor(YrangeCat)+as.factor(sizeCat),
  spcov_type = "exponential",
  additive = "afvArea"
)
summary(ssn_modG) #no basins sig dif
anova(ssn_modG)
