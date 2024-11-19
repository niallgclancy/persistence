################################################################
# R script for Clancy Dissertation Chapter on Empirical Characteristics of Refugia



################################################################
library(SSNbler)
library(SSN2)
library(caret)
library(tidyverse)
library(sf)


## import the streams, observation sites
streams <- st_read("NSI_fix6.shp")
obs <- st_read("PersistenceMetrics.shp")

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





# ---- Fit spatial stream-network model --------------------------
## Fit a spatial linear model to Summer mean temperature with a
## mixture of TU/TD/EUC covariance models
library(PerformanceAnalytics)


#----------------Identify predictors with multicolinearity 
mdata=obs[,c(15,22,23, 24,25,127)]
mdata=as.data.frame(mdata)
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)
#None of the selected predictors (for data exploration) are highly correlated (max r = 0.28)

#--------------------Beta regression can't use 0 or 1...convert to near 0 and 1--------------------
nssn$obs$persNative=as.numeric(nssn$obs$persNative)
nssn2=nssn

#Native Sp. Persistence
nssn2$obs$persNative=as.numeric(nssn2$obs$persNative)
nssn2$obs$persNative[which(nssn2$obs$persNative==0)]=0.000000001
nssn2$obs$persNative[which(nssn2$obs$persNative==1)]=0.999999999

#Full Community Persistence
nssn2$obs$persComm=as.numeric(nssn2$obs$persComm)
nssn2$obs$persComm[which(nssn2$obs$persComm==0)]=0.000000001
nssn2$obs$persComm[which(nssn2$obs$persComm==1)]=0.999999999
anyNA(nssn2$obs$persComm)


hist(nssn2$obs$persComm)


#Individual Species Test
#nssn2$obs$LKCH=as.numeric(nssn2$obs$LKCH)
#nssn2$obs$LKCH[which(nssn2$obs$LKCH==1)]=NA
#nssn2$obs$LKCH[which(nssn2$obs$LKCH==0)]=1
#nssn2$obs$LKCH[which(nssn2$obs$LKCH==-1)]=0


###############################Beta Regression for proportion metrics not working 11-19-2024

##Scale covariates to determine variable importance
ssn_mod <- ssn_glm(
  formula =   persComm~ scale(S1_93_11)+scale(F_MAUG_HIS)+scale(Yrange)+scale(DSbarrier)+scale(Length_km)+scale(nnPreds),
  family = "beta",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  additive = "afvArea", estmethod = "ml"
)
##Non-scaled version for comparison 
ssn_mod <- ssn_glm(
  formula =   persComm~ S1_93_11+F_MAUG_HIS+Yrange+DSbarrier+Length_km+nnPreds,
  family = "beta",
  ssn.object = nssn2,
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  additive = "afvArea", estmethod = "ml"
)










#reg.glm=glm(persComm ~ S1_93_11+F_MAUG_HIS+Yrange+DSbarrier+Length_km+nnPreds, data = nssn2$obs,
 #   family = "binomial")

#summary(reg.glm)


## Summarise model results
summary(ssn_mod)
loocv_mod <- loocv(ssn_mod)
loocv_mod$RMSPE


tg <- Torgegram(
  formula = sttemp ~ cfs+bfi+air+pr+gw,
  ssn.object = niob_ssn,
  type = c("flowcon", "flowuncon", "euclid")
)

## Visualize Torgegram
plot(tg)


## Summarize the model
summary(ssn_mod)

## Inspect the variance components
varcomp(ssn_mod)

## Tidy the model output
tidy(ssn_mod, conf.int = TRUE)

## Glance at the model fit
glance(ssn_mod)



## Fit a model with just nugget (equivalent to lm)
ssn_mod3 <- ssn_lm(
  formula = meanAug ~ AirTemp_mA+PrecipHIST*ELEV+CUMDRAINAG,
  ssn.object = NEMont_ssn,estmethod = "ml"
)

## Glance at all model fits (look for lowest AIC)
glances(ssn_mod, ssn_mod3)

## leave-one-out cross validation for each model; find RMSPE
loocv_mod <- loocv(ssn_mod)
loocv_mod
loocv_mod$RMSPE

loocv_mod3 <- loocv(ssn_mod3)
loocv_mod3$RMSPE


## glance at model fit and leave-one-out
glances(ml_mod, ml_mod2)
loocv_mod_ml <- loocv(ml_mod)
loocv_mod_ml$RMSPE
loocv_mod_ml2 <- loocv(ml_mod2)
loocv_mod_ml2$RMSPE

## Refit final model using REML
ssn_mod <- ssn_lm(
  formula = sttemp ~ cfs+bfi+air+pr+gw,
  ssn.object = niob_ssn,
  tailup_type = "mariah",
  taildown_type = "mariah",
  euclid_type = "exponential",
  random = ~ as.factor(Year),
  additive = "afvArea", estmethod = "reml"
)


## Summarise model results
summary(ssn_mod)
loocv_mod <- loocv(ssn_mod)
loocv_mod
loocv_mod$RMSPE



#:::::::::::::::::::::::::::::::::::::
# Prediction for historical period (1993-2011)####
#:::::::::::::::::::::::::::::::::::::