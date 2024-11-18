################################################################
# R script for SSNbler_vignette.pdf ####
#
# Includes additional comments and R code not included in the
# SSNbler vignete
################################################################
library(SSNbler)
library(SSN2)
library(caret)
library(tidyverse)
# ------------- Import and view the input data ----------------

## Copy the example dataset to a temporary directory
#copy_streams_to_temp()
#path <- paste0(tempdir(), "/streamsdata")

## List local input files
#list.files(path)

## Load the sf package and import the streams, observation sites, and
## two prediction datasets
library(sf)
streams <- st_read("NSI_fix6.shp")
obs <- st_read("WILDCARD_SNAPPED.shp")

streams=st_cast(streams, to="LINESTRING")

streams <- st_transform(streams, crs = 5070)
obs <- st_transform(obs, crs = 5070)


## Look at column names
names(obs)

## Plot the data using ggplot2
ggplot() +
  geom_sf(data = streams) +
  geom_sf(data = obs, color = "blue", size = 2) +
  coord_sf(datum = st_crs(streams))

# --------------- Build the LSN -------------------------------

## Set the lsn.path variable
lsn.path <- "persist_LSN/"

## Build the LSN and take note of the messages printed to the console

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
View(nssn$obs)

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
View(obs.distmat[[1]])

## Create symmetric hydrologic distance matrix
obs.distmat2 <- obs.distmat[[1]] + t(obs.distmat[[1]])
View(obs.distmat2)








############Create Persistence dataset [Left here 11-17-2024]
















# ---- Fit spatial stream-network model --------------------------

## Fit a spatial linear model to Summer mean temperature with a
## mixture of TU/TD/EUC covariance models

pairs(cfs ~ bfi+air+pr+gw, data = niob_ssn$obs)
#niob_ssn$obs$Agency[which(is.na(niob_ssn$obs$Agency))]="None"
ssn_mod <- ssn_lm(
  formula = sttemp ~ cfs+bfi+air+pr+gw,
  ssn.object = niob_ssn,
  tailup_type = "mariah",
  taildown_type = "mariah",
  euclid_type = "exponential",
  random = ~ as.factor(Year),
  additive = "afvArea", estmethod = "ml"
)

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

## Make predictions at all locations in the pred1km data
library(tidyverse)
names(niob_ssn$preds$pred1km)[names(niob_ssn$preds$pred1km) == 'air9311'] <- 'air'
names(niob_ssn$preds$pred1km)[names(niob_ssn$preds$pred1km) == 'pr9311'] <- 'pr'
## Recreate model with updates names for air and precip
ssn_mod <- ssn_lm(
  formula = sttemp ~ cfs+bfi+air+pr+gw,
  ssn.object = niob_ssn,
  tailup_type = "mariah",
  taildown_type = "mariah",
  euclid_type = "exponential",
  #random = ~ as.factor(Year),
  additive = "afvArea", estmethod = "reml"
)

predict(ssn_mod, newdata = "pred1km")

## Augment prediction data with predictions
aug_preds <- augment(ssn_mod, newdata = "pred1km")
aug_preds[, ".fitted"]


## Visualize predictions on the network
ggplot() +
  geom_sf(data = aug_preds, aes(color = .fitted), size = 2) +
  scale_color_viridis_b()

## Write out augmented prediction data to a geopackage
st_write(aug_preds, paste0("aug_preds_93-11.gpkg"))

