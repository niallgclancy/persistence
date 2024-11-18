################################################################
# R script for SSNbler_vignette.pdf ####
#
# Includes additional comments and R code not included in the
# SSNbler vignete
################################################################

# ------- Install SSNbler (if necessary) and load library ------


## Install package from CRAN repository 
## install.packages("SSNbler")

## Install package from Github
##install.packages("remotes") ## Install remotes if needed
##library(remotes)

## Install the latest version of SSNbler from github
##remotes::install_github("pet221/SSNbler", ref = "main")

## Load SSNbler library
#install.packages("SSNbler")
#install.packages("SSN2")
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
streams <- st_read("nsi.shp")
obs <- st_read("observations.shp")
preds <- st_read("predictions.shp")

streams <- st_transform(streams, crs = 5070)
obs <- st_transform(obs, crs = 5070)
preds <- st_transform(preds, crs = 5070)




## Notice that the imported data are of class sf data.frame
class(obs)

## Look at column names
names(obs)
help(MF_obs)

## Plot the data using ggplot2
ggplot() +
  geom_sf(data = streams) +
  geom_sf(data = preds, colour = "purple", size = 1.7) +
  geom_sf(data = obs, color = "blue", size = 2) +
  coord_sf(datum = st_crs(streams))

# --------------- Build the LSN -------------------------------

## Set the lsn.path variable
lsn.path <- "nebraska_LSN/"

## Build the LSN and take note of the messages printed to the console
edges <- lines_to_lsn(
  streams = streams,
  lsn_path = lsn.path,
  check_topology = TRUE,
  snap_tolerance = 0.05,
  topo_tolerance = 20,
  overwrite = TRUE
)

## list new files saved in lsn.path
list.files(lsn.path)

# ----------- Incorporate sites into LSN -----------------------

## Incorporate observations. Pay attention to output messages in R
## console to ensure all sites are snapped successfully.
obs <- sites_to_lsn(
  sites = obs,
  edges = edges,
  lsn_path = lsn.path,
  file_name = "obs",
  snap_tolerance = 200,
  save_local = TRUE,
  overwrite = TRUE
)

## Predictions 
preds <- sites_to_lsn(
  sites = preds,
  edges = edges,
  save_local = TRUE,
  lsn_path = lsn.path,
  file_name = "preds.gpkg",
  snap_tolerance = 10,
  overwrite = TRUE
)



## list files saved in lsn.path. Notice that geopackages 
## for the sites have been added locally
list.files(lsn.path)

## Notice the new columns rid, ratio, and snapdist
names(obs)

## Summarise snapdist, looking for sites that were snapped large
## distances
summary(obs$snapdist)

## Look at ratio values to ensure they range from 0-1
summary(obs$ratio)

# --------- Calculate upstream distance for edges ------------
edges$GNIS_ID=as.numeric(edges$GNIS_ID)
edges$REACHCODE=as.numeric(edges$REACHCODE)

edges <- updist_edges(
  edges = edges,
  save_local = TRUE,
  lsn_path = lsn.path,
  calc_length = TRUE
)

## View edges column names
names(edges)

## Summarise upDist
summary(edges$upDist)

# ---------- Calculate upstream distance for sites -----------

site.list <- updist_sites(
  sites = list(
    obs = obs,      ## Input is a named list of sf objects
    pred1km = preds),
  edges = edges,
  length_col = "Length",
  save_local = TRUE,
  lsn_path = lsn.path
)

## Notice that length of site.list = 2
length(site.list)

## View output site.list names
names(site.list)

## View column names in obs
names(site.list$obs)       

## Plot the upstream distances for edges and obs
ggplot() +
  geom_sf(data = edges, aes(color = upDist)) +
  geom_sf(data = site.list$obs, aes(color = upDist)) +
  coord_sf(datum = st_crs(MF_streams)) +
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

## Look at edges column names. Notice the two new columns, 
## areaPI and afvArea
names(edges)

## Summarize the AFV column to ensure it ranges from 0-1
summary(edges$afvArea)

# ------- Calculate AFV for sites -------------------------------

site.list <- afv_sites(
  sites = site.list,    ## Input is a named list
  edges = edges,
  afv_col = "afvArea",
  save_local = TRUE,
  lsn_path = lsn.path
)

## View column names in pred1km. Notice the new afvArea column
names(site.list$pred1km)

## Summarise AFVs in pred1km to ensure they range from 0-1
summary(site.list$pred1km$afvArea)

# ---------- Assemble the SSN Object -----------------------------

niob_ssn <- ssn_assemble(
  edges = edges,
  lsn_path = lsn.path,
  obs_sites = site.list$obs,          ## Input is sf data.frame
  preds_list = site.list["pred1km"],
  ssn_path = "niobssn.ssn",
  import = TRUE,
  check = TRUE,
  afv_col = "afvArea",
  overwrite = TRUE
)

## Look at files
list.files("niobssn.ssn")

## Check class for mf04_ssn
class(niob_ssn)

## Get SSN object names
names(NEMont_ssn)

## Print path to local .ssn
NEMont_ssn$path

## Print names of prediction datasets in SSN object
names(NEMont_ssn$preds)

## Look at the observations in the SSN object. Notice the new columns
## added when the SSN was assembled
class(NEMont_ssn$obs)
names(NEMont_ssn$obs)
View(NEMont_ssn$obs)

## Plot ssn
ggplot() +
  geom_sf(
    data = niob_ssn$edges,
    color = "medium blue",
    aes(linewidth = TotDASqKM)
  ) +
  #scale_linewidth(range = c(0.1, 2.5)) +
  geom_sf(
    data = niob_ssn$preds$pred1km,
    size = 1.5,
    shape = 21,
    fill = "white",
    color = "dark grey"
  ) +
  geom_sf(
    data = niob_ssn$obs,
    size = 1.7,
    aes(color = sttemp)
  ) +
  coord_sf(datum = st_crs(streams)) +
  scale_color_viridis_c() +
  labs(color = "Temperature", linewidth = "WS Area") +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )




# ---- Create Distance Matrices ----------------------------------
library(SSN2)

## Generate hydrologic distance matrices
ssn_create_distmat(niob_ssn, predpts = "pred1km")

## Look at distance matrix files for obs
#list.files(paste0(niob_ssn$path, "/distance/obs"))

## Look at distance matrix files for pred1km
#list.files(paste0(mf04_ssn$path, "/distance/pred1km"))

## Output a list of distance matrices for obs
obs.distmat<- ssn_get_stream_distmat(niob_ssn,
                                     name = "obs")
names(obs.distmat)
colnames(obs.distmat[[1]])
View(obs.distmat[[1]])

## Create symmetric hydrologic distance matrix
obs.distmat2 <- obs.distmat[[1]] + t(obs.distmat[[1]])
View(obs.distmat2)

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

                              