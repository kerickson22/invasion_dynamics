library(rgbif)
library(sp)
library(maps)
library(maptools)
library(rgdal)
library(omnibus)
library(enmSdm)
library(bayesLopod)
library(calibrate)
library(scrubr)

setwd("C:/Users/kerickson/Documents/GitHub/invasion_dynamics")
load(".RData")

anacards_FL_XY <- anacards_FL[!is.na(anacards_FL$decimalLongitude), ]
anacards_FL_XY <- anacards_FL_XY[!is.na(anacards_FL_XY$decimalLatitude), ]
anacards_FL_XY <- anacards_FL_XY[anacards_FL_XY$decimalLongitude < -75, ]
anacards_FL_XY<-SpatialPoints(cbind(anacards_FL_XY$decimalLongitude, anacards_FL_XY$decimalLatitude), getCRS('wgs84', TRUE))



schinus_FL_XY <- schinus_FL[!is.na(schinus_FL$decimalLongitude), ]
schinus_FL_XY <- schinus_FL_XY[!is.na(schinus_FL_XY$decimalLatitude), ]
schinus_FL_XY <- schinus_FL_XY[schinus_FL_XY$decimalLongitude < -75, ]
schinus_FL_XY <- SpatialPoints(cbind(schinus_FL_XY$decimalLongitude, schinus_FL_XY$decimalLatitude), getCRS('wgs84', TRUE))

thing_schinus <- crop(schinus_FL_XY, rr)

elev <- getData("worldclim", var = "alt", res=10)
elev_FL <- mask(elev, florida)
elev_FL2 <- crop(elev, florida)
schinus_raster <- xyToRaster(schinus_FL_XY, anacards_FL_XY, basemap = elev_FL2)

spplot(schinus_raster, names.attr	=c("Species Detections","Sampling Effort"))

ld_Raster_adMatrix = rasterLopodData( rasterN = schinus_raster[["samplingEffort"]],
                                      rasterY = schinus_raster[["spDetections"]],
                                      Adjacency = T,
                                      extSample = 1.0,
                                      extDetection = 1.0 )

mLopodRaster = modelLopod( LopodData = ld_Raster_adMatrix,
                           varP = F,
                           q = NULL,
                           pmin = 0,
                           CAR = T,
                           nChains = 3,
                           warmup = 500,
                           sampling = 1000,
                           nCores = 4 )

lopodTrace( mLopodRaster, inc_warmup = T)
lopodDens(mLopodRaster, c("q", "p"))

ppRaster = lopodRaster(mLopodRaster, param = "pp", extrapolate = T, metric = "mean")
psiRaster = lopodRaster(mLopodRaster, param = "psi_i", extrapolate = T, metric = "mean")
qRaster = lopodRaster(mLopodRaster, param = "q", extrapolate = F, metric = "mean")

spplot( raster::stack(psiRaster,ppRaster),  
        names.attr = c("Occupancy (Psi)", "Probability of Presence"))


#getAlbersNA

getCRS("albersNA")
thing <- sp::spTransform(florida, getCRS("albersNA", TRUE))


#Read in records of all Tracheophytes in Florida from GBIF: 
tracheophyta <-read.delim("./Tracheophyta/occurrence.txt")
tracheophyta_XY <- SpatialPoints(cbind(tracheophyta$decimalLongitude, tracheophyta$decimalLatitude), getCRS('wgs84', TRUE))