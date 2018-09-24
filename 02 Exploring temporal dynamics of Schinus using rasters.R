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
tracheophyta2 <- read.csv("./Tracheophyta/occurrence2.csv", header=T)
#To create tracheophyta2, I opened the txt file in Excel and saved it as a csv. 
# It seems like the csv file has retained more of the records without introducing errors 



#Check: 
summary(tracheophyta2$stateProvince)
#When I run this, there are 128,525 with state Florida, 1 with state florida, and 69 with random strings of characters...
# This is data that I downloaded from GBIF directly (I did not use rgbif)

tracheophyta2 <- subset(tracheophyta2, tracheophyta2$stateProvince == "florida" | tracheophyta2$stateProvince=="Florida")

tracheophyta2$stateProvince <- gsub(tracheophyta2$stateProvince, pattern='florida', replacement='Florida', fixed=T)
tracheophyta2$stateProvince <- factor(tracheophyta2$stateProvince)

coords_tracheophyta <- subset(tracheophyta2, tracheophyta2$decimalLatitude != '')
coords_tracheophyta <- subset(coords_tracheophyta, coords_tracheophyta$decimalLongitude != '')

tracheophyta_XY <- SpatialPoints(cbind(coords_tracheophyta$decimalLongitude, coords_tracheophyta$decimalLatitude), getCRS('wgs84', TRUE))




#Add a column for tracheophytes to the florida shapefile 
x <- data.frame(numTracheo=rep(0, nrow(florida)))
florida@data <- omnibus::insertCol(x, into=florida@data, at='numAnacard', before=FALSE)


# intersect coordinates with counties


TracheoInFlorida <- over(tracheophyta_XY, florida) 
TracheoInFlorida <- TracheoInFlorida[!is.na(TracheoInFlorida$NAME_2), ]


# find number of records in each county (do by county)
#For all Tracheophyta
for (i in 1:nrow(TracheoInFlorida)) {
  
  thisRow <- which(florida$NAME_1 == TracheoInFlorida$NAME_1[i] & florida$NAME_2 == TracheoInFlorida$NAME_2[i])
  florida@data$numTracheo[thisRow] <- florida@data$numTracheo[thisRow] + 1
  
}
