library(rgbif)
library(sp)
library(maps)
library(maptools)
library(rgdal)
library(omnibus)
library(enmSdm)
library(bayesLopod)
library(calibrate)

setwd("H:/Global Change Program/Research/Occupancy Modeling of Messy Detections/Kelley Erickson/Schinus spread/Code")
schinus <- read.delim("H:/Global Change Program/Research/Occupancy Modeling of Messy Detections/Kelley Erickson/Schinus spread/Data/Schinus terebinthifolia records from GBIF/S_terebinthifolia records downloaded 2018-09-05 from GBIF/occurrence.txt")
anacards_US <-read.delim("H:/Global Change Program/Research/Occupancy Modeling of Messy Detections/Kelley Erickson/Schinus spread/Data/US Anacardiaceae records from GBIF/USA Anacardiaceae records from GBIF - downloaded 2018-09-05/occurrence.txt")
save(schinus, file="./Saved R Objects/00 Schinus records downloaded from GBIF on 2018-09-05.RData" )
save(anacards_US, file="./Saved R Objects/00 USA Anacardiaceae records downloaded from GBIF on 2018-09-05.RData")

counties <- rgdal::readOGR (
  'H:/Global Change Program/Research/Occupancy Modeling of Messy Detections/Kelley Erickson/Schinus spread/Data/GADM/ver2pt8/WGS84',
  'USA_adm2'
)

hawaii <- counties[counties@data$NAME_1=="Hawaii", ]
alaska <- counties[counties@data$NAME_1=="Alaska", ]
florida <- counties[counties@data$NAME_1=="Florida", ]
contiguousUSA <- counties[counties@data$NAME_1!= "Hawaii" & counties@data$NAME_1!="Alaska", ]

x <- data.frame(numSchinus=rep(0, nrow(contiguousUSA)), numAnacard=rep(0, nrow(contiguousUSA)))
contiguousUSA@data <- omnibus::insertCol(x, into=contiguousUSA@data, at='VARNAME_2', before=FALSE)
x2 <- data.frame(numSchinus=rep(0, nrow(florida)), numAnacard=rep(0, nrow(florida)))
florida@data <- omnibus::insertCol(x2, into=florida@data, at='VARNAME_2', before=F)

#Figure out how many records/time interval there are
dim(anacards_US[anacards_US$stateProvince == "Florida", ])
#There are 1274 records of Anacardiaceae in Florida 
#(Because sum(florida@data$numAnacard) = 1374, this means there are 
# 100 occurrences with latitude/longitude that fall within florida, but perhaps
# have stateProvince=NA)
anacards_FL <- anacards_US[anacards_US$stateProvince == "Florida", ]
anacards_FL_pre2000<-anacards_FL[anacards_FL$year<2000, ] #There are 991 records
anacards_FL_pre2000<-anacards_FL_pre2000[!is.na(anacards_FL_pre2000$year), ] #Now there are 960 records
#For these records, should not rely on reported latitude and longitude, but use state and county instead?

coords_anacards_FL_pre2000 <- anacards_FL_pre2000[!is.na(anacards_FL_pre2000$decimalLongitude) & !is.na(anacards_FL_pre2000$decimalLatitude), ]
coords_anacards_FL_pre2000 <- coords_anacards_FL_pre2000[!is.na(coords_anacards_FL_pre2000$decimalLongitude), ] 
coords_anacards_FL_pre2000 <- coords_anacards_FL_pre2000[!is.na(coords_anacards_FL_pre2000$decimalLatitude), ]
coords_anacards_FL_pre2000 <- SpatialPoints(cbind(coords_anacards_FL_pre2000$decimalLongitude, coords_anacards_FL_pre2000$decimalLatitude), getCRS('wgs84', TRUE))


#Do the same thing for Schinus: 

#Figure out how many records/time interval there are
dim(schinus[schinus$stateProvince == "Florida", ])
#There are 338 records of Schinus in Florida 
#(Because sum(florida@data$numAnacard) = 431, this means there are 
# 93 occurrences with latitude/longitude that fall within florida, but perhaps
# have stateProvince=NA)
schinus_FL <- schinus[schinus$stateProvince == "Florida", ]
schinus_FL_pre2000<-schinus_FL[schinus_FL$year<2000, ] #There are 272 records
schinus_FL_pre2000<-schinus_FL_pre2000[!is.na(schinus_FL_pre2000$year), ] #Now there are 268 records
#For these records, should not rely on reported latitude and longitude, but use state and county instead?

coords_schinus_FL_pre2000 <- schinus_FL_pre2000[!is.na(schinus_FL_pre2000$decimalLongitude) & !is.na(schinus_FL_pre2000$decimalLatitude), ]
coords_schinus_FL_pre2000 <- coords_schinus_FL_pre2000[!is.na(coords_schinus_FL_pre2000$decimalLongitude), ] 
coords_schinus_FL_pre2000 <- coords_schinus_FL_pre2000[!is.na(coords_schinus_FL_pre2000$decimalLatitude), ]
coords_schinus_FL_pre2000 <- SpatialPoints(cbind(coords_schinus_FL_pre2000$decimalLongitude, coords_schinus_FL_pre2000$decimalLatitude), getCRS('wgs84', TRUE))


thing<-schinus_FL_pre2000 #268 items
thing<-thing[!is.na(thing$decimalLongitude),] #24 items 
thing<-thing[!is.na(thing$decimalLatitude), ] #24 items
thing<-thing[!is.na(thing$year), ] #24 items
textxy(thing$decimalLongitude, thing$decimalLatitude, thing$year, cex=.8)

###What about records after 2000? 
dim(schinus_FL[schinus_FL$year>=2000 & !is.na(schinus_FL$year), ]) #66 
dim(schinus_FL[schinus_FL$year>=2010 & !is.na(schinus_FL$year), ]) #24
dim(schinus_FL[schinus_FL$year<2010 & schinus_FL$year>=2000 & !is.na(schinus_FL$year),] ) #42

schinus_FL_2000_2010 <- schinus_FL[schinus_FL$year<2010 & schinus_FL$year>=2000 & !is.na(schinus_FL$year),]
schinus_FL_2000_2010 <- schinus_FL_2000_2010[!is.na(schinus_FL_2000_2010$decimalLatitude), ]
schinus_FL_2000_2010 <- schinus_FL_2000_2010[!is.na(schinus_FL_2000_2010$decimalLongitude), ]
schinus_FL_2000_2010 <- schinus_FL_2000_2010[!is.na(schinus_FL_2000_2010$year), ] 

coords_schinus_FL_2000_2010 <- SpatialPoints(cbind(schinus_FL_2000_2010$decimalLongitude, 
                                                   schinus_FL_2000_2010$decimalLatitude),
                                                   getCRS('wgs84', TRUE))

points(coords_schinus_FL_2000_2010, pch=19, col="blue", cex=1.5)
textxy(schinus_FL_2000_2010$decimalLongitude, schinus_FL_2000_2010$decimalLatitude, 
        schinus_FL_2000_2010$year, col="blue", cex=0.8)


schinus_FL_2010_2018 <- schinus_FL[schinus_FL$year>2010 & !is.na(schinus_FL$year), ]
schinus_FL_2010_2018 <- schinus_FL_2010_2018[!is.na(schinus_FL_2010_2018$decimalLatitude), ]
schinus_FL_2010_2018 <- schinus_FL_2010_2018[!is.na(schinus_FL_2010_2018$decimalLongitude), ]
schinus_FL_2010_2018 <- schinus_FL_2010_2018[!is.na(schinus_FL_2010_2018$year), ]

coords_schinus_FL_2010_2018 <- SpatialPoints(cbind(schinus_FL_2010_2018$decimalLongitude, 
                                                   schinus_FL_2010_2018$decimalLatitude),
                                                   getCRS('wgs84', TRUE))

points(coords_schinus_FL_2010_2018, pch=19, col="darkgreen", cex=1.5)
textxy(schinus_FL_2010_2018$decimalLongitude, schinus_FL_2010_2018$decimalLatitude, 
       schinus_FL_2010_2018$year, col="darkgreen", cex=0.8)

### count number of records WITH coordinates in each county
###########################################################

# clean records for Anacardiaceae 
coordsAnacard <- anacards_US[!is.na(anacards_US$decimalLongitude) & !is.na(anacards_US$decimalLatitude), ]
coordsAnacard <- coordsAnacard[!is.na(coordsAnacard$decimalLongitude), ] 
coordsAnacard <- coordsAnacard[!is.na(coordsAnacard$decimalLatitude), ]

# clean records for just Schinus 
coordsSchinus <- schinus[!is.na(schinus$decimalLongitude) & !is.na(schinus$decimalLatitude), ]
coordsSchinus <- coordsSchinus[!is.na(coordsSchinus$decimalLongitude), ] 
coordsSchinus <- coordsSchinus[!is.na(coordsSchinus$decimalLatitude), ]


# convert to spatial object
coordsAnacard <- SpatialPoints(cbind(coordsAnacard$decimalLongitude, coordsAnacard$decimalLatitude), getCRS('wgs84', TRUE))
coordsSchinus <- SpatialPoints(cbind(coordsSchinus$decimalLongitude, coordsSchinus$decimalLatitude), getCRS('wgs84', TRUE))

# intersect coordinates with counties
AnacardInCounties <- over(coordsAnacard, contiguousUSA)
SchinusInCounties <- over(coordsSchinus, contiguousUSA)
AnacardInCounties <- AnacardInCounties[!is.na(AnacardInCounties$NAME_2), ]
SchinusInCounties <- SchinusInCounties[!is.na(SchinusInCounties$NAME_2), ]

AnacardInFlorida <- over(coordsAnacard, florida) 
AnacardInFlorida <- AnacardInFlorida[!is.na(AnacardInFlorida$NAME_2), ]

SchinusInFlorida <- over(coordsSchinus, florida)
SchinusInFlorida <- SchinusInFlorida[!is.na(SchinusInFlorida$NAME_2), ]
# find number of records in each county (do by county)
#For all Anacardiaceae 
for (i in 1:nrow(AnacardInFlorida)) {
  
  thisRow <- which(florida$NAME_1 == AnacardInFlorida$NAME_1[i] & florida$NAME_2 == AnacardInFlorida$NAME_2[i])
  florida@data$numAnacard[thisRow] <- florida@data$numAnacard[thisRow] + 1
  
}

#For all Schinus
for (i in 1:nrow(SchinusInFlorida)) {
  
  thisRow <- which(florida$NAME_1 == SchinusInFlorida$NAME_1[i] & florida$NAME_2 == SchinusInFlorida$NAME_2[i])
  florida@data$numSchinus[thisRow] <- florida@data$numSchinus[thisRow] + 1
  
}

### count number of records WITHOUT coordinates in each county
#############################################################

# clean records
noCoord_anacard <- anacards_US[is.na(anacards_US$decimalLongitude) & is.na(anacards_US$decimalLatitude) & anacards_US$hasGeospatialIssues=='false' & !is.na(anacards_US$county) & !is.na(anacards_US$stateProvince), ]
noCoord_schinus <- schinus[is.na(schinus$decimalLongitude) & is.na(schinus$decimalLatitude) & schinus$hasGeospatialIssues=='false' & !is.na(schinus$county) & !is.na(schinus$stateProvince), ]

#restrict to US records 
noCoord_anacard <- noCoord_anacard[noCoord_anacard$countryCode=="US", ] #shouldn't change anything
noCoord_schinus <- noCoord_schinus[noCoord_schinus$countryCode=="US", ]



#Anacardiaceae first: 
if (nrow(noCoord_anacard) > 0) {
  
  # ensure county names match GADM
  noCoord_anacard$county <- gsub(noCoord_anacard$county, pattern='(', replacement='', fixed=T)
  noCoord_anacard$county <- gsub(noCoord_anacard$county, pattern=')', replacement='', fixed=T)
  noCoord_anacard$county <- gsub(noCoord_anacard$county, pattern=' County', replacement='')
  noCoord_anacard$county <- gsub(noCoord_anacard$county, pattern='County of ', replacement='')
  noCoord_anacard$county <- gsub(noCoord_anacard$county, pattern=' Cty.', replacement='')
  noCoord_anacard$county <- gsub(noCoord_anacard$county, pattern=' Cty', replacement='')
  noCoord_anacard$county <- gsub(noCoord_anacard$county, pattern=' Co.', replacement='')
  noCoord_anacard$county <- gsub(noCoord_anacard$county, pattern='Ste ', replacement='Saint ')
  noCoord_anacard$county <- gsub(noCoord_anacard$county, pattern='St. ', replacement='Saint ')
  
  
  # add each record to county in which it falls--match state/province AND county

  for (i in 1:nrow(noCoord_anacard)) {
    
    thisRow <- which(as.character(florida$NAME_1) == as.character(noCoord_anacard$stateProvince[i]) & as.character(florida$NAME_2) == as.character(noCoord_anacard$county[i]))
    florida@data$numAnacard[thisRow] <- florida@data$numAnacard[thisRow] + 1
    
  }
  
}




#Now Schinus 
if (nrow(noCoord_schinus) > 0) {
  
  # ensure county names match GADM
  noCoord_schinus$county <- gsub(noCoord_schinus$county, pattern=' County', replacement='')
  noCoord_schinus$county <- gsub(noCoord_schinus$county, pattern='(', replacement='', fixed=T)
  noCoord_schinus$county <- gsub(noCoord_schinus$county, pattern=')', replacement='', fixed=T)
  noCoord_schinus$county <- gsub(noCoord_schinus$county, pattern='County of ', replacement='')
  noCoord_schinus$county <- gsub(noCoord_schinus$county, pattern=' Cty.', replacement='')
  noCoord_schinus$county <- gsub(noCoord_schinus$county, pattern=' Cty', replacement='')
  noCoord_schinus$county <- gsub(noCoord_schinus$county, pattern=' Co.', replacement='')
  noCoord_schinus$county <- gsub(noCoord_schinus$county, pattern='Ste ', replacement='Saint ')
  noCoord_schinus$county <- gsub(noCoord_schinus$county, pattern='St. ', replacement='Saint ')
  
  
  # add each record to county in which it falls--match state/province AND county
  
  for (i in 1:nrow(noCoord_schinus)) {
    
    thisRow <- which(as.character(florida$NAME_1) == as.character(noCoord_schinus$stateProvince[i]) & as.character(florida$NAME_2) == as.character(noCoord_schinus$county[i]))
    florida@data$numSchinus[thisRow] <- florida@data$numSchinus[thisRow] + 1
    
  }
  
}
#Had to put in as.character() as the number of levels of counties and states in the shape file is different
# from the number of levels of counties and states in the GBIF-missing-coordinates data

ld_Shape = shapeLopodData( Shapefile = florida,
                           fieldN = "numAnacard",
                           fieldY = "numSchinus",  
                           Adjacency = T,
                           keepFields = F)

schinusSEff = spplot(florida, zcol = c("numAnacard"), main = "Sampling Effort", border=NA)
schinusDetect = spplot(florida, zcol = c("numSchinus"), main = "Detections", border=NA)
plot(schinusSEff, split = c(1,1,2,1), more = T, border=NA)
plot(schinusDetect,  split = c(2,1,2,1), more = F, border=NA)


mLopodShape = modelLopod(LopodData = ld_Shape,
                          varP = T,
                          q = NULL,
                          pmin = 0,
                          CAR = T,
                          nChains = 3,
                          warmup = 400,
                          sampling = 100,
                          nCores =3)

#When I run it with CAR=T, I get the following error message: 
# "Error in new_CppObject_xp(fields$.module, fields$.pointer, ...) : 
# Exception: mismatch in number dimensions declared and found in context;
# processing stage=data initialization; variable name=notSampledId; 
# dims declared=(1); dims found=()  (in 'model_psiipiq_CAR' at line 44)

# Is it because county #3 has 0 sampling effort? 
ld_Shape@geoDataObject$sampEffort[3]
#(originally it was 0)
#YES. 

lopodTrace( mLopodShape, inc_warmup = T)
lopodDens(mLopodShape, c("q", "pmin", "pmax"))

SchinusShape = mLopodShape@LopodData@geoDataObject

psiShape = lopodShape(mLopodShape, "psi_i", extrapolate = T,  quant = 0.05)
SchinusShape@data[,"psi_05"] = psiShape@data[,"psi_i"]

psiShape = lopodShape(mLopodShape, "psi_i", extrapolate = T,  quant = 0.5)
SchinusShape@data[,"psi_50"] = psiShape@data[,"psi_i"]

psiShape = lopodShape(mLopodShape, "psi_i", extrapolate = T,  quant = 0.95)
SchinusShape@data[,"psi_95"] = psiShape@data[,"psi_i"]


spplot( SchinusShape,
        zcol = c("psi_05", "psi_50", "psi_95"),
        names.attr = c("Psi (5% quantile)", "Psi (median)", "Psi (95% quantile)"), main = "Occupancy (Psi)")

#Next time: figure out where CAR model is going wrong
