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
library(plyr)
library(purrr)
library(magick)

setwd("C:/Users/kerickson/Documents/GitHub/invasion_dynamics")
load(".RData")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
memory.limit(memory.limit() *2^30)


#All US occurrences of anacardiaceae 
anacardeaceae <- read.csv("./US Anacardiaceae records from GBIF/anacardiaceae_occurrences.csv", header=T)
#There are 43, 947 records of Anacardeacea in the United States 
#Global occurrences of schinus
schinus <- read.csv("./Schinus records from GBIF/schinus_records.csv", header=T)
# There are 7,240 global records of S. terebinthifolia
tracheophyta<- read.csv("./Tracheophyta/occurrence2.csv", header=T)
#There are 383,697 Florida records with coordinates of Tracheophytes 

#Anacardiaceae with coordinates: 
anacardeaceae_coords <- anacardeaceae[!is.na(anacardeaceae$decimalLongitude), ]
anacardeaceae_coords <- anacardeaceae_coords[!is.na(anacardeaceae_coords$decimalLatitude), ]
#There are 33,965 US Anacardeaceae records with coordinates 
#Anacardiaceae without coordinates: 
anacardeaceae_nocoords <- anacardeaceae[is.na(anacardeaceae$decimalLongitude) | is.na(anacardeaceae$decimalLatitude), ]
#To restrict to Florida, have to clean up the stateProvince codes first: 
anacardeaceae_nocoords$stateProvince <- gsub(anacardeaceae_nocoords$stateProvince, pattern='Miami, Florida', replacement='Florida')
anacardeaceae_nocoords$stateProvince <- gsub(anacardeaceae_nocoords$stateProvince, pattern='Flórida', replacement='Florida')
anacardeaceae_nocoords <- anacardeaceae_nocoords[anacardeaceae_nocoords$stateProvince =="Florida", ]
#There are 969 FL Anacardeaceae records without coordinates 


#Schinus with coordinates: 
schinus_coords <- schinus[!is.na(schinus$decimalLongitude), ]
schinus_coords <- schinus_coords[!is.na(schinus$decimalLatitude), ]
# There are 4,567 Schinus observations with coordinates 
#Schinus without coordinates:
schinus_nocoords <- schinus[is.na(schinus$decimalLongitude)| is.na(schinus$decimalLatitude), ]
#There are 2673 Schinus observations with no coordinates 
#Restrict to US: 
schinus_nocoords <-schinus_nocoords[schinus_nocoords$countryCode =="US", ]
#There are 396 US Schinus records with no coordinates  
#Restrict to Florida: 
schinus_nocoords <- schinus_nocoords[schinus_nocoords$stateProvince =="Florida", ]
#There are 281 Florida Schinus records with no coordinates 
#Remove those records without a county
schinus_nocoords <- schinus_nocoords[!is.na(schinus_nocoords$county) & schinus_nocoords$county!="", ]
#There are 268 Florida Schinus records with no coordinates but with a county listed 



#Divide into year bins
schinus_coords_years<-schinus_coords[!is.na(schinus_coords$year), ]
schinus_nocoords_years <- schinus_nocoords[!is.na(schinus_nocoords$year), ]
#There are 2565 + 266 records with years (not all are in FL necessarily)
anacardeaceae_coords_years <- anacardeaceae_coords[!is.na(anacardeaceae_coords$year), ]
anacardeaceae_nocoords_years <- anacardeaceae_nocoords[!is.na(anacardeaceae_nocoords$year), ]
years<-c(schinus_coords_years$year, schinus_nocoords_years$year)
#Range of years 1876-2018 

#(1870, 1880]
#(1880, 1890]
#(1890, 1900]



categories <- seq(1880, 2020, by=10)
schinus_coords_1880 <- schinus_coords_years[schinus_coords_years$year<1880, ]
schinus_coords_1890 <- schinus_coords_years[schinus_coords_years$year<1890, ]
schinus_coords_1900 <- schinus_coords_years[schinus_coords_years$year<1900, ]
schinus_coords_1910 <- schinus_coords_years[schinus_coords_years$year<1910, ]
schinus_coords_1920 <- schinus_coords_years[schinus_coords_years$year<1920, ]
schinus_coords_1930 <- schinus_coords_years[schinus_coords_years$year<1930, ]
schinus_coords_1940 <- schinus_coords_years[schinus_coords_years$year<1940, ]
schinus_coords_1950 <- schinus_coords_years[schinus_coords_years$year<1950, ]
schinus_coords_1960 <- schinus_coords_years[schinus_coords_years$year<1960, ]
schinus_coords_1970 <- schinus_coords_years[schinus_coords_years$year<1970, ]
schinus_coords_1980 <- schinus_coords_years[schinus_coords_years$year<1980, ]
schinus_coords_1990 <- schinus_coords_years[schinus_coords_years$year<1990, ]
schinus_coords_2000 <- schinus_coords_years[schinus_coords_years$year<2000, ]
schinus_coords_2010 <- schinus_coords_years[schinus_coords_years$year<2010, ]
schinus_coords_2020 <- schinus_coords_years[schinus_coords_years$year<2020, ]


schinus_XY_1880<-SpatialPoints(cbind(schinus_coords_1880$decimalLongitude, schinus_coords_1880$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_1890<-SpatialPoints(cbind(schinus_coords_1890$decimalLongitude, schinus_coords_1890$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_1900<-SpatialPoints(cbind(schinus_coords_1900$decimalLongitude, schinus_coords_1900$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_1910<-SpatialPoints(cbind(schinus_coords_1910$decimalLongitude, schinus_coords_1910$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_1920<-SpatialPoints(cbind(schinus_coords_1920$decimalLongitude, schinus_coords_1920$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_1930<-SpatialPoints(cbind(schinus_coords_1930$decimalLongitude, schinus_coords_1930$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_1940<-SpatialPoints(cbind(schinus_coords_1940$decimalLongitude, schinus_coords_1940$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_1950<-SpatialPoints(cbind(schinus_coords_1950$decimalLongitude, schinus_coords_1950$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_1960<-SpatialPoints(cbind(schinus_coords_1960$decimalLongitude, schinus_coords_1960$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_1970<-SpatialPoints(cbind(schinus_coords_1970$decimalLongitude, schinus_coords_1970$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_1980<-SpatialPoints(cbind(schinus_coords_1980$decimalLongitude, schinus_coords_1980$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_1990<-SpatialPoints(cbind(schinus_coords_1990$decimalLongitude, schinus_coords_1990$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_2000<-SpatialPoints(cbind(schinus_coords_2000$decimalLongitude, schinus_coords_2000$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_2010<-SpatialPoints(cbind(schinus_coords_2010$decimalLongitude, schinus_coords_2010$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_2020<-SpatialPoints(cbind(schinus_coords_2020$decimalLongitude, schinus_coords_2020$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_df <- list(schinus_XY_1880, schinus_XY_1890, schinus_XY_1900, schinus_XY_1910, schinus_XY_1920,
                      schinus_XY_1930, schinus_XY_1940, schinus_XY_1950, schinus_XY_1960, schinus_XY_1970, 
                      schinus_XY_1980, schinus_XY_1990, schinus_XY_2000, schinus_XY_2010, schinus_XY_2020)


schinus_nocoords_1880 <- schinus_nocoords_years[schinus_nocoords_years$year<1880, ]
schinus_nocoords_1890 <- schinus_nocoords_years[schinus_nocoords_years$year<1890, ]
schinus_nocoords_1900 <- schinus_nocoords_years[schinus_nocoords_years$year<1900, ]
schinus_nocoords_1910 <- schinus_nocoords_years[schinus_nocoords_years$year<1910, ]
schinus_nocoords_1920 <- schinus_nocoords_years[schinus_nocoords_years$year<1920, ]
schinus_nocoords_1930 <- schinus_nocoords_years[schinus_nocoords_years$year<1930, ]
schinus_nocoords_1940 <- schinus_nocoords_years[schinus_nocoords_years$year<1940, ]
schinus_nocoords_1950 <- schinus_nocoords_years[schinus_nocoords_years$year<1950, ]
schinus_nocoords_1960 <- schinus_nocoords_years[schinus_nocoords_years$year<1960, ]
schinus_nocoords_1970 <- schinus_nocoords_years[schinus_nocoords_years$year<1970, ]
schinus_nocoords_1980 <- schinus_nocoords_years[schinus_nocoords_years$year<1980, ]
schinus_nocoords_1990 <- schinus_nocoords_years[schinus_nocoords_years$year<1990, ]
schinus_nocoords_2000 <- schinus_nocoords_years[schinus_nocoords_years$year<2000, ]
schinus_nocoords_2010 <- schinus_nocoords_years[schinus_nocoords_years$year<2010, ]
schinus_nocoords_2020 <- schinus_nocoords_years[schinus_nocoords_years$year<2020, ]

anacardeaceae_coords_1880 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<1880, ]
anacardeaceae_coords_1890 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<1890, ]
anacardeaceae_coords_1900 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<1900, ]
anacardeaceae_coords_1910 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<1910, ]
anacardeaceae_coords_1920 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<1920, ]
anacardeaceae_coords_1930 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<1930, ]
anacardeaceae_coords_1940 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<1940, ]
anacardeaceae_coords_1950 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<1950, ]
anacardeaceae_coords_1960 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<1960, ]
anacardeaceae_coords_1970 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<1970, ]
anacardeaceae_coords_1980 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<1980, ]
anacardeaceae_coords_1990 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<1990, ]
anacardeaceae_coords_2000 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<2000, ]
anacardeaceae_coords_2010 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<2010, ]
anacardeaceae_coords_2020 <- anacardeaceae_coords_years[anacardeaceae_coords_years$year<2020, ]

anacard_XY_1880<-SpatialPoints(cbind(anacardeaceae_coords_1880$decimalLongitude, anacardeaceae_coords_1880$decimalLatitude), getCRS('wgs84', TRUE))
anacard_XY_1890<-SpatialPoints(cbind(anacardeaceae_coords_1890$decimalLongitude, anacardeaceae_coords_1890$decimalLatitude), getCRS('wgs84', TRUE))
anacard_XY_1900<-SpatialPoints(cbind(anacardeaceae_coords_1900$decimalLongitude, anacardeaceae_coords_1900$decimalLatitude), getCRS('wgs84', TRUE))
anacard_XY_1910<-SpatialPoints(cbind(anacardeaceae_coords_1910$decimalLongitude, anacardeaceae_coords_1910$decimalLatitude), getCRS('wgs84', TRUE))
anacard_XY_1920<-SpatialPoints(cbind(anacardeaceae_coords_1920$decimalLongitude, anacardeaceae_coords_1920$decimalLatitude), getCRS('wgs84', TRUE))
anacard_XY_1930<-SpatialPoints(cbind(anacardeaceae_coords_1930$decimalLongitude, anacardeaceae_coords_1930$decimalLatitude), getCRS('wgs84', TRUE))
anacard_XY_1940<-SpatialPoints(cbind(anacardeaceae_coords_1940$decimalLongitude, anacardeaceae_coords_1940$decimalLatitude), getCRS('wgs84', TRUE))
anacard_XY_1950<-SpatialPoints(cbind(anacardeaceae_coords_1950$decimalLongitude, anacardeaceae_coords_1950$decimalLatitude), getCRS('wgs84', TRUE))
anacard_XY_1960<-SpatialPoints(cbind(anacardeaceae_coords_1960$decimalLongitude, anacardeaceae_coords_1960$decimalLatitude), getCRS('wgs84', TRUE))
anacard_XY_1970<-SpatialPoints(cbind(anacardeaceae_coords_1970$decimalLongitude, anacardeaceae_coords_1970$decimalLatitude), getCRS('wgs84', TRUE))
anacard_XY_1980<-SpatialPoints(cbind(anacardeaceae_coords_1980$decimalLongitude, anacardeaceae_coords_1980$decimalLatitude), getCRS('wgs84', TRUE))
anacard_XY_1990<-SpatialPoints(cbind(anacardeaceae_coords_1990$decimalLongitude, anacardeaceae_coords_1990$decimalLatitude), getCRS('wgs84', TRUE))
anacard_XY_2000<-SpatialPoints(cbind(anacardeaceae_coords_2000$decimalLongitude, anacardeaceae_coords_2000$decimalLatitude), getCRS('wgs84', TRUE))
anacard_XY_2010<-SpatialPoints(cbind(anacardeaceae_coords_2010$decimalLongitude, anacardeaceae_coords_2010$decimalLatitude), getCRS('wgs84', TRUE))
anacard_XY_2020<-SpatialPoints(cbind(anacardeaceae_coords_2020$decimalLongitude, anacardeaceae_coords_2020$decimalLatitude), getCRS('wgs84', TRUE))

anacard_XY_df <- list(anacard_XY_1880,anacard_XY_1890, anacard_XY_1900, anacard_XY_1910, anacard_XY_1920, 
                      anacard_XY_1930, anacard_XY_1940, anacard_XY_1950, anacard_XY_1960, anacard_XY_1970,
                      anacard_XY_1980, anacard_XY_1990, anacard_XY_2000, anacard_XY_2010, anacard_XY_2020)


anacardeaceae_nocoords_1880 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<1880, ] #no anacards
anacardeaceae_nocoords_1890 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<1890, ]
anacardeaceae_nocoords_1900 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<1900, ]
anacardeaceae_nocoords_1910 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<1910, ]
anacardeaceae_nocoords_1920 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<1920, ]
anacardeaceae_nocoords_1930 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<1930, ]
anacardeaceae_nocoords_1940 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<1940, ]
anacardeaceae_nocoords_1950 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<1950, ]
anacardeaceae_nocoords_1960 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<1960, ]
anacardeaceae_nocoords_1970 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<1970, ]
anacardeaceae_nocoords_1980 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<1980, ]
anacardeaceae_nocoords_1990 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<1990, ]
anacardeaceae_nocoords_2000 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<2000, ]
anacardeaceae_nocoords_2010 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<2010, ]
anacardeaceae_nocoords_2020 <- anacardeaceae_nocoords_years[anacardeaceae_nocoords_years$year<2020, ]


schinus_coords_df<-list(schinus_coords_1880, schinus_coords_1890, schinus_coords_1900, schinus_coords_1910, 
                        schinus_coords_1920, schinus_coords_1930, schinus_coords_1940, schinus_coords_1950,
                        schinus_coords_1960, schinus_coords_1970, schinus_coords_1980, schinus_coords_1990,
                        schinus_coords_2000, schinus_coords_2010, schinus_coords_2020)
schinus_nocoords_df <- list (schinus_nocoords_1880,schinus_nocoords_1890, schinus_nocoords_1900,
                             schinus_nocoords_1910, schinus_nocoords_1920, schinus_nocoords_1930,
                             schinus_nocoords_1940, schinus_nocoords_1950, schinus_nocoords_1960, 
                             schinus_nocoords_1970, schinus_nocoords_1980, schinus_nocoords_1990,
                             schinus_nocoords_2000, schinus_nocoords_2010, schinus_nocoords_2020)
anacard_coords_df <-list(anacardeaceae_coords_1880, anacardeaceae_coords_1890, anacardeaceae_coords_1900,
                         anacardeaceae_coords_1910, anacardeaceae_coords_1920, anacardeaceae_coords_1930,
                         anacardeaceae_coords_1940, anacardeaceae_coords_1950, anacardeaceae_coords_1960,
                         anacardeaceae_coords_1970, anacardeaceae_coords_1980, anacardeaceae_coords_1990,
                         anacardeaceae_coords_2000, anacardeaceae_coords_2010, anacardeaceae_coords_2020)
anacard_nocoords_df <-list(anacardeaceae_nocoords_1880,anacardeaceae_nocoords_1890, anacardeaceae_nocoords_1900,
                           anacardeaceae_nocoords_1910, anacardeaceae_nocoords_1920, anacardeaceae_nocoords_1930,
                           anacardeaceae_nocoords_1940, anacardeaceae_nocoords_1950, anacardeaceae_nocoords_1960,
                           anacardeaceae_nocoords_1970, anacardeaceae_nocoords_1980, anacardeaceae_nocoords_1990,
                           anacardeaceae_nocoords_2000, anacardeaceae_nocoords_2010, anacardeaceae_nocoords_2020)

#refer to a single df by schinus_coords_df[[i]]


counties <- rgdal::readOGR (
  'H:/Global Change Program/Research/Occupancy Modeling of Messy Detections/Kelley Erickson/Schinus spread/Data/GADM/ver2pt8/WGS84',
  'USA_adm2'
)
florida <- counties[counties@data$NAME_1=="Florida", ]
x <- data.frame(numSchinus=rep(0, nrow(florida)), numAnacard=rep(0, nrow(florida)))
florida@data <- omnibus::insertCol(x, into=florida@data, at='VARNAME_2', before=F)

florida_1880 <- florida
florida_1890 <- florida
florida_1900 <- florida
florida_1910 <- florida
florida_1920 <- florida
florida_1930 <- florida
florida_1940 <- florida
florida_1950 <- florida
florida_1960 <- florida
florida_1970 <- florida
florida_1980 <- florida
florida_1990 <- florida
florida_2000 <- florida
florida_2010 <- florida
florida_2020 <- florida 

florida_df <- list(florida_1880, florida_1890, florida_1900, florida_1910, florida_1920, florida_1930, 
                   florida_1940, florida_1950, florida_1960, florida_1970, florida_1980, florida_1990,
                   florida_2000, florida_2010, florida_2020)



for (i in 1:15) {
#i =1

AnacardInFlorida <- over(anacard_XY_df[[i]], florida_df[[i]]) 
AnacardInFlorida <- AnacardInFlorida[!is.na(AnacardInFlorida$NAME_2), ]

SchinusInFlorida <- over(schinus_XY_df[[i]], florida_df[[i]])
SchinusInFlorida <- SchinusInFlorida[!is.na(SchinusInFlorida$NAME_2), ]


# find number of records in each county (do by county)
#For all Tracheophyta
if(nrow(AnacardInFlorida)>0) {
  for (j in 1:nrow(AnacardInFlorida)) {
  
    thisRow <- which(florida_df[[i]]$NAME_1 == AnacardInFlorida$NAME_1[j] & florida_df[[i]]$NAME_2 == AnacardInFlorida$NAME_2[j])
    florida_df[[i]]@data$numAnacard[thisRow] <- florida_df[[i]]@data$numAnacard[thisRow] + 1
  
}
}

if(nrow(SchinusInFlorida)>0) {
  for (k in 1:nrow(SchinusInFlorida)) {
    thisRow <- which(florida_df[[i]]$NAME_1 == SchinusInFlorida$NAME_1[k] & florida_df[[i]]$NAME_2 == SchinusInFlorida$NAME_2[k])
    florida_df[[i]]@data$numSchinus[thisRow] <- florida_df[[i]]@data$numSchinus[thisRow] + 1
  }
  
}
} 




#Now to tackle all of the records without coordinates: 
#First focus on the Anacardiaceae 

for (i in 1:15) {
  if (nrow(anacard_nocoords_df[[i]]) > 0) {
  
    # ensure county names match GADM
    anacard_nocoords_df[[i]]$county <- gsub(anacard_nocoords_df[[i]]$county, pattern='(', replacement='', fixed=T)
    anacard_nocoords_df[[i]]$county <- gsub(anacard_nocoords_df[[i]]$county, pattern=')', replacement='', fixed=T)
    anacard_nocoords_df[[i]]$county <- gsub(anacard_nocoords_df[[i]]$county, pattern=' County', replacement='')
    anacard_nocoords_df[[i]]$county <- gsub(anacard_nocoords_df[[i]]$county, pattern='County of ', replacement='')
    anacard_nocoords_df[[i]]$county <- gsub(anacard_nocoords_df[[i]]$county, pattern=' Cty.', replacement='')
    anacard_nocoords_df[[i]]$county <- gsub(anacard_nocoords_df[[i]]$county, pattern=' Cty', replacement='')
    anacard_nocoords_df[[i]]$county <- gsub(anacard_nocoords_df[[i]]$county, pattern=' Co.', replacement='')
    anacard_nocoords_df[[i]]$county <- gsub(anacard_nocoords_df[[i]]$county, pattern='Ste ', replacement='Saint ')
    anacard_nocoords_df[[i]]$county <- gsub(anacard_nocoords_df[[i]]$county, pattern='St. ', replacement='Saint ')
  
  
  # add each record to county in which it falls--match state/province AND county
  
      for (j in 1:nrow(anacard_nocoords_df[[i]])) {
    
        thisRow <- which(as.character(florida_df[[i]]$NAME_1) == as.character(anacard_nocoords_df[[i]]$stateProvince[j]) & as.character(florida_df[[i]]$NAME_2) == as.character(anacard_nocoords_df[[i]]$county[j]))
        florida_df[[i]]@data$numAnacard[thisRow] <- florida_df[[i]]@data$numAnacard[thisRow] + 1
    
       }
  
  }
}

#Prior to running the above script (so only including occurences with coordinates, the total number of Anacards
# is 1024 ) Afterwards it is 4545 
sum(florida_df[[1]]$numAnacard,florida_df[[2]]$numAnacard, florida_df[[3]]$numAnacard,florida_df[[4]]$numAnacard,
    florida_df[[5]]$numAnacard, florida_df[[6]]$numAnacard, florida_df[[7]]$numAnacard, florida_df[[8]]$numAnacard,
    florida_df[[9]]$numAnacard, florida_df[[10]]$numAnacard, florida_df[[11]]$numAnacard, florida_df[[12]]$numAnacard,
    florida_df[[13]]$numAnacard, florida_df[[14]]$numAnacard, florida_df[[15]]$numAnacard)


#Now Schinus: 

for (i in 1:15) {
  if (nrow(schinus_nocoords_df[[i]]) > 0) {
    
    # ensure county names match GADM
    schinus_nocoords_df[[i]]$county <- gsub(schinus_nocoords_df[[i]]$county, pattern='(', replacement='', fixed=T)
    schinus_nocoords_df[[i]]$county <- gsub(schinus_nocoords_df[[i]]$county, pattern=')', replacement='', fixed=T)
    schinus_nocoords_df[[i]]$county <- gsub(schinus_nocoords_df[[i]]$county, pattern=' County', replacement='')
    schinus_nocoords_df[[i]]$county <- gsub(schinus_nocoords_df[[i]]$county, pattern='County of ', replacement='')
    schinus_nocoords_df[[i]]$county <- gsub(schinus_nocoords_df[[i]]$county, pattern=' Cty.', replacement='')
    schinus_nocoords_df[[i]]$county <- gsub(schinus_nocoords_df[[i]]$county, pattern=' Cty', replacement='')
    schinus_nocoords_df[[i]]$county <- gsub(schinus_nocoords_df[[i]]$county, pattern=' Co.', replacement='')
    schinus_nocoords_df[[i]]$county <- gsub(schinus_nocoords_df[[i]]$county, pattern='Ste ', replacement='Saint ')
    schinus_nocoords_df[[i]]$county <- gsub(schinus_nocoords_df[[i]]$county, pattern='St. ', replacement='Saint ')
    
    
    # add each record to county in which it falls--match state/province AND county
    
    for (j in 1:nrow(schinus_nocoords_df[[i]])) {
      
      thisRow <- which(as.character(florida_df[[i]]$NAME_1) == as.character(schinus_nocoords_df[[i]]$stateProvince[j]) & as.character(florida_df[[i]]$NAME_2) == as.character(schinus_nocoords_df[[i]]$county[j]))
      florida_df[[i]]@data$numSchinus[thisRow] <- florida_df[[i]]@data$numSchinus[thisRow] + 1
      
    }
    
  }
}


sum(florida_df[[1]]$numSchinus,florida_df[[2]]$numSchinus, florida_df[[3]]$numSchinus,florida_df[[4]]$numSchinus,
    florida_df[[5]]$numSchinus, florida_df[[6]]$numSchinus, florida_df[[7]]$numSchinus, florida_df[[8]]$numSchinus,
    florida_df[[9]]$numSchinus, florida_df[[10]]$numSchinus, florida_df[[11]]$numSchinus, florida_df[[12]]$numSchinus,
    florida_df[[13]]$numSchinus, florida_df[[14]]$numSchinus, florida_df[[15]]$numSchinus)




#### Use these objects in BayesLopod 

ld_Shape_1930 = shapeLopodData( Shapefile = florida_df[[6]],
                                   fieldN = "numAnacard",
                                   fieldY = "numSchinus",  
                                   Adjacency = T,
                                   keepFields = F)
ld_Shape_1940 = shapeLopodData( Shapefile = florida_df[[7]],
                                fieldN = "numAnacard",
                                fieldY = "numSchinus",  
                                Adjacency = T,
                                keepFields = F)
ld_Shape_1950 = shapeLopodData( Shapefile = florida_df[[8]],
                                fieldN = "numAnacard",
                                fieldY = "numSchinus",  
                                Adjacency = T,
                                keepFields = F)
ld_Shape_1960 = shapeLopodData( Shapefile = florida_df[[9]],
                                fieldN = "numAnacard",
                                fieldY = "numSchinus",  
                                Adjacency = T,
                                keepFields = F)
ld_Shape_1970 = shapeLopodData( Shapefile = florida_df[[10]],
                                fieldN = "numAnacard",
                                fieldY = "numSchinus",  
                                Adjacency = T,
                                keepFields = F)
ld_Shape_1980 = shapeLopodData( Shapefile = florida_df[[11]],
                                fieldN = "numAnacard",
                                fieldY = "numSchinus",  
                                Adjacency = T,
                                keepFields = F)
ld_Shape_1990 = shapeLopodData( Shapefile = florida_df[[12]],
                                fieldN = "numAnacard",
                                fieldY = "numSchinus",  
                                Adjacency = T,
                                keepFields = F)
ld_Shape_2000 = shapeLopodData( Shapefile = florida_df[[13]],
                                fieldN = "numAnacard",
                                fieldY = "numSchinus",  
                                Adjacency = T,
                                keepFields = F)
ld_Shape_2010 = shapeLopodData( Shapefile = florida_df[[14]],
                                fieldN = "numAnacard",
                                fieldY = "numSchinus",  
                                Adjacency = T,
                                keepFields = F)
ld_Shape_2020 = shapeLopodData( Shapefile = florida_df[[15]],
                                fieldN = "numAnacard",
                                fieldY = "numSchinus",  
                                Adjacency = T,
                                keepFields = F)
lsd_Shape_df <- list(ld_Shape_1930, ld_Shape_1940, ld_Shape_1950, ld_Shape_1960, ld_Shape_1970, ld_Shape_1980,
                     ld_Shape_1990, ld_Shape_2000, ld_Shape_2010, ld_Shape_2020)

for (i in 1:10) {
schinusSEff_1930 = spplot(florida_df[[6]], zcol = c("numAnacard"), main = "Sampling Effort-1930", border=NA)
schinusDetect_1930 = spplot(florida_df[[6]], zcol = c("numSchinus"), main = "Detections-1930", border=NA)
plot(schinusSEff_1930, split = c(1,1,2,1), xlab="1930", more = T, border=NA)
plot(schinusDetect_1930,  split = c(2,1,2,1), more = F, border=NA)
}
x11()
par(mfrow=c(1, 2))
spplot(florida_df[[6]], zcol=c("numAnacard"), main="Sampling Effort", border=NA, more=T)
spplot(florida_df[[6]], zcol=c("numSchinus"), main = "Detections", border=NA, split=c(2, 1, 2, 1), more=T)




rename <- function(x){
  if (x < 10) {
    return(name <- paste('000',i,'plot.png',sep=''))
  }
  if (x < 100 && i >= 10) {
    return(name <- paste('00',i,'plot.png', sep=''))
  }
  if (x >= 100) {
    return(name <- paste('0', i,'plot.png', sep=''))
  }
}

directory<-paste0('./', bquote(.(j)))

dir.create(directory, recursive=TRUE, showWarnings=FALSE)
setwd(directory)
for(i in 1:10){
  name <- rename(i)
  #saves the plot as a .png file in the working directory
  png(name)
  par(xpd=T)
  schinusSEff= spplot(florida_df[[i+5]], zcol=c("numAnacard"), main = bquote("Sampling Effort-" ~ .(categories[i+5])), border=NA)
  schinusDetect = spplot(florida_df[[i+5]], zcol=c("numSchinus"), main = bquote("Detections-" ~ .(categories[i+5])), border=NA)
  plot(schinusSEff, split=c(1, 1, 2, 1), more=T, border=NA)
  plot(schinusDetect, split=c(2, 1, 2, 1), more=F, border=NA)
  dev.off()
}

fileName<-"SamplingEffort.gif"
#Use the purr and magick packages to combine the png files into a GIF
# This part of the code comes from a post by Ryan Peek on how to make GIF 
# files in R: https://ryanpeek.github.io/2016-10-19-animated-gif_maps_in_R/ 
list.files(path = "./", pattern = "*.png", full.names = T) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=2) %>% # animates, can opt for number of loops
  image_write(fileName) # write to current dir



#### Run the model for each year
#There were only 2 schinus in 1930, which is probably not enough detections 
mLopodShape_1930 = modelLopod(LopodData = ld_Shape_1930,
                                 varP = T,
                                 q = NULL,
                                 pmin = 0,
                                 CAR = T,
                                 nChains = 3,
                                 warmup = 500,
                                 sampling = 1000,
                                 nCores =3)
lopodTrace( mLopodShape_1930, inc_warmup = T)

mLopodShape_1960 = modelLopod(LopodData = ld_Shape_1960,
                              varP = T, 
                              q = NULL, 
                              pmin = 0,
                              CAR = T, 
                              nChains = 3, 
                              warmup = 500, 
                              sampling = 3000, 
                              nCores =3)
lopodTrace(mLopodShape_1960, inc_warmup=T)


#Moving to 1970, where there are 68 detections (i = 10)

mLopodShape_1970 = modelLopod(LopodData = ld_Shape_1970,
                              varP = T, 
                              q = NULL, 
                              pmin = 0,
                              CAR = T, 
                              nChains = 4, 
                              warmup = 500, 
                              sampling = 1000, 
                              nCores =4)
lopodTrace(mLopodShape_1970, inc_warmup=T)

mLopodShape_1980 = modelLopod(LopodData = ld_Shape_1980,
                              varP = T, 
                              q = NULL, 
                              pmin = 0, 
                              CAR = T, 
                              nChains = 3, 
                              warmup = 500, 
                              sampling = 1000, 
                              nCores =3)
lopodTrace(mLopodShape_1980, inc_warmup=T)

mLopodShape_1990 = modelLopod(LopodData = ld_Shape_1990,
                              varP = T, 
                              q = NULL, 
                              pmin = 0, 
                              CAR = T, 
                              nChains = 3, 
                              warmup = 500, 
                              sampling = 1000, 
                              nCores =3)
lopodTrace(mLopodShape_1990, inc_warmup=T)

mLopodShape_2000 = modelLopod(LopodData = ld_Shape_2000,
                              varP = T, 
                              q = NULL, 
                              pmin = 0, 
                              CAR = T, 
                              nChains = 3, 
                              warmup = 500, 
                              sampling = 1000, 
                              nCores =3)
lopodTrace(mLopodShape_2000, inc_warmup=T)

mLopodShape_2010 = modelLopod(LopodData = ld_Shape_2010,
                              varP = T, 
                              q = NULL, 
                              pmin = 0, 
                              CAR = T, 
                              nChains = 3, 
                              warmup = 500, 
                              sampling = 1000, 
                              nCores =3)
lopodTrace(mLopodShape_2010, inc_warmup=T)

#Have the same issue that I ran into before...
mLopodShape_2020 = modelLopod(LopodData = ld_Shape_2020,
                              varP = T, 
                              q = NULL, 
                              pmin = 0, 
                              CAR = T, 
                              nChains = 3, 
                              warmup = 500, 
                              sampling = 1000, 
                              nCores =3)
lopodTrace(mLopodShape_1980, inc_warmup=T)

#####
lopodDens(mLopodShape_1970, c("q", "pmin", "pmax"))




#1970
png("1970.png")
SchinusShape_1970 = mLopodShape_1970@LopodData@geoDataObject
psiShape_1970 = lopodShape(mLopodShape_1970, "psi_i", extrapolate = T,  quant = 0.05)
SchinusShape_1970@data[,"psi_05"] = psiShape_1970@data[,"psi_i"]
psiShape_1970 = lopodShape(mLopodShape_1970, "psi_i", extrapolate = T,  quant = 0.5)
SchinusShape_1970@data[,"psi_50"] = psiShape_1970@data[,"psi_i"]
psiShape_1970 = lopodShape(mLopodShape_1970, "psi_i", extrapolate = T,  quant = 0.95)
SchinusShape_1970@data[,"psi_95"] = psiShape_1970@data[,"psi_i"]
spplot( SchinusShape_1970,
        zcol = c("psi_05", "psi_50", "psi_95"),
        names.attr = c("Psi (5% quantile)", "Psi (median)", "Psi (95% quantile)"), main = "Occupancy (Psi): 1970")
dev.off()

#1980
png("1980.png")
SchinusShape_1980 = mLopodShape_1980@LopodData@geoDataObject
psiShape_1980 = lopodShape(mLopodShape_1980, "psi_i", extrapolate = T,  quant = 0.05)
SchinusShape_1980@data[,"psi_05"] = psiShape_1980@data[,"psi_i"]
psiShape_1980 = lopodShape(mLopodShape_1980, "psi_i", extrapolate = T,  quant = 0.5)
SchinusShape_1980@data[,"psi_50"] = psiShape_1980@data[,"psi_i"]
psiShape_1980 = lopodShape(mLopodShape_1980, "psi_i", extrapolate = T,  quant = 0.95)
SchinusShape_1980@data[,"psi_95"] = psiShape_1980@data[,"psi_i"]
spplot( SchinusShape_1980,
        zcol = c("psi_05", "psi_50", "psi_95"),
        names.attr = c("Psi (5% quantile)", "Psi (median)", "Psi (95% quantile)"), main = "Occupancy (Psi): 1980")
dev.off()


#1990
png("1990.png")
SchinusShape_1990 = mLopodShape_1990@LopodData@geoDataObject
psiShape_1990 = lopodShape(mLopodShape_1990, "psi_i", extrapolate = T,  quant = 0.05)
SchinusShape_1990@data[,"psi_05"] = psiShape_1990@data[,"psi_i"]
psiShape_1990 = lopodShape(mLopodShape_1990, "psi_i", extrapolate = T,  quant = 0.5)
SchinusShape_1990@data[,"psi_50"] = psiShape_1990@data[,"psi_i"]
psiShape_1990 = lopodShape(mLopodShape_1990, "psi_i", extrapolate = T,  quant = 0.95)
SchinusShape_1990@data[,"psi_95"] = psiShape_1990@data[,"psi_i"]
spplot( SchinusShape_1990,
        zcol = c("psi_05", "psi_50", "psi_95"),
        names.attr = c("Psi (5% quantile)", "Psi (median)", "Psi (95% quantile)"), main = "Occupancy (Psi): 1990")
dev.off()

#2000
png("2000.png")
SchinusShape_2000 = mLopodShape_2000@LopodData@geoDataObject
psiShape_2000 = lopodShape(mLopodShape_2000, "psi_i", extrapolate = T,  quant = 0.05)
SchinusShape_2000@data[,"psi_05"] = psiShape_2000@data[,"psi_i"]
psiShape_2000 = lopodShape(mLopodShape_2000, "psi_i", extrapolate = T,  quant = 0.5)
SchinusShape_2000@data[,"psi_50"] = psiShape_2000@data[,"psi_i"]
psiShape_2000 = lopodShape(mLopodShape_2000, "psi_i", extrapolate = T,  quant = 0.95)
SchinusShape_2000@data[,"psi_95"] = psiShape_2000@data[,"psi_i"]
spplot( SchinusShape_2000,
        zcol = c("psi_05", "psi_50", "psi_95"),
        names.attr = c("Psi (5% quantile)", "Psi (median)", "Psi (95% quantile)"), main = "Occupancy (Psi): 2000")
dev.off()

png("2010.png")
SchinusShape_2010 = mLopodShape_2010@LopodData@geoDataObject
psiShape_2010 = lopodShape(mLopodShape_2010, "psi_i", extrapolate = T,  quant = 0.05)
SchinusShape_2010@data[,"psi_05"] = psiShape_2010@data[,"psi_i"]
psiShape_2010 = lopodShape(mLopodShape_2010, "psi_i", extrapolate = T,  quant = 0.5)
SchinusShape_2010@data[,"psi_50"] = psiShape_2010@data[,"psi_i"]
psiShape_2010 = lopodShape(mLopodShape_2010, "psi_i", extrapolate = T,  quant = 0.95)
SchinusShape_2010@data[,"psi_95"] = psiShape_2010@data[,"psi_i"]
spplot( SchinusShape_2010,
        zcol = c("psi_05", "psi_50", "psi_95"),
        names.attr = c("Psi (5% quantile)", "Psi (median)", "Psi (95% quantile)"), main = "Occupancy (Psi): 2010")
dev.off()

#
png("2020.png")
SchinusShape_2020 = mLopodShape_2020@LopodData@geoDataObject
psiShape_2020 = lopodShape(mLopodShape_2020, "psi_i", extrapolate = T,  quant = 0.05)
SchinusShape_2020@data[,"psi_05"] = psiShape_2020@data[,"psi_i"]
psiShape_2020 = lopodShape(mLopodShape_2020, "psi_i", extrapolate = T,  quant = 0.5)
SchinusShape_2020@data[,"psi_50"] = psiShape_2020@data[,"psi_i"]


lopodDens(mLopodShape_1970, c("q", "pmin", "pmax"))
psiShape_2020 = lopodShape(mLopodShape_2020, "psi_i", extrapolate = T,  quant = 0.95)
SchinusShape_2020@data[,"psi_95"] = psiShape_2020@data[,"psi_i"]
spplot( SchinusShape_2020,
        zcol = c("psi_05", "psi_50", "psi_95"),
        names.attr = c("Psi (5% quantile)", "Psi (median)", "Psi (95% quantile)"), main = "Occupancy (Psi): 2020")
dev.off()



fileName<-"Occupancy.gif"
#Use the purr and magick packages to combine the png files into a GIF
# This part of the code comes from a post by Ryan Peek on how to make GIF 
# files in R: https://ryanpeek.github.io/2016-10-19-animated-gif_maps_in_R/ 
list.files(path = "./", pattern = "*.png", full.names = T) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=2) %>% # animates, can opt for number of loops
  image_write(fileName) # write to current dir