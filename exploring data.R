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
library(classInt)
library(gridExtra)

setwd("C:/Users/kerickson/Documents/GitHub/invasion_dynamics")
#load(".RData")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
memory.limit(memory.limit() *2^30)


#All US occurrences of anacardiaceae 
anacardiaceae <- read.csv("./US Anacardiaceae records from GBIF/anacardiaceae_occurrences.csv", header=T)
#There are 43, 947 records of Anacardeacea in the United States 
#Global occurrences of schinus
schinus <- read.csv("./Schinus records from GBIF/schinus_records.csv", header=T)

tracheophyta<- read.csv("./Tracheophyta/occurrence2.csv", header=T)
#There are 383,697 Florida records of Tracheophytes 

#For those records without coordinates, want to only keep those
# that are in Florida, and have a county code 

#For those records that /do/ have coordinates, the metadata
# on countrycode, stateProvince, etc. doesn't seem to have
# always been filled out, so subsetting using those
# columns drastically reduces the amount of data available. 

#Instead, we will clip coordinates to only select those
# that fall within the boundaries of the florida shapefile
# (from gadm2)

#actually what I need is a state polygon

load("florida.RData")


states <- rgdal::readOGR (
  'C:/Users/kerickson/Documents/GitHub/invasion_dynamics/United States Shape',
  'USA_adm1'
)

florida_state <- states[states@data$NAME_1=="Florida", ]

schinus_coord<-schinus[!is.na(schinus$decimalLatitude), ]
schinus_coord <- schinus_coord[!is.na(schinus_coord$decimalLongitude), ]
#schinus_XY<-SpatialPoints(cbind(schinus_coord$decimalLongitude, schinus_coord$decimalLatitude), getCRS('wgs84', TRUE))
schinus_XY_df<-schinus_coord
coordinates(schinus_XY_df) <- c(134, 133)
proj4string(schinus_XY_df) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
proj4string(florida_state) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

anacard_coord<-anacardiaceae[!is.na(anacardiaceae$decimalLatitude), ]
anacard_coord<-anacard_coord[!is.na(anacard_coord$decimalLongitude), ]
anacard_XY_df<-anacard_coord
coordinates(anacard_XY_df) <- c(134, 133)
proj4string(anacard_XY_df) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

tracheophyte_coord<-tracheophyta[!is.na(tracheophyta$decimalLatitude),]
tracheophyte_coord<-tracheophyte_coord[!is.na(tracheophyte_coord$decimalLongitude), ]
tracheophyte_XY_df<-tracheophyte_coord
coordinates(tracheophyte_XY_df) <- c(134, 133)
proj4string(tracheophyte_XY_df) <- CRS("+proj=longlat +datum=WGS84 +no_defs + ellps=WGS84 +towgs84=0,0,0")
# thing<- schinus_XY_df[florida_state,]
# dim(thing@data)
# 
# summary(thing@data$date) #all NA
# summary(thing@data$recordedBy) #4003 levels, 99 in Florida S. Koptur
# summary(thing@data$startDayOfYear) #196NAs (out of 200)
# summary(thing@data$endDayOfYear) #all NA
# summary(thing@data$year) #1936-2018
# summary(thing@data$month) #most often collected in December and January, but also collected throughout
# summary(thing@data$day)
# summary(thing@data$eventDate)
# summary(thing@data$verbatimEventDate)
# summary(thing@data$coordinatePrecision) #NA
# unique(thing@data$georeferencedBy) # ""
# unique(thing@data$georeferenceProtocol) #""
# unique(thing@data$identifiedBy) #8 people + ""
# summary(thing@data$dateIdentified) #140 dates of the form 2014-12-18T08:30Z
# unique(thing@data$issue) #10 categories of issue 
# # [1] TAXON_MATCH_FUZZY;GEODETIC_DATUM_ASSUMED_WGS84:24                        
# # [2] TAXON_MATCH_FUZZY;GEODETIC_DATUM_ASSUMED_WGS84;COORDINATE_ROUNDED : 139   
# # [3] "" : 2                                                         
# # [4] GEODETIC_DATUM_ASSUMED_WGS84 :27                                           
# # [5] GEODETIC_DATUM_ASSUMED_WGS84;PRESUMED_NEGATED_LONGITUDE: 1                
# # [6] GEODETIC_DATUM_ASSUMED_WGS84;GEODETIC_DATUM_INVALID: 2                 
# # [7] GEODETIC_DATUM_ASSUMED_WGS84;GEODETIC_DATUM_INVALID;COORDINATE_ROUNDED: 1
# # [8] GEODETIC_DATUM_ASSUMED_WGS84;IDENTIFIED_DATE_UNLIKELY  :2                
# # [9] TAXON_MATCH_FUZZY;GEODETIC_DATUM_ASSUMED_WGS84;IDENTIFIED_DATE_UNLIKELY: 1
# # [10] RECORDED_DATE_MISMATCH;TAXON_MATCH_FUZZY;GEODETIC_DATUM_ASSUMED_WGS84: 1 
# summary(thing@data$hasGeospatialIssues) #allfalse
# 
# 
# #Assign counties to each of the records 
# thing@data$county<-factor(thing@data$county) #remove extra levels 
# 
# column<-rep(0, 200)
# colnames(column) <- "test2"
# test<-cbind(thing@data, b=rep(0,200))
# 
# column<-data.frame(column)
# test<-thing
# 
# #test2<-omnibus::insertCol(column, into=test@data, at='repatriated', before=FALSE)
# 
# test[["z3"]]=rep(0,200)

#Select only those records that fall within the state polygon 
schinus_coords <- schinus_XY_df[florida_state,]
anacard_coords <- anacard_XY_df[florida_state,]
tracheophyte_coords <- tracheophyte_XY_df[florida_state, ]

#Assign counties to each of the records 

#First, remove extra levels of 'county':
schinus_coords@data$county <- factor(schinus_coords@data$county)
anacard_coords@data$county <- factor(anacard_coords@data$county)
tracheophyte_coords@data$county <- factor(tracheophyte_coords@data$county)

#create a column called geo_county for storing the county of each record
# as determined by overlaying county level polygons 

schinus_coords[["geo_county"]]<-rep("NA", 200)
anacard_coords[["geo_county"]]<-rep("NA", 623)
tracheophyte_coords[["geo_county"]]<-rep("NA", 65689)

#test<-schinus_coords

counties<-unique(florida@data$NAME_2) #There are 67 counties




schinus_coords@data$geo_county <- as.character(over(schinus_coords, florida)$NAME_2)
anacard_coords@data$geo_county <- as.character(over(anacard_coords, florida)$NAME_2)
tracheophyte_coords@data$geo_county <- as.character(over(tracheophyte_coords, florida)$NAME_2)

#Now deal with records without coordinates: 

schinus_nocoords <- schinus[is.na(schinus$decimalLongitude)| is.na(schinus$decimalLatitude), ]
anacard_nocoords <- anacardiaceae[is.na(anacardiaceae$decimalLongitude)|is.na(anacardiaceae$decimalLatitude), ]
tracheophyte_nocoords <- tracheophyta[is.na(tracheophyta$decimalLongitude)|is.na(tracheophyta$decimalLatitude), ]

#Select only those records within Florida
schinus_nocoords$stateProvince <- gsub(schinus_nocoords$stateProvince, pattern='Miami, Florida', replacement='Florida')
schinus_nocoords$stateProvince <- gsub(schinus_nocoords$stateProvince, pattern='Flórida', replacement='Florida')
schinus_nocoords <- schinus_nocoords[schinus_nocoords$stateProvince =="Florida", ]

#Select only those records within Florida
anacard_nocoords$stateProvince <- gsub(anacard_nocoords$stateProvince, pattern='Miami, Florida', replacement='Florida')
anacard_nocoords$stateProvince <- gsub(anacard_nocoords$stateProvince, pattern='Flórida', replacement='Florida')
anacard_nocoords <- anacard_nocoords[anacard_nocoords$stateProvince =="Florida", ]

#Select only those records within Florida
tracheophyte_nocoords$stateProvince <- gsub(tracheophyte_nocoords$stateProvince, pattern='Miami, Florida', replacement='Florida')
tracheophyte_nocoords$stateProvince <- gsub(tracheophyte_nocoords$stateProvince, pattern='Flórida', replacement='Florida')
tracheophyte_nocoords <- tracheophyte_nocoords[tracheophyte_nocoords$stateProvince =="Florida", ]




#Then clean up counties:
schinus_nocoords$county <- gsub(schinus_nocoords$county, pattern='(', replacement='', fixed=T)
schinus_nocoords$county <- gsub(schinus_nocoords$county, pattern=')', replacement='', fixed=T)
schinus_nocoords$county <- gsub(schinus_nocoords$county, pattern=' County', replacement='')
schinus_nocoords$county <- gsub(schinus_nocoords$county, pattern='County of ', replacement='')
schinus_nocoords$county <- gsub(schinus_nocoords$county, pattern=' Cty.', replacement='')
schinus_nocoords$county <- gsub(schinus_nocoords$county, pattern=' Cty', replacement='')
schinus_nocoords$county <- gsub(schinus_nocoords$county, pattern=' Co.', replacement='')
schinus_nocoords$county <- gsub(schinus_nocoords$county, pattern='Ste ', replacement='Saint ')
schinus_nocoords$county <- gsub(schinus_nocoords$county, pattern='St. ', replacement='Saint ')
schinus_nocoords$county <- gsub(schinus_nocoords$county, pattern='Miami-Dade', replacement='Miami')
schinus_nocoords$county <- gsub(schinus_nocoords$county, pattern='Dade', replacement='Miami-Dade')
schinus_nocoords$county <- gsub(schinus_nocoords$county, pattern='Miami', replacement='Miami-Dade')
schinus_nocoords$county <- gsub(schinus_nocoords$county, pattern='Miami-Dade-Dade', replacement='Miami-Dade')

schinus_nocoords <- subset(schinus_nocoords, schinus_nocoords$county !="")
#There are 268 records without coordinates 

#Then clean up counties:
anacard_nocoords$county <- gsub(anacard_nocoords$county, pattern='(', replacement='', fixed=T)
anacard_nocoords$county <- gsub(anacard_nocoords$county, pattern=')', replacement='', fixed=T)
anacard_nocoords$county <- gsub(anacard_nocoords$county, pattern=' County', replacement='')
anacard_nocoords$county <- gsub(anacard_nocoords$county, pattern='County of ', replacement='')
anacard_nocoords$county <- gsub(anacard_nocoords$county, pattern=' Cty.', replacement='')
anacard_nocoords$county <- gsub(anacard_nocoords$county, pattern=' Cty', replacement='')
anacard_nocoords$county <- gsub(anacard_nocoords$county, pattern=' Co.', replacement='')
anacard_nocoords$county <- gsub(anacard_nocoords$county, pattern='Ste ', replacement='Saint ')
anacard_nocoords$county <- gsub(anacard_nocoords$county, pattern='St. ', replacement='Saint ')
anacard_nocoords$county <- gsub(anacard_nocoords$county, pattern='Miami-Dade', replacement='Miami')
anacard_nocoords$county <- gsub(anacard_nocoords$county, pattern='Dade', replacement='Miami-Dade')
anacard_nocoords$county <- gsub(anacard_nocoords$county, pattern='Miami', replacement='Miami-Dade')
anacard_nocoords$county <- gsub(anacard_nocoords$county, pattern='Miami-Dade-Dade', replacement='Miami-Dade')

anacard_nocoords <- subset(anacard_nocoords, anacard_nocoords$county !="")


#Then clean up counties:
tracheophyte_nocoords$county <- gsub(tracheophyte_nocoords$county, pattern='(', replacement='', fixed=T)
tracheophyte_nocoords$county <- gsub(tracheophyte_nocoords$county, pattern=')', replacement='', fixed=T)
tracheophyte_nocoords$county <- gsub(tracheophyte_nocoords$county, pattern=' County', replacement='')
tracheophyte_nocoords$county <- gsub(tracheophyte_nocoords$county, pattern='County of ', replacement='')
tracheophyte_nocoords$county <- gsub(tracheophyte_nocoords$county, pattern=' Cty.', replacement='')
tracheophyte_nocoords$county <- gsub(tracheophyte_nocoords$county, pattern=' Cty', replacement='')
tracheophyte_nocoords$county <- gsub(tracheophyte_nocoords$county, pattern=' Co.', replacement='')
tracheophyte_nocoords$county <- gsub(tracheophyte_nocoords$county, pattern='Ste ', replacement='Saint ')
tracheophyte_nocoords$county <- gsub(tracheophyte_nocoords$county, pattern='St. ', replacement='Saint ')
tracheophyte_nocoords$county <- gsub(tracheophyte_nocoords$county, pattern='Miami-Dade', replacement='Miami')
tracheophyte_nocoords$county <- gsub(tracheophyte_nocoords$county, pattern='Dade', replacement='Miami-Dade')
tracheophyte_nocoords$county <- gsub(tracheophyte_nocoords$county, pattern='Miami', replacement='Miami-Dade')
tracheophyte_nocoords$county <- gsub(tracheophyte_nocoords$county, pattern='Miami-Dade-Dade', replacement='Miami-Dade')

tracheophyte_nocoords <- subset(tracheophyte_nocoords, tracheophyte_nocoords$county !="")








schinus_coords@data$geo_county <- as.character(over(schinus_coords, florida)$NAME_2)






#unique(schinus_nocoords$recordedBy, schinus_coords$recordedBy)
#There are 188 unique collectors .... not quite because J. Richard Abbott is there a lot
  
collectors <- c(as.character(schinus_nocoords$recordedBy), as.character(schinus_coords$recordedBy),
                as.character(anacard_nocoords$recordedBy), as.character(schinus_coords$recordedBy))
collectors_tracheophytes <- c(as.character(tracheophyte_coords$recordedBy), as.character(tracheophyte_nocoords$recordedBy))
collectors_tracheophytes <- unique(collectors_tracheophytes)
collectors_tracheophytes<-data.frame(collectors_tracheophytes)
collectors_list_tracheophytes <- c(as.character(tracheophyte_coords$recordedBy), as.character(tracheophyte_nocoords$recordedBy))
collectors_list_tracheophytes <- unique(collectors_list_tracheophytes)

collectors <-unique(collectors)
collectors2<-data.frame(collectors)
collectors2[["numSchinusCoords"]]<-rep(NA, length(collectors))
collectors2[["numSchinusNoCoords"]] <-rep(NA, length(collectors))
collectors2[["numAnacardCoords"]]<-rep(NA, length(collectors))
collectors2[["numAnacardNoCoords"]]<-rep(NA, length(collectors))
for (i in 1:length(collectors)){
  collectors2$numSchinusCoords[i]<-dim(schinus_coords[schinus_coords@data$recordedBy==collectors[i], ] )[1]
  collectors2$numSchinusNoCoords[i]<-dim(schinus_nocoords[schinus_nocoords$recordedBy==collectors[i], ])[1]
  collectors2$numAnacardCoords[i]<-dim(anacard_coords[anacard_coords@data$recordedBy==collectors[i], ])[1]
  collectors2$numAnacardNoCoords[i]<-dim(anacard_nocoords[anacard_nocoords$recordedBy==collectors[i], ])[1]
}
collectors2[["totalCollections"]] <-rep(NA, length(collectors))
collectors2[["totalSchinus"]] <-rep(NA, length(collectors))
#Only sum up Anacardiaceae because they include the Schinus records (I hope)
for(i in 1:length(collectors)) {
  collectors2$totalCollections[i] <- sum(collectors2[i, 4:5 ])
  collectors2$totalSchinus[i] <-sum(collectors2[i, 2:3])
}

collectors2[, c(1, 6)]

hist(collectors2$totalCollections)

plot(collectors2$totalCollections)
points(collectors2$totalSchinus, col="blue")


ordered_collectors <- collectors2[order(-collectors2$totalCollections, -collectors2$totalSchinus),]

x11()
par(las=2)
par(mar=c(18, 4.1, 4.1, 2.1))
barplot(log(ordered_collectors$totalCollections[1:205]), names.arg=ordered_collectors$collectors[1:205], cex.names=.6)
barplot(log(ordered_collectors$totalSchinus[1:205]), add=T, col="blue")



#Plot everyone 
x11()
par(las=2)
par(mfrow=c(2,1))
par(mar=c(12, 4.1, 4.1, 2.1))
barplot(log(ordered_collectors$totalCollections[1:205])+1, ylim=c(0,4), names.arg=ordered_collectors$collectors[1:205], cex.names=.6)
barplot(log(ordered_collectors$totalSchinus[1:205])+1, add=T, col="blue")

barplot(log(ordered_collectors$totalCollections[206:517])+1, ylim=c(0, 4), names.arg=ordered_collectors$collectors[206:517], cex.names=0.6)
barplot(log(ordered_collectors$totalSchinus[206:517])+1, add=T, col="blue")




x11()
par(mfrow=c(1,1))
par(las=2)
par(mar=c(11, 4.1, 4.1, 2.1))
barplot(ordered_collectors$totalCollections[1:22], names.arg=ordered_collectors$collectors[1:22])
barplot(ordered_collectors$totalSchinus[1:22], names.arg=ordered_collectors$collectors[1:22], add=T, col="blue")



length(ordered_collectors)
issues <- array(NA, 1)
collectors <- subset(collectors, collectors!="")
for (i in 1:548){
  thing<-agrep(collectors[i], collectors)
        if (length(thing)>1) {
          issues <- c(issues, collectors[i])
        }
}

#Now do by hand

agrep(issues[49], collectors, value=T)
#Keep going with 41 next time

#collectors <- c(as.character(schinus_nocoords$recordedBy), as.character(schinus_coords$recordedBy),
#as.character(anacard_nocoords$recordedBy), as.character(schinus_coords$recordedBy))
#Actually, think I have to do this for each source of data 

rename_collectors<-function(dataframe_varname) {
  dataframe_varname <- gsub(dataframe_varname, pattern='Wm. L. McCart', replacement='W. L. McCart', fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="George R. Cooley", replacement="G. R. Cooley", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="D. P. Cole", replacement="D. Cole", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="J. Dan Skean, Jr.", replacement="J. D. Skean, Jr.", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="B. F. Hansen, D.R. Richardson", replacement="B. F. Hansen, D. Richardson", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="R. P. Wunderlin, A.N. Arcuri, B.F. Hansen", replacement="R. P. Wunderlin, A. Arcuri, B. Hansen", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="D.B. Ward & E.S. Ford", replacement="D. B. Ward, E.S. Ford", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="D. B. Ward", replacement="D.B. Ward", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="O. K. Lakela", replacement="O. Lakela", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="L.M. Baltzell", replacement="L. Baltzell", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="James B. McFarlin", replacement="J.B. McFarlin", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="S. Mori; C. Gracie", replacement="S. Mori & C. Gracie", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="W. C. Brumbach", replacement="W.C. Brumbach", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="S. A. Mori, C. Gracie", replacement="S. Mori & C. Gracie", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="S. Mori; C. Gracie", replacement="S. Mori & C. Gracie", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="W.M. Buswell", replacement="W. M. Buswell", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="William T. Gillis", replacement="W. T. Gillis", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="B. F. Hansen, D.R. Richardson", replacement="B. F. Hansen, D. Richardson", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="Darin S. Penneys", replacement="D. S. Penneys", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="Jame L. Amoroso", replacement="J. L. Amoroso", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="S. A. Mori, C. Gracie", replacement="S. Mori & C. Gracie", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="S. Mori; C. Gracie", replacement="S. Mori & C. Gracie", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="P. K. Fawcett", replacement="P. Fawcett", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="R K Godfrey", replacement="R. K. Godfrey", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="R.K. Godfrey", replacement="R. K. Godfrey", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="F.C. Craighead", replacement="F. C. Craighead", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="Loran C. Anderson", replacement="L. C. Anderson", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="B. F. Hansen, G.L. Robinson", replacement="B. F. Hansen, G. Robinson", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="R. P. Wunderlin, A.N. Arcuri, B.F. Hansen", replacement="R. P. Wunderlin, A. Arcuri, B. Hansen", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="Robert R. Smith", replacement="R.R. Smith", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="Samuel F. Brockington", replacement="S.F. Brockington", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="W.L. Stern", replacement="W. L. Stern", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="D. B. Ward, E.S. Ford", replacement="D.B. Ward & E.S. Ford", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="Lucia McCulloch", replacement="L. McCulloch", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="A. S. Hitchcock", replacement="A.S. Hitchcock", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="Brenda J. Herring;", replacement="Brenda J. Herring", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="P. Genelle, G. Fleming", replacement="G. Fleming, P. Genelle", fixed=T)
  dataframe_varname <- gsub(dataframe_varname, pattern="C. V. Shroads", replacement="C. Shroads", fixed=T)
  dataframe_varname <- dataframe_varname
return(dataframe_varname)
}

#collectors <- c(as.character(schinus_nocoords$recordedBy), as.character(schinus_coords$recordedBy),
#as.character(anacard_nocoords$recordedBy), as.character(schinus_coords$recordedBy))
#Actually, think I have to do this for each source of data 
schinus_nocoords$recordedBy<-rename_collectors(schinus_nocoords$recordedBy)
schinus_coords$recordedBy<-rename_collectors(schinus_coords$recordedBy)
anacard_nocoords$recordedBy<-rename_collectors(anacard_nocoords$recordedBy)
anacard_coords$recordedBy<-rename_collectors(anacard_coords$recordedBy)


#Now let's look at the more active collectors: 
summary(ordered_collectors$totalCollections)
#3rd quantile = 2

abundant_collectors <- subset(ordered_collectors, ordered_collectors$totalCollections > 2)

temp<-subset(anacard_nocoords, anacard_nocoords$recordedBy==abundant_collectors$collectors[1])
temp2<-subset(anacard_coords, anacard_coords$recordedBy==abundant_collectors$collectors[1])

plot(florida)
points(temp2, col="darkgreen", pch=19)

plot(temp$year)

########Compare to Tracheophyte collectors
#First fix known issues: 

tracheophyte_nocoords$recordedBy<-rename_collectors(tracheophyte_nocoords$recordedBy)
tracheophyte_coords$recordedBy<-rename_collectors(tracheophyte_coords$recordedBy)

issues2 <- array(NA, 1)
collectors_tracheophytes <- subset(collectors_tracheophytes, collectors_tracheophytes!="")


for (i in 1:15126){
  thing<-agrep(collectors_tracheophytes$collectors_tracheophytes[i], collectors)
  if (length(thing)>1) {
    issues2 <- c(issues2, collectors_tracheophytes$collectors_tracheophytes[i])
  }
}

agrep(issues2[3], collectors_tracheophytes$collectors_tracheophytes, value=T)

#So there are still some records with known Anacardiaceae collectors (R. K. Godfrey) with messy "RecordedBy" fields. 
#Potentially 440 of them....

collected_by_known_collectors <- tracheophyte_coords$recordedBy %in% collectors
start_with_M <- nm %in% grep("^M", nm, value = TRUE)
subset(state.x77, start_with_M, Illiteracy:Murder)

tracheophyte_coords_knownCollectors<-subset(tracheophyte_coords, tracheophyte_coords$recordedBy %in% collectors)
tracheophyte_nocoords_knownCollectors<-subset(tracheophyte_nocoords, tracheophyte_nocoords$recordedBy %in% collectors)

nonAnacard_coords <-subset(tracheophyte_coords_knownCollectors, tracheophyte_coords_knownCollectors$family!="Anacardiaceae")
points(nonAnacard_coords, col="green")
points(schinus_coords, col="black", pch=19)


plot(florida, col="lightgrey")

points(tracheophyte_coords_knownCollectors, col="black")
#Red = tracheophyte records of known collectors 
points(anacard_coords, col="green", pch=19)
#blue= Anacardiaceae  collected by known Anacardiaceae collectors
#green = non-Anacardiaceae collected by known Anacardiaceae collectors
#points(nonAnacard_coords, col="green")
#black = records of Schinus terebinthifolia 
points(schinus_coords, col="red", pch=19)
legend("bottomleft", pch=c(1, 19, 19), col=c("black", "green", "red"), c("Tracheophytes", "Anacardiaceae", "Schinus"))




#Make plots of individual collector's efforts 

# find number of records in each county (do by county)
#For all Tracheophyta
if(nrow(AnacardInFlorida)>0) {
  for (j in 1:nrow(AnacardInFlorida)) {
    
    thisRow <- which(florida_df[[i]]$NAME_1 == AnacardInFlorida$NAME_1[j] & florida_df[[i]]$NAME_2 == AnacardInFlorida$NAME_2[j])
    florida_df[[i]]@data$numAnacard[thisRow] <- florida_df[[i]]@data$numAnacard[thisRow] + 1
    
  }
}
