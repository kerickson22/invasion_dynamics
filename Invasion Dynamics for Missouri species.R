##### Using historical records to plot the spread of invasive species through time 

library(rgbif)
library(purrr)
library(magick)

#increase memory limit to 3 GB
memory.limit(memory.limit() *2^30)

#A selection of invasive species from the MDC website (https://nature.mdc.mo.gov/status/invasive)
# that are possibly dispersed by animals 
speciesList<-c("Elaeagnus umbellata", "Lonicera maackii", "Pyrus calleryana", "Rhomnus cathartica",
                "Nandina domestica", "Lonicera japonica", "Rusa multiflora", "Elaeagnus angustifolia",
                "Euonymus fortunei", "Duchesnea indica")

E_umbellata <- occ_search(scientificName = "Elaeagnus umbellata", return =  "data", limit=20000)
gbifmap(E_umbellata)
save(E_umbellata, file='00 Original data downloaded from GBIF.rData')
#There are 7037 observations of 160 variables. 


E_umbellata<- data.frame( 
  idNum = E_umbellata$gbifID,
  rawSpecies=E_umbellata$scientificName,
  longitude=E_umbellata$decimalLongitude,
  latitude=E_umbellata$decimalLatitude,
  coordUncer=E_umbellata$coordinateUncertaintyInMeters,
  recordType=E_umbellata$basisOfRecord,
  country=E_umbellata$countryCode,
  state=E_umbellata$stateProvince,
  county=E_umbellata$county,
  locality=E_umbellata$locality,
  date=E_umbellata$eventDate,
  year=E_umbellata$year,
  institution=E_umbellata$institutionCode,
  identifiedBy=E_umbellata$identifiedBy
)
#Now there are 7037 observations of 14 variables. 
save(E_umbellata,  file='01 Saved as dataframe.rData')


#Remove observations that are missing coordinates 
E_umbellata_missingCoords<-E_umbellata[is.na(E_umbellata$longitude), ]
#There are 3860 observations that are missing coordinates (over half!)
save(E_umbellata_missingCoords, file='02b Observations missing coordinates.rData')
summary(E_umbellata_missingCoords$country)
# Of these 3860 observations that are missing coordinates, 902 have no country 'none'
# 244 records are listed as US
summary(E_umbellata_missingCoords$state)
# 21 recrods from Massachusetts
# 21 records from South Carolina
# 15 records from New York
# 13 records from Michigan
# 13 records from Missouri
# 12 records from Florida
# 11 records from North Carolina
# 9 records from Hawaii 
# 6 records from Virginia 
# 5 records from Georgia 
# 5 records from Oregon
# 4 records from Connecticut
# 4 records from New Hampshire
# 4 records from Pennsylvania
# 3 records from Louisiana
# 3 records from Maine
# 3 records from Rhode Island
# # 3 records from Wisconsin
# and 95 records from (Other)

summary(E_umbellata_missingCoords$county)






E_umbellata<-E_umbellata[!is.na(E_umbellata$longitude), ]
save(E_umbellata, file='02 Remove missing coordinates.rData')

#Remove observations without a year
E_umbellata_coordsMissingYear<-E_umbellata[is.na(E_umbellata$year), ]
#There are 218 observations that have coordinates but are missing a year 
save(E_umbellata_coordsMissingYear, file = '03b Observations with coordinates missing year.rData')


E_umbellata<- E_umbellata[!is.na(E_umbellata$year), ]
save(E_umbellata, file='03 Remove data points that are missing year.rData')

#Remove observations that are not in the US
#Is this step removing valid points?
E_umbellata_missingCountry<-E_umbellata[is.na(E_umbellata$country), ]
#There are no observations that lack a country 


E_umbellata <- E_umbellata[E_umbellata$country=="US", ]
save(E_umbellata, file='04 Remove data points not in US.rData')

#Remove coordinates with longitude outside of US 
E_umbellata <- E_umbellata[E_umbellata$longitude < -40, ]
save(E_umbellata, file='05 Remove data points with longitude outside of US.rData')

hist(E_umbellata$year,
     xlab='Collection year',
     ylab='Number of records',
     main=''
)


years<-unique(E_umbellata$year)
##Rearrange in ascending order 
years2<-years[order(years)]



#Load in shapefile of United States (edited in ArcGIS)
unitedStates <- rgdal::readOGR (
  './United States Shape',
  'United_states'
)

par(xpd=T)
sp::plot(unitedStates, col='gray80', xlim=c(-120, -70), ylim=c(26, 54))
points(E_umbellata$longitude, E_umbellata$latitude, bg='red', pch=21)
mtext("E. umbellata", side=3, cex=2)


#number of frames or plots
frames <- length(years2)

# function for creating file name with leading zeros
# makes it easier to process them sequentially
# This function was written by Dave Tang
# https://davetang.org/muse/2015/02/12/animated-plots-using-r/ 
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

setwd("C:/Users/kerickson/Documents/GitHub/invasion_dynamics/E_umbellata_plots")

#Create the individual frames of the GIF
for(i in 1:frames){
  name <- rename(i)
  
  #saves the plot as a .png file in the working directory
  png(name)
  par(xpd=T)
  dat<-subset(E_umbellata, (E_umbellata$year <=years2[i]))
  sp::plot(unitedStates, col='gray80', xlim=c(-120, -70), ylim=c(26, 54))
  points(dat$longitude, dat$latitude,
         bg='red', pch=21)
  mtext(bquote("Species: E. umbellata \n Year" ~ .(years2[i])), side=3, cex=2)
  #text(-100, 54, bquote("Year:" ~ .(years2[i])), cex=2)
  dev.off()
}

#Use the purr and magick packages to combine the png files into a GIF
# This part of the code comes from a post by Ryan Peek on how to make GIF 
# files in R: https://ryanpeek.github.io/2016-10-19-animated-gif_maps_in_R/ 
list.files(path = "./", pattern = "*.png", full.names = T) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=2) %>% # animates, can opt for number of loops
  image_write("E_umbellata_invasion.gif") # write to current dir

