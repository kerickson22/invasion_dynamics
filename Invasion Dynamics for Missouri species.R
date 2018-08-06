##### Using historical records to plot the spread of invasive species through time 

library(rgbif)
library(purrr)
library(magick)

#increase memory limit to 3 GB
memory.limit(memory.limit() *2^30)

#A selection of invasive species from the MDC website (https://nature.mdc.mo.gov/status/invasive)
# that are possibly dispersed by animals 
speciesList<-c("Schinus terebinthifolia", "Elaeagnus umbellata", "Lonicera maackii", "Pyrus calleryana", "Rhamnus cathartica",
                "Nandina domestica", "Lonicera japonica", "Rosa multiflora", "Elaeagnus angustifolia",
                "Euonymus fortunei", "Duchesnea indica")

dat2<-data.frame()
#speciesList<-c("Elaeagnus umbellata", "Lonicera maackii")
for (i in speciesList) {
  cat ("starting: ", i, "\n")
  dat <- occ_search(scientificName = i, return = "data", limit=20000)
  dat<- data.frame( 
    idNum = dat$gbifID,
    rawSpecies=dat$scientificName,
    species = i,
    longitude=dat$decimalLongitude,
    latitude=dat$decimalLatitude,
    coordUncer=dat$coordinateUncertaintyInMeters,
    recordType=dat$basisOfRecord,
    country=dat$countryCode,
    state=dat$stateProvince,
    county=dat$county,
    locality=dat$locality,
    date=dat$eventDate,
    year=dat$year,
    institution=dat$institutionCode,
    identifiedBy=dat$identifiedBy
  )
  dat2 <- rbind(dat, dat2)
  cat("finished with: ", i,"\n");
}

save(dat2, file='00 Original data downloaded from GBIF.rData')

#With just two species there are 9689 observations

### Remove observations that are missing coordinates: 
dat2_missingCoords <- dat2[is.na(dat2$longitude), ]
dat2 <- dat2[!is.na(dat2$longitude), ]
#Removing observations that are missing a longitude coordinate gets rid of 4959
# records (more than half of the records!)

save(dat2, file='02 Remove missing coordinates.rData')
save(dat2_missingCoords, file= '02b Records that are missing coordinates.rData')

#### Remove observations that are missing year
dat2_missingYear <- dat2[is.na(dat2$year), ]
dat2 <- dat2[!is.na(dat2$year), ]

save(dat2, file = '03 Remove datapoints missing year.rData')
save(dat2_missingYear, file='03b Records that have coords but are missing year.rData')

### Remove coordinates that are not in the US? 
dat2_missingCountry <- subset (dat2, (dat2$country == 'none' | dat2$country == 'ZZ') )
dat2 <- dat2[dat2$country == "US", ]



save(dat2, file='04 Remove data points not in US.rData')
save(dat2_missingCountry, file = '04b Records that have coords and year but are missing country code.rData')

### Draw histograms of the number of records by year for each species 

summary(dat2$year)
for (i in speciesList) {
    species_i <- subset(dat2, dat2$species == i)
    hist(species_i$year,
    xlim=c(1873, 2018),
    xlab='Collection year',
    ylab='Number of records',
    main=bquote("Species: " ~ .(i) )
)
}

### Make yearly plots for each species


#Load in shapefile of United States (edited in ArcGIS)
unitedStates <- rgdal::readOGR (
  './United States Shape',
  'United_states'
)

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


#par(xpd=T)
#sp::plot(unitedStates, col='gray80', xlim=c(-120, -70), ylim=c(26, 54))
#points(dat2$longitude, dat2$latitude, bg='red', pch=21)
#mtext("All species, All Years", side=3, cex=2)


for (j in speciesList) {
  cat ("starting: ", j, "\n")
  species_j <- subset(dat2, dat2$species == j)
  years <- unique(species_j$year)
  #rearrange so that years are in ascending order 
  years2 <- years[order(years)]
  
  directory<-paste0('./', bquote(.(j)))
  
  dir.create(directory, recursive=TRUE, showWarnings=FALSE)
  
  #number of frames or plots
  frames <- length(years2)
  
  
  setwd(directory)
  #Create the individual frames of the GIF
  for(i in 1:frames){
    name <- rename(i)
    
    #saves the plot as a .png file in the working directory
    png(name)
    par(xpd=T)
    dat3<-subset(species_j, (species_j$year <=years2[i]))
    sp::plot(unitedStates, col='gray80', xlim=c(-120, -70), ylim=c(26, 54))
    points(dat3$longitude, dat3$latitude,
           bg='red', pch=21)
    mtext(bquote("Species:" ~  .(j)), side=3, cex=2)
    mtext(bquote("\n Year:" ~ .(years2[i])), side=1, cex=2)
    #text(-100, 54, bquote("Year:" ~ .(years2[i])), cex=2)
    dev.off()
  }
  
  fileName<-paste0(bquote(.(j)), ".gif")
  #Use the purr and magick packages to combine the png files into a GIF
  # This part of the code comes from a post by Ryan Peek on how to make GIF 
  # files in R: https://ryanpeek.github.io/2016-10-19-animated-gif_maps_in_R/ 
  list.files(path = "./", pattern = "*.png", full.names = T) %>% 
    map(image_read) %>% # reads each path file
    image_join() %>% # joins image
    image_animate(fps=2) %>% # animates, can opt for number of loops
    image_write(fileName) # write to current dir
  setwd("C:/Users/kerickson/Documents/GitHub/invasion_dynamics")
  
}













