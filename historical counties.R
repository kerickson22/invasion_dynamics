library(USAboundaries)
library(dplyr)
library(sf)
library(leaflet)


#Have to run the md_votes section first to prompt for an install... 
# Check that USAboundariesData is available, since it is not on CRAN
avail <- requireNamespace("USAboundariesData", quietly = TRUE)



if (avail) {
  md_votes <- read.csv(system.file("extdata", "md-115.csv",
                                   package = "USAboundaries"),
                       colClasses = c("character", "integer", "character", "integer",
                                      "character", "integer", "integer", "integer",
                                      "integer", "integer", "integer",  "numeric",
                                      "numeric", "numeric", "numeric", "numeric"))
  md_districts <- us_congressional(states = "MD")
}


fl_districts <- us_congressional(states="FL")
fl_1840 <- us_counties("1840-03-12", states="FL")
fl_1930 <- us_counties("1930-03-12", states = "FL")
plot(st_geometry(fl_1930))

fl_2000 <- us_counties("2000-03-12", states="FL")