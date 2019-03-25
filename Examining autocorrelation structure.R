library(stringdist)
library(spdep)
#Examine autocorrelation in data 
adjacency<-read.table("./Data/county_adjacency.txt", sep = "\t")
# Column 1: County Name: 2010 State and County Name
# Column 2: County GEOID: 2010 State and County FIPS Codes
# Column 3: Neighbor Name: 2010 State and County name of each neighboring county or county equivalent
# Column 4: Neighbor GEOID: 2010 State nad County FIPS codes of each neighboring county or county equivalent


colnames(adjacency) <- c("County_Name", "County_GEOID", "Neighbor_Name", "Neighbor_GEOID")


###Different approach using spdep package
flnbs <- poly2nb(florida)
summary(flnbs)

florida@polygons[[1]]@Polygons[[1]]@coords

florida@polygons[[3]]@Polygons[[1]]@coords
florida@polygons[[3]]@area
florida@polygons[[3]]@ID
florida@polygons[[3]]@labpt

florida@polygons[[66]]@ID


florida@polygons@ID
