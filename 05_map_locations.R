
## devtools::install_github("hadley/ggplot2@v2.2.0") #the currentv2.2.1 does not work
library(dplyr)
library(ggplot2)
library(ggmap)
library(gridExtra)
library(gdata)

RACE_data      <- read.csv("~/Dropbox (Personal)/Dartmouth/ActiveStepRCT_data/Race Data.csv")
table(RACE_data$SiteName)

## Retrieve GEO coordinates for 5 major facilities:
myLocations <- c(
  "Cheshire Medical Center, Keene, NH", 
  "Concord Hospital, Concord, NH",
  "Dartmouth-Hitchcock Medical Center, Lebanon, NH",
  "Elliot Health System, Manchester, NH",
  "White River Junction Medical Center, White River Junction, VT"
)

mySpatialDF <- geocode(myLocations)
rownames(mySpatialDF) <- myLocations

myMap <- get_googlemap(
  center  = "Dartmouth-Hitchcock Medical Center, Lebanon, NH", 
  markers = mySpatialDF, 
  maptype = "roadmap",
  zoom    = 8
)
ggmap(myMap, extent="device")
