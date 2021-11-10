# 03_Occ_filtering

# Joe Bellis & Sarah Dalrymple
# June/July 2021

# Thin occurrence records to correct for spatial bias #

# as before, use the following code to tidy up 'environment' if necessary
rm(list=ls())

# set wd if not already done so (and at least check you're in the right place with getwd)
setwd(choose.dir())
getwd()

# paste the full file pathway from the console into the setwd() function for future reference
setwd("C:/Users/Public/Documents/Modelling/...")

# DATA PREPARATION
ProjW = "+proj=longlat +ellps=WGS84 +no_defs" #map projection (EPSG: 4326)

library(sp)

# Load cleaned occurrence dataset and convert to spatial points df
species_occ <- read.table("occurrences.txt", header=TRUE)
xy <- species_occ [,2:3]
df <- species_occ 
Sp_occ_SPDF <- SpatialPointsDataFrame(coords=xy, data=df, proj4string = CRS(ProjW)) # convert to SPDF


library(rangeBuilder)

# Thin occurrences - note that 20 = 20km and may need to be reduced for rarer species (with prevalence of <0.1)
# be prepared for this to take a long time!!!
thinned_occ <- filterByProximity(xy, dist=20, mapUnits = FALSE, returnIndex = FALSE)
str(thinned_occ)

# Plot on SDM raster to check how many records will be used in modelling
# Load raster
PC1 <- raster("PCA_Out_PC1.tif")
plot(PC1)
# Plot points
Points_on_rast <- extract(PC1, thinned_occ) #?!# changed to 'thinned_occ' from 'Occ'
Removed_NAs <- as.data.frame(na.omit(Points_on_rast)) # this is the total number that will be used in the SDM

# IF NUMBER OF REMAINING PRESENCES IS LESS THAN 100, RE-THIN AT 10KM RESOLUTION, AS PER METHODOLOGY
# thinned_10km_occ <- filterByProximity(xy, dist=10, mapUnits = FALSE, returnIndex = FALSE)

# IF NUMBER OF REMAINING PRESENCES IS STILL LESS THAN 100, RE-THIN AT 5KM RESOLUTION, AS PER METHODOLOGY
# thinned_5km_occ <- filterByProximity(xy, dist=5, mapUnits = FALSE, returnIndex = FALSE)

# view points on map
library(maps)
world <- map("world", fill = FALSE)
points(thinned_occ, pch = 16, col = "blue")

write.csv(thinned_occ, file=paste("presences",".csv",sep=""))

# Create index showing which records were removed
# may also take a long time!!
removed_points <- filterByProximity(xy, dist=20, mapUnits = FALSE, returnIndex = TRUE)
# Add row number column called ID
xy$ID <- seq.int(nrow(xy))


# Get Long Lat of removed occurrences (these will be used for 'independent' evaluation in continuous boyce index)
library(dplyr)
Removed_occs <- filter(xy, ID %in% removed_points)

write.csv(Removed_occs, file=paste("removed_occs", ".csv",sep=""))

