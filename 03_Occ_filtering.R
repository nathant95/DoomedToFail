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
# needs to end with species name separated by underscore
setwd("C:/Users/Public/Documents/Modelling/...")

wd <-paste0(getwd())
species <- basename(wd)

# DATA PREPARATION
ProjW = "+proj=longlat +ellps=WGS84 +no_defs" #map projection (EPSG: 4326)

# Load cleaned occurrence dataset and convert to spatial points df

library(sp)
species_occ <- read.table("occurrences.txt", header=TRUE)
xy <- species_occ [,2:3]
df <- species_occ 
Sp_occ_SPDF <- SpatialPointsDataFrame(coords=xy, data=df, proj4string = CRS(ProjW)) # convert to SPDF


library(rangeBuilder)

# Thin occurrences - by changing the value for distance, you can set the resolution for thinning
# note that 20 = 20km and may need to be reduced for rarer species (with prevalence of <0.1)
# You should only change the value in the code which defines the value 'distance'

# IMPORTANT:
# start with 20km but if this reduces the number of remaining presences to <100, try 10km or 5km resolution

# be prepared for this to take a long time!!!

distance <- as.numeric(20)
thinned_occ <- filterByProximity(xy, dist=paste0(distance), mapUnits = FALSE, returnIndex = FALSE)
str(thinned_occ)

# Plot on SDM raster to check how many records will be used in modelling
# Load raster
PC1 <- raster("PCA_out_current/__PC1.tif")
plot(PC1)
Points_on_rast <- extract(PC1, thinned_occ) #?!# changed to 'thinned_occ' from 'Occ'
Removed_NAs <- as.data.frame(na.omit(Points_on_rast)) # this is the total number that will be used in the SDM


# IF NUMBER OF REMAINING PRESENCES IS STILL LESS THAN 100 at 1km resolution, 
# consider omitting this species from the dataset

# view points on world map
library(maps)
world <- map("world", fill = FALSE)
points(thinned_occ, pch = 16, col = "blue")

# view points on ecoregions map
ecoregion <- readOGR("ecoreg.shp", p4s = ProjW)
plot(ecoregion)
points(thinned_occ, pch = 16, col = "blue")

write.csv(thinned_occ, file=paste("presences",".csv",sep=""))

# Create index showing which records were removed
# may also take a long time!!
removed_points <- filterByProximity(xy, dist=paste0(distance), mapUnits = FALSE, returnIndex = TRUE)
# Add row number column called ID
xy$ID <- seq.int(nrow(xy))


# Get Long Lat of removed occurrences (these will be used for 'independent' evaluation in continuous boyce index)
library(dplyr)
Removed_occs <- filter(xy, ID %in% removed_points)

write.csv(Removed_occs, file=paste("removed_occs", ".csv",sep=""))

# following assigns the values to the number of original occurrences and 
# number remaining after thinning as csv file in working directory
# needed for reporting

pre_thinned_occs<- length(species_occ$species_name)
thinned_occ_df <- as.data.frame(thinned_occ)
post_thinned_occs <- length(thinned_occ_df$decimalLongitude)
removed_occ <- length(Removed_occs$ID)

# following saves minimum allowable distance, number of occurrences and 
# number remaining after thinning as csv file in working directory

write.table(distance, file="thinning_variables.csv",sep=",", 
            row.names = c("Distance_km"), col.names = paste0(species))
write.table(pre_thinned_occs, file="thinning_variables.csv",sep=",", 
            row.names = c("All_occurrences"), col.names = FALSE, append = TRUE)
write.table(post_thinned_occs, file="thinning_variables.csv",sep=",", 
            row.names = c("Remaining_occurrences"), col.names = FALSE, append = TRUE)
write.table(removed_occ, file="thinning_variables.csv",sep=",", 
            row.names = c("Removed_occurrences"), col.names = FALSE, append = TRUE)
