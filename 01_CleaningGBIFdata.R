# 01 Cleaning GBIF data

# Sarah Dalrymple
# July 2021


# following code removes all objects from the environment - 
# good practice to run this if you've been working on other projects this session
rm(list=ls()) #Removes all objects from working environment


### Manually download the GBIF occurrence data

###       Download the occurrence data as 'simple' csv from GBIF
###       save to a folder assigned the species name in the format 'Genus_species'
###       this folder is the working directory and was created using script '00 Setting up files and folders...'
###       save the citation in a txt file named 'GBIFcitation.txt' also in the working directory

# use the following line of code to navigate to your species working directory (if you haven't already)
setwd(choose.dir())
# use the getwd() function to display the file pathway in the console and confirm you're in the correct wd
getwd()
# paste the full file pathway from the console into the setwd() function for future reference
setwd("C:/Users/Public/Documents/Modelling/...")

# create object called 'species' to paste species name into code rather than retyping it
wd <-paste0(getwd())
species <- basename(wd)

csv_file_name <- list.files(pattern=".csv")

GBIFdata <- read.table(csv_file_name, sep="\t",fill=TRUE,header=TRUE, quote="",encoding="UTF-8")


### clean up the GBIF data

# remove rows without coordinates
GBIFdata_GeoRef <- GBIFdata[!(is.na(GBIFdata$decimalLatitude) | GBIFdata$decimalLatitude==""),]

# remove records pre-1970
GBIFdata_GeoRef_post1970 <- GBIFdata_GeoRef[GBIFdata_GeoRef$year >= 1970,]

# remove coordinates with greater uncertainty than 1km
GBIFdata_GeoRef_post1970_1kmcertainty <- 
  GBIFdata_GeoRef_post1970[GBIFdata_GeoRef_post1970$coordinateUncertaintyInMeters <= 1000,]

# remove coordinates reported to a precision of less than 2 d.p.
GBIFdata_GeoRef_post1970_1kmcertainty_2dp <- 
  GBIFdata_GeoRef_post1970_1kmcertainty[grep("\\.[0-9][0-9]", GBIFdata_GeoRef_post1970_1kmcertainty$decimalLatitude), ] 

# create a dataframe consisting of three columns: species, decimalLongiude, decimalLatitude
occ <- as.data.frame(cbind(GBIFdata_GeoRef_post1970_1kmcertainty_2dp$decimalLongitude, 
                           GBIFdata_GeoRef_post1970_1kmcertainty_2dp$decimalLatitude))
names(occ) <- c("decimalLongitude", "decimalLatitude")

species_name <- rep(species, length.out = nrow(occ))
coords <-cbind(species_name, occ)
head(coords)

library(CoordinateCleaner)

# help available at https://ropensci.github.io/CoordinateCleaner/

flags <- clean_coordinates(coords, 
                           lon = "decimalLongitude", lat = "decimalLatitude",
                           species = "species_name", countries = NULL, 
                           tests = c("capitals", "centroids", "equal", "gbif", "institutions", "outliers", "seas",
                                     "zeros"), capitals_rad = 10000, centroids_rad = 1000,
                           centroids_detail = "both", inst_rad = 100,
                           outliers_method = "quantile", outliers_mtp = 5, outliers_td = 1000,
                           outliers_size = 7, range_rad = 0, zeros_rad = 0.5,
                           capitals_ref = NULL, centroids_ref = NULL, country_ref = NULL,
                           inst_ref = NULL, range_ref = NULL, seas_ref = NULL,
                           seas_scale = 50, urban_ref = NULL, value = "spatialvalid",
                           verbose = TRUE, report = FALSE)
Occ_cleaned <- coords[flags$.summary,]

write.table(Occ_cleaned, "occurrences.txt", 
            row.names = FALSE, 
            quote = FALSE, 
            col.names = TRUE)

# Load cleaned occurrence dataset and convert to spatial points df
ProjW = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" #map projection (EPSG: 4326)
xy <- Occ_cleaned [,2:3]
df <- Occ_cleaned 
Sp_occ_SPDF <- SpatialPointsDataFrame(coords=xy, data=df, proj4string = CRS(ProjW)) # convert to SPDF

# view points on map
library(maps)

map("world", fill = FALSE)
points(Sp_occ_SPDF, pch = 16, col = "red")
