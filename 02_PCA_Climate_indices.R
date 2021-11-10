# 02 Deriving PCA climate indices

# Joe Bellis
# June 2021

# THIS SCRIPT REQUIRES A LOT OF MEMORY #
# IF YOU RUN INTO PROBLEMS TRY INCREASING THE MEMORY ALLOCATION to R #
# you can do this using the last bit of code in script '00 setting up files and folders...'

library(raster)
library(sp)
library(maptools)
library(rgdal)
library(RStoolbox)

rm(list=ls())

# use the following line of code to navigate to your species working directory (if you haven't already)
setwd(choose.dir())
# use the getwd() function to display the file pathway in the console and confirm you're in the correct wd
getwd()
# paste the full file pathway from the console into the setwd() function for future reference
setwd("C:/Users/Public/Documents/Modelling/Melampyrum_sylvaticum")

wd <-paste0(getwd())
species <- basename(wd)

# Load in rasters
bioclim_all <- stack(list.files("C:\\Users\\Public\\Documents\\Modelling\\wc2.1_2.5m_bioc_current", full.names = T))

# Set data projection to the standard WGS84
ProjW = "+proj=longlat +ellps=WGS84 +no_defs" #map projection (EPSG: 4326)

# Load terrestrial ecoregions shapefile with a 10km buffer to ensure coastal cells are preserved
geogrextent<-readOGR("C:\\Users\\Public\\Documents\\Modelling\\wwf_ecoregions_shp\\wwf_terr_ecos_10kmbuffered.shp", p4s = ProjW)


# Load cleaned occurrence dataset and convert to spatial points df
species_occ <- read.table("occurrences.txt", header=TRUE)
xy <- species_occ [,2:3]
df <- species_occ 
Sp_occ_SPDF <- SpatialPointsDataFrame(coords=xy, data=df, proj4string = CRS(ProjW))

# Select ecoregions where species occurs
ecoreg_sp<-over(Sp_occ_SPDF,geogrextent) #extract ecoregions of each locality point. 
uniqueEcoreg<-as.character(na.omit(unique(ecoreg_sp$ECO_NAME)))# unique list of occupied ecoregions
CurrentProj<- subset(geogrextent, geogrextent$ECO_NAME %in% uniqueEcoreg)#polygons of occupied ecoregions
plot(CurrentProj)

# Write the shapefile to save for future reference
writeOGR(CurrentProj, layer = "ecoreg", "ecoreg.shp" ,driver="ESRI Shapefile")

# Crop and mask
CropClim=crop(bioclim_all, CurrentProj) #Crop before mask to speed up processing time
MaskClim=mask(CropClim, CurrentProj)

# Standardise each one of the masked rasters based on Z standardisation (both center and scale need to be TRUE for this)
bioclim_stand <- scale(MaskClim, center = TRUE, scale = TRUE)

# Do the PCA using the RStoolbox package
PCA_Out <- rasterPCA(img = bioclim_stand, nSamples = NULL, ncomp = 19, spca = FALSE)
PCA_Out

dir.create("PCA_out_current")
setwd("PCA_out_current")
getwd()
# Write all principal components (a subset will be used for current conditions in SDMs)
writeRaster(PCA_Out$map, "_", format = 'GTiff', bylayer = TRUE, 
            c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18', 'PC19'))

# Variance explained
variance_explained <- summary(PCA_Out$model)
variance_explained

# Function for writing variance explained (see https://stackoverflow.com/questions/58395013/how-to-print-write-or-export-princomp-results-to-a-csv-file)
as.data.frame.summary.princomp <- function(x, ...) {
  vars <- x$sdev^2
  vars <- vars/sum(vars)
  type.convert(
    as.data.frame(
      rbind(`Standard deviation` = x$sdev, `Proportion of Variance` = vars, 
            `Cumulative Proportion` = cumsum(vars))
    )
  )
}

setwd(paste0("C:/Users/Public/Documents/Modelling/", species))
# Write the SD and proportion of variance output
write.csv(summary(PCA_Out$model),file ="variance.csv")


####################### Now prepare PCA outputs for obtaining future climate projection rasters ########################################
# Get the Loadings
PCA_load <- loadings(PCA_Out$model)
PCA_Loadings <- print(PCA_load, cutoff = 0.0001) # this number just ensures all values are retained

# Now we want to multiply the loadings for each PC variable by the future values of the climate variables
# Convert the atomic vector to a df
PCA_Loadings_df <- as.data.frame(unclass(PCA_Loadings))

# Rearrange order of row names to go from 1-19
row_order <- c("wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_2", "wc2.1_2.5m_bio_3", "wc2.1_2.5m_bio_4", "wc2.1_2.5m_bio_5", "wc2.1_2.5m_bio_6", "wc2.1_2.5m_bio_7", "wc2.1_2.5m_bio_8", "wc2.1_2.5m_bio_9", "wc2.1_2.5m_bio_10", "wc2.1_2.5m_bio_11",
               "wc2.1_2.5m_bio_12", "wc2.1_2.5m_bio_13", "wc2.1_2.5m_bio_14", "wc2.1_2.5m_bio_15", "wc2.1_2.5m_bio_16", "wc2.1_2.5m_bio_17", "wc2.1_2.5m_bio_18", "wc2.1_2.5m_bio_19")
PCA_Loadings_df <- PCA_Loadings_df[row_order, ]

# Write the loadings for potential future interpretation/use
write.csv(PCA_Loadings_df, file ="Loadings.csv")

# Extract loadings for each PC (first 5 should be enough)
PC1_load <- PCA_Loadings_df[1]
PC2_load <- PCA_Loadings_df[2]
PC3_load <- PCA_Loadings_df[3]
PC4_load <- PCA_Loadings_df[4]
PC5_load <- PCA_Loadings_df[5]

# Then get mean and SD of the masked area
x <- stack(MaskClim)
( x.stats <- data.frame(x.mean=cellStats(x, "mean")) )

# Combine mean and SD statistics as rows
( x.stats <- t(data.frame(x.mean=cellStats(x, "mean"))) )
( x.stats <- rbind(x.stats, t(data.frame(x.sd=cellStats(x, "sd")))) )
current.stats <- data.frame(x.stats)

# Rearrange order of column names to go from 1-19
col_order <- c("wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_2", "wc2.1_2.5m_bio_3", 
               "wc2.1_2.5m_bio_4", "wc2.1_2.5m_bio_5", "wc2.1_2.5m_bio_6", 
               "wc2.1_2.5m_bio_7", "wc2.1_2.5m_bio_8", "wc2.1_2.5m_bio_9", 
               "wc2.1_2.5m_bio_10", "wc2.1_2.5m_bio_11",
               "wc2.1_2.5m_bio_12", "wc2.1_2.5m_bio_13", "wc2.1_2.5m_bio_14", 
               "wc2.1_2.5m_bio_15", "wc2.1_2.5m_bio_16", "wc2.1_2.5m_bio_17", 
               "wc2.1_2.5m_bio_18", "wc2.1_2.5m_bio_19")
current.stats <- current.stats[, col_order]

# Load future projections
########################################## CanESM5 ##########################################
########################################## 2030 ############################################
########################################## 126 ##############################################

bioclim_Can_126_2030 <- stack("C:\\Users\\Public\\Documents\\Modelling\\wc2.1_2.5m_bioc_2021-2040\\wc2.1_2.5m_bioc_CanESM5_ssp126_2021-2040.tif")

# Crop and mask them to species-specific extent
CropClim_Can_126_2030=crop(bioclim_Can_126_2030, CurrentProj) #Crop before mask to speed up processing time
MaskClim_Can_126_2030=mask(CropClim_Can_126_2030, CurrentProj)

# Then standardise them based on the current mean and SD values
val1 <- MaskClim_Can_126_2030[[1]] - current.stats[1,1] 
Can_126_2030_sbio1 <- val1 / current.stats[2,1]

val2 <- MaskClim_Can_126_2030[[2]] - current.stats[1,2] 
Can_126_2030_sbio2 <- val2 / current.stats[2,2]

val3 <- MaskClim_Can_126_2030[[3]] - current.stats[1,3] 
Can_126_2030_sbio3 <- val3 / current.stats[2,3]

val4 <- MaskClim_Can_126_2030[[4]] - current.stats[1,4] 
Can_126_2030_sbio4 <- val4 / current.stats[2,4]

val5 <- MaskClim_Can_126_2030[[5]] - current.stats[1,5] 
Can_126_2030_sbio5 <- val5 / current.stats[2,5]

val6 <- MaskClim_Can_126_2030[[6]] - current.stats[1,6] 
Can_126_2030_sbio6 <- val6 / current.stats[2,6]

val7 <- MaskClim_Can_126_2030[[7]] - current.stats[1,7] 
Can_126_2030_sbio7 <- val7 / current.stats[2,7]

val8 <- MaskClim_Can_126_2030[[8]] - current.stats[1,8] 
Can_126_2030_sbio8 <- val8 / current.stats[2,8]

val9 <- MaskClim_Can_126_2030[[9]] - current.stats[1,9] 
Can_126_2030_sbio9 <- val9 / current.stats[2,9]

val10 <- MaskClim_Can_126_2030[[10]] - current.stats[1,10] 
Can_126_2030_sbio10 <- val10 / current.stats[2,10]

val11 <- MaskClim_Can_126_2030[[11]] - current.stats[1,11] 
Can_126_2030_sbio11 <- val11 / current.stats[2,11]

val12 <- MaskClim_Can_126_2030[[12]] - current.stats[1,12] 
Can_126_2030_sbio12 <- val12 / current.stats[2,12]

val13 <- MaskClim_Can_126_2030[[13]] - current.stats[1,13] 
Can_126_2030_sbio13 <- val13 / current.stats[2,13]

val14 <- MaskClim_Can_126_2030[[14]] - current.stats[1,14] 
Can_126_2030_sbio14 <- val14 / current.stats[2,14]

val15 <- MaskClim_Can_126_2030[[15]] - current.stats[1,15] 
Can_126_2030_sbio15 <- val15 / current.stats[2,15]

val16 <- MaskClim_Can_126_2030[[16]] - current.stats[1,16] 
Can_126_2030_sbio16 <- val16 / current.stats[2,16]

val17 <- MaskClim_Can_126_2030[[17]] - current.stats[1,17] 
Can_126_2030_sbio17 <- val17 / current.stats[2,17]

val18 <- MaskClim_Can_126_2030[[18]] - current.stats[1,18] 
Can_126_2030_sbio18 <- val18 / current.stats[2,18]

val19 <- MaskClim_Can_126_2030[[19]] - current.stats[1,19] 
Can_126_2030_sbio19 <- val19 / current.stats[2,19]

# Then multiply PC1 loadings for each variable by the standardised future rasters
can_126_2030_stand <- stack(Can_126_2030_sbio1, Can_126_2030_sbio2, Can_126_2030_sbio3, Can_126_2030_sbio4, Can_126_2030_sbio5, Can_126_2030_sbio6,
                            Can_126_2030_sbio7, Can_126_2030_sbio8, Can_126_2030_sbio9, Can_126_2030_sbio10, Can_126_2030_sbio11, Can_126_2030_sbio12,
                            Can_126_2030_sbio13, Can_126_2030_sbio14, Can_126_2030_sbio15, Can_126_2030_sbio16, Can_126_2030_sbio17, Can_126_2030_sbio18,
                            Can_126_2030_sbio19)

# PC1
PC1_bio1 <- can_126_2030_stand[[1]] * PC1_load[1,] 
PC1_bio2 <- can_126_2030_stand[[2]] * PC1_load[2,] 
PC1_bio3 <- can_126_2030_stand[[3]] * PC1_load[3,] 
PC1_bio4 <- can_126_2030_stand[[4]] * PC1_load[4,] 
PC1_bio5 <- can_126_2030_stand[[5]] * PC1_load[5,] 
PC1_bio6 <- can_126_2030_stand[[6]] * PC1_load[6,] 
PC1_bio7 <- can_126_2030_stand[[7]] * PC1_load[7,] 
PC1_bio8 <- can_126_2030_stand[[8]] * PC1_load[8,] 
PC1_bio9 <- can_126_2030_stand[[9]] * PC1_load[9,] 
PC1_bio10 <- can_126_2030_stand[[10]] * PC1_load[10,] 
PC1_bio11 <- can_126_2030_stand[[11]] * PC1_load[11,] 
PC1_bio12 <- can_126_2030_stand[[12]] * PC1_load[12,] 
PC1_bio13 <- can_126_2030_stand[[13]] * PC1_load[13,] 
PC1_bio14 <- can_126_2030_stand[[14]] * PC1_load[14,] 
PC1_bio15 <- can_126_2030_stand[[15]] * PC1_load[15,] 
PC1_bio16 <- can_126_2030_stand[[16]] * PC1_load[16,] 
PC1_bio17 <- can_126_2030_stand[[17]] * PC1_load[17,] 
PC1_bio18 <- can_126_2030_stand[[18]] * PC1_load[18,] 
PC1_bio19 <- can_126_2030_stand[[19]] * PC1_load[19,] 

PC1_bio_can_126_2030 <- PC1_bio1 + PC1_bio2 + PC1_bio3 + PC1_bio4 + PC1_bio5 + PC1_bio6 + PC1_bio7 + PC1_bio8 + PC1_bio9 + PC1_bio10 + PC1_bio11 + PC1_bio12 + PC1_bio13 + PC1_bio14 + PC1_bio15 + PC1_bio16 + PC1_bio17 + PC1_bio18 + PC1_bio19

# PC2
PC2_bio1 <- can_126_2030_stand[[1]] * PC2_load[1,] 
PC2_bio2 <- can_126_2030_stand[[2]] * PC2_load[2,] 
PC2_bio3 <- can_126_2030_stand[[3]] * PC2_load[3,] 
PC2_bio4 <- can_126_2030_stand[[4]] * PC2_load[4,] 
PC2_bio5 <- can_126_2030_stand[[5]] * PC2_load[5,] 
PC2_bio6 <- can_126_2030_stand[[6]] * PC2_load[6,] 
PC2_bio7 <- can_126_2030_stand[[7]] * PC2_load[7,] 
PC2_bio8 <- can_126_2030_stand[[8]] * PC2_load[8,] 
PC2_bio9 <- can_126_2030_stand[[9]] * PC2_load[9,] 
PC2_bio10 <- can_126_2030_stand[[10]] * PC2_load[10,] 
PC2_bio11 <- can_126_2030_stand[[11]] * PC2_load[11,] 
PC2_bio12 <- can_126_2030_stand[[12]] * PC2_load[12,] 
PC2_bio13 <- can_126_2030_stand[[13]] * PC2_load[13,] 
PC2_bio14 <- can_126_2030_stand[[14]] * PC2_load[14,] 
PC2_bio15 <- can_126_2030_stand[[15]] * PC2_load[15,] 
PC2_bio16 <- can_126_2030_stand[[16]] * PC2_load[16,] 
PC2_bio17 <- can_126_2030_stand[[17]] * PC2_load[17,] 
PC2_bio18 <- can_126_2030_stand[[18]] * PC2_load[18,] 
PC2_bio19 <- can_126_2030_stand[[19]] * PC2_load[19,] 

PC2_bio_can_126_2030 <- PC2_bio1 + PC2_bio2 + PC2_bio3 + PC2_bio4 + PC2_bio5 + PC2_bio6 + PC2_bio7 + PC2_bio8 + PC2_bio9 + PC2_bio10 + PC2_bio11 + PC2_bio12 + PC2_bio13 + PC2_bio14 + PC2_bio15 + PC2_bio16 + PC2_bio17 + PC2_bio18 + PC2_bio19

# PC3
PC3_bio1 <- can_126_2030_stand[[1]] * PC3_load[1,] 
PC3_bio2 <- can_126_2030_stand[[2]] * PC3_load[2,] 
PC3_bio3 <- can_126_2030_stand[[3]] * PC3_load[3,] 
PC3_bio4 <- can_126_2030_stand[[4]] * PC3_load[4,] 
PC3_bio5 <- can_126_2030_stand[[5]] * PC3_load[5,] 
PC3_bio6 <- can_126_2030_stand[[6]] * PC3_load[6,] 
PC3_bio7 <- can_126_2030_stand[[7]] * PC3_load[7,] 
PC3_bio8 <- can_126_2030_stand[[8]] * PC3_load[8,] 
PC3_bio9 <- can_126_2030_stand[[9]] * PC3_load[9,] 
PC3_bio10 <- can_126_2030_stand[[10]] * PC3_load[10,] 
PC3_bio11 <- can_126_2030_stand[[11]] * PC3_load[11,] 
PC3_bio12 <- can_126_2030_stand[[12]] * PC3_load[12,] 
PC3_bio13 <- can_126_2030_stand[[13]] * PC3_load[13,] 
PC3_bio14 <- can_126_2030_stand[[14]] * PC3_load[14,] 
PC3_bio15 <- can_126_2030_stand[[15]] * PC3_load[15,] 
PC3_bio16 <- can_126_2030_stand[[16]] * PC3_load[16,] 
PC3_bio17 <- can_126_2030_stand[[17]] * PC3_load[17,] 
PC3_bio18 <- can_126_2030_stand[[18]] * PC3_load[18,] 
PC3_bio19 <- can_126_2030_stand[[19]] * PC3_load[19,] 

PC3_bio_can_126_2030 <- PC3_bio1 + PC3_bio2 + PC3_bio3 + PC3_bio4 + PC3_bio5 + PC3_bio6 + PC3_bio7 + PC3_bio8 + PC3_bio9 + PC3_bio10 + PC3_bio11 + PC3_bio12 + PC3_bio13 + PC3_bio14 + PC3_bio15 + PC3_bio16 + PC3_bio17 + PC3_bio18 + PC3_bio19

# PC4
PC4_bio1 <- can_126_2030_stand[[1]] * PC4_load[1,] 
PC4_bio2 <- can_126_2030_stand[[2]] * PC4_load[2,] 
PC4_bio3 <- can_126_2030_stand[[3]] * PC4_load[3,] 
PC4_bio4 <- can_126_2030_stand[[4]] * PC4_load[4,] 
PC4_bio5 <- can_126_2030_stand[[5]] * PC4_load[5,] 
PC4_bio6 <- can_126_2030_stand[[6]] * PC4_load[6,] 
PC4_bio7 <- can_126_2030_stand[[7]] * PC4_load[7,] 
PC4_bio8 <- can_126_2030_stand[[8]] * PC4_load[8,] 
PC4_bio9 <- can_126_2030_stand[[9]] * PC4_load[9,] 
PC4_bio10 <- can_126_2030_stand[[10]] * PC4_load[10,] 
PC4_bio11 <- can_126_2030_stand[[11]] * PC4_load[11,] 
PC4_bio12 <- can_126_2030_stand[[12]] * PC4_load[12,] 
PC4_bio13 <- can_126_2030_stand[[13]] * PC4_load[13,] 
PC4_bio14 <- can_126_2030_stand[[14]] * PC4_load[14,] 
PC4_bio15 <- can_126_2030_stand[[15]] * PC4_load[15,] 
PC4_bio16 <- can_126_2030_stand[[16]] * PC4_load[16,] 
PC4_bio17 <- can_126_2030_stand[[17]] * PC4_load[17,] 
PC4_bio18 <- can_126_2030_stand[[18]] * PC4_load[18,] 
PC4_bio19 <- can_126_2030_stand[[19]] * PC4_load[19,] 

PC4_bio_can_126_2030 <- PC4_bio1 + PC4_bio2 + PC4_bio3 + PC4_bio4 + PC4_bio5 + PC4_bio6 + PC4_bio7 + PC4_bio8 + PC4_bio9 + PC4_bio10 + PC4_bio11 + PC4_bio12 + PC4_bio13 + PC4_bio14 + PC4_bio15 + PC4_bio16 + PC4_bio17 + PC4_bio18 + PC4_bio19

# PC5
PC5_bio1 <- can_126_2030_stand[[1]] * PC5_load[1,] 
PC5_bio2 <- can_126_2030_stand[[2]] * PC5_load[2,] 
PC5_bio3 <- can_126_2030_stand[[3]] * PC5_load[3,] 
PC5_bio4 <- can_126_2030_stand[[4]] * PC5_load[4,] 
PC5_bio5 <- can_126_2030_stand[[5]] * PC5_load[5,] 
PC5_bio6 <- can_126_2030_stand[[6]] * PC5_load[6,] 
PC5_bio7 <- can_126_2030_stand[[7]] * PC5_load[7,] 
PC5_bio8 <- can_126_2030_stand[[8]] * PC5_load[8,] 
PC5_bio9 <- can_126_2030_stand[[9]] * PC5_load[9,] 
PC5_bio10 <- can_126_2030_stand[[10]] * PC5_load[10,] 
PC5_bio11 <- can_126_2030_stand[[11]] * PC5_load[11,] 
PC5_bio12 <- can_126_2030_stand[[12]] * PC5_load[12,] 
PC5_bio13 <- can_126_2030_stand[[13]] * PC5_load[13,] 
PC5_bio14 <- can_126_2030_stand[[14]] * PC5_load[14,] 
PC5_bio15 <- can_126_2030_stand[[15]] * PC5_load[15,] 
PC5_bio16 <- can_126_2030_stand[[16]] * PC5_load[16,] 
PC5_bio17 <- can_126_2030_stand[[17]] * PC5_load[17,] 
PC5_bio18 <- can_126_2030_stand[[18]] * PC5_load[18,] 
PC5_bio19 <- can_126_2030_stand[[19]] * PC5_load[19,] 

PC5_bio_can_126_2030 <- PC5_bio1 + PC5_bio2 + PC5_bio3 + PC5_bio4 + PC5_bio5 + PC5_bio6 + PC5_bio7 + PC5_bio8 + PC5_bio9 + PC5_bio10 + PC5_bio11 + PC5_bio12 + PC5_bio13 + PC5_bio14 + PC5_bio15 + PC5_bio16 + PC5_bio17 + PC5_bio18 + PC5_bio19
getwd()
dir.create("PCA_out_bio_can_126_2030")
setwd("PCA_out_bio_can_126_2030")
getwd()

# Write the rasters for the 5 PCs
writeRaster(PC1_bio_can_126_2030, "__PC1.tif")
writeRaster(PC2_bio_can_126_2030, "__PC2.tif")
writeRaster(PC3_bio_can_126_2030, "__PC3.tif")
writeRaster(PC4_bio_can_126_2030, "__PC4.tif")
writeRaster(PC5_bio_can_126_2030, "__PC5.tif")

setwd(paste0("C:/Users/Public/Documents/Modelling/", species))

########################################## CanESM5 ##########################################
########################################## 2030 ############################################
########################################## 370 ##############################################
bioclim_Can_370_2030 <- stack("C:\\Users\\Public\\Documents\\Modelling\\wc2.1_2.5m_bioc_2021-2040\\wc2.1_2.5m_bioc_CanESM5_ssp370_2021-2040.tif")

# Crop and mask them to species-specific extent
CropClim_Can_370_2030=crop(bioclim_Can_370_2030, CurrentProj) #Crop before mask to speed up processing time
MaskClim_Can_370_2030=mask(CropClim_Can_370_2030, CurrentProj)

# Then standardise them based on the current mean and SD values
val1 <- MaskClim_Can_370_2030[[1]] - current.stats[1,1] 
Can_370_2030_sbio1 <- val1 / current.stats[2,1]

val2 <- MaskClim_Can_370_2030[[2]] - current.stats[1,2] 
Can_370_2030_sbio2 <- val2 / current.stats[2,2]

val3 <- MaskClim_Can_370_2030[[3]] - current.stats[1,3] 
Can_370_2030_sbio3 <- val3 / current.stats[2,3]

val4 <- MaskClim_Can_370_2030[[4]] - current.stats[1,4] 
Can_370_2030_sbio4 <- val4 / current.stats[2,4]

val5 <- MaskClim_Can_370_2030[[5]] - current.stats[1,5] 
Can_370_2030_sbio5 <- val5 / current.stats[2,5]

val6 <- MaskClim_Can_370_2030[[6]] - current.stats[1,6] 
Can_370_2030_sbio6 <- val6 / current.stats[2,6]

val7 <- MaskClim_Can_370_2030[[7]] - current.stats[1,7] 
Can_370_2030_sbio7 <- val7 / current.stats[2,7]

val8 <- MaskClim_Can_370_2030[[8]] - current.stats[1,8] 
Can_370_2030_sbio8 <- val8 / current.stats[2,8]

val9 <- MaskClim_Can_370_2030[[9]] - current.stats[1,9] 
Can_370_2030_sbio9 <- val9 / current.stats[2,9]

val10 <- MaskClim_Can_370_2030[[10]] - current.stats[1,10] 
Can_370_2030_sbio10 <- val10 / current.stats[2,10]

val11 <- MaskClim_Can_370_2030[[11]] - current.stats[1,11] 
Can_370_2030_sbio11 <- val11 / current.stats[2,11]

val12 <- MaskClim_Can_370_2030[[12]] - current.stats[1,12] 
Can_370_2030_sbio12 <- val12 / current.stats[2,12]

val13 <- MaskClim_Can_370_2030[[13]] - current.stats[1,13] 
Can_370_2030_sbio13 <- val13 / current.stats[2,13]

val14 <- MaskClim_Can_370_2030[[14]] - current.stats[1,14] 
Can_370_2030_sbio14 <- val14 / current.stats[2,14]

val15 <- MaskClim_Can_370_2030[[15]] - current.stats[1,15] 
Can_370_2030_sbio15 <- val15 / current.stats[2,15]

val16 <- MaskClim_Can_370_2030[[16]] - current.stats[1,16] 
Can_370_2030_sbio16 <- val16 / current.stats[2,16]

val17 <- MaskClim_Can_370_2030[[17]] - current.stats[1,17] 
Can_370_2030_sbio17 <- val17 / current.stats[2,17]

val18 <- MaskClim_Can_370_2030[[18]] - current.stats[1,18] 
Can_370_2030_sbio18 <- val18 / current.stats[2,18]

val19 <- MaskClim_Can_370_2030[[19]] - current.stats[1,19] 
Can_370_2030_sbio19 <- val19 / current.stats[2,19]

# Then multiply PC1 loadings for each variable by the standardised future rasters
can_370_2030_stand <- stack(Can_370_2030_sbio1, Can_370_2030_sbio2, Can_370_2030_sbio3, Can_370_2030_sbio4, Can_370_2030_sbio5, Can_370_2030_sbio6,
                            Can_370_2030_sbio7, Can_370_2030_sbio8, Can_370_2030_sbio9, Can_370_2030_sbio10, Can_370_2030_sbio11, Can_370_2030_sbio12,
                            Can_370_2030_sbio13, Can_370_2030_sbio14, Can_370_2030_sbio15, Can_370_2030_sbio16, Can_370_2030_sbio17, Can_370_2030_sbio18,
                            Can_370_2030_sbio19)

# PC1
PC1_bio1 <- can_370_2030_stand[[1]] * PC1_load[1,] 
PC1_bio2 <- can_370_2030_stand[[2]] * PC1_load[2,] 
PC1_bio3 <- can_370_2030_stand[[3]] * PC1_load[3,] 
PC1_bio4 <- can_370_2030_stand[[4]] * PC1_load[4,] 
PC1_bio5 <- can_370_2030_stand[[5]] * PC1_load[5,] 
PC1_bio6 <- can_370_2030_stand[[6]] * PC1_load[6,] 
PC1_bio7 <- can_370_2030_stand[[7]] * PC1_load[7,] 
PC1_bio8 <- can_370_2030_stand[[8]] * PC1_load[8,] 
PC1_bio9 <- can_370_2030_stand[[9]] * PC1_load[9,] 
PC1_bio10 <- can_370_2030_stand[[10]] * PC1_load[10,] 
PC1_bio11 <- can_370_2030_stand[[11]] * PC1_load[11,] 
PC1_bio12 <- can_370_2030_stand[[12]] * PC1_load[12,] 
PC1_bio13 <- can_370_2030_stand[[13]] * PC1_load[13,] 
PC1_bio14 <- can_370_2030_stand[[14]] * PC1_load[14,] 
PC1_bio15 <- can_370_2030_stand[[15]] * PC1_load[15,] 
PC1_bio16 <- can_370_2030_stand[[16]] * PC1_load[16,] 
PC1_bio17 <- can_370_2030_stand[[17]] * PC1_load[17,] 
PC1_bio18 <- can_370_2030_stand[[18]] * PC1_load[18,] 
PC1_bio19 <- can_370_2030_stand[[19]] * PC1_load[19,] 

PC1_bio_can_370_2030 <- PC1_bio1 + PC1_bio2 + PC1_bio3 + PC1_bio4 + PC1_bio5 + PC1_bio6 + PC1_bio7 + PC1_bio8 + PC1_bio9 + PC1_bio10 + PC1_bio11 + PC1_bio12 + PC1_bio13 + PC1_bio14 + PC1_bio15 + PC1_bio16 + PC1_bio17 + PC1_bio18 + PC1_bio19

# PC2
PC2_bio1 <- can_370_2030_stand[[1]] * PC2_load[1,] 
PC2_bio2 <- can_370_2030_stand[[2]] * PC2_load[2,] 
PC2_bio3 <- can_370_2030_stand[[3]] * PC2_load[3,] 
PC2_bio4 <- can_370_2030_stand[[4]] * PC2_load[4,] 
PC2_bio5 <- can_370_2030_stand[[5]] * PC2_load[5,] 
PC2_bio6 <- can_370_2030_stand[[6]] * PC2_load[6,] 
PC2_bio7 <- can_370_2030_stand[[7]] * PC2_load[7,] 
PC2_bio8 <- can_370_2030_stand[[8]] * PC2_load[8,] 
PC2_bio9 <- can_370_2030_stand[[9]] * PC2_load[9,] 
PC2_bio10 <- can_370_2030_stand[[10]] * PC2_load[10,] 
PC2_bio11 <- can_370_2030_stand[[11]] * PC2_load[11,] 
PC2_bio12 <- can_370_2030_stand[[12]] * PC2_load[12,] 
PC2_bio13 <- can_370_2030_stand[[13]] * PC2_load[13,] 
PC2_bio14 <- can_370_2030_stand[[14]] * PC2_load[14,] 
PC2_bio15 <- can_370_2030_stand[[15]] * PC2_load[15,] 
PC2_bio16 <- can_370_2030_stand[[16]] * PC2_load[16,] 
PC2_bio17 <- can_370_2030_stand[[17]] * PC2_load[17,] 
PC2_bio18 <- can_370_2030_stand[[18]] * PC2_load[18,] 
PC2_bio19 <- can_370_2030_stand[[19]] * PC2_load[19,] 

PC2_bio_can_370_2030 <- PC2_bio1 + PC2_bio2 + PC2_bio3 + PC2_bio4 + PC2_bio5 + PC2_bio6 + PC2_bio7 + PC2_bio8 + PC2_bio9 + PC2_bio10 + PC2_bio11 + PC2_bio12 + PC2_bio13 + PC2_bio14 + PC2_bio15 + PC2_bio16 + PC2_bio17 + PC2_bio18 + PC2_bio19

# PC3
PC3_bio1 <- can_370_2030_stand[[1]] * PC3_load[1,] 
PC3_bio2 <- can_370_2030_stand[[2]] * PC3_load[2,] 
PC3_bio3 <- can_370_2030_stand[[3]] * PC3_load[3,] 
PC3_bio4 <- can_370_2030_stand[[4]] * PC3_load[4,] 
PC3_bio5 <- can_370_2030_stand[[5]] * PC3_load[5,] 
PC3_bio6 <- can_370_2030_stand[[6]] * PC3_load[6,] 
PC3_bio7 <- can_370_2030_stand[[7]] * PC3_load[7,] 
PC3_bio8 <- can_370_2030_stand[[8]] * PC3_load[8,] 
PC3_bio9 <- can_370_2030_stand[[9]] * PC3_load[9,] 
PC3_bio10 <- can_370_2030_stand[[10]] * PC3_load[10,] 
PC3_bio11 <- can_370_2030_stand[[11]] * PC3_load[11,] 
PC3_bio12 <- can_370_2030_stand[[12]] * PC3_load[12,] 
PC3_bio13 <- can_370_2030_stand[[13]] * PC3_load[13,] 
PC3_bio14 <- can_370_2030_stand[[14]] * PC3_load[14,] 
PC3_bio15 <- can_370_2030_stand[[15]] * PC3_load[15,] 
PC3_bio16 <- can_370_2030_stand[[16]] * PC3_load[16,] 
PC3_bio17 <- can_370_2030_stand[[17]] * PC3_load[17,] 
PC3_bio18 <- can_370_2030_stand[[18]] * PC3_load[18,] 
PC3_bio19 <- can_370_2030_stand[[19]] * PC3_load[19,] 

PC3_bio_can_370_2030 <- PC3_bio1 + PC3_bio2 + PC3_bio3 + PC3_bio4 + PC3_bio5 + PC3_bio6 + PC3_bio7 + PC3_bio8 + PC3_bio9 + PC3_bio10 + PC3_bio11 + PC3_bio12 + PC3_bio13 + PC3_bio14 + PC3_bio15 + PC3_bio16 + PC3_bio17 + PC3_bio18 + PC3_bio19

# PC4
PC4_bio1 <- can_370_2030_stand[[1]] * PC4_load[1,] 
PC4_bio2 <- can_370_2030_stand[[2]] * PC4_load[2,] 
PC4_bio3 <- can_370_2030_stand[[3]] * PC4_load[3,] 
PC4_bio4 <- can_370_2030_stand[[4]] * PC4_load[4,] 
PC4_bio5 <- can_370_2030_stand[[5]] * PC4_load[5,] 
PC4_bio6 <- can_370_2030_stand[[6]] * PC4_load[6,] 
PC4_bio7 <- can_370_2030_stand[[7]] * PC4_load[7,] 
PC4_bio8 <- can_370_2030_stand[[8]] * PC4_load[8,] 
PC4_bio9 <- can_370_2030_stand[[9]] * PC4_load[9,] 
PC4_bio10 <- can_370_2030_stand[[10]] * PC4_load[10,] 
PC4_bio11 <- can_370_2030_stand[[11]] * PC4_load[11,] 
PC4_bio12 <- can_370_2030_stand[[12]] * PC4_load[12,] 
PC4_bio13 <- can_370_2030_stand[[13]] * PC4_load[13,] 
PC4_bio14 <- can_370_2030_stand[[14]] * PC4_load[14,] 
PC4_bio15 <- can_370_2030_stand[[15]] * PC4_load[15,] 
PC4_bio16 <- can_370_2030_stand[[16]] * PC4_load[16,] 
PC4_bio17 <- can_370_2030_stand[[17]] * PC4_load[17,] 
PC4_bio18 <- can_370_2030_stand[[18]] * PC4_load[18,] 
PC4_bio19 <- can_370_2030_stand[[19]] * PC4_load[19,] 

PC4_bio_can_370_2030 <- PC4_bio1 + PC4_bio2 + PC4_bio3 + PC4_bio4 + PC4_bio5 + PC4_bio6 + PC4_bio7 + PC4_bio8 + PC4_bio9 + PC4_bio10 + PC4_bio11 + PC4_bio12 + PC4_bio13 + PC4_bio14 + PC4_bio15 + PC4_bio16 + PC4_bio17 + PC4_bio18 + PC4_bio19

# PC5
PC5_bio1 <- can_370_2030_stand[[1]] * PC5_load[1,] 
PC5_bio2 <- can_370_2030_stand[[2]] * PC5_load[2,] 
PC5_bio3 <- can_370_2030_stand[[3]] * PC5_load[3,] 
PC5_bio4 <- can_370_2030_stand[[4]] * PC5_load[4,] 
PC5_bio5 <- can_370_2030_stand[[5]] * PC5_load[5,] 
PC5_bio6 <- can_370_2030_stand[[6]] * PC5_load[6,] 
PC5_bio7 <- can_370_2030_stand[[7]] * PC5_load[7,] 
PC5_bio8 <- can_370_2030_stand[[8]] * PC5_load[8,] 
PC5_bio9 <- can_370_2030_stand[[9]] * PC5_load[9,] 
PC5_bio10 <- can_370_2030_stand[[10]] * PC5_load[10,] 
PC5_bio11 <- can_370_2030_stand[[11]] * PC5_load[11,] 
PC5_bio12 <- can_370_2030_stand[[12]] * PC5_load[12,] 
PC5_bio13 <- can_370_2030_stand[[13]] * PC5_load[13,] 
PC5_bio14 <- can_370_2030_stand[[14]] * PC5_load[14,] 
PC5_bio15 <- can_370_2030_stand[[15]] * PC5_load[15,] 
PC5_bio16 <- can_370_2030_stand[[16]] * PC5_load[16,] 
PC5_bio17 <- can_370_2030_stand[[17]] * PC5_load[17,] 
PC5_bio18 <- can_370_2030_stand[[18]] * PC5_load[18,] 
PC5_bio19 <- can_370_2030_stand[[19]] * PC5_load[19,] 

PC5_bio_can_370_2030 <- PC5_bio1 + PC5_bio2 + PC5_bio3 + PC5_bio4 + PC5_bio5 + PC5_bio6 + PC5_bio7 + PC5_bio8 + PC5_bio9 + PC5_bio10 + PC5_bio11 + PC5_bio12 + PC5_bio13 + PC5_bio14 + PC5_bio15 + PC5_bio16 + PC5_bio17 + PC5_bio18 + PC5_bio19

dir.create("PCA_out_bio_can_370_2030")
setwd("PCA_out_bio_can_370_2030")
getwd()

# Write the rasters for the 5 PCs
writeRaster(PC1_bio_can_370_2030, "__PC1.tif")
writeRaster(PC2_bio_can_370_2030, "__PC2.tif")
writeRaster(PC3_bio_can_370_2030, "__PC3.tif")
writeRaster(PC4_bio_can_370_2030, "__PC4.tif")
writeRaster(PC5_bio_can_370_2030, "__PC5.tif")

setwd(paste0("C:/Users/Public/Documents/Modelling/", species))

########################################## CNRM-CM6-1 ##########################################
########################################## 2030 ############################################
########################################## 126 ##############################################
bioclim_Cnrm_126_2030 <- stack("C:\\Users\\Public\\Documents\\Modelling\\wc2.1_2.5m_bioc_2021-2040\\wc2.1_2.5m_bioc_CNRM-CM6-1_ssp126_2021-2040.tif")

# Crop and mask them to species-specific extent
CropClim_Cnrm_126_2030=crop(bioclim_Cnrm_126_2030, CurrentProj) #Crop before mask to speed up processing time
MaskClim_Cnrm_126_2030=mask(CropClim_Cnrm_126_2030, CurrentProj)

# Then standardise them based on the current mean and SD values
val1 <- MaskClim_Cnrm_126_2030[[1]] - current.stats[1,1] 
Cnrm_126_2030_sbio1 <- val1 / current.stats[2,1]

val2 <- MaskClim_Cnrm_126_2030[[2]] - current.stats[1,2] 
Cnrm_126_2030_sbio2 <- val2 / current.stats[2,2]

val3 <- MaskClim_Cnrm_126_2030[[3]] - current.stats[1,3] 
Cnrm_126_2030_sbio3 <- val3 / current.stats[2,3]

val4 <- MaskClim_Cnrm_126_2030[[4]] - current.stats[1,4] 
Cnrm_126_2030_sbio4 <- val4 / current.stats[2,4]

val5 <- MaskClim_Cnrm_126_2030[[5]] - current.stats[1,5] 
Cnrm_126_2030_sbio5 <- val5 / current.stats[2,5]

val6 <- MaskClim_Cnrm_126_2030[[6]] - current.stats[1,6] 
Cnrm_126_2030_sbio6 <- val6 / current.stats[2,6]

val7 <- MaskClim_Cnrm_126_2030[[7]] - current.stats[1,7] 
Cnrm_126_2030_sbio7 <- val7 / current.stats[2,7]

val8 <- MaskClim_Cnrm_126_2030[[8]] - current.stats[1,8] 
Cnrm_126_2030_sbio8 <- val8 / current.stats[2,8]

val9 <- MaskClim_Cnrm_126_2030[[9]] - current.stats[1,9] 
Cnrm_126_2030_sbio9 <- val9 / current.stats[2,9]

val10 <- MaskClim_Cnrm_126_2030[[10]] - current.stats[1,10] 
Cnrm_126_2030_sbio10 <- val10 / current.stats[2,10]

val11 <- MaskClim_Cnrm_126_2030[[11]] - current.stats[1,11] 
Cnrm_126_2030_sbio11 <- val11 / current.stats[2,11]

val12 <- MaskClim_Cnrm_126_2030[[12]] - current.stats[1,12] 
Cnrm_126_2030_sbio12 <- val12 / current.stats[2,12]

val13 <- MaskClim_Cnrm_126_2030[[13]] - current.stats[1,13] 
Cnrm_126_2030_sbio13 <- val13 / current.stats[2,13]

val14 <- MaskClim_Cnrm_126_2030[[14]] - current.stats[1,14] 
Cnrm_126_2030_sbio14 <- val14 / current.stats[2,14]

val15 <- MaskClim_Cnrm_126_2030[[15]] - current.stats[1,15] 
Cnrm_126_2030_sbio15 <- val15 / current.stats[2,15]

val16 <- MaskClim_Cnrm_126_2030[[16]] - current.stats[1,16] 
Cnrm_126_2030_sbio16 <- val16 / current.stats[2,16]

val17 <- MaskClim_Cnrm_126_2030[[17]] - current.stats[1,17] 
Cnrm_126_2030_sbio17 <- val17 / current.stats[2,17]

val18 <- MaskClim_Cnrm_126_2030[[18]] - current.stats[1,18] 
Cnrm_126_2030_sbio18 <- val18 / current.stats[2,18]

val19 <- MaskClim_Cnrm_126_2030[[19]] - current.stats[1,19] 
Cnrm_126_2030_sbio19 <- val19 / current.stats[2,19]

# Then multiply PC1 loadings for each variable by the standardised future rasters
Cnrm_126_2030_stand <- stack(Cnrm_126_2030_sbio1, Cnrm_126_2030_sbio2, Cnrm_126_2030_sbio3, Cnrm_126_2030_sbio4, Cnrm_126_2030_sbio5, Cnrm_126_2030_sbio6,
                             Cnrm_126_2030_sbio7, Cnrm_126_2030_sbio8, Cnrm_126_2030_sbio9, Cnrm_126_2030_sbio10, Cnrm_126_2030_sbio11, Cnrm_126_2030_sbio12,
                             Cnrm_126_2030_sbio13, Cnrm_126_2030_sbio14, Cnrm_126_2030_sbio15, Cnrm_126_2030_sbio16, Cnrm_126_2030_sbio17, Cnrm_126_2030_sbio18,
                             Cnrm_126_2030_sbio19)

# PC1
PC1_bio1 <- Cnrm_126_2030_stand[[1]] * PC1_load[1,] 
PC1_bio2 <- Cnrm_126_2030_stand[[2]] * PC1_load[2,] 
PC1_bio3 <- Cnrm_126_2030_stand[[3]] * PC1_load[3,] 
PC1_bio4 <- Cnrm_126_2030_stand[[4]] * PC1_load[4,] 
PC1_bio5 <- Cnrm_126_2030_stand[[5]] * PC1_load[5,] 
PC1_bio6 <- Cnrm_126_2030_stand[[6]] * PC1_load[6,] 
PC1_bio7 <- Cnrm_126_2030_stand[[7]] * PC1_load[7,] 
PC1_bio8 <- Cnrm_126_2030_stand[[8]] * PC1_load[8,] 
PC1_bio9 <- Cnrm_126_2030_stand[[9]] * PC1_load[9,] 
PC1_bio10 <- Cnrm_126_2030_stand[[10]] * PC1_load[10,] 
PC1_bio11 <- Cnrm_126_2030_stand[[11]] * PC1_load[11,] 
PC1_bio12 <- Cnrm_126_2030_stand[[12]] * PC1_load[12,] 
PC1_bio13 <- Cnrm_126_2030_stand[[13]] * PC1_load[13,] 
PC1_bio14 <- Cnrm_126_2030_stand[[14]] * PC1_load[14,] 
PC1_bio15 <- Cnrm_126_2030_stand[[15]] * PC1_load[15,] 
PC1_bio16 <- Cnrm_126_2030_stand[[16]] * PC1_load[16,] 
PC1_bio17 <- Cnrm_126_2030_stand[[17]] * PC1_load[17,] 
PC1_bio18 <- Cnrm_126_2030_stand[[18]] * PC1_load[18,] 
PC1_bio19 <- Cnrm_126_2030_stand[[19]] * PC1_load[19,] 

PC1_bio_Cnrm_126_2030 <- PC1_bio1 + PC1_bio2 + PC1_bio3 + PC1_bio4 + PC1_bio5 + PC1_bio6 + PC1_bio7 + PC1_bio8 + PC1_bio9 + PC1_bio10 + PC1_bio11 + PC1_bio12 + PC1_bio13 + PC1_bio14 + PC1_bio15 + PC1_bio16 + PC1_bio17 + PC1_bio18 + PC1_bio19

# PC2
PC2_bio1 <- Cnrm_126_2030_stand[[1]] * PC2_load[1,] 
PC2_bio2 <- Cnrm_126_2030_stand[[2]] * PC2_load[2,] 
PC2_bio3 <- Cnrm_126_2030_stand[[3]] * PC2_load[3,] 
PC2_bio4 <- Cnrm_126_2030_stand[[4]] * PC2_load[4,] 
PC2_bio5 <- Cnrm_126_2030_stand[[5]] * PC2_load[5,] 
PC2_bio6 <- Cnrm_126_2030_stand[[6]] * PC2_load[6,] 
PC2_bio7 <- Cnrm_126_2030_stand[[7]] * PC2_load[7,] 
PC2_bio8 <- Cnrm_126_2030_stand[[8]] * PC2_load[8,] 
PC2_bio9 <- Cnrm_126_2030_stand[[9]] * PC2_load[9,] 
PC2_bio10 <- Cnrm_126_2030_stand[[10]] * PC2_load[10,] 
PC2_bio11 <- Cnrm_126_2030_stand[[11]] * PC2_load[11,] 
PC2_bio12 <- Cnrm_126_2030_stand[[12]] * PC2_load[12,] 
PC2_bio13 <- Cnrm_126_2030_stand[[13]] * PC2_load[13,] 
PC2_bio14 <- Cnrm_126_2030_stand[[14]] * PC2_load[14,] 
PC2_bio15 <- Cnrm_126_2030_stand[[15]] * PC2_load[15,] 
PC2_bio16 <- Cnrm_126_2030_stand[[16]] * PC2_load[16,] 
PC2_bio17 <- Cnrm_126_2030_stand[[17]] * PC2_load[17,] 
PC2_bio18 <- Cnrm_126_2030_stand[[18]] * PC2_load[18,] 
PC2_bio19 <- Cnrm_126_2030_stand[[19]] * PC2_load[19,] 

PC2_bio_Cnrm_126_2030 <- PC2_bio1 + PC2_bio2 + PC2_bio3 + PC2_bio4 + PC2_bio5 + PC2_bio6 + PC2_bio7 + PC2_bio8 + PC2_bio9 + PC2_bio10 + PC2_bio11 + PC2_bio12 + PC2_bio13 + PC2_bio14 + PC2_bio15 + PC2_bio16 + PC2_bio17 + PC2_bio18 + PC2_bio19

# PC3
PC3_bio1 <- Cnrm_126_2030_stand[[1]] * PC3_load[1,] 
PC3_bio2 <- Cnrm_126_2030_stand[[2]] * PC3_load[2,] 
PC3_bio3 <- Cnrm_126_2030_stand[[3]] * PC3_load[3,] 
PC3_bio4 <- Cnrm_126_2030_stand[[4]] * PC3_load[4,] 
PC3_bio5 <- Cnrm_126_2030_stand[[5]] * PC3_load[5,] 
PC3_bio6 <- Cnrm_126_2030_stand[[6]] * PC3_load[6,] 
PC3_bio7 <- Cnrm_126_2030_stand[[7]] * PC3_load[7,] 
PC3_bio8 <- Cnrm_126_2030_stand[[8]] * PC3_load[8,] 
PC3_bio9 <- Cnrm_126_2030_stand[[9]] * PC3_load[9,] 
PC3_bio10 <- Cnrm_126_2030_stand[[10]] * PC3_load[10,] 
PC3_bio11 <- Cnrm_126_2030_stand[[11]] * PC3_load[11,] 
PC3_bio12 <- Cnrm_126_2030_stand[[12]] * PC3_load[12,] 
PC3_bio13 <- Cnrm_126_2030_stand[[13]] * PC3_load[13,] 
PC3_bio14 <- Cnrm_126_2030_stand[[14]] * PC3_load[14,] 
PC3_bio15 <- Cnrm_126_2030_stand[[15]] * PC3_load[15,] 
PC3_bio16 <- Cnrm_126_2030_stand[[16]] * PC3_load[16,] 
PC3_bio17 <- Cnrm_126_2030_stand[[17]] * PC3_load[17,] 
PC3_bio18 <- Cnrm_126_2030_stand[[18]] * PC3_load[18,] 
PC3_bio19 <- Cnrm_126_2030_stand[[19]] * PC3_load[19,] 

PC3_bio_Cnrm_126_2030 <- PC3_bio1 + PC3_bio2 + PC3_bio3 + PC3_bio4 + PC3_bio5 + PC3_bio6 + PC3_bio7 + PC3_bio8 + PC3_bio9 + PC3_bio10 + PC3_bio11 + PC3_bio12 + PC3_bio13 + PC3_bio14 + PC3_bio15 + PC3_bio16 + PC3_bio17 + PC3_bio18 + PC3_bio19

# PC4
PC4_bio1 <- Cnrm_126_2030_stand[[1]] * PC4_load[1,] 
PC4_bio2 <- Cnrm_126_2030_stand[[2]] * PC4_load[2,] 
PC4_bio3 <- Cnrm_126_2030_stand[[3]] * PC4_load[3,] 
PC4_bio4 <- Cnrm_126_2030_stand[[4]] * PC4_load[4,] 
PC4_bio5 <- Cnrm_126_2030_stand[[5]] * PC4_load[5,] 
PC4_bio6 <- Cnrm_126_2030_stand[[6]] * PC4_load[6,] 
PC4_bio7 <- Cnrm_126_2030_stand[[7]] * PC4_load[7,] 
PC4_bio8 <- Cnrm_126_2030_stand[[8]] * PC4_load[8,] 
PC4_bio9 <- Cnrm_126_2030_stand[[9]] * PC4_load[9,] 
PC4_bio10 <- Cnrm_126_2030_stand[[10]] * PC4_load[10,] 
PC4_bio11 <- Cnrm_126_2030_stand[[11]] * PC4_load[11,] 
PC4_bio12 <- Cnrm_126_2030_stand[[12]] * PC4_load[12,] 
PC4_bio13 <- Cnrm_126_2030_stand[[13]] * PC4_load[13,] 
PC4_bio14 <- Cnrm_126_2030_stand[[14]] * PC4_load[14,] 
PC4_bio15 <- Cnrm_126_2030_stand[[15]] * PC4_load[15,] 
PC4_bio16 <- Cnrm_126_2030_stand[[16]] * PC4_load[16,] 
PC4_bio17 <- Cnrm_126_2030_stand[[17]] * PC4_load[17,] 
PC4_bio18 <- Cnrm_126_2030_stand[[18]] * PC4_load[18,] 
PC4_bio19 <- Cnrm_126_2030_stand[[19]] * PC4_load[19,] 

PC4_bio_Cnrm_126_2030 <- PC4_bio1 + PC4_bio2 + PC4_bio3 + PC4_bio4 + PC4_bio5 + PC4_bio6 + PC4_bio7 + PC4_bio8 + PC4_bio9 + PC4_bio10 + PC4_bio11 + PC4_bio12 + PC4_bio13 + PC4_bio14 + PC4_bio15 + PC4_bio16 + PC4_bio17 + PC4_bio18 + PC4_bio19

# PC5
PC5_bio1 <- Cnrm_126_2030_stand[[1]] * PC5_load[1,] 
PC5_bio2 <- Cnrm_126_2030_stand[[2]] * PC5_load[2,] 
PC5_bio3 <- Cnrm_126_2030_stand[[3]] * PC5_load[3,] 
PC5_bio4 <- Cnrm_126_2030_stand[[4]] * PC5_load[4,] 
PC5_bio5 <- Cnrm_126_2030_stand[[5]] * PC5_load[5,] 
PC5_bio6 <- Cnrm_126_2030_stand[[6]] * PC5_load[6,] 
PC5_bio7 <- Cnrm_126_2030_stand[[7]] * PC5_load[7,] 
PC5_bio8 <- Cnrm_126_2030_stand[[8]] * PC5_load[8,] 
PC5_bio9 <- Cnrm_126_2030_stand[[9]] * PC5_load[9,] 
PC5_bio10 <- Cnrm_126_2030_stand[[10]] * PC5_load[10,] 
PC5_bio11 <- Cnrm_126_2030_stand[[11]] * PC5_load[11,] 
PC5_bio12 <- Cnrm_126_2030_stand[[12]] * PC5_load[12,] 
PC5_bio13 <- Cnrm_126_2030_stand[[13]] * PC5_load[13,] 
PC5_bio14 <- Cnrm_126_2030_stand[[14]] * PC5_load[14,] 
PC5_bio15 <- Cnrm_126_2030_stand[[15]] * PC5_load[15,] 
PC5_bio16 <- Cnrm_126_2030_stand[[16]] * PC5_load[16,] 
PC5_bio17 <- Cnrm_126_2030_stand[[17]] * PC5_load[17,] 
PC5_bio18 <- Cnrm_126_2030_stand[[18]] * PC5_load[18,] 
PC5_bio19 <- Cnrm_126_2030_stand[[19]] * PC5_load[19,] 

PC5_bio_Cnrm_126_2030 <- PC5_bio1 + PC5_bio2 + PC5_bio3 + PC5_bio4 + PC5_bio5 + PC5_bio6 + PC5_bio7 + PC5_bio8 + PC5_bio9 + PC5_bio10 + PC5_bio11 + PC5_bio12 + PC5_bio13 + PC5_bio14 + PC5_bio15 + PC5_bio16 + PC5_bio17 + PC5_bio18 + PC5_bio19

dir.create("PCA_out_bio_Cnrm_126_2030")
setwd("PCA_out_bio_Cnrm_126_2030")
getwd()

# Write the rasters for the 5 PCs
writeRaster(PC1_bio_Cnrm_126_2030, "__PC1.tif")
writeRaster(PC2_bio_Cnrm_126_2030, "__PC2.tif")
writeRaster(PC3_bio_Cnrm_126_2030, "__PC3.tif")
writeRaster(PC4_bio_Cnrm_126_2030, "__PC4.tif")
writeRaster(PC5_bio_Cnrm_126_2030, "__PC5.tif")

setwd(paste0("C:/Users/Public/Documents/Modelling/", species))

########################################## CNRM-CM6-1 ##########################################
########################################## 2030 ############################################
########################################## 370 ##############################################
bioclim_Cnrm_370_2030 <- stack("C:\\Users\\Public\\Documents\\Modelling\\wc2.1_2.5m_bioc_2021-2040\\wc2.1_2.5m_bioc_CNRM-CM6-1_ssp370_2021-2040.tif")

# Crop and mask them to species-specific extent
CropClim_Cnrm_370_2030=crop(bioclim_Cnrm_370_2030, CurrentProj) #Crop before mask to speed up processing time
MaskClim_Cnrm_370_2030=mask(CropClim_Cnrm_370_2030, CurrentProj)

# Then standardise them based on the current mean and SD values
val1 <- MaskClim_Cnrm_370_2030[[1]] - current.stats[1,1] 
Cnrm_370_2030_sbio1 <- val1 / current.stats[2,1]

val2 <- MaskClim_Cnrm_370_2030[[2]] - current.stats[1,2] 
Cnrm_370_2030_sbio2 <- val2 / current.stats[2,2]

val3 <- MaskClim_Cnrm_370_2030[[3]] - current.stats[1,3] 
Cnrm_370_2030_sbio3 <- val3 / current.stats[2,3]

val4 <- MaskClim_Cnrm_370_2030[[4]] - current.stats[1,4] 
Cnrm_370_2030_sbio4 <- val4 / current.stats[2,4]

val5 <- MaskClim_Cnrm_370_2030[[5]] - current.stats[1,5] 
Cnrm_370_2030_sbio5 <- val5 / current.stats[2,5]

val6 <- MaskClim_Cnrm_370_2030[[6]] - current.stats[1,6] 
Cnrm_370_2030_sbio6 <- val6 / current.stats[2,6]

val7 <- MaskClim_Cnrm_370_2030[[7]] - current.stats[1,7] 
Cnrm_370_2030_sbio7 <- val7 / current.stats[2,7]

val8 <- MaskClim_Cnrm_370_2030[[8]] - current.stats[1,8] 
Cnrm_370_2030_sbio8 <- val8 / current.stats[2,8]

val9 <- MaskClim_Cnrm_370_2030[[9]] - current.stats[1,9] 
Cnrm_370_2030_sbio9 <- val9 / current.stats[2,9]

val10 <- MaskClim_Cnrm_370_2030[[10]] - current.stats[1,10] 
Cnrm_370_2030_sbio10 <- val10 / current.stats[2,10]

val11 <- MaskClim_Cnrm_370_2030[[11]] - current.stats[1,11] 
Cnrm_370_2030_sbio11 <- val11 / current.stats[2,11]

val12 <- MaskClim_Cnrm_370_2030[[12]] - current.stats[1,12] 
Cnrm_370_2030_sbio12 <- val12 / current.stats[2,12]

val13 <- MaskClim_Cnrm_370_2030[[13]] - current.stats[1,13] 
Cnrm_370_2030_sbio13 <- val13 / current.stats[2,13]

val14 <- MaskClim_Cnrm_370_2030[[14]] - current.stats[1,14] 
Cnrm_370_2030_sbio14 <- val14 / current.stats[2,14]

val15 <- MaskClim_Cnrm_370_2030[[15]] - current.stats[1,15] 
Cnrm_370_2030_sbio15 <- val15 / current.stats[2,15]

val16 <- MaskClim_Cnrm_370_2030[[16]] - current.stats[1,16] 
Cnrm_370_2030_sbio16 <- val16 / current.stats[2,16]

val17 <- MaskClim_Cnrm_370_2030[[17]] - current.stats[1,17] 
Cnrm_370_2030_sbio17 <- val17 / current.stats[2,17]

val18 <- MaskClim_Cnrm_370_2030[[18]] - current.stats[1,18] 
Cnrm_370_2030_sbio18 <- val18 / current.stats[2,18]

val19 <- MaskClim_Cnrm_370_2030[[19]] - current.stats[1,19] 
Cnrm_370_2030_sbio19 <- val19 / current.stats[2,19]

# Then multiply PC1 loadings for each variable by the standardised future rasters
Cnrm_370_2030_stand <- stack(Cnrm_370_2030_sbio1, Cnrm_370_2030_sbio2, Cnrm_370_2030_sbio3, Cnrm_370_2030_sbio4, Cnrm_370_2030_sbio5, Cnrm_370_2030_sbio6,
                             Cnrm_370_2030_sbio7, Cnrm_370_2030_sbio8, Cnrm_370_2030_sbio9, Cnrm_370_2030_sbio10, Cnrm_370_2030_sbio11, Cnrm_370_2030_sbio12,
                             Cnrm_370_2030_sbio13, Cnrm_370_2030_sbio14, Cnrm_370_2030_sbio15, Cnrm_370_2030_sbio16, Cnrm_370_2030_sbio17, Cnrm_370_2030_sbio18,
                             Cnrm_370_2030_sbio19)

# PC1
PC1_bio1 <- Cnrm_370_2030_stand[[1]] * PC1_load[1,] 
PC1_bio2 <- Cnrm_370_2030_stand[[2]] * PC1_load[2,] 
PC1_bio3 <- Cnrm_370_2030_stand[[3]] * PC1_load[3,] 
PC1_bio4 <- Cnrm_370_2030_stand[[4]] * PC1_load[4,] 
PC1_bio5 <- Cnrm_370_2030_stand[[5]] * PC1_load[5,] 
PC1_bio6 <- Cnrm_370_2030_stand[[6]] * PC1_load[6,] 
PC1_bio7 <- Cnrm_370_2030_stand[[7]] * PC1_load[7,] 
PC1_bio8 <- Cnrm_370_2030_stand[[8]] * PC1_load[8,] 
PC1_bio9 <- Cnrm_370_2030_stand[[9]] * PC1_load[9,] 
PC1_bio10 <- Cnrm_370_2030_stand[[10]] * PC1_load[10,] 
PC1_bio11 <- Cnrm_370_2030_stand[[11]] * PC1_load[11,] 
PC1_bio12 <- Cnrm_370_2030_stand[[12]] * PC1_load[12,] 
PC1_bio13 <- Cnrm_370_2030_stand[[13]] * PC1_load[13,] 
PC1_bio14 <- Cnrm_370_2030_stand[[14]] * PC1_load[14,] 
PC1_bio15 <- Cnrm_370_2030_stand[[15]] * PC1_load[15,] 
PC1_bio16 <- Cnrm_370_2030_stand[[16]] * PC1_load[16,] 
PC1_bio17 <- Cnrm_370_2030_stand[[17]] * PC1_load[17,] 
PC1_bio18 <- Cnrm_370_2030_stand[[18]] * PC1_load[18,] 
PC1_bio19 <- Cnrm_370_2030_stand[[19]] * PC1_load[19,] 

PC1_bio_Cnrm_370_2030 <- PC1_bio1 + PC1_bio2 + PC1_bio3 + PC1_bio4 + PC1_bio5 + PC1_bio6 + PC1_bio7 + PC1_bio8 + PC1_bio9 + PC1_bio10 + PC1_bio11 + PC1_bio12 + PC1_bio13 + PC1_bio14 + PC1_bio15 + PC1_bio16 + PC1_bio17 + PC1_bio18 + PC1_bio19

# PC2
PC2_bio1 <- Cnrm_370_2030_stand[[1]] * PC2_load[1,] 
PC2_bio2 <- Cnrm_370_2030_stand[[2]] * PC2_load[2,] 
PC2_bio3 <- Cnrm_370_2030_stand[[3]] * PC2_load[3,] 
PC2_bio4 <- Cnrm_370_2030_stand[[4]] * PC2_load[4,] 
PC2_bio5 <- Cnrm_370_2030_stand[[5]] * PC2_load[5,] 
PC2_bio6 <- Cnrm_370_2030_stand[[6]] * PC2_load[6,] 
PC2_bio7 <- Cnrm_370_2030_stand[[7]] * PC2_load[7,] 
PC2_bio8 <- Cnrm_370_2030_stand[[8]] * PC2_load[8,] 
PC2_bio9 <- Cnrm_370_2030_stand[[9]] * PC2_load[9,] 
PC2_bio10 <- Cnrm_370_2030_stand[[10]] * PC2_load[10,] 
PC2_bio11 <- Cnrm_370_2030_stand[[11]] * PC2_load[11,] 
PC2_bio12 <- Cnrm_370_2030_stand[[12]] * PC2_load[12,] 
PC2_bio13 <- Cnrm_370_2030_stand[[13]] * PC2_load[13,] 
PC2_bio14 <- Cnrm_370_2030_stand[[14]] * PC2_load[14,] 
PC2_bio15 <- Cnrm_370_2030_stand[[15]] * PC2_load[15,] 
PC2_bio16 <- Cnrm_370_2030_stand[[16]] * PC2_load[16,] 
PC2_bio17 <- Cnrm_370_2030_stand[[17]] * PC2_load[17,] 
PC2_bio18 <- Cnrm_370_2030_stand[[18]] * PC2_load[18,] 
PC2_bio19 <- Cnrm_370_2030_stand[[19]] * PC2_load[19,] 

PC2_bio_Cnrm_370_2030 <- PC2_bio1 + PC2_bio2 + PC2_bio3 + PC2_bio4 + PC2_bio5 + PC2_bio6 + PC2_bio7 + PC2_bio8 + PC2_bio9 + PC2_bio10 + PC2_bio11 + PC2_bio12 + PC2_bio13 + PC2_bio14 + PC2_bio15 + PC2_bio16 + PC2_bio17 + PC2_bio18 + PC2_bio19

# PC3
PC3_bio1 <- Cnrm_370_2030_stand[[1]] * PC3_load[1,] 
PC3_bio2 <- Cnrm_370_2030_stand[[2]] * PC3_load[2,] 
PC3_bio3 <- Cnrm_370_2030_stand[[3]] * PC3_load[3,] 
PC3_bio4 <- Cnrm_370_2030_stand[[4]] * PC3_load[4,] 
PC3_bio5 <- Cnrm_370_2030_stand[[5]] * PC3_load[5,] 
PC3_bio6 <- Cnrm_370_2030_stand[[6]] * PC3_load[6,] 
PC3_bio7 <- Cnrm_370_2030_stand[[7]] * PC3_load[7,] 
PC3_bio8 <- Cnrm_370_2030_stand[[8]] * PC3_load[8,] 
PC3_bio9 <- Cnrm_370_2030_stand[[9]] * PC3_load[9,] 
PC3_bio10 <- Cnrm_370_2030_stand[[10]] * PC3_load[10,] 
PC3_bio11 <- Cnrm_370_2030_stand[[11]] * PC3_load[11,] 
PC3_bio12 <- Cnrm_370_2030_stand[[12]] * PC3_load[12,] 
PC3_bio13 <- Cnrm_370_2030_stand[[13]] * PC3_load[13,] 
PC3_bio14 <- Cnrm_370_2030_stand[[14]] * PC3_load[14,] 
PC3_bio15 <- Cnrm_370_2030_stand[[15]] * PC3_load[15,] 
PC3_bio16 <- Cnrm_370_2030_stand[[16]] * PC3_load[16,] 
PC3_bio17 <- Cnrm_370_2030_stand[[17]] * PC3_load[17,] 
PC3_bio18 <- Cnrm_370_2030_stand[[18]] * PC3_load[18,] 
PC3_bio19 <- Cnrm_370_2030_stand[[19]] * PC3_load[19,] 

PC3_bio_Cnrm_370_2030 <- PC3_bio1 + PC3_bio2 + PC3_bio3 + PC3_bio4 + PC3_bio5 + PC3_bio6 + PC3_bio7 + PC3_bio8 + PC3_bio9 + PC3_bio10 + PC3_bio11 + PC3_bio12 + PC3_bio13 + PC3_bio14 + PC3_bio15 + PC3_bio16 + PC3_bio17 + PC3_bio18 + PC3_bio19

# PC4
PC4_bio1 <- Cnrm_370_2030_stand[[1]] * PC4_load[1,] 
PC4_bio2 <- Cnrm_370_2030_stand[[2]] * PC4_load[2,] 
PC4_bio3 <- Cnrm_370_2030_stand[[3]] * PC4_load[3,] 
PC4_bio4 <- Cnrm_370_2030_stand[[4]] * PC4_load[4,] 
PC4_bio5 <- Cnrm_370_2030_stand[[5]] * PC4_load[5,] 
PC4_bio6 <- Cnrm_370_2030_stand[[6]] * PC4_load[6,] 
PC4_bio7 <- Cnrm_370_2030_stand[[7]] * PC4_load[7,] 
PC4_bio8 <- Cnrm_370_2030_stand[[8]] * PC4_load[8,] 
PC4_bio9 <- Cnrm_370_2030_stand[[9]] * PC4_load[9,] 
PC4_bio10 <- Cnrm_370_2030_stand[[10]] * PC4_load[10,] 
PC4_bio11 <- Cnrm_370_2030_stand[[11]] * PC4_load[11,] 
PC4_bio12 <- Cnrm_370_2030_stand[[12]] * PC4_load[12,] 
PC4_bio13 <- Cnrm_370_2030_stand[[13]] * PC4_load[13,] 
PC4_bio14 <- Cnrm_370_2030_stand[[14]] * PC4_load[14,] 
PC4_bio15 <- Cnrm_370_2030_stand[[15]] * PC4_load[15,] 
PC4_bio16 <- Cnrm_370_2030_stand[[16]] * PC4_load[16,] 
PC4_bio17 <- Cnrm_370_2030_stand[[17]] * PC4_load[17,] 
PC4_bio18 <- Cnrm_370_2030_stand[[18]] * PC4_load[18,] 
PC4_bio19 <- Cnrm_370_2030_stand[[19]] * PC4_load[19,] 

PC4_bio_Cnrm_370_2030 <- PC4_bio1 + PC4_bio2 + PC4_bio3 + PC4_bio4 + PC4_bio5 + PC4_bio6 + PC4_bio7 + PC4_bio8 + PC4_bio9 + PC4_bio10 + PC4_bio11 + PC4_bio12 + PC4_bio13 + PC4_bio14 + PC4_bio15 + PC4_bio16 + PC4_bio17 + PC4_bio18 + PC4_bio19

# PC5
PC5_bio1 <- Cnrm_370_2030_stand[[1]] * PC5_load[1,] 
PC5_bio2 <- Cnrm_370_2030_stand[[2]] * PC5_load[2,] 
PC5_bio3 <- Cnrm_370_2030_stand[[3]] * PC5_load[3,] 
PC5_bio4 <- Cnrm_370_2030_stand[[4]] * PC5_load[4,] 
PC5_bio5 <- Cnrm_370_2030_stand[[5]] * PC5_load[5,] 
PC5_bio6 <- Cnrm_370_2030_stand[[6]] * PC5_load[6,] 
PC5_bio7 <- Cnrm_370_2030_stand[[7]] * PC5_load[7,] 
PC5_bio8 <- Cnrm_370_2030_stand[[8]] * PC5_load[8,] 
PC5_bio9 <- Cnrm_370_2030_stand[[9]] * PC5_load[9,] 
PC5_bio10 <- Cnrm_370_2030_stand[[10]] * PC5_load[10,] 
PC5_bio11 <- Cnrm_370_2030_stand[[11]] * PC5_load[11,] 
PC5_bio12 <- Cnrm_370_2030_stand[[12]] * PC5_load[12,] 
PC5_bio13 <- Cnrm_370_2030_stand[[13]] * PC5_load[13,] 
PC5_bio14 <- Cnrm_370_2030_stand[[14]] * PC5_load[14,] 
PC5_bio15 <- Cnrm_370_2030_stand[[15]] * PC5_load[15,] 
PC5_bio16 <- Cnrm_370_2030_stand[[16]] * PC5_load[16,] 
PC5_bio17 <- Cnrm_370_2030_stand[[17]] * PC5_load[17,] 
PC5_bio18 <- Cnrm_370_2030_stand[[18]] * PC5_load[18,] 
PC5_bio19 <- Cnrm_370_2030_stand[[19]] * PC5_load[19,] 

PC5_bio_Cnrm_370_2030 <- PC5_bio1 + PC5_bio2 + PC5_bio3 + PC5_bio4 + PC5_bio5 + PC5_bio6 + PC5_bio7 + PC5_bio8 + PC5_bio9 + PC5_bio10 + PC5_bio11 + PC5_bio12 + PC5_bio13 + PC5_bio14 + PC5_bio15 + PC5_bio16 + PC5_bio17 + PC5_bio18 + PC5_bio19

dir.create("PCA_out_bio_Cnrm_370_2030")
setwd("PCA_out_bio_Cnrm_370_2030")
getwd()

# Write the rasters for the 5 PCs
writeRaster(PC1_bio_Cnrm_370_2030, "__PC1.tif")
writeRaster(PC2_bio_Cnrm_370_2030, "__PC2.tif")
writeRaster(PC3_bio_Cnrm_370_2030, "__PC3.tif")
writeRaster(PC4_bio_Cnrm_370_2030, "__PC4.tif")
writeRaster(PC5_bio_Cnrm_370_2030, "__PC5.tif")

setwd(paste0("C:/Users/Public/Documents/Modelling/", species))

########################################## IPSL-CM6A-LR ##########################################
########################################## 2030 ############################################
########################################## 126 ##############################################
bioclim_Ipsl_126_2030 <- stack("C:\\Users\\Public\\Documents\\Modelling\\wc2.1_2.5m_bioc_2021-2040\\wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp126_2021-2040.tif")

# Crop and mask them to species-specific extent
CropClim_Ipsl_126_2030=crop(bioclim_Ipsl_126_2030, CurrentProj) #Crop before mask to speed up processing time
MaskClim_Ipsl_126_2030=mask(CropClim_Ipsl_126_2030, CurrentProj)

# Then standardise them based on the current mean and SD values
val1 <- MaskClim_Ipsl_126_2030[[1]] - current.stats[1,1] 
Ipsl_126_2030_sbio1 <- val1 / current.stats[2,1]

val2 <- MaskClim_Ipsl_126_2030[[2]] - current.stats[1,2] 
Ipsl_126_2030_sbio2 <- val2 / current.stats[2,2]

val3 <- MaskClim_Ipsl_126_2030[[3]] - current.stats[1,3] 
Ipsl_126_2030_sbio3 <- val3 / current.stats[2,3]

val4 <- MaskClim_Ipsl_126_2030[[4]] - current.stats[1,4] 
Ipsl_126_2030_sbio4 <- val4 / current.stats[2,4]

val5 <- MaskClim_Ipsl_126_2030[[5]] - current.stats[1,5] 
Ipsl_126_2030_sbio5 <- val5 / current.stats[2,5]

val6 <- MaskClim_Ipsl_126_2030[[6]] - current.stats[1,6] 
Ipsl_126_2030_sbio6 <- val6 / current.stats[2,6]

val7 <- MaskClim_Ipsl_126_2030[[7]] - current.stats[1,7] 
Ipsl_126_2030_sbio7 <- val7 / current.stats[2,7]

val8 <- MaskClim_Ipsl_126_2030[[8]] - current.stats[1,8] 
Ipsl_126_2030_sbio8 <- val8 / current.stats[2,8]

val9 <- MaskClim_Ipsl_126_2030[[9]] - current.stats[1,9] 
Ipsl_126_2030_sbio9 <- val9 / current.stats[2,9]

val10 <- MaskClim_Ipsl_126_2030[[10]] - current.stats[1,10] 
Ipsl_126_2030_sbio10 <- val10 / current.stats[2,10]

val11 <- MaskClim_Ipsl_126_2030[[11]] - current.stats[1,11] 
Ipsl_126_2030_sbio11 <- val11 / current.stats[2,11]

val12 <- MaskClim_Ipsl_126_2030[[12]] - current.stats[1,12] 
Ipsl_126_2030_sbio12 <- val12 / current.stats[2,12]

val13 <- MaskClim_Ipsl_126_2030[[13]] - current.stats[1,13] 
Ipsl_126_2030_sbio13 <- val13 / current.stats[2,13]

val14 <- MaskClim_Ipsl_126_2030[[14]] - current.stats[1,14] 
Ipsl_126_2030_sbio14 <- val14 / current.stats[2,14]

val15 <- MaskClim_Ipsl_126_2030[[15]] - current.stats[1,15] 
Ipsl_126_2030_sbio15 <- val15 / current.stats[2,15]

val16 <- MaskClim_Ipsl_126_2030[[16]] - current.stats[1,16] 
Ipsl_126_2030_sbio16 <- val16 / current.stats[2,16]

val17 <- MaskClim_Ipsl_126_2030[[17]] - current.stats[1,17] 
Ipsl_126_2030_sbio17 <- val17 / current.stats[2,17]

val18 <- MaskClim_Ipsl_126_2030[[18]] - current.stats[1,18] 
Ipsl_126_2030_sbio18 <- val18 / current.stats[2,18]

val19 <- MaskClim_Ipsl_126_2030[[19]] - current.stats[1,19] 
Ipsl_126_2030_sbio19 <- val19 / current.stats[2,19]

# Then multiply PC1 loadings for each variable by the standardised future rasters
Ipsl_126_2030_stand <- stack(Ipsl_126_2030_sbio1, Ipsl_126_2030_sbio2, Ipsl_126_2030_sbio3, Ipsl_126_2030_sbio4, Ipsl_126_2030_sbio5, Ipsl_126_2030_sbio6,
                             Ipsl_126_2030_sbio7, Ipsl_126_2030_sbio8, Ipsl_126_2030_sbio9, Ipsl_126_2030_sbio10, Ipsl_126_2030_sbio11, Ipsl_126_2030_sbio12,
                             Ipsl_126_2030_sbio13, Ipsl_126_2030_sbio14, Ipsl_126_2030_sbio15, Ipsl_126_2030_sbio16, Ipsl_126_2030_sbio17, Ipsl_126_2030_sbio18,
                             Ipsl_126_2030_sbio19)

# PC1
PC1_bio1 <- Ipsl_126_2030_stand[[1]] * PC1_load[1,] 
PC1_bio2 <- Ipsl_126_2030_stand[[2]] * PC1_load[2,] 
PC1_bio3 <- Ipsl_126_2030_stand[[3]] * PC1_load[3,] 
PC1_bio4 <- Ipsl_126_2030_stand[[4]] * PC1_load[4,] 
PC1_bio5 <- Ipsl_126_2030_stand[[5]] * PC1_load[5,] 
PC1_bio6 <- Ipsl_126_2030_stand[[6]] * PC1_load[6,] 
PC1_bio7 <- Ipsl_126_2030_stand[[7]] * PC1_load[7,] 
PC1_bio8 <- Ipsl_126_2030_stand[[8]] * PC1_load[8,] 
PC1_bio9 <- Ipsl_126_2030_stand[[9]] * PC1_load[9,] 
PC1_bio10 <- Ipsl_126_2030_stand[[10]] * PC1_load[10,] 
PC1_bio11 <- Ipsl_126_2030_stand[[11]] * PC1_load[11,] 
PC1_bio12 <- Ipsl_126_2030_stand[[12]] * PC1_load[12,] 
PC1_bio13 <- Ipsl_126_2030_stand[[13]] * PC1_load[13,] 
PC1_bio14 <- Ipsl_126_2030_stand[[14]] * PC1_load[14,] 
PC1_bio15 <- Ipsl_126_2030_stand[[15]] * PC1_load[15,] 
PC1_bio16 <- Ipsl_126_2030_stand[[16]] * PC1_load[16,] 
PC1_bio17 <- Ipsl_126_2030_stand[[17]] * PC1_load[17,] 
PC1_bio18 <- Ipsl_126_2030_stand[[18]] * PC1_load[18,] 
PC1_bio19 <- Ipsl_126_2030_stand[[19]] * PC1_load[19,] 

PC1_bio_Ipsl_126_2030 <- PC1_bio1 + PC1_bio2 + PC1_bio3 + PC1_bio4 + PC1_bio5 + PC1_bio6 + PC1_bio7 + PC1_bio8 + PC1_bio9 + PC1_bio10 + PC1_bio11 + PC1_bio12 + PC1_bio13 + PC1_bio14 + PC1_bio15 + PC1_bio16 + PC1_bio17 + PC1_bio18 + PC1_bio19

# PC2
PC2_bio1 <- Ipsl_126_2030_stand[[1]] * PC2_load[1,] 
PC2_bio2 <- Ipsl_126_2030_stand[[2]] * PC2_load[2,] 
PC2_bio3 <- Ipsl_126_2030_stand[[3]] * PC2_load[3,] 
PC2_bio4 <- Ipsl_126_2030_stand[[4]] * PC2_load[4,] 
PC2_bio5 <- Ipsl_126_2030_stand[[5]] * PC2_load[5,] 
PC2_bio6 <- Ipsl_126_2030_stand[[6]] * PC2_load[6,] 
PC2_bio7 <- Ipsl_126_2030_stand[[7]] * PC2_load[7,] 
PC2_bio8 <- Ipsl_126_2030_stand[[8]] * PC2_load[8,] 
PC2_bio9 <- Ipsl_126_2030_stand[[9]] * PC2_load[9,] 
PC2_bio10 <- Ipsl_126_2030_stand[[10]] * PC2_load[10,] 
PC2_bio11 <- Ipsl_126_2030_stand[[11]] * PC2_load[11,] 
PC2_bio12 <- Ipsl_126_2030_stand[[12]] * PC2_load[12,] 
PC2_bio13 <- Ipsl_126_2030_stand[[13]] * PC2_load[13,] 
PC2_bio14 <- Ipsl_126_2030_stand[[14]] * PC2_load[14,] 
PC2_bio15 <- Ipsl_126_2030_stand[[15]] * PC2_load[15,] 
PC2_bio16 <- Ipsl_126_2030_stand[[16]] * PC2_load[16,] 
PC2_bio17 <- Ipsl_126_2030_stand[[17]] * PC2_load[17,] 
PC2_bio18 <- Ipsl_126_2030_stand[[18]] * PC2_load[18,] 
PC2_bio19 <- Ipsl_126_2030_stand[[19]] * PC2_load[19,] 

PC2_bio_Ipsl_126_2030 <- PC2_bio1 + PC2_bio2 + PC2_bio3 + PC2_bio4 + PC2_bio5 + PC2_bio6 + PC2_bio7 + PC2_bio8 + PC2_bio9 + PC2_bio10 + PC2_bio11 + PC2_bio12 + PC2_bio13 + PC2_bio14 + PC2_bio15 + PC2_bio16 + PC2_bio17 + PC2_bio18 + PC2_bio19

# PC3
PC3_bio1 <- Ipsl_126_2030_stand[[1]] * PC3_load[1,] 
PC3_bio2 <- Ipsl_126_2030_stand[[2]] * PC3_load[2,] 
PC3_bio3 <- Ipsl_126_2030_stand[[3]] * PC3_load[3,] 
PC3_bio4 <- Ipsl_126_2030_stand[[4]] * PC3_load[4,] 
PC3_bio5 <- Ipsl_126_2030_stand[[5]] * PC3_load[5,] 
PC3_bio6 <- Ipsl_126_2030_stand[[6]] * PC3_load[6,] 
PC3_bio7 <- Ipsl_126_2030_stand[[7]] * PC3_load[7,] 
PC3_bio8 <- Ipsl_126_2030_stand[[8]] * PC3_load[8,] 
PC3_bio9 <- Ipsl_126_2030_stand[[9]] * PC3_load[9,] 
PC3_bio10 <- Ipsl_126_2030_stand[[10]] * PC3_load[10,] 
PC3_bio11 <- Ipsl_126_2030_stand[[11]] * PC3_load[11,] 
PC3_bio12 <- Ipsl_126_2030_stand[[12]] * PC3_load[12,] 
PC3_bio13 <- Ipsl_126_2030_stand[[13]] * PC3_load[13,] 
PC3_bio14 <- Ipsl_126_2030_stand[[14]] * PC3_load[14,] 
PC3_bio15 <- Ipsl_126_2030_stand[[15]] * PC3_load[15,] 
PC3_bio16 <- Ipsl_126_2030_stand[[16]] * PC3_load[16,] 
PC3_bio17 <- Ipsl_126_2030_stand[[17]] * PC3_load[17,] 
PC3_bio18 <- Ipsl_126_2030_stand[[18]] * PC3_load[18,] 
PC3_bio19 <- Ipsl_126_2030_stand[[19]] * PC3_load[19,] 

PC3_bio_Ipsl_126_2030 <- PC3_bio1 + PC3_bio2 + PC3_bio3 + PC3_bio4 + PC3_bio5 + PC3_bio6 + PC3_bio7 + PC3_bio8 + PC3_bio9 + PC3_bio10 + PC3_bio11 + PC3_bio12 + PC3_bio13 + PC3_bio14 + PC3_bio15 + PC3_bio16 + PC3_bio17 + PC3_bio18 + PC3_bio19

# PC4
PC4_bio1 <- Ipsl_126_2030_stand[[1]] * PC4_load[1,] 
PC4_bio2 <- Ipsl_126_2030_stand[[2]] * PC4_load[2,] 
PC4_bio3 <- Ipsl_126_2030_stand[[3]] * PC4_load[3,] 
PC4_bio4 <- Ipsl_126_2030_stand[[4]] * PC4_load[4,] 
PC4_bio5 <- Ipsl_126_2030_stand[[5]] * PC4_load[5,] 
PC4_bio6 <- Ipsl_126_2030_stand[[6]] * PC4_load[6,] 
PC4_bio7 <- Ipsl_126_2030_stand[[7]] * PC4_load[7,] 
PC4_bio8 <- Ipsl_126_2030_stand[[8]] * PC4_load[8,] 
PC4_bio9 <- Ipsl_126_2030_stand[[9]] * PC4_load[9,] 
PC4_bio10 <- Ipsl_126_2030_stand[[10]] * PC4_load[10,] 
PC4_bio11 <- Ipsl_126_2030_stand[[11]] * PC4_load[11,] 
PC4_bio12 <- Ipsl_126_2030_stand[[12]] * PC4_load[12,] 
PC4_bio13 <- Ipsl_126_2030_stand[[13]] * PC4_load[13,] 
PC4_bio14 <- Ipsl_126_2030_stand[[14]] * PC4_load[14,] 
PC4_bio15 <- Ipsl_126_2030_stand[[15]] * PC4_load[15,] 
PC4_bio16 <- Ipsl_126_2030_stand[[16]] * PC4_load[16,] 
PC4_bio17 <- Ipsl_126_2030_stand[[17]] * PC4_load[17,] 
PC4_bio18 <- Ipsl_126_2030_stand[[18]] * PC4_load[18,] 
PC4_bio19 <- Ipsl_126_2030_stand[[19]] * PC4_load[19,] 

PC4_bio_Ipsl_126_2030 <- PC4_bio1 + PC4_bio2 + PC4_bio3 + PC4_bio4 + PC4_bio5 + PC4_bio6 + PC4_bio7 + PC4_bio8 + PC4_bio9 + PC4_bio10 + PC4_bio11 + PC4_bio12 + PC4_bio13 + PC4_bio14 + PC4_bio15 + PC4_bio16 + PC4_bio17 + PC4_bio18 + PC4_bio19

# PC5
PC5_bio1 <- Ipsl_126_2030_stand[[1]] * PC5_load[1,] 
PC5_bio2 <- Ipsl_126_2030_stand[[2]] * PC5_load[2,] 
PC5_bio3 <- Ipsl_126_2030_stand[[3]] * PC5_load[3,] 
PC5_bio4 <- Ipsl_126_2030_stand[[4]] * PC5_load[4,] 
PC5_bio5 <- Ipsl_126_2030_stand[[5]] * PC5_load[5,] 
PC5_bio6 <- Ipsl_126_2030_stand[[6]] * PC5_load[6,] 
PC5_bio7 <- Ipsl_126_2030_stand[[7]] * PC5_load[7,] 
PC5_bio8 <- Ipsl_126_2030_stand[[8]] * PC5_load[8,] 
PC5_bio9 <- Ipsl_126_2030_stand[[9]] * PC5_load[9,] 
PC5_bio10 <- Ipsl_126_2030_stand[[10]] * PC5_load[10,] 
PC5_bio11 <- Ipsl_126_2030_stand[[11]] * PC5_load[11,] 
PC5_bio12 <- Ipsl_126_2030_stand[[12]] * PC5_load[12,] 
PC5_bio13 <- Ipsl_126_2030_stand[[13]] * PC5_load[13,] 
PC5_bio14 <- Ipsl_126_2030_stand[[14]] * PC5_load[14,] 
PC5_bio15 <- Ipsl_126_2030_stand[[15]] * PC5_load[15,] 
PC5_bio16 <- Ipsl_126_2030_stand[[16]] * PC5_load[16,] 
PC5_bio17 <- Ipsl_126_2030_stand[[17]] * PC5_load[17,] 
PC5_bio18 <- Ipsl_126_2030_stand[[18]] * PC5_load[18,] 
PC5_bio19 <- Ipsl_126_2030_stand[[19]] * PC5_load[19,] 

PC5_bio_Ipsl_126_2030 <- PC5_bio1 + PC5_bio2 + PC5_bio3 + PC5_bio4 + PC5_bio5 + PC5_bio6 + PC5_bio7 + PC5_bio8 + PC5_bio9 + PC5_bio10 + PC5_bio11 + PC5_bio12 + PC5_bio13 + PC5_bio14 + PC5_bio15 + PC5_bio16 + PC5_bio17 + PC5_bio18 + PC5_bio19

dir.create("PCA_out_bio_Ipsl_126_2030")
setwd("PCA_out_bio_Ipsl_126_2030")
getwd()

# Write the rasters for the 5 PCs
writeRaster(PC1_bio_Ipsl_126_2030, "__PC1.tif")
writeRaster(PC2_bio_Ipsl_126_2030, "__PC2.tif")
writeRaster(PC3_bio_Ipsl_126_2030, "__PC3.tif")
writeRaster(PC4_bio_Ipsl_126_2030, "__PC4.tif")
writeRaster(PC5_bio_Ipsl_126_2030, "__PC5.tif")

setwd(paste0("C:/Users/Public/Documents/Modelling/", species))

########################################## IPSL-CM6A-LR ##########################################
########################################## 2030 ############################################
########################################## 370 ##############################################
bioclim_Ipsl_370_2030 <- stack("C:\\Users\\Public\\Documents\\Modelling\\wc2.1_2.5m_bioc_2021-2040\\wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp370_2021-2040.tif")

# Crop and mask them to species-specific extent
CropClim_Ipsl_370_2030=crop(bioclim_Ipsl_370_2030, CurrentProj) #Crop before mask to speed up processing time
MaskClim_Ipsl_370_2030=mask(CropClim_Ipsl_370_2030, CurrentProj)

# Then standardise them based on the current mean and SD values
val1 <- MaskClim_Ipsl_370_2030[[1]] - current.stats[1,1] 
Ipsl_370_2030_sbio1 <- val1 / current.stats[2,1]

val2 <- MaskClim_Ipsl_370_2030[[2]] - current.stats[1,2] 
Ipsl_370_2030_sbio2 <- val2 / current.stats[2,2]

val3 <- MaskClim_Ipsl_370_2030[[3]] - current.stats[1,3] 
Ipsl_370_2030_sbio3 <- val3 / current.stats[2,3]

val4 <- MaskClim_Ipsl_370_2030[[4]] - current.stats[1,4] 
Ipsl_370_2030_sbio4 <- val4 / current.stats[2,4]

val5 <- MaskClim_Ipsl_370_2030[[5]] - current.stats[1,5] 
Ipsl_370_2030_sbio5 <- val5 / current.stats[2,5]

val6 <- MaskClim_Ipsl_370_2030[[6]] - current.stats[1,6] 
Ipsl_370_2030_sbio6 <- val6 / current.stats[2,6]

val7 <- MaskClim_Ipsl_370_2030[[7]] - current.stats[1,7] 
Ipsl_370_2030_sbio7 <- val7 / current.stats[2,7]

val8 <- MaskClim_Ipsl_370_2030[[8]] - current.stats[1,8] 
Ipsl_370_2030_sbio8 <- val8 / current.stats[2,8]

val9 <- MaskClim_Ipsl_370_2030[[9]] - current.stats[1,9] 
Ipsl_370_2030_sbio9 <- val9 / current.stats[2,9]

val10 <- MaskClim_Ipsl_370_2030[[10]] - current.stats[1,10] 
Ipsl_370_2030_sbio10 <- val10 / current.stats[2,10]

val11 <- MaskClim_Ipsl_370_2030[[11]] - current.stats[1,11] 
Ipsl_370_2030_sbio11 <- val11 / current.stats[2,11]

val12 <- MaskClim_Ipsl_370_2030[[12]] - current.stats[1,12] 
Ipsl_370_2030_sbio12 <- val12 / current.stats[2,12]

val13 <- MaskClim_Ipsl_370_2030[[13]] - current.stats[1,13] 
Ipsl_370_2030_sbio13 <- val13 / current.stats[2,13]

val14 <- MaskClim_Ipsl_370_2030[[14]] - current.stats[1,14] 
Ipsl_370_2030_sbio14 <- val14 / current.stats[2,14]

val15 <- MaskClim_Ipsl_370_2030[[15]] - current.stats[1,15] 
Ipsl_370_2030_sbio15 <- val15 / current.stats[2,15]

val16 <- MaskClim_Ipsl_370_2030[[16]] - current.stats[1,16] 
Ipsl_370_2030_sbio16 <- val16 / current.stats[2,16]

val17 <- MaskClim_Ipsl_370_2030[[17]] - current.stats[1,17] 
Ipsl_370_2030_sbio17 <- val17 / current.stats[2,17]

val18 <- MaskClim_Ipsl_370_2030[[18]] - current.stats[1,18] 
Ipsl_370_2030_sbio18 <- val18 / current.stats[2,18]

val19 <- MaskClim_Ipsl_370_2030[[19]] - current.stats[1,19] 
Ipsl_370_2030_sbio19 <- val19 / current.stats[2,19]

# Then multiply PC1 loadings for each variable by the standardised future rasters
Ipsl_370_2030_stand <- stack(Ipsl_370_2030_sbio1, Ipsl_370_2030_sbio2, Ipsl_370_2030_sbio3, Ipsl_370_2030_sbio4, Ipsl_370_2030_sbio5, Ipsl_370_2030_sbio6,
                             Ipsl_370_2030_sbio7, Ipsl_370_2030_sbio8, Ipsl_370_2030_sbio9, Ipsl_370_2030_sbio10, Ipsl_370_2030_sbio11, Ipsl_370_2030_sbio12,
                             Ipsl_370_2030_sbio13, Ipsl_370_2030_sbio14, Ipsl_370_2030_sbio15, Ipsl_370_2030_sbio16, Ipsl_370_2030_sbio17, Ipsl_370_2030_sbio18,
                             Ipsl_370_2030_sbio19)

# PC1
PC1_bio1 <- Ipsl_370_2030_stand[[1]] * PC1_load[1,] 
PC1_bio2 <- Ipsl_370_2030_stand[[2]] * PC1_load[2,] 
PC1_bio3 <- Ipsl_370_2030_stand[[3]] * PC1_load[3,] 
PC1_bio4 <- Ipsl_370_2030_stand[[4]] * PC1_load[4,] 
PC1_bio5 <- Ipsl_370_2030_stand[[5]] * PC1_load[5,] 
PC1_bio6 <- Ipsl_370_2030_stand[[6]] * PC1_load[6,] 
PC1_bio7 <- Ipsl_370_2030_stand[[7]] * PC1_load[7,] 
PC1_bio8 <- Ipsl_370_2030_stand[[8]] * PC1_load[8,] 
PC1_bio9 <- Ipsl_370_2030_stand[[9]] * PC1_load[9,] 
PC1_bio10 <- Ipsl_370_2030_stand[[10]] * PC1_load[10,] 
PC1_bio11 <- Ipsl_370_2030_stand[[11]] * PC1_load[11,] 
PC1_bio12 <- Ipsl_370_2030_stand[[12]] * PC1_load[12,] 
PC1_bio13 <- Ipsl_370_2030_stand[[13]] * PC1_load[13,] 
PC1_bio14 <- Ipsl_370_2030_stand[[14]] * PC1_load[14,] 
PC1_bio15 <- Ipsl_370_2030_stand[[15]] * PC1_load[15,] 
PC1_bio16 <- Ipsl_370_2030_stand[[16]] * PC1_load[16,] 
PC1_bio17 <- Ipsl_370_2030_stand[[17]] * PC1_load[17,] 
PC1_bio18 <- Ipsl_370_2030_stand[[18]] * PC1_load[18,] 
PC1_bio19 <- Ipsl_370_2030_stand[[19]] * PC1_load[19,] 

PC1_bio_Ipsl_370_2030 <- PC1_bio1 + PC1_bio2 + PC1_bio3 + PC1_bio4 + PC1_bio5 + PC1_bio6 + PC1_bio7 + PC1_bio8 + PC1_bio9 + PC1_bio10 + PC1_bio11 + PC1_bio12 + PC1_bio13 + PC1_bio14 + PC1_bio15 + PC1_bio16 + PC1_bio17 + PC1_bio18 + PC1_bio19

# PC2
PC2_bio1 <- Ipsl_370_2030_stand[[1]] * PC2_load[1,] 
PC2_bio2 <- Ipsl_370_2030_stand[[2]] * PC2_load[2,] 
PC2_bio3 <- Ipsl_370_2030_stand[[3]] * PC2_load[3,] 
PC2_bio4 <- Ipsl_370_2030_stand[[4]] * PC2_load[4,] 
PC2_bio5 <- Ipsl_370_2030_stand[[5]] * PC2_load[5,] 
PC2_bio6 <- Ipsl_370_2030_stand[[6]] * PC2_load[6,] 
PC2_bio7 <- Ipsl_370_2030_stand[[7]] * PC2_load[7,] 
PC2_bio8 <- Ipsl_370_2030_stand[[8]] * PC2_load[8,] 
PC2_bio9 <- Ipsl_370_2030_stand[[9]] * PC2_load[9,] 
PC2_bio10 <- Ipsl_370_2030_stand[[10]] * PC2_load[10,] 
PC2_bio11 <- Ipsl_370_2030_stand[[11]] * PC2_load[11,] 
PC2_bio12 <- Ipsl_370_2030_stand[[12]] * PC2_load[12,] 
PC2_bio13 <- Ipsl_370_2030_stand[[13]] * PC2_load[13,] 
PC2_bio14 <- Ipsl_370_2030_stand[[14]] * PC2_load[14,] 
PC2_bio15 <- Ipsl_370_2030_stand[[15]] * PC2_load[15,] 
PC2_bio16 <- Ipsl_370_2030_stand[[16]] * PC2_load[16,] 
PC2_bio17 <- Ipsl_370_2030_stand[[17]] * PC2_load[17,] 
PC2_bio18 <- Ipsl_370_2030_stand[[18]] * PC2_load[18,] 
PC2_bio19 <- Ipsl_370_2030_stand[[19]] * PC2_load[19,] 

PC2_bio_Ipsl_370_2030 <- PC2_bio1 + PC2_bio2 + PC2_bio3 + PC2_bio4 + PC2_bio5 + PC2_bio6 + PC2_bio7 + PC2_bio8 + PC2_bio9 + PC2_bio10 + PC2_bio11 + PC2_bio12 + PC2_bio13 + PC2_bio14 + PC2_bio15 + PC2_bio16 + PC2_bio17 + PC2_bio18 + PC2_bio19

# PC3
PC3_bio1 <- Ipsl_370_2030_stand[[1]] * PC3_load[1,] 
PC3_bio2 <- Ipsl_370_2030_stand[[2]] * PC3_load[2,] 
PC3_bio3 <- Ipsl_370_2030_stand[[3]] * PC3_load[3,] 
PC3_bio4 <- Ipsl_370_2030_stand[[4]] * PC3_load[4,] 
PC3_bio5 <- Ipsl_370_2030_stand[[5]] * PC3_load[5,] 
PC3_bio6 <- Ipsl_370_2030_stand[[6]] * PC3_load[6,] 
PC3_bio7 <- Ipsl_370_2030_stand[[7]] * PC3_load[7,] 
PC3_bio8 <- Ipsl_370_2030_stand[[8]] * PC3_load[8,] 
PC3_bio9 <- Ipsl_370_2030_stand[[9]] * PC3_load[9,] 
PC3_bio10 <- Ipsl_370_2030_stand[[10]] * PC3_load[10,] 
PC3_bio11 <- Ipsl_370_2030_stand[[11]] * PC3_load[11,] 
PC3_bio12 <- Ipsl_370_2030_stand[[12]] * PC3_load[12,] 
PC3_bio13 <- Ipsl_370_2030_stand[[13]] * PC3_load[13,] 
PC3_bio14 <- Ipsl_370_2030_stand[[14]] * PC3_load[14,] 
PC3_bio15 <- Ipsl_370_2030_stand[[15]] * PC3_load[15,] 
PC3_bio16 <- Ipsl_370_2030_stand[[16]] * PC3_load[16,] 
PC3_bio17 <- Ipsl_370_2030_stand[[17]] * PC3_load[17,] 
PC3_bio18 <- Ipsl_370_2030_stand[[18]] * PC3_load[18,] 
PC3_bio19 <- Ipsl_370_2030_stand[[19]] * PC3_load[19,] 

PC3_bio_Ipsl_370_2030 <- PC3_bio1 + PC3_bio2 + PC3_bio3 + PC3_bio4 + PC3_bio5 + PC3_bio6 + PC3_bio7 + PC3_bio8 + PC3_bio9 + PC3_bio10 + PC3_bio11 + PC3_bio12 + PC3_bio13 + PC3_bio14 + PC3_bio15 + PC3_bio16 + PC3_bio17 + PC3_bio18 + PC3_bio19

# PC4
PC4_bio1 <- Ipsl_370_2030_stand[[1]] * PC4_load[1,] 
PC4_bio2 <- Ipsl_370_2030_stand[[2]] * PC4_load[2,] 
PC4_bio3 <- Ipsl_370_2030_stand[[3]] * PC4_load[3,] 
PC4_bio4 <- Ipsl_370_2030_stand[[4]] * PC4_load[4,] 
PC4_bio5 <- Ipsl_370_2030_stand[[5]] * PC4_load[5,] 
PC4_bio6 <- Ipsl_370_2030_stand[[6]] * PC4_load[6,] 
PC4_bio7 <- Ipsl_370_2030_stand[[7]] * PC4_load[7,] 
PC4_bio8 <- Ipsl_370_2030_stand[[8]] * PC4_load[8,] 
PC4_bio9 <- Ipsl_370_2030_stand[[9]] * PC4_load[9,] 
PC4_bio10 <- Ipsl_370_2030_stand[[10]] * PC4_load[10,] 
PC4_bio11 <- Ipsl_370_2030_stand[[11]] * PC4_load[11,] 
PC4_bio12 <- Ipsl_370_2030_stand[[12]] * PC4_load[12,] 
PC4_bio13 <- Ipsl_370_2030_stand[[13]] * PC4_load[13,] 
PC4_bio14 <- Ipsl_370_2030_stand[[14]] * PC4_load[14,] 
PC4_bio15 <- Ipsl_370_2030_stand[[15]] * PC4_load[15,] 
PC4_bio16 <- Ipsl_370_2030_stand[[16]] * PC4_load[16,] 
PC4_bio17 <- Ipsl_370_2030_stand[[17]] * PC4_load[17,] 
PC4_bio18 <- Ipsl_370_2030_stand[[18]] * PC4_load[18,] 
PC4_bio19 <- Ipsl_370_2030_stand[[19]] * PC4_load[19,] 

PC4_bio_Ipsl_370_2030 <- PC4_bio1 + PC4_bio2 + PC4_bio3 + PC4_bio4 + PC4_bio5 + PC4_bio6 + PC4_bio7 + PC4_bio8 + PC4_bio9 + PC4_bio10 + PC4_bio11 + PC4_bio12 + PC4_bio13 + PC4_bio14 + PC4_bio15 + PC4_bio16 + PC4_bio17 + PC4_bio18 + PC4_bio19

# PC5
PC5_bio1 <- Ipsl_370_2030_stand[[1]] * PC5_load[1,] 
PC5_bio2 <- Ipsl_370_2030_stand[[2]] * PC5_load[2,] 
PC5_bio3 <- Ipsl_370_2030_stand[[3]] * PC5_load[3,] 
PC5_bio4 <- Ipsl_370_2030_stand[[4]] * PC5_load[4,] 
PC5_bio5 <- Ipsl_370_2030_stand[[5]] * PC5_load[5,] 
PC5_bio6 <- Ipsl_370_2030_stand[[6]] * PC5_load[6,] 
PC5_bio7 <- Ipsl_370_2030_stand[[7]] * PC5_load[7,] 
PC5_bio8 <- Ipsl_370_2030_stand[[8]] * PC5_load[8,] 
PC5_bio9 <- Ipsl_370_2030_stand[[9]] * PC5_load[9,] 
PC5_bio10 <- Ipsl_370_2030_stand[[10]] * PC5_load[10,] 
PC5_bio11 <- Ipsl_370_2030_stand[[11]] * PC5_load[11,] 
PC5_bio12 <- Ipsl_370_2030_stand[[12]] * PC5_load[12,] 
PC5_bio13 <- Ipsl_370_2030_stand[[13]] * PC5_load[13,] 
PC5_bio14 <- Ipsl_370_2030_stand[[14]] * PC5_load[14,] 
PC5_bio15 <- Ipsl_370_2030_stand[[15]] * PC5_load[15,] 
PC5_bio16 <- Ipsl_370_2030_stand[[16]] * PC5_load[16,] 
PC5_bio17 <- Ipsl_370_2030_stand[[17]] * PC5_load[17,] 
PC5_bio18 <- Ipsl_370_2030_stand[[18]] * PC5_load[18,] 
PC5_bio19 <- Ipsl_370_2030_stand[[19]] * PC5_load[19,] 

PC5_bio_Ipsl_370_2030 <- PC5_bio1 + PC5_bio2 + PC5_bio3 + PC5_bio4 + PC5_bio5 + PC5_bio6 + PC5_bio7 + PC5_bio8 + PC5_bio9 + PC5_bio10 + PC5_bio11 + PC5_bio12 + PC5_bio13 + PC5_bio14 + PC5_bio15 + PC5_bio16 + PC5_bio17 + PC5_bio18 + PC5_bio19

dir.create("PCA_out_bio_Ipsl_370_2030")
setwd("PCA_out_bio_Ipsl_370_2030")
getwd()

# Write the rasters for the 5 PCs
writeRaster(PC1_bio_Ipsl_370_2030, "__PC1.tif")
writeRaster(PC2_bio_Ipsl_370_2030, "__PC2.tif")
writeRaster(PC3_bio_Ipsl_370_2030, "__PC3.tif")
writeRaster(PC4_bio_Ipsl_370_2030, "__PC4.tif")
writeRaster(PC5_bio_Ipsl_370_2030, "__PC5.tif")

setwd(paste0("C:/Users/Public/Documents/Modelling/", species))

########################################## MIROC6 ##########################################
########################################## 2030 ############################################
########################################## 126 ##############################################
bioclim_Miroc_126_2030 <- stack("C:\\Users\\Public\\Documents\\Modelling\\wc2.1_2.5m_bioc_2021-2040\\wc2.1_2.5m_bioc_MIROC6_ssp126_2021-2040.tif")

# Crop and mask them to species-specific extent
CropClim_Miroc_126_2030=crop(bioclim_Miroc_126_2030, CurrentProj) #Crop before mask to speed up processing time
MaskClim_Miroc_126_2030=mask(CropClim_Miroc_126_2030, CurrentProj)

# Then standardise them based on the current mean and SD values
val1 <- MaskClim_Miroc_126_2030[[1]] - current.stats[1,1] 
Miroc_126_2030_sbio1 <- val1 / current.stats[2,1]

val2 <- MaskClim_Miroc_126_2030[[2]] - current.stats[1,2] 
Miroc_126_2030_sbio2 <- val2 / current.stats[2,2]

val3 <- MaskClim_Miroc_126_2030[[3]] - current.stats[1,3] 
Miroc_126_2030_sbio3 <- val3 / current.stats[2,3]

val4 <- MaskClim_Miroc_126_2030[[4]] - current.stats[1,4] 
Miroc_126_2030_sbio4 <- val4 / current.stats[2,4]

val5 <- MaskClim_Miroc_126_2030[[5]] - current.stats[1,5] 
Miroc_126_2030_sbio5 <- val5 / current.stats[2,5]

val6 <- MaskClim_Miroc_126_2030[[6]] - current.stats[1,6] 
Miroc_126_2030_sbio6 <- val6 / current.stats[2,6]

val7 <- MaskClim_Miroc_126_2030[[7]] - current.stats[1,7] 
Miroc_126_2030_sbio7 <- val7 / current.stats[2,7]

val8 <- MaskClim_Miroc_126_2030[[8]] - current.stats[1,8] 
Miroc_126_2030_sbio8 <- val8 / current.stats[2,8]

val9 <- MaskClim_Miroc_126_2030[[9]] - current.stats[1,9] 
Miroc_126_2030_sbio9 <- val9 / current.stats[2,9]

val10 <- MaskClim_Miroc_126_2030[[10]] - current.stats[1,10] 
Miroc_126_2030_sbio10 <- val10 / current.stats[2,10]

val11 <- MaskClim_Miroc_126_2030[[11]] - current.stats[1,11] 
Miroc_126_2030_sbio11 <- val11 / current.stats[2,11]

val12 <- MaskClim_Miroc_126_2030[[12]] - current.stats[1,12] 
Miroc_126_2030_sbio12 <- val12 / current.stats[2,12]

val13 <- MaskClim_Miroc_126_2030[[13]] - current.stats[1,13] 
Miroc_126_2030_sbio13 <- val13 / current.stats[2,13]

val14 <- MaskClim_Miroc_126_2030[[14]] - current.stats[1,14] 
Miroc_126_2030_sbio14 <- val14 / current.stats[2,14]

val15 <- MaskClim_Miroc_126_2030[[15]] - current.stats[1,15] 
Miroc_126_2030_sbio15 <- val15 / current.stats[2,15]

val16 <- MaskClim_Miroc_126_2030[[16]] - current.stats[1,16] 
Miroc_126_2030_sbio16 <- val16 / current.stats[2,16]

val17 <- MaskClim_Miroc_126_2030[[17]] - current.stats[1,17] 
Miroc_126_2030_sbio17 <- val17 / current.stats[2,17]

val18 <- MaskClim_Miroc_126_2030[[18]] - current.stats[1,18] 
Miroc_126_2030_sbio18 <- val18 / current.stats[2,18]

val19 <- MaskClim_Miroc_126_2030[[19]] - current.stats[1,19] 
Miroc_126_2030_sbio19 <- val19 / current.stats[2,19]

# Then multiply PC1 loadings for each variable by the standardised future rasters
Miroc_126_2030_stand <- stack(Miroc_126_2030_sbio1, Miroc_126_2030_sbio2, Miroc_126_2030_sbio3, Miroc_126_2030_sbio4, Miroc_126_2030_sbio5, Miroc_126_2030_sbio6,
                              Miroc_126_2030_sbio7, Miroc_126_2030_sbio8, Miroc_126_2030_sbio9, Miroc_126_2030_sbio10, Miroc_126_2030_sbio11, Miroc_126_2030_sbio12,
                              Miroc_126_2030_sbio13, Miroc_126_2030_sbio14, Miroc_126_2030_sbio15, Miroc_126_2030_sbio16, Miroc_126_2030_sbio17, Miroc_126_2030_sbio18,
                              Miroc_126_2030_sbio19)

# PC1
PC1_bio1 <- Miroc_126_2030_stand[[1]] * PC1_load[1,] 
PC1_bio2 <- Miroc_126_2030_stand[[2]] * PC1_load[2,] 
PC1_bio3 <- Miroc_126_2030_stand[[3]] * PC1_load[3,] 
PC1_bio4 <- Miroc_126_2030_stand[[4]] * PC1_load[4,] 
PC1_bio5 <- Miroc_126_2030_stand[[5]] * PC1_load[5,] 
PC1_bio6 <- Miroc_126_2030_stand[[6]] * PC1_load[6,] 
PC1_bio7 <- Miroc_126_2030_stand[[7]] * PC1_load[7,] 
PC1_bio8 <- Miroc_126_2030_stand[[8]] * PC1_load[8,] 
PC1_bio9 <- Miroc_126_2030_stand[[9]] * PC1_load[9,] 
PC1_bio10 <- Miroc_126_2030_stand[[10]] * PC1_load[10,] 
PC1_bio11 <- Miroc_126_2030_stand[[11]] * PC1_load[11,] 
PC1_bio12 <- Miroc_126_2030_stand[[12]] * PC1_load[12,] 
PC1_bio13 <- Miroc_126_2030_stand[[13]] * PC1_load[13,] 
PC1_bio14 <- Miroc_126_2030_stand[[14]] * PC1_load[14,] 
PC1_bio15 <- Miroc_126_2030_stand[[15]] * PC1_load[15,] 
PC1_bio16 <- Miroc_126_2030_stand[[16]] * PC1_load[16,] 
PC1_bio17 <- Miroc_126_2030_stand[[17]] * PC1_load[17,] 
PC1_bio18 <- Miroc_126_2030_stand[[18]] * PC1_load[18,] 
PC1_bio19 <- Miroc_126_2030_stand[[19]] * PC1_load[19,] 

PC1_bio_Miroc_126_2030 <- PC1_bio1 + PC1_bio2 + PC1_bio3 + PC1_bio4 + PC1_bio5 + PC1_bio6 + PC1_bio7 + PC1_bio8 + PC1_bio9 + PC1_bio10 + PC1_bio11 + PC1_bio12 + PC1_bio13 + PC1_bio14 + PC1_bio15 + PC1_bio16 + PC1_bio17 + PC1_bio18 + PC1_bio19

# PC2
PC2_bio1 <- Miroc_126_2030_stand[[1]] * PC2_load[1,] 
PC2_bio2 <- Miroc_126_2030_stand[[2]] * PC2_load[2,] 
PC2_bio3 <- Miroc_126_2030_stand[[3]] * PC2_load[3,] 
PC2_bio4 <- Miroc_126_2030_stand[[4]] * PC2_load[4,] 
PC2_bio5 <- Miroc_126_2030_stand[[5]] * PC2_load[5,] 
PC2_bio6 <- Miroc_126_2030_stand[[6]] * PC2_load[6,] 
PC2_bio7 <- Miroc_126_2030_stand[[7]] * PC2_load[7,] 
PC2_bio8 <- Miroc_126_2030_stand[[8]] * PC2_load[8,] 
PC2_bio9 <- Miroc_126_2030_stand[[9]] * PC2_load[9,] 
PC2_bio10 <- Miroc_126_2030_stand[[10]] * PC2_load[10,] 
PC2_bio11 <- Miroc_126_2030_stand[[11]] * PC2_load[11,] 
PC2_bio12 <- Miroc_126_2030_stand[[12]] * PC2_load[12,] 
PC2_bio13 <- Miroc_126_2030_stand[[13]] * PC2_load[13,] 
PC2_bio14 <- Miroc_126_2030_stand[[14]] * PC2_load[14,] 
PC2_bio15 <- Miroc_126_2030_stand[[15]] * PC2_load[15,] 
PC2_bio16 <- Miroc_126_2030_stand[[16]] * PC2_load[16,] 
PC2_bio17 <- Miroc_126_2030_stand[[17]] * PC2_load[17,] 
PC2_bio18 <- Miroc_126_2030_stand[[18]] * PC2_load[18,] 
PC2_bio19 <- Miroc_126_2030_stand[[19]] * PC2_load[19,] 

PC2_bio_Miroc_126_2030 <- PC2_bio1 + PC2_bio2 + PC2_bio3 + PC2_bio4 + PC2_bio5 + PC2_bio6 + PC2_bio7 + PC2_bio8 + PC2_bio9 + PC2_bio10 + PC2_bio11 + PC2_bio12 + PC2_bio13 + PC2_bio14 + PC2_bio15 + PC2_bio16 + PC2_bio17 + PC2_bio18 + PC2_bio19

# PC3
PC3_bio1 <- Miroc_126_2030_stand[[1]] * PC3_load[1,] 
PC3_bio2 <- Miroc_126_2030_stand[[2]] * PC3_load[2,] 
PC3_bio3 <- Miroc_126_2030_stand[[3]] * PC3_load[3,] 
PC3_bio4 <- Miroc_126_2030_stand[[4]] * PC3_load[4,] 
PC3_bio5 <- Miroc_126_2030_stand[[5]] * PC3_load[5,] 
PC3_bio6 <- Miroc_126_2030_stand[[6]] * PC3_load[6,] 
PC3_bio7 <- Miroc_126_2030_stand[[7]] * PC3_load[7,] 
PC3_bio8 <- Miroc_126_2030_stand[[8]] * PC3_load[8,] 
PC3_bio9 <- Miroc_126_2030_stand[[9]] * PC3_load[9,] 
PC3_bio10 <- Miroc_126_2030_stand[[10]] * PC3_load[10,] 
PC3_bio11 <- Miroc_126_2030_stand[[11]] * PC3_load[11,] 
PC3_bio12 <- Miroc_126_2030_stand[[12]] * PC3_load[12,] 
PC3_bio13 <- Miroc_126_2030_stand[[13]] * PC3_load[13,] 
PC3_bio14 <- Miroc_126_2030_stand[[14]] * PC3_load[14,] 
PC3_bio15 <- Miroc_126_2030_stand[[15]] * PC3_load[15,] 
PC3_bio16 <- Miroc_126_2030_stand[[16]] * PC3_load[16,] 
PC3_bio17 <- Miroc_126_2030_stand[[17]] * PC3_load[17,] 
PC3_bio18 <- Miroc_126_2030_stand[[18]] * PC3_load[18,] 
PC3_bio19 <- Miroc_126_2030_stand[[19]] * PC3_load[19,] 

PC3_bio_Miroc_126_2030 <- PC3_bio1 + PC3_bio2 + PC3_bio3 + PC3_bio4 + PC3_bio5 + PC3_bio6 + PC3_bio7 + PC3_bio8 + PC3_bio9 + PC3_bio10 + PC3_bio11 + PC3_bio12 + PC3_bio13 + PC3_bio14 + PC3_bio15 + PC3_bio16 + PC3_bio17 + PC3_bio18 + PC3_bio19

# PC4
PC4_bio1 <- Miroc_126_2030_stand[[1]] * PC4_load[1,] 
PC4_bio2 <- Miroc_126_2030_stand[[2]] * PC4_load[2,] 
PC4_bio3 <- Miroc_126_2030_stand[[3]] * PC4_load[3,] 
PC4_bio4 <- Miroc_126_2030_stand[[4]] * PC4_load[4,] 
PC4_bio5 <- Miroc_126_2030_stand[[5]] * PC4_load[5,] 
PC4_bio6 <- Miroc_126_2030_stand[[6]] * PC4_load[6,] 
PC4_bio7 <- Miroc_126_2030_stand[[7]] * PC4_load[7,] 
PC4_bio8 <- Miroc_126_2030_stand[[8]] * PC4_load[8,] 
PC4_bio9 <- Miroc_126_2030_stand[[9]] * PC4_load[9,] 
PC4_bio10 <- Miroc_126_2030_stand[[10]] * PC4_load[10,] 
PC4_bio11 <- Miroc_126_2030_stand[[11]] * PC4_load[11,] 
PC4_bio12 <- Miroc_126_2030_stand[[12]] * PC4_load[12,] 
PC4_bio13 <- Miroc_126_2030_stand[[13]] * PC4_load[13,] 
PC4_bio14 <- Miroc_126_2030_stand[[14]] * PC4_load[14,] 
PC4_bio15 <- Miroc_126_2030_stand[[15]] * PC4_load[15,] 
PC4_bio16 <- Miroc_126_2030_stand[[16]] * PC4_load[16,] 
PC4_bio17 <- Miroc_126_2030_stand[[17]] * PC4_load[17,] 
PC4_bio18 <- Miroc_126_2030_stand[[18]] * PC4_load[18,] 
PC4_bio19 <- Miroc_126_2030_stand[[19]] * PC4_load[19,] 

PC4_bio_Miroc_126_2030 <- PC4_bio1 + PC4_bio2 + PC4_bio3 + PC4_bio4 + PC4_bio5 + PC4_bio6 + PC4_bio7 + PC4_bio8 + PC4_bio9 + PC4_bio10 + PC4_bio11 + PC4_bio12 + PC4_bio13 + PC4_bio14 + PC4_bio15 + PC4_bio16 + PC4_bio17 + PC4_bio18 + PC4_bio19

# PC5
PC5_bio1 <- Miroc_126_2030_stand[[1]] * PC5_load[1,] 
PC5_bio2 <- Miroc_126_2030_stand[[2]] * PC5_load[2,] 
PC5_bio3 <- Miroc_126_2030_stand[[3]] * PC5_load[3,] 
PC5_bio4 <- Miroc_126_2030_stand[[4]] * PC5_load[4,] 
PC5_bio5 <- Miroc_126_2030_stand[[5]] * PC5_load[5,] 
PC5_bio6 <- Miroc_126_2030_stand[[6]] * PC5_load[6,] 
PC5_bio7 <- Miroc_126_2030_stand[[7]] * PC5_load[7,] 
PC5_bio8 <- Miroc_126_2030_stand[[8]] * PC5_load[8,] 
PC5_bio9 <- Miroc_126_2030_stand[[9]] * PC5_load[9,] 
PC5_bio10 <- Miroc_126_2030_stand[[10]] * PC5_load[10,] 
PC5_bio11 <- Miroc_126_2030_stand[[11]] * PC5_load[11,] 
PC5_bio12 <- Miroc_126_2030_stand[[12]] * PC5_load[12,] 
PC5_bio13 <- Miroc_126_2030_stand[[13]] * PC5_load[13,] 
PC5_bio14 <- Miroc_126_2030_stand[[14]] * PC5_load[14,] 
PC5_bio15 <- Miroc_126_2030_stand[[15]] * PC5_load[15,] 
PC5_bio16 <- Miroc_126_2030_stand[[16]] * PC5_load[16,] 
PC5_bio17 <- Miroc_126_2030_stand[[17]] * PC5_load[17,] 
PC5_bio18 <- Miroc_126_2030_stand[[18]] * PC5_load[18,] 
PC5_bio19 <- Miroc_126_2030_stand[[19]] * PC5_load[19,] 

PC5_bio_Miroc_126_2030 <- PC5_bio1 + PC5_bio2 + PC5_bio3 + PC5_bio4 + PC5_bio5 + PC5_bio6 + PC5_bio7 + PC5_bio8 + PC5_bio9 + PC5_bio10 + PC5_bio11 + PC5_bio12 + PC5_bio13 + PC5_bio14 + PC5_bio15 + PC5_bio16 + PC5_bio17 + PC5_bio18 + PC5_bio19

dir.create("PCA_out_bio_Miroc_126_2030")
setwd("PCA_out_bio_Miroc_126_2030")
getwd()

# Write the rasters for the 5 PCs
writeRaster(PC1_bio_Miroc_126_2030, "__PC1.tif")
writeRaster(PC2_bio_Miroc_126_2030, "__PC2.tif")
writeRaster(PC3_bio_Miroc_126_2030, "__PC3.tif")
writeRaster(PC4_bio_Miroc_126_2030, "__PC4.tif")
writeRaster(PC5_bio_Miroc_126_2030, "__PC5.tif")

setwd(paste0("C:/Users/Public/Documents/Modelling/", species))

########################################## MIROC6 ##########################################
########################################## 2030 ############################################
########################################## 370 ##############################################
bioclim_Miroc_370_2030 <- stack("C:\\Users\\Public\\Documents\\Modelling\\wc2.1_2.5m_bioc_2021-2040\\wc2.1_2.5m_bioc_MIROC6_ssp370_2021-2040.tif")

# Crop and mask them to species-specific extent
CropClim_Miroc_370_2030=crop(bioclim_Miroc_370_2030, CurrentProj) #Crop before mask to speed up processing time
MaskClim_Miroc_370_2030=mask(CropClim_Miroc_370_2030, CurrentProj)

# Then standardise them based on the current mean and SD values
val1 <- MaskClim_Miroc_370_2030[[1]] - current.stats[1,1] 
Miroc_370_2030_sbio1 <- val1 / current.stats[2,1]

val2 <- MaskClim_Miroc_370_2030[[2]] - current.stats[1,2] 
Miroc_370_2030_sbio2 <- val2 / current.stats[2,2]

val3 <- MaskClim_Miroc_370_2030[[3]] - current.stats[1,3] 
Miroc_370_2030_sbio3 <- val3 / current.stats[2,3]

val4 <- MaskClim_Miroc_370_2030[[4]] - current.stats[1,4] 
Miroc_370_2030_sbio4 <- val4 / current.stats[2,4]

val5 <- MaskClim_Miroc_370_2030[[5]] - current.stats[1,5] 
Miroc_370_2030_sbio5 <- val5 / current.stats[2,5]

val6 <- MaskClim_Miroc_370_2030[[6]] - current.stats[1,6] 
Miroc_370_2030_sbio6 <- val6 / current.stats[2,6]

val7 <- MaskClim_Miroc_370_2030[[7]] - current.stats[1,7] 
Miroc_370_2030_sbio7 <- val7 / current.stats[2,7]

val8 <- MaskClim_Miroc_370_2030[[8]] - current.stats[1,8] 
Miroc_370_2030_sbio8 <- val8 / current.stats[2,8]

val9 <- MaskClim_Miroc_370_2030[[9]] - current.stats[1,9] 
Miroc_370_2030_sbio9 <- val9 / current.stats[2,9]

val10 <- MaskClim_Miroc_370_2030[[10]] - current.stats[1,10] 
Miroc_370_2030_sbio10 <- val10 / current.stats[2,10]

val11 <- MaskClim_Miroc_370_2030[[11]] - current.stats[1,11] 
Miroc_370_2030_sbio11 <- val11 / current.stats[2,11]

val12 <- MaskClim_Miroc_370_2030[[12]] - current.stats[1,12] 
Miroc_370_2030_sbio12 <- val12 / current.stats[2,12]

val13 <- MaskClim_Miroc_370_2030[[13]] - current.stats[1,13] 
Miroc_370_2030_sbio13 <- val13 / current.stats[2,13]

val14 <- MaskClim_Miroc_370_2030[[14]] - current.stats[1,14] 
Miroc_370_2030_sbio14 <- val14 / current.stats[2,14]

val15 <- MaskClim_Miroc_370_2030[[15]] - current.stats[1,15] 
Miroc_370_2030_sbio15 <- val15 / current.stats[2,15]

val16 <- MaskClim_Miroc_370_2030[[16]] - current.stats[1,16] 
Miroc_370_2030_sbio16 <- val16 / current.stats[2,16]

val17 <- MaskClim_Miroc_370_2030[[17]] - current.stats[1,17] 
Miroc_370_2030_sbio17 <- val17 / current.stats[2,17]

val18 <- MaskClim_Miroc_370_2030[[18]] - current.stats[1,18] 
Miroc_370_2030_sbio18 <- val18 / current.stats[2,18]

val19 <- MaskClim_Miroc_370_2030[[19]] - current.stats[1,19] 
Miroc_370_2030_sbio19 <- val19 / current.stats[2,19]

# Then multiply PC1 loadings for each variable by the standardised future rasters
Miroc_370_2030_stand <- stack(Miroc_370_2030_sbio1, Miroc_370_2030_sbio2, Miroc_370_2030_sbio3, Miroc_370_2030_sbio4, Miroc_370_2030_sbio5, Miroc_370_2030_sbio6,
                              Miroc_370_2030_sbio7, Miroc_370_2030_sbio8, Miroc_370_2030_sbio9, Miroc_370_2030_sbio10, Miroc_370_2030_sbio11, Miroc_370_2030_sbio12,
                              Miroc_370_2030_sbio13, Miroc_370_2030_sbio14, Miroc_370_2030_sbio15, Miroc_370_2030_sbio16, Miroc_370_2030_sbio17, Miroc_370_2030_sbio18,
                              Miroc_370_2030_sbio19)

# PC1
PC1_bio1 <- Miroc_370_2030_stand[[1]] * PC1_load[1,] 
PC1_bio2 <- Miroc_370_2030_stand[[2]] * PC1_load[2,] 
PC1_bio3 <- Miroc_370_2030_stand[[3]] * PC1_load[3,] 
PC1_bio4 <- Miroc_370_2030_stand[[4]] * PC1_load[4,] 
PC1_bio5 <- Miroc_370_2030_stand[[5]] * PC1_load[5,] 
PC1_bio6 <- Miroc_370_2030_stand[[6]] * PC1_load[6,] 
PC1_bio7 <- Miroc_370_2030_stand[[7]] * PC1_load[7,] 
PC1_bio8 <- Miroc_370_2030_stand[[8]] * PC1_load[8,] 
PC1_bio9 <- Miroc_370_2030_stand[[9]] * PC1_load[9,] 
PC1_bio10 <- Miroc_370_2030_stand[[10]] * PC1_load[10,] 
PC1_bio11 <- Miroc_370_2030_stand[[11]] * PC1_load[11,] 
PC1_bio12 <- Miroc_370_2030_stand[[12]] * PC1_load[12,] 
PC1_bio13 <- Miroc_370_2030_stand[[13]] * PC1_load[13,] 
PC1_bio14 <- Miroc_370_2030_stand[[14]] * PC1_load[14,] 
PC1_bio15 <- Miroc_370_2030_stand[[15]] * PC1_load[15,] 
PC1_bio16 <- Miroc_370_2030_stand[[16]] * PC1_load[16,] 
PC1_bio17 <- Miroc_370_2030_stand[[17]] * PC1_load[17,] 
PC1_bio18 <- Miroc_370_2030_stand[[18]] * PC1_load[18,] 
PC1_bio19 <- Miroc_370_2030_stand[[19]] * PC1_load[19,] 

PC1_bio_Miroc_370_2030 <- PC1_bio1 + PC1_bio2 + PC1_bio3 + PC1_bio4 + PC1_bio5 + PC1_bio6 + PC1_bio7 + PC1_bio8 + PC1_bio9 + PC1_bio10 + PC1_bio11 + PC1_bio12 + PC1_bio13 + PC1_bio14 + PC1_bio15 + PC1_bio16 + PC1_bio17 + PC1_bio18 + PC1_bio19

# PC2
PC2_bio1 <- Miroc_370_2030_stand[[1]] * PC2_load[1,] 
PC2_bio2 <- Miroc_370_2030_stand[[2]] * PC2_load[2,] 
PC2_bio3 <- Miroc_370_2030_stand[[3]] * PC2_load[3,] 
PC2_bio4 <- Miroc_370_2030_stand[[4]] * PC2_load[4,] 
PC2_bio5 <- Miroc_370_2030_stand[[5]] * PC2_load[5,] 
PC2_bio6 <- Miroc_370_2030_stand[[6]] * PC2_load[6,] 
PC2_bio7 <- Miroc_370_2030_stand[[7]] * PC2_load[7,] 
PC2_bio8 <- Miroc_370_2030_stand[[8]] * PC2_load[8,] 
PC2_bio9 <- Miroc_370_2030_stand[[9]] * PC2_load[9,] 
PC2_bio10 <- Miroc_370_2030_stand[[10]] * PC2_load[10,] 
PC2_bio11 <- Miroc_370_2030_stand[[11]] * PC2_load[11,] 
PC2_bio12 <- Miroc_370_2030_stand[[12]] * PC2_load[12,] 
PC2_bio13 <- Miroc_370_2030_stand[[13]] * PC2_load[13,] 
PC2_bio14 <- Miroc_370_2030_stand[[14]] * PC2_load[14,] 
PC2_bio15 <- Miroc_370_2030_stand[[15]] * PC2_load[15,] 
PC2_bio16 <- Miroc_370_2030_stand[[16]] * PC2_load[16,] 
PC2_bio17 <- Miroc_370_2030_stand[[17]] * PC2_load[17,] 
PC2_bio18 <- Miroc_370_2030_stand[[18]] * PC2_load[18,] 
PC2_bio19 <- Miroc_370_2030_stand[[19]] * PC2_load[19,] 

PC2_bio_Miroc_370_2030 <- PC2_bio1 + PC2_bio2 + PC2_bio3 + PC2_bio4 + PC2_bio5 + PC2_bio6 + PC2_bio7 + PC2_bio8 + PC2_bio9 + PC2_bio10 + PC2_bio11 + PC2_bio12 + PC2_bio13 + PC2_bio14 + PC2_bio15 + PC2_bio16 + PC2_bio17 + PC2_bio18 + PC2_bio19

# PC3
PC3_bio1 <- Miroc_370_2030_stand[[1]] * PC3_load[1,] 
PC3_bio2 <- Miroc_370_2030_stand[[2]] * PC3_load[2,] 
PC3_bio3 <- Miroc_370_2030_stand[[3]] * PC3_load[3,] 
PC3_bio4 <- Miroc_370_2030_stand[[4]] * PC3_load[4,] 
PC3_bio5 <- Miroc_370_2030_stand[[5]] * PC3_load[5,] 
PC3_bio6 <- Miroc_370_2030_stand[[6]] * PC3_load[6,] 
PC3_bio7 <- Miroc_370_2030_stand[[7]] * PC3_load[7,] 
PC3_bio8 <- Miroc_370_2030_stand[[8]] * PC3_load[8,] 
PC3_bio9 <- Miroc_370_2030_stand[[9]] * PC3_load[9,] 
PC3_bio10 <- Miroc_370_2030_stand[[10]] * PC3_load[10,] 
PC3_bio11 <- Miroc_370_2030_stand[[11]] * PC3_load[11,] 
PC3_bio12 <- Miroc_370_2030_stand[[12]] * PC3_load[12,] 
PC3_bio13 <- Miroc_370_2030_stand[[13]] * PC3_load[13,] 
PC3_bio14 <- Miroc_370_2030_stand[[14]] * PC3_load[14,] 
PC3_bio15 <- Miroc_370_2030_stand[[15]] * PC3_load[15,] 
PC3_bio16 <- Miroc_370_2030_stand[[16]] * PC3_load[16,] 
PC3_bio17 <- Miroc_370_2030_stand[[17]] * PC3_load[17,] 
PC3_bio18 <- Miroc_370_2030_stand[[18]] * PC3_load[18,] 
PC3_bio19 <- Miroc_370_2030_stand[[19]] * PC3_load[19,] 

PC3_bio_Miroc_370_2030 <- PC3_bio1 + PC3_bio2 + PC3_bio3 + PC3_bio4 + PC3_bio5 + PC3_bio6 + PC3_bio7 + PC3_bio8 + PC3_bio9 + PC3_bio10 + PC3_bio11 + PC3_bio12 + PC3_bio13 + PC3_bio14 + PC3_bio15 + PC3_bio16 + PC3_bio17 + PC3_bio18 + PC3_bio19

# PC4
PC4_bio1 <- Miroc_370_2030_stand[[1]] * PC4_load[1,] 
PC4_bio2 <- Miroc_370_2030_stand[[2]] * PC4_load[2,] 
PC4_bio3 <- Miroc_370_2030_stand[[3]] * PC4_load[3,] 
PC4_bio4 <- Miroc_370_2030_stand[[4]] * PC4_load[4,] 
PC4_bio5 <- Miroc_370_2030_stand[[5]] * PC4_load[5,] 
PC4_bio6 <- Miroc_370_2030_stand[[6]] * PC4_load[6,] 
PC4_bio7 <- Miroc_370_2030_stand[[7]] * PC4_load[7,] 
PC4_bio8 <- Miroc_370_2030_stand[[8]] * PC4_load[8,] 
PC4_bio9 <- Miroc_370_2030_stand[[9]] * PC4_load[9,] 
PC4_bio10 <- Miroc_370_2030_stand[[10]] * PC4_load[10,] 
PC4_bio11 <- Miroc_370_2030_stand[[11]] * PC4_load[11,] 
PC4_bio12 <- Miroc_370_2030_stand[[12]] * PC4_load[12,] 
PC4_bio13 <- Miroc_370_2030_stand[[13]] * PC4_load[13,] 
PC4_bio14 <- Miroc_370_2030_stand[[14]] * PC4_load[14,] 
PC4_bio15 <- Miroc_370_2030_stand[[15]] * PC4_load[15,] 
PC4_bio16 <- Miroc_370_2030_stand[[16]] * PC4_load[16,] 
PC4_bio17 <- Miroc_370_2030_stand[[17]] * PC4_load[17,] 
PC4_bio18 <- Miroc_370_2030_stand[[18]] * PC4_load[18,] 
PC4_bio19 <- Miroc_370_2030_stand[[19]] * PC4_load[19,] 

PC4_bio_Miroc_370_2030 <- PC4_bio1 + PC4_bio2 + PC4_bio3 + PC4_bio4 + PC4_bio5 + PC4_bio6 + PC4_bio7 + PC4_bio8 + PC4_bio9 + PC4_bio10 + PC4_bio11 + PC4_bio12 + PC4_bio13 + PC4_bio14 + PC4_bio15 + PC4_bio16 + PC4_bio17 + PC4_bio18 + PC4_bio19

# PC5
PC5_bio1 <- Miroc_370_2030_stand[[1]] * PC5_load[1,] 
PC5_bio2 <- Miroc_370_2030_stand[[2]] * PC5_load[2,] 
PC5_bio3 <- Miroc_370_2030_stand[[3]] * PC5_load[3,] 
PC5_bio4 <- Miroc_370_2030_stand[[4]] * PC5_load[4,] 
PC5_bio5 <- Miroc_370_2030_stand[[5]] * PC5_load[5,] 
PC5_bio6 <- Miroc_370_2030_stand[[6]] * PC5_load[6,] 
PC5_bio7 <- Miroc_370_2030_stand[[7]] * PC5_load[7,] 
PC5_bio8 <- Miroc_370_2030_stand[[8]] * PC5_load[8,] 
PC5_bio9 <- Miroc_370_2030_stand[[9]] * PC5_load[9,] 
PC5_bio10 <- Miroc_370_2030_stand[[10]] * PC5_load[10,] 
PC5_bio11 <- Miroc_370_2030_stand[[11]] * PC5_load[11,] 
PC5_bio12 <- Miroc_370_2030_stand[[12]] * PC5_load[12,] 
PC5_bio13 <- Miroc_370_2030_stand[[13]] * PC5_load[13,] 
PC5_bio14 <- Miroc_370_2030_stand[[14]] * PC5_load[14,] 
PC5_bio15 <- Miroc_370_2030_stand[[15]] * PC5_load[15,] 
PC5_bio16 <- Miroc_370_2030_stand[[16]] * PC5_load[16,] 
PC5_bio17 <- Miroc_370_2030_stand[[17]] * PC5_load[17,] 
PC5_bio18 <- Miroc_370_2030_stand[[18]] * PC5_load[18,] 
PC5_bio19 <- Miroc_370_2030_stand[[19]] * PC5_load[19,] 

PC5_bio_Miroc_370_2030 <- PC5_bio1 + PC5_bio2 + PC5_bio3 + PC5_bio4 + PC5_bio5 + PC5_bio6 + PC5_bio7 + PC5_bio8 + PC5_bio9 + PC5_bio10 + PC5_bio11 + PC5_bio12 + PC5_bio13 + PC5_bio14 + PC5_bio15 + PC5_bio16 + PC5_bio17 + PC5_bio18 + PC5_bio19

dir.create("PCA_out_bio_Miroc_370_2030")
setwd("PCA_out_bio_Miroc_370_2030")
getwd()

# Write the rasters for the 5 PCs
writeRaster(PC1_bio_Miroc_370_2030, "__PC1.tif")
writeRaster(PC2_bio_Miroc_370_2030, "__PC2.tif")
writeRaster(PC3_bio_Miroc_370_2030, "__PC3.tif")
writeRaster(PC4_bio_Miroc_370_2030, "__PC4.tif")
writeRaster(PC5_bio_Miroc_370_2030, "__PC5.tif")

setwd(paste0("C:/Users/Public/Documents/Modelling/", species))

########################################## MRI-ESM2-0 ##########################################
########################################## 2030 ############################################
########################################## 126 ##############################################
bioclim_Mri_126_2030 <- stack("C:\\Users\\Public\\Documents\\Modelling\\wc2.1_2.5m_bioc_2021-2040\\wc2.1_2.5m_bioc_MRI-ESM2-0_ssp126_2021-2040.tif")

# Crop and mask them to species-specific extent
CropClim_Mri_126_2030=crop(bioclim_Mri_126_2030, CurrentProj) #Crop before mask to speed up processing time
MaskClim_Mri_126_2030=mask(CropClim_Mri_126_2030, CurrentProj)

# Then standardise them based on the current mean and SD values
val1 <- MaskClim_Mri_126_2030[[1]] - current.stats[1,1] 
Mri_126_2030_sbio1 <- val1 / current.stats[2,1]

val2 <- MaskClim_Mri_126_2030[[2]] - current.stats[1,2] 
Mri_126_2030_sbio2 <- val2 / current.stats[2,2]

val3 <- MaskClim_Mri_126_2030[[3]] - current.stats[1,3] 
Mri_126_2030_sbio3 <- val3 / current.stats[2,3]

val4 <- MaskClim_Mri_126_2030[[4]] - current.stats[1,4] 
Mri_126_2030_sbio4 <- val4 / current.stats[2,4]

val5 <- MaskClim_Mri_126_2030[[5]] - current.stats[1,5] 
Mri_126_2030_sbio5 <- val5 / current.stats[2,5]

val6 <- MaskClim_Mri_126_2030[[6]] - current.stats[1,6] 
Mri_126_2030_sbio6 <- val6 / current.stats[2,6]

val7 <- MaskClim_Mri_126_2030[[7]] - current.stats[1,7] 
Mri_126_2030_sbio7 <- val7 / current.stats[2,7]

val8 <- MaskClim_Mri_126_2030[[8]] - current.stats[1,8] 
Mri_126_2030_sbio8 <- val8 / current.stats[2,8]

val9 <- MaskClim_Mri_126_2030[[9]] - current.stats[1,9] 
Mri_126_2030_sbio9 <- val9 / current.stats[2,9]

val10 <- MaskClim_Mri_126_2030[[10]] - current.stats[1,10] 
Mri_126_2030_sbio10 <- val10 / current.stats[2,10]

val11 <- MaskClim_Mri_126_2030[[11]] - current.stats[1,11] 
Mri_126_2030_sbio11 <- val11 / current.stats[2,11]

val12 <- MaskClim_Mri_126_2030[[12]] - current.stats[1,12] 
Mri_126_2030_sbio12 <- val12 / current.stats[2,12]

val13 <- MaskClim_Mri_126_2030[[13]] - current.stats[1,13] 
Mri_126_2030_sbio13 <- val13 / current.stats[2,13]

val14 <- MaskClim_Mri_126_2030[[14]] - current.stats[1,14] 
Mri_126_2030_sbio14 <- val14 / current.stats[2,14]

val15 <- MaskClim_Mri_126_2030[[15]] - current.stats[1,15] 
Mri_126_2030_sbio15 <- val15 / current.stats[2,15]

val16 <- MaskClim_Mri_126_2030[[16]] - current.stats[1,16] 
Mri_126_2030_sbio16 <- val16 / current.stats[2,16]

val17 <- MaskClim_Mri_126_2030[[17]] - current.stats[1,17] 
Mri_126_2030_sbio17 <- val17 / current.stats[2,17]

val18 <- MaskClim_Mri_126_2030[[18]] - current.stats[1,18] 
Mri_126_2030_sbio18 <- val18 / current.stats[2,18]

val19 <- MaskClim_Mri_126_2030[[19]] - current.stats[1,19] 
Mri_126_2030_sbio19 <- val19 / current.stats[2,19]

# Then multiply PC1 loadings for each variable by the standardised future rasters
Mri_126_2030_stand <- stack(Mri_126_2030_sbio1, Mri_126_2030_sbio2, Mri_126_2030_sbio3, Mri_126_2030_sbio4, Mri_126_2030_sbio5, Mri_126_2030_sbio6,
                            Mri_126_2030_sbio7, Mri_126_2030_sbio8, Mri_126_2030_sbio9, Mri_126_2030_sbio10, Mri_126_2030_sbio11, Mri_126_2030_sbio12,
                            Mri_126_2030_sbio13, Mri_126_2030_sbio14, Mri_126_2030_sbio15, Mri_126_2030_sbio16, Mri_126_2030_sbio17, Mri_126_2030_sbio18,
                            Mri_126_2030_sbio19)

# PC1
PC1_bio1 <- Mri_126_2030_stand[[1]] * PC1_load[1,] 
PC1_bio2 <- Mri_126_2030_stand[[2]] * PC1_load[2,] 
PC1_bio3 <- Mri_126_2030_stand[[3]] * PC1_load[3,] 
PC1_bio4 <- Mri_126_2030_stand[[4]] * PC1_load[4,] 
PC1_bio5 <- Mri_126_2030_stand[[5]] * PC1_load[5,] 
PC1_bio6 <- Mri_126_2030_stand[[6]] * PC1_load[6,] 
PC1_bio7 <- Mri_126_2030_stand[[7]] * PC1_load[7,] 
PC1_bio8 <- Mri_126_2030_stand[[8]] * PC1_load[8,] 
PC1_bio9 <- Mri_126_2030_stand[[9]] * PC1_load[9,] 
PC1_bio10 <- Mri_126_2030_stand[[10]] * PC1_load[10,] 
PC1_bio11 <- Mri_126_2030_stand[[11]] * PC1_load[11,] 
PC1_bio12 <- Mri_126_2030_stand[[12]] * PC1_load[12,] 
PC1_bio13 <- Mri_126_2030_stand[[13]] * PC1_load[13,] 
PC1_bio14 <- Mri_126_2030_stand[[14]] * PC1_load[14,] 
PC1_bio15 <- Mri_126_2030_stand[[15]] * PC1_load[15,] 
PC1_bio16 <- Mri_126_2030_stand[[16]] * PC1_load[16,] 
PC1_bio17 <- Mri_126_2030_stand[[17]] * PC1_load[17,] 
PC1_bio18 <- Mri_126_2030_stand[[18]] * PC1_load[18,] 
PC1_bio19 <- Mri_126_2030_stand[[19]] * PC1_load[19,] 

PC1_bio_Mri_126_2030 <- PC1_bio1 + PC1_bio2 + PC1_bio3 + PC1_bio4 + PC1_bio5 + PC1_bio6 + PC1_bio7 + PC1_bio8 + PC1_bio9 + PC1_bio10 + PC1_bio11 + PC1_bio12 + PC1_bio13 + PC1_bio14 + PC1_bio15 + PC1_bio16 + PC1_bio17 + PC1_bio18 + PC1_bio19

# PC2
PC2_bio1 <- Mri_126_2030_stand[[1]] * PC2_load[1,] 
PC2_bio2 <- Mri_126_2030_stand[[2]] * PC2_load[2,] 
PC2_bio3 <- Mri_126_2030_stand[[3]] * PC2_load[3,] 
PC2_bio4 <- Mri_126_2030_stand[[4]] * PC2_load[4,] 
PC2_bio5 <- Mri_126_2030_stand[[5]] * PC2_load[5,] 
PC2_bio6 <- Mri_126_2030_stand[[6]] * PC2_load[6,] 
PC2_bio7 <- Mri_126_2030_stand[[7]] * PC2_load[7,] 
PC2_bio8 <- Mri_126_2030_stand[[8]] * PC2_load[8,] 
PC2_bio9 <- Mri_126_2030_stand[[9]] * PC2_load[9,] 
PC2_bio10 <- Mri_126_2030_stand[[10]] * PC2_load[10,] 
PC2_bio11 <- Mri_126_2030_stand[[11]] * PC2_load[11,] 
PC2_bio12 <- Mri_126_2030_stand[[12]] * PC2_load[12,] 
PC2_bio13 <- Mri_126_2030_stand[[13]] * PC2_load[13,] 
PC2_bio14 <- Mri_126_2030_stand[[14]] * PC2_load[14,] 
PC2_bio15 <- Mri_126_2030_stand[[15]] * PC2_load[15,] 
PC2_bio16 <- Mri_126_2030_stand[[16]] * PC2_load[16,] 
PC2_bio17 <- Mri_126_2030_stand[[17]] * PC2_load[17,] 
PC2_bio18 <- Mri_126_2030_stand[[18]] * PC2_load[18,] 
PC2_bio19 <- Mri_126_2030_stand[[19]] * PC2_load[19,] 

PC2_bio_Mri_126_2030 <- PC2_bio1 + PC2_bio2 + PC2_bio3 + PC2_bio4 + PC2_bio5 + PC2_bio6 + PC2_bio7 + PC2_bio8 + PC2_bio9 + PC2_bio10 + PC2_bio11 + PC2_bio12 + PC2_bio13 + PC2_bio14 + PC2_bio15 + PC2_bio16 + PC2_bio17 + PC2_bio18 + PC2_bio19

# PC3
PC3_bio1 <- Mri_126_2030_stand[[1]] * PC3_load[1,] 
PC3_bio2 <- Mri_126_2030_stand[[2]] * PC3_load[2,] 
PC3_bio3 <- Mri_126_2030_stand[[3]] * PC3_load[3,] 
PC3_bio4 <- Mri_126_2030_stand[[4]] * PC3_load[4,] 
PC3_bio5 <- Mri_126_2030_stand[[5]] * PC3_load[5,] 
PC3_bio6 <- Mri_126_2030_stand[[6]] * PC3_load[6,] 
PC3_bio7 <- Mri_126_2030_stand[[7]] * PC3_load[7,] 
PC3_bio8 <- Mri_126_2030_stand[[8]] * PC3_load[8,] 
PC3_bio9 <- Mri_126_2030_stand[[9]] * PC3_load[9,] 
PC3_bio10 <- Mri_126_2030_stand[[10]] * PC3_load[10,] 
PC3_bio11 <- Mri_126_2030_stand[[11]] * PC3_load[11,] 
PC3_bio12 <- Mri_126_2030_stand[[12]] * PC3_load[12,] 
PC3_bio13 <- Mri_126_2030_stand[[13]] * PC3_load[13,] 
PC3_bio14 <- Mri_126_2030_stand[[14]] * PC3_load[14,] 
PC3_bio15 <- Mri_126_2030_stand[[15]] * PC3_load[15,] 
PC3_bio16 <- Mri_126_2030_stand[[16]] * PC3_load[16,] 
PC3_bio17 <- Mri_126_2030_stand[[17]] * PC3_load[17,] 
PC3_bio18 <- Mri_126_2030_stand[[18]] * PC3_load[18,] 
PC3_bio19 <- Mri_126_2030_stand[[19]] * PC3_load[19,] 

PC3_bio_Mri_126_2030 <- PC3_bio1 + PC3_bio2 + PC3_bio3 + PC3_bio4 + PC3_bio5 + PC3_bio6 + PC3_bio7 + PC3_bio8 + PC3_bio9 + PC3_bio10 + PC3_bio11 + PC3_bio12 + PC3_bio13 + PC3_bio14 + PC3_bio15 + PC3_bio16 + PC3_bio17 + PC3_bio18 + PC3_bio19

# PC4
PC4_bio1 <- Mri_126_2030_stand[[1]] * PC4_load[1,] 
PC4_bio2 <- Mri_126_2030_stand[[2]] * PC4_load[2,] 
PC4_bio3 <- Mri_126_2030_stand[[3]] * PC4_load[3,] 
PC4_bio4 <- Mri_126_2030_stand[[4]] * PC4_load[4,] 
PC4_bio5 <- Mri_126_2030_stand[[5]] * PC4_load[5,] 
PC4_bio6 <- Mri_126_2030_stand[[6]] * PC4_load[6,] 
PC4_bio7 <- Mri_126_2030_stand[[7]] * PC4_load[7,] 
PC4_bio8 <- Mri_126_2030_stand[[8]] * PC4_load[8,] 
PC4_bio9 <- Mri_126_2030_stand[[9]] * PC4_load[9,] 
PC4_bio10 <- Mri_126_2030_stand[[10]] * PC4_load[10,] 
PC4_bio11 <- Mri_126_2030_stand[[11]] * PC4_load[11,] 
PC4_bio12 <- Mri_126_2030_stand[[12]] * PC4_load[12,] 
PC4_bio13 <- Mri_126_2030_stand[[13]] * PC4_load[13,] 
PC4_bio14 <- Mri_126_2030_stand[[14]] * PC4_load[14,] 
PC4_bio15 <- Mri_126_2030_stand[[15]] * PC4_load[15,] 
PC4_bio16 <- Mri_126_2030_stand[[16]] * PC4_load[16,] 
PC4_bio17 <- Mri_126_2030_stand[[17]] * PC4_load[17,] 
PC4_bio18 <- Mri_126_2030_stand[[18]] * PC4_load[18,] 
PC4_bio19 <- Mri_126_2030_stand[[19]] * PC4_load[19,] 

PC4_bio_Mri_126_2030 <- PC4_bio1 + PC4_bio2 + PC4_bio3 + PC4_bio4 + PC4_bio5 + PC4_bio6 + PC4_bio7 + PC4_bio8 + PC4_bio9 + PC4_bio10 + PC4_bio11 + PC4_bio12 + PC4_bio13 + PC4_bio14 + PC4_bio15 + PC4_bio16 + PC4_bio17 + PC4_bio18 + PC4_bio19

# PC5
PC5_bio1 <- Mri_126_2030_stand[[1]] * PC5_load[1,] 
PC5_bio2 <- Mri_126_2030_stand[[2]] * PC5_load[2,] 
PC5_bio3 <- Mri_126_2030_stand[[3]] * PC5_load[3,] 
PC5_bio4 <- Mri_126_2030_stand[[4]] * PC5_load[4,] 
PC5_bio5 <- Mri_126_2030_stand[[5]] * PC5_load[5,] 
PC5_bio6 <- Mri_126_2030_stand[[6]] * PC5_load[6,] 
PC5_bio7 <- Mri_126_2030_stand[[7]] * PC5_load[7,] 
PC5_bio8 <- Mri_126_2030_stand[[8]] * PC5_load[8,] 
PC5_bio9 <- Mri_126_2030_stand[[9]] * PC5_load[9,] 
PC5_bio10 <- Mri_126_2030_stand[[10]] * PC5_load[10,] 
PC5_bio11 <- Mri_126_2030_stand[[11]] * PC5_load[11,] 
PC5_bio12 <- Mri_126_2030_stand[[12]] * PC5_load[12,] 
PC5_bio13 <- Mri_126_2030_stand[[13]] * PC5_load[13,] 
PC5_bio14 <- Mri_126_2030_stand[[14]] * PC5_load[14,] 
PC5_bio15 <- Mri_126_2030_stand[[15]] * PC5_load[15,] 
PC5_bio16 <- Mri_126_2030_stand[[16]] * PC5_load[16,] 
PC5_bio17 <- Mri_126_2030_stand[[17]] * PC5_load[17,] 
PC5_bio18 <- Mri_126_2030_stand[[18]] * PC5_load[18,] 
PC5_bio19 <- Mri_126_2030_stand[[19]] * PC5_load[19,] 

PC5_bio_Mri_126_2030 <- PC5_bio1 + PC5_bio2 + PC5_bio3 + PC5_bio4 + PC5_bio5 + PC5_bio6 + PC5_bio7 + PC5_bio8 + PC5_bio9 + PC5_bio10 + PC5_bio11 + PC5_bio12 + PC5_bio13 + PC5_bio14 + PC5_bio15 + PC5_bio16 + PC5_bio17 + PC5_bio18 + PC5_bio19

dir.create("PCA_out_bio_Mri_126_2030")
setwd("PCA_out_bio_Mri_126_2030")
getwd()

# Write the rasters for the 5 PCs
writeRaster(PC1_bio_Mri_126_2030, "__PC1.tif")
writeRaster(PC2_bio_Mri_126_2030, "__PC2.tif")
writeRaster(PC3_bio_Mri_126_2030, "__PC3.tif")
writeRaster(PC4_bio_Mri_126_2030, "__PC4.tif")
writeRaster(PC5_bio_Mri_126_2030, "__PC5.tif")

setwd(paste0("C:/Users/Public/Documents/Modelling/", species))

########################################## MRI-ESM2-0 ##########################################
########################################## 2030 ############################################
########################################## 370 ##############################################
bioclim_Mri_370_2030 <- stack("C:\\Users\\Public\\Documents\\Modelling\\wc2.1_2.5m_bioc_2021-2040\\wc2.1_2.5m_bioc_MRI-ESM2-0_ssp370_2021-2040.tif")

# Crop and mask them to species-specific extent
CropClim_Mri_370_2030=crop(bioclim_Mri_370_2030, CurrentProj) #Crop before mask to speed up processing time
MaskClim_Mri_370_2030=mask(CropClim_Mri_370_2030, CurrentProj)

# Then standardise them based on the current mean and SD values
val1 <- MaskClim_Mri_370_2030[[1]] - current.stats[1,1] 
Mri_370_2030_sbio1 <- val1 / current.stats[2,1]

val2 <- MaskClim_Mri_370_2030[[2]] - current.stats[1,2] 
Mri_370_2030_sbio2 <- val2 / current.stats[2,2]

val3 <- MaskClim_Mri_370_2030[[3]] - current.stats[1,3] 
Mri_370_2030_sbio3 <- val3 / current.stats[2,3]

val4 <- MaskClim_Mri_370_2030[[4]] - current.stats[1,4] 
Mri_370_2030_sbio4 <- val4 / current.stats[2,4]

val5 <- MaskClim_Mri_370_2030[[5]] - current.stats[1,5] 
Mri_370_2030_sbio5 <- val5 / current.stats[2,5]

val6 <- MaskClim_Mri_370_2030[[6]] - current.stats[1,6] 
Mri_370_2030_sbio6 <- val6 / current.stats[2,6]

val7 <- MaskClim_Mri_370_2030[[7]] - current.stats[1,7] 
Mri_370_2030_sbio7 <- val7 / current.stats[2,7]

val8 <- MaskClim_Mri_370_2030[[8]] - current.stats[1,8] 
Mri_370_2030_sbio8 <- val8 / current.stats[2,8]

val9 <- MaskClim_Mri_370_2030[[9]] - current.stats[1,9] 
Mri_370_2030_sbio9 <- val9 / current.stats[2,9]

val10 <- MaskClim_Mri_370_2030[[10]] - current.stats[1,10] 
Mri_370_2030_sbio10 <- val10 / current.stats[2,10]

val11 <- MaskClim_Mri_370_2030[[11]] - current.stats[1,11] 
Mri_370_2030_sbio11 <- val11 / current.stats[2,11]

val12 <- MaskClim_Mri_370_2030[[12]] - current.stats[1,12] 
Mri_370_2030_sbio12 <- val12 / current.stats[2,12]

val13 <- MaskClim_Mri_370_2030[[13]] - current.stats[1,13] 
Mri_370_2030_sbio13 <- val13 / current.stats[2,13]

val14 <- MaskClim_Mri_370_2030[[14]] - current.stats[1,14] 
Mri_370_2030_sbio14 <- val14 / current.stats[2,14]

val15 <- MaskClim_Mri_370_2030[[15]] - current.stats[1,15] 
Mri_370_2030_sbio15 <- val15 / current.stats[2,15]

val16 <- MaskClim_Mri_370_2030[[16]] - current.stats[1,16] 
Mri_370_2030_sbio16 <- val16 / current.stats[2,16]

val17 <- MaskClim_Mri_370_2030[[17]] - current.stats[1,17] 
Mri_370_2030_sbio17 <- val17 / current.stats[2,17]

val18 <- MaskClim_Mri_370_2030[[18]] - current.stats[1,18] 
Mri_370_2030_sbio18 <- val18 / current.stats[2,18]

val19 <- MaskClim_Mri_370_2030[[19]] - current.stats[1,19] 
Mri_370_2030_sbio19 <- val19 / current.stats[2,19]

# Then multiply PC1 loadings for each variable by the standardised future rasters
Mri_370_2030_stand <- stack(Mri_370_2030_sbio1, Mri_370_2030_sbio2, Mri_370_2030_sbio3, Mri_370_2030_sbio4, Mri_370_2030_sbio5, Mri_370_2030_sbio6,
                            Mri_370_2030_sbio7, Mri_370_2030_sbio8, Mri_370_2030_sbio9, Mri_370_2030_sbio10, Mri_370_2030_sbio11, Mri_370_2030_sbio12,
                            Mri_370_2030_sbio13, Mri_370_2030_sbio14, Mri_370_2030_sbio15, Mri_370_2030_sbio16, Mri_370_2030_sbio17, Mri_370_2030_sbio18,
                            Mri_370_2030_sbio19)

# PC1
PC1_bio1 <- Mri_370_2030_stand[[1]] * PC1_load[1,] 
PC1_bio2 <- Mri_370_2030_stand[[2]] * PC1_load[2,] 
PC1_bio3 <- Mri_370_2030_stand[[3]] * PC1_load[3,] 
PC1_bio4 <- Mri_370_2030_stand[[4]] * PC1_load[4,] 
PC1_bio5 <- Mri_370_2030_stand[[5]] * PC1_load[5,] 
PC1_bio6 <- Mri_370_2030_stand[[6]] * PC1_load[6,] 
PC1_bio7 <- Mri_370_2030_stand[[7]] * PC1_load[7,] 
PC1_bio8 <- Mri_370_2030_stand[[8]] * PC1_load[8,] 
PC1_bio9 <- Mri_370_2030_stand[[9]] * PC1_load[9,] 
PC1_bio10 <- Mri_370_2030_stand[[10]] * PC1_load[10,] 
PC1_bio11 <- Mri_370_2030_stand[[11]] * PC1_load[11,] 
PC1_bio12 <- Mri_370_2030_stand[[12]] * PC1_load[12,] 
PC1_bio13 <- Mri_370_2030_stand[[13]] * PC1_load[13,] 
PC1_bio14 <- Mri_370_2030_stand[[14]] * PC1_load[14,] 
PC1_bio15 <- Mri_370_2030_stand[[15]] * PC1_load[15,] 
PC1_bio16 <- Mri_370_2030_stand[[16]] * PC1_load[16,] 
PC1_bio17 <- Mri_370_2030_stand[[17]] * PC1_load[17,] 
PC1_bio18 <- Mri_370_2030_stand[[18]] * PC1_load[18,] 
PC1_bio19 <- Mri_370_2030_stand[[19]] * PC1_load[19,] 

PC1_bio_Mri_370_2030 <- PC1_bio1 + PC1_bio2 + PC1_bio3 + PC1_bio4 + PC1_bio5 + PC1_bio6 + PC1_bio7 + PC1_bio8 + PC1_bio9 + PC1_bio10 + PC1_bio11 + PC1_bio12 + PC1_bio13 + PC1_bio14 + PC1_bio15 + PC1_bio16 + PC1_bio17 + PC1_bio18 + PC1_bio19

# PC2
PC2_bio1 <- Mri_370_2030_stand[[1]] * PC2_load[1,] 
PC2_bio2 <- Mri_370_2030_stand[[2]] * PC2_load[2,] 
PC2_bio3 <- Mri_370_2030_stand[[3]] * PC2_load[3,] 
PC2_bio4 <- Mri_370_2030_stand[[4]] * PC2_load[4,] 
PC2_bio5 <- Mri_370_2030_stand[[5]] * PC2_load[5,] 
PC2_bio6 <- Mri_370_2030_stand[[6]] * PC2_load[6,] 
PC2_bio7 <- Mri_370_2030_stand[[7]] * PC2_load[7,] 
PC2_bio8 <- Mri_370_2030_stand[[8]] * PC2_load[8,] 
PC2_bio9 <- Mri_370_2030_stand[[9]] * PC2_load[9,] 
PC2_bio10 <- Mri_370_2030_stand[[10]] * PC2_load[10,] 
PC2_bio11 <- Mri_370_2030_stand[[11]] * PC2_load[11,] 
PC2_bio12 <- Mri_370_2030_stand[[12]] * PC2_load[12,] 
PC2_bio13 <- Mri_370_2030_stand[[13]] * PC2_load[13,] 
PC2_bio14 <- Mri_370_2030_stand[[14]] * PC2_load[14,] 
PC2_bio15 <- Mri_370_2030_stand[[15]] * PC2_load[15,] 
PC2_bio16 <- Mri_370_2030_stand[[16]] * PC2_load[16,] 
PC2_bio17 <- Mri_370_2030_stand[[17]] * PC2_load[17,] 
PC2_bio18 <- Mri_370_2030_stand[[18]] * PC2_load[18,] 
PC2_bio19 <- Mri_370_2030_stand[[19]] * PC2_load[19,] 

PC2_bio_Mri_370_2030 <- PC2_bio1 + PC2_bio2 + PC2_bio3 + PC2_bio4 + PC2_bio5 + PC2_bio6 + PC2_bio7 + PC2_bio8 + PC2_bio9 + PC2_bio10 + PC2_bio11 + PC2_bio12 + PC2_bio13 + PC2_bio14 + PC2_bio15 + PC2_bio16 + PC2_bio17 + PC2_bio18 + PC2_bio19

# PC3
PC3_bio1 <- Mri_370_2030_stand[[1]] * PC3_load[1,] 
PC3_bio2 <- Mri_370_2030_stand[[2]] * PC3_load[2,] 
PC3_bio3 <- Mri_370_2030_stand[[3]] * PC3_load[3,] 
PC3_bio4 <- Mri_370_2030_stand[[4]] * PC3_load[4,] 
PC3_bio5 <- Mri_370_2030_stand[[5]] * PC3_load[5,] 
PC3_bio6 <- Mri_370_2030_stand[[6]] * PC3_load[6,] 
PC3_bio7 <- Mri_370_2030_stand[[7]] * PC3_load[7,] 
PC3_bio8 <- Mri_370_2030_stand[[8]] * PC3_load[8,] 
PC3_bio9 <- Mri_370_2030_stand[[9]] * PC3_load[9,] 
PC3_bio10 <- Mri_370_2030_stand[[10]] * PC3_load[10,] 
PC3_bio11 <- Mri_370_2030_stand[[11]] * PC3_load[11,] 
PC3_bio12 <- Mri_370_2030_stand[[12]] * PC3_load[12,] 
PC3_bio13 <- Mri_370_2030_stand[[13]] * PC3_load[13,] 
PC3_bio14 <- Mri_370_2030_stand[[14]] * PC3_load[14,] 
PC3_bio15 <- Mri_370_2030_stand[[15]] * PC3_load[15,] 
PC3_bio16 <- Mri_370_2030_stand[[16]] * PC3_load[16,] 
PC3_bio17 <- Mri_370_2030_stand[[17]] * PC3_load[17,] 
PC3_bio18 <- Mri_370_2030_stand[[18]] * PC3_load[18,] 
PC3_bio19 <- Mri_370_2030_stand[[19]] * PC3_load[19,] 

PC3_bio_Mri_370_2030 <- PC3_bio1 + PC3_bio2 + PC3_bio3 + PC3_bio4 + PC3_bio5 + PC3_bio6 + PC3_bio7 + PC3_bio8 + PC3_bio9 + PC3_bio10 + PC3_bio11 + PC3_bio12 + PC3_bio13 + PC3_bio14 + PC3_bio15 + PC3_bio16 + PC3_bio17 + PC3_bio18 + PC3_bio19

# PC4
PC4_bio1 <- Mri_370_2030_stand[[1]] * PC4_load[1,] 
PC4_bio2 <- Mri_370_2030_stand[[2]] * PC4_load[2,] 
PC4_bio3 <- Mri_370_2030_stand[[3]] * PC4_load[3,] 
PC4_bio4 <- Mri_370_2030_stand[[4]] * PC4_load[4,] 
PC4_bio5 <- Mri_370_2030_stand[[5]] * PC4_load[5,] 
PC4_bio6 <- Mri_370_2030_stand[[6]] * PC4_load[6,] 
PC4_bio7 <- Mri_370_2030_stand[[7]] * PC4_load[7,] 
PC4_bio8 <- Mri_370_2030_stand[[8]] * PC4_load[8,] 
PC4_bio9 <- Mri_370_2030_stand[[9]] * PC4_load[9,] 
PC4_bio10 <- Mri_370_2030_stand[[10]] * PC4_load[10,] 
PC4_bio11 <- Mri_370_2030_stand[[11]] * PC4_load[11,] 
PC4_bio12 <- Mri_370_2030_stand[[12]] * PC4_load[12,] 
PC4_bio13 <- Mri_370_2030_stand[[13]] * PC4_load[13,] 
PC4_bio14 <- Mri_370_2030_stand[[14]] * PC4_load[14,] 
PC4_bio15 <- Mri_370_2030_stand[[15]] * PC4_load[15,] 
PC4_bio16 <- Mri_370_2030_stand[[16]] * PC4_load[16,] 
PC4_bio17 <- Mri_370_2030_stand[[17]] * PC4_load[17,] 
PC4_bio18 <- Mri_370_2030_stand[[18]] * PC4_load[18,] 
PC4_bio19 <- Mri_370_2030_stand[[19]] * PC4_load[19,] 

PC4_bio_Mri_370_2030 <- PC4_bio1 + PC4_bio2 + PC4_bio3 + PC4_bio4 + PC4_bio5 + PC4_bio6 + PC4_bio7 + PC4_bio8 + PC4_bio9 + PC4_bio10 + PC4_bio11 + PC4_bio12 + PC4_bio13 + PC4_bio14 + PC4_bio15 + PC4_bio16 + PC4_bio17 + PC4_bio18 + PC4_bio19

# PC5
PC5_bio1 <- Mri_370_2030_stand[[1]] * PC5_load[1,] 
PC5_bio2 <- Mri_370_2030_stand[[2]] * PC5_load[2,] 
PC5_bio3 <- Mri_370_2030_stand[[3]] * PC5_load[3,] 
PC5_bio4 <- Mri_370_2030_stand[[4]] * PC5_load[4,] 
PC5_bio5 <- Mri_370_2030_stand[[5]] * PC5_load[5,] 
PC5_bio6 <- Mri_370_2030_stand[[6]] * PC5_load[6,] 
PC5_bio7 <- Mri_370_2030_stand[[7]] * PC5_load[7,] 
PC5_bio8 <- Mri_370_2030_stand[[8]] * PC5_load[8,] 
PC5_bio9 <- Mri_370_2030_stand[[9]] * PC5_load[9,] 
PC5_bio10 <- Mri_370_2030_stand[[10]] * PC5_load[10,] 
PC5_bio11 <- Mri_370_2030_stand[[11]] * PC5_load[11,] 
PC5_bio12 <- Mri_370_2030_stand[[12]] * PC5_load[12,] 
PC5_bio13 <- Mri_370_2030_stand[[13]] * PC5_load[13,] 
PC5_bio14 <- Mri_370_2030_stand[[14]] * PC5_load[14,] 
PC5_bio15 <- Mri_370_2030_stand[[15]] * PC5_load[15,] 
PC5_bio16 <- Mri_370_2030_stand[[16]] * PC5_load[16,] 
PC5_bio17 <- Mri_370_2030_stand[[17]] * PC5_load[17,] 
PC5_bio18 <- Mri_370_2030_stand[[18]] * PC5_load[18,] 
PC5_bio19 <- Mri_370_2030_stand[[19]] * PC5_load[19,] 

PC5_bio_Mri_370_2030 <- PC5_bio1 + PC5_bio2 + PC5_bio3 + PC5_bio4 + PC5_bio5 + PC5_bio6 + PC5_bio7 + PC5_bio8 + PC5_bio9 + PC5_bio10 + PC5_bio11 + PC5_bio12 + PC5_bio13 + PC5_bio14 + PC5_bio15 + PC5_bio16 + PC5_bio17 + PC5_bio18 + PC5_bio19

dir.create("PCA_out_bio_Mri_370_2030")
setwd("PCA_out_bio_Mri_370_2030")
getwd()

# Write the rasters for the 5 PCs
writeRaster(PC1_bio_Mri_370_2030, "PC1.tif")
writeRaster(PC2_bio_Mri_370_2030, "PC2.tif")
writeRaster(PC3_bio_Mri_370_2030, "PC3.tif")
writeRaster(PC4_bio_Mri_370_2030, "PC4.tif")
writeRaster(PC5_bio_Mri_370_2030, "PC5.tif")

setwd(paste0("C:/Users/Public/Documents/Modelling/", species))
