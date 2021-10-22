# 00 Setting up files and folders for modelling work

# Sarah Dalrymple
# July 2021 updated 29 Sep 2021

# Script for setting up the correct working directory and populating wd with scripts and other necessary files

### INSERT SPECIES NAME BETWEEN QUOTES ###
species <- ("Melampyrum_sylvaticum")

# set working directory
setwd("C:/Users/Public/Documents/Modelling")

# create new folder 
speciesfolder <- dir.create(species)

# copy scripts and maxent file into working directory
file.copy(from = "C:/Users/Public/Documents/Modelling/maxent.jar", 
          to = paste0("C:/Users/Public/Documents/Modelling/",species,"/maxent.jar"))
file.copy(from = "C:/Users/Public/Documents/Modelling/01_CleaningGBIFdata.R", 
          to = paste0("C:/Users/Public/Documents/Modelling/",species,"/01_CleaningGBIFdata_",species,".R"))
file.copy(from = "C:/Users/Public/Documents/Modelling/02_PCA_Climate_indices.R", 
          to = paste0("C:/Users/Public/Documents/Modelling/",species,"/02_PCA_Climate_indices_",species,".R"))
file.copy(from = "C:/Users/Public/Documents/Modelling/03_Occ_filtering.R", 
          to = paste0("C:/Users/Public/Documents/Modelling/",species,"/03_Occ_filtering_",species,".R"))
file.copy(from = "C:/Users/Public/Documents/Modelling/04_biomod.R", 
          to = paste0("C:/Users/Public/Documents/Modelling/",species,"/04_biomod_",species,".R"))
file.copy(from = "C:/Users/Public/Documents/Modelling/05_extracting_climate_suitability.R", 
          to = paste0("C:/Users/Public/Documents/Modelling/",species,"/05_extracting_climate_suitability_",species,".R"))

# after the wd is created and set up with the necessary files, do the manual occurrence download from GBIF
# this is described in script entitled '01_CleaningGBIFdata' which has now been moved to the wd and the
# species name appended to the script

# packages used by this analysis can be installed here but doesn't need to be done every time
install.packages("CoordinateCleaner")
install.packages("rangeBuilder")
install.packages("maptools")
install.packages("raster")
install.packages("sp")
install.packages("rgdal")
install.packages("RStoolbox")
install.packages("dplyr")
install.packages("adehabitatHR")
install.packages("maps")
install.packages("dismo")
install.packages("ecodist")
install.packages("biomod2")
install.packages("rgeos")
install.packages("abind")
install.packages("gridExtra")
install.packages("lattice")
install.packages("ecospat")


# if you get memory problems (error messages will tell you that vectors of certain sizes can't be handled),
# try using the code below to work out what the current memory allocation is, and then increase it
# if you've asked for more than is available, it will return an error
memory.limit()
memory.limit(size = 112000)
