#Species Distribution Modeling with biomod2 for BES funded project

# Joe Bellis
# June 2021

# as before, use the following code to tidy up 'environment' if necessary
rm(list=ls())

setwd("C:/Users/Public/Documents/Modelling/...")
getwd()

library(adehabitatHR)
library(sp)
library(maps)
library(maptools)
library(rgdal)
library(raster)
library(dismo)
library(ecodist)
library(biomod2)
library(rgeos)
library(abind)
library(gridExtra)
library(lattice)
library(ecospat)
library(rangeBuilder)
library(dplyr)


# DATA PREPARATION
ProjW = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" #map projection (EPSG: 4326)

# Load Principal components for use as expl variables
# NOTE: THE NUMBER OF PC's WILL DEPEND ON THE AMOUNT OF VARIATION EXPLAINED BY THE ACCUMULATED PCs, i.e. >90%
# will need to add/remove PCs accordingly, using "#" if less
# use the 'variance.csv' and look at the 'cumulative proportion' row to work out how many components should be included.
PC1 <- raster("PCA_Out_current/__PC1.tif")
PC2 <- raster("PCA_Out_current/__PC2.tif")
PC3 <- raster("PCA_Out_current/__PC3.tif")
PC4 <- raster("PCA_Out_current/__PC4.tif")
PC5 <- raster("PCA_Out_current/__PC5.tif")

# Stack them
clim <-stack(PC1, PC2, PC3, PC4, PC5)

# Assign species name/code
wd <-paste0(getwd())
mySpc <- basename(wd)

# Load presences
Presences <- read.csv("presences.csv")

# Get data into biomod2 format
my.resp.xy <- Presences[,2:3]

# Add a '1' to presences column for biomod2 formatting (my.resp.var argument)
Presences['Binary'] = '1'
my.resp.var <- as.data.frame(Presences[,4])
colnames(my.resp.var) <- "Binary"

# Get the number of occurrences (this will be slightly less than sp occ file due to NAs)
xy <- Presences[,2:3]
df <- Presences
Occ <- SpatialPointsDataFrame(coords=xy, data=df, proj4string = CRS(ProjW))
Points_on_rast <- extract(PC1, Occ)
Removed_NAs <- as.data.frame(na.omit(Points_on_rast)) # this is the total number that will be used in the SDM
Removed_NAs['Random'] = '1'
Column_1 <- as.data.frame(Removed_NAs$Random)
Column_2 <- as.data.frame(Removed_NAs$Random)
colnames(Column_1) <- "Binary"
colnames(Column_2) <- "Binary"

# Merge the two columns to get the value that's double the number of pseudo-absences as presences, i.e. based on Liu et al. (2019)
df <- data.frame(a=Column_1, b=Column_2)
Number_of_PAs <- data.frame(Binary = c(df[,"Binary"], df[,"Binary.1"]))
nrow(Number_of_PAs)

############################################### Biodmod formatting ##############################################
myBiomodData <- BIOMOD_FormatingData(resp.var = my.resp.var,
                                       resp.xy = my.resp.xy,
                                       expl.var = clim, 
                                       resp.name = mySpc,
                                       PA.nb.rep = 5,
                                       PA.nb.absences = nrow(Number_of_PAs),
                                       PA.strategy = 'random',
                                       na.rm = TRUE)

myBiomodData # check everything looks okay

# Extract PA1 - PA5 for MESS
## function to get presence and PA datasets
get_PAtab <- function(bfd){
  dplyr::bind_cols(
    x = bfd@coord[, 1],
    y = bfd@coord[, 2],
    status = bfd@data.species,
    bfd@PA
  )
}

# Note that you need to filter by status in excel to see the xy coordinates for the PAs (1 = presence, NA = pseudo-absence)
PA_datasets <- get_PAtab(myBiomodData)
write.csv(PA_datasets, file=paste("pseudoabsences", ".csv",sep=""))


############################################## Run Species Distribution Models ##############################################
myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,models = c('GAM', 'GBM', 'RF', 'MARS', 'MAXENT.Phillips'),
                                      models.options = BIOMOD_ModelingOptions(),
                                      NbRunEval=5, 
                                      DataSplit=70,
                                      Yweights=NULL,
                                      VarImport=5, 
                                      models.eval.meth = c('TSS', 'ROC'),
                                      SaveObj = TRUE,
                                      rescal.all.models = FALSE, 
                                      do.full.models = FALSE)

# --------------------- BIOMOD Ensemble models ------------------------------------- #

# Ensemble model for each algorithm 
myBiomodEM_algo <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut, 
                                             chosen.models = 'all', 
                                             em.by='algo', 
                                             eval.metric = c('TSS', 'ROC'), 
                                             eval.metric.quality.threshold = NULL, 
                                             prob.mean = F, 
                                             prob.cv = T, 
                                             prob.ci = F, 
                                             prob.ci.alpha = 0.05, 
                                             prob.median = F, 
                                             committee.averaging = F, 
                                             prob.mean.weight = T, 
                                             prob.mean.weight.decay = 'proportional')#ensemble model for each algorithm 



# Ensemble model for consensus across algorithms 
myBiomodEM_all <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                            chosen.models = 'all', 
                                            em.by='all', 
                                            eval.metric = c('TSS', 'ROC'), 
                                            eval.metric.quality.threshold = c(0.4, 0.7), 
                                            prob.mean = F, 
                                            prob.cv = T, 
                                            prob.ci = F, 
                                            prob.ci.alpha = 0.05, 
                                            prob.median = F, 
                                            committee.averaging = F, 
                                            prob.mean.weight = T, 
                                            prob.mean.weight.decay = 'proportional')#general ensemble model

 # ---------------------------- Evaluation and variable (principal component) importance ----------------------------------- #

# Variable importance 
myModelsVarImport<-get_variables_importance(myBiomodModelOut) # importance of explanatory variables in each run

#evaluation of the ensemble model
myBiomodEMEval_all<- get_evaluations(myBiomodEM_all) #get evaluation scores for general ensemble model ('all')
#get overall variable importance (mean across all runs for each algorithm and variable)
varimp<-matrix(nrow=dim(myModelsVarImport)[1],ncol=dim(myModelsVarImport)[2],NA)# create an empty matrix
row.names(varimp)=dimnames(myModelsVarImport)[[1]]#assign names to rows
colnames(varimp)=dimnames(myModelsVarImport)[[2]]#assign names to columns                     
for(j in 1:dim(myModelsVarImport)[1]){
  for(k in 1:dim(myModelsVarImport)[2]){
    varimpmean<-mean(myModelsVarImport[j,k,,1])
    varimp[j,k]=varimpmean
  }
} 

myBiomodModelEval <- get_evaluations(myBiomodModelOut) # get evaluation scores for each run
myBiomodEMEval_algo<- get_evaluations(myBiomodEM_algo) #get evaluation scores for each technique across runs

# Write variable/PC importance
write.csv(varimp,file=paste("Variable_Importance", ".csv",sep="")) # write csv variable importance doc

# Write evaluation files
# Ensemble (all)
write.csv(myBiomodEMEval_all, file=paste("Model_Evaluation_Ensemble",".csv",sep="")) # write csv ensemble model evaluation doc 
# Each algorithm in the ensemble (all)
write.csv(myBiomodModelEval, file=paste("Model_Evaluation_all",".csv",sep="")) # write csv all individual model evaluation doc 

# Individual algorithms (algo)
write.csv(myBiomodEMEval_algo, file=paste("Model_Evaluation_algo",".csv",sep="")) # write csv algo model evaluation doc 


################################################ PREDICTIONS ################################################
current_stack <- clim # Changing name of masked current raster stack

# Project ensemble SDM onto current PCA'd climate
myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                          new.env = current_stack, 
                                          proj.name = 'current',
                                          binary.meth = 'TSS', 'ROC')
myCurrentProjection <-  get_predictions(myBiomodProjection)#predictions for each run (maps)

# Ensemble forecast
myBiomodEF_algo <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_algo, 
                                                projection.output = myBiomodProjection,
                                                proj.name = 'algo_current')
myBiomodEF_all <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all, 
                                               projection.output = myBiomodProjection,
                                               proj.name = 'all_current')

myCurrentProjectionEM_algo<-get_predictions(myBiomodEF_algo)#predictions for each algorithm (maps)
myCurrentProjectionEM_all<-get_predictions(myBiomodEF_all)#predictions for ensemble model (maps)

#Write current suitability for ensemble model with all algorithms 
All_finalmap<-myCurrentProjectionEM_all[[2]] # Extract from the ensemble ('all') map
writeRaster(All_finalmap, filename="All_finalmap_Current", format = 'GTiff')

# Write current suitability outputs for each algorithm (algo)
GBM_algo_<-subset(myCurrentProjectionEM_algo,names(myCurrentProjectionEM_algo)[grep("_GBM",names(myCurrentProjectionEM_algo))])#select GBM raster
RF_algo_<-subset(myCurrentProjectionEM_algo,names(myCurrentProjectionEM_algo)[grep("_RF",names(myCurrentProjectionEM_algo))])#select RF raster
MAXENT_algo_<-subset(myCurrentProjectionEM_algo,names(myCurrentProjectionEM_algo)[grep("_MAXENT.Phillips",names(myCurrentProjectionEM_algo))])#select MAXENT raster
MARS_algo_ <-subset(myCurrentProjectionEM_algo,names(myCurrentProjectionEM_algo)[grep("_MARS",names(myCurrentProjectionEM_algo))])
GAM_algo_ <-subset(myCurrentProjectionEM_algo,names(myCurrentProjectionEM_algo)[grep("_GAM",names(myCurrentProjectionEM_algo))])

writeRaster(GBM_algo_, filename='Continuous_Outputs_algo_GBM', format = 'GTiff')
writeRaster(RF_algo_, filename='Continuous_Outputs_algo_RF', format = 'GTiff')
writeRaster(MAXENT_algo_, filename='Continuous_Outputs_algo_MAXENT', format = 'GTiff')
writeRaster(MARS_algo_, filename='Continuous_Outputs_algo_MARS', format = 'GTiff')
writeRaster(GAM_algo_, filename='Continuous_Outputs_algo_GAM', format = 'GTiff')

##################################### Calculate continuous boyce index on ensemble model with ecospat #####################################
# Thin the removed occurrences to match the thinning resolution/structure of the presences
removed_occs <- read.csv("removed_occs.csv", header = TRUE)
Lat <- removed_occs$decimalLatitude
Long <- removed_occs$decimalLongitude
LongLat <- as.data.frame(cbind(Long,Lat))

Thinned_eval <- filterByProximity(LongLat, dist=20, mapUnits = FALSE, returnIndex = FALSE)

CBI_out <- ecospat.boyce(fit = All_finalmap, obs = Thinned_eval, nclass= 0, PEplot = TRUE)

CBI_val <- CBI_out$Spearman.cor
write.csv(CBI_val, file=paste("CBI",".csv",sep=""))

############################################################ PROJECTIONS ##########################################################
############################## Loading climate change projection rasters ##############################
# Make sure names perfectly match with the names of variables used to build the model (even the names of the files)
##############################################################################################
################################## CanESM5 126 ###################################################
# can126 scenario
# load future scenario rasters
PC1 <- raster("PCA_out_bio_can_126_2030/__PC1.tif")
PC2 <- raster("PCA_out_bio_can_126_2030/__PC2.tif")
PC3 <- raster("PCA_out_bio_can_126_2030/__PC3.tif")
PC4 <- raster("PCA_out_bio_can_126_2030/__PC4.tif")
PC5 <- raster("PCA_out_bio_can_126_2030/__PC5.tif")


# stack rasters
can126_stack <- stack(PC1, PC2, PC3, PC4, PC5) # making a raster object
myBiomodProjcan126 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                           new.env = can126_stack, 
                                           proj.name = 'can126_scen',
                                           build.clamping.mask = T,
                                           binary.meth = 'TSS', 'ROC')
can126_scen <-  get_predictions(myBiomodProjcan126) 

# Ensemble forecast
myBiomodEF_can126 <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all, 
                                                   projection.output = myBiomodProjcan126,
                                                   proj.name = 'can126' )
myfuture_can126_Proj_ <- get_predictions(myBiomodEF_can126)#predictions for ensemble model (maps)

################################## CanESM5 370 ###################################################
# can370 scenario
# load future scenario rasters
PC1 <- raster("PCA_out_bio_can_370_2030/__PC1.tif")
PC2 <- raster("PCA_out_bio_can_370_2030/__PC2.tif")
PC3 <- raster("PCA_out_bio_can_370_2030/__PC3.tif")
PC4 <- raster("PCA_out_bio_can_370_2030/__PC4.tif")
PC5 <- raster("PCA_out_bio_can_370_2030/__PC5.tif")


# stack rasters
can370_stack <- stack(PC1, PC2, PC3, PC4, PC5) # making a raster object
myBiomodProjcan370 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = can370_stack, 
                                        proj.name = 'can370_scen',
                                        build.clamping.mask = T,
                                        binary.meth = 'TSS', 'ROC')
can370_scen <-  get_predictions(myBiomodProjcan370) 

# Ensemble forecast
myBiomodEF_can370 <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all, 
                                                projection.output = myBiomodProjcan370,
                                                proj.name = 'can370' )
myfuture_can370_Proj_ <- get_predictions(myBiomodEF_can370)#predictions for ensemble model (maps)

################################## CNRM-CM6-1 126 ###################################################
# cnrm126 scenario
# load future scenario rasters
PC1 <- raster("PCA_out_bio_Cnrm_126_2030/__PC1.tif")
PC2 <- raster("PCA_out_bio_Cnrm_126_2030/__PC2.tif")
PC3 <- raster("PCA_out_bio_Cnrm_126_2030/__PC3.tif")
PC4 <- raster("PCA_out_bio_Cnrm_126_2030/__PC4.tif")
PC5 <- raster("PCA_out_bio_Cnrm_126_2030/__PC5.tif")


# stack rasters
cnrm126_stack <- stack(PC1, PC2, PC3, PC4, PC5) # making a raster object
myBiomodProjcnrm126 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                         new.env = cnrm126_stack, 
                                         proj.name = 'cnrm126_scen',
                                         build.clamping.mask = T,
                                         binary.meth = 'TSS', 'ROC')
cnrm126_scen <-  get_predictions(myBiomodProjcnrm126) 

# Ensemble forecast
myBiomodEF_cnrm126 <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all, 
                                                 projection.output = myBiomodProjcnrm126,
                                                 proj.name = 'cnrm126' )
myfuture_cnrm126_Proj_ <- get_predictions(myBiomodEF_cnrm126)#predictions for ensemble model (maps)

################################## CNRM-CM6-1 370 ###################################################
# cnrm370 scenario
# load future scenario rasters
PC1 <- raster("PCA_out_bio_Cnrm_370_2030/__PC1.tif")
PC2 <- raster("PCA_out_bio_Cnrm_370_2030/__PC2.tif")
PC3 <- raster("PCA_out_bio_Cnrm_370_2030/__PC3.tif")
PC4 <- raster("PCA_out_bio_Cnrm_370_2030/__PC4.tif")
PC5 <- raster("PCA_out_bio_Cnrm_370_2030/__PC5.tif")


# stack rasters
cnrm370_stack <- stack(PC1, PC2, PC3, PC4, PC5) # making a raster object
myBiomodProjcnrm370 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                         new.env = cnrm370_stack, 
                                         proj.name = 'cnrm370_scen',
                                         build.clamping.mask = T,
                                         binary.meth = 'TSS', 'ROC')
cnrm370_scen <-  get_predictions(myBiomodProjcnrm370) 

# Ensemble forecast
myBiomodEF_cnrm370 <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all, 
                                                 projection.output = myBiomodProjcnrm370,
                                                 proj.name = 'cnrm370' )
myfuture_cnrm370_Proj_ <- get_predictions(myBiomodEF_cnrm370)#predictions for ensemble model (maps)

################################## IPSL-CM6A-LR 126 ###################################################
# ipsl126 scenario
# load future scenario rasters
PC1 <- raster("PCA_out_bio_Ipsl_126_2030/__PC1.tif")
PC2 <- raster("PCA_out_bio_Ipsl_126_2030/__PC2.tif")
PC3 <- raster("PCA_out_bio_Ipsl_126_2030/__PC3.tif")
PC4 <- raster("PCA_out_bio_Ipsl_126_2030/__PC4.tif")
PC5 <- raster("PCA_out_bio_Ipsl_126_2030/__PC5.tif")


# stack rasters
ipsl126_stack <- stack(PC1, PC2, PC3, PC4, PC5) # making a raster object
myBiomodProjipsl126 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                         new.env = ipsl126_stack, 
                                         proj.name = 'ipsl126_scen',
                                         build.clamping.mask = T,
                                         binary.meth = 'TSS', 'ROC')
ipsl126_scen <-  get_predictions(myBiomodProjipsl126) 

# Ensemble forecast
myBiomodEF_ipsl126 <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all, 
                                                 projection.output = myBiomodProjipsl126,
                                                 proj.name = 'ipsl126' )
myfuture_ipsl126_Proj_ <- get_predictions(myBiomodEF_ipsl126)#predictions for ensemble model (maps)

################################## IPSL-CM6A-LR 370 ###################################################
# ipsl370 scenario
# load future scenario rasters
PC1 <- raster("PCA_out_bio_Ipsl_370_2030/__PC1.tif")
PC2 <- raster("PCA_out_bio_Ipsl_370_2030/__PC2.tif")
PC3 <- raster("PCA_out_bio_Ipsl_370_2030/__PC3.tif")
PC4 <- raster("PCA_out_bio_Ipsl_370_2030/__PC4.tif")
PC5 <- raster("PCA_out_bio_Ipsl_370_2030/__PC5.tif")


# stack rasters
ipsl370_stack <- stack(PC1, PC2, PC3, PC4, PC5) # making a raster object
myBiomodProjipsl370 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                         new.env = ipsl370_stack, 
                                         proj.name = 'ipsl370_scen',
                                         build.clamping.mask = T,
                                         binary.meth = 'TSS', 'ROC')
ipsl370_scen <-  get_predictions(myBiomodProjipsl370) 

# Ensemble forecast
myBiomodEF_ipsl370 <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all, 
                                                 projection.output = myBiomodProjipsl370,
                                                 proj.name = 'ipsl370' )
myfuture_ipsl370_Proj_ <- get_predictions(myBiomodEF_ipsl370)#predictions for ensemble model (maps)

################################## MIROC6 126 ###################################################
# miroc126 scenario
# load future scenario rasters
PC1 <- raster("PCA_out_bio_Miroc_126_2030/__PC1.tif")
PC2 <- raster("PCA_out_bio_Miroc_126_2030/__PC2.tif")
PC3 <- raster("PCA_out_bio_Miroc_126_2030/__PC3.tif")
PC4 <- raster("PCA_out_bio_Miroc_126_2030/__PC4.tif")
PC5 <- raster("PCA_out_bio_Miroc_126_2030/__PC5.tif")


# stack rasters
miroc126_stack <- stack(PC1, PC2, PC3, PC4, PC5) # making a raster object
myBiomodProjmiroc126 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                          new.env = miroc126_stack, 
                                          proj.name = 'miroc126_scen',
                                          build.clamping.mask = T,
                                          binary.meth = 'TSS', 'ROC')
miroc126_scen <-  get_predictions(myBiomodProjmiroc126) 

# Ensemble forecast
myBiomodEF_miroc126 <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all, 
                                                  projection.output = myBiomodProjmiroc126,
                                                  proj.name = 'miroc126' )
myfuture_miroc126_Proj_ <- get_predictions(myBiomodEF_miroc126)#predictions for ensemble model (maps)

################################## MIROC6 370 ###################################################
# miroc370 scenario
# load future scenario rasters
PC1 <- raster("PCA_out_bio_Miroc_370_2030/__PC1.tif")
PC2 <- raster("PCA_out_bio_Miroc_370_2030/__PC2.tif")
PC3 <- raster("PCA_out_bio_Miroc_370_2030/__PC3.tif")
PC4 <- raster("PCA_out_bio_Miroc_370_2030/__PC4.tif")
PC5 <- raster("PCA_out_bio_Miroc_370_2030/__PC5.tif")


# stack rasters
miroc370_stack <- stack(PC1, PC2, PC3, PC4, PC5) # making a raster object
myBiomodProjmiroc370 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                          new.env = miroc370_stack, 
                                          proj.name = 'miroc370_scen',
                                          build.clamping.mask = T,
                                          binary.meth = 'TSS', 'ROC')
miroc370_scen <-  get_predictions(myBiomodProjmiroc370) 

# Ensemble forecast
myBiomodEF_miroc370 <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all, 
                                                  projection.output = myBiomodProjmiroc370,
                                                  proj.name = 'miroc370' )
myfuture_miroc370_Proj_ <- get_predictions(myBiomodEF_miroc370)#predictions for ensemble model (maps)

################################## MRI-ESM2-0 126 ###################################################
# mri126 scenario
# load future scenario rasters
PC1 <- raster("PCA_out_bio_Mri_126_2030/__PC1.tif")
PC2 <- raster("PCA_out_bio_Mri_126_2030/__PC2.tif")
PC3 <- raster("PCA_out_bio_Mri_126_2030/__PC3.tif")
PC4 <- raster("PCA_out_bio_Mri_126_2030/__PC4.tif")
PC5 <- raster("PCA_out_bio_Mri_126_2030/__PC5.tif")


# stack rasters
mri126_stack <- stack(PC1, PC2, PC3, PC4, PC5) # making a raster object
myBiomodProjmri126 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = mri126_stack, 
                                        proj.name = 'mri126_scen',
                                        build.clamping.mask = T,
                                        binary.meth = 'TSS', 'ROC')
mri126_scen <-  get_predictions(myBiomodProjmri126) 

# Ensemble forecast
myBiomodEF_mri126 <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all, 
                                                projection.output = myBiomodProjmri126,
                                                proj.name = 'mri126' )
myfuture_mri126_Proj_ <- get_predictions(myBiomodEF_mri126)#predictions for ensemble model (maps)

################################## MRI-ESM2-0 370 ###################################################
# mri370 scenario
# load future scenario rasters
PC1 <- raster("PCA_out_bio_Mri_370_2030/__PC1.tif")
PC2 <- raster("PCA_out_bio_Mri_370_2030/__PC2.tif")
PC3 <- raster("PCA_out_bio_Mri_370_2030/__PC3.tif")
PC4 <- raster("PCA_out_bio_Mri_370_2030/__PC4.tif")
PC5 <- raster("PCA_out_bio_Mri_370_2030/__PC5.tif")


# stack rasters
mri370_stack <- stack(PC1, PC2, PC3, PC4, PC5) # making a raster object
myBiomodProjmri370 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = mri370_stack, 
                                        proj.name = 'mri370_scen',
                                        build.clamping.mask = T,
                                        binary.meth = 'TSS', 'ROC')
mri370_scen <-  get_predictions(myBiomodProjmri370) 

# Ensemble forecast
myBiomodEF_mri370 <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all, 
                                                projection.output = myBiomodProjmri370,
                                                proj.name = 'mri370' )
myfuture_mri370_Proj_ <- get_predictions(myBiomodEF_mri370)#predictions for ensemble model (maps)