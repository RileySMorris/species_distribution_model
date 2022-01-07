#Author: Riley S. Morris
#Date: May 10, 2021
#Latest Revision: Sep 30, 2021
#Description: Species distribution model using various modelling methodologies.
#----

#Data Resources-------------
# Raster data: https://www.gis-blog.com/r-raster-data-acquisition/ - Catalogue of raster data imports using raster package and getData() function
# VertNet: http://vertnet.org/ - Species occurrence data
# VertNet download package: https://cran.r-project.org/web/packages/rvertnet/rvertnet.pdf - Package for accessing VertNet data
# Worldclim: https://www.worldclim.org/data/bioclim.html - Bioclimatic data


#Script Improvement Plan-------------
# Add rvertnet data access capabilities
# Add user data input section
# Add output report with data diagnostics results
# Develop own historical bioclim data procedure
# Develop projection bioclim data procedure
# Improve mapping aesthetic
# Add BIOMOD2 SDM analysis for alternative approach or comparison methodology
# Add capability to generate MaxEnt background data for alternative apporach or comparison methodology

#https://damariszurell.github.io/SDM-Intro/#1_Background
#https://www.seascapemodels.org/SDM-fish-course-notes/

#Species Distribution Algorithms ----------------
# BioClimatic envelope (e.g. BioClim)
# Ordinary regression (e.g. incl. in Arc-SDM)
# Generalized additive models (GAM)
# Generalized linear models (GLM)
# Ordination (e.g. Canonical Correspondence Analysis, CCA)
# Classification and regression trees (CART, randomForest)
# Genetic algorithm (e.g. GARP)
# Local weighted regression (Loess)
# Boosting (generalized boosted models)
# Bagging average (randomForest)
# Artificial neural networks
# Bayseian (e.g. WinBUGS)
# Discriminant analysis
# Maxent (point-process model)


#Species Data Import----
  #Load packages
    install.packages(c("maptools","rgeos","raster","rgdal"))
    library(maptools)
    library(rgeos)
    library(raster)
    library(rgdal)

  #Data Download
    occurrance <- read.delim("P:/1609 - Environmental Services/MorrisR/Scripts/sample_data/blanding_turtle_occurrence.txt")
  #Quick plot observations
    data("wrld_simpl")
    plot(wrld_simpl, xlim=c(-110,-60),ylim=c(35,50),axes=T)
    points(occurrance$decimallongitude,occurrance$decimallatitude,col="red",pch=0.75)
  #Clean data
    occurrance <- occurrance[complete.cases(occurrance$decimallongitude)]  #removes obs with missing longitudes
    occurrance <- occurrance[complete.cases(occurrance$decimallatitude)]  #removes obs with missing latitudes
    occurrance <- occurrance[complete.cases(occurrance$country)]  #removes obs with no location data
    dups <- duplicated(occurrance[,c("decimallatitude","decimallongitude")])  #identifies duplicate obs (i.e. with same lat/long)
    sum(dups)  #number of duplicates
    occurrance <- occurrance[!dups,]  #remvoves obs with duplicate lat/longs
    occurrance <- occurrance[occurrance$isfossil==0]  #removes fossil record obs
    occurrance <- occurrance[occurrance$wascaptive==0]  #removes captive record obs
    plot(wrld_simpl, xlim=c(-110,-60),ylim=c(35,50),axes=T)
    points(occurrance$decimallongitude,occurrance$decimallatitude,col="red",pch=0.75)
  #Align coordinate systems among dataset (WGS84 is standard for most GPS systems)
    #Coordinate system catagorization
      Occur.NAD27a <- occurrance[occurrance$geodeticdatum=="NAD27",]
      Occur.NAD27b <- occurrance[occurrance$geodeticdatum=="North American Datum 1927",]
      Occur.NAD27 <- rbind(Occur.NAD27a,Occur.NAD27b)
      Occur.NAD83a <- occurrance[occurrance$geodeticdatum=="NAD83",]
      Occur.NAD83b <- occurrance[occurrance$geodeticdatum=="North American Datum 1983",]
      Occur.NAD83 <- rbind(Occur.NAD83a,Occur.NAD83b)
      Occur.WGS84a <- occurrance[occurrance$geodeticdatum=="WGS84",]
      Occur.WGS84b <- occurrance[occurrance$geodeticdatum=="World Gendetic System 1984",]
      Occur.nr <- occurrance[occurrance$geodeticdatum=="not recorded",]
      Occur.nr83 <- occurrance[occurrance$geodeticdatum=="not recorded (forced WGS84)",]
      Occur.unk <- occurrance[occurrance$geodeticdatum=="unknown",]
      Occur.WGS84 <- rbind(Occur.WGS84a,Occur.WGS84b,Occur.nr,Occur.nr83,Occur.unk)
    #Coorindate system covnersion to align to WGS84
      coordinates(Occur.NAD27) <- ~decimallongitude+decimallatitude
      projection(Occur.NAD27) <- crs('+proj=longlat +datum=NAD27')
      NAD27.converted <- spTransform(Occur.NAD27, crs('+proj=longlat +datum=WGS84'))
      coordinates(Occur.NAD83) <- ~decimallongitude+decimallatitude
      projection(Occur.NAD83) <- crs('+proj=longlat +datum=NAD83')
      NAD83.converted <- spTransform(Occur.NAD83, crs('+proj=longlat +datum=WGS84'))
      coordinates(Occur.WGS84) <- ~decimallongitude+decimallatitude
      projection(Occur.WGS84) <- crs('+proj=longlat +datum=WGS84')
      occur <- rbind(NAD27.converted,NAD83.converted,Occur.WGS84)
  #Plot points with aligned coordinate system
    plot(wrld_simpl, xlim=c(-110,-60),ylim=c(35,50),axes=T)
    points(occur$decimallongitude,occur$decimallatitude,col="red",pch=0.75)

#Climate Data
  #WorldClim Method
    install.packages("dismo")
    library(dismo)
    climate <- getData('worldclim',download=T,var='bio',res=2.5)
    plot(climate)
  #Climdex Method
  
  #Projections
    
###BIOCLIM APPROACH####----
#Load packages
install.packages("dismo")
library(dismo)

#Ecosystem modelling
obsclim <- extract(climate, occur)  #extracts bioclimatic data for each obs point
bioclim.model <- bioclim(obsclim)
pairs(bioclim.model,pa='p')
predictors <- stack(climate$bio1,climate$bio2,climate$bio3,
                    climate$bio4,climate$bio5,climate$bio6,
                    climate$bio7,climate$bio8,climate$bio9,
                    climate$bio10,climate$bio11,climate$bio12,
                    climate$bio13,climate$bio14,climate$bio15,
                    climate$bio16,climate$bio17,climate$bio18,
                    climate$bio19)
predictions <- predict(predictors,bioclim.model)
plot(predictions,xlim=c(-110,-60),ylim=c(35,50),axes=T)

###BIOMOD2 APPROACH###----------------
  #Reference:
    #https://www.youtube.com/watch?v=-IAdf8Vh6uY
    #https://www.youtube.com/watch?v=QrwqhJgRbnY&t=189s
    #https://www.youtube.com/watch?v=ChxIdJBLXE0
#Load packages
devtools::install_github("biomodhub/biomod2")
install.packages(c("biomod2","ggplot2","gridExtra","raster"))
library(c(biomod2,ggplot2,gridExtra,raster))

#Import data
occurrance <- read.delim("P:/1609 - Environmental Services/MorrisR/Scripts/sample_data/blanding_turtle_occurrence.txt")
summary(occurance)