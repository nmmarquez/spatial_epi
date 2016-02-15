rm(list=ls())
library(SpatialEpi)
library(spdep)
library(INLA)
library(RColorBrewer)
penn.dat <- read.table("penn-data-3.txt", header=T, sep=",")

penn.dat$Expected <- penn.dat$pop*(sum(penn.dat$cases)/sum(penn.dat$pop))
penn.dat$Observed <- penn.dat$cases
penn.dat$Population <- penn.dat$pop


## Get Pennsylvania map file
load("USA_adm2.RData") # From http://gadm.org/
penn <- gadm[which(gadm$NAME_1=="Pennsylvania"), ]

## Make a polygon with the same coordinates as the penn_data has 
tmp <- as(penn, "SpatialPolygons")
penn2  <- latlong2grid(tmp)

# Read in nb file.
nb <- read.gal("penn_nb")
