set.seed(123)
rm(list=ls())
library(SpatialEpi)
library(spdep)
library(INLA)
library(RColorBrewer)
library(maptools)
library(maps)
library(ggplot2)
library(sp)
library(lattice)
library(surveillance)
library(parallel)

## Get Pennsylvania map file
setwd("~/Documents/Classes/spatial_epi/project/")
load("../homework/hw3/USA_adm2.RData") # From http://gadm.org/

# define the 50 states remove DC
states <- unique(gadm$NAME_1)
states <- as.character(states[!(states %in% c("District of Columbia"))])

# each data unit will be a state with different data
state_data <- lapply(states, function(x) list(name=x))
names(state_data) <- states

# while testing lets only keep three states
state_data <- state_data[c("California", "Washington")]

apply2 <- function(f, new_key, state_data, cores=1){
    # apply a function in parallel to an old key of state data and then
    # add the newely generated data to the key list
    states <- names(state_data)
    temp <- mclapply(states, function(x) 
        f(state_data[[x]]), mc.cores=cores)
    names(temp) <- states
    for(state in states){
        state_data[[state]][[new_key]] <- temp[[state]]
    }
    state_data
}

state_data <- apply2(function(x) as(gadm[which(gadm$NAME_1==x$name),], "SpatialPolygons"), 
                     "SpatialPolygons", state_data, cores=3)

state_data <- apply2(function(x) poly2adjmat(x$SpatialPolygons), 
                     "adjmat", state_data, cores=3)

distance_matrix <- function(tmat){
    loops <- 0
    islands <- names(colSums(tmat)[colSums(tmat) == 0])
    sans_island <- tmat
    sans_island[,] <- TRUE
    sans_island[islands,] <- FALSE
    sans_island[,islands] <- FALSE
    while(sum(tmat[lower.tri(tmat, diag = FALSE) & sans_island] == 0) != 0){
        loops <- loops + 1
        print(paste0("Running through loop number ", loops))
        for(i in 1:(nrow(tmat) - 1)){
            for(j in (i+1):(nrow(tmat))){
                if(tmat[i,j] == 0){
                    i_neighbors <- tmat[i,][tmat[i,] != 0]
                    j_neighbors <- tmat[j,][tmat[j,] != 0]
                    shared <- intersect(names(i_neighbors), names(j_neighbors))
                    min_ <- nrow(tmat)
                    for(s in shared){
                        distance <- j_neighbors[s] + i_neighbors[s]
                        if(distance < min_){
                            min_ <- distance
                        }
                    }
                    if(length(s) > 0 & min_ <= (loops + 1)){
                        tmat[i,j] <- min_
                        tmat[j,i] <- min_
                    }
                }
            }
        }
    }
    tmat[islands,] <- nrow(tmat) * 3
    tmat[,islands] <- nrow(tmat) * 3
    tmat[diag(tmat)] <- 0
    return (tmat)
}

rmvn <- function(n, mu = 0, V = matrix(1)) {
    # http://rstudio-pubs-static.s3.amazonaws.com/9688_a49c681fab974bbca889e3eae9fbb837.html
    p <- length(mu)
    if (any(is.na(match(dim(V), p)))) 
        stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

simulate_sre <- function(dist_mat){
    N <- nrow(dist_mat)
    phi <- 0.38
    rmvn(1, rep(0, N), exp(-phi * dist_mat))
} 

delta_simulatar <- function(unit){
    # Given a spatial polygons data set applies
    # a delta to each unit
    sp <- unit$SpatialPolygons
    dist_mat <- distance_matrix(unit$adjmat)
    N <- length(sp) # number of counties in a state unit
    sp$ID <- 1:N
    
    datur_od <- sapply(1:100, function(x) rnorm(length(sp), 0, 1))
    datur_se <- sapply(1:100, function(x) simulate_sre(dist_mat))
    x1 <- sapply(1:100, function(x)rnorm(length(sp)))
    x2 <- sapply(1:100, function(x)rnorm(length(sp)))
    beta1 <- 2
    beta2 <- -3
    mean_ <- beta1 * x1 + beta2 * x2
    random <- sapply(1:ncol(datur_od), function(x) 
        rpois(length(datur_od[,x]), exp(datur_od[,x])))
    spatial <- sapply(1:ncol(datur_se), function(x) 
        rpois(length(datur_se[,x]), exp(datur_se[,x])))
    both <- sapply(1:ncol(datur_se), function(x) 
        rpois(length(datur_se[,x]), exp(datur_se[,x] + datur_od[,x])))
    random_covs <- sapply(1:ncol(datur_od), function(x) 
        rpois(length(datur_od[,x]), exp(mean_[,x] + datur_od[,x])))
    spatial_covs <- sapply(1:ncol(datur_se), function(x) 
        rpois(length(datur_se[,x]), exp(mean_[,x] + datur_se[,x])))
    both_covs <- sapply(1:ncol(datur_se), function(x) 
        rpois(length(datur_se[,x]), exp(mean_[,x] + datur_se[,x] + datur_od[,x])))
    
    sp@data[,paste0("over_dispersion", 1:100)] <- datur_od
    sp@data[,paste0("spatial_random", 1:100)] <- datur_se
    sp@data[,paste0("x1", 1:100)] <- x1
    sp@data[,paste0("x2", 1:100)] <- x2
    sp@data[,paste0("random", 1:100)] <- random
    sp@data[,paste0("spatial", 1:100)] <- spatial
    sp@data[,paste0("both", 1:100)] <- both
    sp@data[,paste0("random_covs", 1:100)] <- random_covs
    sp@data[,paste0("spatial_covs", 1:100)] <- spatial_covs
    sp@data[,paste0("both_covs", 1:100)] <- both_covs

    sp
}
                                               
state_data <- apply2(delta_simulatar,"SpatialPolygons", state_data)
par(mfcol=c(1,2))
p1 <- spplot(state_data$California$SpatialPolygons, "spatial100", 
             main="Simulated Spatial Random Effects")
p2 <- spplot(state_data$California$SpatialPolygons, "random100",
             main="Simulated Overdispersion Random Effects")
p3 <- spplot(state_data$Washington$SpatialPolygons, "spatial100", 
             main="Simulated Spatial Random Effects")
p4 <- spplot(state_data$Washington$SpatialPolygons, "random100",
             main="Simulated Overdispersion Random Effects")

print(p1, position=c(0, .5, .5, 1), more=T)
print(p2, position=c(.5, .5, 1, 1), more=T)
print(p3, position=c(0, 0, .5, .5), more=T)
print(p4, position=c(.5, 0, 1, .5))

# specify the latent structure using a formula object
formula <- y ~ 1 + f ( region ,
                      model = " bym2 " ,
                      graph =g ,
                      scale.model = TRUE ,
                      constr = TRUE ,
                      hyper = list (
                          phi = list (
                              prior = "pc" ,
                              param = c (0.5 , 2 / 3) ,
                              initial = -3) ,
                          prec = list (
                              prior = " pc.prec " ,
                              param = c (0.2 / 0.31 , 0.01) ,
                              initial = 5) ) )
# call the inla function
result = inla ( formula , data = sardinia , family = " poisson " , E =E ,
                control.predictor = list ( compute = TRUE ) )
# get improved estimates for the h y p e r p a r a m e t e r s
result = inla.hyperpar (r , dz = 0.2 , diff . logdens =20)