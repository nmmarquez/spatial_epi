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
library(MASS)

phi <- seq(0.01,1,.01)

## Get Pennsylvania map file
setwd("~/Documents/Classes/spatial_epi/project/")
load("../homework/hw3/USA_adm2.RData") # From http://gadm.org/

cont_usa_locs <- c("Texas", "Louisiana", "Arkansas", "Oklahoma")

cont_usa <- gadm[(gadm@data$NAME_1 %in% cont_usa_locs),]

# define the 50 states remove DC
states <- unique(gadm$NAME_1)
states <- as.character(states[!(states %in% c("District of Columbia"))])

# each data unit will be a state with different data
state_data <- lapply(states, function(x) list(name=x))
names(state_data) <- states

# while testing lets only keep three states
state_data <- state_data[c("California", "Tennessee")]

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

# state_data$cont_usa <- list(name="cont_usa", "SpatialPolygons" = cont_usa)

state_data <- apply2(function(x) poly2adjmat(x$SpatialPolygons), 
                     "adjmat", state_data, cores=3)

simulate_sre <- function(Q){
    n <- nrow(Q)
    Qs <- inla.scale.model(Q, constr=list(A=matrix(1, nrow=1, ncol=n), e=0))
    Q_star <- ginv(as.matrix(Qs))
    t(mvrnorm(n=100, 0 * 1:n, Q_star))
}

delta_simulatar <- function(unit){
    # Given a spatial polygons data set applies
    # a delta to each unit
    sp <- unit$SpatialPolygons
    N <- length(sp) # number of counties in a state unit
    sp$ID <- 1:N
    
    n_delta_i <- rowSums(unit$adjmat)
    Q <- unit$adjmat * -1
    diag(Q) <- n_delta_i
    
    v <- sapply(1:100, function(x) rnorm(length(sp), 0, 1))
    u <- sapply(1:100, function(x) simulate_sre(Q))
    x1 <- sapply(1:100, function(x)rnorm(length(sp)))
    x2 <- sapply(1:100, function(x)rnorm(length(sp)))
    beta1 <- 2
    beta2 <- -3
    prec <- 1/sqrt(.5)
    mean_ <- beta1 * x1 + beta2 * x2
    random <- sapply(1:ncol(v), function(x) 
        rpois(length(v[,x]), exp(prec * v[,x])))
    spatial <- sapply(1:ncol(u), function(x) 
        rpois(length(u[,x]), exp(prec * u[,x])))
    both <- sapply(1:ncol(u), function(x) 
        rpois(length(u[,x]), 
              exp(prec * (sqrt(1 - phi[x]) * v[,x] + sqrt(phi[x]) * u[,x]))))
    random_covs <- sapply(1:ncol(v), function(x) 
        rpois(length(v[,x]), exp(mean_[,x] + prec * v[,x])))
    spatial_covs <- sapply(1:ncol(u), function(x) 
        rpois(length(u[,x]), exp(mean_[,x] + prec * u[,x])))
    both_covs <- sapply(1:ncol(u), function(x) 
        rpois(length(u[,x]), 
              exp(mean_[,x] + prec * (sqrt(1 - phi[x]) * v[,x] + sqrt(phi[x]) * u[,x]))))
    
    sp@data[,paste0("v", 1:100)] <- v
    sp@data[,paste0("u", 1:100)] <- u
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

plot_ca_tn <- function(sim_num){
    ca_datur <- state_data$California$SpatialPolygons
    tn_datur <- state_data$Tennessee$SpatialPolygons
    spat_text <- "Simulated Spatial Random Effects"
    od_text <- "Simulated Overdispersion Random Effects"
    p1 <- spplot(ca_datur, paste0("spatial", sim_num), main=spat_text)
    p2 <- spplot(ca_datur, paste0("random", sim_num), main=od_text)
    p3 <- spplot(tn_datur, paste0("spatial", sim_num), main=spat_text)
    p4 <- spplot(tn_datur, paste0("random", sim_num), main=od_text)
    print(p1, position=c(0, .5, .5, 1), more=T)
    print(p2, position=c(.5, .5, 1, 1), more=T)
    print(p3, position=c(0, 0, .5, .5), more=T)
    print(p4, position=c(.5, 0, 1, .5))
}

plot_ca_tn(100)

df <- state_data$Tennessee$SpatialPolygons@data
adj <- state_data$Tennessee$adjmat 

summary(inla(spatial_covs100 ~ x1100 + x2100 + 
                 f(ID, model = "iid", param = c(1, 0.014)), 
             data=df, family="poisson"))

summary(inla(spatial_covs100 ~ x1100 + x2100 +
                 f(ID, model="besag", graph=adj,
                   param = c(1, 0.68)), 
             data=df, family="poisson"))

summary(inla(spatial_covs1 ~ x11 + x21 + f(ID, model="bym2", graph=adj), 
             data=df, family="poisson"))

summary(inla(spatial100 ~ 1 + f(ID, model="bym2", graph=adj), 
             data=df, family="poisson"))

# specify the latent structure using a formula object
formula <- both100 ~ 1 + 
    f(ID, model = "bym2", scale.model=T, constr=T, graph=adj, 
      hyper=list(phi=list(prior ="pc", param=c(0.5, 2/3), initial=-3),
                 prec=list(prior="pc.prec", initial=5, param=c(0.2/0.31, 0.01))))

# call the inla function
result <- inla(formula , data=df,
               family="poisson", control.predictor=list(compute=T))
# get improved estimates for the h y p e r p a r a m e t e r s
result <- inla.hyperpar(result , dz=0.2 , diff.logdens=20)
summary(result)
