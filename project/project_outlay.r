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

M <- 150 # number of simulations
phi <- c(rep(.25, M/3), rep(.5, M/3), rep(.75, M/3))
tau <- 4
true_betas <- c(1, -1)

## Get Pennsylvania map file
setwd("~/Documents/Classes/spatial_epi/project/")
load("../homework/hw3/USA_adm2.RData") # From http://gadm.org/

cont_usa_locs <- c("Texas", "Louisiana")

cont_usa <- gadm[(gadm@data$NAME_1 %in% cont_usa_locs),]

cont_usa <- list(name="cont_usa", "SpatialPolygons" = cont_usa)

cont_usa$adjmat <- poly2adjmat(cont_usa$SpatialPolygons)

simulate_sre <- function(Q, M){
    n <- nrow(Q)
    Qs <- inla.scale.model(Q, constr=list(A=matrix(1, nrow=1, ncol=n), e=0))
    Q_star <- ginv(as.matrix(Qs))
    t(mvrnorm(n=M, 0 * 1:n, Q_star))
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
    print("qflipped")
    
    v <- sapply(1:M, function(x) rnorm(length(sp), 0, 1))
    print("v created")
    u <- simulate_sre(Q, M=M)
    print("u created")
    
    x1 <- sapply(1:M, function(x)rnorm(length(sp)))
    x2 <- sapply(1:M, function(x)rnorm(length(sp)))
    print("x created")
    beta1 <- true_betas[1]
    beta2 <- true_betas[2]
    prec <- 1/sqrt(tau)
    mean_ <- beta1 * x1 + beta2 * x2
    random <- sapply(1:ncol(v), function(x) 
        rpois(length(v[,x]), exp(prec * v[,x])))
    spatial <- sapply(1:ncol(u), function(x) 
        rpois(length(u[,x]), exp(prec * u[,x])))
    both <- sapply(1:ncol(u), function(x) 
        rpois(length(u[,x]), 
              exp(prec * (sqrt(1 - phi[x]) * v[,x] + sqrt(phi[x]) * u[,x]))))
    print("constant created")
    random_covs <- sapply(1:ncol(v), function(x) 
        rpois(length(v[,x]), exp(mean_[,x] + prec * v[,x])))
    spatial_covs <- sapply(1:ncol(u), function(x) 
        rpois(length(u[,x]), exp(mean_[,x] + prec * u[,x])))
    both_covs <- sapply(1:ncol(u), function(x) 
        rpois(length(u[,x]), 
              exp(mean_[,x] + prec * (sqrt(1 - phi[x]) * v[,x] + sqrt(phi[x]) * u[,x]))))
    print("cov adjusted created")
    
    sp@data[,paste0("v", 1:M)] <- v
    sp@data[,paste0("u", 1:M)] <- u
    sp@data[,paste0("x1", 1:M)] <- x1
    sp@data[,paste0("x2", 1:M)] <- x2
    sp@data[,paste0("random", 1:M)] <- random
    sp@data[,paste0("spatial", 1:M)] <- spatial
    sp@data[,paste0("both", 1:M)] <- both
    sp@data[,paste0("random_covs", 1:M)] <- random_covs
    sp@data[,paste0("spatial_covs", 1:M)] <- spatial_covs
    sp@data[,paste0("both_covs", 1:M)] <- both_covs
    
    sp
}

run_bym2 <- function(i, covs, pc, model=F){
    base <- 'f(ID, model = "bym2", scale.model=T, constr=T, graph=cont_usa$adjmat '
    phi_str <- 'phi=list(prior ="pc", param=c(0.5, 2/3), initial=-3), \n\t\t'
    prec_str <- 'prec=list(prior="pc.prec", initial=5, param=c(0.2/0.31, 0.01))'
    hyper <- paste0(',\n hyper=list(', phi_str, prec_str, ")")
    if(pc){
        re <- paste0(base, hyper, ")")
    }
    else{
        re <- paste0(base, ")")
    }
    if(covs){
        mean_ <- paste0("both_covs", i, " ~ x1", i," + x2", i," + ")
    }
    else{
        mean_ <- paste0("both", i, "~ 1 + ")
    }
    form <- formula(paste0(mean_, re))
    result <- inla(form , data=cont_usa$SpatialPolygons@data,
                   family="poisson", control.predictor=list(compute=T))
    if (model){
        return(result)
    }
    beta_hat <- result$summary.fixed[,"0.5quant"]
    beta_sd <- result$summary.fixed[,"sd"]
    if (length(beta_hat) < 3){
        beta_hat <- c(beta_hat, NA, NA)
        beta_sd <- c(beta_sd, NA, NA)
    }
    betas_ <- ((beta_hat - c(0, true_betas))**2)**.5
    tau_hat <- result$summary.hyperpar["Precision for ID", "0.5quant"]
    tau_ <- ((tau_hat - tau)**2)**.5
    tau_sd <- result$summary.hyperpar["Precision for ID", "sd"]
    phi_hat <- result$summary.hyperpar["Phi for ID","0.5quant"]
    phi_ <- ((phi_hat- phi[i])**2)**.5
    phi_sd <- result$summary.hyperpar["Phi for ID","sd"]
    c(beta_hat, tau_hat, phi_hat, betas_, tau_, phi_, beta_sd, tau_sd, phi_sd)
}

run_bym2_safe <- function(i, covs, pc){
    out <- tryCatch({run_bym2(i, covs, pc)}, 
                    error=function(cond){rep(NA, 15)})
    out
}

cont_usa$SpatialPolygons <- delta_simulatar(cont_usa)

timestamp()
no_covs_no_pc <- mclapply(1:M, function(i) run_bym2_safe(i, F, F), mc.cores=4)
no_covs_no_pc <- do.call(rbind, no_covs_no_pc)
print("25% complete")
timestamp()
no_covs_yes_pc <- mclapply(1:M, function(i) run_bym2_safe(i, F, T), mc.cores=4)
no_covs_yes_pc <- do.call(rbind, no_covs_yes_pc)
print("50% complete")
timestamp()
yes_covs_no_pc <- mclapply(1:M, function(i) run_bym2_safe(i, T, F), mc.cores=4)
yes_covs_no_pc <- do.call(rbind, yes_covs_no_pc)
print("75% complete")
timestamp()
yes_covs_yes_pc <- mclapply(1:M, function(i) run_bym2_safe(i, T, T), mc.cores=4)
yes_covs_yes_pc <- do.call(rbind, yes_covs_yes_pc)
print("100% complete")
timestamp()

f <- "./results.Rdata"

save(no_covs_no_pc, no_covs_yes_pc, yes_covs_no_pc, yes_covs_yes_pc, file=f)
load(f)

sim_num=10
spat_text="spatial"
od_text="overdispersion"
spplot(cont_usa$SpatialPolygons, paste0("u", sim_num), main=spat_text)
spplot(cont_usa$SpatialPolygons, paste0("v", sim_num), main=od_text)

results <- list(no_covs_no_pc, no_covs_yes_pc, yes_covs_no_pc, yes_covs_yes_pc)

results_df <- lapply(results, as.data.frame) 
names(results_df) <- c("no_covs_no_pc", "no_covs_yes_pc", 
                       "yes_covs_no_pc", "yes_covs_yes_pc")

for(i in 1:4){
    names(results_df[[i]]) <- c("beta0_hat", "beta1_hat", "beta2_hat", 
                                "tau_hat", "phi_hat",
                                "beta0_rse", "beta1_rse", "beta2_rse",
                                "tau_rse", "phi_rse", 
                                "beta0_sd", "beta1_sd", "beta2_sd", 
                                "tau_sd", "phi_sd")
    results_df[[i]]$phi <- phi
    if(grepl("no_covs", names(results_df)[i])){
        results_df[[i]]$covs <- F
    }
    else{
        results_df[[i]]$covs <- T
    }
    if(grepl("no_pc", names(results_df)[i])){
        results_df[[i]]$pc <- F
    }
    else{
        results_df[[i]]$pc <- T
    }
}

df <- do.call(rbind, results_df)

ggplot(data=subset(df, !covs), aes(x=log(phi_rse), fill=pc)) + 
    geom_density(alpha=.2)
ggplot(data=subset(df, !covs), aes(y=phi_rse, x=factor(phi))) + 
    geom_boxplot(aes(fill = factor(pc)))

ggplot(data=subset(df, covs), aes(x=log(phi_rse), fill=pc)) + 
    geom_density(alpha=.2)
ggplot(data=subset(df, covs), aes(y=phi_rse, x=factor(phi))) + 
    geom_boxplot(aes(fill = factor(pc)))

ggplot(data=subset(df, covs), aes(y=beta1_rse, x=factor(phi))) + 
    geom_boxplot(aes(fill = factor(pc)))

ggplot(data=subset(df, covs), aes(y=beta2_rse, x=factor(phi))) + 
    geom_boxplot(aes(fill = factor(pc)))
