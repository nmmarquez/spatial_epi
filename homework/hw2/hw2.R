rm(list=ls())
library(SpatialEpi)
library(maps)
library(sp)
library(RColorBrewer)
library(plyr)
library(INLA)

source("~/Documents/Classes/spatial_epi/homework/hw1/ohio.dat")

strata <- ddply(ohio, ~age+sex+race, summarise, 
                deaths=sum(deaths), popn=sum(popn))

# the probability of death from any strata is approximately deaths / population
strata <- ddply(ohio, ~age+sex+race, summarise, 
                deaths=sum(deaths), popn=sum(popn))
strata$prob <- strata$deaths / strata$popn

# merge on the strata data to the original ohio data
ohio <- merge(ohio, strata[,c("age", "sex", "race", "prob")])

# expected number is prob for strata times population number
ohio$expected <- ohio$prob * ohio$popn

# group by county the sum of the expected number of deaths 
# and the observed number of deaths
county <- ddply(ohio, ~fips, summarise, expected=sum(expected), 
                deaths=sum(deaths))
county <- county[order(county$fips),]
county$SMR <- county$deaths / county$expected
county$SMR_error <- (county$SMR / county$expected)**.5

plot_ohio <- function(plotvar, title_, nclr=8, brks=NULL){
    # plots a map of Ohio counties assuimg the plotvar argument is of
    # length 88 where each value represents an Ohio county sorted by 
    # fips number.
    if(length(plotvar) != 88){
        stop("'plotvar' argument must be of length 88")
    }
    # next few lines set up the color scheme for plotting
    plotclr <- brewer.pal(nclr,"BuPu")
    if (is.null(brks)){
        brks <- round(quantile(plotvar,probs=seq(0,1,1/(nclr))),digits=1)
    }
    colornum <- findInterval(plotvar,brks,all.inside=T)
    colcode <- plotclr[colornum]
    # Note order of data in file is in terms of increasing FIPS codes, 
    # which is the same as in the map function (see county.fips)
    map("county", "ohio",col=colcode,fill=T)
    title(title_)
    leg.txt <- paste("[",brks[nclr],",",brks[nclr+1],"]",sep="")
    for(i in (nclr-1):1){
        leg.txt <- append(leg.txt,paste("[",brks[i],",",brks[i+1],")",sep=""))
    }
    leg.txt <- rev(leg.txt)
    legend("bottomright",legend=leg.txt,fill=plotclr,bty="n",cex=.8)
}

