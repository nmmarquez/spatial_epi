source("~/Documents/Classes/spatial_epi/homework/hw1/ohio.dat")
# As an example, map the total deaths, across all strata, by county
Y <- NULL
ncounty <- length(unique(ohio$fips))
for (i in 1:ncounty){
  indi <- (i-1)*16+seq(1,16)
  Y[i] <- sum(ohio$deaths[indi])
}
library(maps)
library(sp)
library(RColorBrewer)
plotvar <- Y # variable we want to map
nclr <- 8 # next few lines set up the color scheme for plotting
plotclr <- brewer.pal(nclr,"BuPu")
brks <- round(quantile(plotvar,probs=seq(0,1,1/(nclr))),digits=1)
colornum <- findInterval(plotvar,brks,all.inside=T)
colcode <- plotclr[colornum]
# Note order of data in file is in terms of increasing FIPS codes, which is the same
# as in the map function (see county.fips)
# map("county", "ohio",col=colcode,fill=T)
# title("Lung cancer deaths in Ohio in 1988")
# leg.txt <- paste("[",brks[nclr],",",brks[nclr+1],"]",sep="")
#   for(i in (nclr-1):1){
#     leg.txt <- append(leg.txt,paste("[",brks[i],",",brks[i+1],")",sep=""))
# }
# leg.txt <- rev(leg.txt)
# legend("bottomright",legend=leg.txt,fill=plotclr,bty="n",cex=.8)
 