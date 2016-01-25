x <- read.csv("alaska.csv")
xcoord <- x[,3]
ycoord <- x[,4]
gender <- x[,8]
age <- x[,9]
race1 <- x[,11]
race2 <- 2*x[,12]
race3 <- 3*x[,13]
race4 <- 4*x[,14]
race5 <- 5*x[,15]
race6 <- 6*x[,16]
race7 <- 7*x[,17]
asthma <- x[,23]
smoke <- x[,28]
income <- x[,54]
school <- x[,55]
parental <- x[,22]
exposure <- x[,60]
illness <- x[,21]

x2 <- cbind(gender,age,illness,smoke,income,school,exposure,parental)
x2[x2== -999] <- NA
x2 <- na.omit(x2) # Form a new dataset including complete cases only

gender <- factor(x2[,1])
age <- x2[,2]
illness <- x2[,3]; illness2 <- rep(0,length(illness))

illness2[illness==1] <- 0
illness2[illness==2] <- 0
illness2[illness==3] <- 1
illness2[illness==4] <- 1

smoke <- factor(x2[,4])
income <- factor(x2[,5])
school <- factor(x2[,6])
exposure <- 10*x2[,7]/max(x2[,7])
parental <- factor(x2[,8])
