rm(list = ls()); gc()
setwd('~/Work/CoolingTrend/')
library(mgcv)
library(parallel)
source('code/functions2.R')

stations <- read.table('data/IonoStations.txt', header = T, sep = '|', stringsAsFactors = F)
data <- read.csv('data/DATA_PROC.csv', header = T, stringsAsFactors = F)
#stations for research
stations<-stations[stations$LAT <= 60 & stations$LAT >= 30,] 
nrow(stations)
data = data[data$CODE %in% stations$CODE,]
stations = stations[stations$CODE %in% data$CODE,]
length(unique(data$CODE))
length(stations$CODE)
data$M = sin(data$I/180 * pi)*cos(data$I/180 * pi)


data <- data[!(data$YEAR %in% c(2014, 2015)),]

max.range = 2
fitted <- vector('list', 20)
data$INT <- 1

variables <- c("TREND", "F107", "H_NT", "INT")
fm <- paste0("FREQ~0+",paste0(variables, collapse = "+"))

fit.gam = glm(as.formula(fm), data = data)
data$RES = data$FREQ - predict(fit.gam, type = 'response')

for(iteration in 1:20){

  byyear = by(data, data$YEAR, function(x) x)
  covList <- mclapply(byyear, function(x) {getSTCovariance(x,stations)}, mc.cores = detectCores())
  
  ln <- length(byyear)
  ZMZ <- 0
  ZY <- 0
  for(i in 1:ln){
    Sig       <- as.matrix(covList[[i]])
    diag(Sig) <- mean(byyear[[i]]$RES^2)
    SigI      <- solve(Sig)
    Z <- byyear[[i]][,variables]
    Z <- as.matrix(Z)
    Y <- as.matrix(byyear[[i]][,c("FREQ")]) 
    ZMZ <- ZMZ + t(Z) %*% SigI %*% Z
    ZY  <- ZY  + t(Z) %*% SigI %*% Y
  }
  
  beta <- solve(ZMZ) %*% ZY
  cat(beta, "\n")
  fitted[[iteration]] <- beta
  XX <- data[,variables]
  data$RES <- data$FREQ - (as.matrix(XX) %*% beta)[,1]
  cat(iteration,'\n')
  gc()
  for(i in 1:length(covList))
    cat(i, " " ,try(is.na(covList[[i]][1,1])), "\n")
}

beta[1,1] / sqrt(diag(solve(ZMZ)))[1]
plot.ts(sapply(fitted[1:20], function(x) x[1,1]))

#Models
#fit.gam = gam(FREQ ~ F107 + TREND, data = data)
#fit.gam = gam(FREQ ~ F107 + TREND + F_NT, data = data.w)
#fit.gam = gam(FREQ ~ F107 + TREND + X_NT + Y_NT + Z_NT, weights = WEIGHTS, data = data.w)
#fit.gam = gam(FREQ ~ F107 + TREND + H_NT, weights = WEIGHTS, data = data.w)
#fit.gam = gam(FREQ ~ F107 + TREND + H_NT + F_NT, weights = WEIGHTS, data = data.w)
#fit.gam = gam(FREQ ~ F107 + TREND + M, weights = WEIGHTS, data = data.w)
#fit.gam = gam(FREQ ~ F107 + TREND + M +  F_NT, weights = WEIGHTS, data = data.w)
#fit.gam = gam(FREQ ~ F107 + TREND + M +  X_NT + Y_NT + Z_NT, weights = WEIGHTS, data = data.w)
#fit.gam = gam(FREQ ~ F107 + TREND + M + H_NT, weights = WEIGHTS, data = data.w)
#fit.gam = gam(FREQ ~ F107 + TREND + M + H_NT + F_NT, weights = WEIGHTS, data = data.w)