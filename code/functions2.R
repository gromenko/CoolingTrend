library(mgcv)
library(parallel)

estimate_tempr <- function(xss, covar){
  sx <- by(xss, xss$MONTH, function(x) x)
  dates = names(sx)
  if(length(dates)!=12) stop('Matrix is not 12 X 12')
  surf = matrix(nc = 12, nr = 12, 0)
  weights = matrix(nc = 12, nr = 12, 0)
  for(i in 1:12){
    for(j in i:12){
      if(i!=j){
        res1 = data.frame(RES1 = sx[[dates[i]]]$RES, CODE = sx[[dates[i]]]$CODE, stringsAsFactors = F)
        res2 = data.frame(RES2 = sx[[dates[j]]]$RES, CODE = sx[[dates[j]]]$CODE, stringsAsFactors = F)
        psi = merge(res1, res2, by.x = 'CODE', by.y = 'CODE', all = F)
        if(nrow(psi) > 0){
          nu = solve(covar[psi$CODE,psi$CODE]^2, rep(1, nrow(psi)))
          nu[nu < 0] = 0; nu = nu/sum(nu)
          #surf[i,j] =  sum(psi$RES1 * psi$RES2 * nu)
          surf[i,j] =  mean(psi$RES1 * psi$RES2)
          surf[j,i] = surf[i,j];
          weights[i,j] <- nrow(psi)
          weights[j,i] <- nrow(psi)
        }
      } 
    }
  }
  id <- expand.grid(1:12, 1:12)
  tmpDF <- data.frame(Z = as.vector(surf), Y = id[,1], X = id[,2], W = as.vector(weights))
  fit <- gam(Z ~ te(X,Y, k =c(12,12)), weights = W, data = tmpDF)
  tmpDF$SM <- predict(fit, newdata = tmpDF, type = 'response')
  surfS <- tmpDF$SM
  dim(surfS) <- c(12,12)
  
  return(surfS / sum(svd(surfS)$d))
}

#names(byyear)
#tmp = byyear[[33]]
getSTCovariance <- function(tmp, stations){
  tryCatch(
    {
      prelim <- variog_est(tmp, stations, max.range = 2)
      covarS <- fit_covariance(prelim$VARIOG, prelim$DIST.M)
      covarT <- estimate_tempr(tmp, covarS$COVM)
      len <- nrow(tmp)
      covarST <- matrix(nc = len, nr = len, NA)
      indc <- tmp$CODE
      indt <- as.numeric(tmp$MONTH)
      sg2 <- covarS$CO
      for(i in 1:len){
        for(j in i:len){
          covarST[i,j] <- covarS$COVM[indc[i], indc[j]] * covarT[indt[i], indt[j]]
          covarST[j,i] <- covarST[i,j] 
        }
      }
      diag(covarST) <- sg2
      return(covarST)
    }, error = function(e){return(NA)}
  )}

variog_est <- function(tmp, stations, max.range){

  covariance = c()
  distance = c()
  sx    = by(tmp, tmp$CODE, function(x) x)
  ds.l  = by(stations, stations$CODE, function(x) x[,c(2,3)])
  dist.matr = matrix(nc = length(names(sx)), nr = length(names(sx)), 0)
  colnames(dist.matr) = names(sx)
  rownames(dist.matr) = names(sx)
  for(c1 in names(sx)){
    for(c2 in names(sx)){
      dd = ch_dist(ds.l[[c1]]$LON, ds.l[[c2]]$LON,
                   ds.l[[c1]]$LAT, ds.l[[c2]]$LAT)
      dist.matr[c1,c2] = dd
      if(c1 != c2 & dim(sx[[c1]])[1] >= 0 & dim(sx[[c2]])[1] >= 0){
        sxx = merge(sx[[c1]][,c('DATE','RES')], sx[[c2]][,c('DATE','RES')], 
                    by.x = 'DATE', by.y = 'DATE', all.x = T)
        cv = na.omit(sxx[,2] - sxx[,3])
        if(length(cv) > 1){
          covariance = c(covariance, cv)
          distance   = c(distance,   rep(dd, length(cv)))
        }
      }
    }
  }
  by = 0.1
  covariog <- c()
  bins <- seq(0,max.range - by, by = by)
  for(i in 1:length(bins)){
    i1 <- distance >= bins[i]
    i2 <- distance < bins[i] + by
    xs <- covariance[i1&i2]
    if(length(xs) > 0){
      #xss <- c() 
      #l <- length(xs)
      #for(j in 1:(l-1)){
      #  xss <- c(xss, abs(xs[(j+1):l] - xs[j]))
      #}
      #quantile(xss, prob = 0.25)
      #covariog <- rbind(covariog, c( bins[i] + by/2, mean(xs^2), median(xs^2), (2.2191 * quantile(xss, prob = 0.25))^2, length(xss)))
      covariog <- rbind(covariog, c( bins[i] + by/2, mean(xs^2)))
    }
  }
  return(list(VARIOG = data.frame(DIST = covariog[,1], COV = covariog[,2]), 
              DIST.M = dist.matr))
}



fit_covariance <- function(variog, dist.matr){
  #Gaussian model
  fit.nls.g <- NULL
  modelG    <- NULL
  RSS.g     <- NULL
  
  try(fit.nls.g <- nls(COV ~ (c1 + c3) - c1 * exp(- DIST^2 * c2^2), 
                       start = list(c1 = 2, c2 = 5, c3 = 0.1), 
                       lower = c(0,0.5,0), upper = c(Inf, Inf, Inf), 
                       algorithm = 'port', data = variog))
  
  if(is.null(fit.nls.g)){
    try(fit.nls.g <- nls(COV ~ c1 - c1 * exp(- DIST^2 * c2^2), 
                         start = list(c1 = 2, c2 = 5), 
                         lower = c(0,0.5), upper = c(Inf, Inf), 
                         algorithm = 'port', data = variog))
    if(is.null(fit.nls.g)){
      RSS.g <- NA
    } else { modelG <- 'RED'; RSS.g = var(summary(fit.nls.g)$residuals) }
  } else { modelG <- 'FULL'; RSS.g = var(summary(fit.nls.g)$residuals) }
  
  #Exponential model
  fit.nls.e <- NULL 
  modelE    <- NULL
  RSS.e     <- NULL
  
  try(fit.nls.e <- nls(COV ~ (c1 + c3) - c1 * exp(- DIST * c2), 
                       start = list(c1 = 2, c2 = 5, c3 = 0.1), 
                       lower = c(0,0.5,0), upper = c(Inf, Inf, Inf), 
                       algorithm = 'port', data = variog))
  
  if(is.null(fit.nls.e)){
    try(fit.nls.e <- nls(COV ~ c1 - c1 * exp(- DIST * c2), 
                         start = list(c1 = 2, c2 = 5), 
                         lower = c(0,0.5), upper = c(Inf, Inf), 
                         algorithm = 'port', data = variog))
    if(is.null(fit.nls.e)){
      RSS.e <- NA
    } else {modelE <- 'RED'; RSS.e = var(summary(fit.nls.e)$residuals) }
  } else {modelE <- 'FULL'; RSS.e = var(summary(fit.nls.e)$residuals) }
  
  RSS.n = NA
  fit = c(RSS.n, RSS.e, RSS.g)
  names(fit) = c('N', 'E', 'G')
  fit = na.omit(fit)
  if("G" %in% names(fit)) { 
    model = "G"
  } else {model = names(which(min(fit) == fit))}
  
  if(model == 'E'){
    cff  <- coef(fit.nls.e)
    covm <- cff[1] / 2 * exp(-dist.matr * cff[2])
    co   <- cff[1] / 2
    co   <- ifelse(modelE == 'FULL', co + cff[3] / 2, co)
    mdl      <- model
    mdlType  <- modelE
  }
  if(model == 'G'){
    cff  <- coef(fit.nls.g)
    covm <- cff[1] / 2 * exp(-dist.matr^2 * cff[2]^2)
    co   <- cff[1] / 2
    co   <- ifelse(modelG == 'FULL', co + cff[3] / 2, co)
    mdl      <- model
    mdlType  <- modelG
  }
  if(model == 'N'){
    covm = diag(dim(dist.matr)[1])
    #diag(covm) = mean(cv.df$COV)
    diag(covm) = 0
    colnames(covm) = colnames(dist.matr)
    rownames(covm) = colnames(dist.matr)
  }
  return(list(COVM = covm, CO = co, MODEL = model, TYPE = mdlType))
}


ch_dist = function(ln1,ln2, lt1, lt2){
  lon1 = ln1 * pi / 180; lon2 = ln2 * pi / 180
  lat1 = lt1 * pi / 180; lat2 = lt2 * pi / 180
  2*sqrt(sin((lon1 - lon2)/2)^2 + cos(lon1) * cos(lon2) * sin((lat1 - lat2)/2)^2)
}

make_surface<- function(fit.bm.te, main){
  xs = seq(1,12, by = 0.5)
  ys = seq(50,250, by = 10)
  coord = do.call(rbind, sapply(ys, function(x) lapply(xs, function(y) c(y, x)) ))
  new.data = data.frame(MONTH = coord[,1], F107 = coord[,2], TREND = 0)
  IM = predict(fit.bm.te, newdata = new.data)
  dim(IM) = c(23, length(IM)/23)
  persp(x = xs, y = ys, z = IM, theta = 45, phi = 25, ticktype="detailed",
        shade=.2, col = 'gray90', xlab = 'Month', ylab = 'SRF', zlab = 'f', main = main)
}


modelRun <- function(variables, filename, niter, ys){
  
  write("Starting iteration",file=filename,append=TRUE)
  
  stations <- read.table('data_noon/IonoStations.txt', header = T, sep = '|', stringsAsFactors = F)
  data <- read.csv('data_noon/DATA_PROC.csv', header = T, stringsAsFactors = F)
  stations<-stations[stations$LAT <= 60 & stations$LAT >= 30,] 
  data = data[data$CODE %in% stations$CODE,]
  stations = stations[stations$CODE %in% data$CODE,]
  data$M = sin(data$I/180 * pi)*cos(data$I/180 * pi)
  data <- data[!(data$YEAR %in% c(2014, 2015)),]
  data <- data[data$YEAR %in% ys,]
  data$INT <- 1
  
  fm <- paste0("FREQ~0+",paste0(variables, collapse = "+"))
  fit.gam = glm(as.formula(fm), data = data)
  data$RES = data$FREQ - predict(fit.gam, type = 'response')
  
  write(paste0("E0:",paste0(coef(summary(fit.gam))[,1], collapse = ",")),file=filename,append=TRUE)
  write(paste0("V0:",paste0(coef(summary(fit.gam))[,2], collapse = ",")),file=filename,append=TRUE)
  write(paste0("T0:",paste0(coef(summary(fit.gam))[,3], collapse = ",")),file=filename,append=TRUE)
  write(paste0("P0:",paste0(coef(summary(fit.gam))[,4], collapse = ",")),file=filename,append=TRUE)
  #cat(coef(summary(fit.gam))[,1], "\n")
  
  for(iteration in 1:niter){
    write(paste0("Iteration:",iteration),file=filename,append=TRUE)
    
    byyear = by(data, data$YEAR, function(x) x)
    covList <- mclapply(byyear, function(x) getSTCovariance(x, stations), mc.cores = detectCores())
    
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
    stde <- sqrt(diag(solve(ZMZ)))
    
    write(paste0("E:",paste0(beta, collapse = ",")),file=filename,append=TRUE)
    write(paste0("V:",paste0(stde, collapse = ",")),file=filename,append=TRUE)
    write(paste0("T:",paste0(beta/stde, collapse = ",")),file=filename,append=TRUE)
    #write(paste0("P0:",paste0(coef(summary(fit.gam))[,4], collapse = ",")),file=filename,append=TRUE)
    XX <- data[,variables]
    data$RES <- data$FREQ - (as.matrix(XX) %*% beta)[,1]
    gc()
  }
}
