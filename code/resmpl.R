rm(list = ls()); gc()
setwd('~/Work/CoolingTrend')
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
data$INT <- 1

variables <- c("TREND", "F107", "H_NT", "INT")
fm <- paste0("FREQ~0+",paste0(variables, collapse = "+"))

xr <- unique(data$YEAR)
hx <- head(xr,1)
tx <- tail(xr,1)

result <- c()
centers <- xr[3:(length(xr)-2)]
for(cn in centers){
  ml <- min(abs(cn - hx), abs(cn - tx))
  for(i in 2:ml){
    ys <- (cn - i) : (cn + i)
    tmp <- data[data$YEAR %in% ys ,]
    fit <- glm(as.formula(fm), data = tmp)
    result <- rbind(result, c(cn, i, coef(summary(fit))[1,c(1,4)]))
  }
  cat(cn, '\n')
}

dir <- 'resmplExp/'
fnames <- list.files(dir)
xs <- gsub(pattern = 'log.txt', '', fnames)
year   <- as.numeric(substr(xs, 2, 5)) 
length <- as.numeric(substr(xs, 7, 10)) 

resultST <- c()
for(i in 1:length(fnames)){
  tmp <- readLines(paste0(dir, fnames[i]))
  if(length(grep('Iteration:', tmp)) > 15){
    trEst <- na.omit(as.numeric(sapply(strsplit(gsub('E:','',tmp[grep('Iteration:', tmp) + 1]), '[,]'), function(x) x[1])))
    trSde <- na.omit(as.numeric(sapply(strsplit(gsub('V:','',tmp[grep('Iteration:', tmp) + 2]), '[,]'), function(x) x[1])))
    mn <- mean(trEst[10:length(trEst)])
    md <- median(trEst[10:length(trEst)])
    
    mnSe <- mean(trSde[10:length(trSde)])
    mdSe <- median(trSde[10:length(trSde)])
    ind <- 0
    ind <- ifelse(length(which(trEst > 0)), 1, 0)
    
    if(abs(mn - md) < tail(trSde,1)){
      e <- mn * 365 * 1000
      s <- mnSe * 365 * 1000
      p <- 2 * pnorm(abs(e / s),lower.tail = F)
      resultST <- rbind(resultST, c(year[i], length[i], e, p))
    }
  }
}

xlim <- range(result[,1])
ylim <- range(result[,2])
atv = seq(1960, 2010, by = 10)
ath = seq(5, 25, by = 5)
resps <- result[result[,3] >  0 & result[,4] <= 0.05,]
respu <- result[result[,3] >  0 & result[,4] >  0.05,]
resns <- result[result[,3] <= 0 & result[,4] <= 0.05,]
resnu <- result[result[,3] <= 0 & result[,4] >  0.05,]

respsST <- resultST[resultST[,3] >  0 & resultST[,4] <= 0.05,]
respuST <- resultST[resultST[,3] >  0 & resultST[,4] >  0.05,]
resnsST <- resultST[resultST[,3] <= 0 & resultST[,4] <= 0.05,]
resnuST <- resultST[resultST[,3] <= 0 & resultST[,4] >  0.05,]
cex = 0.8

postscript('figs/F6.eps', width = 9, height = 3)

  par(mar = c(2,2,1,1))
  layout(matrix(nc = 2, nr = 1, 1:2))
  plot(1,1, xlim = xlim, ylim = ylim, type = "n", xlab = NA, ylab = NA)
  abline(v = atv, lty = 'dotted')
  abline(h = ath, lty = 'dotted')
  
  points(resps[,1], resps[,2], pch = 3, col = "red", cex = cex)
  points(respu[,1], respu[,2], pch = 5, col = "red", cex = cex)
  points(resns[,1], resns[,2], pch = 19, col = "blue", cex = cex)
  points(resnu[,1], resnu[,2], pch = 1, col = "blue", cex = cex)
  
  plot(1,1, xlim = xlim, ylim = ylim, type = "n", xlab = NA, ylab = NA)
  abline(v = atv, lty = 'dotted')
  abline(h = ath, lty = 'dotted')
  
  points(respsST[,1], respsST[,2], pch = 3, col = "red", cex = cex)
  points(respuST[,1], respuST[,2], pch = 5, col = "red", cex = cex)
  points(resnsST[,1], resnsST[,2], pch = 19, col = "blue", cex = cex)
  points(resnuST[,1], resnuST[,2], pch = 1, col = "blue", cex = cex)

dev.off()
