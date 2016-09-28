rm(list = ls())
setwd('~/Work/CoolingTrend/')
library(maps)
library(mapproj)

make_vert <- function(lat, y){
  x <- c(rep(lat, length(y))); 
  lines(mapproject(list(x=c(x), y=c(y))), col=col, lwd=lwd)
}

make_hor <- function(lon, x){
  y <- c(rep(lon, length(x)));
  lines(mapproject(list(x=c(x), y=c(y))), col=col, lwd=lwd)
}


stations <- read.table('data/IonoStations.txt', header = T, sep = '|', stringsAsFactors = F)
data <- read.csv('data/DATA_PROC.csv', header = T, stringsAsFactors = F)
#stations for research
stations<-stations[stations$LAT <= 60 & stations$LAT >= 30,] 
nrow(stations)
data = data[data$CODE %in% stations$CODE,]
stations = stations[stations$CODE %in% data$CODE,]
loc <- stations
data(worldMapEnv)
nrow(stations)

dts <- c()
for(i in 1958:2015){
  for(j in 1:12){
    dts <- c(dts, paste(i,ifelse(j < 10, paste0(0,j),j),'01', sep = '-'))
  }}

count <- c()
for(ds in dts){
  count <- c(count, length(which(data$DATE == ds)))
}

cnt.df <- data.frame(DATE = as.Date(dts), COUNT = count)

postscript("figs/Miss_U.eps", height=4, width=8)
  par(mar = c(2,4,1,4))
  plot(cnt.df$DATE, cnt.df$COUNT, type = 'l', xlab = NA, ylab = "Number of stations")
dev.off()

data$DATE <- as.Date(data$DATE)
range(data$FREQ)


postscript("figs/Data_U.eps", height=4, width=8)
  par(mar = c(2,4,1,4))
  plot(cnt.df$DATE, cnt.df$COUNT , ylim = c(3,17), type = 'n', ylab = "foF2, MHz")
  fit <- lm(FREQ ~ F107, data = data)
  data$PR <- predict(fit)
  srf <- c()
  for(dt in cnt.df$DATE){
    tmp <- unique(data[data$DATE == dt,'PR'])
    tmp <- ifelse(length(tmp) > 0, tmp, NA)
    srf <- c(srf, tmp)
  }
  
  for(cd in unique( data$CODE )){
    ind <- data$CODE == cd
    tmp <- data[data$CODE == cd, c('DATE', 'FREQ')]
    tmp <- merge(cnt.df, tmp, by = 'DATE', all.x = T)
    lines(tmp$DATE, tmp$FREQ, col = 'gray60')
  }
  
  lines(cnt.df$DATE, srf, lwd = 1.5)
  
  a <- coef(summary(fit))[1,1]
  b <- coef(summary(fit))[2,1]
  
  x <- seq(0, 320, by = 50)
  at = x * b + a
  axis(4, at = at, labels = x)
  mtext(expression("SRF, W/"~"m"^2~"/Hz"), side = 4, line = 2.5)
dev.off()


est <- read.table('data/sample_est.txt', header = T, sep = ",")
esd <- read.table('data/sample_sd.txt', header = T, sep = ",")
data$RES <- data$FREQ - (data$F107 * est$F107 + est$INT + data$H_NT * est$H_NT)

postscript("figs/Data_Trend.eps", height=4, width=8)
  par(mar = c(2,4,1,4))
  plot(cnt.df$DATE, cnt.df$COUNT , ylim = c(-5,5), type = 'n', ylab = "foF2, MHz")
  
  for(cd in unique( data$CODE )){
    ind <- data$CODE == cd
    tmp <- data[data$CODE == cd, c('DATE', 'RES')]
    tmp <- merge(cnt.df, tmp, by = 'DATE', all.x = T)
    lines(tmp$DATE, tmp$RES, col = 'gray60')
  }
  
  abline(h = 0, lty = 'dashed')
  tm <- as.numeric(cnt.df$DATE - min(cnt.df$DATE))
  lines(cnt.df$DATE, tm * est$TREND, col = 1, lwd = 1.5)
  lines(cnt.df$DATE, tm * est$TREND + 3 * tm * esd$TREND , col = 1, lwd = 1, lty = 'dotted')
  lines(cnt.df$DATE, tm * est$TREND - 3 * tm * esd$TREND , col = 1, lwd = 1, lty = 'dotted')
dev.off()


col1="gray60"
col="gray20"
lwd=0.2
postscript(file="figs/Stations.eps", height=8, width=8, 
           horizontal = F, onefile = FALSE, paper = "special")
  
  map(type="n", projection="albers", par=c(89,90), xlim=c(-180,180), ylim=c(20, 90), mar = c(4.1, 4.1, 4.1, 4.1))
  
  symbols(0, 0, 0.812, bg = col1, inches = FALSE, add=T, lwd=0.00000001)
  
  map("world", fill = T, interior = FALSE, boundary = TRUE, projection="albers", 
      par=c(90,90), xlim=c(-180,180), ylim=c(21.5, 90), col="white", add=T, lwd=0.9, border = "white")
  
  ######################################
  x <- c(seq(-180, 180, by = 0.1)); 
  for(lon in seq(20,80,by=20)){
    make_hor(lon, x)
  }
  ######################################
  y <- c(seq(20, 80, by=0.1));
  for(lat in seq(-160, 180, by = 20)){
    make_vert(lat,y)
  }
  #lines(mapproject(list(x=c(line75[,1]), y=c(line75[,2]))),  lwd=1, lty=5, col="blue")
  #lines(mapproject(list(x=c(line80[,1]), y=c(line80[,2]))),  lwd=1, lty=5, col="blue")
  
  x <- c(seq(0, 180, by = 20)); 
  y <- c(rep(13, length(x)));
  text(mapproject(list(x=c(x), y=c(y))), 
  labels= c(expression(paste("0"^o)  ),
            expression(paste("20"^o) ),
            expression(paste("40"^o) ),
            expression(paste("60"^o) ),
            expression(paste("80"^o) ),
            expression(paste("100"^o)),
            expression(paste("120"^o)),
            expression(paste("140"^o)),
            expression(paste("160"^o)),
            expression(paste("180"^o))), col=1, cex=0.7)
  
  x <- c(seq(-20, -160, by = -20)); 
  y <- c(rep(13, length(x)));
  text(mapproject(list(x=c(x), y=c(y))), 
  labels= c(expression(paste("-20"^o) ),
            expression(paste("-40"^o) ),
            expression(paste("-60"^o) ),
            expression(paste("-80"^o) ),
            expression(paste("-100"^o)),
            expression(paste("-120"^o)),
            expression(paste("-140"^o)),
            expression(paste("-160"^o))), col=1, cex=0.7)
  ######################################
  y <- c(seq(20,80, by=20))-3;
  x <- c(rep(-170, length(y))); 
  text(mapproject(list(x=c(x), y=c(y))), 
  labels=c(expression(paste("20"^o)),
  		expression(paste("40"^o)),
  		expression(paste("60"^o)),
  		expression(paste("80"^o))), col=1, cex=0.7)
  
  points(mapproject(list(x=loc$LON, y=loc$LAT)), col="black", pch=19, cex=0.5)
dev.off()
