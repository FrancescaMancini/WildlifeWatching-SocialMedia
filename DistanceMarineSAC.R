#######################################
##Calculate distance to marine SAC
#######################################
library(rgeos)
library(rgdal)
library(classInt)

data<-read.table(".//data//CombinedData_v9.txt",stringsAsFactors = F,header=T)

SAC<-readOGR("H://National Scale Data//ChapterIII//Chapter_3//data//SAC//SAC_Marine.shp",layer="SAC_Marine")

coordinates(data)<-c("Longitude","Latitude")

wgs84 <- '+proj=longlat +datum=WGS84'
data@proj4string <- CRS(wgs84)
data_bng <- spTransform(data, SAC@proj4string)

Fdist <- list()
for(i in 1:dim(data_bng)[1]) {
  pDist <- vector()
  for(j in 1:dim(SAC)[1]) { 
    pDist <- append(pDist, gDistance(data_bng[i,],SAC[j,])) 
  }
  Fdist[[i]] <- pDist
} 

PolyDist <- unlist(lapply(Fdist, FUN=function(x) min(x)[1]))

summary(PolyDist)

data_bng@data <- data.frame(data_bng@data, PDist=PolyDist)


cuts <- classIntervals(data_bng@data$PDist, 10, style="quantile") 
plotclr <- colorRampPalette(c("cyan", "yellow", "red"))( 20 )
colcode <- findColours(cuts, plotclr)
plot(SAC,col="black")
plot(data_bng, col=colcode, pch=19, add=TRUE)

write.table(data_bng@data,".//data//CombinedData_v10.txt")
