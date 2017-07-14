#########################################
##########Google maps API
#########################################
setwd("H:/National Scale Data/ChapterIII/ChapterIII")

library(httr)
library(jsonlite)
require(maptools)
library(rgdal)

location<-read.table("H:\\National Scale Data\\ChapterIII\\Locations.txt",header=T)

radius<-2500

coordinates(location)<-c("X","Y")
proj4string(location)=CRS("+init=epsg:27700")
location<-spTransform(location,CRS("+init=epsg:4326"))
coords<-as.data.frame(unlist(location@coords))

#day1<-coords[1:1000,]
#day2<-coords[1001:2000,]
#day3<-coords[2001:3000,]
#day4<-coords[3001:4000,]
#day5<-coords[4001:5000,]
#day6<-coords[5001:6000,]
#day7<-coords[6001:7000,]
#day8<-coords[7001:8000,]
#day9<-coords[9001:10000,]
#day10<-coords[10001:11000,]
#day11<-coords[11001:12000,]
#day12<-coords[12001:12845,]


###############################
######Hotels
###############################

hotels<-NULL

for(i in 1:dim(coords)[1]){

res<-GET(paste( "https://maps.googleapis.com/maps/api/place/radarsearch/json?location=",
                paste(coords$Y[i],coords$X[i],sep=","),"&radius=",radius,
                "&type=lodging&key=AIzaSyDRWYQO3gTUNXOQa0t1FerdyOIyfUU8vKQ",sep=""))

myHotels<- fromJSON(content(res,"text"),flatten=T)$results
if (is.null(dim(myHotels))==FALSE){
hotels<-rbind(hotels,myHotels[,c("place_id","geometry.location.lat","geometry.location.lng")])
}}

write.table(hotels,"H:\\National Scale Data\\ChapterIII\\Hotels.txt",col.names=T,row.names=F,sep="\t")

details<-array(NA,c(length(hotels$place_id),3))
reviews<-data.frame(NA,NA,NA,NA)
colnames(reviews)<-c("rating","text","time","hotel")


for (i in 1:length(hotels$place_id)){
res<-GET( paste("https://maps.googleapis.com/maps/api/place/details/json?placeid=",hotels$place_id[i], "&key=AIzaSyDRWYQO3gTUNXOQa0t1FerdyOIyfUU8vKQ",sep=""))
myCost<- fromJSON(content(res,"text"),flatten=T)$result
if (is.null(myCost$name)==FALSE){
details[i,1]<-myCost$name}
if(is.null(myCost$price_level)==FALSE){
details[i,2]<-myCost$price_level}
if(is.null(myCost$rating)==FALSE){
details[i,3]<-myCost$rating}
if(is.null(myCost$reviews)==FALSE){
tmp<-cbind(myCost$reviews[which(names(myCost$reviews)=="rating" |names(myCost$reviews)=="text" | names(myCost$reviews)=="time" )],rep(myCost$name,dim(myCost$review)[1]))
colnames(tmp)<-c("rating","text","time","hotel")
reviews<-rbind(reviews,tmp)
}
}

details_df<-data.frame(Name=details[,1],Price=details[,2],Rating=details[,3])

hotels_det<-cbind(hotels,details_df)


write.table(hotels_det,"H:\\National Scale Data\\ChapterIII\\HotelsDetail.txt",col.names=T,row.names=F,sep="\t")

write.table(reviews,"H:\\National Scale Data\\ChapterIII\\Reviews.txt",col.names=T,row.names=F,sep="\t")


##############################
########Tour Operators
##############################

TO<-NULL

for(i in 1:dim(coords)[1]){

res<-GET(paste( "https://maps.googleapis.com/maps/api/place/radarsearch/json?location=",
                paste(coords$Y[i],coords$X[i],sep=","),"&radius=",radius,
                "&keyword=wildlife&type=tour_operator&key=AIzaSyDRWYQO3gTUNXOQa0t1FerdyOIyfUU8vKQ",sep=""))

myTO<- fromJSON(content(res,"text"),flatten=T)$results
if (is.null(dim(myTO))==FALSE){
TO<-rbind(TO,myTO[,c("place_id","geometry.location.lat","geometry.location.lng")])
}}

write.table(TO,"H:\\National Scale Data\\ChapterIII\\TourOperators.txt",col.names=T,row.names=F,sep="\t")

details<-array(NA,c(length(TO$place_id),4))
reviews<-data.frame(NA,NA,NA,NA)
colnames(reviews)<-c("rating","text","time","name")


for (i in 1:length(TO$place_id)){
res<-GET( paste("https://maps.googleapis.com/maps/api/place/details/json?placeid=",TO$place_id[i], "&key=AIzaSyDRWYQO3gTUNXOQa0t1FerdyOIyfUU8vKQ",sep=""))
myDetails<- fromJSON(content(res,"text"),flatten=T)$result
if (is.null(myDetails$name)==FALSE){
details[i,1]<-myDetails$name}
if(is.null(myDetails$price_level)==FALSE){
details[i,2]<-myDetails$price_level}
if(is.null(myDetails$rating)==FALSE){
details[i,3]<-myDetails$rating}
if(is.null(myDetails$types)==FALSE){
details[i,4]<-paste(myDetails$types,collapse="/")}
if(is.null(myDetails$reviews)==FALSE){
tmp<-cbind(myDetails$reviews[which(names(myDetails$reviews)=="rating" |names(myDetails$reviews)=="text" | names(myDetails$reviews)=="time" )],rep(myDetails$name,dim(myDetails$review)[1]))
colnames(tmp)<-c("rating","text","time","name")
reviews<-rbind(reviews,tmp)
}
}

details_df<-data.frame(Name=details[,1],Price=details[,2],Rating=details[,3],Type=details[,4])

TO_det<-cbind(TO,details_df)


write.table(TO_det,"H:\\National Scale Data\\ChapterIII\\TODetail.txt",col.names=T,row.names=F,sep="\t")

write.table(reviews,"H:\\National Scale Data\\ChapterIII\\TOReviews.txt",col.names=T,row.names=F,sep="\t")


###############################
######Airports
###############################

airports<-NULL

for(i in 1:dim(coords)[1]){

res<-GET(paste( "https://maps.googleapis.com/maps/api/place/radarsearch/json?location=",
                paste(coords$Y[i],coords$X[i],sep=","),"&radius=",radius,
                "&type=airport&key=AIzaSyDRWYQO3gTUNXOQa0t1FerdyOIyfUU8vKQ",sep=""))

myAirports<- fromJSON(content(res,"text"),flatten=T)$results
if (is.null(dim(myAirports))==FALSE){
airports<-rbind(airports,myAirports[,c("place_id","geometry.location.lat","geometry.location.lng")])
}}

write.table(airports,"H:\\National Scale Data\\ChapterIII\\Airports.txt",col.names=T,row.names=F,sep="\t")

details<-array(NA,c(length(airports$place_id),3))
reviews<-data.frame(NA,NA,NA,NA)
colnames(reviews)<-c("rating","text","time","name")


for (i in 1:length(airports$place_id)){
res<-GET( paste("https://maps.googleapis.com/maps/api/place/details/json?placeid=",airports$place_id[i], "&key=AIzaSyDRWYQO3gTUNXOQa0t1FerdyOIyfUU8vKQ",sep=""))
myDetails<- fromJSON(content(res,"text"),flatten=T)$result
if (is.null(myDetails$name)==FALSE){
details[i,1]<-myDetails$name}
if(is.null(myDetails$price_level)==FALSE){
details[i,2]<-myDetails$price_level}
if(is.null(myDetails$rating)==FALSE){
details[i,3]<-myDetails$rating}
if(is.null(myDetails$reviews)==FALSE){
tmp<-cbind(myDetails$reviews[which(names(myDetails$reviews)=="rating" |names(myDetails$reviews)=="text" | names(myDetails$reviews)=="time" )],rep(myDetails$name,dim(myDetails$review)[1]))
colnames(tmp)<-c("rating","text","time","name")
reviews<-rbind(reviews,tmp)
}
}

details_df<-data.frame(Name=details[,1],Price=details[,2],Rating=details[,3])

airports_det<-cbind(airports,details_df)


write.table(airports_det,"H:\\National Scale Data\\ChapterIII\\AirportsDetail.txt",col.names=T,row.names=F,sep="\t")

write.table(reviews,"H:\\National Scale Data\\ChapterIII\\AirportsReviews.txt",col.names=T,row.names=F,sep="\t")


###############################
######Bus station
###############################

bus<-NULL

for(i in 1:dim(coords)[1]){

res<-GET(paste( "https://maps.googleapis.com/maps/api/place/radarsearch/json?location=",
                paste(coords$Y[i],coords$X[i],sep=","),"&radius=",radius,
                "&type=bus_station&key=AIzaSyDRWYQO3gTUNXOQa0t1FerdyOIyfUU8vKQ",sep=""))

myBus<- fromJSON(content(res,"text"),flatten=T)$results
if (is.null(dim(myBus))==FALSE){
bus<-rbind(bus,myBus[,c("place_id","geometry.location.lat","geometry.location.lng")])
}}



details<-array(NA,c(length(bus$place_id),3))


for (i in 1:length(bus$place_id)){
res<-GET( paste("https://maps.googleapis.com/maps/api/place/details/json?placeid=",bus$place_id[i], "&key=AIzaSyDRWYQO3gTUNXOQa0t1FerdyOIyfUU8vKQ",sep=""))
myDetails<- fromJSON(content(res,"text"),flatten=T)$result
if (is.null(myDetails$name)==FALSE){
details[i,1]<-myDetails$name}
if(is.null(myDetails$price_level)==FALSE){
details[i,2]<-myDetails$price_level}
if(is.null(myDetails$rating)==FALSE){
details[i,3]<-myDetails$rating}
}

details_df<-data.frame(Name=details[,1],Price=details[,2],Rating=details[,3])

bus_det<-cbind(bus,details_df)


write.table(bus_det,"H:\\National Scale Data\\ChapterIII\\BusStationsDetail.txt",col.names=T,row.names=F,sep="\t")


###############################
######Subway station
###############################

subway<-NULL

for(i in 1:dim(coords)[1]){

res<-GET(paste( "https://maps.googleapis.com/maps/api/place/radarsearch/json?location=",
                paste(coords$Y[i],coords$X[i],sep=","),"&radius=",radius,
                "&type=subway_station&key=AIzaSyDRWYQO3gTUNXOQa0t1FerdyOIyfUU8vKQ",sep=""))

mySubway<- fromJSON(content(res,"text"),flatten=T)$results
if (is.null(dim(mySubway))==FALSE){
subway<-rbind(subway,mySubway[,c("place_id","geometry.location.lat","geometry.location.lng")])
}}


details<-array(NA,c(length(subway$place_id),3))


for (i in 1:length(subway$place_id)){
res<-GET( paste("https://maps.googleapis.com/maps/api/place/details/json?placeid=",subway$place_id[i], "&key=AIzaSyDRWYQO3gTUNXOQa0t1FerdyOIyfUU8vKQ",sep=""))
myDetails<- fromJSON(content(res,"text"),flatten=T)$result
if (is.null(myDetails$name)==FALSE){
details[i,1]<-myDetails$name}
if(is.null(myDetails$price_level)==FALSE){
details[i,2]<-myDetails$price_level}
if(is.null(myDetails$rating)==FALSE){
details[i,3]<-myDetails$rating}
}

details_df<-data.frame(Name=details[,1],Price=details[,2],Rating=details[,3])

subway_det<-cbind(subway,details_df)


write.table(subway_det,"H:\\National Scale Data\\ChapterIII\\SubwayStationsDetail.txt",col.names=T,row.names=F,sep="\t")


###############################
######Train station
###############################

train<-NULL

for(i in 1:dim(coords)[1]){

res<-GET(paste( "https://maps.googleapis.com/maps/api/place/radarsearch/json?location=",
                paste(coords$Y[i],coords$X[i],sep=","),"&radius=",radius,
                "&type=train_station&key=AIzaSyDRWYQO3gTUNXOQa0t1FerdyOIyfUU8vKQ",sep=""))

myTrain<- fromJSON(content(res,"text"),flatten=T)$results
if (is.null(dim(myTrain))==FALSE){
train<-rbind(train,myTrain[,c("place_id","geometry.location.lat","geometry.location.lng")])
}}


d<-duplicated(train$place_id)
train<-train[-which(d=="TRUE"),]

details<-array(NA,c(length(train$place_id),3))


for (i in 1:length(train$place_id)){
res<-GET( paste("https://maps.googleapis.com/maps/api/place/details/json?placeid=",train$place_id[i], "&key=AIzaSyDRWYQO3gTUNXOQa0t1FerdyOIyfUU8vKQ",sep=""))
myDetails<- fromJSON(content(res,"text"),flatten=T)$result
if (is.null(myDetails$name)==FALSE){
details[i,1]<-myDetails$name}
if(is.null(myDetails$price_level)==FALSE){
details[i,2]<-myDetails$price_level}
if(is.null(myDetails$rating)==FALSE){
details[i,3]<-myDetails$rating}
}

details_df<-data.frame(Name=details[,1],Price=details[,2],Rating=details[,3])

train_det<-cbind(train,details_df)


write.table(train_det,"H:\\National Scale Data\\ChapterIII\\TrainStationsDetail.txt",col.names=T,row.names=F,sep="\t")



###############################
######Car parks
###############################

parking<-NULL

for(i in 1:dim(coords)[1]){

res<-GET(paste( "https://maps.googleapis.com/maps/api/place/radarsearch/json?location=",
                paste(coords$Y[i],coords$X[i],sep=","),"&radius=",radius,
                "&type=parking&key=AIzaSyDRWYQO3gTUNXOQa0t1FerdyOIyfUU8vKQ",sep=""))

myParking<- fromJSON(content(res,"text"),flatten=T)$results
if (is.null(dim(myParking))==FALSE){
parking<-rbind(parking,myParking[,c("place_id","geometry.location.lat","geometry.location.lng")])
}}


d<-duplicated(parking$place_id)
parking<-parking[-which(d=="TRUE"),]

details<-array(NA,c(length(parking$place_id),3))


for (i in 1:length(parking$place_id)){
res<-GET( paste("https://maps.googleapis.com/maps/api/place/details/json?placeid=",parking$place_id[i], "&key=AIzaSyDRWYQO3gTUNXOQa0t1FerdyOIyfUU8vKQ",sep=""))
myDetails<- fromJSON(content(res,"text"),flatten=T)$result
if (is.null(myDetails$name)==FALSE){
details[i,1]<-myDetails$name}
if(is.null(myDetails$price_level)==FALSE){
details[i,2]<-myDetails$price_level}
if(is.null(myDetails$rating)==FALSE){
details[i,3]<-myDetails$rating}
}

details_df<-data.frame(Name=details[,1],Price=details[,2],Rating=details[,3])

parking_det<-cbind(parking,details_df)


write.table(parking_det,"H:\\National Scale Data\\ChapterIII\\ParkingDetail.txt",col.names=T,row.names=F,sep="\t")



####################
#Delete duplicates hotels

hotels<-read.table("HotelsDetail.txt",header=T,stringsAsFactors=F)
hotels_dupl<-duplicated(hotels[,"place_id"])
hotels_unique<-hotels[-which(hotels_dupl=="TRUE"),]
write.table(hotels_unique,"HotelsDetail.txt", row.names=F,sep="\t", quote=F)


#########delete duplicates Tour Operators

TO<-read.table("TODetail.txt",header=T)

d<-duplicated(TO$place_id)
TOnew<-TO[-which(d=="TRUE"),]

TOnew$Type<-as.character(TOnew$Type)

TOclean<-TOnew[-grep("[[:<:]]lodging[[:>:]]",TOnew$Type,value=FALSE,perl=TRUE),]

write.table(TOclean,"TODetail.txt",col.names=T,row.names=F,sep="\t")


######################
#Delete duplicates airports

airports<-read.table("AirportsDetail.txt",header=T,stringsAsFactors=F)
airports_dupl<-duplicated(airports[,"place_id"])
airports_unique<-airports[-which(airports_dupl=="TRUE"),]
write.table(airports_unique,"AirportsDetail.txt", row.names=F,sep="\t", quote=F)


######################
#Delete duplicates bus stations

bus<-read.table("BusStationsDetail.txt",header=T,stringsAsFactors=F)
bus_dupl<-duplicated(bus[,"place_id"])
bus_unique<-bus[-which(bus_dupl=="TRUE"),]
write.table(bus_unique,"BusStationsDetail.txt", row.names=F,sep="\t", quote=F)

######################
#Delete duplicates subway stations

subway<-read.table("SubwayStationsDetail.txt",header=T,stringsAsFactors=F)
subway_dupl<-duplicated(subway[,"place_id"])
subway_unique<-subway[-which(subway_dupl=="TRUE"),]
write.table(subway_unique,"SubwayStationsDetail.txt", row.names=F,sep="\t", quote=F)

######################
#Delete duplicates train stations

train<-read.table("TrainStationsDetail.txt",header=T,stringsAsFactors=F)
train_dupl<-duplicated(train[,"place_id"])
table(train_dupl)
#train_dupl
#FALSE 
#  371 

#no duplicates

######################
#Delete duplicates car parks

parking<-read.table("ParkingDetail.txt",header=T,stringsAsFactors=F)
parking_dupl<-duplicated(parking[,"place_id"])
table(parking_dupl)
#parking_dupl
#FALSE 
#  692 

#no duplicates



 
 