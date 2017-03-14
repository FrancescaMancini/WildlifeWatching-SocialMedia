########################################
########Combine data
########Francesca Mancini
########last modified 14/03/2017
########################################

#load the column for distance from marine SAC from the subfolder //data
Dist_MSAC<-read.table(".//data//CombinedData_v10.txt",stringsAsFactors = F,header=T)[,43]

#load columns for distance from marine SAC, MPA and MCA and area of terrestrial SAC
new<-read.table(".//data//Area_SAC_L.txt",stringsAsFactors = F,header=T)

#load the rest of the data
data<-read.table(".//data//CombinedData_v9.txt",stringsAsFactors = F,header=T)

#put data together
data_new<-cbind(data,new,Dist_MSAC)

#include biodiversity records and species richness
bio<-read.table(".//data//SpeciesRichness.txt",header=T)
bio_data<-cbind(data_new,bio[,2:3])

str(bio_data)

#subset the data to exclude all the observations without natural values
data_sub<-bio_data[is.na(bio_data$Mean_Nat)==F,]

write.table(data_sub,".//data//CombinedData_v11.txt", sep="\t", row.names=F)




