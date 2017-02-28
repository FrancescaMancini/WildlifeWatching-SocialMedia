########################################
########Combine data
########Francesca Mancini
########last modified 28/02/2017
########################################

#load the column for distance from marine SAC from the subfolder //data
Dist_MSAC<-read.table(".//data//CombinedData_v10.txt",stringsAsFactors = F,header=T)[,43]

#load columns for distance from marine SAC, MPA and MCA and area of terrestrial SAC
new<-read.table(".//data//Area_SAC_L.txt",stringsAsFactors = F,header=T)

#load the rest of the data
data<-read.table(".//data//CombinedData_v9.txt",stringsAsFactors = F,header=T)

#put data together
data_new<-cbind(data,new,Dist_MSAC)

#subset the data to exclude all the observations without natural values
data_sub<-data_new[is.na(data_new$Mean_Nat)==F,]

write.table(data_sub,".//data//CombinedData_v11.txt", sep="\t", row.names=F)




