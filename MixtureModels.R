#########################################
####Mixture model to classify data:
####low vs high number of records
####David Lusseau & Francesca Mancini
####last modified 20/03/2017
#########################################

library(mclust)

#create variable log of number of records
data_sub$logrec<-log(data_sub$Records+1)

#Produces a density estimate for each data point using a Gaussian finite mixture model from Mclust.
clusrec<-densityMclust(data_sub$logrec,G=2)     #G=2 we know that there are 2 groups (low and high)

#plot the distribution
plot(clusrec, what = "density", data = data_sub$logrec, breaks = 15)

#Very high nuber of 0s 

#remove 0 records
sublogrec<-data_sub$logrec[data_sub$logrec>0]

#refit mixture model without 0s
clusrecsub<-densityMclust(sublogrec,G=2)

#plot
plot(clusrecsub, what = "density", data = sublogrec, breaks = 15)

#returns a summary of the model
classification<- summary(clusrecsub,classification=T)

#classification of each datapoint into the two groups is then in:
membership<-classification$classification

table(membership)








