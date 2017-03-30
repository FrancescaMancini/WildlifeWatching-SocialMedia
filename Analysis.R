####################################################
###########Francesca Mancini
###########last modified 14/03/2017
####################################################

#load required libraries
library(gstat)
library(MASS)
library(sp)

data_sub<-read.table(".//data//CombinedData_v11.txt",stringsAsFactors = F,header=T)


#check for collinearity between linear predictors
#first create a dataframe containing only the predictors of interest
preds<-data_sub[,c("Area_WHS","Area_SSSI","Area_SPA","Area_SAC_L","Area_RAMSA","Area_NR","Area_NNR",
                   "Dist_MPA","Dist_MCA","Area_LNR","Area_COUNE","Area_CNTRY","Area_BIOSP","Area_BIOGE",
                   "Area_NP","Dist_Air","Count_Bus","Count_Hotel","Dist_CarPark","Dist_TourOp",
                   "Dist_Train","Dist_road","Mean_Nat","Dist_MSAC")]

#then calculate variance inflation factor

#calculate VIF
source("Collinearity")
VIF<-corvif(preds)

#fit the full model to the logged count data with scaled predictors

full.lm<-lm(log10(Count_WW+1)~scale(Area_WHS,scale=T) + scale(Area_SPA,scale=T) + 
              scale(Area_SAC_L,scale=T) +scale(Area_RAMSA,scale=T) + scale(Area_NR,scale=T) +
              scale(Area_NNR,scale=T) + scale(Dist_MPA,scale=T) +scale(Dist_MSAC,scale=T) +
              scale(Dist_MCA,scale=T) + scale(Area_LNR,scale=T) + scale(Area_COUNE,scale=T) +
              scale(Area_CNTRY,scale=T) + scale(Area_BIOSP,scale=T) +
              scale(Area_BIOGE,scale=T) + scale(Area_NP,scale=T) + scale(Dist_Air,scale=T) + 
              scale(Count_Bus,scale=T) + scale(Count_Hotel,scale=T) + scale(Dist_CarPark,scale=T) + 
              scale(Dist_TourOp,scale=T) +  scale(Dist_Train,scale=T)  + scale(Dist_road,scale=T) + 
              scale(Mean_Nat,scale=T) + offset(log10(Pop_dens+1)),data=data_sub)


#look at the summary
summary(full.lm)

#and check the residuals
par(mfrow=c(2,2))
plot(full.lm)

#extract the standardised residuals
S.res.lm<-rstandard(full.lm)

#create a spatial dataframe
mydata_sp<-data_sub
coordinates(mydata_sp)<-c("Longitude", "Latitude")

#calculate and plot variograms
Vario<-variogram(S.res.lm~1,data=mydata_sp)
Variodir<-variogram(S.res.lm~1,data=mydata_sp,alpha=c(0,45,90,135))

plot(Vario)
plot(Variodir)

#make a bubble plot of the residuals to check for spatial patterns
bubble.data<-data.frame(S.res.lm,data_sub$Longitude,data_sub$Latitude)
coordinates(bubble.data)<-c("data_sub.Longitude","data_sub.Latitude")

bubble(bubble.data,"S.res.lm",col=c("black","grey"),main="Residuals",xlab="Longitude",ylab="Latitude")

##########################################
#######fix autocorrelation with GLS
##########################################

#load required library
library(nlme)

#create a new dataframe with scaled predictors and logged response
gls.data<-data.frame(Count_WW=log10(data_sub$Count_WW+1),Area_WHS=scale(data_sub$Area_WHS,scale=T),
                     Area_SPA=scale(data_sub$Area_SPA,scale=T), Dist_MSAC=scale(data_sub$Dist_MSAC,scale=T),
                     Area_SAC_L=scale(data_sub$Area_SAC_L,scale=T), Area_RAMSA=scale(data_sub$Area_RAMSA,scale=T),
                     Area_NR=scale(data_sub$Area_NR,scale=T), Area_NNR=scale(data_sub$Area_NNR,scale=T), 
                     Dist_MPA=scale(data_sub$Dist_MPA,scale=T), Dist_MCA=scale(data_sub$Dist_MCA,scale=T), 
                     Area_LNR=scale(data_sub$Area_LNR,scale=T), Area_COUNE=scale(data_sub$Area_COUNE,scale=T),
                     Area_CNTRY=scale(data_sub$Area_CNTRY,scale=T), Area_BIOSP=scale(data_sub$Area_BIOSP,scale=T), 
                     Area_BIOGE=scale(data_sub$Area_BIOGE,scale=T), Area_NP=scale(data_sub$Area_NP,scale=T), 
                     Dist_Air=scale(data_sub$Dist_Air,scale=T), Count_Bus=scale(data_sub$Count_Bus,scale=T),
                     Count_Hotel=log10(data_sub$Count_Hotel+1), Dist_CarPark=scale(data_sub$Dist_CarPark,scale=T), 
                     Dist_TourOp=scale(data_sub$Dist_TourOp,scale=T), Dist_Train=scale(data_sub$Dist_Train,scale=T), 
                     Dist_road=scale(data_sub$Dist_road,scale=T), Mean_Nat=scale(data_sub$Mean_Nat,scale=T),
                     Pop_dens=log10(data_sub$Pop_dens+1),x=data_sub$Longitude,y=data_sub$Latitude)

#transform coordinates from latlong to utm to avoid issues with correlation structures

library(rgdal)

coordinates(gls.data) <- gls.data[,c(26,27)]
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84

proj4string(gls.data) <- crs.geo                             # assign the coordinate system

coords_proj <- spTransform(gls.data, CRS("+proj=utm +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "))

gls.data<-cbind(gls.data@data, coords_proj@coords)

names(gls.data)[c(26,27)] <- c("Long","Lat")

#fit the full model with gls
full.gls<-gls(Count_WW~ Area_WHS +  Area_SPA + Dist_MSAC +Area_SAC_L + Area_RAMSA+ Area_NR +
                Area_NNR+ Dist_MPA + Dist_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                Dist_Train  + Dist_road + Mean_Nat + offset(Pop_dens),data=gls.data)

summary(full.gls)

#check for autocorrelation
Vario.gls<- Variogram(full.gls,form= ~x+y,robust=T,resType = "pearson",maxDist=1)

plot(Vario.gls,smooth=T)

#fit the same model with different autocorrelation structures
full.gls.sph<-gls(Count_WW~ Area_WHS +  Area_SPA + Dist_MSAC +Area_SAC_L + Area_RAMSA+ Area_NR +
                    Area_NNR+ Dist_MPA + Dist_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                    Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                    Dist_Train  + Dist_road + Mean_Nat + offset(Pop_dens),data=gls.data,
                  correlation=corSpher(form=~x+y,nugget=T))

full.gls.Lin<-gls(Count_WW~ Area_WHS +  Area_SPA + Dist_MSAC +Area_SAC_L + Area_RAMSA+ Area_NR +
                    Area_NNR+ Dist_MPA + Dist_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                    Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                    Dist_Train  + Dist_road + Mean_Nat + offset(Pop_dens),data=gls.data,
                  correlation=corLin(form=~x+y,nugget=T))

full.gls.ratio<-gls(Count_WW~ Area_WHS +  Area_SPA + Dist_MSAC +Area_SAC_L + Area_RAMSA+ Area_NR +
                      Area_NNR+ Dist_MPA + Dist_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                      Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                      Dist_Train  + Dist_road + Mean_Nat + offset(Pop_dens),data=gls.data,
                    correlation=corRatio(form=~x+y,nugget=T))

full.gls.Gaus<-gls(Count_WW~ Area_WHS +  Area_SPA + Dist_MSAC +Area_SAC_L + Area_RAMSA+ Area_NR +
                     Area_NNR+ Dist_MPA + Dist_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                     Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                     Dist_Train  + Dist_road + Mean_Nat + offset(Pop_dens),data=gls.data,
                   correlation=corGaus(form=~x+y,nugget=T))

full.gls.exp<-gls(Count_WW~ Area_WHS +  Area_SPA + Dist_MSAC +Area_SAC_L + Area_RAMSA+ Area_NR +
                    Area_NNR+ Dist_MPA + Dist_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                    Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                    Dist_Train  + Dist_road + Mean_Nat + offset(Pop_dens),data=gls.data,
                  correlation=corExp(form=~x+y,nugget=T))

#use AIC to choose the best one
AIC(full.gls,full.gls.sph,full.gls.Lin,full.gls.ratio,full.gls.Gaus,full.gls.exp)

#model with ratio correlation structure is the best model



###########################################
##model selection

#Variables Area_LSAC and Area_SSSI are collinear 
#put each in a different model
#use AIC to pick best model and do model validation

gls.data<-data.frame(Count_WW=log10(data_sub$Count_WW+1),Area_WHS=scale(data_sub$Area_WHS,scale=T),
                       Dist_MSAC=scale(data_sub$Dist_MSAC,scale=T), Area_SPA=scale(data_sub$Area_SPA,scale=T),
                       Area_LSAC=scale(data_sub$Area_SAC_L,scale=T), Area_RAMSA=scale(data_sub$Area_RAMSA,scale=T),
                       Area_NR=scale(data_sub$Area_NR,scale=T), Area_NNR=scale(data_sub$Area_NNR,scale=T),
                       Dist_MPA=scale(data_sub$Dist_MPA,scale=T), Dist_MCA=scale(data_sub$Dist_MCA,scale=T),
                       Area_LNR=scale(data_sub$Area_LNR,scale=T), Area_COUNE=scale(data_sub$Area_COUNE,scale=T),
                       Area_CNTRY=scale(data_sub$Area_CNTRY,scale=T), Area_BIOSP=scale(data_sub$Area_BIOSP,scale=T),
                       Area_BIOGE=scale(data_sub$Area_BIOGE,scale=T), Area_NP=scale(data_sub$Area_NP,scale=T),
                       Dist_Air=scale(data_sub$Dist_Air,scale=T), Count_Bus=scale(data_sub$Count_Bus,scale=T),
                       Count_Hotel=log10(data_sub$Count_Hotel+1), Dist_CarPark=scale(data_sub$Dist_CarPark,scale=T),
                       Dist_TourOp=scale(data_sub$Dist_TourOp,scale=T), Dist_Train=scale(data_sub$Dist_Train,scale=T),
                       Dist_road=scale(data_sub$Dist_road,scale=T), Mean_Nat=scale(data_sub$Mean_Nat,scale=T),
                       Area_PA=scale(data_sub$Area_PA,scale=T),Count_Inf=log10(data_sub$Count_Inf+1),
                       Dist_MSAC=scale(data_sub$Dist_MSAC,scale=T),Area_SSSI=scale(data_sub$Area_SSSI,scale=T),
                       Pop_dens=log10(data_sub$Pop_dens+1),x=gls.data$x,y=gls.data$y)



#fit the full model
full.gls.exp.2<-gls(Count_WW~ Area_WHS + Area_LSAC +  Dist_MSAC +Dist_MPA +Dist_MCA +
                    Area_SPA +Area_RAMSA+ Area_NR + Area_NNR+   Area_LNR + Area_COUNE + 
                    Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Dist_Air + Count_Bus + 
                    Count_Hotel + Dist_CarPark + Dist_TourOp + Dist_Train  + Dist_road + 
                    Mean_Nat + offset(Pop_dens),data=gls.data,
                    correlation=corRatio(form=~x+y,nugget=T),method="ML")

full.gls.exp.3<-gls(Count_WW~ Area_WHS + Area_SSSI +  Dist_MSAC +Dist_MPA +Dist_MCA +
                    Area_SPA +Area_RAMSA+ Area_NR + Area_NNR+   Area_LNR + Area_COUNE + 
                    Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Dist_Air + Count_Bus + 
                    Count_Hotel + Dist_CarPark + Dist_TourOp + Dist_Train  + Dist_road + 
                    Mean_Nat + offset(Pop_dens),data=gls.data,
                    correlation=corRatio(form=~x+y,nugget=T),method="ML")

AIC(full.gls.exp.2,full.gls.exp.3)

summary(full.gls.exp.3)


res.gls<-residuals(full.gls.exp.3,type="normalized")
fit.gls<-fitted(full.gls.exp.3)

plot(res.gls~fit.gls)

qqnorm(res.gls)
qqline(res.gls)


#check that autocorrelation is not an issue anymore
Vario.gls.exp<-Variogram(full.gls.exp.3,form= ~x+y,robust=T,resType = "normalized")

plot(Vario.gls.exp, smooth=F)

#fit a model with aggregated variables

preds.agg<-data_sub[,c("Area_PA","Mean_Nat","Count_Inf")]

#then calculate variance inflation factor

#calculate VIF
source("Collinearity")
VIF<-corvif(preds.agg)

agg.gls.ratio<-gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf + offset(Pop_dens),data=gls.data,
             correlation=corRatio(form=~x+y,nugget=T))

agg.gls.sph<-gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf + offset(Pop_dens),data=gls.data,
                   correlation=corSpher(form=~x+y,nugget=T))

agg.gls.Lin<-gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf + offset(Pop_dens),data=gls.data,
                 correlation=corLin(form=~x+y,nugget=T))

agg.gls.Gaus<-gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf + offset(Pop_dens),data=gls.data,
                 correlation=corGaus(form=~x+y,nugget=T))

agg.gls.exp<-gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf + offset(Pop_dens),data=gls.data,
                  correlation=corExp(form=~x+y,nugget=T))

AIC(agg.gls.ratio, agg.gls.sph, agg.gls.Lin, agg.gls.Gaus, agg.gls.exp)

#Spherical autocrrelation structure is the best one

summary(agg.gls.sph)


#check for patterns in the residuals

res.agg.gls<-residuals(agg.gls.sph,type="normalized")
fit.agg.gls<-fitted(agg.gls.sph)

plot(res.agg.gls~fit.agg.gls)

qqnorm(res.agg.gls)
qqline(res.agg.gls)

#check that autocorrelation is not an issue anymore
Vario.agg.gls<-Variogram(agg.gls.sph,form= ~x+y,robust=T,resType = "normalized")

plot(Vario.agg.gls, smooth=F)


# all variables but MeanNat have an effect on the response 
# fit an environmental and an infrastructure model
# to select important variables

#Environment

preds.env<-data_sub[,c("Area_WHS","Area_SPA","Area_SAC_L","Area_SSSI","Area_RAMSA","Area_NR","Area_NNR",
                   "Dist_MPA","Dist_MCA","Area_LNR","Area_COUNE","Area_CNTRY","Area_BIOSP","Area_BIOGE",
                   "Area_NP","Mean_Nat","Dist_MSAC")]

#then calculate variance inflation factor

#calculate VIF
source("Collinearity.R")
VIF<-corvif(preds.env)


env.gls.exp<-gls(Count_WW~ Area_WHS + Area_LSAC + Area_SSSI + Dist_MSAC +Dist_MPA +Dist_MCA +
               Area_SPA +Area_RAMSA+ Area_NR + Area_NNR+   Area_LNR + Area_COUNE + 
               Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat + 
               offset(Pop_dens),data=gls.data, correlation=corExp(form=~x+y,nugget=T))

env.gls.ratio<-gls(Count_WW~ Area_WHS + Area_LSAC + Area_SSSI + Dist_MSAC +Dist_MPA +Dist_MCA +
               Area_SPA +Area_RAMSA+ Area_NR + Area_NNR+   Area_LNR + Area_COUNE + 
               Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat + 
               offset(Pop_dens),data=gls.data, correlation=corRatio(form=~x+y,nugget=T))

env.gls.sph<-gls(Count_WW~ Area_WHS + Area_LSAC + Area_SSSI + Dist_MSAC +Dist_MPA +Dist_MCA +
                 Area_SPA +Area_RAMSA+ Area_NR + Area_NNR+   Area_LNR + Area_COUNE + 
                 Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat + 
                 offset(Pop_dens),data=gls.data, correlation=corSpher(form=~x+y,nugget=T))

env.gls.Lin<-gls(Count_WW~ Area_WHS + Area_LSAC + Area_SSSI + Dist_MSAC +Dist_MPA +Dist_MCA +
                 Area_SPA +Area_RAMSA+ Area_NR + Area_NNR+   Area_LNR + Area_COUNE + 
                 Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat + 
                 offset(Pop_dens),data=gls.data,correlation=corLin(form=~x+y,nugget=T))

env.gls.Gaus<-gls(Count_WW~ Area_WHS + Area_LSAC + Area_SSSI + Dist_MSAC +Dist_MPA +Dist_MCA +
                  Area_SPA +Area_RAMSA+ Area_NR + Area_NNR+   Area_LNR + Area_COUNE + 
                  Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat + 
                  offset(Pop_dens),data=gls.data,correlation=corGaus(form=~x+y,nugget=T))

AIC(env.gls.exp, env.gls.ratio, env.gls.sph, env.gls.Lin, env.gls.Gaus)


#check for patterns in the residuals

res.env.gls<-residuals(env.gls.exp,type="normalized")
fit.env.gls<-fitted(env.gls.exp)

plot(res.env.gls~fit.env.gls)

qqnorm(res.env.gls)
qqline(res.env.gls)

#check that autocorrelation is not an issue anymore
Vario.env.gls<-Variogram(env.gls.exp,form= ~x+y,robust=T,resType = "normalized")

plot(Vario.env.gls, smooth=F)


####Infrastructure

#check for collinearity between linear predictors
#first create a dataframe containing only the predictors of interest
preds.inf<-data_sub[,c("Dist_Air","Count_Bus","Count_Hotel","Dist_CarPark","Dist_TourOp",
                   "Dist_Train","Dist_road")]

#then calculate variance inflation factor

#calculate VIF
source("Collinearity.R")
VIF<-corvif(preds.inf)


inf.gls.exp<-gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                 Dist_Train  + Dist_road +offset(Pop_dens),data=gls.data,
                 correlation=corExp(form=~x+y,nugget=T))

inf.gls.sph<-gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                 Dist_Train  + Dist_road +offset(Pop_dens),data=gls.data,
                 correlation=corSpher(form=~x+y,nugget=T))

inf.gls.ratio<-gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                   Dist_Train  + Dist_road +offset(Pop_dens),data=gls.data,
                   correlation=corRatio(form=~x+y,nugget=T))

inf.gls.Lin<-gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                 Dist_Train  + Dist_road +offset(Pop_dens),data=gls.data,
                 correlation=corLin(form=~x+y,nugget=T))

inf.gls.Gaus<-gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                  Dist_Train  + Dist_road +offset(Pop_dens),data=gls.data,
                  correlation=corGaus(form=~x+y,nugget=T))

AIC(inf.gls.exp, inf.gls.sph, inf.gls.ratio, inf.gls.Lin, inf.gls.Gaus)


#check for patterns in the residuals

res.inf.gls<-residuals(inf.gls.exp,type="normalized")
fit.inf.gls<-fitted(inf.gls.exp)

plot(res.inf.gls~fit.inf.gls)

qqnorm(res.inf.gls)
qqline(res.inf.gls)

#check that autocorrelation is not an issue anymore
Vario.inf.gls<-Variogram(inf.gls.exp,form= ~x+y,robust=T,resType = "normalized")

plot(Vario.inf.gls, smooth=F)


#now use dredge to find best combination of variables for both env and infr models

#use pdredge to use parallell computing
library(parallel)
library(doParallel)
library(MuMIn)

cl <- makeCluster(3)            #split into 3 cores
registerDoParallel(cl)          #register the parallel backend
clusterExport(cl,"gls.data")  #export the dataframe to the cluster
clusterEvalQ(cl,library(nlme))  #load the required package onto the cluster

#model selection for infrastructure model
inf.sel<-pdredge(inf.gls.exp,cluster=cl,rank = "AICc",trace=2, REML = FALSE,
                 fixed =c("offset(Pop_dens)"))   # include offset(Pop_dens) in all models
  

saveRDS(inf.sel,"Infrastructure_sel.rds")

#model selection for environmental model
env.sel<-pdredge(env.gls,cluster=cl,rank = "AICc",trace=2, REML = FALSE,
                 fixed =c("offset(Pop_dens)"),# include offset(Pop_dens) in all models
                 subset["Area_SSSI", "Area_LSAC"] == FALSE) #do not include collinear variables in the same model
str(env.sel)

saveRDS(env.sel,"Environment_sel.rds")

stopCluster(cl)

inf.sel<-readRDS("Infrastructure_sel.rds")
inf.sel

inf.sel.sub<-subset(inf.sel, delta <5)

inf.var.imp<-importance(inf.sel.sub)

#refit subset of models with REML
inf.sel.REML<-get.models(inf.sel, subset = delta < 5, method = "REML")

#model averaging
inf.avg<-model.avg(inf.sel.REML, revised.var = TRUE) 
summary(inf.avg) 

inf.confint<-confint(inf.avg)

#inf.pred.parms<-get.models(inf.sel.sub,subset=T)

newdata.Bus<-data.frame(Count_WW= gls.data$Count_WW,
                        Count_Bus=seq(min(gls.data$Count_Bus),max(gls.data$Count_Bus),
                        length.out=length(gls.data$Count_Bus)), Count_Hotel=rep(mean(gls.data$Count_Hotel),
                        length(gls.data$Count_Bus)),Dist_Air=rep(mean(gls.data$Dist_Air),
                        length(gls.data$Dist_Air)),Dist_Train=rep(mean(gls.data$Dist_Train),
                        length(gls.data$Dist_Train)),Pop_dens=rep(mean(gls.data$Pop_dens),
                        length(gls.data$Pop_dens)),Dist_TourOp=rep(mean(gls.data$Dist_TourOp),
                        length(gls.data$Dist_TourOp)),Dist_CarPark=rep(mean(gls.data$Dist_CarPark),
                        length(gls.data$Dist_CarPark)),Dist_road=rep(mean(gls.data$Dist_road),
                        length(gls.data$Dist_road)))

newdata.Hotel<-data.frame(Count_WW= gls.data$Count_WW,
                          Count_Bus=rep(mean(gls.data$Count_Bus),length(gls.data$Count_Bus)),
                          Count_Hotel=seq(from=min(gls.data$Count_Hotel),to=max(gls.data$Count_Hotel),
                          length.out=length(gls.data$Count_Hotel)),Dist_Air=rep(mean(gls.data$Dist_Air),
                        length(gls.data$Dist_Air)),Dist_Train=rep(mean(gls.data$Dist_Train),
                        length(gls.data$Dist_Train)),Pop_dens=rep(mean(gls.data$Pop_dens),
                        length(gls.data$Pop_dens)),Dist_TourOp=rep(mean(gls.data$Dist_TourOp),
                        length(gls.data$Dist_TourOp)),Dist_CarPark=rep(mean(gls.data$Dist_CarPark),
                        length(gls.data$Dist_CarPark)),Dist_road=rep(mean(gls.data$Dist_road),
                        length(gls.data$Dist_road)))

newdata.Air<-data.frame(Count_WW= gls.data$Count_WW,
                        Count_Bus=rep(mean(gls.data$Count_Bus),length(gls.data$Count_Bus)),
                          Count_Hotel=rep(mean(gls.data$Count_Hotel),length(gls.data$Count_Bus)),
                          Dist_Air=seq(from=min(gls.data$Dist_Air),to=max(gls.data$Dist_Air),
                         length.out=length(gls.data$Dist_Air)),Dist_Train=rep(mean(gls.data$Dist_Train),
                          length(gls.data$Dist_Train)),Pop_dens=rep(mean(gls.data$Pop_dens),
                          length(gls.data$Pop_dens)),Dist_TourOp=rep(mean(gls.data$Dist_TourOp),
                          length(gls.data$Dist_TourOp)),Dist_CarPark=rep(mean(gls.data$Dist_CarPark),
                          length(gls.data$Dist_CarPark)),Dist_road=rep(mean(gls.data$Dist_road),
                          length(gls.data$Dist_road)))


Bus.preds <- sapply(inf.pred.parms, predict, newdata = newdata.Bus) 
Bus.ave4plot<-Bus.preds %*% Weights(inf.sel.sub) 

plot(gls.data$Count_WW~newdata.Bus$Count_Bus)
lines(Bus.ave4plot~newdata.Bus$Count_Bus)

Hotel.preds <- sapply(inf.pred.parms, predict, newdata = newdata.Hotel)
Hotel.ave4plot<-Hotel.preds %*% Weights(inf.sel.sub) 

plot(gls.data$Count_WW~newdata.Hotel$Count_Hotel)
lines(Hotel.ave4plot~newdata.Hotel$Count_Hotel)


Air.preds <- sapply(inf.pred.parms, predict, newdata = newdata.Air)
Air.ave4plot<-Air.preds %*% Weights(inf.sel.sub)

plot(gls.data$Count_WW~newdata.Air$Dist_Air)
lines(Air.ave4plot~newdata.Air$Dist_Air)

######################################
#####Biodiversity

data_sub<-read.table(".//data//CombinedData_v11.txt",stringsAsFactors = F,header=T)

#look at distribution of biodiversity records
par(mfrow=c(2,1))
hist(log(data_sub$Records+1),main="Histogram of species records",xlab="Species Records (log)")
plot(density(log(data_sub$Records+1)),main="Density of species records")

#look at distribution of biodiversity richness
par(mfrow=c(2,1))
hist(log(data_sub$Species+1),main="Histogram of species richness",xlab="Species Richness (log)")
plot(density(log(data_sub$Species+1)),main="Density of species richness")

#look at distribution of species richness where Records are above 7
par(mfrow=c(2,1))
hist(log(data_sub$Species[which(data_sub$Records>7)]),main="Histogram of species richness",xlab="Species Richness (log)")
plot(density(log(data_sub$Species[which(data_sub$Records>7)])),main="Density of species richness")

###############################
####Effect of biodiversity 
###############################

#read classification from mixture model
classes<-read.table("data/classification.txt")

#calculate log of records
data_sub$logrec<-log(data_sub$Records+1)

#only select non 0 records
data.pos<-data_sub[data_sub$logrec>0,]

#now select only records belonging to group 2
data.gr2<-cbind(data.pos,classes)
data.gr2<-data.gr2[data.gr2$x==2,]
str(data.gr2)

##model effect of biodiversity

bio.1<-lm(log10(Count_WW+1)~log10(Species),data=data.gr2)
summary(bio.1)

par(mfrow=c(2,2))
plot(bio.1)

#extract the standardised residuals
S.res.bio.1<-rstandard(bio.1)

#create a spatial dataframe
mydata_sp<-data.gr2
coordinates(mydata_sp)<-c("Longitude", "Latitude")

#calculate and plot variograms
Vario<-variogram(S.res.bio.1~1,data=mydata_sp)
Variodir<-variogram(S.res.bio.1~1,data=mydata_sp,alpha=c(0,45,90,135))

plot(Vario)
plot(Variodir)

#make a bubble plot of the residuals to check for spatial patterns
bubble.data<-data.frame(S.res.bio.1,data.gr2$Longitude,data.gr2$Latitude)
coordinates(bubble.data)<-c("data.gr2.Longitude","data.gr2.Latitude")

bubble(bubble.data,"S.res.bio.1",col=c("black","grey"),main="Residuals",xlab="Longitude",ylab="Latitude")

#####gls

data.bio.gls<-data.frame(Count_WW=log10(data.gr2$Count_WW+1),Species=log10(data.gr2$Species),
                         x=data.gr2$Longitude,y=data.gr2$Latitude)


bio.gls.exp<-gls(Count_WW~ Species,data=data.bio.gls, correlation=corExp(form=~x+y,nugget=T))
summary(bio.gls.exp)

res.bio.gls<-residuals(bio.gls.exp, type="normalized")

fit.bio.gls<-fitted(bio.gls.exp)

par(mfrow=c(1,2))
plot(res.bio.gls ~ fit.bio.gls)

qqnorm(res.bio.gls)

qqline(res.bio.gls)


#calculate and plot variograms
Vario<-variogram(res.bio.gls~1,data=mydata_sp)

plot(Vario)

#####calculate and plot predictions

pred_data<-data.frame(Species=seq(from=min(data.bio.gls$Species), 
                                  to=max(data.bio.gls$Species),by=0.001))

preds<-predict(bio.gls.exp,pred_data,type="link")
preds<-as.data.frame(preds)

plot(Count_WW~Species,data=data.bio.gls)
lines(preds$preds~pred_data$Species)

#not a very good fit
#maybe non linear relationship
#try gam

library(mgcv)

bio.gam<-gam(log10( Count_WW +1 )~s(Species),data=data.gr2)
summary(bio.gam)

plot(bio.gam)
gam.check(bio.gam)

pred_data<-data.frame(Species=seq(from=min(data.gr2$Species), 
                                  to=max(data.gr2$Species),by=0.001))

preds<-predict(bio.gam,pred_data,type="response",se=T)
preds<-as.data.frame(preds)
CIup<-preds$fit + 1.96 *preds$se.fit
CIlow<-preds$fit - 1.96 *preds$se.fit


library(scales)

plot(log10(Count_WW +1)~Species,data=data.gr2, pch=20,col=alpha("cadetblue",0.5),xlab="Number of species",ylab="Number of Flickr users")
lines(preds$fit~pred_data$Species, lwd=3, col="cadetblue")
lines(CIup~pred_data$Species,lty=2,lwd=2,col="cadetblue")
lines(CIlow~pred_data$Species,lty=2,lwd=2,col="cadetblue")


