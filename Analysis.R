# Script for statistical analysis of "Modelling the distribution of wildlife watchers using social media"
# Francesca Mancini
# last modified 06/06/2017

# Packages ######
library(gstat)
library(MASS)
library(sp)

# Linear models #####

data_sub <- read.table(".//data//CombinedData_v11.txt", stringsAsFactors = F, header=T)


#check for collinearity between linear predictors
#first create a dataframe containing only the predictors of interest
preds_env <- data_sub[,c("Area_WHS", "Area_SSSI", "Area_SPA", "Area_SAC_L", "Area_RAMSA", "Area_NR", "Area_NNR",
                   "Dist_MPA", "Dist_MCA", "Area_LNR", "Area_COUNE", "Area_CNTRY", "Area_BIOSP", "Area_BIOGE",
                   "Area_NP", "Mean_Nat", "Dist_MSAC")]

preds_inf <- data_sub[,c("Dist_Air", "Count_Bus", "Count_Hotel", "Dist_CarPark",
                         "Dist_TourOp", "Dist_Train", "Dist_road")]

# visually inspect correlations
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(preds_env, lower.panel = panel.cor, upper.panel = panel.smooth)

pairs(preds_inf, lower.panel = panel.cor, upper.panel = panel.smooth)


#then calculate variance inflation factor

#calculate VIF
source("Collinearity.R")
VIF_env <- corvif(preds_env)
# Area_SSSI and Area_SAC_L are highly collinear (VIF > 3)

preds_env <- preds_env[,-4]

VIF_env <- corvif(preds_env)
# all VIF < 2

VIF_inf <- corvif(preds_inf)
# all VIF < 2


#fit the full model to the logged count data with scaled predictors

#transform coordinates from latlong to utm to avoid issues with correlation structures

library(rgdal)

coordinates(data_sub) <- data_sub[,c(29,30)]
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84

proj4string(data_sub) <- crs.geo                             # assign the coordinate system

coords_proj <- spTransform(data_sub, CRS("+proj=utm +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "))

data_sub<-cbind(data_sub@data, coords_proj@coords)

names(data_sub)[c(51,52)] <- c("x","y")


full.lm<-lm(log10(Count_WW+1)~scale(Area_WHS,scale=T) + scale(Area_SPA,scale=T) + 
              scale(Area_SSSI,scale=T) +scale(Area_RAMSA,scale=T) + scale(Area_NR,scale=T) +
              scale(Area_NNR,scale=T) + scale(Dist_MPA,scale=T) +scale(Dist_MSAC,scale=T) +
              scale(Dist_MCA,scale=T) + scale(Area_LNR,scale=T) + scale(Area_COUNE,scale=T) +
              scale(Area_CNTRY,scale=T) + scale(Area_BIOSP,scale=T) +
              scale(Area_BIOGE,scale=T) + scale(Area_NP,scale=T) + scale(Dist_Air,scale=T) + 
              scale(Count_Bus,scale=T) + scale(Count_Hotel,scale=T) + scale(Dist_CarPark,scale=T) + 
              scale(Dist_TourOp,scale=T) +  scale(Dist_Train,scale=T)  + scale(Dist_road,scale=T) + 
              scale(Mean_Nat,scale=T),data=data_sub)


#look at the summary
summary(full.lm)

#and check the residuals
par(mfrow=c(2,2))
plot(full.lm)

#extract the standardised residuals
S.res.lm<-rstandard(full.lm)

#create a spatial dataframe
mydata_sp<-data_sub
coordinates(mydata_sp)<-c("x", "y")

#calculate and plot variograms
Vario<-variogram(S.res.lm~1,data=mydata_sp)
Variodir<-variogram(S.res.lm~1,data=mydata_sp,alpha=c(0,45,90,135))

plot(Vario)
plot(Variodir)

#make a bubble plot of the residuals to check for spatial patterns
bubble.data<-data.frame(S.res.lm,data_sub$x,data_sub$y)
coordinates(bubble.data)<-c("data_sub.x","data_sub.y")

bubble(bubble.data,"S.res.lm",col=c("black","grey"),main="Residuals",xlab="Longitude",ylab="Latitude")


# obvious spatial autocorrelation in residuals
# Full GLS and selection of autocorrelation structure ######


#load required library
library(nlme)

#create a new dataframe with scaled predictors and logged response
gls.data<-data.frame(Count_WW=log10(data_sub$Count_WW+1),Area_WHS=scale(data_sub$Area_WHS,scale=T),
                     Area_SPA=scale(data_sub$Area_SPA,scale=T), Dist_MSAC=scale(data_sub$Dist_MSAC,scale=T),
                    Area_SSSI=scale(data_sub$Area_SSSI,scale=T), Area_RAMSA=scale(data_sub$Area_RAMSA,scale=T),
                     Area_NR=scale(data_sub$Area_NR,scale=T), Area_NNR=scale(data_sub$Area_NNR,scale=T), 
                     Dist_MPA=scale(data_sub$Dist_MPA,scale=T), Dist_MCA=scale(data_sub$Dist_MCA,scale=T), 
                     Area_LNR=scale(data_sub$Area_LNR,scale=T), Area_COUNE=scale(data_sub$Area_COUNE,scale=T),
                     Area_CNTRY=scale(data_sub$Area_CNTRY,scale=T), Area_BIOSP=scale(data_sub$Area_BIOSP,scale=T), 
                     Area_BIOGE=scale(data_sub$Area_BIOGE,scale=T), Area_NP=scale(data_sub$Area_NP,scale=T), 
                     Dist_Air=scale(data_sub$Dist_Air,scale=T), Count_Bus=scale(data_sub$Count_Bus,scale=T),
                     Count_Hotel=log10(data_sub$Count_Hotel+1), Dist_CarPark=scale(data_sub$Dist_CarPark,scale=T), 
                     Dist_TourOp=scale(data_sub$Dist_TourOp,scale=T), Dist_Train=scale(data_sub$Dist_Train,scale=T), 
                     Dist_road=scale(data_sub$Dist_road,scale=T), Mean_Nat=scale(data_sub$Mean_Nat,scale=T),
                     Area_PA = scale(data_sub$Area_PA, scale = T), Count_Inf = scale(data_sub$Count_Inf,scale = T),
                     x=data_sub$x,y=data_sub$y, Pop_dens = data_sub$Pop_dens, 
                     Area_SAC_L = scale(data_sub$Area_SAC_L, scale = T))



# Aggregated variables #######

preds.agg<-data_sub[,c("Area_PA","Mean_Nat","Count_Inf")]

#then calculate variance inflation factor

#calculate VIF
source("Collinearity")
VIF<-corvif(preds.agg)

agg.gls.ratio<-gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf,data=gls.data,
             correlation=corRatio(form=~x+y,nugget=T))

agg.gls.sph<-gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf,data=gls.data,
                   correlation=corSpher(form=~x+y,nugget=T))

agg.gls.Lin<-gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf,data=gls.data,
                 correlation=corLin(form=~x+y,nugget=T))

agg.gls.Gaus<-gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf,data=gls.data,
                 correlation=corGaus(form=~x+y,nugget=T))

agg.gls.exp<-gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf,data=gls.data,
                  correlation=corExp(form=~x+y,nugget=T))


AIC(agg.gls.ratio, agg.gls.sph, agg.gls.Lin, agg.gls.Gaus, agg.gls.exp)

#Exponential autocrrelation structure is the best one

summary(agg.gls.exp)


#check for patterns in the residuals

res.agg.gls<-residuals(agg.gls.exp,type="normalized")
fit.agg.gls<-fitted(agg.gls.exp)

plot(res.agg.gls~fit.agg.gls)

qqnorm(res.agg.gls)
qqline(res.agg.gls)

#check that autocorrelation is not an issue anymore
Vario.agg.gls<-Variogram(agg.gls.exp,form= ~x+y,robust=T,resType = "normalized")

plot(Vario.agg.gls, smooth=F)


# fit an environmental and an infrastructure model
# to select important variables

# Environmental infrastructure model #######

# Select best correlation structure

env.gls <- gls(Count_WW~ Area_WHS + Area_SAC_L + Area_SSSI + Dist_MSAC +Dist_MPA +Dist_MCA +
                 Area_SPA +Area_RAMSA+ Area_NR + Area_NNR+   Area_LNR + Area_COUNE + 
                 Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat,data=gls.data)

env.gls.exp<-gls(Count_WW~ Area_WHS + Area_SAC_L + Area_SSSI + Dist_MSAC +Dist_MPA +Dist_MCA +
               Area_SPA +Area_RAMSA+ Area_NR + Area_NNR+   Area_LNR + Area_COUNE + 
               Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat,
               data=gls.data, correlation=corExp(form=~x+y,nugget=T))

env.gls.ratio<-gls(Count_WW~ Area_WHS + Area_SAC_L + Area_SSSI + Dist_MSAC +Dist_MPA +Dist_MCA +
               Area_SPA +Area_RAMSA+ Area_NR + Area_NNR+   Area_LNR + Area_COUNE + 
               Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat,
               data=gls.data, correlation=corRatio(form=~x+y,nugget=T))

env.gls.sph<-gls(Count_WW~ Area_WHS + Area_SAC_L + Area_SSSI + Dist_MSAC +Dist_MPA +Dist_MCA +
                 Area_SPA +Area_RAMSA+ Area_NR + Area_NNR+   Area_LNR + Area_COUNE + 
                 Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat,
                 data=gls.data, correlation=corSpher(form=~x+y,nugget=T))

env.gls.Lin<-gls(Count_WW~ Area_WHS + Area_SAC_L + Area_SSSI + Dist_MSAC +Dist_MPA +Dist_MCA +
                 Area_SPA +Area_RAMSA+ Area_NR + Area_NNR+   Area_LNR + Area_COUNE + 
                 Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat,
                 data=gls.data,correlation=corLin(form=~x+y,nugget=T))

env.gls.Gaus<-gls(Count_WW~ Area_WHS + Area_SAC_L + Area_SSSI + Dist_MSAC +Dist_MPA +Dist_MCA +
                  Area_SPA +Area_RAMSA+ Area_NR + Area_NNR+   Area_LNR + Area_COUNE + 
                  Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat,
                  data=gls.data,correlation=corGaus(form=~x+y,nugget=T))

AIC(env.gls,env.gls.exp, env.gls.ratio, env.gls.sph, env.gls.Lin, env.gls.Gaus)


#check for patterns in the residuals

res.env.gls<-residuals(env.gls.exp,type="normalized")
fit.env.gls<-fitted(env.gls.exp)

plot(res.env.gls~fit.env.gls)

qqnorm(res.env.gls)
qqline(res.env.gls)

#check that autocorrelation is not an issue anymore
Vario.env.gls<-Variogram(env.gls.exp,form= ~x+y,robust=T,resType = "normalized")

plot(Vario.env.gls, smooth=F)

summary(env.gls.exp)

# Infrastructure model ######

#check for collinearity between linear predictors
#first create a dataframe containing only the predictors of interest

preds.inf<-data_sub[,c("Dist_Air","Count_Bus","Count_Hotel","Dist_CarPark","Dist_TourOp",
                   "Dist_Train","Dist_road")]

#then calculate variance inflation factor

#calculate VIF
source("Collinearity.R")
VIF<-corvif(preds.inf)

inf.gls <- gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                 Dist_Train  + Dist_road,data=gls.data)

inf.gls.exp<-gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                 Dist_Train  + Dist_road,data=gls.data, correlation=corExp(form=~x+y,nugget=T))

inf.gls.sph<-gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                 Dist_Train  + Dist_road,data=gls.data, correlation=corSpher(form=~x+y,nugget=T))

inf.gls.ratio<-gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                   Dist_Train  + Dist_road,data=gls.data, correlation=corRatio(form=~x+y,nugget=T))

inf.gls.Lin<-gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                 Dist_Train  + Dist_road,data=gls.data, correlation=corLin(form=~x+y,nugget=T))

inf.gls.Gaus<-gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                  Dist_Train  + Dist_road,data=gls.data, correlation=corGaus(form=~x+y,nugget=T))

AIC(inf.gls,inf.gls.exp, inf.gls.sph, inf.gls.ratio, inf.gls.Lin, inf.gls.Gaus)


#check for patterns in the residuals

res.inf.gls<-residuals(inf.gls.exp,type="normalized")
fit.inf.gls<-fitted(inf.gls.exp)

plot(res.inf.gls~fit.inf.gls)

qqnorm(res.inf.gls)
qqline(res.inf.gls)

#check that autocorrelation is not an issue anymore
Vario.inf.gls<-Variogram(inf.gls.exp,form= ~x+y,robust=T,resType = "normalized")

plot(Vario.inf.gls, smooth=F)

summary(inf.gls.exp)

## Variable selection ######

#now use dredge to find best combination of variables for both env and infr models

#use pdredge to use parallell computing
library(parallel)
library(MuMIn)

# Calculate the number of cores
cores <- detectCores()

# Determine cluster type (mine is a PSOCK)
#clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"

# Set up a cluster with number of cores specified as result of detectCores() 
#   and call it "clust" 
# For laptop with 4 cores
clust <- makeCluster(cores)#, type = clusterType)


# Load required packages onto worker nodes
#   (in this example, load packages {MASS} and {MuMIn} to be used by pdredg)
clusterEvalQ(clust,library(nlme))
clusterEvalQ(clust,library(MuMIn))
clusterExport(clust,"gls.data")  #export the dataframe to the cluster
clusterExport(clust,"inf.gls.exp")  #export the model to the cluster

#model selection for infrastructure model
inf.sel<-pdredge(inf.gls.exp, cluster=clust,rank = "AICc",trace=2, REML = FALSE)

saveRDS(inf.sel,"Infrastructure_sel.rds")

stopCluster(clust)

#model selection for environmental model
clust <- makeCluster(cores)#, type = clusterType)


# Load required packages onto worker nodes
#   (in this example, load packages {MASS} and {MuMIn} to be used by pdredg)
clusterEvalQ(clust,library(nlme))
clusterEvalQ(clust,library(MuMIn))
clusterExport(clust,"gls.data")  #export the dataframe to the cluster
clusterExport(clust,"env.gls.exp")  #export the model to the cluster

env.sel<-pdredge(env.gls.exp, cluster=clust, rank = "AICc",trace=2, REML = FALSE,
                 subset = !("Area_SSSI"&& "Area_SAC_L")) #do not include collinear variables in the same model
str(env.sel)

saveRDS(env.sel,"Environment_sel.rds")

stopCluster(clust)

inf.sel<-readRDS("Infrastructure_sel.rds")
inf.sel



#refit subset of models with REML
inf.sel.REML<-get.models(inf.sel, subset = delta < 5)

saveRDS(inf.sel.REML,"Infrastructure_Best.rds")

inf.sel.REML<-readRDS("Infrastructure_Best.rds")

# calculate variable importance 

inf.var.imp<-importance(inf.sel.REML)

# put it into a dataframe format

df <- as.data.frame(inf.var.imp)

inf.var.imp <- cbind(df, attr(inf.var.imp, "names"))

names(inf.var.imp) <- c("Importance", "Var")

# and plot

library(ggplot2)

InfVarImpVis <- ggplot (inf.var.imp, aes(Var, Importance, fill = Importance)) +
  geom_hline(yintercept = seq(0, 1.2, by = 0.5), colour = "grey90", size = 1) +
  geom_vline(aes(xintercept = Var), colour = "grey90", size = 1) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  scale_y_continuous(breaks = 0:nlevels(as.factor(inf.var.imp$Var))) +
  scale_fill_gradient(low = "thistle1", high = "thistle4") +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))

InfVarImpVis +coord_polar()


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


# Biodiversity #####

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


####Effect of biodiversity 

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

#transform coordinates from latlong to utm to avoid issues with correlation structures

library(rgdal)

coordinates(data.gr2) <- data.gr2[,c(29,30)]
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84

proj4string(data.gr2) <- crs.geo                             # assign the coordinate system

coords_proj <- spTransform(data.gr2, CRS("+proj=utm +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "))

data.gr2<-cbind(data.gr2@data[,-52], coords_proj@coords)

str(data.gr2)
names(data.gr2)[c(52,53)] <- c("x","y")


##model effect of biodiversity

bio.1<-lm(log10(Count_WW+1)~log10(Species) +offset(log10(Pop_dens+1)),data=data.gr2)
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
                         Pop_dens = log10(data.gr2$Pop_dens + 1), x=data.gr2$x,y=data.gr2$y)


bio.gls.exp<-gls(Count_WW ~ Species + offset(Pop_dens), 
                 data = data.bio.gls, correlation = corExp(form = ~ x+y,nugget = T))

bio.gls.exp2<-gls(Count_WW ~ Species + Pop_dens,
                 data = data.bio.gls, correlation = corExp(form = ~ x+y,nugget = T))


bio.gls.lin<-gls(Count_WW ~ Species + offset(Pop_dens), 
                 data = data.bio.gls, correlation = corLin(form = ~ x+y,nugget = T))

bio.gls.ratio<-gls(Count_WW ~ Species + offset(Pop_dens), 
                 data = data.bio.gls, correlation = corRatio(form = ~ x+y,nugget = T))

bio.gls.gaus<-gls(Count_WW ~ Species + offset(Pop_dens), 
                 data = data.bio.gls, correlation = corGaus(form = ~ x+y,nugget = T))

bio.gls.sph<-gls(Count_WW ~ Species + offset(Pop_dens), 
                  data = data.bio.gls, correlation = corSpher(form = ~ x+y,nugget = T))

AIC(bio.gls.exp, bio.gls.lin, bio.gls.ratio, bio.gls.gaus, bio.gls.sph)

summary(bio.gls.exp)

res.bio.gls<-residuals(bio.gls.exp2, type="normalized")

fit.bio.gls<-fitted(bio.gls.exp2)

par(mfrow=c(1,2))
plot(res.bio.gls ~ fit.bio.gls)

qqnorm(res.bio.gls)

qqline(res.bio.gls)


#calculate and plot variograms
Vario<-variogram(res.bio.gls~1,data=mydata_sp)

plot(Vario)

#####calculate and plot predictions

pred_data<-data.frame(Species=seq(from=min(data.bio.gls$Species), 
                                  to=max(data.bio.gls$Species),by=0.001),
                      Pop_dens = rep(max(data.bio.gls$Pop_dens), n = 2604))

preds<-predict(bio.gls.exp2,pred_data,type="response")

plot(Count_WW~Species,data=data.bio.gls)
lines(preds~pred_data$Species)

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


