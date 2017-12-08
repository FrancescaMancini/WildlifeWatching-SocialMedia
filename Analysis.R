# Script for statistical analysis of 
# "“Big Data” for quantifying recreational ecosystem services 
# and characterising people’s connection to nature"
# created by Francesca Mancini
# last modified 28/09/2017

# Packages ######

# load required packages
library(gstat)
library(MASS)
library(sp)
library(rgdal)
library(nlme)
library(MuMIn)
library(parallel)
library(ggplot2)
source("Multiplot.R")
source("nseq.R")
source("Collinearity.R")

# Linear models #####

# load the data
data_sub <- read.table(".//data//CombinedData_v11.txt", stringsAsFactors = F, header=T)

# check for collinearity between linear predictors
# first create a dataframe containing only the predictors of interest
preds_env <- data_sub[,c("Area_WHS", "Area_SSSI", "Area_SPA", "Area_SAC_L", 
                         "Area_RAMSA", "Area_NR", "Area_NNR","Dist_MPA", 
                         "Dist_MCA", "Area_LNR", "Area_COUNE", "Area_CNTRY", 
                         "Area_BIOSP", "Area_BIOGE","Area_NP", "Mean_Nat", "Dist_MSAC")]

preds_inf <- data_sub[,c("Dist_Air", "Count_Bus", "Count_Hotel", "Dist_CarPark",
                         "Dist_TourOp", "Dist_Train", "Dist_road")]

# visually inspect correlations
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456759), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.5/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(preds_env, lower.panel = panel.cor, upper.panel = panel.smooth)

pairs(preds_inf, lower.panel = panel.cor, upper.panel = panel.smooth)


# then calculate variance inflation factor

# calculate VIF
VIF_env <- corvif(preds_env)

# Area_SSSI and Area_SAC_L are highly collinear (VIF > 3)

# calculate VIF for environmental variables without the collinear one
preds_env <- preds_env[,-4]

VIF_env <- corvif(preds_env)
# all VIF < 2

VIF_inf <- corvif(preds_inf)
# all VIF < 2

# transform coordinates from latlong to utm to avoid issues with correlation structures
coordinates(data_sub) <- data_sub[,c(29,30)]
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84

proj4string(data_sub) <- crs.geo                             # assign the coordinate system

coords_proj <- spTransform(data_sub, CRS("+proj=utm +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "))

data_sub <- cbind(data_sub@data, coords_proj@coords)

names(data_sub)[c(51,52)] <- c("x","y")


# fit the full model to the logged count data with scaled predictors
full.lm <- lm(log10(Count_WW + 1) ~ scale(Area_WHS, scale=T) + scale(Area_SPA, scale=T) + 
              scale(Area_SSSI, scale=T) + scale(Area_RAMSA, scale=T) + 
              scale(Area_NR, scale=T) + scale(Area_NNR, scale=T) + scale(Dist_MPA, scale=T) +
              scale(Dist_MSAC, scale=T) + scale(Dist_MCA, scale=T) + scale(Area_LNR, scale=T) +
              scale(Area_COUNE, scale=T) + scale(Area_CNTRY, scale=T) + 
              scale(Area_BIOSP, scale=T) + scale(Area_BIOGE, scale=T) + 
              scale(Area_NP, scale=T) + scale(Dist_Air, scale=T) + scale(Count_Bus, scale=T) +
              scale(Count_Hotel, scale=T) + scale(Dist_CarPark, scale=T) + 
              scale(Dist_TourOp, scale=T) + scale(Dist_Train, scale=T)  + 
              scale(Dist_road, scale=T) + scale(Mean_Nat, scale=T), data=data_sub)


# look at the summary
summary(full.lm)

# and check the residuals
par(mfrow=c(2,2))
plot(full.lm)

# extract the standardised residuals
S.res.lm <- rstandard(full.lm)

# create a spatial dataframe
mydata_sp <- data_sub
coordinates(mydata_sp) <- c("x","y")

# calculate and plot variograms
Vario <- variogram(S.res.lm ~ 1, data=mydata_sp)
Variodir <- variogram(S.res.lm ~ 1, data=mydata_sp, alpha=c(0,45,90,135))

plot(Vario)
plot(Variodir)

# make a bubble plot of the residuals to check for spatial patterns
bubble.data <- data.frame(S.res.lm, data_sub$x, data_sub$y)
coordinates(bubble.data) <- c("data_sub.x", "data_sub.y")

bubble(bubble.data, "S.res.lm", col=c("black","grey"), 
       main="Residuals", xlab="Longitude", ylab="Latitude")


# obvious spatial autocorrelation in residuals

# GLS ######

# create a new dataframe with scaled predictors and logged response
gls.data<-data.frame(Count_WW=log10(data_sub$Count_WW + 1), Area_WHS=scale(data_sub$Area_WHS, scale=T),
                     Area_SPA=scale(data_sub$Area_SPA, scale=T), Dist_MSAC=scale(data_sub$Dist_MSAC, scale=T),
                     Area_SSSI=scale(data_sub$Area_SSSI, scale=T), Area_RAMSA=scale(data_sub$Area_RAMSA, scale=T),
                     Area_NR=scale(data_sub$Area_NR, scale=T), Area_NNR=scale(data_sub$Area_NNR, scale=T), 
                     Dist_MPA=scale(data_sub$Dist_MPA, scale=T), Dist_MCA=scale(data_sub$Dist_MCA, scale=T), 
                     Area_LNR=scale(data_sub$Area_LNR, scale=T), Area_COUNE=scale(data_sub$Area_COUNE, scale=T),
                     Area_CNTRY=scale(data_sub$Area_CNTRY, scale=T), Area_BIOSP=scale(data_sub$Area_BIOSP, scale=T), 
                     Area_BIOGE=scale(data_sub$Area_BIOGE, scale=T), Area_NP=scale(data_sub$Area_NP, scale=T), 
                     Dist_Air=scale(data_sub$Dist_Air, scale=T), Count_Bus=scale(data_sub$Count_Bus, scale=T),
                     Count_Hotel=scale(log10(data_sub$Count_Hotel + 1), scale = T), 
                     Dist_CarPark=scale(data_sub$Dist_CarPark, scale=T), 
                     Dist_TourOp=scale(data_sub$Dist_TourOp, scale=T), Dist_Train=scale(data_sub$Dist_Train, scale=T), 
                     Dist_road=scale(data_sub$Dist_road, scale=T), Mean_Nat=scale(data_sub$Mean_Nat, scale=T),
                     Area_PA=scale(data_sub$Area_PA, scale = T), Count_Inf=scale(data_sub$Count_Inf, scale = T),
                     Area_SAC_L=scale(data_sub$Area_SAC_L, scale = T), x=data_sub$x, y=data_sub$y)


# Aggregated variables #######
preds.agg<-data_sub[,c("Area_PA","Mean_Nat","Count_Inf")]

# calculate VIF
VIF<-corvif(preds.agg)

agg.gls.ratio <- gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf, data=gls.data,
                     correlation=corRatio(form=~x+y, nugget=T))

agg.gls.sph <- gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf, data=gls.data,
                   correlation=corSpher(form=~x+y, nugget=T))

agg.gls.Lin <- gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf, data=gls.data,
                   correlation=corLin(form=~x+y, nugget=T))

agg.gls.Gaus <- gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf, data=gls.data,
                    correlation=corGaus(form=~x+y, nugget=T))

agg.gls.exp <- gls(Count_WW ~ Area_PA + Mean_Nat + Count_Inf, data=gls.data,
                   correlation=corExp(form=~x+y, nugget=T))


AIC(agg.gls.ratio, agg.gls.sph, agg.gls.Lin, agg.gls.Gaus, agg.gls.exp)

# Exponential autocrrelation structure is the best one

summary(agg.gls.exp)


# check for patterns in the residuals

res.agg.gls <- residuals(agg.gls.exp, type="normalized")
fit.agg.gls <- fitted(agg.gls.exp)

plot(res.agg.gls ~ fit.agg.gls)

qqnorm(res.agg.gls)
qqline(res.agg.gls)

# check that autocorrelation is not an issue anymore
Vario.agg.gls <- Variogram(agg.gls.exp, form= ~x+y, robust=T, resType="normalized")

plot(Vario.agg.gls, smooth=F)

# plot predictions

# Area_PA
newdata <- as.data.frame(lapply(lapply(gls.data, mean), rep, 1131))
newdata$Area_PA_org <- gls.data$Area_PA 
newdata$Count_WW <- gls.data$Count_WW

newdata$Area_PA <- nseq(gls.data$Area_PA, nrow(newdata))

preds_PA <- predict(agg.gls.exp, se.fit=TRUE, newdata)

pred_PA <- ggplot(newdata, aes(x = Area_PA, y = preds_PA$fit)) +
            geom_ribbon(aes(ymin = preds_PA$fit - 1.96 * preds_PA$se.fit, 
                            ymax = preds_PA$fit + 1.96 * preds_PA$se.fit), 
                        alpha = 0.4, fill = "#afcbff") +
            geom_line(size = 2, col = "#afcbff") +
            geom_rug(aes(x = Area_PA_org, y = Count_WW), alpha = 1/2, position = "jitter", size = 0.2)+
            labs(x = "Protected Areas", y = "Wildlife pictures (log10)") +
            theme_bw() +
            theme(plot.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  axis.line = element_line(color = 'black'),
                  axis.text = element_text(size = 16, face = "bold.italic"),
                  axis.title = element_text(size = 18, face = "bold.italic"))


# Mean_Nat

newdata <- as.data.frame(lapply(lapply(gls.data, mean), rep, 1131))
newdata$Mean_Nat_org <- gls.data$Mean_Nat
newdata$Count_WW <- gls.data$Count_WW

newdata$Mean_Nat <- nseq(gls.data$Mean_Nat, nrow(newdata))

preds_nat <- predict(agg.gls.exp, se.fit=TRUE, full=T, newdata)

pred_nat <- ggplot(newdata, aes(x = Mean_Nat, y = preds_nat$fit)) +
             geom_ribbon(aes(ymin = preds_nat$fit - 1.96 * preds_nat$se.fit, 
                             ymax =  preds_nat$fit + 1.96 * preds_nat$se.fit), 
                         alpha = .4, fill = "#afcbff") +
             geom_line(size = 2, col = "#afcbff") +
             geom_rug(aes(x = Mean_Nat_org, y = Count_WW), alpha = 1/2, position = "jitter", size = 0.2)+
             labs(x = "Naturalness", y = "")+
             theme_bw() +
             theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   axis.line = element_line(color = 'black'),
                   axis.text = element_text(size = 16, face = "bold.italic"),
                   axis.title = element_text(size = 18, face = "bold.italic"))


# Count_Inf

newdata <- as.data.frame(lapply(lapply(gls.data, mean), rep, 1131))
newdata$Count_Inf_org <- gls.data$Count_Inf
newdata$Count_WW <- gls.data$Count_WW

newdata$Count_Inf <- nseq(gls.data$Count_Inf, nrow(newdata))

preds_inf <- predict(agg.gls.exp, se.fit=TRUE, full=T, newdata)

pred_inf <- ggplot(newdata, aes(x = Count_Inf, y = preds_inf$fit)) +
             geom_ribbon(aes(ymin = preds_inf$fit - 1.96 * preds_inf$se.fit, 
                             ymax =  preds_inf$fit + 1.96 * preds_inf$se.fit), 
                         alpha = .4, fill = "#afcbff") +
             geom_line(size = 2, col = "#afcbff") +
             geom_rug(aes(x=Count_Inf_org, y=Count_WW), alpha = 1/2, position = "jitter", size = 0.2)+
             labs(x = "Infrastructures", y = "") +
             theme_bw() +
             theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   axis.line = element_line(color = 'black'),
                   axis.text = element_text(size = 16, face = "bold.italic"),
                   axis.title = element_text(size = 18, face = "bold.italic"))


png("Agg_Effects.png", bg = "transparent", width = 17, height = 15, units = "cm", res = 600)
multiplot(pred_PA, pred_nat, pred_inf, cols = 3)
dev.off()

# spatial plot of residuals
res.agg <- gls.data$Count_WW - predict(agg.gls.exp)

res.agg.sp <- data.frame(Residuals = res.agg, Longitude = gls.data$y, Latitude = gls.data$x)

res.agg.plot <- ggplot(data = res.agg.sp, aes(x = Latitude, y = Longitude, height=10000, width=10000))+ 
                 geom_tile(aes(fill = Residuals)) + 
                 scale_fill_gradient2(low = "darkolivegreen", mid = "gold", high = "orangered") +
                 coord_fixed() + scale_x_continuous(labels = function(x){format(x, scientific=F)})+
                 theme_bw() +
                 theme(plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       axis.line = element_line(color = 'black'),
                       axis.text = element_text(size = 6, face = "bold.italic"),
                       axis.title = element_text(size = 9, face = "bold.italic"),
                       legend.title = element_text(size = 6, face = "bold.italic"),
                       legend.text = element_text(size = 9))


png("Agg.res.png", bg = "transparent", width = 8, height = 12, units = "cm", res = 600)
res.agg.plot
dev.off()


# Environmental model #######

# Select best correlation structure

env.gls <- gls(Count_WW ~ Area_WHS + Area_SAC_L + Area_SSSI + Dist_MSAC + Dist_MPA + Dist_MCA +
               Area_SPA + Area_RAMSA + Area_NR + Area_NNR + Area_LNR + Area_COUNE + 
               Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat, data=gls.data)

env.gls.exp <- gls(Count_WW ~ Area_WHS + Area_SAC_L + Area_SSSI + Dist_MSAC + Dist_MPA + Dist_MCA +
                   Area_SPA + Area_RAMSA + Area_NR + Area_NNR + Area_LNR + Area_COUNE + 
                   Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat, data=gls.data, 
                   correlation=corExp(form=~x+y,nugget=T))

env.gls.ratio <- gls(Count_WW ~ Area_WHS + Area_SAC_L + Area_SSSI + Dist_MSAC + Dist_MPA + Dist_MCA +
                     Area_SPA + Area_RAMSA + Area_NR + Area_NNR + Area_LNR + Area_COUNE + 
                     Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat, data=gls.data,
                     correlation=corRatio(form=~x+y,nugget=T))

env.gls.sph <- gls(Count_WW ~ Area_WHS + Area_SAC_L + Area_SSSI + Dist_MSAC + Dist_MPA + Dist_MCA +
                   Area_SPA + Area_RAMSA + Area_NR + Area_NNR + Area_LNR + Area_COUNE + 
                   Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat, data=gls.data,
                   correlation=corSpher(form=~x+y,nugget=T))

env.gls.Lin <- gls(Count_WW ~ Area_WHS + Area_SAC_L + Area_SSSI + Dist_MSAC + Dist_MPA + Dist_MCA +
                   Area_SPA + Area_RAMSA + Area_NR + Area_NNR + Area_LNR + Area_COUNE + 
                   Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat, data=gls.data,
                   correlation=corLin(form=~x+y,nugget=T))

env.gls.Gaus <- gls(Count_WW ~ Area_WHS + Area_SAC_L + Area_SSSI + Dist_MSAC + Dist_MPA + Dist_MCA +
                   Area_SPA + Area_RAMSA + Area_NR + Area_NNR + Area_LNR + Area_COUNE + 
                   Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat, data=gls.data,
                    correlation=corGaus(form=~x+y,nugget=T))

AIC(env.gls,env.gls.exp, env.gls.ratio, env.gls.sph, env.gls.Lin, env.gls.Gaus)


#check for patterns in the residuals

res.env.gls < -residuals(env.gls.exp, type="normalized")
fit.env.gls <- fitted(env.gls.exp)

plot(res.env.gls ~ fit.env.gls)

qqnorm(res.env.gls)
qqline(res.env.gls)

#check that autocorrelation is not an issue anymore
Vario.env.gls <- Variogram(env.gls.exp, form= ~x+y, robust=T, resType = "normalized")

plot(Vario.env.gls, smooth=F)

summary(env.gls.exp)

# Infrastructure model ######

#check for collinearity between linear predictors
preds.inf<-data_sub[,c("Dist_Air","Count_Bus","Count_Hotel","Dist_CarPark","Dist_TourOp",
                   "Dist_Train","Dist_road")]

#calculate VIF
VIF<-corvif(preds.inf)

inf.gls <- gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
               Dist_Train  + Dist_road, data=gls.data)

inf.gls.exp <- gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                   Dist_Train  + Dist_road, data=gls.data, correlation=corExp(form=~x+y, nugget=T))

inf.gls.sph <- gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                   Dist_Train  + Dist_road, data=gls.data, correlation=corSpher(form=~x+y, nugget=T))
 
inf.gls.ratio <- gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                     Dist_Train  + Dist_road, data=gls.data, correlation=corRatio(form=~x+y, nugget=T))

inf.gls.Lin <- gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                 Dist_Train  + Dist_road, data=gls.data, correlation=corLin(form=~x+y, nugget=T))

inf.gls.Gaus <- gls(Count_WW ~ Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                  Dist_Train  + Dist_road, data=gls.data, correlation=corGaus(form=~x+y, nugget=T))

AIC(inf.gls,inf.gls.exp, inf.gls.sph, inf.gls.ratio, inf.gls.Lin, inf.gls.Gaus)


#check for patterns in the residuals

res.inf.gls <- residuals(inf.gls.exp, type="normalized")
fit.inf.gls <- fitted(inf.gls.exp)

plot(res.inf.gls ~ fit.inf.gls)

qqnorm(res.inf.gls)
qqline(res.inf.gls)

#check that autocorrelation is not an issue anymore
Vario.inf.gls <- Variogram(inf.gls.exp, form= ~x+y, robust=T, resType = "normalized")

plot(Vario.inf.gls, smooth=F)

summary(inf.gls.exp)

# Variable selection ######

# use pdredge to use parallell computing
# calculate the number of cores
cores <- detectCores()

# Set up a cluster with number of cores specified as result of detectCores() 
# and call it "clust" 
clust <- makeCluster(cores)

# Load required packages onto worker nodes
# (in this example, load packages {nlme} and {MuMIn} to be used by pdredge)
clusterEvalQ(clust,library(nlme))
clusterEvalQ(clust,library(MuMIn))
clusterExport(clust,"gls.data")     #export the dataframe to the cluster
clusterExport(clust,"inf.gls.exp")  #export the model to the cluster

# model selection for infrastructure model
inf.sel<-pdredge(inf.gls.exp, cluster=clust, rank="AICc", trace=2, REML=FALSE)

saveRDS(inf.sel,"Infrastructure_sel.rds")

stopCluster(clust)

# model selection for environmental model
# this code was used to initiate a cluster and run the model selection in parallel
# on the University of Aberdeen HPC Maxwell.
# For bash script see Env_ModSel_bash.txt in this repository.
# Do not run

# library(nlme)
# library(MuMIn)
# library(snow)
# library(snowfall)
# library(Rmpi)
# 
# gls.data<-read.table("gls.data.txt",stringsAsFactors = F,header=T)
# 
# env.gls.exp<-gls(Count_WW~ Area_WHS + Area_SAC_L + Area_SSSI + Dist_MSAC +Dist_MPA +Dist_MCA +
#                    Area_SPA +Area_RAMSA+ Area_NR + Area_NNR+   Area_LNR + Area_COUNE + 
#                    Area_CNTRY + Area_BIOSP + Area_BIOGE + Area_NP + Mean_Nat,
#                  data=gls.data, correlation=corExp(form=~x+y,nugget=T))
# 
# clust <- makeCluster(mpi.universe.size())
# 
# # Load required packages onto worker nodes
# #   (in this example, load packages {nlme} and {MuMIn} to be used by pdredge)
# clusterEvalQ(clust,library(nlme))
# clusterEvalQ(clust,library(MuMIn))
# clusterExport(clust,"gls.data")  #export the dataframe to the cluster
# clusterExport(clust,"env.gls.exp")  #export the model to the cluster
# 
# env.sel<-pdredge(env.gls.exp, cluster=clust, rank = "AICc",trace=2, REML = FALSE,
#                  subset = !("Area_SSSI"&& "Area_SAC_L")) #do not include collinear variables in the same model
# 
# saveRDS(env.sel,"Environment_sel.rds")
# 
# stopCluster(clust)


# Variable importance and model averaging #####

# Infrastructures

#inf.sel <- readRDS("Infrastructure_sel.rds")

# refit subset of models with REML
inf.sel.REML <- get.models(inf.sel, subset=delta < 5)

#saveRDS(inf.sel.REML, "Infrastructure_Best.rds")

#inf.sel.REML <- readRDS("Infrastructure_Best.rds")

# calculate variable importance 
inf.var.imp <- importance(inf.sel.REML)

# put it into a dataframe format
df <- as.data.frame(inf.var.imp)
inf.var.imp <- cbind(df, attr(inf.var.imp, "names"))
names(inf.var.imp) <- c("Importance", "Var")

# and plot

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
                       #axis.text.y = element_blank(),
                       panel.grid = element_blank(),
                       #axis.text.x = element_blank(),
                       legend.title = element_text(size = 15),
                       legend.text = element_text(size = 12))

png("Infr_VarImp.png", bg = "transparent", width = 17, height = 17, units = "cm", res = 600)
InfVarImpVis +coord_polar()
dev.off()

# model averaging
inf.avg <- model.avg(inf.sel.REML, revised.var=TRUE) 
summary(inf.avg) 

# plotting predictions

# Count_Hotel
newdata <- as.data.frame(lapply(lapply(gls.data, mean), rep, 1131))
newdata$Count_Hotel_org <- gls.data$Count_Hotel 
newdata$Count_WW <- gls.data$Count_WW

newdata$Count_Hotel <- nseq(gls.data$Count_Hotel, nrow(newdata))

preds_hotel <- predict(inf.avg, se.fit=TRUE, full=T, newdata)
# full=T defines the use of  full model-averaged coefficients

pred_hotel <- ggplot(newdata, aes(x = Count_Hotel, y = preds_hotel$fit)) +
               geom_ribbon(aes(ymin = preds_hotel$fit - 1.96 * preds_hotel$se.fit, 
                               ymax =  preds_hotel$fit+1.96*preds_hotel$se.fit), 
                           alpha = .4, fill = "thistle3") +
               geom_line(size = 2, col = "thistle4") +
               geom_rug(aes(x = Count_Hotel_org, y = Count_WW), alpha = 1/2, position = "jitter", size = 0.2) +
               labs(x = "Hotels", y = "Wildlife pictures (log10)") +
               theme_bw() +
               theme(plot.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(color = 'black'),
                     axis.text = element_text(size = 16, face = "bold.italic"),
                     axis.title = element_text(size = 18, face = "bold.italic"))


# Count_Bus
newdata <- as.data.frame(lapply(lapply(gls.data, mean), rep, 1131))
newdata$Count_Bus_org <- gls.data$Count_Bus 
newdata$Count_WW <- gls.data$Count_WW

newdata$Count_Bus <- nseq(gls.data$Count_Bus, nrow(newdata))

preds_bus <- predict(inf.avg, se.fit=TRUE, full=T, newdata)
# full=T defines the use of  full model-averaged coefficients

pred_bus <- ggplot(newdata, aes(x = Count_Bus, y = preds_bus$fit)) +
             geom_ribbon(aes(ymin = preds_bus$fit - 1.96 * preds_bus$se.fit, 
                             ymax = preds_bus$fit + 1.96 * preds_bus$se.fit), 
                         alpha = .4, fill = "thistle3") +
             geom_line(size = 2, col = "Thistle4") +
             geom_rug(aes(x=Count_Bus_org, y=Count_WW), alpha = 1/2, position = "jitter", size = 0.2)+
             labs(x = "Bus Stations", y = "")+
             theme_bw() +
             theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   axis.line = element_line(color = 'black'),
                   axis.text = element_text(size = 16, face = "bold.italic"),
                   axis.title = element_text(size = 18, face = "bold.italic"))

# Dist_Air
newdata <- as.data.frame(lapply(lapply(gls.data, mean), rep, 1131))
newdata$Dist_Air_org <- gls.data$Dist_Air 
newdata$Count_WW <- gls.data$Count_WW

newdata$Dist_Air <- nseq(gls.data$Dist_Air, nrow(newdata))

preds_air <- predict(inf.avg, se.fit=TRUE, full=T, newdata)
# full=T defines the use of  full model-averaged coefficients 

pred_air <- ggplot(newdata, aes(x = Dist_Air, y = preds_air$fit)) +
             geom_ribbon(aes(ymin = preds_air$fit - 1.96 * preds_air$se.fit, 
                             ymax = preds_air$fit + 1.96 * preds_air$se.fit), 
                         alpha = .5, fill = "thistle3") +
             geom_line(size = 2, col = "Thistle4") +
             geom_rug(aes(x=Dist_Air_org, y=Count_WW), alpha = 1/2, position = "jitter", size = 0.2)+
             labs(x = "Airports (dist)", y = "") +
             theme_bw() +
             theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   axis.line = element_line(color = 'black'),
                   axis.text = element_text(size = 16, face = "bold.italic"),
                   axis.title = element_text(size = 18, face = "bold.italic"))

png("Inf_Effects.png", bg = "transparent", width = 17, height = 15, units = "cm", res = 600)
multiplot(pred_hotel, pred_bus, pred_air, cols = 3)
dev.off()

# spatial plots of residuals 
res.inf <- gls.data$Count_WW - predict(inf.avg, full = T)

res.inf.sp <- data.frame(Residuals = res.inf, Longitude = gls.data$y, Latitude = gls.data$x)


res.plot <- ggplot(data = res.inf.sp, aes(x = Latitude, y = Longitude, height=10000, width=10000))+ 
             geom_tile(aes(fill = Residuals)) + 
             scale_fill_gradient2(low = "darkolivegreen", mid = "gold", high = "orangered") +
             coord_fixed() + 
             scale_x_continuous(labels = function(x){format(x, scientific=F)})+
             theme_bw() +
             theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   axis.line = element_line(color = 'black'),
                   axis.text = element_text(size = 6, face = "bold.italic"),
                   axis.title = element_text(size = 9, face = "bold.italic"),
                   legend.title = element_text(size = 6, face = "bold.italic"),
                   legend.text = element_text(size = 9))


png("Inf.res.png", bg = "transparent", width = 8, height = 12, units = "cm", res = 600)
res.plot
dev.off()

# Environment 
# this was run on on the UoA HPC Maxwell in an interactive session
# subset object to avoid running out of memory

env.sel <- readRDS("Environment_sel.rds")
env.sel.sub1 <- subset(env.sel, delta < 1)
saveRDS(env.sel.sub1, "Environment_sel1.rds")
rm(env.sel.sub1)

env.sel.sub2 <- subset(env.sel, delta>=1 & delta < 2)
saveRDS(env.sel.sub2, "Environment_sel2.rds")
rm(env.sel.sub2)

env.sel.sub3 <- subset(env.sel, delta>=2 & delta < 3)
saveRDS(env.sel.sub3, "Environment_sel3.rds")
rm(env.sel.sub3)

env.sel.sub4 <- subset(env.sel, delta>=3 & delta < 3.5)
saveRDS(env.sel.sub4, "Environment_sel4.rds")
rm(env.sel.sub4)

env.sel.sub5 <- subset(env.sel, delta> 3.5 & delta < 4)
saveRDS(env.sel.sub5, "Environment_sel5.rds")
rm(env.sel.sub5)

env.sel.sub6 <- subset(env.sel, delta> 4 & delta < 4.5)
saveRDS(env.sel.sub6, "Environment_sel6.rds")
rm(env.sel.sub6)

env.sel.sub7 <- subset(env.sel, delta> 4.5 & delta < 5)
saveRDS(env.sel.sub7, "Environment_sel7.rds")
rm(env.sel.sub7)

# now refit each subset of models with REML

env.sel.sub1 <- readRDS("Environment_sel1.rds")
str(env.sel.sub1)
env.sel.REML1 <- get.models(env.sel.sub1, subset = TRUE, method = "REML")
saveRDS(env.sel.REML1,"Environment_Best1.rds")


env.sel.sub2 <- readRDS("Environment_sel2.rds")
str(env.sel.sub2)
env.sel.REML2 <- get.models(env.sel.sub2, subset = TRUE, method = "REML")
saveRDS(env.sel.REML2,"Environment_Best2.rds")


env.sel.sub3 <- readRDS("Environment_sel3.rds")
str(env.sel.sub3)
env.sel.REML3 <- get.models(env.sel.sub3, subset = TRUE, method = "REML")
saveRDS(env.sel.REML3,"Environment_Best3.rds")


env.sel.sub4 <- readRDS("Environment_sel4.rds")
str(env.sel.sub4)
env.sel.REML4 <- get.models(env.sel.sub4, subset = TRUE, method = "REML")
saveRDS(env.sel.REML4,"Environment_Best4.rds")


env.sel.sub5 <- readRDS("Environment_sel5.rds")
str(env.sel.sub5)
env.sel.REML5 <- get.models(env.sel.sub5, subset = TRUE, method = "REML")
saveRDS(env.sel.REML5,"Environment_Best5.rds")



env.sel.sub6 <- readRDS("Environment_sel6.rds")
str(env.sel.sub6)
env.sel.REML6 <- get.models(env.sel.sub6, subset = TRUE, method = "REML")
saveRDS(env.sel.REML6,"Environment_Best6.rds")


env.sel.sub7 <- readRDS("Environment_sel7.rds")
str(env.sel.sub7)
env.sel.REML7 <- get.models(env.sel.sub7, subset = TRUE, method = "REML")
saveRDS(env.sel.REML7,"Environment_Best7.rds")



env.sel.REML1<-readRDS("Environment_Best1.rds")

env.sel.REML2<-readRDS("Environment_Best2.rds")

env.sel.REML3<-readRDS("Environment_Best3.rds")

env.sel.REML4<-readRDS("Environment_Best4.rds")

env.sel.REML5<-readRDS("Environment_Best5.rds")

env.sel.REML6<-readRDS("Environment_Best6.rds")

env.sel.REML7.1<-readRDS("Environment_Best7.1.rds")

env.sel.REML7.2<-readRDS("Environment_Best7.2.rds")


# calculate variable importance 
env.sel.REML <- c(env.sel.REML1, env.sel.REML2, env.sel.REML3, 
                  env.sel.REML4, env.sel.REML5, env.sel.REML6,
                  env.sel.REML7.1, env.sel.REML7.2)

env.var.imp <- importance(env.sel.REML)

# put it into a dataframe format
df <- as.data.frame(env.var.imp)
env.var.imp <- cbind(df, attr(env.var.imp, "names"))
names(env.var.imp) <- c("Importance", "Var")

# and plot

library(ggplot2)

EnvVarImpVis <- ggplot (env.var.imp, aes(Var, Importance, fill = Importance)) +
                 geom_hline(yintercept = seq(0, 1.2, by = 0.5), colour = "grey90", size = 1) +
                 geom_vline(aes(xintercept = Var), colour = "grey90", size = 1) +
                 geom_bar(width = 1, stat = "identity", color = "white") +
                 scale_y_continuous(breaks = 0:nlevels(as.factor(env.var.imp$Var))) +
                 scale_fill_gradient(low = "darkolivegreen1", high = "darkolivegreen4") +
                 theme_bw() +
                 theme(axis.title = element_blank(),
                       panel.border = element_blank(),
                       axis.ticks = element_blank(),
                       axis.text.y = element_blank(),
                       panel.grid = element_blank(),
                       axis.text.x = element_blank(),
                       legend.title = element_text(size = 15),
                       legend.text = element_text(size = 12))

png("Env_VarImp.png", bg = "transparent", width = 17, height = 17, units = "cm", res = 600)
EnvVarImpVis + coord_polar()
dev.off()

# model averaging
env.avg <- model.avg(env.sel.REML, revised.var = TRUE) 
summary(env.avg) 

# plotting predictions

# Area_CNTRY
newdata <- as.data.frame(lapply(lapply(gls.data, mean), rep, 1131))
newdata$Area_CNTRY_org <- gls.data$Area_CNTRY 
newdata$Count_WW <- gls.data$Count_WW

newdata$Area_CNTRY <- nseq(gls.data$Area_CNTRY, nrow(newdata))

preds_CNTRY <- predict(env.avg, se.fit=TRUE, full=T, newdata)

pred_CNTRY <- ggplot(newdata, aes(x = Area_CNTRY, y = preds_CNTRY$fit)) +
               geom_ribbon(aes(ymin = preds_CNTRY$fit - 1.96 * preds_CNTRY$se.fit, 
                               ymax = preds_CNTRY$fit + 1.96 * preds_CNTRY$se.fit), 
                           alpha = .4, fill = "darkolivegreen3") +
               geom_line(size = 2, col = "darkolivegreen4") +
               geom_rug(aes(x=Area_CNTRY_org, y=Count_WW), alpha = 1/2, position = "jitter", size = 0.2)+
               labs(x = "Country Parks", y = "Wildlife pictures (log10)") +
               theme_bw() +
               theme(plot.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(color = 'black'),
                     axis.text = element_text(size = 16, face = "bold.italic"),
                     axis.title = element_text(size = 18, face = "bold.italic"))

# Area_LNR
newdata <- as.data.frame(lapply(lapply(gls.data, mean), rep, 1131))
newdata$Area_LNR_org <- gls.data$Area_LNR 
newdata$Count_WW <- gls.data$Count_WW

newdata$Area_LNR <- nseq(gls.data$Area_LNR, nrow(newdata))

preds_LNR <- predict(env.avg, se.fit=TRUE, full=T, newdata)

pred_LNR <- ggplot(newdata, aes(x = Area_LNR, y = preds_LNR$fit)) +
             geom_ribbon(aes(ymin = preds_LNR$fit - 1.96 * preds_LNR$se.fit, 
                             ymax = preds_LNR$fit + 1.96 * preds_LNR$se.fit), 
                         alpha = .4, fill = "darkolivegreen3") +
             geom_line(size = 2, col = "darkolivegreen4") +
             geom_rug(aes(x=Area_LNR_org, y=Count_WW), alpha = 1/2, position = "jitter", size = 0.2)+
             labs(x = "Local Nature Reserves", y = "Wildlife pictures (log10)")+
             theme_bw() +
             theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   axis.line = element_line(color = 'black'),
                   axis.text = element_text(size = 16, face = "bold.italic"),
                   axis.title = element_text(size = 18, face = "bold.italic"))

# Area_NP
newdata <- as.data.frame(lapply(lapply(gls.data, mean), rep, 1131))
newdata$Area_NP_org <- gls.data$Area_NP 
newdata$Count_WW <- gls.data$Count_WW

newdata$Area_NP <- nseq(gls.data$Area_NP, nrow(newdata))

preds_NP <- predict(env.avg, se.fit=TRUE, full=T, newdata)

pred_NP <- ggplot(newdata, aes(x = Area_NP, y = preds_NP$fit)) +
            geom_ribbon(aes(ymin = preds_NP$fit - 1.96 * preds_NP$se.fit, 
                            ymax = preds_NP$fit + 1.96 * preds_NP$se.fit), 
                        alpha = .4, fill = "darkolivegreen3") +
            geom_line(size = 2, col = "darkolivegreen4") +
            geom_rug(aes(x=Area_NP_org, y=Count_WW), alpha = 1/2, position = "jitter", size = 0.2)+
            labs(x = "National Parks", y = "") +
            theme_bw() +
            theme(plot.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  axis.line = element_line(color = 'black'),
                  axis.text = element_text(size = 16, face = "bold.italic"),
                  axis.title = element_text(size = 18, face = "bold.italic"))

# Dist_MSAC
newdata <- as.data.frame(lapply(lapply(gls.data, mean), rep, 1131))
newdata$Dist_MSAC_org <- gls.data$Dist_MSAC 
newdata$Count_WW <- gls.data$Count_WW

newdata$Dist_MSAC <- nseq(gls.data$Dist_MSAC, nrow(newdata))

preds_MSAC <- predict(env.avg, se.fit=TRUE, full=T, newdata)

pred_MSAC <- ggplot(newdata, aes(x = Dist_MSAC, y = preds_MSAC$fit)) +
              geom_ribbon(aes(ymin = preds_MSAC$fit - 1.96 * preds_MSAC$se.fit, 
                              ymax = preds_MSAC$fit + 1.96 * preds_MSAC$se.fit), 
                          alpha = .4, fill = "darkolivegreen3") +
              geom_line(size = 2, col = "darkolivegreen4") +
              geom_rug(aes(x=Dist_MSAC_org, y=Count_WW), alpha = 1/2, position = "jitter", size = 0.2)+
              labs(x = "Marine SAC (dist)", y = "") +
              theme_bw() +
              theme(plot.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color = 'black'),
                    axis.text = element_text(size = 16, face = "bold.italic"),
                    axis.title = element_text(size = 18, face = "bold.italic"))

png("Env_Effects.png", bg = "transparent", width = 17, height = 25, units = "cm", res = 600)
multiplot(pred_CNTRY, pred_LNR, pred_NP, pred_MSAC, cols = 2)
dev.off()

# spatial plots of residuals 
res.env <- gls.data$Count_WW - predict(env.avg, full = T)

res.env.sp <- data.frame(Residuals = res.env, Longitude = gls.data$y, Latitude = gls.data$x)


res.plot <- ggplot(data = res.env.sp, aes(x = Latitude, y = Longitude, height=10000, width=10000))+ 
             geom_tile(aes(fill = Residuals)) + 
             scale_fill_gradient2(low = "darkolivegreen", mid = "gold", high = "orangered") +
             coord_fixed() + 
             scale_x_continuous(labels = function(x){format(x, scientific=F)})+
             theme_bw() +
             theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   axis.line = element_line(color = 'black'),
                   axis.text = element_text(size = 6, face = "bold.italic"),
                   axis.title = element_text(size = 9, face = "bold.italic"),
                   legend.title = element_text(size = 6, face = "bold.italic"),
                   legend.text = element_text(size = 9))


png("Env_res.png", bg = "transparent", width = 8, height = 12, units = "cm", res = 600)
res.plot
dev.off()


# Biodiversity #####

# look at distribution of biodiversity records
par(mfrow=c(2,1))
hist(log(data_sub$Records + 1), main="Histogram of species records", xlab="Species Records (log)")
plot(density(log(data_sub$Records + 1)), main="Density of species records")

# look at distribution of biodiversity richness
par(mfrow=c(2,1))
hist(log(data_sub$Species + 1), main="Histogram of species richness", xlab="Species Richness (log)")
plot(density(log(data_sub$Species + 1)), main="Density of species richness")

#look at distribution of species richness where Records are above 7
par(mfrow=c(2,1))
hist(log(data_sub$Species[which(data_sub$Records > 7)]), 
     main="Histogram of species richness", xlab="Species Richness (log)")
plot(density(log(data_sub$Species[which(data_sub$Records > 7)])), main="Density of species richness")


# Effect of biodiversity 

#read classification from mixture model
classes <- read.table(".//data//classification.txt")
names(classes) <- "Class"

# calculate log of records
data_sub$logrec <- log(data_sub$Records + 1)

# only select non 0 records
data.pos <- data_sub[data_sub$logrec > 0,]

# now select only records belonging to group 2
data.gr2 <- cbind(data.pos, classes)
data.gr2 <- data.gr2[data.gr2$Class==2,]
str(data.gr2)

# create a new dataframe with scaled predictors and logged response
gls.data.bio <- data.frame(Count_WW=log10(data.gr2$Count_WW + 1), Dist_MSAC=scale(data.gr2$Dist_MSAC, scale=T),
                     Area_LNR=scale(data.gr2$Area_LNR, scale=T), Area_CNTRY=scale(data.gr2$Area_CNTRY, scale=T), 
                     Area_NP=scale(data.gr2$Area_NP, scale=T), Mean_Nat=scale(data.gr2$Mean_Nat, scale=T),
                     Species=log10(data.gr2$Species), x=data.gr2$x, y=data.gr2$y)


bio.gls.exp <- gls(Count_WW ~ Dist_MSAC + Area_LNR + Area_CNTRY + Area_NP + Mean_Nat + Species,
                   data = gls.data.bio, correlation = corExp(form = ~ x+y, nugget = T))

bio.gls.lin <- gls(Count_WW ~ Dist_MSAC + Area_LNR + Area_CNTRY + Area_NP + Mean_Nat + Species,
                   data = gls.data.bio, correlation = corLin(form = ~ x+y, nugget = T))

bio.gls.ratio <- gls(Count_WW ~ Dist_MSAC + Area_LNR + Area_CNTRY + Area_NP + Mean_Nat + Species,
                     data = gls.data.bio, correlation = corRatio(form = ~ x+y, nugget = T))

bio.gls.gaus <- gls(Count_WW ~ Dist_MSAC + Area_LNR + Area_CNTRY + Area_NP + Mean_Nat + Species,
                    data = gls.data.bio, correlation = corGaus(form = ~ x+y, nugget = T))

bio.gls.sph <- gls(Count_WW ~ Dist_MSAC + Area_LNR + Area_CNTRY + Area_NP + Mean_Nat + Species,
                   data = gls.data.bio, correlation = corSpher(form = ~ x+y, nugget = T))

AIC(bio.gls.exp, bio.gls.lin, bio.gls.ratio, bio.gls.gaus, bio.gls.sph)

summary(bio.gls.exp)

# check for patterns in the residuals
res.bio.gls <- residuals(bio.gls.exp, type="normalized")

fit.bio.gls <- fitted(bio.gls.exp)

par(mfrow=c(1,2))
plot(res.bio.gls ~ fit.bio.gls)

qqnorm(res.bio.gls)
qqline(res.bio.gls)


# check for autocorrelation

mydata_sp <- gls.data.bio
coordinates(mydata_sp)<-c("x", "y")

Vario <- variogram(res.bio.gls ~ 1, data=mydata_sp)

plot(Vario)

# plotting predictions

# Species

newdata <- as.data.frame(lapply(lapply(gls.data.bio, mean), rep, 760))
newdata$Species_org <- gls.data.bio$Species 
newdata$Count_WW <- gls.data.bio$Count_WW

newdata$Species <- nseq(gls.data.bio$Species, nrow(newdata))

preds_sp <- predict(bio.gls.exp, se.fit=TRUE, newdata)

pred_sp <- ggplot(newdata, aes(x = Species, y = preds_sp$fit)) +
            geom_ribbon(aes(ymin = preds_sp$fit - 1.96 * preds_sp$se.fit, 
                            ymax =  preds_sp$fit + 1.96 * preds_sp$se.fit), 
                        alpha = 0.4, fill = "#A32765") +
            geom_line(size = 2, col = "#A32765") +
            geom_rug(aes(x = Species_org, y = Count_WW), alpha = 1/2, position = "jitter", size = 0.2) +
            labs(x = "Species (log10)", y = "Wildlife pictures (log10)") +
            theme_bw() +
            theme(plot.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  axis.line = element_line(color = 'black'),
                  axis.text = element_text(size = 16, face = "bold.italic"),
                  axis.title = element_text(size = 18, face = "bold.italic"))

png("Bio_Effect.png", bg = "transparent", width = 8, height = 10, units = "cm", res = 600)
pred_sp
dev.off()

# spatial plot of residuals
res.bio <- gls.data.bio$Count_WW - predict(bio.gls.exp)

res.bio.sp <- data.frame(Residuals = res.bio, Longitude = gls.data.bio$y, Latitude = gls.data.bio$x)

res.bio.plot <- ggplot(data = res.bio.sp, aes(x = Latitude, y = Longitude, height=10000, width=10000))+ 
                 geom_tile(aes(fill = Residuals)) + 
                 scale_fill_gradient2(low = "darkolivegreen", mid = "gold", high = "orangered") +
                 coord_fixed() + 
                 scale_x_continuous(labels = function(x){format(x, scientific=F)})+
                 theme_bw() +
                 theme(plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       axis.line = element_line(color = 'black'),
                       axis.text = element_text(size = 6, face = "bold.italic"),
                       axis.title = element_text(size = 9, face = "bold.italic"),
                       legend.title = element_text(size = 6, face = "bold.italic"),
                       legend.text = element_text(size = 9))


png("Bio.res.png", bg = "transparent", width = 8, height = 12, units = "cm", res = 600)
res.bio.plot
dev.off()


