####################################################
###########Francesca Mancini
###########created 21/02/2017
####################################################

#load required libraries
library(gstat)
library(MASS)
library(sp)

#load the data from the subfolder //data
data<-read.table(".//data//CombinedData_v9.txt",stringsAsFactors = F,header=T)

#subset the data to exclude all the observation without natural values
data_sub<-data[is.na(data$Mean_Nat)==F,]

#check for collinearity between linear predictors
#first create a dataframe containing only the predictors of interest
pred.df<-data_sub[,c(2:15,18:19,29:31,33,36:37,39,41:42)]

#then calculate variance inflation factor
#the following functions are written by Gudrun Carl, 2005-2007 
#and available at http://www.ecography.org/appendix/e5171
#Dormann, F. C., McPherson, J. M., Araújo, M. B., Bivand, R., Bolliger, J., 
#Carl, G., Davies, R. G., Hirzel, A., Jetz, W., Kissling, W. D., Kühn, I., 
#Ohlemüller, R., Peres-Neto, P. R., Reineking, B., Schröder, B., Schurr, F. M. 
#and Wilson, R. 2007. Methods to account for spatial autocorrelation in the analysis 
#of species distributional data: a review. – Ecography 30: 609–628.

#function to calculate VIF for the predictors
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  #correlation part
  cat("Correlations of the variables\n\n")
  tmp_cor <- cor(dataz,use="complete.obs")
  print(tmp_cor)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

#calculate VIF
VIF<-corvif(pred.df)

#fit the full model to the logged count data with scaled predictors

full.lm<-lm(log10(Count_WW+1)~scale(Area_WHS,scale=T) + scale(Area_SSSI,scale=T) + scale(Area_SPA,scale=T) + 
              scale(Area_SAC,scale=T) +scale(Area_RAMSA,scale=T) + scale(Area_NR,scale=T) +
              scale(Area_NNR,scale=T) + scale(Area_MPAdi,scale=T) +
              scale(Area_MCA,scale=T) + scale(Area_LNR,scale=T) + scale(Area_COUNE,scale=T) +
              scale(Area_CNTRY,scale=T) + scale(Area_BIOSP,scale=T) +
              scale(Area_BIOGE,scale=T) + scale(Area_NP,scale=T) + scale(Dist_Air,scale=T) + scale(Count_Bus,scale=T) +
              scale(Count_Hotel,scale=T) + scale(Dist_CarPark,scale=T) + scale(Dist_TourOp,scale=T) +
              scale(Dist_Train,scale=T)  + scale(Dist_road,scale=T) + scale(Mean_Nat,scale=T) +
              offset(log10(Pop_dens+1)),data=data_sub)

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

#create a new dataframe with scaled predictorsand logged response
gls.data<-data.frame(Count_WW=log10(data_sub$Count_WW+1),Area_WHS=scale(data_sub$Area_WHS,scale=T),
                     Area_SSSI=scale(data_sub$Area_SSSI,scale=T), Area_SPA=scale(data_sub$Area_SPA,scale=T), 
                     Area_SAC=scale(data_sub$Area_SAC,scale=T), Area_RAMSA=scale(data_sub$Area_RAMSA,scale=T),
                     Area_NR=scale(data_sub$Area_NR,scale=T), Area_NNR=scale(data_sub$Area_NNR,scale=T), 
                     Area_MPAdi=scale(data_sub$Area_MPAdi,scale=T), Area_MCA=scale(data_sub$Area_MCA,scale=T), 
                     Area_LNR=scale(data_sub$Area_LNR,scale=T), Area_COUNE=scale(data_sub$Area_COUNE,scale=T),
                     Area_CNTRY=scale(data_sub$Area_CNTRY,scale=T), Area_BIOSP=scale(data_sub$Area_BIOSP,scale=T), 
                     Area_BIOGE=scale(data_sub$Area_BIOGE,scale=T), Area_NP=scale(data_sub$Area_NP,scale=T), 
                     Dist_Air=scale(data_sub$Dist_Air,scale=T), Count_Bus=scale(data_sub$Count_Bus,scale=T),
                     Count_Hotel=scale(data_sub$Count_Hotel,scale=T), Dist_CarPark=scale(data_sub$Dist_CarPark,scale=T), 
                     Dist_TourOp=scale(data_sub$Dist_TourOp,scale=T), Dist_Train=scale(data_sub$Dist_Train,scale=T), 
                     Dist_road=scale(data_sub$Dist_road,scale=T), Mean_Nat=scale(data_sub$Mean_Nat,scale=T),
                     Pop_dens=log10(data_sub$Pop_dens+1),x=data_sub$Longitude,y=data_sub$Latitude)


#fit the full model with gls
full.gls<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +Area_RAMSA+ Area_NR +
                Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                Dist_Train  + Dist_road + Mean_Nat + offset(Pop_dens),data=gls.data)

summary(full.gls)

#check for autocorrelation
Vario.gls<- Variogram(full.gls,form= ~x+y,robust=T,resType = "pearson",maxDist=1)

plot(Vario.gls,smooth=T)

#fit the same model with different autocorrelation structures
full.gls.sph<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +Area_RAMSA+ Area_NR +
                    Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                    Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                    Dist_Train  + Dist_road + Mean_Nat + offset(Pop_dens),data=gls.data,
                  correlation=corSpher(form=~x+y,nugget=T))

full.gls.Lin<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +Area_RAMSA+ Area_NR +
                    Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                    Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                    Dist_Train  + Dist_road + Mean_Nat + offset(Pop_dens),data=gls.data,
                  correlation=corLin(form=~x+y,nugget=T))

full.gls.ratio<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +Area_RAMSA+ Area_NR +
                      Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                      Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                      Dist_Train  + Dist_road + Mean_Nat + offset(Pop_dens),data=gls.data,
                    correlation=corRatio(form=~x+y,nugget=T))

full.gls.Gaus<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +Area_RAMSA+ Area_NR +
                     Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                     Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                     Dist_Train  + Dist_road + Mean_Nat + offset(Pop_dens),data=gls.data,
                   correlation=corGaus(form=~x+y,nugget=T))

full.gls.exp<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +Area_RAMSA+ Area_NR +
                    Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                    Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                    Dist_Train  + Dist_road + Mean_Nat + offset(Pop_dens),data=gls.data,
                  correlation=corExp(form=~x+y,nugget=T))

#use AIC to choose the best one
AIC(full.gls,full.gls.sph,full.gls.Lin,full.gls.ratio,full.gls.Gaus,full.gls.exp)

#model with exponential correlation structure is the best model

summary(full.gls.exp)


#check for patterns in the residuals
plot(full.gls.exp)

res.gls<-residuals(full.gls.exp,type="normalized")
fit.gls<-fitted(full.gls.exp)

plot(res.gls~fit.gls)

qqnorm(res.gls)
qqline(res.gls)

#check that autocorrelation is not an issue anymore
Vario.gls.exp<-Variogram(full.gls.exp,form= ~x+y,robust=T,resType = "normalized")

plot(Vario.gls.exp, smooth=F)

bubble.data.gls<-data.frame(res.gls,data_sub$Longitude,data_sub$Latitude)
coordinates(bubble.data.gls)<-c("data_sub.Longitude","data_sub.Latitude")

bubble(bubble.data.gls,"res.gls",col=c("black","grey"),main="Residuals",xlab="Longitude",ylab="Latitude")

#######################
##model selection
#######################

#fit every alternative model with maximum likelihood

full.gls.exp<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +Area_RAMSA+ Area_NR +
                    Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                    Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                    Dist_Train  + Dist_road + Mean_Nat + offset(Pop_dens),data=gls.data,
                  correlation=corExp(form=~x+y,nugget=T),method="ML")

gls.exp1<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +Area_RAMSA+ Area_NR +
                Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                Dist_Train  + Mean_Nat + offset(Pop_dens),data=gls.data,
              correlation=corExp(form=~x+y,nugget=T),method="ML")


gls.exp2<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +Area_RAMSA+ Area_NR +
                Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                Area_BIOGE + Area_NP + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                Mean_Nat + offset(Pop_dens),data=gls.data,
              correlation=corExp(form=~x+y,nugget=T),method="ML")




gls.exp3<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +Area_RAMSA+ Area_NR +
                Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                Area_BIOGE + Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                Mean_Nat + offset(Pop_dens),data=gls.data,
              correlation=corExp(form=~x+y,nugget=T),method="ML")



gls.exp4<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +Area_RAMSA+ Area_NR +
                Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_COUNE + Area_CNTRY + Area_BIOSP +
                Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                Mean_Nat + offset(Pop_dens),data=gls.data,
              correlation=corExp(form=~x+y,nugget=T),method="ML")



gls.exp5<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +Area_RAMSA+ Area_NR +
                Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_COUNE + Area_CNTRY +
                Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                Mean_Nat + offset(Pop_dens),data=gls.data,
              correlation=corExp(form=~x+y,nugget=T),method="ML")



gls.exp6<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +Area_RAMSA+ Area_NR +
                Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_CNTRY +
                Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                Mean_Nat + offset(Pop_dens),data=gls.data,
              correlation=corExp(form=~x+y,nugget=T),method="ML")



gls.exp7<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +Area_RAMSA+
                Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_CNTRY +
                Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                Mean_Nat + offset(Pop_dens),data=gls.data,
              correlation=corExp(form=~x+y,nugget=T),method="ML")



gls.exp8<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +
                Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_CNTRY +
                Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                Mean_Nat + offset(Pop_dens),data=gls.data,
              correlation=corExp(form=~x+y,nugget=T),method="ML")

gls.exp9<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SAC +
                Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_CNTRY +
                Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                Mean_Nat + offset(Pop_dens),data=gls.data,
              correlation=corExp(form=~x+y,nugget=T),method="ML")


#use AIC to determine best model

AIC(full.gls.exp, gls.exp1, gls.exp2, gls.exp3, gls.exp4, gls.exp5, gls.exp6, gls.exp7, gls.exp8,gls.exp9)



gls.exp8<-gls(Count_WW~ Area_WHS + Area_WHS + Area_SPA + Area_SAC +
                Area_NNR+ Area_MPAdi + Area_MCA + Area_LNR + Area_CNTRY +
                Dist_Air + Count_Bus + Count_Hotel + Dist_CarPark + Dist_TourOp +
                Mean_Nat + offset(Pop_dens),data=gls.data,
              correlation=corExp(form=~x+y,nugget=T),method="REML")

########################
##model averaging!


