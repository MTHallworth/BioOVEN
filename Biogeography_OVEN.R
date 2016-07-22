
.libPaths("D:/Google_Drive/Hallworth_R_Library")

library(raster)
library(sp)
library(rgdal)
library(rgeos)

######################################################################
#### Weighted Correlation Mean Min Temperature and Latitude ####

# Breeding temperature gradient #

OVENbbs<-raster("Spatial_Layers/bbsoven.txt")

clim<-getData('worldclim',var='tmin',res=10)

BreedMin<-(sum(clim[[5:8]])/4)/10

BreedMin<-crop(BreedMin,extent(OVENbbs))

BreedMin<-resample(BreedMin,OVENbbs)

BreedMin<-mask(BreedMin,OVENbbs)

Breedtemps<-rasterToPoints(BreedMin)

BreedAbund<-extract(OVENbbs,cbind(Breedtemps[,1],Breedtemps[,2]))


# Non-breeding temperature gradient #

Carib<-shapefile("Spatial_Layers/OVEN_Caribbean.shp")
CentAm<-shapefile("Spatial_Layers/OVEN_CentralAmerica.shp")

OVENwinter<-gUnion(Carib,CentAm)

OVENebird<-raster("Spatial_Layers/ebirdoven.txt")

OVENebird<-mask(OVENebird,OVENwinter)

WinterMin<-(sum(clim[[c(1:4,11:12)]]/6))/10
WinterMin<-crop(WinterMin,OVENwinter)
WinterMin<-mask(WinterMin,OVENwinter)
WinterMin<-resample(WinterMin,OVENebird)

WinterTemps<-rasterToPoints(WinterMin)

WinterAbund<-extract(OVENebird,cbind(WinterTemps[,1],WinterTemps[,2]))

######################################################################
#### FIGURE  ####
par(bty="l")
plot(Breedtemps[,3]~Breedtemps[,2],
     xlim=c(7,60),ylim=c(0,25),
     pch=19,col="gray",cex=BreedAbund/10,
     ylab="Mean minimum temperature",xlab="Latitude",
     yaxt="n")
axis(2,las=2)
par(new=TRUE)
plot(WinterTemps[,3]~WinterTemps[,2],
     xlim=c(7,60),ylim=c(0,25),
     pch=19,col="black",cex=WinterAbund/2,
     ylab="",xlab="",axes=F)
legend(12,26, legend="Non-breeding",bty="n")
legend(43,26, legend="Breeding",bty="n")
abline(v=31)
abline(lm(c(WinterTemps[,3],Breedtemps[,3])~c(WinterTemps[,2],Breedtemps[,2])))
abline(lm(WinterTemps[,3]~WinterTemps[,2],weight=WinterAbund),col="black")
abline(lm(Breedtemps[,3]~Breedtemps[,2],weight=BreedAbund),col="darkgray")

#####

precip<-getData('worldclim',var='prec',res=10)


BreedPrecip<-(sum(precip[[5:8]])/4)

BreedPrecip<-crop(BreedPrecip,extent(OVENbbs))

BreedPrecip<-resample(BreedPrecip,OVENbbs)

BreedPrecip<-mask(BreedPrecip,OVENbbs)

BreedPrecipAmt<-rasterToPoints(BreedPrecip)

BreedAbund<-extract(OVENbbs,cbind(BreedPrecipAmt[,1],BreedPrecipAmt[,2]))


# Non-breeding temperature gradient #

Carib<-shapefile("Spatial_Layers/OVEN_Caribbean.shp")
CentAm<-shapefile("Spatial_Layers/OVEN_CentralAmerica.shp")

OVENwinter<-gUnion(Carib,CentAm)

OVENebird<-raster("Spatial_Layers/ebirdoven.txt")

OVENebird<-mask(OVENebird,OVENwinter)

WinterPrecip<-(sum(precip[[c(1:4,11:12)]]/6))
WinterPrecip<-crop(WinterPrecip,OVENwinter)
WinterPrecip<-mask(WinterPrecip,OVENwinter)
WinterPrecip<-resample(WinterPrecip,OVENebird)

WinterPrecipAmt<-rasterToPoints(WinterPrecip)

WinterAbund<-extract(OVENebird,cbind(WinterPrecipAmt[,1],WinterPrecipAmt[,2]))

#### FIGURE Latitude  ####
par(bty="l")
plot(BreedPrecipAmt[,3]~BreedPrecipAmt[,2],
     xlim=c(7,60),ylim=c(0,410),
     pch=19,col="gray",cex=BreedAbund/10,
     ylab="Mean precipitation",xlab="Latitude",
     yaxt="n")
axis(2,las=2)
par(new=TRUE)
plot(WinterPrecipAmt[,3]~WinterPrecipAmt[,2],
     xlim=c(7,60),ylim=c(0,410),
     pch=19,col="black",cex=WinterAbund/2,
     ylab="",xlab="",axes=F)
legend(12,410, legend="Non-breeding",bty="n")
legend(43,410, legend="Breeding",bty="n")
abline(v=31)

#### FIGURE Long Winter  ####
par(bty="l",mfrow=c(1,2),mar=c(4,4,4,1))
plot(WinterPrecipAmt[,3]~WinterPrecipAmt[,1],
    ylim=c(0,410),
    ylab="Mean Precipitation",xlab="Longitude",
     pch=19,col="black",cex=WinterAbund/2)
legend(-110,410, legend="Non-breeding",bty="n")
par(mar=c(4,0,4,4))
plot(BreedPrecipAmt[,3]~BreedPrecipAmt[,1],
    ylim=c(0,410),yaxt="n",
    ylab="",xlab="Longitude",
     pch=19,col="gray",cex=BreedAbund/10)
legend(-110,410, legend="Breeding",bty="n")


##############################################################################
#
# Tarsus in relation to rainfall during non-breeding/breeding
#
##############################################################################
OVENnb<-read.csv("Data/Oven_non-breeding_LatLong.csv")
plot(Carib)
plot(WinterPrecip,add=TRUE)
points(-OVENnb$long,OVENnb$lat,cex=OVENnb$tarsus-min(OVENnb$tarsus,na.rm=TRUE)+1,pch=19)

Year<-1998:2015
nYears<-length(Year)

YearList<-RainRasters<-vector('list',nYears)
for(i in 1:nYears){
  YearList[[i]]<-list.files(paste0("D:/MC_Demography/Data/TRMM/",Year[i],"/"),full.names=TRUE)
}

for(i in 1:(nYears-1)){
  temp<-vector('list',12)
  for(m in 1:12){
    temp[[m]]<-raster(YearList[[i]][m])
    temp[[m]][temp[[m]]>10000]<-NA
  }
  RainRasters[[i]]<-stack(temp)
}

temp<-vector('list',8)
for(m in 1:12){
  temp[[m]]<-raster(YearList[[18]][m])
  temp[[m]][temp[[m]]>10000]<-NA
}
RainRasters[[18]]<-stack(temp)

Jan<-Feb<-Mar<-Apr<-May<-Jun<-Jul<-Aug<-Sep<-Oct<-Nov<-Dec<-list()
for(i in 1:nYears){
  Jan[[i]]<-RainRasters[[i]][[1]]
  Feb[[i]]<-RainRasters[[i]][[2]]
  Mar[[i]]<-RainRasters[[i]][[3]]
  Apr[[i]]<-RainRasters[[i]][[4]]
  May[[i]]<-RainRasters[[i]][[5]]
  Jun[[i]]<-RainRasters[[i]][[6]]
  Jul[[i]]<-RainRasters[[i]][[7]]
  Aug[[i]]<-RainRasters[[i]][[8]]
  Sep[[i]]<-RainRasters[[i]][[9]]
  Oct[[i]]<-RainRasters[[i]][[10]]
  Nov[[i]]<-RainRasters[[i]][[11]]
  Dec[[i]]<-RainRasters[[i]][[12]]
} 

winterRain<-vector('list',nYears-1)
for(i in 1:(nYears-1)){
winterRain[[i]]<-sum(Oct[[i]],Nov[[i]],Dec[[i]],Jan[[i+1]],Feb[[i+1]],Mar[[i+1]],Apr[[i+1]])
}



n<-array(NA,c(8,18))
colnames(n)<-c(1999:2015,"mean")
for(i in 1:17){
n[,i]<-extract(winterRain[[i]],cbind(-OVENnb$long,OVENnb$lat))
}
n[,18]<-apply(n[,1:17],1,median)



par(bty="l")
plot(OVENnb$tarsus[1:7]~log(n[1:7,18]),yaxt="n",ylab="Mean Tarsus Length",xlab="Log(Precipitation)",pch=19,cex=1.5)
text(log(n[,18]),OVENnb$tarsus-0.025,label=OVENnb[,1])
axis(2,las=2)
abline(lm((OVENnb$tarsus[1:7]~log(n[1:7,18]))))





























