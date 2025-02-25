
#set home and data directories if not already existing from Preprocess_goes_nldas.R
if(!exists("homedir")){homedir<-"/Users/bethanyblakely/Desktop/Analysis/CHEESEHEAD_local/CHEESEHEAD/scripts/LST_R"}
if(!exists("datadir")){datadir<-"/Users/bethanyblakely/Desktop/Analysis/LST/Data"}

#Run prior scripts if needed####
setwd(homedir)
if(!exists("goes.filled")){
source("Preprocess_goes_nldas_3.R")
source("GoesGF.R")
}

#Get stuff you need to run this script####
setwd(homedir)
library(terra)
library(lubridate)
library(smoothr)
library(spatialEco)


#Set extents. Already done in preprocess_goes_nldas but you can do it here again if you want#####

# #spatial
# und<-ext(-92,-88, 45, 47) #UNDERC
# big.wi<-ext(-93,-87, 43, 47); med.wi<-ext(-90,-89, 45, 46) #Wisconsin larger and smaller portion respectively
# swth<-vect("ancillary/stei_tree.kml"); swth.utm<-project(swth, "epsg:32616") #Steigerwalt/Treehaven from NEON IOP
# swth2<-buffer(swth, 12000); swth2.utm<-buffer(swth.utm, 12000) #Steigerwalt/Treehaven with an additional 12km buffer
# 
# wi<-swth2.utm #assign the one you'll use
# 
# #temporal
# timespan<-c(166:274) #DOYs to do all this for; must have data for these obviously



#Get ecostress QC Lookup table####
setwd(datadir)
setwd("lut")

qc.lut<-read.csv("ECO2LSTE-001-SDS-QC-lookup.csv")
vals.qc<-qc.lut$Value[qc.lut$Mandatory.QA.flag%in%c("Pixel produced, best quality")]#,"Pixel produced, nominal quality", "Pixel produced, but cloud detected")]
vals.nom<-qc.lut$Value[qc.lut$Mandatory.QA.flag%in%c("Pixel produced, nominal quality")]#,"Pixel produced, nominal quality", "Pixel produced, but cloud detected")]
#vals.cld<-qc.lut$Value[qc.lut$Mandatory.QA.flag%in%c("Pixel produced, but cloud detected")]#,"Pixel produced, nominal quality", "Pixel produced, but cloud detected")]

cld.lut<-read.csv("ECO2CLD-001-SDS-CloudMask-lookup.csv")
vals.cld<-cld.lut$Value[cld.lut$Final.Cloud..either.one.of.bits.2..3..or.4.set =="no" & cld.lut$Cloud.Mask.Flag=="determined"]


#Get ecostress LST and QC files listed and temporally aligned####

setwd(datadir)
setwd("ecostress/lst")

files.ecs.all<-list.files()
rdate.lst<-as.POSIXlt(substr(files.ecs.all, 25,37), format="%Y%j%H%M%S", tz="UTC"); rdoy<-as.numeric(format(rdate.lst, "%j"))
files.sub<-which(rdoy%in%timespan)
files.ecs<-files.ecs.all[files.sub]

setwd(datadir)
setwd("ecostress/qc")

files.qc.all<-list.files()
rdate.qc<-as.POSIXlt(substr(files.qc.all, 24,36), format="%Y%j%H%M%S", tz="UTC"); rdoy<-as.numeric(format(rdate.qc, "%j"))
files.sub<-which(rdoy%in%timespan)
files.qc<-files.qc.all[files.sub]

setwd(datadir)
setwd("ecostress/cld")

files.cld.all<-list.files()
rdate.cld<-as.POSIXlt(substr(files.cld.all, 30,42), format="%Y%j%H%M%S", tz="UTC"); rdoy<-as.numeric(format(rdate.cld, "%j"))
files.sub<-which(rdoy%in%timespan)
files.cld<-files.cld.all[files.sub]

#Extract dates from filenames
rdate.1<-as.POSIXlt(substr(files.qc, 24,36), format="%Y%j%H%M%S", tz="UTC") 
rdate.2<-as.POSIXlt(substr(files.ecs, 25,37), format="%Y%j%H%M%S", tz="UTC")
rdate.3<-as.POSIXlt(substr(files.cld, 30,42), format="%Y%j%H%M%S", tz="UTC")

#get a standard set of scenes across products by date
sub.qc<-which(rdate.1%in%rdate.2)
files.qc<-files.qc[sub.qc]

sub.cld<-which(rdate.3%in%rdate.2)
files.cld<-files.cld[sub.cld]

length(files.qc)==length(files.ecs)
length(files.cld)==length(files.ecs) #time alignment check. should both eval to true


#Read in ecostress and QC it, and align it to the GOES grid, also retaining an ecostress-gridded one for later####

setwd(homedir)

allrast.ecs<-list(); allrast.fullres<-list()
for (i in 1:length(files.ecs)){
  
  setwd(datadir)
  setwd("ecostress/lst")
  erst<-terra::rast(files.ecs[i])
  
  ts<-substr(files.ecs[i], 25, 37)
  
  setwd(datadir)
  setwd("ecostress/qc")
  qc.match<-which(substr(files.qc, 24, 36)==ts)
  qc<-terra::rast(files.qc[qc.match])
  
  setwd(datadir)
  setwd("ecostress/cld")
  cld.match<-which(substr(files.qc, 24, 36)==ts)
  cld<-terra::rast(files.cld[cld.match])
  
  setwd(homedir)
  
  #this one for regressions - good data only
  elst<-erst
  elst[!qc%in%vals.qc&!qc%in%vals.nom]<-NA
  
  #this one for final; less strict about cloud
  elst.gen<-erst
  elst.gen[!cld%in%vals.cld|elst.gen<273]<-NA

  elst.utm<-project(elst,"epsg:32616"); elst.gen.utm<-project(elst.gen,"epsg:32616")
  
  
  #tests for data availability in ROI:
  
  crop.attempt<-try(elst.wi<-crop(elst.utm, swth2.utm))

  filetime<-substr(files.ecs[i], 25, 37)
  #filetime<-substr(files.goes[i], 48, 61)
  ts.posix<-as.POSIXct(filetime, format="%Y%j%H%M%S", tz="UTC")
  
  if(class(crop.attempt)=="try-error"){
    print(paste("no ecostress in domain for", ts.posix))
          next}
  
  #this "cov" is the proportion of cells that are NA
  cov<-1-(freq(is.na(elst.wi))$count[freq(is.na(elst.wi))$value==1]/ncell(elst.wi))
  if(length(cov)==0){cov<-1} # if length (cov) is zero it means no NA's could be found. Thus, it is assigned full coverage
  print(paste(cov, "coverage"))
  
  if(cov<0.05){
    print("less than 5% coverage; skipping")
    next}
  
  #once it's passed ROI tests, make a copy of the generous ("gen") version for more analyses
  elst.gen.wi<-crop(elst.gen.utm, swth2.utm)


  #full resolution output
  time(elst.wi)<-time(elst.gen.wi)<-ts.posix
  allrast.fullres[[i]]<-elst.gen.wi #elst.wi
  #plot(elst.gen.wi, main=paste(time(elst.gen.wi)))
  
  #GOES-resolution output
  elst.align<-project(elst, goesdat, method="average")
  plot(elst.align); time(elst.align)<-ts.posix
  allrast.ecs[[i]]<-elst.align
}

allrast.ecs<-Filter(Negate(is.null), allrast.ecs)
brick.ecs<-rast(allrast.ecs)

fullres.ecs<-Filter(Negate(is.null), allrast.fullres)

#Find matching GOES times. This will be ugly.. ####

goeslist<-rep(0, length(time(brick.ecs)))
for(t in 1:length(time(brick.ecs))){
  time.ecs<-time(brick.ecs)[t]
  ind<-which(abs(time(goes.filled)-time.ecs)==min(abs(time(goes.filled)-time.ecs)))
  goeslist[t]<-ind
}

goes.regress<-goes.filled[[goeslist]]

#Make QC-step GOES-ecostress regressions (on GOES grid, within-scene across space) ####

pf<-rep(0,dim(brick.ecs)[3])
par(mfrow=c(1,3))

for(i in 1:dim(brick.ecs)[3]){ #for each available ecostress scene...

  
  plot(goes.regress[[i]], main=paste("goes",time(goes.regress)[i])); plot(brick.ecs[[i]], main=paste("ecostress",time(brick.ecs)[i]))
  
  goeslyr<-as.array(goes.regress[[i]])
  ecoslyr<-as.array(brick.ecs[[i]])
  
  ecovec<-c(ecoslyr); goesvec<-c(goeslyr)
  mod<-lm(ecovec~goesvec)
  
  print(summary(mod)$r.squared); print(summary(mod)$coefficients[2,4])
  
  plot(ecovec~goesvec, main=paste("R2:", round(summary(mod)$r.squared,2), "p:", round(summary(mod)$coefficients[2,4],2)))
  
  if(summary(mod)$r.squared>0.25 & summary(mod)$coefficients[2,4]<0.001){pf[i]<-1}
  
}

ecs.good<-brick.ecs[[pf==1]] #lose about half to bad regressions


#Resample goes and ecostress to common 50m grid #####

##make a grid
ext(fullres.ecs[[1]]) #fullres.ecs is the NAN-filtered original ecostress data

#get extents that are multiples of 50m
coord<-rep(NA, 4)
for(i in 1:4){ #for each extent
  coord[i]<-floor((ext(fullres.ecs[[1]])[i])/50)*50
}
cellsx<-(coord[2]-coord[1])/50
cellsy<-(coord[4]-coord[3])/50
grid<-rast(resolution=50,nrows=cellsy,ncols=cellsx, ext=coord, crs=crs(fullres.ecs[[1]]), vals=rnorm(n=cellsy*cellsx))

#resample ecostress stack to this grid

ecs.gridalign.l<-list()
for(e in 1:length(fullres.ecs)){
  
  r<-fullres.ecs[[e]]
  r.align<-resample(r, grid)
  
  ecs.gridalign.l[[e]]<-r.align
  
}

ecs.gridalign<-rast(ecs.gridalign.l)

#resample goes to the grid

goesdat.resamp<-resample(goes.filled, grid, method="bilinear") #slow! don't plot this.


#Pair retained Ecostress scenes with goes scenes####

goodind<-which(time(ecs.gridalign)%in%time(ecs.good)) #index of kept scenes
ecs.grid<-ecs.gridalign[[goodind]]#subset ecostress

g<-goesdat.resamp[[goeslist]] #subset goes to all ecostress scenes...
goes.grid<-g[[goodind]]#and again to kept ecostress scenes

#Ankur's debias: subtracting mean of values####

#plot(ecs.grid); plot(goes.grid)
tdiff<-ecs.grid-goes.grid
mu<-mean(ecs.grid, na.rm=TRUE)
ecs.mudiff<-ecs.grid-mu
goes.mudiff<-goes.grid-mu
#mudiff<-mean(tdiff, na.rm=TRUE)
#goes.deb<-goes.grid+mudiff

#Ecostress-goes regressions: downscaling (on common 50m grid, pixel by pixel regression through time).####
par(mfrow=c(3,3))

goes.array<-as.array(goes.mudiff); ecs.array<-as.array(ecs.mudiff)

par(mfrow=c(2,2))

d1<-dim(goes.grid)[1]; d2=dim(goes.grid)[2]
nlyr<-dim(goes.array)[3]

slopes<-matrix(nrow=d1, ncol=d2, data=NA)
ints<-matrix(nrow=d1, ncol=d2, data=NA)

start<-Sys.time()


for(i in c(1:d1)){
  if(i%%10==0){print(i)}
  for (j in c(1:d2)){
    
    g<-goes.array[i,j,] #entire "stack" at x/y location i,j. i.e. this is a regression across time
    e<-ecs.array[i,j,]
    
    #if both have values for more than half the available timepoints (=ecostress scenes)...
    if(length(which(is.finite(e)))>(nlyr/2) & length(which(is.finite(g)))>(nlyr/2)){
    
    #keep all to see which ones get removed
    e.all<-e
    
    #remove all points >4c different
    diffs<-abs(g-e)
    e[diffs>4]<-NaN
    
    if(length(which(is.finite(e)))<2){next} #skip to next if <2 points remain
     
    m<-lm(g~e)
    #m.noint<-lm(goes.array[i,j,]~0+ecs.array[i,j,]); m<-m.noint
    
    if(i%%10==0 &j%%100==0){
    plot(g~e.all,col='red', main=paste("slope:",round(coef(m)[2], 2), "int:",round(coef(m)[1], 2)), ylim=c(-12,12), xlim=c(-12,12))
    points(g~e)
    abline(coef(m), col='blue');
    abline(0,1)
    }
    
    slopes[i,j]<-coef(m)[2]
    ints[i,j]<-coef(m)[1]
    
    rm('m')
  
    }#else {(print("x"))}
    
  }
}
end<-Sys.time()
print(end-start)

#Investigate and QC slopes and ints

template<-goesdat.resamp[[1]]
sloperast<-template; values(sloperast)<-slopes
intrast<-template; values(intrast)<-ints

sloperast[sloperast<0.7 |sloperast>1.3]<-NA
intrast[abs(intrast)>5]<-NA


#Apply regression and debias (rebias lol?) to downscale GOES based on patterns in ecostress####

bk<-goesdat.resamp #backup with full precision
values(goesdat.resamp)<-round(values(goesdat.resamp), 2) #reduce precision for faster processing

gstack<-list()
for(i in 1:nlyr(goesdat.resamp)){ 
  if(i%%50==0){print(i)}
  #"Take the GOES-NLDAS 1x1 km scene, resample to 50x50 m, then apply a 1x1 km (20x20 pixel) smoother" (goesdat_resamp is this)
  gdat<-goesdat.resamp[[i]]
  #"Subtract the multi-scene pixel bypixel mean Ecostress LST"
  gdat.deb<-gdat-mu
  #"multiple by slope and add intercept"
  gdat.scl<-(gdat.deb*sloperast)+intrast
  
  gdat.out<-gdat.scl+mu
  gstack[[i]]<-gdat.out
}

lst.downscale<-rast(gstack)
plot(lst.downscale, range=c(285,305), maxnl=48)

lst.downscale.c<-crop(lst.downscale, swth.utm)
par(mfrow=c(1,1))
animate(lst.downscale.c[[491:731]], pause=0.01, main=paste(time(lst.downscale.c[[491:731]])))

writeRaster(lst.downscale.c, filename="STEI_TREE_LST", filetype="GTiff", overwrite=TRUE)


#####

