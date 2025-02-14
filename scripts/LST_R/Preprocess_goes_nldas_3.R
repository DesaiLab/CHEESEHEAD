#preprocessing for LST downscaling

library(terra)
library(lubridate)
library(ncdf4)
library(ncdf4.helpers)

setwd("/Users/bethanyblakely/Desktop/Analysis/LST")

#Set extents #####
#spatial
und<-ext(-92,-88, 45, 47)
big.wi<-ext(-93,-87, 43, 47); med.wi<-ext(-90,-89, 45, 46)
swth<-vect("ancillary/stei_tree.kml"); swth.utm<-project(swth, "epsg:32616")
swth2<-buffer(swth, 12000); swth2.utm<-buffer(swth.utm, 12000)

wi<-swth2.utm #the one you'll use

#temporal
timespan<-c(166:274)#304 #245 #DOYs to do all this for
#timespan<-c(182:212) #July

#####

#read in NLDAS: mosaic ####
setwd("/Users/bethanyblakely/Desktop/Analysis/LST")
setwd("Data/NLDAS/Incoming2/mosaic")

files.mo.all<-allfiles<-list.files()

#subset files in date range of interest
rdate<-as.POSIXlt(substr(files.mo.all, 22,25), format="%m%d", tz="UTC"); rdoy<-as.numeric(format(rdate, "%j"))
files.sub<-which(rdoy%in%timespan)

files.mo<-files.mo.all[files.sub]
allrast.mo<-list()

for (i in 1:length(files.mo)){
  
  if(i%%24==0){print(paste("DOY", timespan[i/24]))} #day counter
  
  file<-files.mo[i]
  
  rst<-terra::rast(file)
  lst<-project(rst,"epsg:32616")
  
  lst.wisc<-crop(lst, wi,snap="out")
  plot(lst.wisc)
  
  datestr<-substr(file, 18,25);timestr<-substr(file, 27,30)
  ts<-as.POSIXlt(paste(datestr, timestr, sep=""), format="%Y%m%d%H%M",tz="UTC")
  
  time(lst.wisc)<-ts
  
  allrast.mo[[i]]<-lst.wisc
  
}



brick.mosaic<-rast(allrast.mo)
#####

#read in NLDAS: noah #####
setwd("/Users/bethanyblakely/Desktop/Analysis/LST")
setwd("Data/NLDAS/Incoming2/noah")

files.no.all<-allfiles<-list.files()

#subset files in date range of interest
rdate<-as.POSIXlt(substr(files.no.all, 23,26), format="%m%d", tz="UTC"); rdoy<-as.numeric(format(rdate, "%j"))
files.sub<-which(rdoy%in%timespan)

files.no<-files.no.all[files.sub]
allrast.no<-list()

for (i in 1:length(files.no)){
  
  if(i%%24==0){print(paste("DOY", timespan[i/24]))} #day counter
  
  file<-files.no[i]
  
  rst<-terra::rast(file)
  lst<-project(rst,"epsg:32616")

  lst.wisc<-crop(lst, wi, snap="out")
  #plot(lst.wisc)

  datestr<-substr(file, 19,26);timestr<-substr(file, 28,31)
  ts<-as.POSIXlt(paste(datestr, timestr, sep=""), format="%Y%m%d%H%M",tz="UTC")

  time(lst.wisc)<-ts
  
  allrast.no[[i]]<-lst.wisc
  
}

brick.noah<-rast(allrast.no)

#####

#read in NLDAS: vic #####
setwd("/Users/bethanyblakely/Desktop/Analysis/LST")
setwd("Data/NLDAS/Incoming2/vic")

files.vic.all<-allfiles<-list.files()

#subset files in date range of interest
rdate<-as.POSIXlt(substr(files.vic.all, 22,25), format="%m%d", tz="UTC"); rdoy<-as.numeric(format(rdate, "%j"))
files.sub<-which(rdoy%in%timespan)

files.vic<-files.vic.all[files.sub]
allrast.vic<-list()

for (i in 1:length(files.vic)){
  
  if(i%%24==0){print(paste("DOY", timespan[i/24]))} #day counter
  
  file<-files.vic[i]
  
  rst<-terra::rast(file)
  lst<-project(rst,"epsg:32616")
  
  lst.wisc<-crop(lst, wi,snap="out")
  #plot(lst.wisc)
  
  datestr<-substr(file, 18,25);timestr<-substr(file, 27,30)
  ts<-as.POSIXlt(paste(datestr, timestr, sep=""), format="%Y%m%d%H%M",tz="UTC")
  
  time(lst.wisc)<-ts
  
  allrast.vic[[i]]<-lst.wisc
  
}

brick.vic<-rast(allrast.vic)

#####



#merge NLDAS datasets
nldas.brick<-mean(brick.vic, brick.mosaic, brick.noah)

#nldas.sum<-(brick.vic+brick.mosaic+brick.noah)
#nldas.mean<-nldas.sum/3

nldas.sd<-stdev(brick.vic, brick.mosaic, brick.noah)
#####


#read in goes####
  

  setwd("/Users/bethanyblakely/Desktop/Analysis/LST")
  setwd("Data/goesday")
  
  #missing a couple files; should have either 720 (for 30 days) or 744 (for 31)
  #actually have 718
  
  #function to correct values for unsigned integer problem: removes stretch,
  #"converts" raw values to what they would be if read correctlty, then reapplies stretch
  cv<-function(x){
    nval<-((((x-190)/0.0025)+2^16)*0.0025)+190
  }
  
  hdrs<-list.files(pattern = "*.aux.xml")
  
  
  allfiles<-list.files()
  
  files.goes.all<-setdiff(allfiles, hdrs)
  
  rdate<-as.POSIXlt(substr(files.goes.all, 24,37), format="%Y%j%H%M%S", tz="UTC"); rdoy<-as.numeric(format(rdate, "%j"))
  
  files.sub<-which(rdoy%in%timespan)
  
  files.goes<-files.goes.all[files.sub]
  
  #read in sample file to create target projection for initial clip
  pset<-terra::rast(files.goes[1])
  bbox<-vect(med.wi)
  crs(bbox)<-"EPSG:4326"
  bbox.proj<-project(bbox, pset)
  rm(pset)
  
  #main loop
  allrast.goes<-list()
  par(mfrow=c(2,2))
  for (i in 1:length(files.goes)){
    
  if(i%%24==0){print(paste("DOY", timespan[i/24]))} #day counter
  
  grst<-terra::rast(files.goes[i])
  

  
  grst.crop<-crop(grst, bbox.proj)
  grst.utm<-project(grst.crop,"epsg:32616", method="near") #project to UTM zone 16
  grst.wi<-crop(grst.utm, swth2.utm)
  
  
  
  glst.wi<-grst.wi$LST
  gqc.wi<-grst.wi$DQF
  
  
  #apply function fixing unsigned integer problem
  glst.wi<-cv(glst.wi)
  
  glst.wi[gqc.wi!=0]<-NA
  glst.wi[glst.wi<280]<-NA
  
  filetime<-substr(files.goes[i], 24, 37)
  #filetime<-substr(files.goes[i], 48, 61)
  ts.posix<-as.POSIXct(filetime, format="%Y%j%H%M%S", tz="UTC")
 
  time(glst.wi)<-ts.posix
  
  
  #if(!all(is.nan(glst.wi))){
  #plot(glst.wi, main=as.character(time(glst.wi)))
  #plot(gqc.wi, main="qc")
  #}
  
  allrast.goes[[i]]<-glst.wi
  }
  
  brick.goes<-rast(allrast.goes)

#####
    

#Resample and realign#####
goesdat<-disagg(brick.goes, fact=2, method="near"); #why does Ankur use NN here?
nldasdat<-project(nldas.brick, goesdat, method="near")

##This currently breaks R  
#setwd("/Users/bethanyblakely/Desktop/Analysis/LST")
#writeRaster(goesdat, "GOES_July_2019.nc", overwrite=TRUE)
#writeRaster(nldasdat, "NLDAS_July_2019.nc", overwrite=TRUE)

rm("allrast.mo", "allrast.no", "allrast.goes", "allrast.vic","brick.goes", "brick.mosaic", "brick.noah", "brick.vic")
#####
  