#preprocessing for LST downscaling

library(terra)
library(lubridate)
library(ncdf4)
library(ncdf4.helpers)

homedir<-"/Users/bethanyblakely/Desktop/Analysis/CHEESEHEAD_local/CHEESEHEAD/scripts/LST_R" #working directory
datadir<-"/Users/bethanyblakely/Desktop/Analysis/LST/Data" #directory where the input data live
setwd(homedir)

#Set extents #####
#list of possible spatial extents
und<-ext(-92,-88, 45, 47) #underc
big.wi<-ext(-93,-87, 43, 47); med.wi<-ext(-90,-89, 45, 46) #large bounding boxes for initial cropping (to decrease runtime)
swth<-vect("stei_tree.kml"); swth.utm<-project(swth, "epsg:32616") #steigerwaldt/treehaven neon AOP, projected to UTM16
swth2<-buffer(swth, 12000); swth2.utm<-buffer(swth.utm, 12000) #steigerwaldt/treehaven AOP with a 12K buffer, projected to UTM16

wi<-swth2.utm #the area subset you want the ultimate output in. Must be UTM16

#temporal
timespan<-c(182:201) #DOYs to do all this for


#####

#read in, subset, and assign time to NLDAS: mosaic ####
setwd(datadir)
setwd("nldas/mosaic") #directory where the nldas data are

files.mo.all<-allfiles<-list.files()

#subset files in date range of interest based on file names
rdate<-as.POSIXlt(substr(files.mo.all, 22,25), format="%m%d", tz="UTC"); rdoy<-as.numeric(format(rdate, "%j"))
files.sub<-which(rdoy%in%timespan)

files.mo<-files.mo.all[files.sub]
allrast.mo<-list()

for (i in 1:length(files.mo)){
  
  if(i%%24==0){print(paste("DOY", timespan[i/24]))} #day counter for monitoring loop. No other function
  
  file<-files.mo[i]
  
  rst<-terra::rast(file)
  lst<-project(rst,"epsg:32616") #reproject to UTM zone 16
  
  lst.wisc<-crop(lst, wi,snap="out") #crop to area of intrest
  plot(lst.wisc)
  
  datestr<-substr(file, 18,25);timestr<-substr(file, 27,30) #subset out date and time from file name
  ts<-as.POSIXlt(paste(datestr, timestr, sep=""), format="%Y%m%d%H%M",tz="UTC") #...convert to posix
  
  time(lst.wisc)<-ts #...and assign as the time dimension for that image"
  
  allrast.mo[[i]]<-lst.wisc
  
}



brick.mosaic<-rast(allrast.mo) #aggregate into multilayer raster
#####

#read in, subset, and assign time to NLDAS: noah. Process identical to mosaic except for different file names#####
setwd(datadir)
setwd("nldas/noah")

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

#read in, subset, and assign time to NLDAS: vic. Same thing again #####
setwd(datadir)
setwd("nldas/vic")

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


#merge NLDAS datasets - mean and sd
nldas.brick<-mean(brick.vic, brick.mosaic, brick.noah)
nldas.sd<-stdev(brick.vic, brick.mosaic, brick.noah)

#####


#read in goes####
  

  setwd(datadir)
  setwd("goes")
  
  #Note: missing a couple files; should have either 720 (for 30 days) or 744 (for 31)
  #actually have 718
  
  #function to correct values for unsigned integer problem: removes stretch,
  #"converts" raw values to what they would be if read correctly, then reapplies stretch
  cv<-function(x){
    nval<-((((x-190)/0.0025)+2^16)*0.0025)+190
  }
  
  hdrs<-list.files(pattern = "*.aux.xml") #get list of header files so we can exclude them later
  
  allfiles<-list.files()
  
  files.goes.all<-setdiff(allfiles, hdrs) #exclude hdr files
  
  rdate<-as.POSIXlt(substr(files.goes.all, 24,37), format="%Y%j%H%M%S", tz="UTC"); rdoy<-as.numeric(format(rdate, "%j"))
  
  files.sub<-which(rdoy%in%timespan)
  
  files.goes<-files.goes.all[files.sub]
  
  #read in sample file to reproject bounding box to GOES CRS
  pset<-terra::rast(files.goes[1])
  bbox<-vect(med.wi) #med.wi is a user-designated (in ~line 15) large lat/lon bounding box that acts as an initial clip for (natively CONUS) goes data
  crs(bbox)<-"EPSG:4326" #assign WGS84 projection to med.wi
  bbox.proj<-project(bbox, pset) #reproject this box into goes CRS
  rm(pset)
  
  #main loop
  allrast.goes<-list()
  par(mfrow=c(2,2))
  for (i in 1:length(files.goes)){
    
  if(i%%24==0){print(paste("DOY", timespan[i/24]))} #Day counter. Not necessary
  
  grst<-terra::rast(files.goes[i])
  #plot(grst) #entire CONUS, LST and QC ("DQF"). Slows down loop a lot.
  
  grst.crop<-crop(grst, bbox.proj)
  grst.utm<-project(grst.crop,"epsg:32616", method="near") #project to UTM zone 16
  grst.wi<-crop(grst.utm, swth2.utm)
  #plot(grst.wi) #cropped to user-selected large bounding box
  
  #separate out LST and QC
  glst.wi<-grst.wi$LST
  gqc.wi<-grst.wi$DQF
  
  #apply function fixing unsigned integer problem
  glst.wi<-cv(glst.wi)
  
  #NaN pixels with QC != best quality (leaves only cloud-free pixels)
  glst.wi[gqc.wi!=0]<-NA

  
  #extract and assign time dimension
  filetime<-substr(files.goes[i], 24, 37)
  ts.posix<-as.POSIXct(filetime, format="%Y%j%H%M%S", tz="UTC")
  time(glst.wi)<-ts.posix
  
  ##little plotting loop if you want to see how it's working
  #if(!all(is.nan(glst.wi))){
  #plot(glst.wi, main=as.character(time(glst.wi)))
  #plot(gqc.wi, main="qc")
  #}
  
  allrast.goes[[i]]<-glst.wi
  }
  
  #combine into multilayer raster
  brick.goes<-rast(allrast.goes)

#####
    

#Resample and realign#####
goesdat<-disagg(brick.goes, fact=2, method="near")
nldasdat<-project(nldas.brick, goesdat, method="near")

#cleanup
rm("allrast.mo", "allrast.no", "allrast.goes", "allrast.vic","brick.goes", "brick.mosaic", "brick.noah", "brick.vic")

#####
  