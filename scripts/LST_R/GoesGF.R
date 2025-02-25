#Gapfilling regressions for GOES/NLDAS

#set home and data directories if not already existing from Preprocess_goes_nldas.R
if(!exists("homedir")){homedir<-"/Users/bethanyblakely/Desktop/Analysis/CHEESEHEAD_local/CHEESEHEAD/scripts/LST_R"}
if(!exists("datadir")){datadir<-"/Users/bethanyblakely/Desktop/Analysis/LST/Data"}


setwd("homedir")

goesdat.bk<-goesdat #make a backup of ungapfilled goes data; mostly for development purposes

#NaN LSTs out of reasonable range
goesdat[goesdat>310|goesdat<273]<-NA


# #an example timeseries plot
# nts<-as.numeric(unname(nldasdat[10,10,]))
# gts<-as.numeric(unname(goesdat[10,10,]))
# timex<-as.numeric(format(time(nldasdat), "%j"))+(as.numeric(format(time(nldasdat), "%H"))/24)
# timex2<-as.numeric(format(time(goesdat), "%j"))+(as.numeric(format(time(goesdat), "%H"))/24)
# 
# par(mfrow=c(1,1))
# plot(nts~timex, ylab="LST", type="l", xlab="DOY"); points(gts~timex2, col='red')
# legend(max(timespan)-5, min(gts, na.rm=TRUE)+10, legend=c("GOES", "NLDAS"), pch=c("o", "---"), col=c("red", "black"), bty="n")


#get rasters into same temporal dimensions

both<-which(timex%in%timex2) #timex: nldas tme; timex2: goes time; "which nldas times are found in goes times"
nldasdat2<-nldasdat[[both]]

timex3<-as.numeric(format(time(nldasdat2), "%j"))+(as.numeric(format(time(nldasdat2), "%H"))/24) #convert to decimal DOY
both2<-which(timex2%in%timex3) #which goestime is in the trimmed nldas time
goesdat2<-goesdat[[both2]]

#Check: should evaluate true (rounded time because GOES is timestamped ~1min after each hour)
identical(date(time(goesdat2)), date(time(nldasdat2))) & identical(hour(time(goesdat2)), hour(time(nldasdat2)))
goesdat<-goesdat2 #hacky way to not have to change names later. 


##QC and simple interpolation steps. This makes a lot of plots and console output you may or may not want

for (i in 1:dim(nldasdat2)[3]){
  
  g<-goesdat[[i]]
  n<-nldasdat2[[i]]
  
  #coverage: proportion goes pixels with data
  cov<-length(which(!is.na(as.array(g))))/length(as.array(g))
  
  
  if(cov>0.1){ #only attempt for scenes with coverage >10%
  
  if(i!=1&i!=dim(nldasdat2)[3]){ #if not first or last
  
    #temporal interpolation: interpolate across one hour gaps
    
    #get previous, next, and to-be filled scenes
    prev<-goesdat[[i-1]]
    nxt<-goesdat[[i+1]]
    fill<-goesdat[[i]]
    gap<-which(is.na(values(fill)))
      
    #if there are some data in both the previous and next scene, and the to-be-filled scene has some but not all pixels...
    if(!all(is.na(values(prev)))&!all(is.na(values(nxt)))&!all(is.na(values(fill)))&length(gap)!=0){
        
  
      fill[gap]<-(prev[gap]+nxt[gap])/2 #mean of prev and next timestep anywhere to-be-filled has missing data
      
      #see how it's going
      par(mfrow=c(1,2))
      plot(goesdat[[i]], main=paste(time(g))); plot(fill, main=paste("timeinterp",time(goesdat)[i]))
      
      values(goesdat[[i]])<-values(fill)
      g<-goesdat[[i]]
        
      }else{print(paste("skipping temporal interp for", time(goesdat)[i]))}
      
    }
    
  rm("fill", "prev", "nxt")
  
  
  
  #Spatial interpolation

  #get the coverage post temporal interpolation
  cov2<-length(which(!is.na(as.array(g))))/length(as.array(g))
  
  if(cov2>0.9 & cov2<1){ #if we're missing 10% of the data or less (don't want to be filling big holes this way)...
    
    par(mfrow=c(1,2))
    plot(g, main=paste(time(g))); #plot(n)
    
    #moving window smoother applied only to NA values
    plot(focal(g, na.policy="only", fun="mean"), main=paste("spatinterp", time(goesdat)[i]))
    values(goesdat[[i]])<-values(focal(g, na.policy="only", fun="mean"))
    
    g<-goesdat[[i]]
    
  }else{print(paste("skipping spatial interp for", time(goesdat)[i]))}
  
  }else{
    
    print(paste("goes data coverage too low for", time(g), ":", round(cov,2), "NANing any remaining values"))
  
    values(goesdat[[i]])<-NA
    
    }
  
}




##Build regressions#####

allslopes<-array(data=NA, dim=c(dim(goesdat)[1:2], 24)) #same spatial dimensions as goes, 24 hours
allints<-array(data=NA, dim=c(dim(goesdat)[1:2], 24)) #same spatial dimensions as goes, 24 hours
alldebias<-array(data=NA, dim=c(dim(goesdat)[1:2], 24)) #same spatial dimensions as goes, 24 hours
allslopes.noint<-array(data=NA, dim=c(dim(goesdat)[1:2], 24)) #same spatial dimensions as goes, 24 hours

#Loop over hours
for (h in 1:24){
 
#set hour 
hr<-h-1 #because hours are starting at 0 
print(paste("doing regressions for hour", hr, "UTC"))

#subset by that hour across days
goeshr.1<-goesdat[[hour(time(goesdat))==hr]] 
nldashr.1<-nldasdat[[hour(time(nldasdat))==hr]]

#plot(goeshr)

both<-which(date(time(nldashr.1))%in%date(time(goeshr.1)))
nldashr<-nldashr.1[[both]]; #rm(nldashr.1)

both2<-which(date(time(goeshr.1))%in%date(time(nldashr)))
goeshr<-goeshr.1[[both2]]

if(all(date(time(goeshr))==date(time((nldashr))))){print(paste("times match for hour", h-1))} #should eval to true


#plot(nldashr)

terra::plot(goeshr-nldashr)

#make arrays for easier subsetting
goes<-as.array(goeshr); nldas<-as.array(nldashr)

#make holders for regression params
slopes<-(matrix(nrow=dim(goes)[1], ncol=dim(goes)[2])); slopes.noint<-(matrix(nrow=dim(goes)[1], ncol=dim(goes)[2]))
ints<-(matrix(nrow=dim(goes)[1], ncol=dim(goes)[2]))
debias<-(matrix(nrow=dim(goes)[1], ncol=dim(goes)[2]))


par(mfrow=c(2,2))
#main loop
for (i in 1:nrow(goes)){
  #print(paste("row",i))
  print(i)
  for (j in 1:ncol(goes)){
    gpix<-goes[i,j,]; npix<-nldas[i,j,]
    #print(paste("there are", (length(which(!is.na(gpix)))), "good goes values of", length(gpix), "possible"))
   
     if(!all(is.nan(gpix))&!all(is.nan(npix))){
      
      #do a difference of the means
       
      gmean<-mean(gpix, na.rm=TRUE); nmean<-mean(npix[!is.na(gpix)], na.rm=TRUE) #now only usese nldas pix where goes pix exist
      d.means<-nmean-gmean; debias[i,j]<-d.means #calc and extract difference in the means
      
      #print(paste("debias is", d.means))
      
      #appply the difference so NLDAS has the same mean as GOES
      npix.mod<-npix-(d.means) #appply the difference so NLDAS has the same mean as GOES
      #npix.mod<-npix
      
      #calc and extract coefficients
      cf<-coef(lm(gpix~npix.mod))
      slopes[i,j]<-cf[2]; ints[i,j]<-cf[1] 
      
      cf.noint<-coef(lm(gpix~0+npix.mod))
      slopes.noint[i,j]<-cf.noint[1]
      
      #plot
      if(j%%2==0){  #All plots are doubled so we only plot every other.
        plot(gpix~npix.mod, main=paste("hour ", hr, "; (", i,",", j,")", sep=""), ylim=c(280, 310), xlim=c(280,310)); 
        abline(0,1)
        abline(cf,col="blue")
        abline(0, cf.noint[1], col="gray")}
      
    }else{
    
    slopes[i,j]<-NA; ints[i,j]<-NA
    }
    
  }
}

allslopes[,,h]<-slopes; allslopes.noint[,,h]<-slopes.noint
allints[,,h]<-ints
alldebias[,,h]<-debias

}


#####

##Apply regressions#####

#To convert: 
#goes = slope*(nldas-debias)+int

goes.restack<-array(dim=dim(goesdat))

for(i in 1:dim(goesdat)[3]){
  par(mfrow=c(2,2))
  
  if(hour(time(goesdat[[i]]))==0){print(paste("starting", date(time(goesdat[[i]]))))}

  goes.slice.og<-as.array(goesdat[[i]]); goes.slice.og<-goes.slice.og[,,1]
  nldas.slice<-as.array(nldasdat2[[i]]); nldas.slice<-nldas.slice[,,1]
  

  hour<-hour(time(goesdat[[i]]))+1
  
  
  slope<-allslopes.noint[,,hour];#slope<-allslopes[,,hour];
  debias<-alldebias[,,hour]
  int<-0#int<-allints[,,hour]
  
  ind<-which(is.na(goes.slice.og))
  goes.fill<-slope*(nldas.slice-debias)+int
  
  goes.slice<-goes.slice.og
  
  #Filling the gaps
  goes.slice[ind]<-goes.fill[ind]
  
  goes.filld<-goes.slice
  
  #if there's not much goes and the means are way different, go with fill
  subflag<-0
  cov<-1-(length(which(is.na((goes.slice.og))))/length(goes.slice.og))
  ogmean<-mean(goes.slice.og, na.rm=TRUE); #ogsd<-sd(goes.slice.og, na.rm=TRUE)
  nlmean<-mean(goes.fill, na.rm=TRUE)
  if(cov>0 & ((cov<0.3&abs(nlmean-ogmean)>3)|(abs(nlmean-ogmean)>5))){goes.slice<-goes.fill;  subflag<-1} #if cov=0, all inds will be filled anyway

  
  
  all.lims<-c(range(goes.slice.og, na.rm=TRUE),range(goes.slice, na.rm=TRUE), range(goes.fill, na.rm=TRUE), range(goes.filld, na.rm=TRUE))
  all.lims<-all.lims[is.finite(all.lims)]
  range<-c(min(all.lims)-0.1, max(all.lims)+0.1)
  
  #range<-c(280,305)
  
  
  if(!all(is.na(goes.slice.og))){plot(rast(goes.slice.og), main=paste("og",i), range=range)}else{plot(1:10, col="white", main="NO GOES")}
  
  if(subflag==0){
  plot(rast(debias), main="debias")
  plot(rast(goes.fill), main="fill", range=range); 
  plot(rast(goes.slice), range=range, main=paste("final:", date(time(goesdat)[i]), paste(hour(time(goesdat)[i]), ":00", sep="")))
  }
  
  if(subflag==1){
    plot(rast(goes.slice), main="final")
    plot(rast(goes.fill), main="fill; used for final data", range=range); 
    plot(rast(goes.filld), range=range, main=paste("if not filled", date(time(goesdat)[i]), paste(hour(time(goesdat)[i]), ":00", sep="")))
  }
    
  
  goes.restack[,,i]<-goes.slice
  
}

goes.filled<-rast(goes.restack)

time(goes.filled)<-time(goesdat);crs(goes.filled)<-crs(goesdat);ext(goes.filled)<-ext(goesdat)
#####


#example plot showing fill
#an example plot
par(mfrow=c(1,1))
nts<-as.numeric(unname(nldasdat[10,10,]))
gts<-as.numeric(unname(goesdat[10,10,]))
gtf<-goes.restack[10,10,]
timex<-as.numeric(format(time(nldasdat), "%j"))+(as.numeric(format(time(nldasdat), "%H"))/24)
timex2<-as.numeric(format(time(goesdat), "%j"))+(as.numeric(format(time(goesdat), "%H"))/24)

plot(nts~timex, type="l", ylab="lst", xlab="DOY"); 
points(gts~timex2, col='red')
points(gtf~timex2, col='blue', pch="+")

plot(gts~timex2, col='red', ylim=c(280,310))
lines(gtf~timex2, col='blue')



