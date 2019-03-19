% --------- flight_track.m ----------- 
% Script to generate and save (kml & csv) flight tracks for 8 wind directions.
% Created - February 17, 2019
% Updated - February 28, 2019
% Author - Brian Butterworth
% ------------------------------------ 

clear
wd={'w','nw','n','ne','e','se','s','sw'};
bear1=[270,315,0,45,90,135,180,225];
bear2     =[180,45,90,135,180,45,90,135];

incr=1000;%number of iterations used to calculate path between two points.
         %It is necessary because equations used follow great circles, and
         %will cause some offset if not iterated.
         
for n=1:8

%Center of Study Area:
clat = 45.945350;
clon = -90.280300;

R = 6367.137;% radius of Earth at center point latitude
d = 4.5; % distance to move laterally in finding initial point

    for k=1:incr;%increment distance for calculation 
        if k==1
            [lat_tmp,lon_tmp] = coordcalc(clat,clon,d/incr,bear1(n));
        else
           	[lat_tmp,lon_tmp] = coordcalc(lat_tmp,lon_tmp,d/incr,bear1(n));
        end
    end
        
d=15; % distance to move down for finding initial point

    for k=1:incr;%increment distance for calculation 
           	[lat_tmp,lon_tmp] = coordcalc(lat_tmp,lon_tmp,d/incr,bear2(n));
    end
    
    lat(1)=lat_tmp;
    lon(1)=lon_tmp;

% distance to travel each leg
d=[NaN 2 2 2 2 2 30 2 2 2 2 2];

% altitude at each waypoint
order=[1 4 5 8 9 12 11 10 7 6 3 2];

waypoint=[1 2 2 1 3 4 4 3 5 6 6 5 7 8 8 7 9 10 10 9 11 12 12 11];  
alt=    [400 400 100 100 400 400 100 100 400 400 100 100 400 400 100 100 400 400 100 100 400 400 100 100];

% change in bearing between subsequent wind directions
bear_change=[NaN 225 45 45 45 225 45 45];


if n==1
   bear=[NaN 90 90 90 90 90 0 -90 -90 -90 -90 -90];
    ind1=find(bear==-90 | bear==90);
elseif n==2 | n==6
%     
     bear=bear+bear_change(n);
     bear(ind1)=bear(ind1)+180;
     bear(bear>360)=bear(bear>360)-360;
     bear(bear<0)=bear(bear<0)+360;
 else
     bear=bear+bear_change(n);
     bear(bear>360)=bear(bear>360)-360;
     bear(bear<0)=bear(bear<0)+360;

end


for j=2:length(d);
    
    for k=1:incr;
        lat1=lat(j-1);
        lon1=lon(j-1);
        if k==1
            [lat_tmp,lon_tmp] = coordcalc(lat1,lon1,d(j)/incr,bear(j));
        else
           	[lat_tmp,lon_tmp] = coordcalc(lat_tmp,lon_tmp,d(j)/incr,bear(j));
        end
        
    end
    
    lat(j)=lat_tmp;
    lon(j)=lon_tmp;

end


for k=1:length(waypoint)
   ind=find( order==waypoint(k));
tlat(k)=lat(ind);
tlon(k)=lon(ind);
end


%-------------------------------------------------------------------------
% If we want exactly 10 x 10 km on all sides we'll need to come up with a
% more precise way to do this. But for now, for the ease of the pilots I am
% keeping the longitudes of opposite waypoints fixed for both north/south
% points (only applicable for west and east wind sets of waypoints). 
% I am doing so because I assume it will be easier for the pilots,
% when switching between waypoints. Also, I'm not sure, but it's possible that 
% it is easier to follow due north, compared to some slight fraction off of
% due north. It's going to mean about 135m shorter top leg across. 
if n==1 | n==5
    i1=find(waypoint==1);i1=i1(1);
    i2=find(waypoint==2);
    tlon(i2)=tlon(i1);

    i1=find(waypoint==4);i1=i1(1);
    i2=find(waypoint==3);
    tlon(i2)=tlon(i1);
    
    i1=find(waypoint==5);i1=i1(1);
    i2=find(waypoint==6);
    tlon(i2)=tlon(i1);
 
    i1=find(waypoint==8);i1=i1(1);
    i2=find(waypoint==7);
    tlon(i2)=tlon(i1);

    i1=find(waypoint==9);i1=i1(1);
    i2=find(waypoint==10);
    tlon(i2)=tlon(i1);
end
%-------------------------------------------------------------------------

airport_lat= 45.630244;
airport_lon=-89.486242;

[dist_arcend,bear_end]=distance(lat(end),lon(end),airport_lat,airport_lon);
[dist_arcbeg,bear_beg]=distance(airport_lat,airport_lon,lat(1),lon(1));
dist_end=dist_arcend*111.15128505806946;
dist_beg=dist_arcbeg*111.15128505806946;

D=[dist_beg,d,dist_end];
ALT=[0,alt,0];
LAT=[airport_lat,tlat,airport_lat];
LON=[airport_lon,tlon,airport_lon];
TIME=round((D/(160*1.852))*60,1);
WAYPOINT=[0,waypoint,0]';

LAT_DM=degrees2dm(LAT);LAT_DM(:,2)=round(LAT_DM(:,2),2);
LON_DM=degrees2dm(LON);LON_DM(:,2)=round(LON_DM(:,2),2);

% SAVE AS KML
filepath=['/Users/brian/Documents/cheesehead/waypoints/way_' wd{n} '.kml'];
kmlwritepoint(filepath,LAT,LON,ALT);

% SAVE AS CSV
for k=1:length(LAT_DM)
    dec=num2str(LAT_DM(k,2)-floor(LAT_DM(k,2))) ;
    if length(dec)>=4
        dec=dec(3:4);
    else
        dec=[dec(3) '0'];
    end
    lat_temp{k}=['N ' sprintf('%02d',LAT_DM(k,1)) ' ' sprintf('%02d',floor(LAT_DM(k,2))) '.' dec ''''];

    dec=num2str(LON_DM(k,2)-floor(LON_DM(k,2))) ;
    if length(dec)>=4
        dec=dec(3:4);
    elseif length(dec)==1
        dec=['00'];
    else
        dec=[dec(3) '0'];
    end
    lon_temp{k}=['W ' sprintf('%03d',abs(LON_DM(k,1))) ' ' sprintf('%02d',floor(LON_DM(k,2))) '.' dec ''''];

end
fname=['/Users/brian/Documents/cheesehead/waypoints/waypoints_' wd{n} '.csv'];
LAT=lat_temp';
LON=lon_temp';
ALT=[num2cell(ALT)]';
T = table(WAYPOINT,LAT,LON,ALT);
writetable(T,fname)

end

copyname=which(mfilename)
copyfile(copyname,['/Users/brian/Documents/cheesehead/waypoints/flight_track_' date '.m']);
coordpath=which('coordcalc')
copyfile(coordpath,['/Users/brian/Documents/cheesehead/waypoints/coordcalc_' date '.m']);
