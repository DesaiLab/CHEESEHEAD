function [lat_tmp,lon_tmp] = coordcalc(LAT1,LON1,d,bear)

R=6367.137;%radius at clat

lat_tmp = asind( sind(LAT1) * cosd((d/R)*180/pi) + cosd(LAT1) * sind((d/R)*180/pi) * cosd(bear));
lon_tmp = LON1 + atan2d(sind(bear)*sind((d/R)*180/pi)*cosd(LAT1), cosd((d/R)*180/pi)- sind(LAT1)*sind(lat_tmp));

end