function result = geodetic(C,value,direction)
%GEODETIC WGS84 forward/inverse; columns are individual points.
a=C.Constants.EarthA_m; e2=C.Constants.EarthEccentricity^2;
if strcmp(direction,'to_ecef')
    lat=value(1,:); lon=value(2,:); h=value(3,:);
    N=a./sqrt(1-e2*sind(lat).^2);
    result=[(N+h).*cosd(lat).*cosd(lon);(N+h).*cosd(lat).*sind(lon); ...
        (N*(1-e2)+h).*sind(lat)];
else
    x=value(1,:); y=value(2,:); z=value(3,:); p=hypot(x,y);
    lat=atan2(z,p*(1-e2)); lon=atan2(y,x);
    for k=1:15
        N=a./sqrt(1-e2*sin(lat).^2); latNew=atan2(z+e2*N.*sin(lat),p);
        if max(abs(latNew-lat),[],'all')<1e-13, lat=latNew; break; end
        lat=latNew;
    end
    N=a./sqrt(1-e2*sin(lat).^2); h=p./max(cos(lat),eps)-N;
    pole=p<1e-7; h(pole)=abs(z(pole))-a*sqrt(1-e2);
    result=[rad2deg(lat);rad2deg(lon);h];
end
end
