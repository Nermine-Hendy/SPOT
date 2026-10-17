function r=sun(C,t)
%SUN Low-precision geometric solar ephemeris; MOD equatorial -> GCRF.
% Approximate solar longitude/obliquity; no aberration/refraction, metre AU scale.
d=juliandate(t(:))-2451545; L=mod(280.460+.9856474*d,360); g=mod(357.528+.9856003*d,360);
lambda=L+1.915*sind(g)+.020*sind(2*g); obliq=23.439-.0000004*d;
range=C.Constants.AU_m*(1.00014-.01671*cosd(g)-.00014*cosd(2*g));
v=[cosd(lambda)';(cosd(obliq).*sind(lambda))';(sind(obliq).*sind(lambda))'];
Q=SPOT.precession(t); r=zeros(3,numel(t));
for k=1:numel(t), r(:,k)=Q(:,:,k)'*v(:,k)*range(k); end
end
