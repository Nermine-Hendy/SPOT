function Q = precession(t)
%PRECESSION IAU-1976 precession approximation, J2000 -> mean of date.
% Used only by the explicitly approximate frame/Sun models, not TLE SGP4.
T=(juliandate(t(:))-2451545)/36525;
zeta=(2306.2181*T+.30188*T.^2+.017998*T.^3)/3600;
z=(2306.2181*T+1.09468*T.^2+.018203*T.^3)/3600;
theta=(2004.3109*T-.42665*T.^2-.041833*T.^3)/3600;
Rz=@(a)[cosd(a) -sind(a) 0;sind(a) cosd(a) 0;0 0 1];
Ry=@(a)[cosd(a) 0 sind(a);0 1 0;-sind(a) 0 cosd(a)];
Q=zeros(3,3,numel(t));
for k=1:numel(t), Q(:,:,k)=Rz(z(k))*Ry(-theta(k))*Rz(zeta(k)); end
end
