function B=surfaceReflection(S,R,sunDirection,observerDirection,solarIrradiance,range)
%SURFACEREFLECTION Unresolved convex facets, Lambert + reciprocal modified Phong.
% All directions point OUTWARD from the target. Surface area is physical face
% area, not projected area. Both cosines appear exactly once in irradiance.
% f = rho_d/pi + rho_s*(p+2)/(2*pi)*max(r.v,0)^p, only on front hemisphere.
% rho_d+rho_s<=1 gives an energy-conserving empirical BRDF. rho_s is a
% lobe coefficient, not an angle-independent hemispherical albedo.
n=size(sunDirection,2); nf=size(S.NormalBody,2);
B.DiffuseIrradiance_W_m2=zeros(1,n); B.SpecularIrradiance_W_m2=zeros(1,n);
B.IncidenceCosine=zeros(nf,n); B.EmissionCosine=zeros(nf,n);
B.SpecularAlignment=zeros(nf,n);
for k=1:n
    N=R(:,:,k)*S.NormalBody;
    ci=N'*sunDirection(:,k); co=N'*observerDirection(:,k);
    reflected=2*N.*ci'-sunDirection(:,k);
    alignment=max(0,min(1,reflected'*observerDirection(:,k)));
    support=ci>0 & co>0;
    geometric=S.Area_m2(:).*max(0,ci).*max(0,co).*support;
    fd=S.DiffuseReflectance(:)/pi;
    fs=S.SpecularReflectance(:).*(S.SpecularExponent(:)+2)/(2*pi).*alignment.^S.SpecularExponent(:);
    factor=solarIrradiance(k)/range(k)^2;
    B.DiffuseIrradiance_W_m2(k)=factor*sum(geometric.*fd);
    B.SpecularIrradiance_W_m2(k)=factor*sum(geometric.*fs);
    B.IncidenceCosine(:,k)=ci; B.EmissionCosine(:,k)=co; B.SpecularAlignment(:,k)=alignment;
end
B.Irradiance_W_m2=B.DiffuseIrradiance_W_m2+B.SpecularIrradiance_W_m2;
end
