function [trans,gain]=atmosphere(C,t,el,isGround,sourceKey)
%ATMOSPHERE Extinction and reproducible time-correlated, unit-mean scintillation.
trans=ones(size(el)); gain=trans;
if ~C.Atmosphere.Enabled || ~any(isGround(:)), return; end
eff=max(el,C.Atmosphere.MinElevation_deg);
switch C.Atmosphere.AirmassModel
    case 'secant', X=1./max(sind(eff),1e-6);
    case 'kasten_young', X=1./(sind(eff)+.50572*(eff+6.07995).^(-1.6364));
end
trans(isGround)=10.^(-.4*C.Atmosphere.Extinction_mag_airmass*X(isGround));
S=C.Atmosphere.Scintillation;
if ~S.Enabled || S.Sigma0==0, return; end
q=seconds(t(:)-C.Time.StartUTC); h=S.KnotStep_s;
assert(all(q>=-1e-7),'SPOT:Scintillation','Scintillation times precede the acquisition epoch.');
q=max(0,q); nk=ceil(max(q)/h)+2;
rs=RandStream('mt19937ar','Seed',mod(double(C.RandomSeed)+104729*double(sourceKey),2^32));
rho=exp(-h/S.CorrelationTime_s); z=randn(rs,nk,1);
for i=2:nk, z(i)=rho*z(i-1)+sqrt(1-rho*rho)*z(i); end
i=floor(q/h)+1; f=q/h-floor(q/h);
% Normalize interpolant variance, preserving E[exp(sigma*z-sigma^2/2)]=1.
zq=((1-f).*z(i)+f.*z(i+1))./sqrt((1-f).^2+f.^2+2*rho*f.*(1-f));
sigma=S.Sigma0./max(sind(eff),S.MinSin);
if isrow(el), zq=zq'; end
g=exp(sigma.*zq-.5*sigma.^2); gain(isGround)=g(isGround);
end
