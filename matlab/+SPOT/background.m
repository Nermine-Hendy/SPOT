function B=background(C,startUTC,stars)
%BACKGROUND Exposure-average star and diffuse rates with their own time grid.
dur=C.Time.Exposure_s; mid=startUTC+seconds(dur/2);
extra=[];
if ~isempty(C.Pointing.JitterUTC), extra=seconds(C.Pointing.JitterUTC(:)-startUTC); end
if ~isempty(C.Pointing.AttitudeUTC), extra=[extra;seconds(C.Pointing.AttitudeUTC(:)-startUTC)]; end
[t,w]=SPOT.timeGrid(startUTC,dur,C.Integration.Step_s,extra);
for iteration=1:15
    assert(numel(t)<=C.Integration.MaxBackgroundSamples,'SPOT:Sampling','Star/background motion exceeds MaxBackgroundSamples. Increase the limit or shorten the exposure.');
    G=SPOT.sample(C,t,mid); step=zeros(1,numel(t)-1);
    scale=max(([C.Camera.Nx C.Camera.Ny]-1)./C.Camera.FOV_deg)/cosd(max(C.Camera.FOV_deg)/2)^2;
    for k=1:numel(step)
        Q=G.CameraToGCRF(:,:,k)'*G.CameraToGCRF(:,:,k+1);
        a=acosd(max(-1,min(1,(trace(Q)-1)/2)));
        % Observer motion also moves the Earth limb, even with fixed inertial attitude.
        r1=G.Observer.Position_m(:,k); r2=G.Observer.Position_m(:,k+1);
        h=max(1000,min(norm(r1),norm(r2))-C.Constants.EarthRadius_m);
        earthMotion=0; if ~all(G.IsGround(k:k+1)), earthMotion=rad2deg(norm(r2-r1)/h); end
        step(k)=(a+earthMotion)*scale;
    end
    if ~C.Integration.Adaptive || all(step<=C.Integration.MaxStep_pix*1.001), break; end
    q=seconds(t-startUTC)'; add=[];
    for k=find(step>C.Integration.MaxStep_pix)
        n=max(2,ceil(step(k)/C.Integration.MaxStep_pix));
        assert(n+numel(t)<=C.Integration.MaxBackgroundSamples,'SPOT:Sampling','Background sampling limit exceeded.');
        a=linspace(q(k),q(k+1),n+1); add=[add a(2:end-1)]; %#ok<AGROW>
    end
    [t,w]=SPOT.timeGrid(startUTC,dur,dur,[q add]);
end
if iteration==15 && any(step>C.Integration.MaxStep_pix*1.001) && C.Integration.Adaptive, error('SPOT:Sampling','Background grid failed to converge.'); end
% Do not allocate time histories for an entire all-sky catalogue. Retain the
% union of conservative FoV cones over ALL samples, including moving cameras.
keep=SPOT.starSupport(C,stars.Direction_GCRF,G.CameraToGCRF);
stars.Direction_GCRF=stars.Direction_GCRF(:,keep);
stars.Magnitude=stars.Magnitude(keep); stars.ID=stars.ID(keep);
stars.Metadata.WorkingSetCount=nnz(keep);
shape=[C.Camera.Ny C.Camera.Nx]; B.StarRate=zeros(shape); B.DiffuseRate=zeros(shape); B.SkyRate=zeros(shape); B.EarthRate=zeros(shape);
B.StarTruth=struct('x',[],'y',[],'Xi_deg',[],'Eta_deg',[],'RA_deg',[],'DEC_deg',[],'Magnitude',[]);
ns=numel(stars.Magnitude); pos=cell(1,numel(t)); rate=zeros(numel(t),ns); elevations=rate; ground=false(size(rate));
for k=1:numel(t)
    [D,diffuse]=SPOT.rayScene(C,G,k); B.DiffuseRate=B.DiffuseRate+w(k)*diffuse;
    B.SkyRate=B.SkyRate+w(k)*D.SkyRadiance.*D.SolidAngle_sr;
    B.EarthRate=B.EarthRate+w(k)*D.EarthRadiance.*D.SolidAngle_sr;
    if abs(seconds(t(k)-mid))<1e-7, B.Truth=D; end
    if ns==0, continue; end
    p=SPOT.project(C,G.CameraToGCRF(:,:,k),stars.Direction_GCRF); pos{k}=p;
    el=asind(max(-1,min(1,G.Zenith_GCRF(:,k)'*stars.Direction_GCRF))); elevations(k,:)=el; ground(k,:)=G.IsGround(k);
    visible=p.InFOV;
    if G.IsGround(k), visible=visible & el>=C.Visibility.MinElevation_deg;
    elseif C.Visibility.EarthOccultation
        r=G.Observer.Position_m(:,k); d=r'*stars.Direction_GCRF;
        visible=visible & ~(d<0 & d.^2>=dot(r,r)-C.Constants.EarthRadius_m^2);
    end
    su=G.SunPosition_m(:,k)-G.Observer.Position_m(:,k); su=su/norm(su);
    visible=visible & acosd(max(-1,min(1,su'*stars.Direction_GCRF)))>=C.Visibility.SunExclusion_deg;
    rate(k,:)=10.^(-.4*stars.Magnitude').*visible;
    if abs(seconds(t(k)-mid))<1e-7
        u=stars.Direction_GCRF(:,visible);
        B.StarTruth=struct('x',p.x(visible),'y',p.y(visible),'Xi_deg',p.Xi_deg(visible), ...
            'Eta_deg',p.Eta_deg(visible),'RA_deg',mod(atan2d(u(2,:),u(1,:)),360), ...
            'DEC_deg',asind(max(-1,min(1,u(3,:)))),'Magnitude',stars.Magnitude(visible));
    end
end
% Evaluate only stars intersecting this exposure; each has a stable ID stream.
for j=find(any(rate>0,1))
    key=hash(stars.ID(j)); [at,sc]=SPOT.atmosphere(C,t,elevations(:,j),ground(:,j),key);
    rate(:,j)=rate(:,j).*at.*sc;
end
for k=1:numel(t)
    if ns>0, B.StarRate=B.StarRate+w(k)*SPOT.deposit(C,pos{k}.x,pos{k}.y,rate(k,:),B.Truth.Active); end
end
B.TimeUTC=t; B.Weights=w; B.StarsUsed=nnz(any(rate>0,1));
B.StarCatalogueMetadata=stars.Metadata;
end
function key=hash(id)
key=7919; s=double(char(id));
for k=1:numel(s), key=mod(key*31+s(k),2^32); end
key=mod(key+100000,2^32);
end
