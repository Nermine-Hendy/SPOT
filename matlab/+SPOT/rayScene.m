function [D,rate]=rayScene(C,G,k)
%RAYSCENE Per-pixel support and diffuse radiance from physical detector rays.
D=SPOT.detectorGrid(C,G.CameraToGCRF(:,:,k)); shape=size(D.Active);
r=G.Observer.Position_m(:,k); u=D.LOS_GCRF; RE=C.Constants.EarthRadius_m;
su=G.SunPosition_m(:,k); su=su/norm(su);
el=asind(max(-1,min(1,G.Zenith_GCRF(:,k)'*u)));
D.Elevation_deg=reshape(el,shape); D.EarthMask=false(shape); D.LimbDistance_deg=nan(shape);
disk=zeros(shape); rim=disk; sky=disk;
if G.IsGround(k)
    D.SkyMask=D.Active & D.Elevation_deg>=C.Visibility.MinElevation_deg;
    D.GeometryRegime='ground_sky'; D.ConservativeRegime='ground_sky';
    if C.Background.Enabled
        sunEl=asind(dot(su,G.Zenith_GCRF(:,k)));
        twilight=C.Background.TwilightSolarElevation_deg;
        f=max(0,min(1,(sunEl-twilight(1))/diff(twilight))); f=f*f*(3-2*f);
        base=C.Background.DarkRadiance+f*(C.Background.SunlitRadiance-C.Background.DarkRadiance);
        sky(:)=base;
        if C.Background.LightPollutionEnabled
            sky=sky+C.Background.LightPollutionRadiance*(1+C.Background.HorizonBoost*(1-sind(max(D.Elevation_deg,0))));
        end
    end
else
    dotru=r'*u; disc=dotru.^2-(dot(r,r)-RE*RE);
    near=-dotru-sqrt(max(disc,0)); far=-dotru+sqrt(max(disc,0)); near(near<=0)=far(near<=0);
    hit=disc>=0 & near>0; D.EarthMask=reshape(hit,shape)&D.Active;
    D.SkyMask=D.Active & ~D.EarthMask;
    alpha=asind(min(1,RE/norm(r))); nadir=-r/norm(r);
    beta=acosd(max(-1,min(1,nadir'*u))); D.LimbDistance_deg=reshape(beta-alpha,shape);
    if all(D.EarthMask(D.Active)), D.GeometryRegime='space_earth';
    elseif ~any(D.EarthMask,'all'), D.GeometryRegime='space_sky'; else, D.GeometryRegime='space_limb'; end
    b=G.CameraToGCRF(:,3,k); beta0=acosd(max(-1,min(1,dot(nadir,b))));
    if strcmp(C.Camera.Shape,'circular'), cone=C.Camera.FOV_deg(1)/2;
    else, cone=atand(hypot(tand(C.Camera.FOV_deg(1)/2),tand(C.Camera.FOV_deg(2)/2))); end
    if beta0+cone<alpha, D.ConservativeRegime='space_earth';
    elseif beta0-cone>alpha, D.ConservativeRegime='space_sky'; else, D.ConservativeRegime='space_limb'; end
    if C.Background.Enabled, sky(:)=C.Background.DarkRadiance; end
    if C.Earth.Enabled
        p=r+u.*near; normal=p/RE; mu=max(0,su'*normal); view=max(0,-sum(normal.*u,1));
        value=(C.Earth.NightRadiance+C.Earth.DayRadiance*mu.^C.Earth.LambertGamma).*(C.Earth.LimbMinWeight+(1-C.Earth.LimbMinWeight)*sqrt(view));
        value(~hit)=0; disk=reshape(value,shape).*D.Active;
        p=r+u.*max(0,-dotru); p=p./max(vecnorm(p),realmin); mu=max(0,su'*p);
        rim=C.Earth.RimGain*C.Earth.DayRadiance*exp(-.5*(D.LimbDistance_deg/C.Earth.RimSigma_deg).^2).*reshape(sqrt(mu),shape).*(abs(D.LimbDistance_deg)<=C.Earth.RimHalfWidth_deg).*D.Active;
    end
end
sunDirection=G.SunPosition_m(:,k)-r; sunDirection=sunDirection/norm(sunDirection);
D.SunExcluded=reshape(acosd(max(-1,min(1,sunDirection'*u)))<C.Visibility.SunExclusion_deg,shape)&D.Active;
D.SkyMask=D.SkyMask & ~D.SunExcluded;
D.SkyRadiance=sky.*D.SkyMask;
D.EarthRadiance=(disk+rim).*~D.SunExcluded;
D.DiffuseRadiance=D.SkyRadiance+D.EarthRadiance;
rate=D.DiffuseRadiance.*D.SolidAngle_sr;
D.EarthFraction=nnz(D.EarthMask)/nnz(D.Active);
end
