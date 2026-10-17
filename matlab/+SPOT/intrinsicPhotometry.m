function H=intrinsicPhotometry(C,G,j)
%INTRINSICPHOTOMETRY Deterministic target reflection, before atmospheric gains.
H=G.Target{j}; T=C.Targets(j); P=T.Photometry;
H.LambertPhase=ones(size(H.Range_m));
if strcmp(P.Model,'surface_brdf')
    sun=G.SunPosition_m-H.Position_m; ds=vecnorm(sun); sun=sun./ds;
    solar=C.Radiometry.SolarIrradianceAtAU_W_m2*(C.Constants.AU_m./ds).^2;
    R=SPOT.targetAttitude(T.Attitude,G.TimeUTC);
    H.Reflection=SPOT.surfaceReflection(T.Surface,R,sun,-H.LOS_GCRF,solar,max(H.Range_m,realmin));
    H.BodyToGCRF=R;
    H.UnattenuatedRate=H.Reflection.Irradiance_W_m2/C.Radiometry.ZeroMagnitudeIrradiance_W_m2;
else
    % Explicit legacy comparison path; area/albedo do not multiply magnitude.
    a=pi-H.PhaseAngle_rad;
    H.LambertPhase=max(0,(sin(a)-a.*cos(a))/pi);
    small=a<1e-3; H.LambertPhase(small)=(a(small).^3/3-a(small).^5/30+a(small).^7/840)/pi;
    if strcmp(P.Model,'constant_magnitude'), H.LambertPhase=ones(size(H.Range_m)); end
    H.UnattenuatedRate=10^(-.4*P.ReferenceMagnitude)*(P.ReferenceRange_m./max(H.Range_m,realmin)).^2;
end
H.LightCurveGain=ones(size(H.Range_m));
if ~isempty(P.LightCurveUTC)
    x=seconds(P.LightCurveUTC(:)-P.LightCurveUTC(1)); q=seconds(G.TimeUTC-P.LightCurveUTC(1));
    assert(all(q>=x(1)-1e-7 & q<=x(end)+1e-7),'SPOT:Photometry','Light curve does not cover exposure.');
    H.LightCurveGain=reshape(interp1(x,P.LightCurveGain,q,'linear'),1,[]);
end
H.IntrinsicRate=H.UnattenuatedRate.*H.LambertPhase.*H.Illumination.*H.LightCurveGain;
end
