function C=validate(C)
%VALIDATE Fail before rendering on ambiguous units, frames or unsupported modes.
assert(strcmp(C.SchemaVersion,'7.0.0'),'SPOT:Config','Use SPOT_Config from this release.');
C=SPOT.resolvePools(C);
C=SPOT.applyExampleOverrides(C);
pos(C.Time.Exposure_s,'Exposure_s'); pos(C.Time.Cadence_s,'Cadence_s');
integer(C.Time.NumExposures,1,'NumExposures'); utc(C.Time.StartUTC,'StartUTC',true);
C.Time.StartUTC.TimeZone='UTC';
integer(C.Camera.Nx,2,'Nx'); integer(C.Camera.Ny,2,'Ny');
assert(numel(C.Camera.FOV_deg)==2 && all(C.Camera.FOV_deg>0 & C.Camera.FOV_deg<180),'SPOT:Config','FOV must be two angles in (0,180) degrees.');
C.Camera.FOV_deg=reshape(C.Camera.FOV_deg,1,2);
choice(C.Camera.Shape,{'rectangular','circular'},'Camera.Shape');
if strcmp(C.Camera.Shape,'circular'), assert(diff(C.Camera.FOV_deg)==0,'SPOT:Config','Circular FOV requires equal angular diameters.'); end
choice(C.Camera.SolidAngle,{'mean','per_pixel'},'SolidAngle');
choice(C.Frames.Model,{'iau2006','gmst_precession'},'Frames.Model');
assert(strcmp(C.Frames.CelestialFrame,'GCRF'),'SPOT:Frame','Canonical celestial frame must be GCRF.');
pos(C.Integration.Step_s,'Step_s'); pos(C.Integration.MaxStep_pix,'MaxStep_pix');
integer(C.Integration.MaxSamples,3,'MaxSamples'); integer(C.Integration.MaxBackgroundSamples,3,'MaxBackgroundSamples');
choice(C.Integration.AttitudeMode,{'dynamic','midpoint'},'AttitudeMode');
pos(C.Morphology.MinSmear_pix,'Morphology.MinSmear_pix'); pos(C.Morphology.PSFFWHMFactor,'Morphology.PSFFWHMFactor');
choice(C.Seeding.Propagation,{'linear','two_body'},'Seeding.Propagation');
choice(C.Pointing.Mode,{'az_el_roll','nadir_offset','inertial','target_tracking','attitude'},'Pointing.Mode');
assert(isscalar(C.Pointing.Roll_deg) && isfinite(C.Pointing.Roll_deg),'SPOT:Pointing','Roll must be finite.');
if strcmp(C.Pointing.Mode,'az_el_roll')
    a=C.Pointing.AzEl_deg; assert(numel(a)==2 && all(isfinite(a)) && abs(a(2))<=90,'SPOT:Pointing','AzEl requires finite azimuth and elevation in [-90,90] degrees.');
elseif strcmp(C.Pointing.Mode,'nadir_offset')
    assert(isscalar(C.Pointing.OffNadir_deg) && isfinite(C.Pointing.OffNadir_deg) && C.Pointing.OffNadir_deg>=0 && C.Pointing.OffNadir_deg<=180 && isfinite(C.Pointing.SpaceAzimuth_deg),'SPOT:Pointing','Invalid off-nadir or orbital azimuth.');
end
if strcmp(C.Pointing.Mode,'attitude')
    series(C.Pointing.AttitudeUTC,C.Pointing.CameraToGCRF,3,'Attitude');
    assert(size(C.Pointing.CameraToGCRF,1)==3 && size(C.Pointing.CameraToGCRF,2)==3,'SPOT:Config','CameraToGCRF must be 3-by-3-by-N.');
    for k=1:size(C.Pointing.CameraToGCRF,3)
        A=C.Pointing.CameraToGCRF(:,:,k);
        assert(norm(A'*A-eye(3),'fro')<1e-6 && det(A)>0,'SPOT:Attitude','Attitude matrices must be proper orthonormal rotations.');
        [U,~,V]=svd(A); C.Pointing.CameraToGCRF(:,:,k)=U*V'; % remove accepted roundoff
    end
end
if ~isempty(C.Pointing.JitterUTC)
    series(C.Pointing.JitterUTC,C.Pointing.Jitter_arcsec,1,'Jitter');
    assert(size(C.Pointing.Jitter_arcsec,2)==3,'SPOT:Config','Jitter is N-by-3 [column,row,roll] arcseconds.');
end
pos(C.Integration.PhotometryRelTol,'PhotometryRelTol'); pos(C.Integration.PhotometryAbsTol,'PhotometryAbsTol');
pos(C.Integration.MaxBodyRotation_deg,'MaxBodyRotation_deg'); pos(C.Integration.SpecularSamplesPerHWHM,'SpecularSamplesPerHWHM');
for field={'SolarIrradianceAtAU_W_m2','ZeroMagnitudeIrradiance_W_m2','EffectiveWavelength_m','ApertureDiameter_m'}
    pos(C.Radiometry.(field{1}),field{1});
end
unit(C.Radiometry.QuantumEfficiency,'QuantumEfficiency');
choice(C.Radiometry.DetectorCalibration,{'band_photon','legacy_gain'},'DetectorCalibration');
if strcmp(C.Radiometry.DetectorCalibration,'band_photon')
    photonEnergy=6.62607015e-34*299792458/C.Radiometry.EffectiveWavelength_m;
    C.Detector.RadiometricGain_e_per_unit=C.Radiometry.ZeroMagnitudeIrradiance_W_m2* ...
        pi*(C.Radiometry.ApertureDiameter_m/2)^2*C.Radiometry.QuantumEfficiency/photonEnergy;
end
ids=[C.Targets.ID]; assert(all(ids>=1 & ids==fix(ids) & ids<=double(intmax('uint32'))) && numel(unique(ids))==numel(ids),'SPOT:Config','Target IDs must be unique positive uint32 integers.');
if strcmp(C.Pointing.Mode,'target_tracking')
    j=find(ids==C.Pointing.TrackTargetID,1); assert(~isempty(j),'SPOT:Pointing','Tracked target ID is absent.');
    assert(~strcmpi(C.Targets(j).Trajectory.Type,'demo'),'SPOT:Pointing','Camera-relative demo trajectories cannot define tracking attitude.');
end
choice(C.Visibility.EclipseModel,{'cylinder','finite_sun'},'EclipseModel');
assert(numel(C.Visibility.RangeLimits_m)==2 && C.Visibility.RangeLimits_m(1)>0 && diff(C.Visibility.RangeLimits_m)>0,'SPOT:Config','Range limits must be increasing and positive.');
assert(C.Visibility.MinElevation_deg>=0 && C.Visibility.MinElevation_deg<90,'SPOT:Config','MinElevation must lie in [0,90).');
choice(C.Atmosphere.AirmassModel,{'secant','kasten_young'},'AirmassModel');
pos(C.Atmosphere.Scintillation.CorrelationTime_s,'Scintillation correlation time');
pos(C.Atmosphere.Scintillation.KnotStep_s,'Scintillation knot step');
assert(C.Atmosphere.Scintillation.Sigma0>=0 && C.Atmosphere.Scintillation.MinSin>0 && C.Atmosphere.Scintillation.MinSin<=1,'SPOT:Config','Invalid scintillation parameters.');
unit(C.Optics.Transmission,'Transmission'); unit(C.Optics.EdgeTransmission,'EdgeTransmission');
assert(all(C.Optics.ElementTransmission>=0 & C.Optics.ElementTransmission<=1),'SPOT:Config','Element transmissions must lie in [0,1].');
choice(C.Optics.PSFModel,{'field_dependent_gaussian','gaussian'},'PSFModel');
assert(C.Optics.SigmaCenter_pix>=0 && C.Optics.SigmaEdgeAdd_pix>=0,'SPOT:Config','PSF widths must be nonnegative.');
integer(C.Optics.PSFBanks,2,'PSFBanks'); pos(C.Optics.RadialExponent,'RadialExponent');
pos(C.Masks.RelativeThreshold,'Masks.RelativeThreshold'); unit(C.Masks.RelativeThreshold,'Masks.RelativeThreshold');
pos(C.Detector.RadiometricGain_e_per_unit,'RadiometricGain'); pos(C.Detector.ConversionGain_e_per_DN,'ConversionGain');
pos(C.Detector.FullWell_e,'FullWell'); integer(C.Detector.ADCBits,1,'ADCBits');
assert(C.Detector.ADCBits<=32 && C.Detector.DarkCurrent_e_s>=0 && C.Detector.ReadNoise_e>=0,'SPOT:Config','Invalid detector parameters.');
for name={'PRNUMap','BadPixelMask'}
    a=C.Detector.(name{1}); assert(isempty(a) || isequal(size(a),[C.Camera.Ny C.Camera.Nx]),'SPOT:Config','Detector maps must match Ny-by-Nx.');
    assert(all(isfinite(a(:))) && all(a(:)>=0),'SPOT:Config','Detector maps must be finite and nonnegative.');
end
choice(C.Stars.Source,{'csv','synthetic'},'Stars.Source'); assert(strcmp(C.Stars.Frame,'GCRF'),'SPOT:Frame','Catalogue must use GCRF/ICRS-aligned RA and DEC.');
integer(C.Stars.SyntheticCount,0,'SyntheticCount'); integer(C.Stars.MaxStars,1,'MaxStars');
assert(numel(C.Stars.SyntheticMagnitudeRange)==2 && all(isfinite(C.Stars.SyntheticMagnitudeRange)) && diff(C.Stars.SyntheticMagnitudeRange)>=0,'SPOT:Config','Invalid synthetic magnitude range.');
nonnegative=[C.Background.DarkRadiance C.Background.SunlitRadiance C.Background.LightPollutionRadiance C.Background.HorizonBoost C.Earth.DayRadiance C.Earth.NightRadiance C.Earth.RimGain C.Atmosphere.Extinction_mag_airmass];
twilight=C.Background.TwilightSolarElevation_deg;
assert(numel(twilight)==2 && all(isfinite(twilight)) && diff(twilight)>0 && all(abs(twilight)<=90),'SPOT:Config','Twilight solar-elevation bounds must increase within [-90,90] degrees.');
assert(all(isfinite(nonnegative)) && all(nonnegative>=0),'SPOT:Config','Radiance and extinction settings must be finite and nonnegative.');
pos(C.Earth.RimSigma_deg,'RimSigma_deg'); pos(C.Earth.RimHalfWidth_deg,'RimHalfWidth_deg'); pos(C.Earth.LambertGamma,'LambertGamma'); unit(C.Earth.LimbMinWeight,'LimbMinWeight');
pos(C.Camera.DiagnosticGridStep_deg,'DiagnosticGridStep_deg');
assert(C.Visibility.SunExclusion_deg>=0 && C.Visibility.SunExclusion_deg<180,'SPOT:Config','Sun exclusion must lie in [0,180) degrees.');
assert(C.Atmosphere.MinElevation_deg>0 && C.Atmosphere.MinElevation_deg<=90,'SPOT:Config','Atmospheric elevation floor must lie in (0,90].');
assert(numel(C.Output.PNGPercentiles)==2 && all(C.Output.PNGPercentiles>=0 & C.Output.PNGPercentiles<=100) && diff(C.Output.PNGPercentiles)>0,'SPOT:Config','PNG percentiles must increase within [0,100].');
integer(C.RandomSeed,0,'RandomSeed'); assert(C.RandomSeed<2^32,'SPOT:Config','RandomSeed must be below 2^32.');
C.Observer.Trajectory=traj(C,C.Observer.Trajectory,false);
C=SPOT.seedTargets(C);
for j=1:numel(C.Targets)
    assert(~isfield(C.Targets,'Rendering'),'SPOT:Config','Remove Target.Rendering: morphology is derived from the integrated motion.');
    choice(C.Targets(j).Photometry.Model,{'surface_brdf','lambert_sphere','constant_magnitude'},'Photometry.Model');
    pos(C.Targets(j).Photometry.ReferenceRange_m,'ReferenceRange_m');
    assert(isfinite(C.Targets(j).Photometry.ReferenceMagnitude),'SPOT:Config','ReferenceMagnitude must be finite.');
    p=C.Targets(j).Photometry;
    if strcmp(p.Model,'surface_brdf'), SPOT.validateSurface(C.Targets(j)); end
    if ~isempty(p.LightCurveUTC)
        series(p.LightCurveUTC,p.LightCurveGain,1,'LightCurve');
        assert(all(p.LightCurveGain>=0),'SPOT:Config','Light-curve gain must be nonnegative.');
    end
    C.Targets(j).Trajectory=traj(C,C.Targets(j).Trajectory,true);
end
end
function T=traj(C,T,isTarget)
T.Type=lower(char(T.Type)); T.Frame=upper(char(T.Frame));
choice(T.Type,{'geodetic','ecef','state','ephemeris','keplerian','tle','demo'},'Trajectory.Type');
utc(T.EpochUTC,'Trajectory.EpochUTC',true);
T.EpochUTC.TimeZone='UTC';
choice(T.Frame,{'GCRF','ECEF'},'Trajectory.Frame');
choice(T.Motion,{'fixed','linear','two_body'},'Trajectory.Motion');
choice(T.Interpolation,{'hermite','linear'},'Trajectory.Interpolation');
choice(T.Extrapolation,{'error'},'Trajectory.Extrapolation');
switch T.Type
    case 'geodetic'
        a=T.Geodetic_deg_m; assert(numel(a)==3 && all(isfinite(a)) && abs(a(1))<=90 && abs(a(2))<=180,'SPOT:Config','Geodetic state is [latitude_deg longitude_deg ellipsoid_height_m].');
        assert(~strcmp(T.Motion,'two_body'),'SPOT:Config','Use a GCRF state for two-body propagation.');
    case {'ecef','state'}
        assert(numel(T.Position_m)==3 && norm(T.Position_m)>0 && all(isfinite(T.Position_m)) && numel(T.Velocity_mps)==3 && all(isfinite(T.Velocity_mps)),'SPOT:Config','State needs finite position and velocity in metres and m/s.');
        if strcmp(T.Type,'ecef'), T.Frame='ECEF'; end
        assert(~strcmp(T.Motion,'two_body') || strcmp(T.Frame,'GCRF'),'SPOT:Config','two_body requires GCRF.');
    case 'keplerian'
        a=T.Elements; assert(numel(a)==6 && all(isfinite(a)) && a(1)>0 && a(2)>=0 && a(2)<1 && a(3)>=0 && a(3)<=180,'SPOT:Config','Elements are [a_m e i_deg RAAN_deg argument_perigee_deg true_anomaly_deg], bound elliptic orbits only.');
        T.Frame='GCRF';
    case 'tle'
        integer(T.TLEIndex,1,'TLEIndex'); pos(T.MaxTLEAge_days,'MaxTLEAge_days');
    case 'demo'
        assert(isTarget && C.Dataset.AllowDemoTrajectories,'SPOT:Demo','Demo trajectories are allowed only for explicitly enabled artificial targets.');
        choice(T.Demo.AngleMode,{'az_el_offset','theta_phi','camera_angles'},'Demo.AngleMode');
        choice(T.Demo.MotionModel,{'linear_angles','great_circle'},'Demo.MotionModel');
        pos(T.Demo.Range_m,'Demo.Range_m'); utc(T.Demo.ReferenceUTC,'Demo.ReferenceUTC',true);
        assert(numel(T.Demo.Offset_deg)==2 && numel(T.Demo.AngularRate_deg_s)==2 && all(isfinite([T.Demo.Offset_deg T.Demo.AngularRate_deg_s])),'SPOT:Demo','Invalid demo angles/rates.');
        if strcmp(T.Demo.AngleMode,'az_el_offset'), assert(strcmp(C.Pointing.Mode,'az_el_roll'),'SPOT:Demo','Az/el demo offsets require az_el_roll pointing.'); end
end
T=SPOT.prepareTrajectory(C,T);
end
function utc(t,label,scalar)
assert(isdatetime(t) && ~isempty(t.TimeZone) && all(~isnat(t)) && (~scalar || isscalar(t)),'SPOT:Time','%s must be timezone-aware datetime.',label);
end
function series(t,a,dim,label)
utc(t,label,false); assert(numel(t)>=2 && all(seconds(diff(t(:)))>0) && size(a,dim)==numel(t) && all(isfinite(a(:))),'SPOT:Config','%s must be finite and strictly time ordered, with matching samples.',label);
end
function choice(a,b,label)
assert(any(strcmp(a,b)),'SPOT:Config','Unsupported %s: %s.',label,a);
end
function pos(a,label)
assert(isscalar(a) && isfinite(a) && a>0,'SPOT:Config','%s must be positive.',label);
end
function unit(a,label)
assert(isscalar(a) && isfinite(a) && a>=0 && a<=1,'SPOT:Config','%s must lie in [0,1].',label);
end
function integer(a,lo,label)
assert(isscalar(a) && isfinite(a) && a==fix(a) && a>=lo,'SPOT:Config','%s must be an integer >= %g.',label,lo);
end
