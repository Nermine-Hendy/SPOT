function C = SPOT_Config(preset,quality)
%SPOT_CONFIG Shared defaults and presets for scenes, targets and plots.
% RUN_SPOT                         runs SelectedExample below and saves plots.
% SPOT_Config('list')              lists every example.
% C=SPOT_Config('figure8'); Out=RUN_SPOT(C); % edit C before running if desired.
%% 01 | QUICK START: select an example and acquisition quality
% Edit these two settings for RUN_SPOT with no arguments.
% Explicit SPOT_Config(preset,quality) arguments take precedence.
SelectedExample='ground'; % ground | space_sky | space_limb | space_earth | figure8
SelectedQuality='full';   % full: 1024px, catalogue, IAU2006, Poisson noise
                         % preview: 256px, synthetic stars, approximate frames
%% 02 | ARGUMENT HANDLING AND PRESET NAMES
% Internal selection logic. Keep this block unless adding an example.
% The current_* profiles use linear seed propagation and midpoint attitude sampling.
if nargin<1 || isempty(preset), preset=SelectedExample; end
if nargin<2 || isempty(quality), quality=SelectedQuality; end
if islogical(quality), if quality, quality='preview'; else, quality='full'; end, end
quality=validatestring(char(quality),{'full','preview'}); preset=char(preset);
requestedPreset=preset;
userProfile=startsWith(preset,'current_');
%% 30 | CURRENT_* ACQUISITION PROFILES
% Applied late: replaces epoch, exposure/cadence, catalogue cap, propagation,
% camera sampling, detector calibration, face area and seed velocities.
% For these profiles, edit this block or edit the returned C structure.
% The retained 1 s exposure and 0.8 s cadence overlap by 0.2 s.
if userProfile, preset=extractAfter(preset,'current_'); end
aliases={'paper_ground','ground';'paper_space_sky','space_sky'; ...
    'paper_space_limb','space_limb';'paper_space_earth','space_earth';'paper_mixed','figure8'};
k=find(strcmp(preset,aliases(:,1)),1); if ~isempty(k), preset=aliases{k,2}; end
names={'ground','space_sky','space_limb','space_earth','figure8','figure8_short','figure8_long', ...
    'ground_mixed','space_limb_mixed','space_earth_mixed', ...
    'multi_exposure','orbit_tracking','trajectory_geodetic','trajectory_ecef','trajectory_state', ...
    'trajectory_ephemeris','trajectory_ephemeris_ecef','trajectory_ephemeris_csv', ...
    'trajectory_keplerian','trajectory_tle','trajectory_tle_deep_space'};
if strcmp(preset,'list'), C=[names {'current_ground','current_space_sky','current_space_limb','current_space_earth'}]; return; end
assert(any(strcmp(preset,names)),'SPOT:Config','Unknown example %s. Use SPOT_Config(''list'').',preset);
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root); SPOT_Setup;
%% 03 | RUN IDENTITY AND REPRODUCIBILITY
% General metadata and the master random seed.
C.SchemaVersion='7.0.0';
C.SoftwareVersion='7.2.0'; 
C.Name=char(preset); 
C.RootDir=root; 
% Master seed for repeatable stochastic runs with the same configuration/backend.
C.RandomSeed=1;
%% 04 | EXPOSURE SCHEDULE
% UTC start of the first exposure, integration duration, start-to-start
% cadence and exposure count. Cadence smaller than exposure means overlap.
C.Time.StartUTC=datetime(2026,1,1,12,0,0,'TimeZone','UTC');
C.Time.Exposure_s=1.5; 
C.Time.Cadence_s=1.5; 
C.Time.NumExposures=1;
%% 05 | TEMPORAL INTEGRATION AND ADAPTIVE SAMPLING
% Controls numerical integration within each exposure. Smaller steps and
% tighter tolerances improve resolution at greater computation cost.
% AttitudeMode here refers to the CAMERA, not the target body attitude.
% Initial maximum time step [s]; 0.8/0.02 gives 40 intervals, 41 endpoints.
C.Integration.Step_s=0.02; % Initial grid; adaptive refinement may add more samples.
C.Integration.Adaptive=true; 
% Maximum requested image-plane movement between temporal samples [pixels].
C.Integration.MaxStep_pix=0.35;
% Upper target-sampling budget; MaxBackgroundSamples separately caps background sampling.
C.Integration.MaxSamples=8193; 
C.Integration.MaxBackgroundSamples=257;
C.Integration.AttitudeMode='dynamic'; % 'midpoint' is the paper's approximation.
% Sampling tolerances apply to intrinsic reflection, before scintillation.
C.Integration.PhotometryAdaptive=true;
% Relative tolerance for adaptive intrinsic brightness changes [dimensionless].
C.Integration.PhotometryRelTol=0.01;
C.Integration.PhotometryAbsTol=1e-12; % normalized zero-magnitude irradiance units
% Angular sampling safeguard for target body motion [degrees].
C.Integration.MaxBodyRotation_deg=1;
% Sampling density per specular-lobe half width at half maximum (HWHM).
C.Integration.SpecularSamplesPerHWHM=8;
% Target irradiance / ZeroMagnitudeIrradiance gives scene rate units.
% These are explicit effective-band assumptions, not a measured instrument.
%% 06 | EFFECTIVE-BAND RADIOMETRY AND ELECTRON CALIBRATION
% Surface irradiance is divided by ZeroMagnitudeIrradiance_W_m2 to obtain
% the scene-rate normalization. These effective-band assumptions should
% be replaced by calibrated instrument values for quantitative validation.
C.Radiometry.Band='effective visible band';
% Effective-band incident solar irradiance at one astronomical unit [W/m^2].
C.Radiometry.SolarIrradianceAtAU_W_m2=500;
% Reference irradiance for magnitude zero in the adopted effective band [W/m^2].
C.Radiometry.ZeroMagnitudeIrradiance_W_m2=3.6e-8;
% Representative wavelength used to convert energy to photon count [m].
C.Radiometry.EffectiveWavelength_m=550e-9;
% Detected-electron probability per incident photon [0 to 1].
C.Radiometry.QuantumEfficiency=0.65;
% Collecting aperture diameter [m]; sets the collecting area.
C.Radiometry.ApertureDiameter_m=0.05;
% band_photon derives conversion from aperture/wavelength/QE; legacy_gain uses the explicit gain.
C.Radiometry.DetectorCalibration='band_photon'; % or explicit 'legacy_gain'
%% 07 | OPTIONAL SENSOR AND TARGET CSV POOLS
% false uses the configured observer/targets below; true requests CSV input
% resolution. NumTargets controls selection from the target pool.
% Selection modes: ordered, random, ids. StartRow is one-based; Seed controls
% repeatable random selection; IDs selects named rows. Empty filters disable
% filtering. Sensor and target pools have independent selection settings.
C.Input.UseSensorCSV=false; 
C.Input.UseTargetCSV=false;
C.Input.SensorCSV=fullfile(root,'Inputs','SensorPool.csv');
C.Input.TargetCSV=fullfile(root,'Inputs','TargetPool_Seeded.csv');
C.Input.NumTargets=4;
selection=struct('Mode','ordered','StartRow',1,'Seed',202,'IDs',strings(0,1), ...
    'FilterColumn','','FilterValue','');
C.Input.SensorSelection=selection; 
C.Input.SensorSelection.Seed=101;
C.Input.SensorSelection.FilterColumn='profile'; 
C.Input.SensorSelection.FilterValue='ground';
C.Input.TargetSelection=selection;
%% 08 | APPEARANCE CLASSIFICATION AND SEED PROPAGATION
% All targets use exposure-integrated rendering. These thresholds label
% the resulting appearance; they do not force an object or trail renderer.
% Compact threshold = max(MinSmear_pix, PSFFWHMFactor * local PSF FWHM).
C.Morphology.MinSmear_pix=1; % Label only: rendering always integrates motion.
C.Morphology.PSFFWHMFactor=1; % Compact if equivalent smear <= max(1px, local FWHM).
% Propagation model after the initial seed is resolved: two_body or linear.
C.Seeding.Propagation='two_body'; % Initial seeds become Cartesian states once.
%% 09 | PHYSICAL CONSTANTS
% EarthRadius_m is the spherical radius; EarthA_m and EarthEccentricity
% define the geodetic ellipsoid. Mu is the Earth gravitational parameter.
C.Constants.EarthRadius_m=6371008.8; 
C.Constants.EarthA_m=6378137;
C.Constants.EarthEccentricity=8.1819190842622e-2;
C.Constants.Mu_m3_s2=3.986004418e14; 
C.Constants.AU_m=149597870700;
C.Constants.SunRadius_m=695700e3;
%% 10 | REFERENCE FRAMES AND EARTH ORIENTATION
% GCRF is the adopted inertial frame. DeltaAT = TAI - UTC; DeltaUT1 =
% UT1 - UTC. PolarMotion_rad contains [xp yp] in radians. Use epoch-appropriate
% Earth-orientation inputs when accurate absolute pointing is required.
C.Frames.Model='iau2006'; % Aerospace Toolbox; explicit 'gmst_precession' alternative.
C.Frames.DeltaAT_s=37; 
C.Frames.DeltaUT1_s=0; 
C.Frames.PolarMotion_rad=[0 0];
C.Frames.CelestialFrame='GCRF';
%% 11 | VISIBILITY, OCCULTATION AND ILLUMINATION SUPPORT
% Geometric eligibility and eclipse settings. Atmospheric airmass has its
% own positive minimum elevation in Section 17; it is a separate setting.
% Altitude threshold used for ground/space platform classification [m].
C.Visibility.PlatformAltitude_m=100e3; 
% Ground-horizon visibility cutoff [deg], separate from the airmass floor.
C.Visibility.MinElevation_deg=0;
% Allowed observer-to-target range [minimum maximum] in metres.
C.Visibility.RangeLimits_m=[1 Inf]; 
C.Visibility.EarthOccultation=true;
C.Visibility.SunExclusion_deg=0; % disabled in paper; when >0, mark excluded.
% cylinder is the simple shadow model; finite_sun allows partial illumination.
C.Visibility.EclipseModel='cylinder'; % also 'finite_sun' for partial eclipses.
%% 12 | DETECTOR SAMPLING AND FIELD OF VIEW
% Nx = columns; Ny = rows. FOV_deg = [horizontal vertical] full angular
% extents in degrees. DiagnosticGridStep_deg controls angular diagnostics.
C.Camera.Nx=1024; C.Camera.Ny=1024; 
C.Camera.FOV_deg=[20 20];
C.Camera.Shape='rectangular'; 
% mean uses one average pixel solid angle; per_pixel retains spatial variation.
C.Camera.SolidAngle='mean'; % Optional alternative: 'per_pixel'.
C.Camera.DiagnosticGridStep_deg=0.25;
%% 13 | SENSOR POINTING, CAMERA ATTITUDE AND JITTER
% Only fields used by the selected pointing mode are active. Ground
% pointing uses azimuth/elevation/roll; space examples set nadir_offset;
% target_tracking uses TrackTargetID. GCRF vectors describe directions.
% CameraToGCRF contains 3-by-3 rotation matrices at AttitudeUTC epochs.
% This camera orientation is distinct from each target body attitude.
C.Pointing.Mode='az_el_roll'; 
% Ground boresight [azimuth elevation] in degrees; Roll_deg rotates about the boresight.
C.Pointing.AzEl_deg=[180 60]; 
C.Pointing.Roll_deg=0;
% Space boresight offset from nadir [deg]; SpaceAzimuth_deg selects its orbital azimuth.
C.Pointing.OffNadir_deg=70; 
C.Pointing.SpaceAzimuth_deg=0;
C.Pointing.Boresight_GCRF=[1;0;0]; 
C.Pointing.ReferenceUp_GCRF=[0;0;1];
% Target identifier followed when Mode is target_tracking.
C.Pointing.TrackTargetID=1; 
C.Pointing.AttitudeUTC=datetime.empty(0,1);
C.Pointing.CameraToGCRF=zeros(3,3,0); % optional measured matrices, interpolated on SO(3).
% Optional jitter time history; Jitter_arcsec is N-by-3 [column row roll] in arcsec.
C.Pointing.JitterUTC=datetime.empty(0,1); 
C.Pointing.Jitter_arcsec=zeros(0,3);
%% 14 | OBSERVER TRAJECTORY AND COMMON STATE TEMPLATE
% The trajectory structure carries alternative input representations; the
% Type field selects which representation is read. The same template is
% reused for targets. Positions are metres and velocities metres/second.
% Geodetic_deg_m = [latitude_deg longitude_deg altitude_m].
% Elements = [a_m eccentricity inclination_deg RAAN_deg argumentOfPerigee_deg
%             trueAnomaly_deg]. EpochUTC is the reference state epoch.
% Sampled Positions_m/Velocities_mps are 3-by-N arrays at TimeUTC epochs.
% File identifies an external ephemeris/TLE when selected. TLEIndex selects
% a TLE record; MaxTLEAge_days bounds age. Demo is an explicit demonstration
% trajectory description, not the default physically propagated seed.
C.Observer.Name='Observer';
tr=struct('Type','geodetic','Frame','ECEF','EpochUTC',C.Time.StartUTC, ...
    'Geodetic_deg_m',[-37.8136 144.9631 50],'Position_m',[0;0;0], ...
    'Velocity_mps',[0;0;0],'VelocityENU_mps',[0;0;0],'Motion','fixed', ...
    'Elements',[6878137 0.001 97.5 20 0 0], ...
    'File','','TimeUTC',datetime.empty(0,1),'Positions_m',zeros(3,0), ...
    'Velocities_mps',zeros(3,0),'Interpolation','hermite','Extrapolation','error', ...
    'TLELines',{{}},'TLEIndex',1,'MaxTLEAge_days',14, ...
    'Demo',struct('AngleMode','az_el_offset','MotionModel','linear_angles','Offset_deg',[-6 5], ...
    'Range_m',900e3,'AngularRate_deg_s',[0 0],'ReferenceUTC',C.Time.StartUTC+seconds(.1)));
C.Observer.Trajectory=tr;
% RAAN=20 deg keeps seeded space targets sunlit at the selected epoch.
% The paper's RAAN=40 deg places this near-observer scene in Earth's shadow.
% Eclipse masking remains enabled; brightness is never forced for a figure.
%% 15 | TARGET PHOTOMETRY, MATERIALS, ATTITUDE AND SEEDS
% Target ID/Name/Class are identifiers and metadata. Trajectory defines
% motion; Seed supplies initial geometry/velocity for seeded_state targets.
% The active default Photometry.Model is surface_brdf. ReferenceMagnitude
% and ReferenceRange_m are inactive for that model. Optional light-curve
% gains are dimensionless samples at LightCurveUTC epochs.
ph=struct('Model','surface_brdf','ReferenceMagnitude',4, ...
    'ReferenceRange_m',1e6,'LightCurveUTC',datetime.empty(0,1), ...
    'LightCurveGain',zeros(0,1)); % ReferenceMagnitude/Range are used ONLY by explicit legacy models.
tgt=struct('ID',1,'Name','Target_1','Class','RSO', ...
    'Trajectory',tr,'Photometry',ph,'Metadata',struct(), ...
    'Seed',struct('AngleMode','az_el_offset','Angles_deg',[-6 5], ...
    'Range_m',900e3,'TransverseSpeed_mps',6000,'Direction_deg',-35,'RadialSpeed_mps',30));
% Convex box: six one-sided faces, each 0.1 m^2 (side ~0.316 m).
% Face normals are body-fixed, never aligned automatically with Sun/observer.
% Each column of NormalBody is one unit face normal in target body coordinates.
% Area_m2, DiffuseReflectance, SpecularReflectance and SpecularExponent
% have one entry per face. Reflectances are dimensionless, nonnegative,
% and their sum per face must not exceed one. Higher exponent narrows
% the specular lobe. Changing body attitude rotates the reflection direction.
tgt.Surface=struct('NormalBody',[eye(3) -eye(3)],'Area_m2',0.1*ones(1,6), ...
    'DiffuseReflectance',0.5*ones(1,6),'SpecularReflectance',0.1*ones(1,6), ...
    'SpecularExponent',40*ones(1,6));
% Target attitude modes: fixed, spin, and samples (sampled rotations). BodyToGCRF0 is
% the initial body-to-inertial rotation. SpinAxisBody is a body-frame unit
% axis; SpinRate_deg_s is degrees/second. Sampled BodyToGCRF is 3-by-3-by-N
% at TimeUTC epochs. These fields describe target orientation, not orbit.
tgt.Attitude=struct('Mode','spin','EpochUTC',C.Time.StartUTC, ...
    'BodyToGCRF0',eye(3),'SpinAxisBody',[0;0;1],'SpinRate_deg_s',5, ...
    'TimeUTC',datetime.empty(0,1),'BodyToGCRF',zeros(3,3,0));
C.Targets=repmat(tgt,1,4);
% Default velocities were designed for non-Earth-intersecting orbits.
% They remain editable physical inputs; no runtime adjustment forces a class.
% Ground seeds: angle offset [az,el], range [m], transverse speed [m/s],
% direction [deg from +xi toward +eta], radial speed [m/s]. No rendering type.
groundSeeds=[-6 5 900000 7578.659567 -35 1527.922877; ...
    -2 -4 1.1e+06 7720.538584 -10 630.500087; ...
    3 2 1.2e+06 7589.480648 20 -999.521913; ...
    7 -6 1.5e+06 7234.788812 35 -1998.735390];
for j=1:4
    C.Targets(j).ID=j; C.Targets(j).Name=sprintf('Target_%d',j);
    C.Targets(j).Trajectory.Type='seeded_state';
    C.Targets(j).Seed=SPOT_ConfigSeed('az_el_offset',groundSeeds(j,:));
end
%% 16 | STELLAR BACKGROUND
% Choose csv or synthetic sources. Catalogue selection is bounded by
% MaxMagnitude and MaxStars. Synthetic settings apply to synthetic sources.
% Gaia G magnitudes and the effective instrument band are not identical.
C.Stars.Enabled=true; 
C.Stars.Source='csv';
C.Stars.File=fullfile(root,'Inputs','Stars_648000.csv');
C.Stars.Frame='GCRF'; 
C.Stars.Catalogue='Gaia DR3 (supplied subset)'; C.Stars.Band='Gaia G';
C.Stars.EpochUTC=datetime(2016,1,1,'TimeZone','UTC');
C.Stars.MaxMagnitude=12; 
C.Stars.MaxStars=648000; % Entire supplied catalogue, subject to magnitude cut.
C.Stars.SyntheticCount=20000; 
C.Stars.SyntheticMagnitudeRange=[4 12];
%% 17 | ATMOSPHERIC EXTINCTION AND SCINTILLATION
% Ground-path attenuation and intensity fluctuations. Secant airmass
% requires a strictly positive elevation floor to avoid the horizon
% singularity. Scintillation is multiplicative and mean-preserving.
C.Atmosphere.Enabled=true; 
C.Atmosphere.AirmassModel='secant'; % or 'kasten_young'
% Extinction coefficient [magnitudes per airmass].
C.Atmosphere.Extinction_mag_airmass=.2; 
C.Atmosphere.MinElevation_deg=1;
C.Atmosphere.Scintillation.Enabled=true; 
% Baseline lognormal fluctuation strength [dimensionless].
C.Atmosphere.Scintillation.Sigma0=.08;
% Floor on elevation sine in the scintillation strength law [dimensionless].
C.Atmosphere.Scintillation.MinSin=.2;
% Characteristic correlation time of intensity fluctuations [s].
C.Atmosphere.Scintillation.CorrelationTime_s=.05;
C.Atmosphere.Scintillation.KnotStep_s=.01; % independent of quadrature grid.
%% 18 | DIFFUSE SKY BACKGROUND
% Phenomenological normalized radiance terms for dark sky, twilight,
% light pollution and horizon enhancement; not calibrated W/m^2/sr values.
C.Background.Enabled=true; 
C.Background.DarkRadiance=1e-7;
% Solar elevation limits [dark endpoint daylight endpoint] in degrees.
C.Background.TwilightSolarElevation_deg=[-18 0];
C.Background.SunlitRadiance=1e-4; C.Background.LightPollutionEnabled=true;
C.Background.LightPollutionRadiance=5e-7; 
C.Background.HorizonBoost=3;
%% 19 | EARTH DISK AND LIMB BACKGROUND
% Phenomenological Earth/rim radiance in the explicitly named convention.
% These normalized amplitudes are illustrative; physical radiance requires
% an explicit calibration and conversion model.
C.Earth.Enabled=true; 
C.Earth.DayRadiance=350; 
C.Earth.NightRadiance=2;
% Day/night amplitudes use the normalized convention recorded below.
C.Earth.RadianceConvention='original_figure_code_normalized';
% Exponent controlling Earth illumination dependence; LimbMinWeight is its limb weighting floor.
C.Earth.LambertGamma=1; 
C.Earth.LimbMinWeight=.25;
% Rim amplitude relative to DayRadiance; evaluated at run time.
C.Earth.RimGain=3; % Relative to DayRadiance.
% Angular width of the limb-rim profile [deg]; RimHalfWidth_deg bounds its support.
C.Earth.RimSigma_deg=1.6; 
C.Earth.RimHalfWidth_deg=8;
%% 20 | OPTICAL THROUGHPUT, VIGNETTING AND PSF
% Dimensionless transmission and spatial spreading on the detector grid.
% Gaussian widths are standard deviations in pixels, not FWHM.
% PSFBanks controls the discrete approximation to field-varying blur.
C.Optics.Enabled=true; 
C.Optics.Transmission=.85; 
% Optional per-element throughput factors; their product replaces Transmission.
% Empty leaves the scalar Transmission setting active.
C.Optics.ElementTransmission=[];
C.Optics.VignettingEnabled=true; 
% Relative edge transmission; RadialExponent controls radial falloff.
C.Optics.EdgeTransmission=.6; 
C.Optics.RadialExponent=2;
C.Optics.PSFEnabled=true; 
C.Optics.PSFModel='field_dependent_gaussian';
C.Optics.SigmaCenter_pix=.5; 
% Additional Gaussian width toward the field edge [pixels].
C.Optics.SigmaEdgeAdd_pix=.8; 
C.Optics.PSFBanks=9;
%% 21 | DETECTOR ELECTRON STATISTICS AND DIGITAL READOUT
% Photon/dark noise, full-well saturation, read noise, conversion gain,
% bias and optional ADC clipping/quantization. Full well is an electron
% capacity; ADC range is a separate digital limit. PRNUMap/BadPixelMask
% are optional detector-sized calibration maps; [] means not supplied.
% Explicit legacy conversion factor; band_photon derives its own effective conversion.
C.Detector.RadiometricGain_e_per_unit=1e6; 
% Electrons per digital number [e-/DN]; this is not quantum efficiency.
C.Detector.ConversionGain_e_per_DN=1;
C.Detector.PoissonEnabled=true; 
C.Detector.DarkEnabled=true;
% Mean thermally generated charge rate per pixel [electrons/s].
C.Detector.DarkCurrent_e_s=.02; 
C.Detector.ReadEnabled=true; 
% RMS read noise per pixel [electrons].
C.Detector.ReadNoise_e=1.5;
C.Detector.BiasEnabled=true; 
% Electronic output offset [DN].
C.Detector.Bias_DN=100;
C.Detector.FullWellEnabled=true; 
% Maximum stored charge per pixel [electrons], when full-well clipping is enabled.
C.Detector.FullWell_e=5e5;
C.Detector.ADCEnabled=false; 
% ADC bit depth; ADCEnabled=false retains continuous-valued DN output.
C.Detector.ADCBits=12; % Active only when ADCEnabled is true.
C.Detector.PRNUMap=[]; 
C.Detector.BadPixelMask=[]; % optional measured detector maps.
%% 22 | TARGET ANNOTATION SUPPORT
% Per-target post-optics support threshold relative to that target peak.
% This defines truth annotations, not a noisy-image detection threshold.
C.Masks.RelativeThreshold=.01; % Same post-optics threshold for every target.
%% 23 | ACQUISITION EXPORT AND DATASET METADATA
% Raw products and display exports. The final output directory is assigned
% a timestamp at the end of this function. Display percentiles affect PNG
% appearance only, not detector DN or saved numerical truth.
C.Output.Directory=fullfile(root,'Outputs'); 
C.Output.SavePNG=true;
C.Output.SaveMAT=true; 
C.Output.SaveCSV=true; 
% Retain individual target maps for inspection/export; increases memory and storage.
C.Output.KeepTargetMaps=true;
C.Output.Plot=true; % Show visualizations from RUN_SPOT.
C.Output.SaveResults=true; 
C.Output.PNGPercentiles=[1 99.7]; 
C.Output.Verbose=true;
C.Dataset.Split='unspecified'; 
C.Dataset.SampleID='sample_0001';
% Optional requested regime metadata/validation field; actual regime follows geometry.
C.Dataset.RequestedRegime=''; 
% Explicit opt-in for demonstration trajectories in dataset generation.
C.Dataset.AllowDemoTrajectories=false;
% Plot settings affect visualizations without changing detector arrays.
%% 24 | VISUALIZATION STYLE AND PANELS
% Plot/export appearance and selected diagnostics. These settings do not
% change the physical acquisition.
C.Plot.Visible='on'; 
C.Plot.Save=true; 
C.Plot.DPI=300;
% MATLAB figure export formats; fig preserves editable MATLAB graphics.
C.Plot.Formats={'png','fig'}; 
C.Plot.FontName='Arial'; 
C.Plot.FontSize=11;
C.Plot.TargetFontSize=13; 
C.Plot.TruthLineWidth=2; 
C.Plot.MarkerSize=5;
C.Plot.SingleSize_px=[900 800]; 
C.Plot.StageSize_px=[1800 1200];
C.Plot.DatasetSize_px=[1500 950]; 
% One RGB row per target colour; values are in [0,1].
C.Plot.Colors=[0 .447 .741;.850 .325 .098;.929 .694 .125;.494 .184 .556;.466 .674 .188];
C.Plot.DisplayPercentiles=[1 99.7]; 
C.Plot.StarPointSize=8;
% Subsampling stride for RA/DEC diagnostic plotting; does not change the image.
C.Plot.RADECStride=8; 
C.Plot.RADECFullScale=true;
C.Plot.Camera=true; 
C.Plot.Detector=true; 
C.Plot.RADEC=true;
C.Plot.Stages=true; 
C.Plot.Dataset=true; 
C.Plot.ExportRawDetector=true;
C.Plot.Exposures=1; % Use 'all' to visualize every exposure.
C.Plot.Photometry=false;


%% 25 | SPACE PRESET OVERRIDES
% Applied AFTER the defaults above. Space presets replace observer type,
% FoV, pointing and seed rows. All seed tables use six columns:
% [angle1_deg angle2_deg range_m transverseSpeed_mps direction_deg radialSpeed_mps].
% For theta_phi seeds the first two entries are theta and phi.
if startsWith(preset,'space_') || startsWith(preset,'figure8') || strcmp(preset,'orbit_tracking') || startsWith(preset,'trajectory_')
    C.Observer.Trajectory.Type='keplerian'; C.Observer.Trajectory.Frame='GCRF';
    C.Camera.FOV_deg=[40 40]; C.Pointing.Mode='nadir_offset';
    switch preset
        case {'space_earth','space_earth_mixed'}, C.Pointing.OffNadir_deg=30;
        case {'space_limb','space_limb_mixed'}, C.Pointing.OffNadir_deg=70;
        otherwise, C.Pointing.OffNadir_deg=115;
    end
    % theta, phi, range, transverse speed, direction, radial speed.
    spaceSkySeeds=[8 20 120000 8557.386612 -35 -11082.705217; ...
    14 90 160000 6107.411529 -10 -13184.915339; ...
    20 160 200000 1964.959110 20 -15057.111580; ...
    25 230 240000 2426.520107 35 -13688.043774];
    spaceLimbSeeds=[8 20 120000 3833.620618 145 -14683.177388; ...
    14 90 160000 4990.552745 170 -14189.086109; ...
    20 160 200000 9141.945790 200 -11285.696768; ...
    25 230 240000 10458.418708 215 -11062.571364];
    spaceEarthSeeds=[8 20 120000 11578.670357 145 -7024.402262; ...
    14 90 160000 13416.799564 170 -7358.032989; ...
    20 160 200000 13873.101560 200 -3202.171618; ...
    25 230 240000 13434.782140 215 -391.550556];
    switch C.Pointing.OffNadir_deg
        case 30, spaceSeeds=spaceEarthSeeds;
        case 70, spaceSeeds=spaceLimbSeeds;
        otherwise, spaceSeeds=spaceSkySeeds;
    end
    for j=1:4, C.Targets(j).Seed=SPOT_ConfigSeed('theta_phi',spaceSeeds(j,:)); end
end

%% 26 | MIXED-APPEARANCE PRESET OVERRIDES
% Replaces physical seed rows for targets 2 and 4. Appearance still follows
% projected motion, exposure, brightness and PSF; no rendering class is set.
if endsWith(preset,'_mixed')
    % Ground compact examples are distant orbiting targets; a near-stationary
    % Cartesian body above a fixed ground observer would not be a valid orbit.
    groundMixedSeeds=[-2 -4 3e+07 3714.320301 -10 48.780975; ...
    7 -6 3.6e+07 3396.016025 35 -128.845952];
    limbMixedSeeds=[14 90 160000 148.464146 170 86.591358; ...
    25 230 240000 150.726688 215 215.879733];
    earthMixedSeeds=[14 90 160000 17.340533 350 115.010815; ...
    25 230 240000 122.067305 395 69.625379];
    switch preset
        case 'ground_mixed', mixedRows=groundMixedSeeds; mode='az_el_offset';
        case 'space_limb_mixed', mixedRows=limbMixedSeeds; mode='theta_phi';
        otherwise, mixedRows=earthMixedSeeds; mode='theta_phi';
    end
    C.Targets(2).Seed=SPOT_ConfigSeed(mode,mixedRows(1,:));
    C.Targets(4).Seed=SPOT_ConfigSeed(mode,mixedRows(2,:));
end

%% 27 | FIGURE 8 AND EXPOSURE-VARIANT OVERRIDES
% Camera-angle seeds: first two columns are xi and eta in degrees.
% The short/long variants change exposure; the ordinary figure8 preset
% retains the exposure value defined above.
if startsWith(preset,'figure8')
    % SAME target model, all speeds nonzero. T2 is mainly radial and appears
    % compact; T1/T3 have larger angular displacement. T4 tests an empty mask.
    % xi, eta [deg], range [m], transverse speed [m/s], direction [deg], radial [m/s].
    mixedSeeds=[5.6 4.1 40000 8075.916564 -35 -11461.568789; ...
    0 14 150000 2.000000 60 300.000000; ...
    -19.4 7 40000 1693.103494 25 -15106.411540; ...
    -23 -23 230000 1116.353805 -45 -14441.810356];
    for j=1:4, C.Targets(j).Seed=SPOT_ConfigSeed('camera_angles',mixedSeeds(j,:)); end
    if strcmp(preset,'figure8_short'), C.Time.Exposure_s=.001; end
    if strcmp(preset,'figure8_long')
        C.Time.Exposure_s=3;
        C.Integration.MaxBackgroundSamples=513; % Longer interval needs more limb/star samples.
    end
    C.Time.Cadence_s=max(C.Time.Cadence_s,C.Time.Exposure_s);
end

%% 28 | SEQUENCE AND TRACKING EXAMPLE OVERRIDES
% Examples for repeated acquisition, target tracking and trajectory inputs.
if strcmp(preset,'multi_exposure')
    C.Time.NumExposures=3; C.Time.Cadence_s=1;
end
if strcmp(preset,'orbit_tracking') || startsWith(preset,'trajectory_')
    C.Targets=C.Targets(1:2); C.Pointing.Mode='target_tracking';
    for j=1:2
        C.Targets(j).Trajectory.Type='keplerian'; C.Targets(j).Trajectory.Frame='GCRF';
        C.Targets(j).Trajectory.Elements=[7000000 .001 97.5 20 0 1+j];
    end
end
% Quality is explicit; preview is a different acquisition, not paper evidence.

%% 29 | QUALITY OVERRIDES AND TRAJECTORY INPUT EXAMPLES
% Preview overrides resolution, frames, stars, photon shot noise and DPI.
% It is a different acquisition, not evidence for full-quality validation.
% trajectoryExample then substitutes the requested input representation.
C.Run.Quality=quality;
if strcmp(quality,'preview')
    C.Frames.Model='gmst_precession'; C.Camera.Nx=256; C.Camera.Ny=256;
    C.Stars.Source='synthetic'; C.Stars.SyntheticCount=20000;
    C.Detector.PoissonEnabled=false; C.Plot.DPI=120;
end
C=SPOT_ConfigTrajectoryExample(C,preset);
if startsWith(preset,'space_'), C.Input.SensorSelection.FilterValue=erase(preset,'_mixed');
elseif startsWith(preset,'figure8'), C.Input.SensorSelection.FilterValue='space_sky'; end
if userProfile || strcmp(preset,'space_limb')
    % Sunlit Earth-limb scene and current_* acquisition parameters.
    % space_limb uses two-body propagation and dynamic camera integration;
    % current_* profiles use linear propagation and midpoint attitude.
    C.Name=requestedPreset; C.Time.StartUTC=datetime(2026,5,1,12,0,0,'TimeZone','UTC');
    if strcmp(preset,'ground'), C.Time.StartUTC=datetime(2026,5,1,8,0,0,'TimeZone','UTC'); end
    C.Time.Exposure_s=1; C.Time.Cadence_s=0.8; C.Stars.MaxStars=200000;
    C.Observer.Trajectory.EpochUTC=C.Time.StartUTC;
    C.Observer.Trajectory.Elements=[6878137 .001 97.5 40 0 0];
    if userProfile
        C.Seeding.Propagation='linear'; C.Integration.AttitudeMode='midpoint';
    end
    C.Radiometry.DetectorCalibration='legacy_gain'; C.Detector.RadiometricGain_e_per_unit=1e6;
    ranges=linspace(700e3,1500e3,4); speeds=[7500 7500 6000 8000];
    angles=[90 .2 90 90]; directions=[-35 0 35 -20];
    az=[-6 -2 3 7]; el=[5 -4 2 -6]; theta=[8 10 14 12]; phi=[20 90 160 230];
    for j=1:4
        C.Targets(j).Trajectory.EpochUTC=C.Time.StartUTC; 
        C.Targets(j).Attitude.EpochUTC=C.Time.StartUTC;
        C.Targets(j).Surface.Area_m2=10/6*ones(1,6);
        if strcmp(preset,'ground')
            d=SPOT_ConfigSeed('az_el_offset',[az(j) el(j) ranges(j) 0 directions(j) 0]);
        else
            d=SPOT_ConfigSeed('theta_phi',[theta(j) phi(j) (80+40*j)*1e3 0 directions(j) 0]);
        end
        d.TransverseSpeed_mps=speeds(j)*sind(angles(j)); d.RadialSpeed_mps=speeds(j)*cosd(angles(j));
        C.Targets(j).Seed=d;
    end
    preset=requestedPreset;
end
if any(strcmp(preset,{'space_limb','paper_space_limb','current_space_limb'}))
    C.Plot.DetectorColorbar=true; % Standalone detector figure: native DN scale.
end

%% 31 | FINAL SAMPLE NAME AND TIMESTAMPED OUTPUT DIRECTORY
% Last assignments before returning C. To guarantee a custom override,
% modify the returned structure before RUN_SPOT(C), e.g.
% C = SPOT_Config('figure8','full'); C.Time.Exposure_s = 0.5; RUN_SPOT(C);
C.Dataset.SampleID=preset;
C.Output.Directory=fullfile(root,'Outputs',char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
end
