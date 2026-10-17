"""Unified SPOT 7 configuration. Field names follow MATLAB 7.0.0."""
from copy import deepcopy
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
class Struct(dict):
    def __getattr__(self,k):
        try: return self[k]
        except KeyError as e: raise AttributeError(k) from e
    __setattr__=dict.__setitem__
def struct(x):
    if isinstance(x,dict): return Struct({k:struct(v) for k,v in x.items()})
    if isinstance(x,list): return [struct(v) for v in x]
    return x
def plain(x):
    if isinstance(x,dict): return {k:plain(v) for k,v in x.items() if not k.startswith('_')}
    if isinstance(x,(list,tuple)): return [plain(v) for v in x]
    if isinstance(x,np.ndarray): return plain(x.tolist())
    if isinstance(x,np.generic): return x.item()
    if isinstance(x,Path): return str(x)
    return x
def save_config(c,path): Path(path).write_text(json.dumps(plain(c),indent=2))
def load_config(path):
    c=struct(json.loads(Path(path).read_text()));old_root=Path(c.get('RootDir',ROOT));c.RootDir=str(ROOT)
    if isinstance(c.Targets,dict):c.Targets=[c.Targets]
    if c.Visibility.RangeLimits_m[1] is None:c.Visibility.RangeLimits_m[1]=float('inf')
    # Rebase packaged paths even when the old installation still exists.
    locations=[(c.Stars,'File'),(c.Input,'SensorCSV'),(c.Input,'TargetCSV'),(c.Output,'Directory')]
    locations += [(t,'File') for t in [c.Observer.Trajectory]+[v.Trajectory for v in c.Targets]]
    for parent,key in locations:
        if not parent.get(key):continue
        old=Path(parent[key])
        try:parent[key]=str(ROOT/old.relative_to(old_root))
        except ValueError:
            if not old.is_absolute():parent[key]=str(ROOT/old)
    return c
NAMES=['ground','space_sky','space_limb','space_earth','figure8','figure8_short','figure8_long','ground_mixed','space_limb_mixed','space_earth_mixed','multi_exposure','orbit_tracking','trajectory_geodetic','trajectory_ecef','trajectory_state','trajectory_ephemeris','trajectory_ephemeris_ecef','trajectory_ephemeris_csv','trajectory_keplerian','trajectory_tle','trajectory_tle_deep_space','current_ground','current_space_sky','current_space_limb','current_space_earth']
def default_config(preset='ground',quality='full'):
    if preset=='list': return NAMES.copy()
    aliases={'ground_default':'ground','space_sky_only':'space_sky','space_limb_crossing':'space_limb','space_earth_dominated':'space_earth','paper_mixed':'figure8'}
    preset=aliases.get(preset,preset.removeprefix('paper_'))
    if preset not in NAMES: raise ValueError(f'Unknown preset: {preset}')
    if quality not in ('preview','full'): raise ValueError('quality must be full or preview')
    original=preset;current=preset.startswith('current_');preset=preset.removeprefix('current_')
    epoch='2026-01-01T12:00:00Z'
    tr=dict(Type='geodetic',Frame='ECEF',EpochUTC=epoch,Geodetic_deg_m=[-37.8136,144.9631,50],Position_m=[0,0,0],Velocity_mps=[0,0,0],VelocityENU_mps=[0,0,0],Motion='fixed',Elements=[6878137,.001,97.5,20,0,0],File='',TimeUTC=[],Positions_m=[],Velocities_mps=[],Interpolation='hermite',Extrapolation='error',TLELines=[],TLEIndex=1,MaxTLEAge_days=14,Demo=dict(AngleMode='az_el_offset',MotionModel='linear_angles',Offset_deg=[-6,5],Range_m=900e3,AngularRate_deg_s=[0,0],ReferenceUTC='2026-01-01T12:00:00.100000Z'))
    selection=dict(Mode='ordered',StartRow=1,Seed=202,IDs=[],FilterColumn='',FilterValue='')
    tgt=dict(ID=1,Name='Target_1',Class='RSO',Trajectory=deepcopy(tr),Metadata={},Seed={},Photometry=dict(Model='surface_brdf',ReferenceMagnitude=4,ReferenceRange_m=1e6,LightCurveUTC=[],LightCurveGain=[]),Surface=dict(NormalBody=np.column_stack([np.eye(3),-np.eye(3)]).tolist(),Area_m2=[.1]*6,DiffuseReflectance=[.5]*6,SpecularReflectance=[.1]*6,SpecularExponent=[40]*6),Attitude=dict(Mode='spin',EpochUTC=epoch,BodyToGCRF0=np.eye(3).tolist(),SpinAxisBody=[0,0,1],SpinRate_deg_s=5,TimeUTC=[],BodyToGCRF=[]))
    c=struct(dict(SchemaVersion='7.0.0',SoftwareVersion='7.2.0-py.1',Name=original,RootDir=str(ROOT),RandomSeed=1,
      Time=dict(StartUTC=epoch,Exposure_s=1.5,Cadence_s=1.5,NumExposures=1),
      Integration=dict(Step_s=.02,Adaptive=True,MaxStep_pix=.35,MaxSamples=8193,MaxBackgroundSamples=257,AttitudeMode='dynamic',PhotometryAdaptive=True,PhotometryRelTol=.01,PhotometryAbsTol=1e-12,MaxBodyRotation_deg=1,SpecularSamplesPerHWHM=8),
      Radiometry=dict(Band='effective visible band',SolarIrradianceAtAU_W_m2=500,ZeroMagnitudeIrradiance_W_m2=3.6e-8,EffectiveWavelength_m=550e-9,QuantumEfficiency=.65,ApertureDiameter_m=.05,DetectorCalibration='band_photon'),
      Input=dict(UseSensorCSV=False,UseTargetCSV=False,SensorCSV=str(ROOT/'Inputs/SensorPool.csv'),TargetCSV=str(ROOT/'Inputs/TargetPool_Seeded.csv'),NumTargets=4,SensorSelection=deepcopy(selection),TargetSelection=deepcopy(selection)),
      Morphology=dict(MinSmear_pix=1,PSFFWHMFactor=1),Seeding=dict(Propagation='two_body'),
      Constants=dict(EarthRadius_m=6371008.8,EarthA_m=6378137,EarthEccentricity=.081819190842622,Mu_m3_s2=3.986004418e14,AU_m=149597870700.,SunRadius_m=695700e3),
      Frames=dict(Model='iau2006',DeltaAT_s=37,DeltaUT1_s=0,PolarMotion_rad=[0,0],CelestialFrame='GCRF'),
      Visibility=dict(PlatformAltitude_m=100e3,MinElevation_deg=0,RangeLimits_m=[1,float('inf')],EarthOccultation=True,SunExclusion_deg=0,EclipseModel='cylinder'),
      Camera=dict(Nx=1024,Ny=1024,FOV_deg=[20,20],Shape='rectangular',SolidAngle='mean',DiagnosticGridStep_deg=.25),
      Pointing=dict(Mode='az_el_roll',AzEl_deg=[180,60],Roll_deg=0,OffNadir_deg=70,SpaceAzimuth_deg=0,Boresight_GCRF=[1,0,0],ReferenceUp_GCRF=[0,0,1],TrackTargetID=1,AttitudeUTC=[],CameraToGCRF=[],JitterUTC=[],Jitter_arcsec=[]),
      Observer=dict(Name='Observer',Trajectory=tr),Targets=[deepcopy(tgt) for _ in range(4)],
      Stars=dict(Enabled=True,Source='csv',File=str(ROOT/'Inputs/Stars_648000.csv'),Frame='GCRF',Catalogue='Gaia DR3 (supplied subset)',Band='Gaia G',EpochUTC='2016-01-01T00:00:00Z',MaxMagnitude=12,MaxStars=648000,SyntheticCount=20000,SyntheticMagnitudeRange=[4,12]),
      Atmosphere=dict(Enabled=True,AirmassModel='secant',Extinction_mag_airmass=.2,MinElevation_deg=1,Scintillation=dict(Enabled=True,Sigma0=.08,MinSin=.2,CorrelationTime_s=.05,KnotStep_s=.01)),
      Background=dict(Enabled=True,DarkRadiance=1e-7,TwilightSolarElevation_deg=[-18,0],SunlitRadiance=1e-4,LightPollutionEnabled=True,LightPollutionRadiance=5e-7,HorizonBoost=3),
      Earth=dict(Enabled=True,DayRadiance=350,NightRadiance=2,RadianceConvention='original_figure_code_normalized',LambertGamma=1,LimbMinWeight=.25,RimGain=3,RimSigma_deg=1.6,RimHalfWidth_deg=8),
      Optics=dict(Enabled=True,Transmission=.85,ElementTransmission=[],VignettingEnabled=True,EdgeTransmission=.6,RadialExponent=2,PSFEnabled=True,PSFModel='field_dependent_gaussian',SigmaCenter_pix=.5,SigmaEdgeAdd_pix=.8,PSFBanks=9),
      Detector=dict(RadiometricGain_e_per_unit=1e6,ConversionGain_e_per_DN=1,PoissonEnabled=True,DarkEnabled=True,DarkCurrent_e_s=.02,ReadEnabled=True,ReadNoise_e=1.5,BiasEnabled=True,Bias_DN=100,FullWellEnabled=True,FullWell_e=5e5,ADCEnabled=False,ADCBits=12,PRNUMap=[],BadPixelMask=[]),
      Masks=dict(RelativeThreshold=.01),Output=dict(Directory=str(ROOT/'Outputs'),SaveResults=True,SavePNG=True,SaveNPZ=True,SaveCSV=True,KeepTargetMaps=True,Plot=True,PNGPercentiles=[1,99.7],Verbose=True),
      Dataset=dict(Split='unspecified',SampleID=original,RequestedRegime='',AllowDemoTrajectories=False),
      Plot=dict(DPI=300,Save=True,Formats=['png','pdf'],Exposures=1,Products=['camera','detector','radec','stages','dataset'],Photometry=False,DisplayPercentiles=[1,99.7],Colors=['#0072bd','#d95319','#edb120','#7e2f8e','#77ac30']),Run=dict(Quality=quality)))
    c.Input.SensorSelection.update(Seed=101,FilterColumn='profile',FilterValue='ground')
    def seed(mode,row): return struct(dict(AngleMode=mode,Angles_deg=row[:2],Range_m=row[2],TransverseSpeed_mps=row[3],Direction_deg=row[4],RadialSpeed_mps=row[5]))
    data=json.loads((ROOT/'Inputs/preset_seeds.json').read_text())
    for j,T in enumerate(c.Targets): T.ID=j+1;T.Name=f'Target_{j+1}';T.Trajectory.Type='seeded_state';T.Seed=seed('az_el_offset',data['groundSeeds'][j])
    if preset.startswith(('space_','figure8','trajectory_')) or preset=='orbit_tracking':
        c.Observer.Trajectory.Type='keplerian';c.Observer.Trajectory.Frame='GCRF';c.Camera.FOV_deg=[40,40];c.Pointing.Mode='nadir_offset'
        off=30 if preset.startswith('space_earth') else 70 if preset.startswith('space_limb') else 115;c.Pointing.OffNadir_deg=off
        key={30:'spaceEarthSeeds',70:'spaceLimbSeeds',115:'spaceSkySeeds'}[off]
        for j,T in enumerate(c.Targets): T.Seed=seed('theta_phi',data[key][j])
        c.Input.SensorSelection.FilterValue='space_earth' if off==30 else 'space_limb' if off==70 else 'space_sky'
    if preset.endswith('_mixed'):
        key={'ground_mixed':'groundMixedSeeds','space_limb_mixed':'limbMixedSeeds','space_earth_mixed':'earthMixedSeeds'}[preset]
        for j,row in zip([1,3],data[key]): c.Targets[j].Seed=seed('az_el_offset' if preset=='ground_mixed' else 'theta_phi',row)
    if preset.startswith('figure8'):
        for j,T in enumerate(c.Targets): T.Seed=seed('camera_angles',data['mixedSeeds'][j])
        if preset=='figure8_short': c.Time.Exposure_s=.001
        if preset=='figure8_long': c.Time.Exposure_s=3;c.Integration.MaxBackgroundSamples=513
        c.Time.Cadence_s=max(c.Time.Cadence_s,c.Time.Exposure_s)
    if preset=='multi_exposure': c.Time.NumExposures=3;c.Time.Cadence_s=1
    if preset=='orbit_tracking' or preset.startswith('trajectory_'):
        c.Targets=c.Targets[:2];c.Pointing.Mode='target_tracking'
        for j,T in enumerate(c.Targets): T.Trajectory.Type='keplerian';T.Trajectory.Frame='GCRF';T.Trajectory.Elements=[7e6,.001,97.5,20,0,j+2]
    if quality=='preview': c.Frames.Model='gmst_precession';c.Camera.Nx=c.Camera.Ny=256;c.Stars.Source='synthetic';c.Detector.PoissonEnabled=False;c.Plot.DPI=120
    # Sunlit limb scene: two-body propagation and dynamic integration.
    # current_* profiles use linear propagation and midpoint attitude.
    if current or preset=='space_limb':
        epoch='2026-05-01T08:00:00Z' if preset=='ground' else '2026-05-01T12:00:00Z'
        c.Time.StartUTC=epoch;c.Time.Exposure_s=1;c.Time.Cadence_s=.8;c.Stars.MaxStars=200000;c.Observer.Trajectory.EpochUTC=epoch;c.Observer.Trajectory.Elements=[6878137,.001,97.5,40,0,0]
        if current:
            c.Seeding.Propagation='linear';c.Integration.AttitudeMode='midpoint'
        c.Radiometry.DetectorCalibration='legacy_gain'
        for j,T in enumerate(c.Targets):
            T.Trajectory.EpochUTC=epoch;T.Attitude.EpochUTC=epoch;T.Surface.Area_m2=[10/6]*6
            a=np.deg2rad([90,.2,90,90][j]);speed=[7500,7500,6000,8000][j];direction=[-35,0,35,-20][j]
            row=([-6,-2,3,7][j],[5,-4,2,-6][j],np.linspace(700e3,1500e3,4)[j]) if preset=='ground' else ([8,10,14,12][j],[20,90,160,230][j],(120+40*j)*1e3)
            T.Seed=seed('az_el_offset' if preset=='ground' else 'theta_phi',list(row)+[speed*np.sin(a),direction,speed*np.cos(a)])
    if preset=='space_limb':
        c.Plot.DetectorColorbar=True
    if preset.startswith('trajectory_'):
        from .geometry import trajectory_example
        c=trajectory_example(c,preset)
    return c
# Configuration type alias.
Config=Struct
