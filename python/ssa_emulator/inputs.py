from copy import deepcopy
from pathlib import Path
import numpy as np
import pandas as pd
from .config import struct
from .geometry import euler,prepare_trajectory,seed_targets,utc
from .photometry import rng

def select_pool_rows(table,s,count,idcolumn):
    if count<0 or int(count)!=count:raise ValueError('Count must be nonnegative integer')
    if table[idcolumn].isna().any() or table[idcolumn].duplicated().any():raise ValueError('Pool IDs must be unique and nonempty')
    idx=np.arange(len(table))
    if s.FilterColumn:idx=idx[table[s.FilterColumn].astype(str).to_numpy()==str(s.FilterValue)]
    if s.Mode=='ordered':
        lo=int(s.StartRow)-1
        if lo<0 or lo+count>len(idx):raise ValueError('Insufficient ordered pool rows')
        idx=idx[lo:lo+count]
    elif s.Mode=='random':
        if count>len(idx):raise ValueError('Insufficient rows without replacement')
        idx=idx[rng(s.Seed).choice(len(idx),size=count,replace=False)]
    elif s.Mode=='ids':
        if len(s.IDs)!=count or len(set(s.IDs))!=count:raise ValueError('IDs must be unique and match count')
        lookup={str(table.iloc[i][idcolumn]):i for i in idx}
        try:idx=np.array([lookup[str(k)] for k in s.IDs],int)
        except KeyError as e:raise ValueError('Pool ID missing from filtered pool') from e
    else:raise ValueError('Selection must be ordered, random or ids')
    return table.iloc[idx],idx+1

def resolve_pools(c):
    if c.Input.UseSensorCSV:
        rows,idx=select_pool_rows(pd.read_csv(c.Input.SensorCSV),c.Input.SensorSelection,1,'sensor_id');q=rows.iloc[0];tr=c.Observer.Trajectory;tr.Type=q.position_type;tr.EpochUTC=c.Time.StartUTC
        if tr.Type=='geodetic':tr.Frame='ECEF';tr.Motion='fixed';tr.Geodetic_deg_m=[q.lat_deg,q.lon_deg,q.alt_m]
        elif tr.Type=='ecef':tr.Frame='ECEF';tr.Motion='fixed';tr.Position_m=[q.x_ecef_m,q.y_ecef_m,q.z_ecef_m]
        elif tr.Type=='keplerian':tr.Frame='GCRF';tr.Elements=[q.a_m,q.e,q.i_deg,q.raan_deg,q.argp_deg,q.ta_deg]
        elif tr.Type=='ephemeris_csv':
            tr.Type='ephemeris';f=Path(q.ephemeris_file);f=f if f.is_file() else Path(c.Input.SensorCSV).parent/f;tr.File=str(f);cols=pd.read_csv(f,nrows=0).columns;tr.Frame='GCRF' if 'x_gcrf_m' in cols else 'ECEF'
        else:raise ValueError('Unsupported sensor pool trajectory')
        c.Observer.Name=q['name'];c.Pointing.Mode='inertial' if q.orientation_mode=='boresight_eci' else q.orientation_mode
        c.Pointing.AzEl_deg=[q.azimuth_deg,q.elevation_deg];c.Pointing.Roll_deg=q.roll_deg;c.Pointing.OffNadir_deg=q.offnadir_deg;c.Pointing.SpaceAzimuth_deg=q.space_azimuth_deg;c.Pointing.Boresight_GCRF=[q.bx_eci,q.by_eci,q.bz_eci]
        c.Camera.Shape=q.fov_shape;c.Camera.FOV_deg=[q.fov_diameter_deg]*2 if q.fov_shape=='circular' else [q.fov_azimuth_deg,q.fov_elevation_deg];c.Input.SelectedSensorRow=int(idx[0]);c.Input.SelectedSensor=q.to_dict()
    if c.Input.UseTargetCSV:
        rows,idx=select_pool_rows(pd.read_csv(c.Input.TargetCSV),c.Input.TargetSelection,c.Input.NumTargets,'target_id')
        if len(rows) and not c.Targets:raise ValueError('Keep one target template when selecting a nonempty pool')
        template=deepcopy(c.Targets[0]) if len(rows) else None;c.Targets=[]
        for (_,q),row in zip(rows.iterrows(),idx):
            if q.photometry_model!='surface_brdf':raise ValueError('Pool must use surface_brdf')
            T=deepcopy(template);T.ID=int(row);T.Name=q.target_id;T.Class=q.target_class;tr=T.Trajectory;tr.EpochUTC=c.Time.StartUTC
            if q.source=='seeded':
                tr.Type='seeded_state';s=T.Seed
                if c.Pointing.Mode=='az_el_roll':s.AngleMode='az_el_offset';s.Angles_deg=[q.az_offset_deg,q.el_offset_deg];s.Range_m=q.range_m
                else:s.AngleMode='theta_phi';s.Angles_deg=[q.space_theta_deg,q.space_phi_deg];s.Range_m=q.space_range_m
                a=np.deg2rad(q.angle_to_los_deg);s.TransverseSpeed_mps=q.relative_speed_mps*np.sin(a);s.RadialSpeed_mps=q.relative_speed_mps*np.cos(a);s.Direction_deg=q.position_angle_deg
            elif q.source=='keplerian':tr.Type='keplerian';tr.Frame='GCRF';tr.Elements=[q.a_m,q.e,q.i_deg,q.raan_deg,q.argp_deg,q.ta_deg]
            else:raise ValueError('Unsupported target pool source')
            T.Photometry.Model='surface_brdf';T.Surface.NormalBody=np.column_stack([np.eye(3),-np.eye(3)]).tolist();T.Surface.Area_m2=[q.area_m2/6]*6;T.Surface.DiffuseReflectance=[q.albedo]*6;T.Surface.SpecularReflectance=[q.get('specular_reflectance',.1)]*6;T.Surface.SpecularExponent=[q.get('specular_exponent',40)]*6
            T.Attitude.Mode='spin';T.Attitude.EpochUTC=c.Time.StartUTC;T.Attitude.SpinRate_deg_s=q.get('spin_rate_deg_s',5);T.Attitude.SpinAxisBody=[q.get('spin_axis_x',0),q.get('spin_axis_y',0),q.get('spin_axis_z',1)];T.Attitude.BodyToGCRF0=euler([q.get('body_roll_deg',0),q.get('body_pitch_deg',0),q.get('body_yaw_deg',0)]).tolist();T.Metadata.PoolID=q.target_id;T.Metadata.PoolRow=int(row);T.Metadata.PoolRecord=q.to_dict();c.Targets.append(T)
        c.Input.SelectedTargetRows=idx.tolist();c.Input.SelectedTargets=rows.to_dict('records')
    c.Input.Resolved=True
    return c

def validate(c):
    c=struct(deepcopy(c));resolve_pools(c)
    from .examples import apply_example_overrides
    apply_example_overrides(c)
    if c.SchemaVersion!='7.0.0':raise ValueError('Use the SPOT 7 configuration schema')
    for value in [c.Time.Exposure_s,c.Time.Cadence_s,c.Integration.Step_s,c.Integration.MaxStep_pix,c.Integration.PhotometryRelTol,c.Integration.PhotometryAbsTol,c.Integration.MaxBodyRotation_deg,c.Integration.SpecularSamplesPerHWHM]:
        if not np.isfinite(value) or value<=0:raise ValueError('Timing and integration controls must be positive finite')
    utc(c.Time.StartUTC)
    for value,allowed,label in [(c.Integration.AttitudeMode,['dynamic','midpoint'],'camera attitude'),(c.Visibility.EclipseModel,['cylinder','finite_sun'],'eclipse'),(c.Atmosphere.AirmassModel,['secant','kasten_young'],'airmass'),(c.Optics.PSFModel,['gaussian','field_dependent_gaussian'],'PSF'),(c.Pointing.Mode,['az_el_roll','nadir_offset','inertial','target_tracking','attitude'],'pointing')]:
        if value not in allowed:raise ValueError('Unsupported '+label)
    if c.Integration.MaxSamples<3 or c.Integration.MaxBackgroundSamples<3:raise ValueError('Sample budgets must be >=3')
    if c.Pointing.Mode=='attitude':check_rotations(c,c.Pointing.AttitudeUTC,c.Pointing.CameraToGCRF)
    for tr in [c.Observer.Trajectory]+[T.Trajectory for T in c.Targets]:
        if tr.Type not in ('geodetic','ecef','state','ephemeris','keplerian','tle','demo','seeded_state'):raise ValueError('Unsupported trajectory')
        if tr.Frame not in ('ECEF','GCRF'):raise ValueError('Input frame must be ECEF or GCRF')
        if tr.Type=='keplerian':
            a=np.asarray(tr.Elements)
            if a.shape!=(6,) or not np.isfinite(a).all() or a[0]<=0 or not 0<=a[1]<1 or not 0<=a[2]<=180:raise ValueError('Invalid elliptic elements')
        if tr.Type in ('state','ecef'):
            r=np.asarray(tr.Position_m);v=np.asarray(tr.Velocity_mps)
            if r.shape!=(3,) or v.shape!=(3,) or not np.isfinite(np.r_[r,v]).all() or np.linalg.norm(r)==0:raise ValueError('Invalid Cartesian state')
    
    if min(c.Camera.Nx,c.Camera.Ny)<2 or any(int(v)!=v for v in [c.Camera.Nx,c.Camera.Ny,c.Time.NumExposures]) or c.Time.NumExposures<1:raise ValueError('Invalid dimensions/exposure count')
    if len(c.Camera.FOV_deg)!=2 or np.any(np.asarray(c.Camera.FOV_deg)<=0) or np.any(np.asarray(c.Camera.FOV_deg)>=180):raise ValueError('Invalid FoV')
    if c.Camera.Shape not in ('rectangular','circular') or c.Camera.SolidAngle not in ('mean','per_pixel'):raise ValueError('Unknown camera model')
    if c.Camera.Shape=='circular' and c.Camera.FOV_deg[0]!=c.Camera.FOV_deg[1]:raise ValueError('Circular FoV needs equal diameters')
    if c.Radiometry.DetectorCalibration=='band_photon':
        r=c.Radiometry
        if min(r.ZeroMagnitudeIrradiance_W_m2,r.EffectiveWavelength_m,r.ApertureDiameter_m,r.SolarIrradianceAtAU_W_m2)<=0 or not 0<r.QuantumEfficiency<=1:raise ValueError('Invalid radiometric calibration')
        c.Detector.RadiometricGain_e_per_unit=r.ZeroMagnitudeIrradiance_W_m2*np.pi*(r.ApertureDiameter_m/2)**2*r.QuantumEfficiency/(6.62607015e-34*299792458/r.EffectiveWavelength_m)
    elif c.Radiometry.DetectorCalibration!='legacy_gain':raise ValueError('Invalid calibration mode')
    if not 0<c.Masks.RelativeThreshold<=1 or c.Atmosphere.MinElevation_deg<=0 or c.Optics.PSFBanks<2:raise ValueError('Invalid threshold, airmass floor or PSF bank count')
    for val in [c.Optics.Transmission,c.Optics.EdgeTransmission]+list(c.Optics.ElementTransmission):
        if not 0<=val<=1:raise ValueError('Transmission must be in [0,1]')
    for val in [c.Detector.RadiometricGain_e_per_unit,c.Detector.ConversionGain_e_per_DN,c.Detector.FullWell_e]:
        if val<=0 or not np.isfinite(val):raise ValueError('Invalid detector gain/capacity')
    if c.Detector.ReadNoise_e<0 or c.Detector.DarkCurrent_e_s<0:raise ValueError('Negative noise')
    ids=[T.ID for T in c.Targets]
    if len(set(ids))!=len(ids) or any(i<=0 or int(i)!=i or i>2**32-1 for i in ids):raise ValueError('IDs must be unique positive uint32 integers')
    if c.Pointing.Mode=='target_tracking' and c.Pointing.TrackTargetID not in ids:raise ValueError('Tracked target is absent')
    prepare_trajectory(c,c.Observer.Trajectory);seed_targets(c)
    for T in c.Targets:
        if 'Rendering' in T:raise ValueError('Appearance is derived, not a render setting')
        prepare_trajectory(c,T.Trajectory)
        if T.Trajectory.Type=='demo' and not c.Dataset.AllowDemoTrajectories:raise ValueError('Artificial demo trajectories require explicit opt-in')
        if T.Photometry.Model not in ('surface_brdf','lambert_sphere','constant_magnitude'):raise ValueError('Unsupported photometry')
        if T.Photometry.Model=='surface_brdf':
            s=T.Surface;N=np.asarray(s.NormalBody);nf=N.shape[1]
            if N.shape[0]!=3 or not np.allclose(np.linalg.norm(N,axis=0),1,atol=1e-8,rtol=0):raise ValueError('Normals must be unit 3-by-N vectors')
            for name in ['Area_m2','DiffuseReflectance','SpecularReflectance','SpecularExponent']:
                a=np.asarray(s[name])
                if a.shape!=(nf,) or not np.all(np.isfinite(a)) or np.any(a<0):raise ValueError('Invalid face parameters')
            if min(s.Area_m2)<=0 or np.any(np.asarray(s.DiffuseReflectance)+s.SpecularReflectance>1) or min(s.SpecularExponent)<1 or max(s.SpecularExponent)>1e5:raise ValueError('Material energy/area/exponent bounds violated')
            a=T.Attitude
            if a.Mode not in ('fixed','spin','samples'):raise ValueError('Invalid target attitude mode')
            R=np.asarray(a.BodyToGCRF0)
            if R.shape!=(3,3) or not np.allclose(R.T@R,np.eye(3),atol=1e-8,rtol=0) or np.linalg.det(R)<0:raise ValueError('Body rotation must be proper orthonormal')
            if np.shape(a.SpinAxisBody)!=(3,) or not np.isfinite(a.SpinAxisBody).all() or np.linalg.norm(a.SpinAxisBody)==0 or not np.isfinite(a.SpinRate_deg_s):raise ValueError('Invalid body spin')
            if a.Mode=='samples':check_rotations(c,a.TimeUTC,a.BodyToGCRF,True)
    return c

def check_rotations(c,times,matrices,target=False):
    from .geometry import offsets
    t=offsets(c,times);a=np.asarray(matrices,float)
    if a.shape==(3,3,len(t)):a=np.moveaxis(a,2,0)
    if len(t)<2 or np.any(np.diff(t)<=0) or a.shape!=(len(t),3,3) or not np.isfinite(a).all():raise ValueError('Invalid attitude series')
    for R in a:
        if not np.allclose(R.T@R,np.eye(3),atol=1e-8,rtol=0) or np.linalg.det(R)<0:raise ValueError('Attitude samples must be proper rotations')
    if target and np.any((np.einsum('nij,nij->n',a[:-1],a[1:])-1)/2<=np.cos(np.deg2rad(179))):raise ValueError('Target attitude knots must differ by <179 degrees')
