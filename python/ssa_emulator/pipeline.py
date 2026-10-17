from copy import deepcopy
from datetime import datetime,timezone
import platform
import numpy as np
from .config import default_config,struct,plain
from .inputs import validate
from .sampling import integrate_geometry
from .geometry import iso,D
from .photometry import photometry
from .background import load_stars,background
from .imaging import deposit,optics,detector,morphology,labels

def run(c=None,export=False,progress=None):
    c=validate(c or default_config());stars=load_stars(c)
    out=struct(dict(SchemaVersion='7.0.0',Configuration=c,Provenance=dict(Python=platform.python_version(),CreatedUTC=datetime.now(timezone.utc).isoformat(),Reference='MATLAB ZPublic-SPOT 7.2.0 (schema 7.0.0)',ParityStatus='equation-aligned; MATLAB/Python numerical parity requires exported fixture comparison',Randomness='NumPy MT19937 streams: repeatable in Python, not sample-identical to MATLAB',Frame=c.Frames.Model,DetectorCalibration=c.Radiometry.DetectorCalibration),Exposures=[]))
    for ie in range(c.Time.NumExposures):
        if progress:progress(f'Exposure {ie+1}/{c.Time.NumExposures}: geometry and background')
        start=ie*c.Time.Cadence_s;G,w,requested=integrate_geometry(c,start);B=background(c,start,stars);shape=(c.Camera.Ny,c.Camera.Nx);target=np.zeros(shape);maps=[];records=[]
        for j,T in enumerate(c.Targets):
            H=photometry(c,G,j);H.RA_deg=np.arctan2(H.LOS_GCRF[1],H.LOS_GCRF[0])/D%360;H.DEC_deg=np.arcsin(np.clip(H.LOS_GCRF[2],-1,1))/D;m=morphology(c,H,w)
            rate=deposit(c,H.x,H.y,w*H.RadiometricRate*H.Visible,B.Truth.Active);signal,_=optics(c,rate,B.Truth.Active);maps.append(signal);target+=rate;mean=float(np.sum(w*H.Visible*H.RadiometricRate));H.Time_s=G.Time_s;H.TimeUTC=G.TimeUTC;H.Weights=w
            r=struct(dict(ID=T.ID,Name=T.Name,Class=T.Class,InputTrajectory=T.Trajectory.Type,Appearance=m.Kind,Morphology=m,DisplacementBound_pix=float(np.hypot(np.ptp(H.x),np.ptp(H.y))),RelativeFluxChange=float(np.ptp(H.RadiometricRate)/max(np.max(H.RadiometricRate),1e-300)),AnyVisible=bool(H.Visible.any()),VisibleFraction=float(w@H.Visible),IlluminatedFraction=float(w@H.Illumination),MeanVisibleRate=mean,MeanReceivedIrradiance_W_m2=mean*c.Radiometry.ZeroMagnitudeIrradiance_W_m2,MeanApparentMagnitude=-2.5*np.log10(mean) if mean>0 else np.inf,InputDepositedRate=float(rate.sum()),IntegratedOpticalSignal=float(signal.sum()),Samples=H));records.append(r)
        scene=B.StarRate+B.DiffuseRate+target;signal,opt=optics(c,scene,B.Truth.Active);det=detector(c,signal,ie+1,B.Truth.Active);lab=labels(c,maps,records)
        E=struct(dict(Index=ie+1,StartUTC=iso(c,[start])[0],MidUTC=iso(c,[start+c.Time.Exposure_s/2])[0],EndUTC=iso(c,[start+c.Time.Exposure_s])[0],RequestedTime_s=requested,TargetTime_s=G.Time_s,TargetWeights=w,BackgroundTime_s=B.Time_s,BackgroundWeights=B.Weights,Observer=G.Observer,CameraToGCRF=G.CameraToGCRF,Truth=B.Truth,StarTruth=B.StarTruth,SkyRate=B.SkyRate,EarthRate=B.EarthRate,StarRate=B.StarRate,DiffuseRate=B.DiffuseRate,TargetRate=target,SceneRate=scene,OpticalStages=opt,DetectorStages=det,Image_DN=det.Image_DN,Labels=lab,Targets=records,StarMetadata=B.StarCatalogueMetadata,StarsUsed=B.StarsUsed,RequestedRegime=c.Dataset.RequestedRegime,TotalTargetOpticalSignal=sum(maps,np.zeros(shape)),TargetOpticalSignal=maps if c.Output.KeepTargetMaps else [],Sampling=dict(Adaptive=c.Integration.Adaptive,PhotometryAdaptive=c.Integration.PhotometryAdaptive,TargetSampleCount=len(w),BackgroundSampleCount=len(B.Weights),TimeConvergenceCertified=False)))
        out.Exposures.append(E)
        if progress:progress(f'Exposure {ie+1}: {B.Truth.GeometryRegime}, {sum(r.HasAnnotation for r in records)} annotated targets')
    if export and c.Output.get('SaveResults',True):
        from .export import export_run
        out.ExportDirectory=str(export_run(out))
    return out

def convergence(c):
    c=deepcopy(c);c.Detector.PoissonEnabled=False;c.Detector.ReadEnabled=False;a=run(c);c.Integration.Step_s/=2;c.Integration.MaxStep_pix/=2;c.Integration.PhotometryRelTol/=2;c.Integration.PhotometryAbsTol/=2;c.Integration.MaxBodyRotation_deg/=2;c.Integration.SpecularSamplesPerHWHM*=2;c.Integration.MaxSamples*=2;c.Integration.MaxBackgroundSamples*=2;b=run(c);rows=[]
    for x,y in zip(a.Exposures,b.Exposures):
        f=np.array([r.IntegratedOpticalSignal for r in x.Targets]);g=np.array([r.IntegratedOpticalSignal for r in y.Targets]);union=x.Labels.Binary|y.Labels.Binary
        rows.append(dict(RelativeL1=float(abs(x.OpticalStages.PSF-y.OpticalStages.PSF).sum()/max(abs(y.OpticalStages.PSF).sum(),1e-300)),TargetFluxRelativeChange=(abs(f-g)/np.maximum(g,1e-300)).tolist(),MaskIoU=float((x.Labels.Binary&y.Labels.Binary).sum()/union.sum()) if union.any() else 1.))
    return rows
