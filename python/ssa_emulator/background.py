import numpy as np
import pandas as pd
from .config import struct
from .geometry import D,unit,sample,offsets,project
from .sampling import time_grid
from .imaging import detector_grid,deposit
from .photometry import rng,atmosphere

def load_stars(c):
    empty=struct(dict(ID=np.array([],str),Magnitude=np.array([]),Direction_GCRF=np.zeros((3,0)),Metadata=dict(Enabled=c.Stars.Enabled,Source=c.Stars.Source,Count=0)))
    if not c.Stars.Enabled:return empty
    if c.Stars.Source=='csv':
        a=pd.read_csv(c.Stars.File,dtype={'source_id':str});ra=a.ra.to_numpy();dec=a.dec.to_numpy();mag=a.magnitude.to_numpy();ids=a.source_id.to_numpy() if 'source_id' in a else np.array([f'row_{i+1}' for i in range(len(a))])
    elif c.Stars.Source=='synthetic':
        r=rng(c.RandomSeed+900001);n=c.Stars.SyntheticCount;ra=360*r.random(n);dec=np.arcsin(2*r.random(n)-1)/D;lo,hi=c.Stars.SyntheticMagnitudeRange;mag=lo+(hi-lo)*r.random(n);ids=np.array([f'synthetic_{i+1}' for i in range(n)])
    else:raise ValueError('Stars.Source must be csv or synthetic; no synthetic fallback')
    if not np.all(np.isfinite([ra,dec,mag])) or np.any(abs(dec)>90):raise ValueError('Invalid star astrometry')
    use=np.where(mag<=c.Stars.MaxMagnitude)[0];count=len(use);use=use[np.argsort(mag[use],kind='stable')[:c.Stars.MaxStars]];ra=ra[use]*D;dec=dec[use]*D
    return struct(dict(ID=ids[use],Magnitude=mag[use],Direction_GCRF=np.array([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)]),Metadata=dict(Enabled=True,Source=c.Stars.Source,Count=len(use),CountBeforeCap=count,Truncated=len(use)<count,Catalogue=c.Stars.Catalogue,Band=c.Stars.Band,EpochUTC=c.Stars.EpochUTC)))

def ray_scene(c,G,k):
    d=detector_grid(c,G.CameraToGCRF[k]);shape=d.Active.shape;r=G.Observer.Position_m[:,k];u=d.LOS_GCRF;re=c.Constants.EarthRadius_m;su=unit(G.SunPosition_m[:,k]);reshape=lambda v:np.asarray(v).reshape(shape,order='F')
    d.Elevation_deg=reshape(np.arcsin(np.clip(G.Zenith_GCRF[:,k]@u,-1,1))/D);d.EarthMask=np.zeros(shape,bool);d.LimbDistance_deg=np.full(shape,np.nan);disk=np.zeros(shape);rim=disk.copy();sky=disk.copy()
    if G.IsGround[k]:
        d.SkyMask=d.Active&(d.Elevation_deg>=c.Visibility.MinElevation_deg);d.GeometryRegime=d.ConservativeRegime='ground_sky'
        if c.Background.Enabled:
            sunel=np.arcsin(np.clip(su@G.Zenith_GCRF[:,k],-1,1))/D;lo,hi=c.Background.TwilightSolarElevation_deg;f=np.clip((sunel-lo)/(hi-lo),0,1);f=f*f*(3-2*f);sky[:]=c.Background.DarkRadiance+f*(c.Background.SunlitRadiance-c.Background.DarkRadiance)
            if c.Background.LightPollutionEnabled:sky+=c.Background.LightPollutionRadiance*(1+c.Background.HorizonBoost*(1-np.sin(np.maximum(d.Elevation_deg,0)*D)))
    else:
        dotru=r@u;disc=dotru*dotru-(r@r-re*re);near=-dotru-np.sqrt(np.maximum(disc,0));far=-dotru+np.sqrt(np.maximum(disc,0));near=np.where(near<=0,far,near);hit=(disc>=0)&(near>0);d.EarthMask=reshape(hit)&d.Active;d.SkyMask=d.Active&~d.EarthMask
        alpha=np.arcsin(min(1,re/np.linalg.norm(r)))/D;nadir=-unit(r);beta=np.arccos(np.clip(nadir@u,-1,1))/D;d.LimbDistance_deg=reshape(beta-alpha)
        d.GeometryRegime='space_earth' if np.all(d.EarthMask[d.Active]) else 'space_sky' if not d.EarthMask.any() else 'space_limb';beta0=np.arccos(np.clip(nadir@G.CameraToGCRF[k,:,2],-1,1))/D
        cone=c.Camera.FOV_deg[0]/2 if c.Camera.Shape=='circular' else np.arctan(np.linalg.norm(np.tan(np.asarray(c.Camera.FOV_deg)*D/2)))/D
        d.ConservativeRegime='space_earth' if beta0+cone<alpha else 'space_sky' if beta0-cone>alpha else 'space_limb'
        if c.Background.Enabled:sky[:]=c.Background.DarkRadiance
        if c.Earth.Enabled:
            p=r[:,None]+u*near;normal=p/re;mu=np.maximum(0,su@normal);view=np.maximum(0,-np.sum(normal*u,axis=0));val=(c.Earth.NightRadiance+c.Earth.DayRadiance*mu**c.Earth.LambertGamma)*(c.Earth.LimbMinWeight+(1-c.Earth.LimbMinWeight)*np.sqrt(view));val[~hit]=0;disk=reshape(val)*d.Active
            p=r[:,None]+u*np.maximum(0,-dotru);p/=np.maximum(np.linalg.norm(p,axis=0),1e-300);mu=np.maximum(0,su@p)
            rim=c.Earth.RimGain*c.Earth.DayRadiance*np.exp(-.5*(d.LimbDistance_deg/c.Earth.RimSigma_deg)**2)*reshape(np.sqrt(mu))*(abs(d.LimbDistance_deg)<=c.Earth.RimHalfWidth_deg)*d.Active
    sd=unit(G.SunPosition_m[:,k]-r);d.SunExcluded=reshape(np.arccos(np.clip(sd@u,-1,1))/D<c.Visibility.SunExclusion_deg)&d.Active;d.SkyMask&=~d.SunExcluded;d.SkyRadiance=sky*d.SkyMask;d.EarthRadiance=(disk+rim)*~d.SunExcluded;d.DiffuseRadiance=d.SkyRadiance+d.EarthRadiance;d.EarthFraction=d.EarthMask.sum()/d.Active.sum()
    return d,d.DiffuseRadiance*d.SolidAngle_sr

def background(c,start,stars):
    dur=c.Time.Exposure_s;mid=start+dur/2;extra=[]
    for series in [c.Pointing.JitterUTC,c.Pointing.AttitudeUTC]:
        if series:extra.extend(offsets(c,series)-start)
    t,w=time_grid(start,dur,c.Integration.Step_s,extra);scale=max((np.array([c.Camera.Nx,c.Camera.Ny])-1)/c.Camera.FOV_deg)/np.cos(max(c.Camera.FOV_deg)*D/2)**2
    for it in range(15):
        if len(t)>c.Integration.MaxBackgroundSamples:raise ValueError('Background sample budget exceeded')
        G=sample(c,t,mid);angles=np.arccos(np.clip((np.einsum('nij,nij->n',G.CameraToGCRF[:-1],G.CameraToGCRF[1:])-1)/2,-1,1))/D;r=G.Observer.Position_m
        h=np.maximum(1000,np.minimum(np.linalg.norm(r[:,:-1],axis=0),np.linalg.norm(r[:,1:],axis=0))-c.Constants.EarthRadius_m);earth=np.linalg.norm(np.diff(r,axis=1),axis=0)/h/D;earth[G.IsGround[:-1]&G.IsGround[1:]]=0;step=(angles+earth)*scale
        if not c.Integration.Adaptive or np.all(step<=c.Integration.MaxStep_pix*1.001):break
        q=t-start;add=[]
        for k in np.where(step>c.Integration.MaxStep_pix)[0]:
            n=max(2,int(np.ceil(step[k]/c.Integration.MaxStep_pix)))
            if n+len(t)>c.Integration.MaxBackgroundSamples:raise ValueError('Background refinement budget exceeded')
            add.extend(np.linspace(q[k],q[k+1],n+1)[1:-1])
        t,w=time_grid(start,dur,dur,np.r_[q,add])
    else:raise ValueError('Background grid did not converge')
    half=c.Camera.FOV_deg[0]/2 if c.Camera.Shape=='circular' else np.arctan(np.linalg.norm(np.tan(np.asarray(c.Camera.FOV_deg)*D/2)))/D
    keep=np.zeros(len(stars.Magnitude),bool)
    for A in G.CameraToGCRF:keep|=A[:,2]@stars.Direction_GCRF>=np.cos(half*D)-1e-10
    dirs=stars.Direction_GCRF[:,keep];mag=stars.Magnitude[keep];ids=stars.ID[keep];ns=len(mag);shape=(c.Camera.Ny,c.Camera.Nx)
    B=struct({k:np.zeros(shape) for k in ['StarRate','DiffuseRate','SkyRate','EarthRate']});B.StarTruth=struct({k:np.array([]) for k in ['x','y','Xi_deg','Eta_deg','RA_deg','DEC_deg','Magnitude']});rate=np.zeros((len(t),ns));els=rate.copy();pos=[]
    for k in range(len(t)):
        d,diffuse=ray_scene(c,G,k);B.DiffuseRate+=w[k]*diffuse;B.SkyRate+=w[k]*d.SkyRadiance*d.SolidAngle_sr;B.EarthRate+=w[k]*d.EarthRadiance*d.SolidAngle_sr
        ismid=abs(t[k]-mid)<1e-7
        if ismid:B.Truth=d
        if not ns:continue
        p=project(c,G.CameraToGCRF[k],dirs);pos.append(p);el=np.arcsin(np.clip(G.Zenith_GCRF[:,k]@dirs,-1,1))/D;els[k]=el;visible=p.InFOV.copy()
        if G.IsGround[k]:visible&=el>=c.Visibility.MinElevation_deg
        elif c.Visibility.EarthOccultation:
            r=G.Observer.Position_m[:,k];dot=r@dirs;visible&=~((dot<0)&(dot*dot>=r@r-c.Constants.EarthRadius_m**2))
        su=unit(G.SunPosition_m[:,k]-G.Observer.Position_m[:,k]);visible&=np.arccos(np.clip(su@dirs,-1,1))/D>=c.Visibility.SunExclusion_deg;rate[k]=10**(-.4*mag)*visible
        if ismid:
            B.StarTruth=struct({name:p[name][visible] for name in ('x','y','Xi_deg','Eta_deg')});u=dirs[:,visible];B.StarTruth.update(RA_deg=np.arctan2(u[1],u[0])/D%360,DEC_deg=np.arcsin(np.clip(u[2],-1,1))/D,Magnitude=mag[visible])
    for j in np.where(np.any(rate>0,axis=0))[0]:
        key=7919
        for ch in str(ids[j]):key=(key*31+ord(ch))%2**32
        at,sc=atmosphere(c,t,els[:,j],G.IsGround,(key+100000)%2**32);rate[:,j]*=at*sc
    for k in range(len(t)):
        if ns:B.StarRate+=w[k]*deposit(c,pos[k].x,pos[k].y,rate[k],B.Truth.Active)
    B.Time_s=t;B.Weights=w;B.StarsUsed=int(np.any(rate>0,axis=0).sum());B.StarCatalogueMetadata=struct(dict(stars.Metadata,WorkingSetCount=ns));return B
