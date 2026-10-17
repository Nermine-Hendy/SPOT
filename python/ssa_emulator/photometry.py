"""Surface reflection and atmosphere: MATLAB SPOT 7 equations and units."""
import numpy as np
from scipy.spatial.transform import Rotation
from .config import struct
from .geometry import offsets,utc,rotation_series,D,unit

def target_attitude(c,A,t):
    t=np.atleast_1d(t)
    if A.Mode=='samples':return rotation_series(offsets(c,A.TimeUTC),A.BodyToGCRF,t)
    dt=t-(utc(A.EpochUTC)-utc(c.Time.StartUTC)).total_seconds()
    angle=dt*A.SpinRate_deg_s*D if A.Mode=='spin' else dt*0
    return np.asarray(A.BodyToGCRF0)@Rotation.from_rotvec(angle[:,None]*unit(A.SpinAxisBody)).as_matrix()
def surface_reflection(S,R,s,v,solar,ran):
    n=np.einsum('nij,jf->nif',R,np.asarray(S.NormalBody));ci=np.einsum('nif,in->nf',n,s);co=np.einsum('nif,in->nf',n,v)
    reflected=2*n*ci[:,None,:]-s.T[:,:,None];alignment=np.clip(np.einsum('nif,in->nf',reflected,v),0,1)
    geometry=np.asarray(S.Area_m2)*np.maximum(ci,0)*np.maximum(co,0)
    fd=np.asarray(S.DiffuseReflectance)/np.pi;p=np.asarray(S.SpecularExponent);fs=np.asarray(S.SpecularReflectance)*(p+2)/(2*np.pi)*alignment**p
    factor=np.asarray(solar)/np.asarray(ran)**2;diffuse=factor*np.sum(geometry*fd,axis=1);spec=factor*np.sum(geometry*fs,axis=1)
    return struct(dict(DiffuseIrradiance_W_m2=diffuse,SpecularIrradiance_W_m2=spec,Irradiance_W_m2=diffuse+spec,IncidenceCosine=ci.T,EmissionCosine=co.T,SpecularAlignment=alignment.T))
def intrinsic(c,G,j):
    h=struct(dict(G.Target[j]));T=c.Targets[j];P=T.Photometry;h.LambertPhase=np.ones_like(h.Range_m)
    if P.Model=='surface_brdf':
        s=G.SunPosition_m-h.Position_m;ds=np.linalg.norm(s,axis=0);R=target_attitude(c,T.Attitude,G.Time_s)
        solar=c.Radiometry.SolarIrradianceAtAU_W_m2*(c.Constants.AU_m/ds)**2
        h.Reflection=surface_reflection(T.Surface,R,s/ds,-h.LOS_GCRF,solar,np.maximum(h.Range_m,1e-300));h.BodyToGCRF=R
        h.UnattenuatedRate=h.Reflection.Irradiance_W_m2/c.Radiometry.ZeroMagnitudeIrradiance_W_m2
    else:
        a=np.pi-h.PhaseAngle_rad;phase=np.maximum(0,(np.sin(a)-a*np.cos(a))/np.pi);small=a<1e-3;phase[small]=(a[small]**3/3-a[small]**5/30+a[small]**7/840)/np.pi
        h.LambertPhase=np.ones_like(a) if P.Model=='constant_magnitude' else phase
        h.UnattenuatedRate=10**(-.4*P.ReferenceMagnitude)*(P.ReferenceRange_m/np.maximum(h.Range_m,1e-300))**2
    h.LightCurveGain=np.ones_like(h.Range_m)
    if P.LightCurveUTC:
        x=offsets(c,P.LightCurveUTC)
        if min(G.Time_s)<x[0]-1e-7 or max(G.Time_s)>x[-1]+1e-7:raise ValueError('Light curve does not cover exposure')
        h.LightCurveGain=np.interp(G.Time_s,x,P.LightCurveGain)
    h.IntrinsicRate=h.UnattenuatedRate*h.LambertPhase*h.Illumination*h.LightCurveGain
    return h

def rng(seed):return np.random.Generator(np.random.MT19937(int(seed)%2**32))
def atmosphere(c,t,el,ground,key):
    el=np.asarray(el);ground=np.asarray(ground,bool);trans=np.ones_like(el);gain=np.ones_like(el)
    if not c.Atmosphere.Enabled or not np.any(ground):return trans,gain
    eff=np.maximum(el,c.Atmosphere.MinElevation_deg);X=1/np.maximum(np.sin(eff*D),1e-6)
    if c.Atmosphere.AirmassModel=='kasten_young':X=1/(np.sin(eff*D)+.50572*(eff+6.07995)**-1.6364)
    trans[ground]=10**(-.4*c.Atmosphere.Extinction_mag_airmass*X[ground]);S=c.Atmosphere.Scintillation
    if not S.Enabled or S.Sigma0==0:return trans,gain
    t=np.asarray(t);q=np.maximum(0,t);h=S.KnotStep_s;nk=int(np.ceil(max(q)/h))+2;rho=np.exp(-h/S.CorrelationTime_s)
    z=rng(c.RandomSeed+104729*key).standard_normal(nk)
    for i in range(1,nk):z[i]=rho*z[i-1]+np.sqrt(1-rho*rho)*z[i]
    i=np.floor(q/h).astype(int);f=q/h-i;zq=((1-f)*z[i]+f*z[i+1])/np.sqrt((1-f)**2+f*f+2*rho*f*(1-f))
    sigma=S.Sigma0/np.maximum(np.sin(eff*D),S.MinSin);g=np.exp(sigma*zq-.5*sigma*sigma);gain[ground]=g[ground]
    return trans,gain

def photometry(c,G,j):
    h=intrinsic(c,G,j);h.AtmosphericTransmission,h.ScintillationGain=atmosphere(c,G.Time_s,h.Elevation_deg,G.IsGround,c.Targets[j].ID)
    h.RadiometricRate=h.IntrinsicRate*h.AtmosphericTransmission*h.ScintillationGain;h.ReceivedIrradiance_W_m2=h.RadiometricRate*c.Radiometry.ZeroMagnitudeIrradiance_W_m2
    h.ApparentMagnitude=np.full_like(h.RadiometricRate,np.inf);p=h.RadiometricRate>0;h.ApparentMagnitude[p]=-2.5*np.log10(h.RadiometricRate[p])
    return h
