import numpy as np
from .geometry import sample,offsets,D,unit
from .photometry import intrinsic,target_attitude

def time_grid(start,duration,step,extra=()):
    q=np.unique(np.r_[np.arange(0,duration,step),duration,duration/2,np.asarray(extra).reshape(-1)])
    q=q[(q>=0)&(q<=duration)];q=q[np.r_[True,np.diff(q)>max(64*np.spacing(duration),1e-12*duration)]];q[0]=0;q[-1]=duration;q[np.argmin(abs(q-duration/2))]=duration/2
    d=np.diff(q);w=np.r_[d[0]/2,(d[:-1]+d[1:])/2,d[-1]/2]/duration
    return start+q,w

def integrate_geometry(c,start):
    duration=c.Time.Exposure_s;mid=start+duration/2;step_limit=c.Integration.Step_s
    series=[c.Pointing.AttitudeUTC,c.Pointing.JitterUTC,c.Observer.Trajectory.TimeUTC]
    for T in c.Targets:
        series.extend([T.Trajectory.TimeUTC,T.Photometry.LightCurveUTC])
        if T.Photometry.Model!='surface_brdf':continue
        A=T.Attitude;series.append(A.TimeUTC);speed=abs(A.SpinRate_deg_s) if A.Mode=='spin' else 0
        if A.Mode=='samples':
            x=offsets(c,A.TimeUTC);R=target_attitude(c,A,x);angle=np.arccos(np.clip((np.einsum('nij,nij->n',R[:-1],R[1:])-1)/2,-1,1))/D;speed=max(angle/np.diff(x))
        if speed>0:
            step_limit=min(step_limit,c.Integration.MaxBodyRotation_deg/speed);use=np.asarray(T.Surface.SpecularReflectance)>0
            if np.any(use):half=np.arccos(2**(-1/np.max(np.asarray(T.Surface.SpecularExponent)[use])))/D;step_limit=min(step_limit,half/(2*speed*c.Integration.SpecularSamplesPerHWHM))
    if np.ceil(duration/step_limit)+1>c.Integration.MaxSamples:raise ValueError('Spin/glint sampling exceeds MaxSamples')
    extra=np.concatenate([offsets(c,s)-start for s in series if len(s)]) if any(len(s) for s in series) else []
    requested,_=time_grid(start,duration,c.Integration.Step_s);t,w=time_grid(start,duration,step_limit,extra)
    scale=max((np.array([c.Camera.Nx,c.Camera.Ny])-1)/c.Camera.FOV_deg)/np.cos(max(c.Camera.FOV_deg)*D/2)**2
    for _ in range(20):
        if len(t)>c.Integration.MaxSamples:raise ValueError('Target sampling budget exceeded')
        G=sample(c,t,mid);step=np.zeros(len(t)-1);specsplit=np.zeros(len(t)-1,bool)
        for j,H in enumerate(G.Target):
            u=np.einsum('nji,jn->in',G.CameraToGCRF,H.LOS_GCRF);ang=np.arccos(np.clip(np.sum(u[:,:-1]*u[:,1:],axis=0),-1,1))/D;d=ang*scale;near=u[2]>np.cos(min(89,max(c.Camera.FOV_deg))*D);d[~(near[:-1]|near[1:])]=0;step=np.maximum(step,d)
            T=c.Targets[j]
            if T.Photometry.Model=='surface_brdf' and np.any(np.asarray(T.Surface.SpecularReflectance)>0):
                su=unit(G.SunPosition_m-H.Position_m);st=np.arccos(np.clip(np.sum(su[:,:-1]*su[:,1:],axis=0),-1,1))/D;vt=np.arccos(np.clip(np.sum(H.LOS_GCRF[:,:-1]*H.LOS_GCRF[:,1:],axis=0),-1,1))/D
                R=target_attitude(c,T.Attitude,t);bt=np.arccos(np.clip((np.einsum('nij,nij->n',R[:-1],R[1:])-1)/2,-1,1))/D
                p=max(np.asarray(T.Surface.SpecularExponent)[np.asarray(T.Surface.SpecularReflectance)>0]);limit=np.arccos(2**(-1/p))/D/c.Integration.SpecularSamplesPerHWHM;specsplit|=2*bt+st+vt>limit
        q=t-start;qm=(q[:-1]+q[1:])/2;split=specsplit.copy()
        if c.Integration.Adaptive:split|=step>c.Integration.MaxStep_pix*1.001
        if c.Integration.PhotometryAdaptive and c.Targets:
            M=sample(c,start+qm,mid)
            for j in range(len(c.Targets)):
                a=intrinsic(c,G,j);b=intrinsic(c,M,j);f=a.IntrinsicRate*a.Visible;fm=b.IntrinsicRate*b.Visible
                coarse=(f[:-1]+f[1:])/2;fine=(f[:-1]+2*fm+f[1:])/4;tol=c.Integration.PhotometryAbsTol+c.Integration.PhotometryRelTol*max(np.max(f),np.max(fm));err=abs(fine-coarse)*np.diff(q)/duration
                need=err>tol*np.diff(q)/duration;need[(err<=tol)&(np.diff(q)<=duration*1e-5)]=False;split|=need
        if not np.any(split):return G,w,requested
        if len(t)+sum(split)>c.Integration.MaxSamples:raise ValueError('Spatial/photometric sampling exceeds MaxSamples')
        t,w=time_grid(start,duration,duration,np.r_[q,qm[split]])
    raise ValueError('Adaptive sampling did not converge')
