"""GCRF state, camera geometry and exposure samples. Vectors use 3-by-N arrays."""
from datetime import datetime,timezone,timedelta
from copy import deepcopy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.spatial.transform import Rotation,Slerp
from .config import struct,Struct
D=np.pi/180
def unit(x,axis=0):
    x=np.asarray(x,float);n=np.linalg.norm(x,axis=axis,keepdims=True)
    if np.any(n==0): raise ValueError('Cannot normalize zero vector')
    return x/n
def rx(d): return Rotation.from_rotvec([d*D,0,0]).as_matrix()
def ry(d): return Rotation.from_rotvec([0,d*D,0]).as_matrix()
def rz(d): return Rotation.from_rotvec([0,0,d*D]).as_matrix()
def euler(a): return rz(a[2])@ry(a[1])@rx(a[0])
def utc(s):
    d=datetime.fromisoformat(str(s).replace('Z','+00:00'))
    if d.tzinfo is None: raise ValueError('UTC time needs timezone')
    return d.astimezone(timezone.utc)
def offsets(c,values): return np.array([(utc(v)-utc(c.Time.StartUTC)).total_seconds() for v in values])
def iso(c,t): return [(utc(c.Time.StartUTC)+timedelta(seconds=float(v))).isoformat().replace('+00:00','Z') for v in np.atleast_1d(t)]
def jd(c,t): return 2440587.5+utc(c.Time.StartUTC).timestamp()/86400+np.asarray(t)/86400

def precession(c,t):
    T=(jd(c,np.atleast_1d(t))-2451545)/36525
    zeta=(2306.2181*T+.30188*T*T+.017998*T**3)/3600
    z=(2306.2181*T+1.09468*T*T+.018203*T**3)/3600;th=(2004.3109*T-.42665*T*T-.041833*T**3)/3600
    return np.array([rz(a)@ry(-b)@rz(g) for a,b,g in zip(z,th,zeta)])
def rotation(c,t):
    t=np.atleast_1d(t);j=jd(c,t)
    if c.Frames.Model=='iau2006':
        import erfa
        # Explicit constant EOP inputs, never fetch IERS or silently substitute.
        tt=(c.Frames.DeltaAT_s+32.184)/86400;ut=c.Frames.DeltaUT1_s/86400
        return np.swapaxes(erfa.c2t06a(np.full(len(t),2451545.),j-2451545+tt,np.full(len(t),2451545.),j-2451545+ut,*c.Frames.PolarMotion_rad),1,2)
    if c.Frames.Model!='gmst_precession': raise ValueError('Unknown frame model')
    T=(j-2451545)/36525;th=(280.46061837+360.98564736629*(j-2451545)+.000387933*T*T-T**3/38710000)%360
    return np.array([q.T@rz(a) for q,a in zip(precession(c,t),th)])
def sun(c,t):
    d=jd(c,np.atleast_1d(t))-2451545;L=(280.460+.9856474*d)*D;g=(357.528+.9856003*d)*D
    lam=L+(1.915*np.sin(g)+.020*np.sin(2*g))*D;ob=(23.439-.0000004*d)*D
    v=np.array([np.cos(lam),np.cos(ob)*np.sin(lam),np.sin(ob)*np.sin(lam)])
    r=c.Constants.AU_m*(1.00014-.01671*np.cos(g)-.00014*np.cos(2*g))
    return np.einsum('nji,jn->in',precession(c,t),v)*r

def geodetic(c,v,inverse=False):
    v=np.asarray(v,float).reshape(3,-1);a=c.Constants.EarthA_m;e2=c.Constants.EarthEccentricity**2
    if not inverse:
        lat,lon,h=v;lat=lat*D;lon=lon*D;N=a/np.sqrt(1-e2*np.sin(lat)**2)
        return np.array([(N+h)*np.cos(lat)*np.cos(lon),(N+h)*np.cos(lat)*np.sin(lon),(N*(1-e2)+h)*np.sin(lat)])
    x,y,z=v;p=np.hypot(x,y);lat=np.arctan2(z,p*(1-e2));lon=np.arctan2(y,x)
    for _ in range(15):
        N=a/np.sqrt(1-e2*np.sin(lat)**2);new=np.arctan2(z+e2*N*np.sin(lat),p)
        delta=np.max(abs(new-lat));lat=new
        if delta<1e-13: break
    N=a/np.sqrt(1-e2*np.sin(lat)**2);h=p/np.maximum(np.cos(lat),np.finfo(float).eps)-N;h[p<1e-7]=abs(z[p<1e-7])-a*np.sqrt(1-e2)
    return np.array([lat/D,lon/D,h])
def enu(lat,lon):
    a,b=lat*D,lon*D
    return np.array([[-np.sin(b),-np.sin(a)*np.cos(b),np.cos(a)*np.cos(b)],[np.cos(b),-np.sin(a)*np.sin(b),np.cos(a)*np.sin(b)],[0,np.cos(a),np.sin(a)]])
def prepare_trajectory(c,tr):
    if tr.Type=='ephemeris' and tr.File:
        import pandas as pd
        a=pd.read_csv(tr.File);tr.TimeUTC=a.time_utc.tolist();f=tr.Frame.lower()
        tr.Positions_m=a[[f'{x}_{f}_m' for x in 'xyz']].to_numpy().T.tolist()
        keys=[f'v{x}_{f}_mps' for x in 'xyz'];tr.Velocities_mps=a[keys].to_numpy().T.tolist() if all(k in a for k in keys) else []
    if tr.Type=='ephemeris':
        x=offsets(c,tr.TimeUTC)
        if len(x)<2 or np.any(np.diff(x)<=0) or np.shape(tr.Positions_m)!=(3,len(x)): raise ValueError('Invalid ephemeris samples')
        if tr.Extrapolation!='error': raise ValueError('Ephemeris extrapolation unsupported')
    if tr.Type=='tle':
        from pathlib import Path
        from sgp4.api import Satrec
        lines=tr.TLELines or Path(tr.File).read_text().splitlines();lines=[s.strip() for s in lines if s.strip()]
        pairs=[(a,lines[i+1]) for i,a in enumerate(lines[:-1]) if a.startswith('1 ')];a,b=pairs[tr.TLEIndex-1]
        for line in (a,b):
            if len(line)!=69 or sum(int(x) if x.isdigit() else 1 if x=='-' else 0 for x in line[:68])%10!=int(line[-1]): raise ValueError('Invalid TLE checksum')
        sat=Satrec.twoline2rv(a,b);finish=(c.Time.NumExposures-1)*c.Time.Cadence_s+c.Time.Exposure_s
        if np.max(abs(jd(c,[0,finish])-(sat.jdsatepoch+sat.jdsatepochF)))>tr.MaxTLEAge_days: raise ValueError('TLE exceeds allowed age')
        tr._sat=sat
    return tr

def state(c,tr,t):
    t=np.atleast_1d(t).astype(float);dt=t-(utc(tr.EpochUTC)-utc(c.Time.StartUTC)).total_seconds();n=len(t);frame=tr.Frame;mu=c.Constants.Mu_m3_s2
    if tr.Type=='geodetic':
        r0=geodetic(c,tr.Geodetic_deg_m);vel=enu(*tr.Geodetic_deg_m[:2])@np.asarray(tr.VelocityENU_mps)
        r=np.repeat(r0,n,axis=1);v=np.zeros((3,n))
        if tr.Motion=='linear':r+=vel[:,None]*dt;v[:]=vel[:,None]
        frame='ECEF'
    elif tr.Type in ('ecef','state'):
        r0=np.asarray(tr.Position_m,float);v0=np.asarray(tr.Velocity_mps,float);r=np.repeat(r0[:,None],n,axis=1);v=np.zeros_like(r)
        if tr.Motion=='linear':r+=v0[:,None]*dt;v[:]=v0[:,None]
        elif tr.Motion=='two_body':
            if frame!='GCRF':raise ValueError('Two-body input must be GCRF')
            fun=lambda _,y: np.r_[y[3:],-mu*y[:3]/np.linalg.norm(y[:3])**3]
            # Integrate each side of epoch once and interpolate (same ODE/tolerances).
            for sign in (-1,1):
                use=dt*sign>0
                if np.any(use):
                    end=dt[use].min() if sign<0 else dt[use].max();sol=solve_ivp(fun,[0,end],np.r_[r0,v0],rtol=1e-11,atol=1e-5,dense_output=True)
                    if not sol.success:raise ValueError(sol.message)
                    y=sol.sol(dt[use]);r[:,use]=y[:3];v[:,use]=y[3:]
            v[:,dt==0]=v0[:,None]
        elif tr.Motion!='fixed':raise ValueError('Unknown motion model')
        if tr.Type=='ecef':frame='ECEF'
    elif tr.Type=='keplerian':
        a,e,inc,node,arg,nu=tr.Elements;nu*=D
        E0=2*np.arctan2(np.sqrt(1-e)*np.sin(nu/2),np.sqrt(1+e)*np.cos(nu/2));M=E0-e*np.sin(E0)+np.sqrt(mu/a**3)*dt;E=M.copy()
        for _ in range(40):
            correction=(E-e*np.sin(E)-M)/(1-e*np.cos(E));E-=correction
            if np.max(abs(correction))<1e-13:break
        q=np.sqrt(1-e*e);edot=np.sqrt(mu/a**3)/(1-e*np.cos(E));Q=rz(node)@rx(inc)@rz(arg)
        r=Q@np.array([a*(np.cos(E)-e),a*q*np.sin(E),np.zeros(n)]);v=Q@np.array([-a*np.sin(E)*edot,a*q*np.cos(E)*edot,np.zeros(n)]);frame='GCRF'
    elif tr.Type=='ephemeris':
        x=offsets(c,tr.TimeUTC);p=np.asarray(tr.Positions_m,float)
        if min(t)<x[0]-1e-7 or max(t)>x[-1]+1e-7:raise ValueError('Ephemeris does not cover exposure')
        if tr.Interpolation=='hermite' and len(tr.Velocities_mps):
            f=CubicHermiteSpline(x,p,np.asarray(tr.Velocities_mps),axis=1);r=f(t);v=f(t,1)
        else:
            i=np.clip(np.searchsorted(x,t,side='right')-1,0,len(x)-2);u=(t-x[i])/(x[i+1]-x[i]);v=(p[:,i+1]-p[:,i])/(x[i+1]-x[i]);r=p[:,i]+u*(p[:,i+1]-p[:,i])
    elif tr.Type=='tle':
        from astropy.time import Time
        from astropy.coordinates import TEME,GCRS,CartesianRepresentation,CartesianDifferential
        from astropy import units as u
        from astropy.utils import iers
        iers.conf.auto_download=False
        j=jd(c,t);i=np.floor(j);err,rk,vk=tr._sat.sgp4_array(i,j-i)
        if np.any(err):raise ValueError(f'SGP4 error: {err}')
        at=Time(j,format='jd',scale='utc');at.delta_ut1_utc=c.Frames.DeltaUT1_s
        rep=CartesianRepresentation(rk.T*u.km,differentials=CartesianDifferential(vk.T*u.km/u.s))
        g=TEME(rep,obstime=at).transform_to(GCRS(obstime=at));r=g.cartesian.xyz.to_value(u.m);v=g.velocity.d_xyz.to_value(u.m/u.s);frame='GCRF'
    else:raise ValueError(f'Unsupported trajectory: {tr.Type}')
    R=rotation(c,t)
    if frame=='ECEF':
        re=r.copy();v=np.einsum('nij,jn->in',R,v)+np.einsum('nij,jn->in',rotation(c,t+.5)-rotation(c,t-.5),r);r=np.einsum('nij,jn->in',R,r)
    elif frame=='GCRF':re=np.einsum('nji,jn->in',R,r)
    else:raise ValueError('State frame must be GCRF or ECEF')
    return struct(dict(Position_m=r,Velocity_mps=v,ECEF_m=re,Geodetic_deg_m=geodetic(c,re,True),Frame='GCRF'))

def rotation_series(times,matrices,query):
    x=np.asarray(times);q=np.asarray(query);a=np.asarray(matrices,float)
    if a.shape==(3,3,len(x)):a=np.moveaxis(a,2,0) # MATLAB layout
    if len(x)<2 or np.any(np.diff(x)<=0) or np.min(q)<x[0]-1e-7 or np.max(q)>x[-1]+1e-7:raise ValueError('Attitude series does not cover exposure')
    return Slerp(x,Rotation.from_matrix(a))(np.clip(q,x[0],x[-1])).as_matrix()
def pointing(c,o,t,target=None):
    t=np.atleast_1d(t);P=c.Pointing
    if P.Mode=='attitude':A=rotation_series(offsets(c,P.AttitudeUTC),P.CameraToGCRF,t)
    else:
        R=rotation(c,t);A=[]
        for k in range(len(t)):
            E=R[k]@enu(*o.Geodetic_deg_m[:2,k])
            if P.Mode=='az_el_roll':
                az,el=np.asarray(P.AzEl_deg)*D;b=E@np.array([np.cos(el)*np.sin(az),np.cos(el)*np.cos(az),np.sin(el)]);a=E[:,2]
            elif P.Mode=='nadir_offset':
                nadir=-unit(o.Position_m[:,k]);normal=unit(np.cross(o.Position_m[:,k],o.Velocity_mps[:,k]));along=unit(np.cross(normal,-nadir));off=P.OffNadir_deg*D;az=P.SpaceAzimuth_deg*D
                b=np.cos(off)*nadir+np.sin(off)*(np.cos(az)*along+np.sin(az)*normal);a=normal
            elif P.Mode=='inertial':b=unit(P.Boresight_GCRF);a=unit(P.ReferenceUp_GCRF)
            elif P.Mode=='target_tracking':b=unit(target[:,k]-o.Position_m[:,k]);a=unit(P.ReferenceUp_GCRF)
            else:raise ValueError('Unsupported pointing mode')
            if np.linalg.norm(np.cross(a,b))<1e-8:
                a=E[:,1]
                if np.linalg.norm(np.cross(a,b))<1e-8:a=np.eye(3)[:,np.argmin(abs(b))]
            v=unit(-a+np.dot(a,b)*b);u=unit(np.cross(v,b));angle=P.Roll_deg*D
            A.append(np.column_stack([np.cos(angle)*u+np.sin(angle)*v,-np.sin(angle)*u+np.cos(angle)*v,b]))
        A=np.array(A)
    if P.JitterUTC:
        x=offsets(c,P.JitterUTC)
        if min(t)<x[0] or max(t)>x[-1]:raise ValueError('Jitter lacks time support')
        j=np.array([np.interp(t,x,np.asarray(P.Jitter_arcsec)[:,k]) for k in range(3)]).T/3600
        for k,(a,b,g) in enumerate(j):A[k]=A[k]@ry(a)@rx(-b)@rz(g)
    return A

def project(c,A,los):
    d=A.T@np.asarray(los).reshape(3,-1);xi=np.arctan2(d[0],d[2])/D;eta=np.arctan2(-d[1],d[2])/D
    f=np.asarray(c.Camera.FOV_deg);inside=(d[2]>0)&(abs(xi)<=f[0]/2+1e-10)&(abs(eta)<=f[1]/2+1e-10)
    if c.Camera.Shape=='circular':inside&=np.arccos(np.clip(d[2],-1,1))/D<=f[0]/2+1e-10
    return struct(dict(x=(c.Camera.Nx+1)/2+(c.Camera.Nx-1)*xi/f[0],y=(c.Camera.Ny+1)/2-(c.Camera.Ny-1)*eta/f[1],Xi_deg=xi,Eta_deg=eta,InFOV=inside))

def seed_targets(c):
    targets=[T for T in c.Targets if T.Trajectory.Type=='seeded_state']
    if not targets:return c
    if c.Pointing.Mode=='target_tracking':raise ValueError('Resolve seeds before tracking')
    o=state(c,c.Observer.Trajectory,[0]);A=pointing(c,o,[0])[0];R=rotation(c,[0])[0]
    for T in targets:
        s=T.Seed;a,b=s.Angles_deg
        if s.AngleMode=='az_el_offset':
            E=R@enu(*o.Geodetic_deg_m[:2,0]);az,el=(np.asarray(c.Pointing.AzEl_deg)+[a,b])*D
            p=project(c,A,E@np.array([np.cos(el)*np.sin(az),np.cos(el)*np.cos(az),np.sin(el)]));xy=np.array([p.Xi_deg[0],p.Eta_deg[0]])
        elif s.AngleMode=='theta_phi':xy=np.arctan(np.tan(a*D)*np.array([np.cos(b*D),np.sin(b*D)]))/D
        else:xy=np.array([a,b])
        q=np.r_[np.tan(xy[0]*D),-np.tan(xy[1]*D),1];u=unit(q)
        v=np.array([np.cos(s.Direction_deg*D)/np.cos(xy[0]*D)**2,-np.sin(s.Direction_deg*D)/np.cos(xy[1]*D)**2,0]);v=unit(v-u*np.dot(u,v))
        tr=T.Trajectory;tr.Type='state';tr.Frame='GCRF';tr.Motion=c.Seeding.Propagation;tr.EpochUTC=c.Time.StartUTC
        tr.Position_m=(o.Position_m[:,0]+s.Range_m*A@u).tolist();tr.Velocity_mps=(o.Velocity_mps[:,0]+A@(s.TransverseSpeed_mps*v+s.RadialSpeed_mps*u)).tolist();T.Metadata.Seed=deepcopy(s)
    return c

def illumination(c,r,sunpos):
    re=c.Constants.EarthRadius_m;rr=np.linalg.norm(r,axis=0);s=unit(sunpos);behind=np.sum(r*s,axis=0)<0;dist=np.linalg.norm(r-s*np.sum(r*s,axis=0),axis=0)
    I=(~(behind&(dist<re)) & (rr>re)).astype(float)
    if c.Visibility.EclipseModel=='cylinder':return I
    for k in range(r.shape[1]):
        if rr[k]<=re:continue
        tosun=sunpos[:,k]-r[:,k];ds=np.linalg.norm(tosun);a=np.arcsin(min(1,c.Constants.SunRadius_m/ds));b=np.arcsin(re/rr[k]);d=np.arccos(np.clip(np.dot(-r[:,k]/rr[k],tosun/ds),-1,1))
        if d>=a+b:overlap=0
        elif d<=abs(b-a):overlap=np.pi*min(a,b)**2
        else:
            x=np.clip((d*d+a*a-b*b)/(2*d*a),-1,1);y=np.clip((d*d+b*b-a*a)/(2*d*b),-1,1)
            overlap=a*a*np.arccos(x)+b*b*np.arccos(y)-.5*np.sqrt(max(0,(-d+a+b)*(d+a-b)*(d-a+b)*(d+a+b)))
        I[k]=np.clip(1-overlap/(np.pi*a*a),0,1)
    return I

def sample(c,t,mid):
    t=np.atleast_1d(t);o=state(c,c.Observer.Trajectory,t);R=rotation(c,t);up=np.column_stack([R[k]@enu(*o.Geodetic_deg_m[:2,k])[:,2] for k in range(len(t))])
    states=[state(c,T.Trajectory,t) if T.Trajectory.Type!='demo' else None for T in c.Targets]
    tracked=next((s.Position_m for T,s in zip(c.Targets,states) if T.ID==c.Pointing.TrackTargetID and s is not None),None)
    A=pointing(c,o,t,tracked)
    if c.Integration.AttitudeMode=='midpoint':
        om=state(c,c.Observer.Trajectory,[mid]);track=None
        if c.Pointing.Mode=='target_tracking':track=state(c,next(T.Trajectory for T in c.Targets if T.ID==c.Pointing.TrackTargetID),[mid]).Position_m
        A=np.repeat(pointing(c,om,[mid],track),len(t),axis=0)
    sp=sun(c,t);ground=o.Geodetic_deg_m[2]<c.Visibility.PlatformAltitude_m;out=[]
    for T,S in zip(c.Targets,states):
        if S is None:
            d=T.Trajectory.Demo;r=[];dt=t-(utc(d.ReferenceUTC)-utc(c.Time.StartUTC)).total_seconds()
            for k in range(len(t)):
                off=np.asarray(d.Offset_deg)
                if d.AngleMode=='az_el_offset':
                    E=R[k]@enu(*o.Geodetic_deg_m[:2,k]);az,el=(np.asarray(c.Pointing.AzEl_deg)+off)*D;p=project(c,A[k],E@np.array([np.cos(el)*np.sin(az),np.cos(el)*np.cos(az),np.sin(el)]));xy=np.r_[p.Xi_deg,p.Eta_deg]
                elif d.AngleMode=='theta_phi':xy=np.arctan(np.tan(off[0]*D)*np.array([np.cos(off[1]*D),np.sin(off[1]*D)]))/D
                else:xy=off.copy()
                rates=np.asarray(d.AngularRate_deg_s)
                if d.MotionModel=='linear_angles':
                    xy=xy+rates*dt[k]
                    if np.any(abs(xy)>=89.9):raise ValueError('Demo crossed forward hemisphere')
                    local=unit(np.r_[np.tan(xy[0]*D),-np.tan(xy[1]*D),1])
                else:
                    raw=np.r_[np.tan(xy[0]*D),-np.tan(xy[1]*D),1];local=unit(raw);dv=np.r_[rates[0]/np.cos(xy[0]*D)**2,-rates[1]/np.cos(xy[1]*D)**2,0]*D;dv=(dv-local*np.dot(local,dv))/np.linalg.norm(raw);speed=np.linalg.norm(dv)
                    if speed>0:local=local*np.cos(speed*dt[k])+dv/speed*np.sin(speed*dt[k])
                r.append(o.Position_m[:,k]+d.Range_m*A[k]@local)
            r=np.array(r).T;vel=np.gradient(r,t,axis=1) if len(t)>1 else np.full_like(r,np.nan)
        else:r=S.Position_m;vel=S.Velocity_mps
        delta=r-o.Position_m;ran=np.linalg.norm(delta,axis=0);u=delta/np.maximum(ran,1e-300)
        projections=[project(c,a,u[:,k]) for k,a in enumerate(A)]
        h=struct({name:np.array([p[name][0] for p in projections]) for name in ('x','y','Xi_deg','Eta_deg','InFOV')})
        h.update(ID=T.ID,Position_m=r,Velocity_mps=vel,Range_m=ran,LOS_GCRF=u)
        h.Elevation_deg=np.arcsin(np.clip(np.sum(up*u,axis=0),-1,1))/D
        closest=o.Position_m+u*np.clip(-np.sum(o.Position_m*u,axis=0),0,ran)
        h.EarthOccluded=(np.linalg.norm(closest,axis=0)<c.Constants.EarthRadius_m)&~ground if c.Visibility.EarthOccultation else np.zeros(len(t),bool)
        h.PlatformVisible=(~ground|(h.Elevation_deg>=c.Visibility.MinElevation_deg))&~h.EarthOccluded
        h.RangeValid=(ran>=c.Visibility.RangeLimits_m[0])&(ran<=c.Visibility.RangeLimits_m[1]);su=unit(sp-o.Position_m)
        h.SunExcluded=np.arccos(np.clip(np.sum(u*su,axis=0),-1,1))/D<c.Visibility.SunExclusion_deg
        h.Visible=h.InFOV&h.PlatformVisible&h.RangeValid&~h.SunExcluded;h.Illumination=illumination(c,r,sp)
        h.PhaseAngle_rad=np.arccos(np.clip(np.sum(unit(sp-r)*(-u),axis=0),-1,1));out.append(h)
    return struct(dict(Time_s=t,TimeUTC=iso(c,t),Observer=o,CameraToGCRF=A,Zenith_GCRF=up,IsGround=ground,SunPosition_m=sp,Target=out))

def trajectory_example(c,name):
    if name in ('trajectory_geodetic','trajectory_ecef'):
        c.Observer.Trajectory.Type='geodetic';c.Observer.Trajectory.Frame='ECEF';c.Targets=c.Targets[:1];c.Targets[0].Trajectory=deepcopy(c.Observer.Trajectory);c.Targets[0].Trajectory.Geodetic_deg_m=[-37.8,144.96,800e3]
        if name.endswith('ecef'):
            for tr in [c.Observer.Trajectory,c.Targets[0].Trajectory]:tr.Position_m=geodetic(c,tr.Geodetic_deg_m)[:,0].tolist();tr.Type='ecef'
    elif name=='trajectory_state':
        for tr in [c.Observer.Trajectory]+[T.Trajectory for T in c.Targets]:
            s=state(c,tr,[0]);tr.Type='state';tr.Motion='two_body';tr.Position_m=s.Position_m[:,0].tolist();tr.Velocity_mps=s.Velocity_mps[:,0].tolist()
    elif name in ('trajectory_ephemeris','trajectory_ephemeris_ecef'):
        for tr in [c.Observer.Trajectory]+[T.Trajectory for T in c.Targets]:
            t=np.array([-1,0,.1,.2,1]);s=state(c,tr,t);tr.Type='ephemeris';tr.TimeUTC=iso(c,t);tr.File='';tr.Frame='GCRF';tr.Positions_m=s.Position_m.tolist();tr.Velocities_mps=s.Velocity_mps.tolist()
            if name.endswith('ecef'):tr.Frame='ECEF';tr.Positions_m=s.ECEF_m.tolist();tr.Velocities_mps=[];tr.Interpolation='linear'
    elif name=='trajectory_ephemeris_csv':
        c.Targets=c.Targets[:1]
        for tr,kind in [(c.Observer.Trajectory,'observer'),(c.Targets[0].Trajectory,'target')]:tr.Type='ephemeris';tr.Frame='GCRF';tr.File=str(__import__('pathlib').Path(c.RootDir)/f'Inputs/Ephemeris/{kind}_gcrf.csv')
    elif name.startswith('trajectory_tle'):
        c.Targets=c.Targets[:1];prefix='geo' if name.endswith('deep_space') else 'leo'
        for tr,kind in [(c.Observer.Trajectory,'observer'),(c.Targets[0].Trajectory,'target')]:tr.Type='tle';tr.File=str(__import__('pathlib').Path(c.RootDir)/f'Inputs/TLE/{prefix}_{kind}.tle')
    return c
