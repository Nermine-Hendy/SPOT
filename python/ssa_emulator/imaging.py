"""Detector-aligned rendering, linear optics, detector physics, annotation masks."""
import numpy as np
from scipy.ndimage import convolve1d
from .config import struct
from .photometry import rng
from .geometry import D

def detector_grid(c,A):
    nx,ny=c.Camera.Nx,c.Camera.Ny;f=np.asarray(c.Camera.FOV_deg)
    xi,eta=np.meshgrid(np.linspace(-f[0]/2,f[0]/2,nx),np.linspace(f[1]/2,-f[1]/2,ny));v=np.array([np.tan(xi.ravel(order='F')*D),-np.tan(eta.ravel(order='F')*D),np.ones(nx*ny)]);v/=np.linalg.norm(v,axis=0)
    active=np.ones((ny,nx),bool)
    if c.Camera.Shape=='circular':active=(np.arccos(np.clip(v[2],-1,1))/D<=f[0]/2+1e-10).reshape((ny,nx),order='F')
    if not active.any():raise ValueError('No active pixel centres')
    omega=2*np.pi*(1-np.cos(f[0]*D/2)) if c.Camera.Shape=='circular' else 4*np.arcsin(np.sin(f[0]*D/2)*np.sin(f[1]*D/2))
    if c.Camera.SolidAngle=='mean':w=active*omega/active.sum()
    else:
        tx,ty=np.tan(xi*D),np.tan(eta*D);w=(1+tx*tx)*(1+ty*ty)/(1+tx*tx+ty*ty)**1.5;w[[0,-1]]/=2;w[:,[0,-1]]/=2;w*=active;w*=omega/w.sum()
    los=A@v;shape=(ny,nx)
    out=struct(dict(Xi_deg=xi,Eta_deg=eta,Active=active,SolidAngle_sr=w,FOVSolidAngle_sr=omega,LOS_GCRF=los,RA_deg=(np.arctan2(los[1],los[0])/D%360).reshape(shape,order='F'),DEC_deg=(np.arcsin(np.clip(los[2],-1,1))/D).reshape(shape,order='F')))
    step=c.Camera.DiagnosticGridStep_deg
    gx=np.unique(np.r_[np.arange(-f[0]/2,f[0]/2,step),f[0]/2]);gy=np.unique(np.r_[np.arange(-f[1]/2,f[1]/2,step),f[1]/2])[::-1];xx,yy=np.meshgrid(gx,gy)
    v=np.array([np.tan(xx.ravel(order='F')*D),-np.tan(yy.ravel(order='F')*D),np.ones(xx.size)]);v/=np.linalg.norm(v,axis=0);u=A@v
    out.AngularDiagnosticGrid=struct(dict(Xi_deg=xx,Eta_deg=yy,RA_deg=(np.arctan2(u[1],u[0])/D%360).reshape(xx.shape,order='F'),DEC_deg=(np.arcsin(np.clip(u[2],-1,1))/D).reshape(xx.shape,order='F')))
    return out

def deposit(c,x,y,value,active):
    shape=(c.Camera.Ny,c.Camera.Nx);a=np.zeros(shape);x,y,value=np.broadcast_arrays(np.asarray(x,float).ravel(),np.asarray(y,float).ravel(),np.asarray(value,float).ravel())
    use=np.isfinite(x)&np.isfinite(y)&np.isfinite(value)&(value!=0)&(x>=1-1e-8)&(x<=shape[1]+1e-8)&(y>=1-1e-8)&(y<=shape[0]+1e-8)
    x=np.clip(x[use],1,shape[1])-1;y=np.clip(y[use],1,shape[0])-1;value=value[use]
    if not len(x):return a
    x0=np.floor(x).astype(int);y0=np.floor(y).astype(int);fx=x-x0;fy=y-y0
    xx=np.array([x0,x0+1,x0,x0+1]).T;yy=np.array([y0,y0,y0+1,y0+1]).T;w=np.array([(1-fx)*(1-fy),fx*(1-fy),(1-fx)*fy,fx*fy]).T
    valid=(xx>=0)&(xx<shape[1])&(yy>=0)&(yy<shape[0]);xx=np.clip(xx,0,shape[1]-1);yy=np.clip(yy,0,shape[0]-1);w*=valid*active[yy,xx];w/=np.maximum(w.sum(axis=1,keepdims=True),np.finfo(float).tiny)
    np.add.at(a,(yy.ravel(),xx.ravel()),(w*value[:,None]).ravel());return a

def blur(a,sigma):
    if sigma<=1e-12:return a.copy()
    q=np.arange(-np.ceil(4*sigma),np.ceil(4*sigma)+1);h=np.exp(-.5*(q/sigma)**2);h/=h.sum()
    return convolve1d(convolve1d(a,h,axis=0,mode='constant',cval=0),h,axis=1,mode='constant',cval=0)
def optics(c,rate,active):
    yy,xx=np.indices(rate.shape);r=np.hypot(xx-(c.Camera.Nx-1)/2,yy-(c.Camera.Ny-1)/2);r=np.minimum(1,r/r[active].max());tau=1;vig=np.ones_like(rate);o=c.Optics
    if o.Enabled:
        tau=np.prod(o.ElementTransmission) if len(o.ElementTransmission) else o.Transmission
        if o.VignettingEnabled:vig=1-(1-o.EdgeTransmission)*r**o.RadialExponent
    integrated=rate*c.Time.Exposure_s;through=integrated*tau;vignetted=through*vig;out=vignetted.copy()
    if o.Enabled and o.PSFEnabled:
        if o.PSFModel=='gaussian' or o.SigmaEdgeAdd_pix==0:out=blur(out,o.SigmaCenter_pix)
        else:
            out=np.zeros_like(rate)
            for node in np.linspace(0,1,o.PSFBanks):
                w=np.maximum(0,1-abs(r-node)*(o.PSFBanks-1));out+=blur(vignetted*w,o.SigmaCenter_pix+o.SigmaEdgeAdd_pix*node)
    out=np.maximum(0,out)*active
    return out,struct(dict(Integrated=integrated,Throughput=through,Vignetted=vignetted,PSF=out,Vignetting=vig,ThroughputFactor=tau))

def detector(c,signal,exposure_index,active):
    s=c.Detector;seed=(c.RandomSeed+700001+7919*exposure_index)%2**32;expected=signal*s.RadiometricGain_e_per_unit
    if len(s.PRNUMap):expected*=np.asarray(s.PRNUMap)
    dark=np.full(signal.shape,s.DarkCurrent_e_s*c.Time.Exposure_s if s.DarkEnabled else 0)
    photo=rng(seed).poisson(expected).astype(float) if s.PoissonEnabled else expected.copy()
    de=rng(seed+1000003).poisson(dark).astype(float) if s.PoissonEnabled else dark.copy();before=photo+de;sat=(before>s.FullWell_e)&s.FullWellEnabled;after=np.clip(before,0,s.FullWell_e) if s.FullWellEnabled else before.copy()
    noise=s.ReadNoise_e*rng(seed+2000003).standard_normal(signal.shape) if s.ReadEnabled else np.zeros_like(signal);read=after+noise;analog=read/s.ConversionGain_e_per_DN+(s.Bias_DN if s.BiasEnabled else 0)
    adcsat=(analog>2**s.ADCBits-1)&s.ADCEnabled
    # MATLAB round ties away from zero, and ADC input has already been clipped >=0.
    image=np.floor(np.clip(analog,0,2**s.ADCBits-1)+.5) if s.ADCEnabled else analog.copy();valid=active.copy()
    if len(s.BadPixelMask):valid&=~np.asarray(s.BadPixelMask,bool)
    return struct(dict(ExpectedSignal_e=expected,ExpectedDark_e=dark,SignalElectrons=photo,DarkElectrons=de,BeforeFullWell_e=before,AfterFullWell_e=after,FullWellSaturated=sat,ReadNoise_e=noise,AfterRead_e=read,Analog_DN=analog,ADCSaturated=adcsat,Image_DN=image,ValidPixel=valid))

def morphology(c,h,w):
    p=np.column_stack([h.x,h.y]);a=w*h.RadiometricRate*h.Visible;ok=np.isfinite(p).all(axis=1)&np.isfinite(a)&(a>0)
    m=struct(dict(Kind='unobserved',Centroid_xy=[np.nan,np.nan],Covariance_pix2=np.zeros((2,2)),EquivalentSmear_pix=0.,PathLength_pix=0.,EndToEnd_pix=0.,PSFFWHM_pix=0.,Threshold_pix=0.,DirectionDetector_deg=np.nan))
    if not ok.any():return m
    q=p[ok];a=a[ok]/a[ok].sum();centre=a@q;d=q-centre;V=d.T@(a[:,None]*d);L,U=np.linalg.eigh(V);smear=np.sqrt(12*max(0,L[-1]));path=np.linalg.norm(np.diff(p,axis=0),axis=1);path=path[ok[:-1]&ok[1:]].sum()
    radius=np.linalg.norm(q-[(c.Camera.Nx+1)/2,(c.Camera.Ny+1)/2],axis=1)
    if c.Camera.Shape=='circular':
        grid=detector_grid(c,np.eye(3));yy,xx=np.indices(grid.Active.shape);rmax=np.hypot(xx-(c.Camera.Nx-1)/2,yy-(c.Camera.Ny-1)/2)[grid.Active].max()
    else:rmax=np.hypot((c.Camera.Nx-1)/2,(c.Camera.Ny-1)/2)
    sigma=np.zeros(len(q))
    if c.Optics.Enabled and c.Optics.PSFEnabled:
        sigma[:]=c.Optics.SigmaCenter_pix
        if c.Optics.PSFModel=='field_dependent_gaussian':sigma+=c.Optics.SigmaEdgeAdd_pix*np.minimum(1,radius/rmax)
    fwhm=2*np.sqrt(2*np.log(2))*np.sqrt(a@(sigma*sigma));threshold=max(c.Morphology.MinSmear_pix,c.Morphology.PSFFWHMFactor*fwhm)
    m.update(Kind='trail' if smear>threshold else 'object',Centroid_xy=centre,Covariance_pix2=V,EquivalentSmear_pix=smear,PathLength_pix=path,EndToEnd_pix=np.linalg.norm(q[-1]-q[0]),PSFFWHM_pix=fwhm,Threshold_pix=threshold,DirectionDetector_deg=np.arctan2(U[1,-1],U[0,-1])/D%180)
    return m

def labels(c,maps,records):
    shape=(c.Camera.Ny,c.Camera.Nx);instance=np.zeros(shape,np.uint32);semantic=np.zeros(shape,np.uint8);best=np.zeros(shape);binary=np.zeros(shape,bool)
    objects=np.zeros(shape,bool);trails=np.zeros(shape,bool);per_target=[]
    for j,(a,r) in enumerate(zip(maps,records)):
        peak=float(a.max());threshold=c.Masks.RelativeThreshold*peak;mask=(a>=threshold)&(a>0) if peak>0 else np.zeros(shape,bool);binary|=mask
        per_target.append(mask)
        if r.Appearance=='object':objects|=mask
        else:trails|=mask
        win=mask&(a>best);best[win]=a[win];instance[win]=r.ID;semantic[win]=1 if r.Appearance=='object' else 2
        yy,xx=np.where(mask);box=[int(xx.min()+1),int(yy.min()+1),int(xx.max()-xx.min()+1),int(yy.max()-yy.min()+1)] if len(xx) else [np.nan,np.nan,0,0]
        r.update(HasAnnotation=bool(mask.any()),PeakOpticalSignal=peak,MaskPixelCount=int(mask.sum()),MaskThreshold=c.Masks.RelativeThreshold,AbsoluteMaskThreshold=threshold,BoundingBox_xywh=box,BoxCenter_xy=[box[0]+(box[2]-1)/2,box[1]+(box[3]-1)/2] if len(xx) else [np.nan,np.nan])
    return struct(dict(Binary=binary,Instance=instance,Semantic=semantic,ObjectUnion=objects,TrailUnion=trails,PerTarget=per_target,Schema=dict(Background=0,Object=1,Trail=2,InstanceMeaning='Target ID',PixelCoordinates='one-based pixel centres; boxes [xmin ymin inclusive_width inclusive_height]',OverlapRule='largest target-only optical signal among candidate masks; ties smallest target array index')))
